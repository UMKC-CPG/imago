## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

"""kaleidoscope.wingbeats -- the pluggable wingbeat seam
(DESIGN 6.2.2; PSEUDOCODE 13.2).

A *wingbeat* is the seam (VISION Principle 8) between
kaleidoscope's dispatch core and how a unit actually executes.
It takes a unit and its prepared run directory and returns a
domain-agnostic ``WingbeatOutcome``; the dispatch core never changes
when a new wingbeat is added.

The default wingbeat, ``ImagoWingbeat``, drives the imago.py
callable API (DESIGN 6.1).  An ASE wingbeat (D12) and future
adapters implement the same one-method protocol.  Wingbeats are
looked up by name through the ``WINGBEATS`` registry so a unit can
name its wingbeat as a string and a worker process can resolve it
after import.
"""

import os
import shutil

from .model import WingbeatOutcome, KaleidoscopeError
from .workspace import emit_scalar, toml_line


def _partition_options(options):
    """Split a unit's options into the makeinput build set and the
    imago run set (DESIGN 6.2.10).  Routing is by each tool's
    recognised-key set, TWO buckets and no third:

    - a key in ``imago.OPTION_KEYS`` (job, edge, scf_basis, ...) is
      an imago run-time selection -> the imago set;
    - every other key -> the makeinput build set, where makeinput's
      strict unknown-key check stays the typo backstop.

    There is deliberately no "dropped before forwarding" bucket.
    ``options`` is a dictionary of *tool inputs*, so a fact that
    reaches neither tool has no business in it; bookkeeping such as
    the engine build identity rides on ``unit.record`` instead
    (DESIGN 6.2.4).  That is what keeps makeinput's strictness
    meaningful here -- an unrecognised key is now always a typo,
    never a deliberate passenger.

    Returns ``(makeinput_options, imago_options)``."""
    import imago
    makeinput_options = {}
    imago_options = {}
    for key, value in options.items():
        if key in imago.OPTION_KEYS:
            imago_options[key] = value
        else:
            makeinput_options[key] = value
    return makeinput_options, imago_options


# ------------------------------------------------------------------
#  Wingbeat protocol and registry
# ------------------------------------------------------------------

class Wingbeat:
    """Base class documenting the wingbeat protocol.  A wingbeat
    implements ``run(unit, wingbeat_dir) -> WingbeatOutcome``; ``run``
    executes the calculation however it likes and reports a
    generic outcome.  ``ok`` means the unit *completed*, not
    that it succeeded scientifically; ``detail`` is an opaque
    string kaleidoscope records but never interprets."""

    def run(self, unit, wingbeat_dir):
        raise NotImplementedError


# name -> Wingbeat instance.  The default wingbeat is registered at
#   import (below), so a freshly imported worker process can
#   resolve "imago" without any flight-side setup.
WINGBEATS = {}


def register_wingbeat(name, wingbeat):
    """Register a wingbeat instance under a name (DESIGN 6.2.2).
    A custom wingbeat must be registered at import time in any
    process that will execute it -- the main process always, and
    each worker process when a multi-process executor is used."""
    WINGBEATS[name] = wingbeat


def resolve_wingbeat(name):
    """Return the registered wingbeat for ``name`` or raise."""
    if name not in WINGBEATS:
        raise KaleidoscopeError(
            f"no wingbeat registered under name {name!r}; "
            f"known wingbeats: {sorted(WINGBEATS)}"
        )
    return WINGBEATS[name]


# ------------------------------------------------------------------
#  The default Imago wingbeat
# ------------------------------------------------------------------

class ImagoWingbeat(Wingbeat):
    """Run a unit through the imago.py callable API (DESIGN
    6.2.2).  When the run directory already holds staged inputs
    it is run as a prepared directory; otherwise it is built from
    the unit's structure and options.  The native ``ImagoResult``
    is persisted into the run directory as ``result.toml`` so the
    client can reload it during harvest -- kaleidoscope itself
    never reads it (VISION Principle 9)."""

    def run(self, unit, wingbeat_dir):
        # Imported lazily so the package imports without imago's
        #   own runtime environment ($IMAGO_RC etc.) being set.
        import imago

        # The imago-side options are RUN-TIME settings (job, edge,
        #   scf_basis, ...) that do NOT live in a staged imago.dat
        #   (DESIGN 6.2.10), so they must be re-applied on EVERY
        #   launch.  The job type and the SCF suppression live only
        #   in these settings, so if they are dropped imago no
        #   longer sees the unit's `-loen -scf no` request and falls
        #   back to its DEFAULT job -- a ground-state SCF.
        #   (`-loen -scf no` never runs an SCF itself; the unwanted
        #   SCF is purely the dropped-settings fallback -- the "SCF
        #   after loen" the seed run hit.)  Build the settings once
        #   and pass them however the inputs get staged (6.2.2/6.1).
        _, imago_options = _partition_options(unit.options)
        settings = imago.ScriptSettings.from_options(imago_options)

        self._stage_inputs(unit, wingbeat_dir)
        result = imago.run_prepared(wingbeat_dir, settings=settings)

        self._persist_result(wingbeat_dir, result, unit.record)

        # Map the Imago-native status onto the generic outcome:
        #   "ran" covers CONVERGED / NOT_CONVERGED / SKIPPED;
        #   only a hard FAILED is not-ok.  The status value
        #   becomes the opaque detail string the flight records.
        ok = result.status in (
            imago.RunStatus.CONVERGED,
            imago.RunStatus.NOT_CONVERGED,
            imago.RunStatus.SKIPPED,
        )
        return WingbeatOutcome(
            ok=ok,
            detail=result.status.value,
            runtime_seconds=result.runtime_seconds,
            message=result.message,
        )

    @staticmethod
    def _stage_inputs(unit, wingbeat_dir):
        """Ensure the run directory holds runnable inputs, so
        ``run_prepared`` can execute them (DESIGN 6.2.5; PSEUDOCODE
        13.2).  Three cases:

        - the driver's prepare step already built this unit's
          inputs into ``unit.prepared_dir`` -- commit that staged
          copy into the run directory (the Model-A producer path);
        - the run directory already holds a staged ``imago.dat`` --
          a re-run of a directory a prior launch built -- so there
          is nothing to do;
        - neither (a client that did not prepare) -- build the deck
          from the unit's structure and its makeinput-side options.

        The wingbeat owns the makeinput/imago option split (DESIGN
        6.2.10) only on the build path: a unit carries ONE options
        dict, but makeinput (strict) and imago (lenient) have
        disjoint key vocabularies, so the wingbeat routes each key
        to the tool that recognises it and drops the cache-only
        build identity.
        """
        if unit.prepared_dir is not None:
            ImagoWingbeat._commit_prepared_inputs(
                unit.prepared_dir, wingbeat_dir)
        elif not ImagoWingbeat._is_prepared(wingbeat_dir):
            import makeinput
            makeinput_options, _ = _partition_options(unit.options)
            makeinput.build_run_dir(
                unit.structure, makeinput_options, wingbeat_dir)

    @staticmethod
    def _commit_prepared_inputs(prepared_dir, wingbeat_dir):
        """Copy the driver-staged inputs (structure.dat, imago.dat,
        scfV, kp files -- DESIGN 6.2.5) from ``prepared_dir`` into
        the run directory so ``run_prepared`` finds them.  The
        staging area is transient (the producer's prepare pass
        rebuilds it each run), so the commit simply copies from it,
        merging into any existing run-directory contents."""
        os.makedirs(wingbeat_dir, exist_ok=True)
        for name in os.listdir(prepared_dir):
            source = os.path.join(prepared_dir, name)
            target = os.path.join(wingbeat_dir, name)
            if os.path.isdir(source):
                shutil.copytree(source, target, dirs_exist_ok=True)
            else:
                shutil.copy2(source, target)

    @staticmethod
    def _is_prepared(wingbeat_dir):
        """A run directory is 'prepared' when it already holds
        the primary imago.dat (directly or under inputs/), so it
        can be run as-is without a makeinput build."""
        for candidate in ("imago.dat",
                          os.path.join("inputs", "imago.dat")):
            if os.path.exists(os.path.join(wingbeat_dir, candidate)):
                return True
        return False

    @staticmethod
    def _persist_result(wingbeat_dir, result, record=None):
        """Write the ImagoResult to ``<wingbeat_dir>/result.toml`` for
        the client's harvest.  Flat scalar fields first, then an
        ``[outputs]`` table of logical-key -> path and a ``[job]``
        table echoing what ran.

        One *recorded* fact rides along with the measured ones: the
        engine build identity out of the unit's ``record`` mapping
        (DESIGN 6.2.4), written as ``imago_commit``.  A guidance
        entry's provenance reads it here (DESIGN 7.8), which keeps
        that harvest on the three per-run sources it already has and
        off the dispatch core's ``status.toml``.  The engine's own
        word wins when it has one: imago does not report its build
        yet (TODO C84), and when it does, preferring
        ``result.imago_commit`` over the recorded value is the whole
        change -- one substitution in one field of one file."""
        os.makedirs(wingbeat_dir, exist_ok=True)
        path = os.path.join(wingbeat_dir, "result.toml")
        build = getattr(result, "imago_commit", None) or (
            (record or {}).get("imago_commit"))
        with open(path, "w") as result_file:
            result_file.write(
                toml_line("status", result.status.value)
            )
            result_file.write(toml_line("imago_commit", build))
            result_file.write(toml_line("success", result.success))
            result_file.write(toml_line("converged",
                                       result.converged))
            result_file.write(
                toml_line("reused_checkpoint",
                          result.reused_checkpoint)
            )
            result_file.write(
                toml_line("scf_iterations", result.scf_iterations)
            )
            result_file.write(
                toml_line("total_energy", result.total_energy)
            )
            result_file.write(
                toml_line("total_magnetization",
                          result.total_magnetization)
            )
            result_file.write(toml_line("gap_ev", result.gap_ev))
            result_file.write(toml_line("gap_kind", result.gap_kind))
            result_file.write(
                toml_line("scf_threshold", result.scf_threshold)
            )
            # Resolved mesh (DESIGN 6.1.2): kpoint_mesh renders as
            #   a TOML array, kpoint_count as a scalar; toml_line
            #   omits either when None (an explicit-list run or an
            #   older imago binary that emits neither).
            result_file.write(
                toml_line("kpoint_mesh", result.kpoint_mesh)
            )
            result_file.write(
                toml_line("kpoint_count", result.kpoint_count)
            )
            result_file.write(
                toml_line("runtime_seconds", result.runtime_seconds)
            )
            result_file.write(toml_line("message", result.message))

            result_file.write("\n[outputs]\n")
            for key in sorted(result.outputs):
                result_file.write(
                    toml_line(key, result.outputs[key])
                )

            result_file.write("\n[job]\n")
            result_file.write(toml_line("edge", result.job.edge))
            result_file.write(
                toml_line("job_name", result.job.job_name)
            )
            result_file.write(
                toml_line("basis_scf", result.job.basis_scf)
            )
            result_file.write(
                toml_line("basis_pscf", result.job.basis_pscf)
            )


# Register the default wingbeat at import time.
register_wingbeat("imago", ImagoWingbeat())
