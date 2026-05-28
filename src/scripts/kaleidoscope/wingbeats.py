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

from .model import WingbeatOutcome, KaleidoscopeError
from .workspace import emit_scalar, toml_line


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
#   resolve "imago" without any campaign-side setup.
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

        if self._is_prepared(wingbeat_dir):
            result = imago.run_prepared(wingbeat_dir)
        else:
            result = imago.run_structure(
                unit.structure, unit.options, wingbeat_dir
            )

        self._persist_result(wingbeat_dir, result)

        # Map the Imago-native status onto the generic outcome:
        #   "ran" covers CONVERGED / NOT_CONVERGED / SKIPPED;
        #   only a hard FAILED is not-ok.  The status value
        #   becomes the opaque detail string the campaign records.
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
    def _persist_result(wingbeat_dir, result):
        """Write the ImagoResult to ``<wingbeat_dir>/result.toml`` for
        the client's harvest.  Flat scalar fields first, then an
        ``[outputs]`` table of logical-key -> path and a ``[job]``
        table echoing what ran."""
        os.makedirs(wingbeat_dir, exist_ok=True)
        path = os.path.join(wingbeat_dir, "result.toml")
        with open(path, "w") as result_file:
            result_file.write(
                toml_line("status", result.status.value)
            )
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
