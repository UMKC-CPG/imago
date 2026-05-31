"""guidance_harvest.py -- turn a finished kaleidoscope flight into
staged historical-guidance entries (DESIGN 7.8 harvest half;
PSEUDOCODE 15.7).

Role in the pipeline
--------------------
This is the *producer* half of the library / producer / consumer
split (DESIGN 7 / ARCHITECTURE 10).  After kaleidoscope has
dispatched a predict-then-verify flight (built by
``kaleidoscope.builders.predict_verify``), each structure has been
run across a small grid of k-point densities.  This script reads
that finished workspace back off disk, finds the converged grid
point for each structure, and writes a rich ``GuidanceEntry`` to
``staging/<system_type>/`` so a curator can later promote it
(``guidance_promote.py``).  The promoted entries are what the
predictor learns from, closing the loop: every converged sweep
sharpens the next prediction (VISION Goal 5).

Where each entry field comes from (the "Model 1" sourcing, so the
information flow stays simple and homogeneous -- three inputs,
each with one clear job):

  * ``flight.toml``  -- the *plan*: the unit list (id, structure,
    calc tags), the ``[flight.sweep]`` block (which axis varied and
    which sub-model axes were held fixed), and the
    ``[flight.prediction]`` block the builder stashed.  The swept
    k-density of each grid point is read out of its calc tag (e.g.
    ``kpt-density-100``) using the sweep's ordered ``varied_axes``.
  * each run's ``result.toml`` -- the *per-run facts*: the final
    SCF total energy (for the convergence test and the
    grid-energies array), the measured ``gap_ev`` / ``gap_kind`` /
    ``total_magnetization``, and the ``scf_threshold`` the run used
    (the entry's ``metric_threshold``).  imago.py writes all of
    these (the result.toml is self-contained, DESIGN 6.1).
  * the structure ``.skl`` -- the *structural facts*: the harvest
    must load it anyway to compute the signature
    (composition + lattice family), and the same load yields
    ``cell_atom_count`` (``num_atoms``) and
    ``cell_volume_per_formula_unit`` (the cell volume, Bohr^3).

The convergence rule (DESIGN 7.8 step 3c) is two-sided: a grid
point counts as converged only when BOTH its neighbours' total
energies are within ``metric_threshold`` of it, so a single
numerical fluke cannot masquerade as a flat region.  A sweep
whose energy is still moving at the top of the grid earns no
entry -- it is logged, the flight is tagged
``prediction_mismatch``, and the structure is skipped (the user
widens the grid and re-runs; DESIGN 7.9).

Two v1 conventions, settled with the programmer 2026-05-30:
  * ``metric_threshold = scf_threshold`` -- the same criterion the
    SCF converged to is reused as the grid-flatness threshold.
  * ``imago_commit`` falls back to ``"unknown"`` when the producer
    did not inject a build identity (the C74 producer wiring will
    supply it; C78 hardens build-identity stamping).
  * ``cell_volume_per_formula_unit`` is the whole-cell volume
    (formula-unit count Z = 1) in v1; it is curator-facing metadata
    only -- the predictor never reads it -- so the formula-unit
    reduction is deferred.

This script is intentionally a thin reader over already-tested
pieces: ``guidance_db`` owns the schema, validation, and the
byte-deterministic emitter (``save_entry``); ``kaleidoscope``
owns the flight.toml format (``read_flight_toml``).
"""

from __future__ import annotations

import argparse
import os
import sys
import tomllib
from datetime import datetime, timezone

from kaleidoscope.workspace import read_flight_toml, flight_id_of
from guidance_db import (
    Context,
    GuidanceEntry,
    Measured,
    Provenance,
    Verification,
    compute_signature,
    save_entry,
    load,
)
from structure_control import StructureControl, BOHR_RAD


# Cubic Bohr per cubic Angstrom: the structure stores its lattice
#   (hence ``real_cell_volume``) in Angstroms, but the guidance
#   schema records volumes in Bohr^3, so a cell volume is divided
#   by BOHR_RAD^3 (BOHR_RAD is the Bohr radius in Angstroms).
_ANGSTROM3_TO_BOHR3 = 1.0 / (BOHR_RAD ** 3)

# The build-identity stand-in used when the producer injected no
#   commit (DESIGN 7.8 / C74).  It is non-empty so a harvested
#   flight entry still satisfies the schema's rule-11 check; the
#   curator can spot it on review.
_UNKNOWN_COMMIT = "unknown"


# ==============================================================
#  Reading the swept value back out of a calc tag (DESIGN 6.2.4)
# ==============================================================

def decode_axis_value(token: str) -> float:
    """Invert ``predict_verify.encode_axis_value``: turn the value
    portion of a calc tag back into a number.  ``"p"`` was the
    stand-in for the decimal point and a leading ``"m"`` for a
    minus sign, so ``"100" -> 100``, ``"1p5" -> 1.5``,
    ``"m2" -> -2``.  An integer-valued result is returned as an
    ``int`` so it round-trips to the same on-disk tag and compares
    cleanly as a grid value."""

    negative = token.startswith("m")
    if negative:
        token = token[1:]
    value = float(token.replace("p", "."))
    if negative:
        value = -value
    return int(value) if value == int(value) else value


def swept_value_of(unit, axis: str) -> float:
    """Read ``unit``'s value for the swept ``axis`` out of its calc
    tag.  Each calc component is ``"<axis>-<encoded-value>"``
    (DESIGN 6.2.4); the component whose prefix matches ``axis``
    carries the value.  Raises when the unit has no component for
    the swept axis -- a flight whose units do not match its own
    SweepRecord is malformed and must fail loudly, not be guessed
    around."""

    prefix = axis + "-"
    for component in unit.calc:
        if component.startswith(prefix):
            return decode_axis_value(component[len(prefix):])
    raise ValueError(
        "unit id=" + repr(unit.id) + " calc=" + repr(unit.calc)
        + " has no component for swept axis " + repr(axis))


# ==============================================================
#  The two-sided convergence rule (DESIGN 7.8 step 3c)
# ==============================================================

def pick_converged(energies, threshold):
    """Return the smallest interior grid index whose total energy
    is within ``threshold`` of BOTH neighbours, or ``None`` when
    the energy is still moving (no flat interior point).

    Two-sided by design: requiring both consecutive-pair deltas to
    be small means a single-grid-point numerical dip cannot be
    mistaken for convergence (DESIGN 7.8 step 3c, stricter than a
    one-sided "delta below threshold" rule).  Endpoints are never
    eligible -- a converged-at-the-edge sweep is suspect (the grid
    may have been too narrow) and is left for the auto-promote
    rule to reject (DESIGN 7.8)."""

    for index in range(1, len(energies) - 1):
        below_up = abs(energies[index] - energies[index + 1]) < threshold
        below_down = abs(energies[index] - energies[index - 1]) < threshold
        if below_up and below_down:
            return index
    return None


# ==============================================================
#  Side-channel marker for a non-converged sweep (DESIGN 7.8 3d)
# ==============================================================

def tag_prediction_mismatch(workspace_root: str, unit_id: str) -> None:
    """Record that ``unit_id``'s sweep did not converge within its
    grid (DESIGN 7.8 step 3d / 7.9).  Written as a small
    ``prediction_mismatch.toml`` at the workspace root, keyed by
    unit id, so a later review (or a re-run script) can see which
    structures need a widened grid without re-deriving it.  The
    write merges into any existing marker file so several skipped
    structures accumulate."""

    path = os.path.join(workspace_root, "prediction_mismatch.toml")
    flagged = {}
    if os.path.exists(path):
        with open(path, "rb") as marker_file:
            flagged = tomllib.load(marker_file)
    flagged[unit_id] = True
    with open(path, "w") as marker_file:
        for key in sorted(flagged):
            # unit ids are slugs ([a-z0-9_-]+), so they are valid
            #   bare TOML keys and need no quoting.
            marker_file.write(key + " = true\n")


# ==============================================================
#  Building a structure's GuidanceEntry from its sweep
# ==============================================================

def _load_structure(path: str):
    """Load the structure ``.skl`` into a StructureControl.  Uses
    ``read_input_file`` (not the bare skeleton reader) so the
    element mapping AND the cell geometry -- hence
    ``real_cell_volume`` -- are computed, which the entry's
    ``cell_volume_per_formula_unit`` needs."""

    structure = StructureControl()
    structure.read_input_file(path)
    return structure


def _now_iso8601_utc() -> str:
    """The harvest timestamp in the schema's ISO-8601 UTC form
    (``2026-05-30T14:03:55Z``)."""

    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_entry(workspace_root, grid, kpds, energies, rts, idx,
                prediction, fixed_axes, dataspace):
    """Assemble the rich :class:`GuidanceEntry` for one converged
    structure sweep (DESIGN 7.8 step 3f; PSEUDOCODE 15.7 step g).

    ``grid``/``kpds``/``energies``/``rts`` are the per-grid-point
    units, swept k-densities, total energies, and parsed
    result.toml dicts (all parallel, sorted by k-density); ``idx``
    is the converged index chosen by :func:`pick_converged`.
    ``prediction`` is the ``[flight.prediction]`` dict (or None for
    a hand-built flight) and ``fixed_axes`` the sweep's held-
    constant sub-model axes.  The ``entry_id`` is left empty here;
    :func:`save_entry` fills it with the deterministic slug."""

    chosen_rt = rts[idx]
    chosen_kpd = kpds[idx]
    threshold = rts[0].get("scf_threshold")

    # system_type rides on the prediction record (it carries it
    #   from the builder, DESIGN 7.7); a hand-built flight with no
    #   prediction defaults to crystalline.
    system_type = (prediction["system_type"]
                   if prediction is not None else "crystalline")

    structure = _load_structure(grid[0].structure)
    signature = compute_signature(
        structure, system_type, dataspace.group_table)

    measured = Measured(
        gap_ev=chosen_rt["gap_ev"],
        gap_kind=chosen_rt["gap_kind"],
        # spin_polarization is not surfaced by imago (the iteration
        #   file carries the magnetic moment, not a polarization),
        #   so it is recorded as 0.0; the predictor's spin character
        #   keys on total_magnetization instead (the C72 decision).
        spin_polarization=0.0,
        total_magnetization=chosen_rt.get("total_magnetization", 0.0),
        kpoint_density=chosen_kpd)

    context = Context(
        basis=fixed_axes["basis"],
        functional=fixed_axes["functional"],
        kpoint_integration=fixed_axes["kpoint_integration"],
        scf_threshold=threshold,
        cell_atom_count=structure.num_atoms,
        cell_volume_per_formula_unit=(
            structure.real_cell_volume * _ANGSTROM3_TO_BOHR3))

    verification = Verification(
        grid_values=tuple(kpds),
        grid_energies=tuple(energies),
        converged_at=chosen_kpd,
        metric="total_energy",
        metric_threshold=threshold,
        predictor_confidence=(
            prediction["confidence"]
            if prediction is not None else 0.0),
        predictor_neighbor_ids=tuple(
            prediction["neighbor_entry_ids"]
            if prediction is not None else ()))

    commit = chosen_rt.get("imago_commit") or _UNKNOWN_COMMIT
    provenance = Provenance(
        flight_id=flight_id_of(workspace_root),
        source_structure=grid[0].structure,
        imago_commit=commit,
        curator="guidance_harvest.py")

    return GuidanceEntry(
        entry_id="",
        generated_at=_now_iso8601_utc(),
        source="flight",
        signature=signature,
        measured=measured,
        context=context,
        verification=verification,
        provenance=provenance)


# ==============================================================
#  The harvest driver (DESIGN 7.8; PSEUDOCODE 15.7)
# ==============================================================

def _read_result_toml(workspace_root, unit):
    """Parse one unit's ``result.toml`` from its run directory
    (``<root>/wingbeats/<id>/<calc...>/result.toml``, the layout
    ``kaleidoscope.workspace.unit_run_dir`` writes)."""

    path = os.path.join(workspace_root, "wingbeats", unit.id,
                        *unit.calc, "result.toml")
    with open(path, "rb") as result_file:
        return tomllib.load(result_file)


def harvest_flight(workspace_root, db_root, dataspace):
    """Walk a finished flight workspace and stage a GuidanceEntry
    for every converged structure sweep (DESIGN 7.8; PSEUDOCODE
    15.7).  Returns a list of one human-readable summary line per
    structure (staged / skipped / trusted) so the CLI -- or a
    caller such as the C48.3 producer -- can report what happened.

    ``dataspace`` supplies only ``group_table`` (for
    ``compute_signature``); the staged files are written under
    ``db_root`` via ``save_entry``.  No real Imago run is needed to
    exercise this: it reads flight.toml, the per-unit result.toml
    files, and each structure .skl."""

    flight = read_flight_toml(
        os.path.join(workspace_root, "flight.toml"))
    prediction = flight.metadata.get("prediction")

    # The swept axis (v1: a single axis, "kpt-density") names which
    #   calc-tag component carries each grid point's value, and the
    #   fixed axes carry the sub-model context (basis/functional/
    #   kpoint_integration).  A hand-built flight with no sweep
    #   cannot be harvested into the k-density dataspace.
    if flight.sweep is None or not flight.sweep.varied_axes:
        raise ValueError(
            workspace_root + ": flight.toml has no [flight.sweep] "
            "with a varied axis -- nothing to harvest")
    axis = flight.sweep.varied_axes[0]
    fixed_axes = flight.sweep.fixed_axes

    # Group the units by structure id: one verification sub-grid
    #   per structure (insertion order preserved for a stable
    #   summary).
    groups: dict[str, list] = {}
    for unit in flight.units:
        groups.setdefault(unit.id, []).append(unit)

    summaries = []
    for unit_id, units in groups.items():
        # a. Sort the sub-grid by swept k-density.
        grid = sorted(units, key=lambda u: swept_value_of(u, axis))

        # b. Parse each run's result.toml for the energy and the
        #    measured quantities.
        kpds, energies, rts = [], [], []
        for unit in grid:
            rt = _read_result_toml(workspace_root, unit)
            kpds.append(swept_value_of(unit, axis))
            energies.append(rt["total_energy"])
            rts.append(rt)

        # c. Trust mode harvests deliverables but stages NO entry
        #    (DESIGN 6.2.1 / 7.7): one converged calc is weaker
        #    evidence than a grid.  This MUST precede pick_converged
        #    -- a trust grid is a single point, so the two-sided
        #    test below would report it as "still moving".
        if prediction is not None \
                and prediction.get("policy") == "trust_no_verify":
            summaries.append(
                unit_id + ": trusted point (not staged)")
            continue

        # d. metric_threshold = the SCF criterion the runs used
        #    (the v1 convention); pick the converged grid point.
        threshold = rts[0].get("scf_threshold")
        if threshold is None:
            raise ValueError(
                unit_id + ": result.toml carries no scf_threshold "
                "(needed as the convergence metric_threshold)")
        idx = pick_converged(energies, threshold)

        # e. Energy still moving at the top of the grid: tag the
        #    flight and skip -- a non-converged sweep earns no entry.
        if idx is None:
            tag_prediction_mismatch(workspace_root, unit_id)
            summaries.append(
                unit_id + ": energy still moving at top of grid "
                "-- skipped (tagged prediction_mismatch)")
            continue

        # f/g. Build the rich entry and stage it.
        entry = build_entry(
            workspace_root, grid, kpds, energies, rts, idx,
            prediction, fixed_axes, dataspace)
        path = save_entry(entry, db_root)
        summaries.append(unit_id + ": staged " + path)

    return summaries


# ==============================================================
#  Command-line interface
# ==============================================================

def _default_db_root() -> str:
    """The dataspace root under $IMAGO_DATA (DESIGN 7 layout:
    ``share/historicalGuidanceDB/``).  Returned empty when
    $IMAGO_DATA is unset so the parser can demand ``--db-root``."""

    data_dir = os.environ.get("IMAGO_DATA", "")
    return (os.path.join(data_dir, "historicalGuidanceDB")
            if data_dir else "")


def main(argv=None):
    """CLI entry point: harvest one finished flight workspace into
    the dataspace's staging area (DESIGN 7.8)."""

    parser = argparse.ArgumentParser(
        description="Harvest a finished kaleidoscope flight into "
                    "staged historical-guidance entries "
                    "(DESIGN 7.8).")
    parser.add_argument(
        "workspace_root",
        help="the flight workspace directory (holding flight.toml "
             "and the wingbeats/ run tree)")
    parser.add_argument(
        "--db-root", default=_default_db_root(),
        help="the historicalGuidanceDB root to stage into "
             "(default: $IMAGO_DATA/historicalGuidanceDB)")
    args = parser.parse_args(argv)

    if not args.db_root:
        parser.error(
            "--db-root not given and $IMAGO_DATA is unset")

    # The dataspace is loaded for its group_table (compute_signature
    #   needs it) and is the same root the staged entries land under.
    dataspace = load(args.db_root)
    summaries = harvest_flight(
        args.workspace_root, args.db_root, dataspace)
    for line in summaries:
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
