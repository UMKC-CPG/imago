#!/usr/bin/env python3
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

"""guidance_harvest.py -- turn a finished kaleidoscope flight into
staged historical-guidance entries (DESIGN 7.8 harvest half;
PSEUDOCODE 15.7).

Role in the pipeline
--------------------
This is the *producer* half of the library / producer / consumer
split (DESIGN 7 / ARCHITECTURE 10).  After kaleidoscope has
dispatched a predict-then-verify flight (built by
``kaleidoscope.builders.kpoint_convergence``), each structure has
been
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
    calc tags), the ``[flight.sweep]`` block (which axis varied; in
    v1 nothing is held fixed), and the per-structure
    ``[flight.predictions.<id>]`` blocks the builder stashed.  Each
    block carries that structure's prediction AND the (basis,
    functional, kpoint_integration) sub-model it ran under -- the
    sole home for the sub-model, never duplicated into
    ``fixed_axes`` (DESIGN 6.2.9).  The swept k-density of each grid
    point is read out of its calc tag (e.g. ``kpt-density-100``)
    using the sweep's ordered ``varied_axes``.
  * each run's ``result.toml`` -- the *per-run facts*: the final
    SCF total energy (for the convergence test and the
    grid-energies array), the measured ``gap_ev`` / ``gap_kind`` /
    ``total_magnetization``, and the ``scf_threshold`` the run used
    (recorded in the entry's context, distinct from the k-point
    ``metric_threshold``).  imago.py writes all of these (the
    result.toml is self-contained, DESIGN 6.1).
  * the structure ``.skl`` -- the *structural facts*: the harvest
    must load it anyway to compute the signature
    (composition + lattice family), and the same load yields
    ``cell_atom_count`` (``num_atoms``) and
    ``cell_volume_per_formula_unit`` (the cell volume, Bohr^3).

The convergence rule (DESIGN 7.8 step 3c) is two-sided: a grid
point counts as converged only when BOTH its neighbours' total
energies are within ``metric_threshold`` of it -- PER ATOM and in
eV, the basis the threshold is stated in -- so a single numerical
fluke cannot masquerade as a flat region.  A sweep
whose energy is still moving at the top of the grid earns no
entry -- it is logged, the flight is tagged
``prediction_mismatch``, and the structure is skipped (the user
widens the grid and re-runs; DESIGN 7.9).

v1 conventions, settled with the programmer:
  * The grid-flatness ``metric_threshold`` is the solid's resolved
    ``kpoint_convergence_threshold`` (per atom, in eV; DESIGN 7.8 /
    5.7), which rides on the structure's prediction record -- a
    manifest/resolved fact, absent from result.toml -- and is
    DISTINCT from the run's own ``scf_threshold``.  ``grid_energies``
    are stored RAW (total-cell hartree, Option B); every comparison
    normalizes per atom at the point of use.
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
from structure_control import StructureControl, BOHR_RAD, HARTREE


# Cubic Bohr per cubic Angstrom: the structure stores its lattice
#   (hence ``real_cell_volume``) in Angstroms, but the guidance
#   schema records volumes in Bohr^3, so a cell volume is divided
#   by BOHR_RAD^3 (BOHR_RAD is the Bohr radius in Angstroms).
_ANGSTROM3_TO_BOHR3 = 1.0 / (BOHR_RAD ** 3)


def per_atom_ev(total_energy_hartree: float,
                cell_atom_count: int) -> float:
    """A raw total-cell energy (hartree) expressed as eV per atom
    -- the basis the k-point flatness threshold is stated in
    (DESIGN 7.8 / 5.7, Option B).  Single-sourced here so the
    convergence pick and the auto-promote flatness test normalize
    identically and cannot drift.  ``HARTREE`` is the hartree->eV
    factor (structure_control)."""

    return total_energy_hartree * HARTREE / cell_atom_count

# The build-identity stand-in used when the producer injected no
#   commit (DESIGN 7.8 / C74).  It is non-empty so a harvested
#   flight entry still satisfies the schema's rule-11 check; the
#   curator can spot it on review.
_UNKNOWN_COMMIT = "unknown"


# ==============================================================
#  Reading the swept value back out of a calc tag (DESIGN 6.2.4)
# ==============================================================

def decode_axis_value(token: str) -> float:
    """Invert ``kpoint_convergence.encode_axis_value``: turn the value
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

def pick_converged_climb(energies, cell_atom_count, threshold,
                         flat_needed):
    """Return the smallest interior index whose per-atom energy is
    flat over ``flat_needed`` consecutive interior rungs, or ``None``
    when no such run exists yet (PSEUDOCODE 4e.2; DESIGN 3.12.3).

    A single rung is "flat" when its per-atom energy is within
    ``threshold`` of BOTH neighbours -- the same two-sided test
    :func:`pick_converged` applies (DESIGN 7.8 step 3c).  The
    adaptive climb generalises it by demanding the flatness PERSIST:
    the returned index and the next ``flat_needed - 1`` interior
    rungs must all pass that two-sided test.  A confident search
    sets ``flat_needed = 1`` (a single flat interior rung is
    enough); a cold or bootstrap search sets ``flat_needed = 2`` so
    one lucky flat step cannot end the climb prematurely (DESIGN
    3.12.3).  With ``flat_needed = 1`` this is exactly
    :func:`pick_converged`.

    ``energies`` are raw total-cell values in hartree (Option B),
    normalized once to eV per atom (:func:`per_atom_ev`) before
    comparing, and ``threshold`` is stated per atom -- so a large
    cell is not held to a tighter bound than a small one (DESIGN
    7.8).  Endpoints are never eligible: a rung needs both a lower
    and an upper neighbour to be judged, and a sweep that only goes
    flat at its top edge is suspect and left for the auto-promote
    rule (DESIGN 7.8)."""

    per_atom = [per_atom_ev(energy, cell_atom_count)
                for energy in energies]

    def two_sided_flat(rung):
        below_down = abs(
            per_atom[rung] - per_atom[rung - 1]) < threshold
        below_up = abs(
            per_atom[rung] - per_atom[rung + 1]) < threshold
        return below_down and below_up

    for first in range(1, len(per_atom) - 1):
        last = first + flat_needed - 1     # last interior to confirm
        if last > len(per_atom) - 2:       # not enough rungs above
            break                          #   to confirm `first` yet
        if all(two_sided_flat(rung)
               for rung in range(first, last + 1)):
            return first
    return None


def pick_converged(energies, cell_atom_count, threshold):
    """Return the smallest interior grid index whose total energy
    is within ``threshold`` of BOTH neighbours, or ``None`` when
    the energy is still moving (no flat interior point).

    ``energies`` are raw total-cell values in hartree (Option B) and
    ``threshold`` is per atom, in eV, so the ladder is normalized to
    that basis once (:func:`per_atom_ev`) before comparing.  The
    per-atom scale keeps a large cell from being held to a tighter
    bound than a small one (DESIGN 7.8).

    Two-sided by design: requiring both consecutive-pair deltas to
    be small means a single-grid-point numerical dip cannot be
    mistaken for convergence (DESIGN 7.8 step 3c, stricter than a
    one-sided "delta below threshold" rule).  Endpoints are never
    eligible -- a converged-at-the-edge sweep is suspect (the grid
    may have been too narrow) and is left for the auto-promote
    rule to reject (DESIGN 7.8).

    This is the ``flat_needed = 1`` case of
    :func:`pick_converged_climb`, kept under its own name as the
    single-grid convergence pick the harvest and auto-promote call
    -- the two share the one two-sided rule and cannot drift."""

    return pick_converged_climb(energies, cell_atom_count,
                                threshold, 1)


def stride_is_flat(lo_rung, hi_rung, cell_atom_count, threshold):
    """Return whether a bracket stride's two endpoints are flat
    (PSEUDOCODE 4e.2; DESIGN 3.12.3).

    The bracket phase of the mesh climb strides across many ladder
    positions at once and computes an energy only at each stride's
    endpoints.  A stride is "flat" when those two endpoints' per-atom
    energies are within the SAME ``threshold`` the two-sided
    convergence pick uses.  Because a stride adds many k-points, a
    small energy change across it is strong evidence the energy has
    settled -- but only evidence: the refine phase VERIFIES a
    proposed bracket with the full two-sided test
    (:func:`pick_converged_climb`), so a coincidentally flat stride
    (an oscillating near-metal energy that dips and returns) is
    caught there, not trusted here.

    Lives here, beside :func:`pick_converged_climb` and
    :func:`per_atom_ev`, so every energy-flatness test in the climb
    normalizes on the one per-atom rule and none can drift -- the
    same reason the two-sided pick is single-sourced here rather than
    in the pure-geometry ``mesh_climb``.

    ``lo_rung`` and ``hi_rung`` are the stride's lower and upper
    endpoints, each an object exposing a raw total-cell ``.energy``
    in hartree (Option B); ``threshold`` is per atom, in eV, so a
    large cell is not held to a tighter bound than a small one
    (DESIGN 7.8)."""

    lo_per_atom = per_atom_ev(lo_rung.energy, cell_atom_count)
    hi_per_atom = per_atom_ev(hi_rung.energy, cell_atom_count)
    return abs(hi_per_atom - lo_per_atom) < threshold


def is_gapless(rung, gap_threshold):
    """Return whether a rung is metallic -- its computed band gap at
    or below ``gap_threshold`` -- the metal test (PSEUDOCODE 4e.2;
    DESIGN 3.12.3).

    A metal has no band gap: its total energy oscillates as the
    k-point mesh crosses the Fermi surface and never settles, so
    chasing k-point convergence on it is futile.  The gap is read
    straight from the rung's own result -- a DIRECT metal signal,
    unlike the retired rising-stride proxy that inferred metallicity
    from a finer mesh raising the energy and so missed the common
    small-amplitude oscillator whose rise never cleared the margin.

    ``gap_threshold`` is an absolute band gap in eV (not a per-atom
    energy): low enough that no real insulator crosses it, high enough
    to catch a true metal's near-zero reading (DESIGN 3.12.6).
    ``rung`` exposes a ``.gap`` in eV, taken from its result's
    ``gap_ev``; a rung whose gap is unknown (``None``) is treated as
    NON-metallic, so a missing reading never spuriously stops a
    climb.

    A thin wrapper over :func:`is_gapless_value`: this side knows how
    to find a gap on a *rung*, that side holds what the gap means.
    The guidance harvest calls the same core on a parsed
    ``result.toml`` (:func:`build_entry`), so the climb's metal
    short-circuit and the harvest's metal skip cannot drift apart
    (DESIGN 7.8)."""
    return is_gapless_value(rung.gap, gap_threshold)


def is_gapless_value(gap_ev, gap_threshold):
    """Return whether a bare band-gap reading is metallic -- at or
    below ``gap_threshold`` (DESIGN 3.12.3 / 7.8).

    ``gap_threshold`` is an ABSOLUTE band gap in eV, not a per-atom
    energy: low enough that no real insulator crosses it, high enough
    to catch a true metal's near-zero reading.

    A MISSING gap (``None``) is NOT metallic.  Both callers depend on
    that, for the same reason from opposite ends: in the climb an
    absent reading must not stop a search that was converging, and in
    the harvest it must not suppress a guidance entry a genuine
    insulator earned.  Defaulting the other way would make an unwired
    gap look like a collection with no insulators in it -- the sort of
    failure nobody thinks to question.

    The scalar core, so the rung-shaped :func:`is_gapless` and the
    result-dict-shaped call in :func:`build_entry` share ONE rule and
    neither caller has to build a shape it does not have."""
    return gap_ev is not None and gap_ev <= gap_threshold


# Two runs of the same resolved mesh are the same calculation and
#   must give the same total energy; ENERGY_MATCH_EPS is the gap
#   below which they count as equal.  It is tight (near float
#   noise, in hartree): a single-threaded imago run is
#   deterministic, so genuine duplicates agree to the last digits
#   and a wider gap means the runs were not in fact identical.
ENERGY_MATCH_EPS = 1e-9


def collapse_by_mesh(kpoint_densities, energies, meshes):
    """Reduce a density-sorted grid to one rung per distinct
    resolved mesh, keeping the lowest-density member (DESIGN 7.8
    step 3c guard; PSEUDOCODE 15.7).

    The three inputs are parallel arrays ordered by ascending
    requested density; ``meshes[i]`` is grid point i's resolved
    ``kpoint_mesh`` (the axial counts, DESIGN 6.1.2), or ``None``.
    Returns ``(collapsed_densities, collapsed_energies, kept)``,
    where ``kept[j]`` is the ORIGINAL index of the j-th surviving
    rung, so a caller can map a collapsed index back to the run it
    names.

    Two rungs that resolved to the same mesh are one calculation
    run twice; their energy delta is exactly zero, which the
    two-sided flatness test would misread as convergence
    (DESIGN 3.11).  Collapsing removes that manufactured zero
    before the test runs.

    If any mesh is ``None`` -- an older result.toml, or an imago
    binary that does not yet emit the mesh -- the guard cannot act
    and returns the grid unchanged with identity indices, so it
    stays inert until result.toml carries the mesh (DESIGN 6.1.2 /
    3.11)."""

    if any(mesh is None for mesh in meshes):
        return (kpoint_densities, energies,
                list(range(len(meshes))))

    kept = []
    for index in range(len(meshes)):
        # Same mesh as the last surviving rung?  Then this rung is
        #   that calculation run again (the density-to-mesh map is
        #   monotone in density, DESIGN 3.7, so equal meshes are
        #   contiguous).  An equal mesh MUST give an equal energy;
        #   a mismatch means the runs were not identical, and is
        #   surfaced rather than averaged away.
        if kept and meshes[index] == meshes[kept[-1]]:
            if abs(energies[index]
                   - energies[kept[-1]]) > ENERGY_MATCH_EPS:
                raise ValueError(
                    "k-density rungs {0} and {1} resolved to the "
                    "same mesh {2} but disagree in total energy "
                    "-- not the same calculation".format(
                        kpoint_densities[kept[-1]],
                        kpoint_densities[index], meshes[index]))
            continue
        kept.append(index)

    return ([kpoint_densities[i] for i in kept],
            [energies[i] for i in kept],
            kept)


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

def load_structure(path: str):
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


def _require_field(result_toml: dict, field: str, unit_id: str):
    """Return ``result_toml[field]`` or raise a clear error naming
    the missing field and the structure (DESIGN 7.8 / group-A #8).

    Used for the electronic-character fields ``gap_ev`` / ``gap_kind``
    that the predictor keys on: unlike ``total_magnetization`` (which
    a closed-shell run legitimately omits, defaulting to 0.0), an
    absent gap means the run did not surface the quantity the entry
    exists to record, so the harvest fails loudly rather than
    fabricating a value."""

    if field not in result_toml:
        raise ValueError(
            unit_id + ": converged run's result.toml is missing "
            "the required field " + repr(field) + " (the harvest "
            "needs the measured electronic character; re-run with "
            "an imago build that surfaces it -- TODO C76)")
    return result_toml[field]


def build_entry(workspace_root, source_structure, prediction,
                dataspace, structure, kpoint_threshold,
                grid_values, grid_energies, converged_density,
                chosen_result, ladder_is_metal=False):
    """Assemble the rich :class:`GuidanceEntry` for one converged
    structure from its ALREADY-CHOSEN facts (DESIGN 7.8 step 3f;
    PSEUDOCODE 15.7).

    This is the ONE entry builder both guidance harvests feed, so a
    schema change touches a single place and the two paths cannot
    drift (the Q1-Q2 shared core, DESIGN 5.7).  Each path picks its
    converged rung its own way -- the standalone density sweep
    (:func:`harvest_flight`) via ``collapse_by_mesh`` +
    ``pick_converged``, the producer's in-memory climb
    (``build_initial_potentials``) via the climb +
    ``record_converged`` -- then hands the identical chosen facts
    here:

      * ``grid_values`` / ``grid_energies`` -- the distinct-mesh
        flatness ladder (ascending), the second as raw total-cell
        hartree (Option B; the consumer normalizes per atom).  Stored
        as the entry's verification grid so the curator's
        ``auto_promote_ok`` re-judges flatness on the same genuinely
        distinct calculations, never re-encountering a duplicate-mesh
        zero (DESIGN 7.8 step 3f).
      * ``converged_density`` -- the chosen rung's k-point density.
      * ``chosen_result`` -- the chosen run's parsed ``result.toml``.

    Everything MEASURED (gap, magnetization), the SCF threshold, the
    exact converged mesh, and the Imago commit come from
    ``chosen_result``; the sub-model and ``system_type`` from
    ``prediction`` (its SOLE source, DESIGN 6.2.9 / 7.8 step 3f -- a
    structure with no record is skipped before reaching here, so no
    None-guards are needed); the cell facts and signature from the
    loaded ``structure`` (passed in, not reloaded here).

    ``source_structure`` is the structure's local path -- recorded in
    the entry's provenance and used to name any loud failure.
    ``kpoint_threshold`` is the resolved per-atom eV flatness
    tolerance, stored as the entry's ``metric_threshold`` (DESIGN
    7.8) -- distinct from the run's ``scf_threshold``.  The
    ``entry_id`` is left empty here; :func:`save_entry` fills it with
    the deterministic slug.

    Returns None for a METAL, which stages no guidance entry at all
    (DESIGN 7.8).  An entry's whole content is the claim "for a
    structure like this, this k-density is converged," and a metal
    cannot make it: its energy does not converge in k-points at any
    mesh worth paying for, and the climb acknowledges this by
    short-circuiting at the first gapless rung and settling there as
    a deliberately rough potential (DESIGN 3.12.3).  That settled
    rung is a stopping point, not a converged density.  Recording it
    as one would feed the predictor a claim nobody made -- and
    disproportionately, since a metal is often the only member of its
    lattice family in a young collection and would then dominate
    every prediction for that family through the distance weighting
    of DESIGN 7.6.  The guard sits HERE, in the one builder both
    harvest paths share, so neither can grow its own version of the
    rule.  The producer's potential harvest is unaffected: a rough
    starting potential is exactly what that database is for.

    TWO readings decide it and EITHER is sufficient (DESIGN 7.8 d').
    ``ladder_is_metal`` is the CALLER's multi-rung reading, passed in
    rather than worked out here because only the caller holds the
    ladder: the producer passes the climb's verdict, the standalone
    sweep passes the any-rung test over its own grid.  It defaults to
    False, which means "no multi-rung evidence offered" and leaves
    the chosen rung as the only witness -- never "known not to be a
    metal".  The chosen rung's own gap is the second reading.

    Neither reading is redundant.  A metal on a discrete mesh shows
    an artificial gap whose size depends on where the mesh points
    fall (DESIGN 1.6), so one rung is weak evidence -- fcc Al reads
    zero at several meshes and 0.124 eV at another -- while the
    ladder taken whole is strong.  Accepting either means the
    stronger reading cannot be overruled by the weaker one, and
    nothing the chosen-rung test already caught stops being caught."""

    # The caller's ladder reading first: it rests on more evidence
    #   than anything available here, and settles the question when
    #   it fires.
    if ladder_is_metal:
        return None

    # The chosen rung's own gap.  The cut is read off the prediction
    #   record, which is how a manifest knob reaches a standalone tool
    #   that never sees a manifest (the same channel kpoint_threshold
    #   uses, DESIGN 7.8 step 3d').  is_gapless_value is the scalar
    #   core the climb's rung-shaped is_gapless also calls, so one
    #   rule serves both -- including its side on missing data: an
    #   UNKNOWN gap is NOT metallic, because a missing reading must
    #   never silently suppress an entry a real insulator earned.
    metal_gap_threshold = prediction.get("metal_gap_threshold")
    if metal_gap_threshold is None:
        raise ValueError(
            source_structure + ": prediction record carries no "
            "metal_gap_threshold (the absolute eV band gap below "
            "which a run counts as metallic, which the producer "
            "resolves from the manifest and stamps on the record, "
            "DESIGN 7.8)")
    if is_gapless_value(chosen_result.get("gap_ev"),
                        metal_gap_threshold):
        return None

    # The SCF threshold is a per-run fact from the chosen run's
    #   result.toml, recorded in the entry's context; it is SEPARATE
    #   from the k-point flatness metric_threshold (kpoint_threshold,
    #   per atom in eV).  A converged run must record it (a required
    #   context fact), so an absent one fails loudly, not as None.
    scf_threshold = chosen_result.get("scf_threshold")
    if scf_threshold is None:
        raise ValueError(
            source_structure + ": converged run's result.toml "
            "carries no scf_threshold (a required context fact, "
            "DESIGN 5.2)")

    # system_type rides on this structure's prediction record (it
    #   carries it from the predictor, DESIGN 7.7).
    system_type = prediction["system_type"]

    signature = compute_signature(
        structure, system_type, dataspace.group_table)

    measured = Measured(
        # gap_ev / gap_kind are REQUIRED measured quantities (the
        #   electronic character the predictor keys on); a run that
        #   did not surface them cannot earn an entry, so this is a
        #   loud failure rather than a silent default.
        gap_ev=_require_field(
            chosen_result, "gap_ev", source_structure),
        gap_kind=_require_field(
            chosen_result, "gap_kind", source_structure),
        # spin_polarization is not surfaced by imago (the iteration
        #   file carries the magnetic moment, not a polarization),
        #   so it is recorded as 0.0; the predictor's spin character
        #   keys on total_magnetization instead (the C72 decision).
        spin_polarization=0.0,
        total_magnetization=chosen_result.get(
            "total_magnetization", 0.0),
        kpoint_density=converged_density)

    # The sub-model is read from THIS structure's record -- never
    #   from sweep.fixed_axes (empty), so a combined mixed-sub-model
    #   flight harvests each structure correctly (DESIGN 6.2.9 / 7.8
    #   step 3f).
    context = Context(
        basis=prediction["basis"],
        functional=prediction["functional"],
        kpoint_integration=prediction["kpoint_integration"],
        scf_threshold=scf_threshold,
        cell_atom_count=structure.num_atoms,
        cell_volume_per_formula_unit=(
            structure.real_cell_volume * _ANGSTROM3_TO_BOHR3))

    verification = Verification(
        # The distinct-mesh flatness ladder the picker judged, not a
        #   raw duplicate-bearing ladder, so the stored evidence is
        #   free of the duplicate-mesh zero (DESIGN 7.8 step 3f).
        grid_values=tuple(grid_values),
        # grid_energies are RAW total-cell hartree (Option B);
        #   consumers (auto_promote_ok) normalize per atom.
        grid_energies=tuple(grid_energies),
        converged_at=converged_density,
        # The chosen rung's resolved mesh, stored exact beside the
        #   density so the calculation is auditable where a density
        #   round-trips only up to rounding (DESIGN 3.12.4 / 7.2).
        #   Read from the chosen run's result.toml (the resolved axial
        #   counts, DESIGN 6.1.2) -- the SAME source in both harvests,
        #   so the mesh cannot differ between them; absent on an older
        #   run -> None.
        converged_mesh=(tuple(chosen_result["kpoint_mesh"])
                        if chosen_result.get("kpoint_mesh")
                        is not None else None),
        metric="total_energy",
        metric_threshold=kpoint_threshold,
        predictor_confidence=prediction["confidence"],
        predictor_neighbor_ids=tuple(
            prediction["neighbor_entry_ids"]))

    # The build behind the run, read from result.toml like every other
    #   per-run fact: the wingbeat echoed it there out of the unit's
    #   `record` (DESIGN 6.2.2/6.2.4), so this harvest stays on its
    #   three sources and never opens the dispatch core's status.toml.
    #   _UNKNOWN_COMMIT remains the floor for a run that recorded
    #   nothing -- non-empty, so the schema's rule-11 check passes and
    #   a curator can spot it on review.
    commit = chosen_result.get("imago_commit") or _UNKNOWN_COMMIT
    provenance = Provenance(
        flight_id=flight_id_of(workspace_root),
        source_structure=source_structure,
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

    # Per-structure predictions, keyed by structure id (DESIGN
    #   6.2.9); a single-structure flight carries a one-entry map.
    #   Each record is the SOLE source of its structure's
    #   system_type and (basis, functional, kpoint_integration)
    #   sub-model (7.8 step 3f), so a structure with no record is
    #   skipped below.
    predictions = flight.metadata.get("predictions", {})

    # The swept axis (v1: a single axis, "kpt-density") names which
    #   calc-tag component carries each grid point's value.  A
    #   flight with no sweep cannot be harvested into the k-density
    #   dataspace.  The sub-model is NOT read from sweep.fixed_axes
    #   (now empty): it rides on the per-structure record (3f).
    if flight.sweep is None or not flight.sweep.varied_axes:
        raise ValueError(
            workspace_root + ": flight.toml has no [flight.sweep] "
            "with a varied axis -- nothing to harvest")
    axis = flight.sweep.varied_axes[0]

    # Keep only the convergence-sweep runs (DESIGN 6.2.9 / 7.8 step
    #   2): a producer flight also carries structure-only
    #   "fingerprint" loen units that share a structure id but
    #   belong to a different harvester, and must not be mistaken
    #   for grid points.  Then group by structure id: one
    #   verification sub-grid per structure (insertion order
    #   preserved for a stable summary).
    groups: dict[str, list] = {}
    for unit in flight.units:
        if unit.kind != "convergence":
            continue
        groups.setdefault(unit.id, []).append(unit)

    summaries = []
    for unit_id, units in groups.items():
        # The prediction this structure was launched under.  No
        #   record -> not guidance-harvestable (the record is the
        #   sole source of system_type + sub-model, 7.8 step 3f);
        #   the builder always attaches one, so this only fires for
        #   a hand-built flight outside the predict-then-verify
        #   path (which seeds guidance by hand instead, 7.9).
        prediction = predictions.get(unit_id)
        if prediction is None:
            summaries.append(
                unit_id + ": no prediction record (not staged)")
            continue

        # a. Sort the sub-grid by swept k-density.
        grid = sorted(units, key=lambda u: swept_value_of(u, axis))

        # b. Parse each run's result.toml for the energy, the
        #    measured quantities, and the resolved mesh.
        #    ``meshes[i]`` feeds the duplicate-mesh guard in step
        #    (d); it is None when result.toml carries no
        #    kpoint_mesh, which the guard treats as "cannot
        #    collapse" (collapse_by_mesh), keeping it inert.
        kpoint_densities, energies, meshes, result_tomls = [], [], [], []
        for unit in grid:
            result_toml = _read_result_toml(workspace_root, unit)
            kpoint_densities.append(swept_value_of(unit, axis))
            energies.append(result_toml["total_energy"])
            meshes.append(result_toml.get("kpoint_mesh"))
            result_tomls.append(result_toml)

        # c. A single-point grid harvests deliverables but stages NO
        #    entry (DESIGN 6.2.1 / 7.7): one converged calc is weaker
        #    evidence than a grid.  Covers both trust mode and a
        #    single-point curator override, and MUST precede
        #    pick_converged -- the two-sided test below needs >= 3
        #    points and would report one point as "still moving".
        if len(grid) == 1:
            summaries.append(
                unit_id + ": single point (not staged)")
            continue

        # d. The k-point flatness tolerance rode in on this
        #    structure's prediction record: per atom, in eV, the
        #    solid's resolved kpoint_convergence_threshold (DESIGN
        #    7.8 / 5.7).  It is a manifest/resolved fact, absent from
        #    result.toml.  Load the structure once (for
        #    cell_atom_count here, reused by build_entry), then pick
        #    the converged grid point.
        kpoint_threshold = prediction.get(
            "kpoint_convergence_threshold")
        if kpoint_threshold is None:
            raise ValueError(
                unit_id + ": prediction record carries no "
                "kpoint_convergence_threshold (the per-atom k-point "
                "flatness tolerance the producer resolves and stamps "
                "on the record, DESIGN 7.8)")
        structure = load_structure(grid[0].structure)

        # Collapse duplicate-mesh rungs, then pick the converged
        #   grid point on the distinct-mesh grid (DESIGN 7.8 step
        #   3c).  Two rungs that resolved to the same mesh are one
        #   calculation run twice; their zero energy delta would
        #   fool the two-sided test, so they are merged first.
        #   ``kept[chosen]`` maps the collapsed index back to its
        #   original grid position.
        collapsed_densities, collapsed_energies, kept = \
            collapse_by_mesh(kpoint_densities, energies, meshes)
        chosen = pick_converged(
            collapsed_energies, structure.num_atoms,
            kpoint_threshold)

        # e. No flat interior point -- energy still moving at the
        #    top of the range, or the grid collapsed below the
        #    three distinct meshes the interior test needs.  Tag
        #    the flight and skip: a non-converged sweep earns no
        #    entry.
        if chosen is None:
            tag_prediction_mismatch(workspace_root, unit_id)
            summaries.append(
                unit_id + ": no converged point (energy still "
                "moving, or too few distinct meshes) -- skipped "
                "(tagged prediction_mismatch)")
            continue
        idx = kept[chosen]

        # f/g. Build the rich entry from the already-chosen facts and
        #   stage it (the shared chosen-facts core, DESIGN 5.7).  The
        #   density sweep picked its rung above (collapse_by_mesh +
        #   pick_converged); it now hands build_entry the collapsed
        #   distinct-mesh ladder, the chosen rung's k-density, and its
        #   result.toml -- the SAME shape the producer's in-memory
        #   climb hands it, so the two paths stage identical entries
        #   and cannot drift on a schema change (7.8 3f).
        # This path has no climb to take a verdict FROM, so it makes
        #   the multi-rung reading itself, over the gaps of every
        #   point in its own grid -- the same any-rung rule the climb
        #   applies to its ladder, on the evidence this path holds.
        #   The gaps were parsed in step (b) and simply never looked
        #   at before: only the chosen rung's was consulted, and one
        #   rung's apparent gap on a discrete mesh is close to a coin
        #   toss (DESIGN 1.6 / 7.8 d').
        #   A record with no threshold at all is left to build_entry,
        #   which names the missing field and what it is for; reading
        #   it here would raise a bare comparison error instead.
        ladder_gap_cut = prediction.get("metal_gap_threshold")
        ladder_is_metal = ladder_gap_cut is not None and any(
            is_gapless_value(result_toml.get("gap_ev"), ladder_gap_cut)
            for result_toml in result_tomls)
        entry = build_entry(
            workspace_root, grid[0].structure, prediction,
            dataspace, structure, kpoint_threshold,
            collapsed_densities, collapsed_energies,
            kpoint_densities[idx], result_tomls[idx],
            ladder_is_metal=ladder_is_metal)

        # g'. A metal builds no entry (DESIGN 7.8): build_entry
        #    returns None and BOTH harvest paths skip on it, so the
        #    one place the rule lives is the one builder they share.
        if entry is None:
            summaries.append(
                unit_id + ": metal -- no guidance entry staged")
            continue

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
                    "staged historical-guidance entries.")
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


def record_command():
    """Append the issued command line to a file named "command" in
    the current directory, so the exact invocation can be recovered
    later.  This is a standing project convention: each run appends
    a dated block, so the file builds up a history of how the script
    was called."""

    with open("command", "a") as cmd:
        now = datetime.now()
        stamp = now.strftime("%b. %d, %Y: %H:%M:%S")
        cmd.write(f"Date: {stamp}\n")
        cmd.write("Cmnd:")
        for argument in sys.argv:
            cmd.write(f" {argument}")
        cmd.write("\n\n")


if __name__ == "__main__":
    record_command()
    sys.exit(main())
