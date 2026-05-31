"""guidance_db.py -- Historical-guidance dataspace file-format
and predictor library.

Role in the pipeline
--------------------
This module is the *library* half of the library / producer /
consumer split documented in DESIGN 7 and ARCHITECTURE 10.  It
is the one place that knows how to read, validate, write, and
predict from the historical-guidance dataspace -- the curated
collection of converged calculations whose chemistry signature,
electronic-structure character, and converged k-point density
let a small predictor turn a new system into a predicted
operating point plus an uncertainty (VISION Goal 5).

The dataspace lives under ``share/historicalGuidanceDB/`` and is
partitioned by ``system_type``:

    share/historicalGuidanceDB/
      SCHEMA_VERSION              bare-integer marker file
      elemental_groups.toml       element -> group table
      entries/<system_type>/*.toml    promoted entries
      staging/<system_type>/*.toml    harvested, awaiting curator

The producers (``guidance_harvest.py``, ``guidance_promote.py``)
and the consumer (the kaleidoscope flight-builder helper) all
go through this module; keeping the file format and the
predictor isolated here means a schema bump touches one file
(ARCHITECTURE 10) and unit tests can exercise the format with
synthetic structures instead of real Imago runs.

Why a dataspace and a predictor (not a categorical lookup) and
the full predict-then-verify motivation are in DESIGN 7.1.

Schema version
--------------
This module implements **schema version 1** (DESIGN 7.2).  The
``SCHEMA_VERSION`` marker file at the dataspace root and every
entry's ``schema_version`` key must equal :data:`SCHEMA_VERSION`;
:func:`load` rejects anything else (a future bump is handled by
``guidance_migrate.py``, not by silent coercion).

Public surface
--------------
Signature, Measured, Context, Verification, Provenance,
GuidanceEntry, Dataspace, PredictionResult
    The in-memory records mirroring DESIGN 7.4 field-for-field.
load_elemental_groups(path)
    Read and invert the checked-in element-to-group table into a
    lower-cased ``symbol -> group`` lookup.
compute_signature(structure, system_type, group_table, label)
    Turn a structure into the predictor's feature input: the
    13-dimensional atom-fraction composition vector plus the
    crystalline lattice-family one-hot.
bravais_family_of(structure)
    Map a structure's space-group number to one of the six
    v1 lattice families.
load(root)
    Read and validate the whole dataspace under ``root`` into a
    :class:`Dataspace`, partitioned by system_type.
load_entry(path, system_type_dir, seen_ids)
    Read and validate one entry file against DESIGN 7.2 rules
    1-12, with file/block/field error messages.
save_entry(entry, root) / format_entry(entry, slug) / slug_for
    Emit an entry to ``staging/<system_type>/<slug>.toml`` with
    the byte-deterministic, %.16e hand-formatter; the slug is
    ``<system_type>-<short_sha>`` over the provenance fields.
predict(dataspace, query, basis, functional, kpoint_integration)
    Predict a converged k-density and confidence for a new
    system: the canonical entry for non-crystalline systems, the
    two-stage k-NN (select_submodel -> stage1 -> stage2) for
    crystalline ones.  Always returns a PredictionResult.

Its only external dependency is ``tomllib`` (Python standard
library) plus the duck-typed ``StructureControl`` attributes
``num_atoms``, ``atom_element_name`` (1-indexed, lower-cased
symbols), and ``space_group_num``.
"""

from __future__ import annotations

import glob
import hashlib
import math
import os
import tomllib
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any


# ============================================================
#  Constants and canonical orderings (DESIGN 7.4 / PSEUDOCODE
#  15.1)
# ============================================================

# Schema + partition constants.  Bumping SCHEMA_VERSION is a
# deliberate, migration-backed event (guidance_migrate.py).
SCHEMA_VERSION = 1
VALID_SYSTEM_TYPES = (
    "crystalline", "amorphous", "nanostructure", "molecular")
NON_CRYSTALLINE_TYPES = ("amorphous", "nanostructure", "molecular")
VALID_BASES = ("mb", "fb", "eb")
VALID_GAP_KINDS = ("direct", "indirect", "none")
METRIC_REGISTRY = ("total_energy",)        # DESIGN 7.2 rule 10

# Canonical slot orderings.  These pin which composition-vector
# slot means which element group, and which one-hot slot means
# which Bravais family, so the reader, the emitter,
# compute_signature, and the predictor all agree.  The
# composition vector sums to 1.0 (DESIGN 7.2 rule 4).
CANONICAL_GROUP_ORDER = (       # 13 element groups
    "alkali", "alkali_earth", "halide", "chalcogen",
    "pnictogen", "group_iv", "group_iii", "transition_metal",
    "lanthanide", "actinide", "metalloid", "noble_gas",
    "hydrogen")
CANONICAL_LATTICE_ORDER = (     # 6 Bravais families
    "cubic", "hex", "tet", "ortho", "mono", "tri")

# Predictor tuning knobs (DESIGN 7.6).  Named here so that a
# post-seed recalibration is a one-file change.  (Consumed by
# the predictor, added in a later increment of this module.)
k_min = 3            # below this a sub-model is refused
k_neighbors = 5      # neighbors used at each k-NN stage
epsilon = 1.0e-6     # numerical floor on a distance
w_comp = 1.0         # composition weight in the stage-1 metric
w_latt = 0.25        # lattice-family weight in the stage-1 metric
w_gap = 1.0          # gap weight in the stage-2 metric
w_spin = 0.5         # magnetization weight in the stage-2 metric
sigma_gap = 1.0      # gap normalization (eV) in the stage-2 metric
sigma_spin = 0.5     # magnetization normalization (Bohr magnetons
#                      per atom) in the stage-2 metric
sigma_gap_ref = 1.0  # gap spread -> confidence_1 reference (eV)
sigma_kpd_ref = 50.0  # kpd spread -> confidence_2 reference


# ============================================================
#  Space-group-number to lattice-family mapping (DESIGN 7.4 /
#  PSEUDOCODE 15.2)
# ============================================================
#
# StructureControl carries the International Tables space-group
# NUMBER (``space_group_num``), not a crystal-system label, so
# the crystal system is read off the standard IT number ranges
# below.  The v1 family mapping then lumps trigonal into ``hex``
# (six families, not seven) -- a simplification flagged in
# DESIGN 7.10; isolating it in one table makes a future split a
# one-line change.
#
# IT space-group-number ranges -> crystal system:
##   1- 2  triclinic       75-142  tetragonal
##   3- 15 monoclinic     143-167  trigonal
##  16- 74 orthorhombic   168-194  hexagonal
##                        195-230  cubic

CRYSTAL_SYSTEM_TO_FAMILY = {
    "triclinic":    "tri",   "monoclinic":   "mono",
    "orthorhombic": "ortho", "tetragonal":   "tet",
    "trigonal":     "hex",   "hexagonal":    "hex",
    "cubic":        "cubic"}


def _crystal_system_of(space_group_num: int) -> str:
    """Return the crystal-system name for an IT space-group
    number, using the standard number ranges.  Raises
    ``ValueError`` if the number is outside 1..230 (which
    includes the ``0`` that StructureControl uses to mean "no
    space group recorded").
    """

    sg = space_group_num
    if not (1 <= sg <= 230):
        raise ValueError(
            "space_group_num " + str(sg) + " is outside the "
            "valid International Tables range 1..230 (0 means "
            "no space group was recorded for the structure)")
    if sg <= 2:
        return "triclinic"
    if sg <= 15:
        return "monoclinic"
    if sg <= 74:
        return "orthorhombic"
    if sg <= 142:
        return "tetragonal"
    if sg <= 167:
        return "trigonal"
    if sg <= 194:
        return "hexagonal"
    return "cubic"


# ============================================================
#  In-memory records (DESIGN 7.4 / PSEUDOCODE 15.1)
# ============================================================
#
# These mirror the on-disk schema block-for-block and field
# order, so the reader (a later increment) and the emitter have
# a single target.  The vector fields are tuples (immutable) and
# the records are frozen, so an entry cannot drift after it is
# validated.

@dataclass(frozen=True)
class Signature:
    """The predictor's feature input for one system."""

    system_type: str                  # one of VALID_SYSTEM_TYPES
    composition_vector: tuple[float, ...]   # 13, CANONICAL_GROUP_ORDER
    lattice_family: str               # "" if non-crystalline
    lattice_onehot: tuple[float, ...]  # 6, CANONICAL_LATTICE_ORDER;
    #                                    all zeros if non-crystalline


@dataclass(frozen=True)
class Measured:
    """The quantities harvested from a converged calculation."""

    gap_ev: float
    gap_kind: str                     # direct | indirect | none
    spin_polarization: float
    total_magnetization: float
    kpoint_density: float             # the predictor's target


@dataclass(frozen=True)
class Context:
    """The calculation context an entry was converged under."""

    basis: str                        # mb | fb | eb
    functional: str                   # e.g. "gga-pbe"
    kpoint_integration: str           # e.g. "gaussian-0.1"
    scf_threshold: float
    cell_atom_count: int
    cell_volume_per_formula_unit: float    # Bohr^3


@dataclass(frozen=True)
class Verification:
    """The convergence grid that backs a campaign-sourced entry.

    ``grid_energies`` is parallel to ``grid_values`` and lets the
    curator's auto-promote rule read flatness from a staging file
    alone (DESIGN 7.8); it is optional for hand-entered records.
    """

    grid_values: tuple[float, ...]
    grid_energies: tuple[float, ...] | None
    converged_at: float
    metric: str                       # "total_energy"
    metric_threshold: float
    predictor_confidence: float       # [0.0, 1.0]
    predictor_neighbor_ids: tuple[str, ...]


@dataclass(frozen=True)
class Provenance:
    """Where an entry came from."""

    flight_id: str
    source_structure: str
    imago_commit: str
    curator: str


@dataclass(frozen=True)
class GuidanceEntry:
    """One dataspace entry: a fully-validated converged record."""

    entry_id: str
    generated_at: str                 # ISO-8601 UTC
    source: str                       # "flight" | "manual"
    signature: Signature
    measured: Measured
    context: Context
    verification: Verification | None  # None only for manual
    provenance: Provenance


@dataclass
class Dataspace:
    """The whole loaded dataspace, partitioned by system_type."""

    schema_version: int
    entries_by_system_type: dict[str, list[GuidanceEntry]]
    group_table: dict[str, str]       # lower-cased symbol -> group


@dataclass(frozen=True)
class PredictionResult:
    """What :func:`predict` returns (a later increment)."""

    predicted_kpoint_density: float
    confidence: float                 # [0.0, 1.0]
    is_under_trained: bool
    neighbor_entry_ids: tuple[str, ...]
    predicted_gap: float | None       # None if non-crystalline
    predicted_magnetization: float | None  # intensive moment,
    #                                  Bohr magnetons per atom;
    #                                  None if non-crystalline


# ============================================================
#  Element-group table (DESIGN 7.4 / PSEUDOCODE 15.2)
# ============================================================

def load_elemental_groups(path: str) -> dict[str, str]:
    """Read ``elemental_groups.toml`` and invert it into a
    lower-cased ``symbol -> group`` lookup.

    The table is checked-in *data*, not code (VISION Principle
    11).  Element symbols are stored in conventional case in the
    file (e.g. ``"Si"``) but matched case-insensitively, because
    StructureControl stores atom element names lower-cased; the
    returned lookup is therefore keyed by the lower-cased symbol.

    A symbol that lands in two groups is a data-file typo and is
    rejected loudly (it must fail, not let the last assignment
    silently win).  Every group in :data:`CANONICAL_GROUP_ORDER`
    must be present -- even ``metalloid``, which ships empty in
    v1 (DESIGN 7.4 / 7.10).
    """

    with open(path, "rb") as handle:
        raw = tomllib.load(handle)

    version = raw.get("schema_version")
    if version != SCHEMA_VERSION:
        raise ValueError(
            path + ": elemental_groups.toml schema_version "
            + repr(version) + " != " + str(SCHEMA_VERSION))
    if "groups" not in raw:
        raise ValueError(
            path + ": elemental_groups.toml missing [groups]")

    groups = raw["groups"]
    table: dict[str, str] = {}        # lower-cased symbol -> group
    for group in CANONICAL_GROUP_ORDER:
        if group not in groups:
            raise ValueError(
                path + ": elemental_groups.toml missing group: "
                + group)
        for symbol in groups[group]:
            key = symbol.lower()
            if key in table:
                raise ValueError(
                    path + ": element " + symbol + " assigned to"
                    " two groups (" + table[key] + ", " + group
                    + ")")
            table[key] = group
    return table


# ============================================================
#  Structure -> Signature (DESIGN 7.4 / PSEUDOCODE 15.2)
# ============================================================

def bravais_family_of(structure: Any) -> str:
    """Map a structure's space-group number to one of the six v1
    lattice families (cubic / hex / tet / ortho / mono / tri).

    Reuses StructureControl's recorded ``space_group_num`` (the
    International Tables number) rather than re-deriving the
    Bravais lattice from cell vectors.  A structure with no
    recorded space group (``space_group_num == 0``) cannot be
    signed as crystalline and raises a clear error.
    """

    system = _crystal_system_of(structure.space_group_num)
    return CRYSTAL_SYSTEM_TO_FAMILY[system]


def compute_signature(structure: Any, system_type: str,
                      group_table: dict[str, str],
                      label: str = "<structure>") -> Signature:
    """Turn a structure into the predictor's feature input.

    The composition vector counts atoms per element group and
    normalizes to atom fractions, laid out in
    :data:`CANONICAL_GROUP_ORDER`.  The lattice family (and its
    one-hot) is read only for crystalline systems; for the three
    non-crystalline ``system_type`` values the family is the
    empty string and the one-hot is all zeros, so the
    predictor's stage-1 distance never sees a lattice term for
    them (DESIGN 7.4).

    ``structure`` is duck-typed: it needs ``num_atoms``,
    ``atom_element_name`` (1-indexed, lower-cased symbols, index
    0 unused), and -- for crystalline systems --
    ``space_group_num``.  ``label`` only sharpens error messages.

    An element symbol missing from ``group_table`` is a hard
    error here, at compute time, naming the offending structure
    -- distinct from a malformed dataspace load.
    """

    if system_type not in VALID_SYSTEM_TYPES:
        raise ValueError(
            "unknown system_type: " + repr(system_type))

    # Composition vector: count atoms per group, then normalize.
    counts = {group: 0 for group in CANONICAL_GROUP_ORDER}
    total_atoms = 0
    for index in range(1, structure.num_atoms + 1):
        symbol = structure.atom_element_name[index]
        if symbol is None:
            raise ValueError(
                "atom " + str(index) + " in structure " + label
                + " has no element name")
        key = symbol.lower()
        if key not in group_table:
            raise ValueError(
                "element " + symbol + " (atom " + str(index)
                + " in structure " + label + ") is not in "
                "elemental_groups.toml")
        counts[group_table[key]] += 1
        total_atoms += 1
    if total_atoms <= 0:
        raise ValueError(
            "structure " + label + " has no atoms")
    composition = tuple(
        counts[group] / total_atoms
        for group in CANONICAL_GROUP_ORDER)

    # Lattice family + one-hot, crystalline only.
    if system_type == "crystalline":
        family = bravais_family_of(structure)
        onehot = tuple(
            1.0 if name == family else 0.0
            for name in CANONICAL_LATTICE_ORDER)
    else:
        family = ""
        onehot = tuple(0.0 for _ in CANONICAL_LATTICE_ORDER)

    return Signature(
        system_type=system_type,
        composition_vector=composition,
        lattice_family=family,
        lattice_onehot=onehot)


# ============================================================
#  TOML reader load() (DESIGN 7.2 rules 1-12 / PSEUDOCODE 15.3)
# ============================================================

def _require(condition: bool, path: str, message: str) -> None:
    """Raise a ``ValueError`` naming the file at fault when a
    validation condition fails.  This mirrors the
    file/block/field error discipline of
    ``initial_potential_db.load`` (DESIGN 5.2): a caller reading
    the message learns exactly which file and field broke a
    rule, without a stack dive.
    """

    if not condition:
        raise ValueError(path + ": " + message)


def load(root: str) -> Dataspace:
    """Read and validate the whole dataspace under ``root``.

    Reads the bare-integer ``SCHEMA_VERSION`` marker, the
    ``elemental_groups.toml`` table, and every entry under
    ``entries/<system_type>/``, validating each file against the
    twelve rules (DESIGN 7.2) and partitioning the result by
    system_type.  A ``system_type`` subdirectory that does not
    exist simply contributes no entries -- an empty dataspace is
    valid (it is exactly the cold-start state, DESIGN 7.9).
    """

    # Rule 1 (marker half): the bare-integer marker file at the
    # dataspace root must equal SCHEMA_VERSION.
    marker_path = os.path.join(root, "SCHEMA_VERSION")
    with open(marker_path, "r") as handle:
        marker = handle.read().strip()
    _require(marker == str(SCHEMA_VERSION), marker_path,
             "marker " + marker + " != " + str(SCHEMA_VERSION))

    group_table = load_elemental_groups(
        os.path.join(root, "elemental_groups.toml"))

    entries_by_type: dict[str, list[GuidanceEntry]] = {
        system_type: [] for system_type in VALID_SYSTEM_TYPES}
    seen_ids: dict[str, str] = {}        # entry_id -> source path
    for system_type in VALID_SYSTEM_TYPES:
        subdir = os.path.join(root, "entries", system_type)
        if not os.path.isdir(subdir):
            continue
        for path in sorted(glob.glob(
                os.path.join(subdir, "*.toml"))):
            entry = load_entry(path, system_type, seen_ids)
            entries_by_type[system_type].append(entry)

    return Dataspace(
        schema_version=SCHEMA_VERSION,
        entries_by_system_type=entries_by_type,
        group_table=group_table)


def load_entry(path: str, system_type_dir: str,
               seen_ids: dict[str, str]) -> GuidanceEntry:
    """Read and validate one entry file against DESIGN 7.2 rules
    1-12, returning a :class:`GuidanceEntry`.

    The schema is fully checked BEFORE any dataclass is built
    (rule 12), so an omission surfaces as a clear validation
    failure rather than a bare constructor ``TypeError`` -- the
    same discipline as DESIGN 5.2.  ``system_type_dir`` is the
    directory the file was found under (rule 3 requires the
    entry's own ``system_type`` to match it); ``seen_ids`` maps
    each entry_id seen so far to its file so a duplicate names
    both files (rule 2).
    """

    with open(path, "rb") as handle:
        raw = tomllib.load(handle)

    # Rule 12 (top-level half): the required top-level keys.
    for field_name in ("schema_version", "entry_id",
                       "generated_at", "source"):
        _require(field_name in raw, path,
                 "missing top-level: " + field_name)

    # Rule 1 (entry half): version agrees with the marker.
    _require(raw["schema_version"] == SCHEMA_VERSION, path,
             "schema_version " + str(raw["schema_version"])
             + " != " + str(SCHEMA_VERSION))

    # Rule 11 (source domain): flight | manual.
    source = raw["source"]
    _require(source in ("flight", "manual"), path,
             "source must be flight|manual, got " + str(source))

    # Rule 2: entry_id unique across the whole entries tree.
    entry_id = raw["entry_id"]
    _require(entry_id not in seen_ids, path,
             "duplicate entry_id " + str(entry_id) + " (also in "
             + seen_ids.get(entry_id, "?") + ")")
    seen_ids[entry_id] = path

    _require("entry" in raw, path, "missing [entry]")
    entry = raw["entry"]

    # --- signature block --------------------------------------
    _require("signature" in entry, path,
             "missing [entry.signature]")
    sig = entry["signature"]
    _require("system_type" in sig, path,
             "missing signature.system_type")
    system_type = sig["system_type"]

    # Rule 3: system_type valid AND matches the directory.
    _require(system_type in VALID_SYSTEM_TYPES, path,
             "invalid system_type: " + str(system_type))
    _require(system_type == system_type_dir, path,
             "system_type " + system_type + " under entries/"
             + system_type_dir + "/")

    # Rule 4: composition_vector has exactly the 13 group keys,
    # each in [0, 1], summing to 1.0 +/- 1e-6.
    _require("composition_vector" in sig, path,
             "missing signature.composition_vector")
    cv = sig["composition_vector"]
    _require(set(cv.keys()) == set(CANONICAL_GROUP_ORDER), path,
             "composition_vector keys != the 13 element groups")
    composition = tuple(cv[g] for g in CANONICAL_GROUP_ORDER)
    for group, fraction in zip(CANONICAL_GROUP_ORDER, composition):
        _require(0.0 <= fraction <= 1.0, path,
                 "composition_vector[" + group + "] out of [0,1]")
    _require(abs(sum(composition) - 1.0) <= 1.0e-6, path,
             "composition_vector sums to " + str(sum(composition))
             + " != 1.0")

    # Rule 5: lattice_family present + valid iff crystalline.
    family = sig.get("lattice_family", "")
    if system_type == "crystalline":
        _require(family in CANONICAL_LATTICE_ORDER, path,
                 "crystalline entry: lattice_family must be one of"
                 " the six families, got '" + str(family) + "'")
        onehot = tuple(
            1.0 if name == family else 0.0
            for name in CANONICAL_LATTICE_ORDER)
    else:
        _require(family == "", path,
                 "non-crystalline entry must not set lattice_family")
        onehot = tuple(0.0 for _ in CANONICAL_LATTICE_ORDER)

    signature = Signature(system_type, composition, family, onehot)

    # --- measured block ---------------------------------------
    _require("measured" in entry, path, "missing [entry.measured]")
    m = entry["measured"]
    for field_name in ("gap_ev", "gap_kind", "spin_polarization",
                       "total_magnetization", "kpoint_density"):
        _require(field_name in m, path,
                 "missing measured." + field_name)

    # Rule 6: gap_ev >= 0; gap_kind valid; "none" iff metal.
    _require(m["gap_ev"] >= 0.0, path, "gap_ev < 0")
    _require(m["gap_kind"] in VALID_GAP_KINDS, path,
             "invalid gap_kind: " + str(m["gap_kind"]))
    is_metal = (m["gap_ev"] == 0.0)
    _require((m["gap_kind"] == "none") == is_metal, path,
             "gap_kind=='none' iff gap_ev==0.0 violated")
    # Rule 7: kpoint_density > 0.
    _require(m["kpoint_density"] > 0.0, path,
             "kpoint_density must be > 0")

    measured = Measured(
        gap_ev=m["gap_ev"],
        gap_kind=m["gap_kind"],
        spin_polarization=m["spin_polarization"],
        total_magnetization=m["total_magnetization"],
        kpoint_density=m["kpoint_density"])

    # --- context block ----------------------------------------
    _require("context" in entry, path, "missing [entry.context]")
    c = entry["context"]
    for field_name in ("basis", "functional",
                       "kpoint_integration", "scf_threshold",
                       "cell_atom_count",
                       "cell_volume_per_formula_unit"):
        _require(field_name in c, path,
                 "missing context." + field_name)
    # Rule 8: basis valid; functional + kpoint_integration
    # non-empty.
    _require(c["basis"] in VALID_BASES, path,
             "invalid basis: " + str(c["basis"]))
    _require(len(c["functional"]) > 0, path,
             "functional must be non-empty")
    _require(len(c["kpoint_integration"]) > 0, path,
             "kpoint_integration must be non-empty")
    # Rule 9: cell counts / volumes positive.
    _require(c["cell_atom_count"] > 0, path,
             "cell_atom_count must be > 0")
    _require(c["cell_volume_per_formula_unit"] > 0.0, path,
             "cell_volume_per_formula_unit must be > 0")
    context = Context(
        basis=c["basis"],
        functional=c["functional"],
        kpoint_integration=c["kpoint_integration"],
        scf_threshold=c["scf_threshold"],
        cell_atom_count=c["cell_atom_count"],
        cell_volume_per_formula_unit=c[
            "cell_volume_per_formula_unit"])

    # --- verification block (required for source=flight) ------
    verification = None
    if "verification" in entry:
        verification = load_verification(
            entry["verification"], measured, path)
    _require(verification is not None or source == "manual", path,
             "source=flight requires [entry.verification]")

    # --- provenance block -------------------------------------
    _require("provenance" in entry, path,
             "missing [entry.provenance]")
    p = entry["provenance"]
    for field_name in ("flight_id", "source_structure",
                       "imago_commit", "curator"):
        _require(field_name in p, path,
                 "missing provenance." + field_name)
    if source == "flight":
        # Rule 11: flight entries need a non-empty flight_id,
        # source_structure, and imago_commit.
        for field_name in ("flight_id", "source_structure",
                           "imago_commit"):
            _require(len(p[field_name]) > 0, path,
                     "source=flight needs non-empty " + field_name)
    provenance = Provenance(
        flight_id=p["flight_id"],
        source_structure=p["source_structure"],
        imago_commit=p["imago_commit"],
        curator=p["curator"])

    return GuidanceEntry(
        entry_id=entry_id,
        generated_at=raw["generated_at"],
        source=source,
        signature=signature,
        measured=measured,
        context=context,
        verification=verification,
        provenance=provenance)


def load_verification(raw_verification: dict[str, Any],
                      measured: Measured,
                      path: str) -> Verification:
    """Validate and build the ``[entry.verification]`` block
    (DESIGN 7.2 rule 10).

    The grid must be sorted ascending, ``converged_at`` must be a
    grid point and must equal ``measured.kpoint_density`` (the
    converged density and the recorded target are the same
    number), the metric must be registered, and the confidence
    must lie in [0, 1].  ``grid_energies`` is optional but, when
    present, must be parallel to ``grid_values`` so the curator's
    auto-promote rule can read flatness from the staging file
    alone (DESIGN 7.8).
    """

    v = raw_verification
    for field_name in ("grid_values", "converged_at", "metric",
                       "metric_threshold", "predictor_confidence",
                       "predictor_neighbor_ids"):
        _require(field_name in v, path,
                 "missing verification." + field_name)

    grid = v["grid_values"]
    _require(list(grid) == sorted(grid), path,
             "grid_values not sorted ascending")
    _require(v["converged_at"] in grid, path,
             "converged_at not present in grid_values")
    _require(v["converged_at"] == measured.kpoint_density, path,
             "converged_at != measured.kpoint_density")
    _require(v["metric"] in METRIC_REGISTRY, path,
             "unknown metric: " + str(v["metric"]))
    _require(0.0 <= v["predictor_confidence"] <= 1.0, path,
             "predictor_confidence out of [0,1]")

    energies = v.get("grid_energies")
    if energies is not None:
        _require(len(energies) == len(grid), path,
                 "grid_energies length != grid_values length")

    return Verification(
        grid_values=tuple(grid),
        grid_energies=(tuple(energies)
                       if energies is not None else None),
        converged_at=v["converged_at"],
        metric=v["metric"],
        metric_threshold=v["metric_threshold"],
        predictor_confidence=v["predictor_confidence"],
        predictor_neighbor_ids=tuple(
            v["predictor_neighbor_ids"]))


# ============================================================
#  Hand-formatted emitter save_entry() (DESIGN 7.5 /
#  PSEUDOCODE 15.4)
# ============================================================
#
# Deterministic hand-formatter, same philosophy as
# initial_potential_db.save: a fixed block sequence, a fixed
# key order, %.16e floats, and float arrays one element per
# line with a trailing comma.  The output is byte-identical for
# a given in-memory entry, so version-control diffs are
# meaningful and a save/reload round trip is lossless for
# IEEE-754 binary64 values.  Within each block the `=` signs are
# aligned one space past the longest key in that block (matching
# the DESIGN 7.3 gold sketch); the float arrays are emitted
# unpadded as multi-line blocks.


def _fmt_float(value: float) -> str:
    """Render a float in the canonical 16-significant-digit form
    (one digit before the point, sixteen after, signed two-digit
    exponent) -- e.g. ``5.0000000000000000e+01``.  This is the
    representation that round-trips exactly through a decimal
    string for binary64.
    """

    return format(value, ".16e")


def _toml_string(text: str) -> str:
    """Quote a string as a TOML basic string, escaping backslash
    and double quote (the only escapes the schema's string
    fields can need)."""

    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _scalar(value: Any) -> str:
    """Render one scalar: floats as %.16e, ints bare, strings
    TOML-quoted.  ``bool`` is guarded before ``int`` because it
    is an ``int`` subclass (no boolean fields exist today, but
    the guard keeps a future one honest)."""

    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return _fmt_float(value)
    if isinstance(value, int):
        return str(value)
    return _toml_string(value)


def _emit_scalars(lines: list[str],
                  pairs: list[tuple[str, Any]],
                  width: int | None = None) -> None:
    """Append ``key = value`` lines with the ``=`` signs aligned.
    By default the alignment width is the longest key among
    ``pairs``; ``width`` overrides it when a block's keys must
    align to a wider key emitted separately (the verification
    block, whose float arrays carry the same scalars as
    ``predictor_neighbor_ids``)."""

    if width is None:
        width = max(len(key) for key, _ in pairs)
    for key, value in pairs:
        lines.append(key.ljust(width) + " = " + _scalar(value))


def _emit_float_array(lines: list[str], key: str,
                      values: tuple[float, ...]) -> None:
    """Append a float array as a multi-line block: the opener
    ``key = [`` (unpadded), one ``    <float>,`` line per
    element, and a closing ``]``."""

    lines.append(key + " = [")
    for value in values:
        lines.append("    " + _fmt_float(value) + ",")
    lines.append("]")


def short_sha(flight_id: str, source_structure: str,
              generated_at: str) -> str:
    """The slug guard (DESIGN 7.5): the first six hex digits of
    the SHA-256 over the three provenance fields concatenated.
    Two simultaneous harvests differ in ``flight_id`` or
    ``source_structure``, so their hashes (and filenames) differ.
    """

    blob = (flight_id + source_structure
            + generated_at).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:6]


def slug_for(entry: GuidanceEntry) -> str:
    """The entry's filename stem and ``entry_id``:
    ``<system_type>-<short_sha>``."""

    return (entry.signature.system_type + "-"
            + short_sha(entry.provenance.flight_id,
                        entry.provenance.source_structure,
                        entry.generated_at))


def save_entry(entry: GuidanceEntry, root: str) -> str:
    """Emit ``entry`` into ``staging/<system_type>/<slug>.toml``
    under ``root`` and return the path written.

    Producers stage; a curator promotes (DESIGN 7.8), so a fresh
    harvest always lands in ``staging/``.  The on-disk
    ``entry_id`` is set to the slug here (the harvest builds the
    entry with an empty ``entry_id`` and lets the slug fill it).
    A pre-existing file at the target is refused (DESIGN 7.2 rule
    2 / 7.5): the caller retries with a fresh ``generated_at`` on
    the rare hash clash, rather than silently overwriting.
    """

    slug = slug_for(entry)
    subdir = os.path.join(root, "staging", entry.signature.system_type)
    os.makedirs(subdir, exist_ok=True)
    path = os.path.join(subdir, slug + ".toml")
    if os.path.exists(path):
        raise ValueError(
            "save_entry: " + path
            + " already exists (entry_id collision)")
    with open(path, "w") as handle:
        handle.write(format_entry(entry, slug))
    return path


def format_entry(entry: GuidanceEntry, slug: str) -> str:
    """Render ``entry`` to the canonical, byte-deterministic TOML
    text (the gold sketch of DESIGN 7.3).  ``slug`` becomes the
    on-disk ``entry_id``."""

    lines: list[str] = []

    # Top-level block.  schema_version is a bare int; entry_id
    # always equals the slug.
    _emit_scalars(lines, [
        ("schema_version", SCHEMA_VERSION),
        ("entry_id", slug),
        ("generated_at", entry.generated_at),
        ("source", entry.source)])
    lines.append("")

    # [entry.signature] -- lattice_family only when crystalline.
    lines.append("[entry.signature]")
    signature_pairs = [("system_type", entry.signature.system_type)]
    if entry.signature.system_type == "crystalline":
        signature_pairs.append(
            ("lattice_family", entry.signature.lattice_family))
    _emit_scalars(lines, signature_pairs)
    lines.append("")

    # The composition vector: all 13 groups, aligned to the
    # widest group name, in canonical order.
    lines.append("[entry.signature.composition_vector]")
    width = max(len(group) for group in CANONICAL_GROUP_ORDER)
    for group, fraction in zip(CANONICAL_GROUP_ORDER,
                               entry.signature.composition_vector):
        lines.append(group.ljust(width) + " = " + _fmt_float(fraction))
    lines.append("")

    # [entry.measured].
    lines.append("[entry.measured]")
    _emit_scalars(lines, [
        ("gap_ev", entry.measured.gap_ev),
        ("gap_kind", entry.measured.gap_kind),
        ("spin_polarization", entry.measured.spin_polarization),
        ("total_magnetization", entry.measured.total_magnetization),
        ("kpoint_density", entry.measured.kpoint_density)])
    lines.append("")

    # [entry.context].
    lines.append("[entry.context]")
    _emit_scalars(lines, [
        ("basis", entry.context.basis),
        ("functional", entry.context.functional),
        ("kpoint_integration", entry.context.kpoint_integration),
        ("scf_threshold", entry.context.scf_threshold),
        ("cell_atom_count", entry.context.cell_atom_count),
        ("cell_volume_per_formula_unit",
         entry.context.cell_volume_per_formula_unit)])
    lines.append("")

    # [entry.verification] when present.  The float arrays are
    # multi-line and unpadded; the remaining scalars align to the
    # widest of their own keys (predictor_neighbor_ids).
    if entry.verification is not None:
        v = entry.verification
        lines.append("[entry.verification]")
        _emit_float_array(lines, "grid_values", v.grid_values)
        if v.grid_energies is not None:
            _emit_float_array(lines, "grid_energies",
                              v.grid_energies)
        scalar_width = max(
            len("converged_at"), len("metric"),
            len("metric_threshold"), len("predictor_confidence"),
            len("predictor_neighbor_ids"))
        _emit_scalars(lines, [
            ("converged_at", v.converged_at),
            ("metric", v.metric),
            ("metric_threshold", v.metric_threshold),
            ("predictor_confidence", v.predictor_confidence)],
            width=scalar_width)
        # predictor_neighbor_ids: an inline string array, aligned
        # to the same width.
        ids = ", ".join(_toml_string(i)
                        for i in v.predictor_neighbor_ids)
        lines.append("predictor_neighbor_ids".ljust(scalar_width)
                     + " = [" + ids + "]")
        lines.append("")

    # [entry.provenance].
    lines.append("[entry.provenance]")
    _emit_scalars(lines, [
        ("flight_id", entry.provenance.flight_id),
        ("source_structure", entry.provenance.source_structure),
        ("imago_commit", entry.provenance.imago_commit),
        ("curator", entry.provenance.curator)])

    return "\n".join(lines) + "\n"


# ============================================================
#  Predictor predict() (DESIGN 7.6 / PSEUDOCODE 15.5)
# ============================================================
#
# Given a new system's Signature, predict its converged k-point
# density and an uncertainty.  Non-crystalline systems return
# the single hand-seeded canonical entry (k-density there is set
# by a cell-volume convention, not chemistry).  Crystalline
# systems run a two-stage k-NN: stage 1 maps chemistry to the
# electronic character (gap, intensive magnetization) the system
# is likely to produce, and stage 2 maps that predicted character
# to a k-density -- "find calcs whose gap and magnetization look
# like what this query will produce, and copy their converged
# density."
# predict always returns a PredictionResult; is_under_trained +
# confidence tell the caller (the flight-builder helper, 15.6)
# how much to widen its verification grid.


def predict(dataspace: Dataspace, query: Signature, basis: str,
            functional: str,
            kpoint_integration: str) -> PredictionResult:
    """Predict the converged k-point density for ``query``.

    ``basis``, ``functional``, and ``kpoint_integration`` select
    the predictor sub-model (DESIGN 7.6 step 2): the dataspace is
    conditioned on the calculation settings the new run will use,
    so a prediction never interpolates across incompatible
    settings -- in particular a tetrahedral-integration density
    is never mixed with a Gaussian-smeared one.
    """

    pool = dataspace.entries_by_system_type.get(
        query.system_type, [])

    if query.system_type in NON_CRYSTALLINE_TYPES:
        return predict_non_crystalline(pool)

    entries, under_trained = select_submodel(
        pool, basis, functional, kpoint_integration)
    if under_trained:
        # No usable sub-model: the caller falls back to the
        # wide-grid default (DESIGN 7.9).  The density field is
        # unused in this branch.
        return PredictionResult(
            predicted_kpoint_density=0.0,
            confidence=0.0,
            is_under_trained=True,
            neighbor_entry_ids=(),
            predicted_gap=None,
            predicted_magnetization=None)

    predicted_gap, predicted_mag, conf1, ids1 = stage1(query, entries)
    predicted_kpd, conf2, ids2 = stage2(
        predicted_gap, predicted_mag, entries)
    return PredictionResult(
        predicted_kpoint_density=predicted_kpd,
        confidence=conf1 * conf2,
        is_under_trained=False,
        neighbor_entry_ids=_dedup(ids1 + ids2),
        predicted_gap=predicted_gap,
        predicted_magnetization=predicted_mag)


def predict_non_crystalline(
        pool: list[GuidanceEntry]) -> PredictionResult:
    """Return the canonical-entry prediction for a non-crystalline
    system_type.  k-density there follows a cell-volume
    convention rather than chemistry, so the single hand-seeded
    ``manual`` entry (DESIGN 7.9) carries essentially all the
    signal.  If several canonical entries accumulate, the most
    recent wins (deterministic by ``generated_at``)."""

    canonical = [entry for entry in pool if entry.source == "manual"]
    if not canonical:
        return PredictionResult(
            predicted_kpoint_density=0.0, confidence=0.0,
            is_under_trained=True, neighbor_entry_ids=(),
            predicted_gap=None, predicted_magnetization=None)
    entry = max(canonical, key=lambda e: e.generated_at)
    return PredictionResult(
        predicted_kpoint_density=entry.measured.kpoint_density,
        confidence=1.0,
        is_under_trained=False,
        neighbor_entry_ids=(entry.entry_id,),
        predicted_gap=None,
        predicted_magnetization=None)


def functional_family(functional: str) -> str:
    """The functional family: the token before the first hyphen
    (``"gga-pbe" -> "gga"``, ``"lda" -> "lda"``).  Whether
    functional/basis are sub-model dimensions or regression
    features is open (DESIGN 7.10); isolating the rule here makes
    that a one-line change."""

    return functional.split("-")[0]


def most_populous_submodel(
        entries: list[GuidanceEntry]) -> list[GuidanceEntry]:
    """Group ``entries`` by ``(basis, functional)`` and return the
    largest group (empty list if there are no entries)."""

    groups: dict[tuple[str, str], list[GuidanceEntry]] = {}
    for entry in entries:
        key = (entry.context.basis, entry.context.functional)
        groups.setdefault(key, []).append(entry)
    if not groups:
        return []
    return max(groups.values(), key=len)


def select_submodel(pool: list[GuidanceEntry], basis: str,
                    functional: str, kpoint_integration: str
                    ) -> tuple[list[GuidanceEntry], bool]:
    """Pick the entries the k-NN runs over, via the
    (basis, functional, kpoint_integration) -> functional-family
    -> overall-pool fallback chain (DESIGN 7.6 step 2).  Returns
    ``(entries, is_under_trained)``; ``is_under_trained`` is True
    only when even the whole system_type pool has fewer than
    :data:`k_min` entries.  The family/pool fallbacks ignore
    basis and integration -- they are the degraded best-effort
    path when no exact sub-model is populous enough."""

    # 1. The exact (basis, functional, kpoint_integration)
    #    sub-model.
    exact = [entry for entry in pool
             if entry.context.basis == basis
             and entry.context.functional == functional
             and entry.context.kpoint_integration
             == kpoint_integration]
    if len(exact) >= k_min:
        return exact, False

    # 2. The most-populous sub-model within the same functional
    #    family (e.g. (mb,gga-pbe) backing off to (fb,gga-pbe)).
    family_token = functional_family(functional)
    family = [entry for entry in pool
              if functional_family(entry.context.functional)
              == family_token]
    best = most_populous_submodel(family)
    if len(best) >= k_min:
        return best, False

    # 3. The whole system_type pool, settings ignored.
    if len(pool) >= k_min:
        return pool, False

    # 4. Too thin everywhere.
    return pool, True


def knn_weights(entries: list[GuidanceEntry],
                distance_of: Callable[[GuidanceEntry], float]
                ) -> list[tuple[GuidanceEntry, float]]:
    """Return the ``k_neighbors`` nearest entries as
    ``(entry, weight)`` pairs.  Weights are inverse-distance,
    ``1 / (d + epsilon)``, normalized to sum to 1.0 -- so an
    exact match (d = 0) dominates without dividing by zero."""

    scored = sorted(entries, key=distance_of)
    nearest = scored[:min(k_neighbors, len(scored))]
    raw = [1.0 / (distance_of(entry) + epsilon) for entry in nearest]
    total = sum(raw)
    return [(entry, weight / total)
            for entry, weight in zip(nearest, raw)]


def intensive_magnetization(entry: GuidanceEntry) -> float:
    """The entry's net magnetic moment per atom (|M| / N_atoms),
    in Bohr magnetons per atom -- the predictor's spin-character
    feature (DESIGN 7.6).

    Three deliberate choices make this a sound similarity feature
    where the raw cell moment is not:

    - **Per-atom (intensive).**  The total moment ``M`` is
      extensive -- a primitive cell and an N-fold supercell of the
      same magnetic material report N-fold-different totals while
      needing the SAME converged k-density (a reciprocal-space
      property of the primitive physics).  Dividing by the atom
      count puts cells of different sizes on one footing for the
      k-NN.
    - **Magnitude.**  The up/down labeling of the two spin
      channels is an arbitrary SCF choice, so only ``|M|`` carries
      physical information about integration difficulty.
    - **From the measured moment, not a polarization.**  imago
      surfaces the magnetic moment (iteration-file column 6), never
      a spin-polarization fraction, so ``Measured.spin_polarization``
      is structurally 0.0 for harvested entries; keying on the
      moment is what keeps this feature alive (the C72 decision).

    v1 limitation: an antiferromagnet has ``M = 0`` total yet is
    locally spin-polarized, so this feature reads it as
    non-magnetic (DESIGN 7.6 / 7.10)."""

    count = entry.context.cell_atom_count
    if count <= 0:
        return 0.0
    return abs(entry.measured.total_magnetization) / count


def stage1(query: Signature, entries: list[GuidanceEntry]
           ) -> tuple[float, float, float, list[str]]:
    """Chemistry -> electronic character.  The distance ``d1``
    combines the composition L2 distance with the lattice-family
    one-hot term, the latter halved so a full lattice mismatch
    maps into the same [0, 1] range as composition (DESIGN 7.6
    step 3).  Returns the inverse-distance-weighted predicted gap
    and predicted intensive magnetization (Bohr magnetons per
    atom), a confidence from the weighted gap variance, and the
    neighbor entry_ids."""

    def d1(entry: GuidanceEntry) -> float:
        comp_sq = sum(
            (a - b) ** 2 for a, b in zip(
                query.composition_vector,
                entry.signature.composition_vector))
        latt_sq = sum(
            (a - b) ** 2 for a, b in zip(
                query.lattice_onehot,
                entry.signature.lattice_onehot))
        return math.sqrt(w_comp * comp_sq + w_latt * latt_sq / 2.0)

    neighbors = knn_weights(entries, d1)
    predicted_gap = sum(
        weight * entry.measured.gap_ev
        for entry, weight in neighbors)
    predicted_mag = sum(
        weight * intensive_magnetization(entry)
        for entry, weight in neighbors)
    variance = sum(
        weight * (entry.measured.gap_ev - predicted_gap) ** 2
        for entry, weight in neighbors)
    confidence = math.exp(-math.sqrt(variance) / sigma_gap_ref)
    neighbor_ids = [entry.entry_id for entry, _ in neighbors]
    return predicted_gap, predicted_mag, confidence, neighbor_ids


def stage2(predicted_gap: float, predicted_mag: float,
           entries: list[GuidanceEntry]
           ) -> tuple[float, float, list[str]]:
    """Electronic character -> k-density.  The distance ``d2``
    uses the PREDICTED character from stage 1, not the query's
    chemistry: find calcs whose gap and intensive magnetization
    resemble what this query is likely to produce, and copy their
    converged density.  Returns the weighted predicted k-density,
    a confidence from the weighted k-density variance, and the
    neighbor ids."""

    def d2(entry: GuidanceEntry) -> float:
        gap_term = (w_gap * (predicted_gap - entry.measured.gap_ev) ** 2
                    / sigma_gap ** 2)
        mag_term = (w_spin
                    * (predicted_mag
                       - intensive_magnetization(entry)) ** 2
                    / sigma_spin ** 2)
        return math.sqrt(gap_term + mag_term)

    neighbors = knn_weights(entries, d2)
    predicted_kpd = sum(
        weight * entry.measured.kpoint_density
        for entry, weight in neighbors)
    variance = sum(
        weight * (entry.measured.kpoint_density - predicted_kpd) ** 2
        for entry, weight in neighbors)
    confidence = math.exp(-math.sqrt(variance) / sigma_kpd_ref)
    neighbor_ids = [entry.entry_id for entry, _ in neighbors]
    return predicted_kpd, confidence, neighbor_ids


def _dedup(ids: list[str]) -> tuple[str, ...]:
    """Order-preserving de-duplication of neighbor entry_ids."""

    seen: set[str] = set()
    out: list[str] = []
    for identifier in ids:
        if identifier not in seen:
            seen.add(identifier)
            out.append(identifier)
    return tuple(out)
