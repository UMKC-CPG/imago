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

Public surface (this module is built in increments; the
foundation below covers shapes + signatures)
----------------------------------------------------------------
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

Its only external dependency is ``tomllib`` (Python standard
library) plus the duck-typed ``StructureControl`` attributes
``num_atoms``, ``atom_element_name`` (1-indexed, lower-cased
symbols), and ``space_group_num``.
"""

from __future__ import annotations

import tomllib
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
w_spin = 0.5         # spin weight in the stage-2 metric
sigma_gap = 1.0      # gap normalization (eV) in the stage-2 metric
sigma_spin = 0.5     # spin normalization in the stage-2 metric
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
    dos_at_fermi: float | None        # None when absent


@dataclass(frozen=True)
class Context:
    """The calculation context an entry was converged under."""

    basis: str                        # mb | fb | eb
    functional: str                   # e.g. "gga-pbe"
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
    predicted_spin_pol: float | None


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
