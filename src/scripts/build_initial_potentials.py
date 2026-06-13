"""build_initial_potentials.py -- the augmented potential
database *producer* (DESIGN 5.7; PSEUDOCODE 11.4).

Role in the pipeline
--------------------
This script is the *producer* half of the library / producer /
consumer split documented in DESIGN 5.4.  It builds the augmented
per-element potential database files
(``share/atomicPDB/<elem>/s_gaussian_pot.toml``) that the consumer
(``makeinput.py``, C47) later reads.  Its inputs are a curated set
of reference solids; for each it refreshes the atomSCF-derived
``"isolated"`` baseline and -- once the SCF machinery lands -- runs
(or reuses a cached) Imago SCF, harvests the converged potential at
named atom sites, and writes the results via
``initial_potential_db.save()``.

The curation manifest -- what it is and why
-------------------------------------------
The pipeline's primary input is one human-readable file: the
**curation manifest** (``share/atomicBDB/manifest.toml``).  It
declares *which solids, which atom sites, under which labels, with
which SCF settings*.  Its three jobs (DESIGN 5.7):

1. **Declare the curation set.**  The manifest *is* the curation
   strategy made explicit.  Adding a reference solid means adding a
   manifest entry; reviewing the curation set means reading the
   manifest -- nothing is hidden inside the script.
2. **Tell the pipeline what to harvest.**  For each reference
   solid, which atom sites' converged potentials enter the
   database, and under what labels.
3. **Record the SCF settings used.**  k-points, convergence
   threshold, etc., copied into the provenance fields of DESIGN 5.2
   so every database entry carries the conditions of its reference
   run.

VISION Principle 5 ("the database must be regeneratable from the
curated set, not a hand-edited artifact") rules out hardcoding the
curation set inside this script and rules out folder-of-files
conventions that lose metadata.  The manifest is the smallest piece
of structured data that closes the gap: every curation choice
captured in one version-controlled file alongside the structure
files it points at, so regeneration becomes a deterministic
function of (manifest, structure files, Imago build).

Build analogy
-------------
The manifest is the **build configuration** for the database.  The
structure files are the **source**.  The Imago build is the
**toolchain**.  ``build_initial_potentials.py`` is the **build
script**.  The augmented database is the **compiled output**.  Same
role ``pyproject.toml`` plays for a Python package, or a Makefile
plays for a binary.

Reproducibility (layered)
-------------------------
- **Emitter determinism (bit-level, strict).**  The TOML emitter
  (DESIGN 5.5, in ``initial_potential_db``) writes byte-identical
  bytes for a fixed in-memory ``ElementDatabase``.
- **Pipeline numerical output (precision-level, loose).**  Given
  the same manifest, ``pot1``/``coeff1`` files, and Imago build, the
  numerical outputs agree at the precision the SCF / fit chain can
  reach -- bit-identity is *not* promised (accumulation order,
  threading, library versions can shift the last bits).
- **Provenance metadata (free).**  Timestamps and commit SHAs
  refresh every run and carry no reproducibility guarantee.

Implementation status
----------------------
This file is being built incrementally (C48):

* **C48.1 (landed) -- the manifest reader.**
  :func:`load_manifest_v2` and the manifest dataclasses.  Enforces
  validation rules 1-8 of DESIGN 5.7; rule 9 (method must be a
  registered matcher) is gated on an optional ``known_methods``
  argument and skipped when it is ``None``, exactly as
  ``initial_potential_db.load`` handles its own rule 9 -- the
  matcher registry lives in ``makeinput.py`` (C54) and does not
  exist yet.
* **C48.2 (landed) -- the ``"isolated"`` baseline refresh.**
  :func:`refresh_isolated_entries` (step 1 of the pipeline) loads
  or creates each element's database and rebuilds its
  ``"isolated"`` entry from the current ``pot1``/``coeff1`` files
  via :func:`build_isolated_entry`; :func:`save_databases` writes
  them back.  This alone forms a working producer->consumer loop
  with no SCF: every element's file gets its rule-6 baseline,
  default-tagged when the manifest curates nothing else for it.
* **C48.3 -- SCF orchestration, the content-keyed cache, COD fetch,
  the run log, and the CLI flags** (next; needs a live Imago
  toolchain to exercise).
* **C60 -- fingerprint harvest** (Phase 2): compute and attach the
  ``[[reference_solid.entry.fingerprint]]`` records.  The manifest
  reader here already *parses* those declarations and enforces
  their uniqueness (rule 8); only the harvest is deferred.
"""

import argparse
import os
import re
import subprocess
import sys
import tomllib
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any

import initial_potential_db as ipdb
import guidance_db
import guidance_harvest
from guidance_db import VALID_SYSTEM_TYPES
from kaleidoscope import Flight, LocalExecutor, SweepRecord, dispatch
from kaleidoscope.builders.kpoint_convergence import (
    build_kpoint_convergence)
from kaleidoscope.workspace import toml_line


# ============================================================
#  Manifest dataclasses (DESIGN 5.7)
# ============================================================

@dataclass
class ManifestFingerprint:
    """One ``[[reference_solid.entry.fingerprint]]`` declaration.

    Tells the producer to compute and harvest one fingerprint
    record alongside the numerical potential (the harvest itself
    is C60).  ``method`` names a matcher (ARCHITECTURE 8.9) and
    ``sub_spec`` is its method-specific parameter table.  Unlike a
    :class:`initial_potential_db.FingerprintRecord`, a declaration
    carries *no* payload: the payload is the SCF/loen output the
    producer will harvest, not something the curator writes.
    """

    method: str
    sub_spec: dict[str, Any]


@dataclass
class ReferenceEntry:
    """One ``[[reference_solid.entry]]`` harvest declaration.

    Names a single atom site whose converged potential becomes a
    labeled entry in the element's database.  ``default`` records
    whether this entry should carry the database file's single
    ``default = true`` tag (per-element-database rule 7); the
    manifest is the single source of truth for that choice.

    ``label`` is optional (DESIGN 5.2.1).  When the curator omits
    it, the producer derives the label at harvest from the run's
    site identity --
    ``<reference_id>-<element><species>-t<type>-a<site>`` -- so the
    species and type numbers (unknown until the grouping pass runs)
    land in the label without being authored ahead of time.  A
    present ``label`` is an explicit override of that derived
    default.
    """

    element: str
    atom_site: int
    default: bool
    description: str
    label: str | None = None
    fingerprints: list[ManifestFingerprint] = field(
        default_factory=list)


@dataclass
class ReferenceSolid:
    """One ``[[reference_solid]]`` -- a curated reference system.

    Carries everything the SCF run needs (the structure source,
    ``kpoint_spec``, ``scf_threshold``) plus the list of
    per-site harvest declarations.  Exactly one of ``cod_id`` (with
    ``cod_revision``) or ``structure_path`` is set (rule 4); the
    unused alternative is ``None``.

    ``system_type`` is one of the four guidance system types
    (``"crystalline"`` / ``"amorphous"`` / ``"nanostructure"`` /
    ``"molecular"``); it drives which guidance sub-model the
    predictor consults (DESIGN 7) and is recorded on the produced
    entry for forensics.  ``basis`` / ``functional`` /
    ``kpoint_integration`` are the (basis, functional, integration)
    sub-model the reference run uses; together they select the
    predictor sub-model (DESIGN 7.6) and are recorded on every
    produced entry's context.  All four are *required* in the
    manifest (rule 2): nothing the producer emits depends on an
    implicit default (VISION Principle 5).
    """

    reference_id: str
    system_type: str
    basis: str
    functional: str
    kpoint_integration: str
    kpoint_spec: dict[str, Any]
    scf_threshold: float
    cod_id: int | None
    cod_revision: str | None
    structure_path: str | None
    entries: list[ReferenceEntry] = field(default_factory=list)


@dataclass
class CurationManifest:
    """The whole curation manifest, parsed and validated.

    ``manifest_path`` is retained so later pipeline steps can
    resolve ``structure_path`` entries relative to the manifest's
    directory and report errors against the source file.
    """

    schema_version: int
    manifest_path: str
    reference_solids: list[ReferenceSolid] = field(
        default_factory=list)


# ============================================================
#  Manifest reader (DESIGN 5.7; PSEUDOCODE 11.4 load_manifest_v2)
# ============================================================

def _require(condition: bool, path: str, message: str) -> None:
    """Strict-refusal guard: raise ``ValueError`` if false.

    Mirrors the validation style of ``initial_potential_db``: every
    failure names the manifest file path and the specific rule and
    offending entry, so the curator can fix the source directly.
    There is no warning-and-continue path anywhere in the loader.
    """

    if not condition:
        raise ValueError(f"{path}: {message}")


def load_manifest_v2(path: str,
                     known_methods: set[str] | None = None
                     ) -> CurationManifest:
    """Read and validate the curation manifest (schema v2).

    Implements ``load_manifest_v2`` of PSEUDOCODE 11.4 and the nine
    validation rules of DESIGN 5.7.  Strict refusal throughout: any
    rule violation raises ``ValueError`` naming the rule, the file,
    and the offending reference solid or entry.

    The manifest file must exist (a deliberate choice -- a missing
    manifest is a hard error, not an empty curation set), so an
    absent ``path`` raises ``FileNotFoundError``.

    ``known_methods`` is the optional set of registered matcher
    names used to enforce rule 9 (a fingerprint ``method`` must be a
    registered matcher).  Callers without a registry pass ``None``
    and rule 9 is skipped -- the matcher registry lives in
    ``makeinput.py`` (ARCHITECTURE 8.9) and is not wired in until
    C54.  Rule 8 (per-entry ``(method, sub_spec)`` uniqueness) is
    enforced regardless, using the same canonical sub-spec equality
    as the per-element-database reader.

    The nine rules (DESIGN 5.7):

    1. ``schema_version`` must equal 2.
    2. Every ``[[reference_solid]]`` carries ``reference_id``,
       ``system_type``, ``basis``, ``functional``,
       ``kpoint_integration``, ``kpoint_spec``, and
       ``scf_threshold``; ``system_type`` must be one of the four
       guidance system types (``crystalline`` / ``amorphous`` /
       ``nanostructure`` / ``molecular``).
    3. Every ``[[reference_solid.entry]]`` carries ``element``,
       ``atom_site``, ``default``, ``description``.  ``label`` is
       optional (DESIGN 5.2.1): when present it overrides the
       derived default; when absent the producer assembles it at
       harvest from the run's site identity.
    4. Exactly one of ``cod_id`` / ``structure_path`` per solid;
       ``cod_id`` must be a positive integer with a non-empty
       ``cod_revision``; ``structure_path`` must resolve to an
       existing file under the manifest's directory.
    5. ``reference_id`` is unique across the manifest and is
       label-safe (lowercase letters, digits, ``-``, ``_``),
       because it is embedded verbatim in every derived entry
       label and typed into ``-pot`` (DESIGN 5.2.1).
    6. For entries with an explicit ``label``, ``(element, label)``
       is unique across the manifest.  For entries with a derived
       label, ``(reference_id, element, atom_site)`` is unique --
       two such entries would derive the identical label (the
       per-construction uniqueness of DESIGN 5.2.1).
    7. Exactly one ``default = true`` entry per element that
       appears anywhere in the manifest.
    8. Within one entry's fingerprint declarations,
       ``(method, sub_spec)`` is unique.
    9. Every fingerprint ``method`` is a registered matcher
       (checked only when ``known_methods`` is supplied).
    """

    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"curation manifest not found: {path} (a manifest "
            f"is required; create it or pass --manifest)")

    with open(path, "rb") as handle:
        raw = tomllib.load(handle)

    # ----- Rule 1: schema_version must equal 2.
    _require(raw.get("schema_version") == 2, path,
             f"manifest rule 1: schema_version must equal 2 "
             f"(found {raw.get('schema_version')!r})")

    manifest_dir = os.path.dirname(os.path.abspath(path))

    seen_ref_ids: set[str] = set()
    seen_element_label: set[tuple[str, str]] = set()
    # (reference_id, element, atom_site) for derived-label entries:
    # two with the same triple would derive the identical label.
    seen_derived_key: set[tuple[str, str, int]] = set()
    # elem -> count of default=true entries (rule 7, post-loop).
    default_per_element: dict[str, int] = {}
    seen_elements: set[str] = set()

    reference_solids: list[ReferenceSolid] = []

    for ref in raw.get("reference_solid", []):
        # ----- Rule 2: required per-solid fields.  The sub-model
        # triple (basis, functional, kpoint_integration) and
        # system_type are required alongside the run settings:
        # they select the guidance sub-model and land on every
        # produced entry, so nothing emitted rides on an implicit
        # default (VISION Principle 5).
        for field_name in ("reference_id", "system_type", "basis",
                           "functional", "kpoint_integration",
                           "kpoint_spec", "scf_threshold"):
            _require(field_name in ref, path,
                     f"manifest rule 2: [[reference_solid]] "
                     f"missing field: {field_name}")

        rid = ref["reference_id"]

        # ----- Rule 5 (charset): reference_id is embedded verbatim
        # in every derived entry label and typed into -pot, so it
        # must be label-safe -- lowercase letters, digits, hyphen,
        # underscore (DESIGN 5.2.1).  Checked before the uniqueness
        # tally so a malformed id fails with the clearer message.
        _require(re.fullmatch(r"[a-z0-9_-]+", rid) is not None, path,
                 f"manifest rule 5: reference_id {rid!r} is not "
                 f"label-safe (use lowercase letters, digits, '-', "
                 f"'_'); it is embedded verbatim in derived labels")

        # ----- Rule 2 (domain): system_type must be one of the four
        # guidance system types -- the predictor switches its
        # sub-model on it (DESIGN 7), so an unknown value is a hard
        # error rather than a silently mis-signed entry.
        _require(ref["system_type"] in VALID_SYSTEM_TYPES, path,
                 f"manifest rule 2: [[reference_solid {rid}]] "
                 f"system_type {ref['system_type']!r} is not one "
                 f"of {VALID_SYSTEM_TYPES}")

        # ----- Rule 4: exactly one structure source.
        has_cod = "cod_id" in ref
        has_path = "structure_path" in ref
        _require(has_cod != has_path, path,
                 f"manifest rule 4: [[reference_solid {rid}]] "
                 f"must set exactly one of cod_id or "
                 f"structure_path")
        if has_cod:
            # cod_id must be a positive integer and cod_revision
            # must be present and non-empty (so the build pins
            # against a specific upstream COD revision).
            cod_id = ref["cod_id"]
            _require(isinstance(cod_id, int)
                     and not isinstance(cod_id, bool)
                     and cod_id > 0, path,
                     f"manifest rule 4: cod_id must be a positive "
                     f"integer ({rid}), got {cod_id!r}")
            _require("cod_revision" in ref
                     and isinstance(ref["cod_revision"], str)
                     and len(ref["cod_revision"]) > 0, path,
                     f"manifest rule 4: cod_revision required "
                     f"(non-empty string) when cod_id is set "
                     f"({rid})")
        else:
            sp = os.path.join(manifest_dir, ref["structure_path"])
            _require(os.path.isfile(sp), path,
                     f"manifest rule 4: structure_path resolves "
                     f"to a missing file: {sp} ({rid})")

        # ----- Rule 5: reference_id uniqueness.
        _require(rid not in seen_ref_ids, path,
                 f"manifest rule 5: duplicate reference_id: {rid}")
        seen_ref_ids.add(rid)

        entries: list[ReferenceEntry] = []
        for entry in ref.get("entry", []):
            # ----- Rule 3: required entry fields.  ``label`` is
            # NOT required (DESIGN 5.2.1): absent => derived at
            # harvest, present => explicit override.
            for field_name in ("element", "atom_site",
                               "default", "description"):
                _require(field_name in entry, path,
                         f"manifest rule 3: "
                         f"[[reference_solid.entry]] in {rid} "
                         f"missing field: {field_name}")

            elem = entry["element"]
            label = entry.get("label")
            seen_elements.add(elem)

            # ----- Rule 6: no two entries may produce the same
            # database entry.  With an explicit label that means
            # (element, label) is unique across the manifest; with
            # a derived label it means (reference_id, element,
            # atom_site) is unique, since those three are exactly
            # what the derived label is built from (DESIGN 5.2.1).
            if label is not None:
                key = (elem, label)
                _require(key not in seen_element_label, path,
                         f"manifest rule 6: duplicate "
                         f"(element, label): {key}")
                seen_element_label.add(key)
            else:
                dkey = (rid, elem, entry["atom_site"])
                _require(dkey not in seen_derived_key, path,
                         f"manifest rule 6: two entries derive the "
                         f"same label (same reference_id, element, "
                         f"atom_site): {dkey}")
                seen_derived_key.add(dkey)

            # ----- Rule 7 tally: count default=true per element;
            # the "exactly one" check runs after the full walk.
            if entry["default"]:
                default_per_element[elem] = (
                    default_per_element.get(elem, 0) + 1)

            # ----- Fingerprint declarations (rules 8 and 9).
            fingerprints: list[ManifestFingerprint] = []
            seen_method_subspec: set = set()
            for fp in entry.get("fingerprint", []):
                _require("method" in fp and "sub_spec" in fp, path,
                         f"manifest rule 8: fingerprint "
                         f"declaration must carry both method and "
                         f"sub_spec ({rid}, label={label})")
                method = fp["method"]
                sub_spec = fp["sub_spec"]

                # Rule 9: method must be a registered matcher.
                # Skipped when the caller passed no registry.
                if known_methods is not None:
                    _require(method in known_methods, path,
                             f"manifest rule 9: unknown matcher "
                             f"method {method!r} ({rid}, "
                             f"label={label})")

                # Rule 8: (method, sub_spec) unique within this
                # entry, using the per-element-database reader's
                # canonical sub-spec equality (so key order and
                # int-vs-float spelling do not matter).
                canon = ipdb.canonicalize_sub_spec(sub_spec)
                fp_key = (method, canon)
                _require(fp_key not in seen_method_subspec, path,
                         f"manifest rule 8: duplicate "
                         f"(method, sub_spec) in entry {label} "
                         f"of {rid}")
                seen_method_subspec.add(fp_key)

                fingerprints.append(ManifestFingerprint(
                    method=method, sub_spec=dict(sub_spec)))

            entries.append(ReferenceEntry(
                element=elem,
                atom_site=entry["atom_site"],
                label=label,
                default=entry["default"],
                description=entry["description"],
                fingerprints=fingerprints))

        reference_solids.append(ReferenceSolid(
            reference_id=rid,
            system_type=ref["system_type"],
            basis=ref["basis"],
            functional=ref["functional"],
            kpoint_integration=ref["kpoint_integration"],
            kpoint_spec=dict(ref["kpoint_spec"]),
            scf_threshold=ref["scf_threshold"],
            cod_id=ref.get("cod_id"),
            cod_revision=ref.get("cod_revision"),
            structure_path=ref.get("structure_path"),
            entries=entries))

    # ----- Rule 7 post-loop: exactly one default-tagged entry
    # per element that appears anywhere in the manifest.
    for elem in sorted(seen_elements):
        count = default_per_element.get(elem, 0)
        _require(count == 1, path,
                 f"manifest rule 7: element {elem} has {count} "
                 f"default-tagged entries across the manifest "
                 f"(need exactly one)")

    return CurationManifest(
        schema_version=raw["schema_version"],
        manifest_path=os.path.abspath(path),
        reference_solids=reference_solids)


# ============================================================
#  Isolated-baseline refresh (DESIGN 5.7 step 2; PSEUDOCODE
#  11.4 build_isolated_entry / read_pot1 / read_coeff1)
# ============================================================
#
# Step 1 of the pipeline always rebuilds every element's
# "isolated" entry from the *current* legacy pot1/coeff1 files,
# so any change in atomSCF output propagates into the augmented
# database on the next run.  The isolated entry is also the
# rule-6 baseline every per-element file must carry, and -- for
# elements the manifest does not curate -- it is the file's
# sole entry and therefore its default-tagged one (rule 7).
#
# The legacy pot1/coeff1 text format is the same fixed-line
# layout the C47 consumer materializes; these parsers are the
# inverse of makeinput._write_legacy_pot_files_from_entry.


@dataclass
class _PotFileData:
    """The scalar fields parsed out of a legacy ``pot`` file.

    The per-term coefficients and alphas live in the companion
    ``coeff`` file (see :func:`_parse_coeff_file`); this holds
    only the per-entry scalars the ``pot`` file carries.
    """

    nuclear_z: float
    nuclear_alpha: float
    covalent_radius: float
    num_gaussians: int
    alpha_min: float
    alpha_max: float


def _parse_pot_file(path: str) -> _PotFileData:
    """Parse a legacy ``pot`` file's fixed eight-line layout.

    The layout (atomSCF output, and what the C47 consumer
    regenerates) is positional::

        0  NUCLEAR_CHARGE__ALPHA
        1  <Z> <nuclear_alpha>
        2  COVALENT_RADIUS
        3  <covalent_radius>
        4  NUM_ALPHAS
        5  <num_gaussians>
        6  ALPHAS
        7  <alpha_min> <alpha_max>

    ``Z`` is written by atomSCF as a float (``79.000000``) and is
    kept as the real ``nuclear_z`` the schema (DESIGN 5.2) uses --
    nominally an integer, but Imago consumes Z as a real number.
    The four tag lines are checked so a malformed file fails with
    a clear message rather than a silent misparse.
    """

    with open(path) as handle:
        lines = [ln.rstrip("\n") for ln in handle
                 if ln.strip() != ""]

    if len(lines) < 8:
        raise ValueError(
            f"{path}: malformed pot file (expected 8 content "
            f"lines, found {len(lines)})")
    for index, tag in ((0, "NUCLEAR_CHARGE__ALPHA"),
                       (2, "COVALENT_RADIUS"),
                       (4, "NUM_ALPHAS"), (6, "ALPHAS")):
        if lines[index].strip() != tag:
            raise ValueError(
                f"{path}: malformed pot file (line {index} "
                f"expected tag {tag!r}, found "
                f"{lines[index].strip()!r})")

    charge_alpha = lines[1].split()
    alpha_range = lines[7].split()
    return _PotFileData(
        nuclear_z=float(charge_alpha[0]),
        nuclear_alpha=float(charge_alpha[1]),
        covalent_radius=float(lines[3].split()[0]),
        num_gaussians=int(lines[5].split()[0]),
        alpha_min=float(alpha_range[0]),
        alpha_max=float(alpha_range[1]))


def _parse_coeff_file(path: str) -> tuple[list[float],
                                          list[float]]:
    """Parse a legacy ``coeff`` file into (coefficients, alphas).

    The file is a count line followed by one line per Gaussian
    term; the term lines carry five whitespace-separated columns
    of which only column 1 (the coefficient) and column 2 (the
    alpha) are meaningful here -- columns 3-5 are the placeholder
    fields Imago ignores (see C47).  The count line is
    cross-checked against the number of term lines so a truncated
    or padded file is caught at parse time.
    """

    with open(path) as handle:
        lines = [ln for ln in handle if ln.strip() != ""]

    declared_count = int(lines[0].split()[0])
    coefficients: list[float] = []
    alphas: list[float] = []
    for term_line in lines[1:]:
        tokens = term_line.split()
        coefficients.append(float(tokens[0]))
        alphas.append(float(tokens[1]))

    if len(coefficients) != declared_count:
        raise ValueError(
            f"{path}: coeff file count line says {declared_count} "
            f"terms but {len(coefficients)} term lines follow")
    return coefficients, alphas


def element_path(pdb_root: str, elem: str) -> str:
    """Return the per-element database path under ``pdb_root``.

    Matches both the consumer's lookup (C47) and PSEUDOCODE
    11.4 ``element_path``: the element directory name is
    lower-cased, so ``Au`` and ``au`` resolve to the same
    ``share/atomicPDB/au/s_gaussian_pot.toml``.
    """

    return os.path.join(pdb_root, elem.lower(),
                        "s_gaussian_pot.toml")


def is_isolated_default_for(elem: str,
                            manifest: CurationManifest) -> bool:
    """True iff the manifest curates no default entry for ``elem``.

    The isolated baseline carries the per-element file's
    ``default = true`` tag exactly when the manifest declares no
    other default-tagged entry for the element (PSEUDOCODE 11.4
    ``is_isolated_default_for``).  Manifest rule 7 forbids zero
    defaults for any element that *appears* in the manifest, so a
    manifest entry with ``default = true`` always wins over the
    baseline; an element absent from the manifest has only its
    isolated entry, which must therefore be the default.  Element
    symbols are compared case-insensitively because the manifest
    uses proper case (``Au``) while directory names are lower
    case (``au``).
    """

    for solid in manifest.reference_solids:
        for entry in solid.entries:
            if (entry.element.lower() == elem.lower()
                    and entry.default):
                return False
    return True


def build_isolated_entry(pdb_root: str, elem: str, commit: str,
                         timestamp: str,
                         manifest: CurationManifest
                         ) -> ipdb.PotentialEntry:
    """Build the ``"isolated"`` entry from current pot1/coeff1.

    Reads ``<pdb_root>/<elem>/pot1`` and ``coeff1``, cross-checks
    that the term counts agree, and returns a fresh
    :class:`initial_potential_db.PotentialEntry` tagged
    ``"isolated"`` with atomSCF-source provenance.  The
    ``default`` flag is computed from the manifest via
    :func:`is_isolated_default_for`.  No fingerprints are
    attached -- the baseline never participates in environment-
    scheme matching.
    """

    elem_dir = os.path.join(pdb_root, elem.lower())
    pot = _parse_pot_file(os.path.join(elem_dir, "pot1"))
    coefficients, alphas = _parse_coeff_file(
        os.path.join(elem_dir, "coeff1"))

    if not (len(coefficients) == len(alphas)
            == pot.num_gaussians):
        raise ValueError(
            f"{elem}: pot1/coeff1 disagree on term count "
            f"(pot num_gaussians={pot.num_gaussians}, "
            f"coeff coefficients={len(coefficients)}, "
            f"alphas={len(alphas)})")

    symbol = elem.capitalize()
    return ipdb.PotentialEntry(
        label="isolated",
        default=is_isolated_default_for(elem, manifest),
        description=(f"Single isolated {symbol} atom "
                     f"(from atomSCF)."),
        num_gaussians=pot.num_gaussians,
        alpha_min=pot.alpha_min,
        alpha_max=pot.alpha_max,
        coefficients=coefficients,
        alphas=alphas,
        provenance={
            "source": "atomSCF",
            "commit": commit,
            "generated_at": timestamp},
        fingerprints=[])


def list_element_dirs(pdb_root: str) -> list[str]:
    """Return the element directory names under ``pdb_root``.

    Only directories that actually carry a ``pot1`` file are
    returned, so non-element siblings (a ``cache`` directory, a
    stray ``manifest.toml``) are skipped.  The list is sorted for
    deterministic processing order.
    """

    names = []
    for name in os.listdir(pdb_root):
        if os.path.isfile(os.path.join(pdb_root, name, "pot1")):
            names.append(name)
    return sorted(names)


def refresh_isolated_entries(pdb_root: str,
                             manifest: CurationManifest,
                             commit: str, timestamp: str,
                             elements: list[str] | None = None
                             ) -> dict[str, ipdb.ElementDatabase]:
    """Step 1 of the pipeline: rebuild every isolated baseline.

    For each element (all element directories under ``pdb_root``,
    or just ``elements`` when given), load the existing
    ``s_gaussian_pot.toml`` if present or create a fresh
    :class:`initial_potential_db.ElementDatabase` from the
    element's ``pot1`` scalars, then drop and re-insert the
    ``"isolated"`` entry rebuilt from the current pot1/coeff1.
    Returns the in-memory databases keyed by element directory
    name; the caller saves them (see :func:`save_databases`).

    Implements PSEUDOCODE 11.4 step 1.  The existing-file load
    passes ``known_methods=None`` (rule 9 is skipped) because the
    matcher registry does not exist until C54; any fingerprint
    records already in the file are preserved untouched -- only
    the isolated entry is rewritten here.
    """

    if elements is None:
        elements = list_element_dirs(pdb_root)

    databases: dict[str, ipdb.ElementDatabase] = {}
    for elem in elements:
        path = element_path(pdb_root, elem)
        if os.path.isfile(path):
            database = ipdb.load(path, known_methods=None)
        else:
            pot = _parse_pot_file(
                os.path.join(pdb_root, elem.lower(), "pot1"))
            database = ipdb.ElementDatabase(
                schema_version=2,
                element_symbol=elem.capitalize(),
                nuclear_z=pot.nuclear_z,
                nuclear_alpha=pot.nuclear_alpha,
                covalent_radius=pot.covalent_radius,
                potentials=[])

        # Drop any prior isolated entry and re-insert a fresh one
        # at the front so atomSCF refreshes always propagate.
        database.potentials = [
            entry for entry in database.potentials
            if entry.label != "isolated"]
        database.potentials.insert(0, build_isolated_entry(
            pdb_root, elem, commit, timestamp, manifest))
        databases[elem] = database

    return databases


def save_databases(databases: dict[str, ipdb.ElementDatabase],
                   pdb_root: str) -> None:
    """Write each element database to its on-disk path.

    Implements PSEUDOCODE 11.4 step 3: every affected
    :class:`initial_potential_db.ElementDatabase` is written via
    the deterministic emitter (DESIGN 5.5) to its
    ``element_path``.  The element directory already exists (it
    held the ``pot1`` we read), so no directory creation is
    needed.
    """

    for elem, database in databases.items():
        ipdb.save(database, element_path(pdb_root, elem))


# ============================================================
#  Build identity, workspace, and structure materialization
#  (DESIGN 5.7; PSEUDOCODE 11.4)
# ============================================================

def _git_sha() -> str:
    """The current HEAD commit, or ``"unknown"`` when git is
    unavailable.  Injected into the run options as ``imago_commit``
    so kaleidoscope's run-reuse cache key (DESIGN 6.2.5) and every
    produced entry's provenance both record which build the
    *producer* believed it ran.  This records the producer's belief,
    which can drift from the binary actually executed; C84 hardens
    it by having Imago stamp its own build commit."""

    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True)
        return completed.stdout.strip() or "unknown"
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _now_iso8601_utc() -> str:
    """The build timestamp in the schema's ISO-8601 UTC form."""

    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def curation_workspace_root(pdb_root: str) -> str:
    """The kaleidoscope workspace the producer dispatches its
    combined flight into (PSEUDOCODE 11.4).  It sits beside the
    databases under the shared data root so its run-reuse cache
    (DESIGN 6.2.5) persists across producer runs and reference
    solids dedupe against earlier builds."""

    data_root = os.path.dirname(pdb_root.rstrip("/"))
    return os.path.join(data_root, "curation", "workspace")


def _cod_extension(ref: ReferenceSolid) -> str:
    """The on-disk extension for a fetched COD structure.  COD
    serves CIF, so v1 always writes ``.cif``; kept as a hook so a
    future format negotiation has a single place to change."""

    return ".cif"


def _fetch_cod_structure(cod_id: int, cod_revision: str,
                         dest: str) -> None:
    """Fetch one pinned COD revision to ``dest`` (DESIGN 5.7,
    Option A).  Strict on failure: a network outage, a COD outage,
    or a missing pinned revision raises -- the fetch NEVER falls
    back to a different revision, because a silent fallback would
    desync the reproducible build from the pinned manifest.

    NOTE (C74 end-to-end): the live COD fetch needs network access
    and the COD per-revision API; it is exercised only on the
    cluster.  A ``structure_path`` manifest avoids it entirely for
    offline / unit-test runs."""

    import urllib.error
    import urllib.request

    url = f"https://www.crystallography.net/cod/{cod_id}.cif"
    try:
        with urllib.request.urlopen(url) as response:
            payload = response.read()
    except urllib.error.URLError as exc:
        raise RuntimeError(
            f"COD fetch failed for cod_id={cod_id} (pinned "
            f"revision {cod_revision!r}): {exc}.  The build pins "
            f"this revision and refuses to fall back.") from exc
    with open(dest, "wb") as handle:
        handle.write(payload)


def materialize_structure(ref: ReferenceSolid, manifest_dir: str,
                          pdb_root: str) -> str:
    """Guarantee the reference solid's structure exists as a local
    file and return its path (Option A; DESIGN 5.7 / PSEUDOCODE
    11.4).  This is the producer's ONLY network access and is
    deliberately decoupled from any run cache: it carries no SCF
    state and makes no hit/miss decision.  Recompute avoidance
    belongs to kaleidoscope's run-reuse cache (DESIGN 6.2.5), which
    keys on this file's bytes.

    A ``structure_path`` ref is a plain disk read, resolved under
    the manifest directory (rule 4 already validated it exists).  A
    ``cod_id`` ref fetches the pinned revision once into a plain
    cache location (decoupled from the SCF cache), then reuses it."""

    if ref.structure_path is not None:
        return os.path.join(manifest_dir, ref.structure_path)

    data_root = os.path.dirname(pdb_root.rstrip("/"))
    local = os.path.join(
        data_root, "atomicBDB", "cache", "structures",
        ref.reference_id + _cod_extension(ref))
    if not os.path.exists(local):
        os.makedirs(os.path.dirname(local), exist_ok=True)
        _fetch_cod_structure(ref.cod_id, ref.cod_revision, local)
    return local


# ============================================================
#  Per-solid options + the (deferred) loen fingerprint units
# ============================================================

def make_producer_options(ref: ReferenceSolid,
                          imago_commit: str) -> dict[str, Any]:
    """The fixed (non-swept) makeinput options for a reference
    solid's convergence flight (PSEUDOCODE 11.4).

    Carries the sub-model triple the builder requires (``basis`` /
    ``functional`` / ``kpoint_integration``), the ``scf_threshold``
    and ``imago_commit`` that form kaleidoscope's run-reuse cache
    key (DESIGN 6.2.5), and the fixed k-point shift.  The swept
    k-density is deliberately absent: the builder adds it per grid
    point (the ``kpd`` option), so a fixed value here would collide
    with the sweep."""

    options: dict[str, Any] = {
        "basis": ref.basis,
        "functional": ref.functional,
        "kpoint_integration": ref.kpoint_integration,
        "scf_threshold": ref.scf_threshold,
        "imago_commit": imago_commit,
    }
    shift = ref.kpoint_spec.get("shift")
    if shift is not None:
        options["kpoint_shift"] = shift
    return options


def build_loen_units(ref: ReferenceSolid, struct_path: str) -> list:
    """Structure-only ``imago -loen -scf no`` units, one per
    Fortran-side fingerprint declaration (PSEUDOCODE 11.4).

    DEFERRED (C54 matcher registry / C60 fingerprint harvest): the
    Python-vs-Fortran split needs the MATCHERS registry to know
    which declarations require a loen run, and the fingerprint
    harvest itself is C60.  Until both land the producer dispatches
    no loen units and harvests no fingerprints -- the convergence
    path (the C74 deliverable) is unaffected.  Returns an empty
    list, tagged here so the wiring point is obvious when C54/C60
    fill it in."""

    return []


def harvest_fingerprints(flight: Flight, ref: ReferenceSolid,
                         spec: ReferenceEntry,
                         struct_path: str) -> list:
    """Build one FingerprintRecord per declared fingerprint
    (PSEUDOCODE 11.4 ``harvestFingerprints``).

    DEFERRED (C60): the Python-side / Fortran-side matcher harvest
    needs the C54 matcher registry and the C60 harvest itself.  For
    C74 it returns an empty list so a produced entry carries no
    fingerprints yet; the convergence potential -- C74's deliverable
    -- is harvested fully."""

    return []


# ============================================================
#  Harvest helpers: pick the converged run, read the potential
#  (DESIGN 5.7 / 7.8; PSEUDOCODE 11.4)
# ============================================================

def _read_unit_result(workspace_root: str, unit) -> dict:
    """Parse one dispatched unit's ``result.toml`` from its run
    directory (the ``kaleidoscope.workspace.unit_run_dir`` layout:
    ``<root>/wingbeats/<id>/<calc...>/result.toml``)."""

    path = os.path.join(workspace_root, "wingbeats", unit.id,
                        *unit.calc, "result.toml")
    with open(path, "rb") as handle:
        return tomllib.load(handle)


def pick_converged_unit(flight: Flight, reference_id: str,
                        workspace_root: str, scf_threshold: float):
    """Return ``(unit, result_toml)`` for the converged grid point
    of one solid's convergence sweep, or ``None`` when the energy is
    still moving at the top of the grid (PSEUDOCODE 11.4).

    Uses the same two-sided delta-below-threshold rule as the
    guidance harvest (DESIGN 7.8 step 3c), reused via
    ``guidance_harvest.pick_converged``: the smallest interior
    k-density whose total energy is within ``scf_threshold`` of both
    neighbours.  A single-point grid (a trust-mode run or a
    single-point curator override) has no interior point to judge,
    so its lone run IS the deliverable -- the producer trusts the
    pinned/predicted point and harvests its potential directly."""

    units = [unit for unit in flight.units
             if unit.id == reference_id
             and unit.kind == "convergence"]
    units = sorted(
        units,
        key=lambda unit: guidance_harvest.swept_value_of(
            unit, "kpt-density"))
    results = [_read_unit_result(workspace_root, unit)
               for unit in units]

    if len(units) == 1:
        return units[0], results[0]

    energies = [result["total_energy"] for result in results]
    index = guidance_harvest.pick_converged(energies, scf_threshold)
    if index is None:
        return None
    return units[index], results[index]


def extract_potential(result_toml: dict, atom_site: int
                      ) -> tuple[list[float], list[float]]:
    """Harvest the converged Gaussian potential for one atom site
    from a dispatched run's result (DESIGN 5.7 / ARCHITECTURE 9.7;
    PSEUDOCODE 11.4).

    The converged ``scfV`` output (``result.outputs["scfV"]``, the
    ``<edge>_initPot-<basis>.dat`` Imago writes from ``fort.8``)
    holds the self-consistent potential as Gaussian terms.  Per the
    5.7 harvest contract -- "converged ``scfV`` matches input
    ``scfV``: coefficients from the output, alphas from the input,
    taken together" -- the producer reads the coefficients out of
    that file; the file shares the legacy ``coeff`` layout, so the
    same parser yields both the coefficients and the alphas.

    NOTE (C74 end-to-end): the per-atom-site selection within a
    multi-site solid's ``scfV`` is validated only against a live
    Imago run; ``atom_site`` is threaded through here for that
    wiring (the single-site case reads the whole file)."""

    scfv_path = result_toml["outputs"]["scfV"]
    coefficients, alphas = _parse_coeff_file(scfv_path)
    return coefficients, alphas


def read_site_identity_map(result_toml: dict
                           ) -> dict[int, tuple[str, int, int]]:
    """Read the run's ``datSkl.map`` into ``{skeleton_atom:
    (element, species, type)}`` (DESIGN 5.2.1 / 5.7; ARCHITECTURE
    9.7).

    makeinput writes this file during input preparation, recording
    -- per atom -- the sorted-dat number, the original skeleton
    number, and the site's element symbol, OLCAO species number,
    and potential-type number (C87).  The producer keys harvest by
    ``atom_site``, which is a *skeleton* numbering index, so this
    returns the map keyed by the skeleton column (``values[1]``).

    The columns are ``DAT#  SKELETON#  ELEMENT  SPECIES  TYPE``; the
    header line is skipped.  Reading the whole file once per
    converged solid lets every site's label be assembled without a
    second pass."""

    path = result_toml["outputs"]["datSkl_map"]
    identity: dict[int, tuple[str, int, int]] = {}
    with open(path) as handle:
        rows = [line for line in handle if line.strip()]
    for row in rows[1:]:                       # skip the header line
        columns = row.split()
        skeleton_atom = int(columns[1])
        element = columns[2]
        species = int(columns[3])
        type_number = int(columns[4])
        identity[skeleton_atom] = (element, species, type_number)
    return identity


def assemble_entry_label(reference_id: str, element: str,
                         species: int, type_number: int,
                         atom_site: int) -> str:
    """Assemble the DESIGN 5.2.1 entry label
    ``<reference_id>-<element><species>-t<type>-a<site>``.

    The element symbol is lowercased so it fuses with the species
    number into the OLCAO species token the CLI speaks (``si1``),
    and the whole label is lowercase.  ``reference_id`` is already
    label-safe (manifest rule 5), so no further escaping is
    needed."""

    return (f"{reference_id}-{element.lower()}{species}"
            f"-t{type_number}-a{atom_site}")


def make_imago_provenance(commit: str, timestamp: str,
                          ref: ReferenceSolid, atom_site: int,
                          scf_iterations) -> dict[str, Any]:
    """The ``[potential.provenance]`` block for a harvested
    Imago-source entry (DESIGN 5.2 / 5.7).

    Carries the ``source = "Imago"`` discriminant and every field
    ``initial_potential_db.require_provenance`` demands of an Imago
    entry (``reference_id``, ``atom_site``, ``kpoint_spec``,
    ``scf_threshold``, ``scf_iterations``) so the 5.8 validation
    harness can re-run the originating SCF, plus ``system_type``
    recorded for forensics (5.7 rule 2)."""

    return {
        "source": "Imago",
        "commit": commit,
        "generated_at": timestamp,
        "reference_id": ref.reference_id,
        "system_type": ref.system_type,
        "atom_site": atom_site,
        "kpoint_spec": dict(ref.kpoint_spec),
        "scf_threshold": ref.scf_threshold,
        "scf_iterations": scf_iterations,
    }


# ============================================================
#  Run log (DESIGN 5.7 / 5.8; PSEUDOCODE 11.4 write_run_log)
# ============================================================

def make_run_log_entry(ref: ReferenceSolid, unit,
                       result_toml: dict) -> dict[str, Any]:
    """One converged-solid row for the run log: the reference id,
    the converged k-density (read off the chosen unit's calc tag),
    and the SCF iteration count the 5.8 harness reads."""

    return {
        "reference_id": ref.reference_id,
        "converged": True,
        "converged_kpoint_density": guidance_harvest.swept_value_of(
            unit, "kpt-density"),
        "scf_iterations": result_toml.get("scf_iterations"),
    }


def make_nonconverged_log_entry(ref: ReferenceSolid
                                ) -> dict[str, Any]:
    """One non-converged-solid row for the run log: the sweep never
    flattened, so no potential was harvested and the curator must
    widen the grid (DESIGN 7.9)."""

    return {"reference_id": ref.reference_id, "converged": False}


def write_run_log(path: str, imago_commit: str, timestamp: str,
                  per_run_log: list[dict[str, Any]]) -> None:
    """Write the producer's run log (PSEUDOCODE 11.4): a manifest
    snapshot header (the Imago commit + timestamp) followed by one
    ``[[run]]`` block per reference solid.  The 5.8 validation
    harness reads this to know which solids converged and in how
    many SCF iterations.  Emitted with the kaleidoscope
    ``toml_line`` helper so scalars/arrays format consistently."""

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        handle.write(toml_line("schema_version", 1))
        handle.write(toml_line("imago_commit", imago_commit))
        handle.write(toml_line("generated_at", timestamp))
        for entry in per_run_log:
            handle.write("\n[[run]]\n")
            for key, value in entry.items():
                if value is not None:
                    handle.write(toml_line(key, value))


def curation_executor(force: bool = False):
    """Build the executor kaleidoscope dispatches the producer's
    combined flight on (PSEUDOCODE 11.4).  v1 uses the in-process
    ``LocalExecutor``; ``force`` is threaded for the planned
    cache-bypass wiring (DESIGN 6.2.5) -- when set, cached runs
    should re-run rather than dedupe.

    NOTE (C74 end-to-end): a cluster producer swaps in a
    ``ParslExecutor``; the force-driven cache bypass is wired
    against the live run-reuse cache only on the cluster."""

    return LocalExecutor()


# ============================================================
#  The producer pipeline (DESIGN 5.7; PSEUDOCODE 11.4)
# ============================================================

def build_initial_potentials(manifest_path: str, pdb_root: str,
                             data_root: str, *, force: bool = False,
                             single_element: str | None = None,
                             dispatch_fn=dispatch,
                             extract_fn=extract_potential,
                             identity_fn=read_site_identity_map
                             ) -> list[dict[str, Any]]:
    """The three-phase producer (DESIGN 5.7; PSEUDOCODE 11.4):
    *build* a single combined flight, *dispatch* it through
    kaleidoscope, then *harvest* each solid's converged potential
    and contribute the same converged point back to the guidance
    dataspace.  Returns the per-run log (also written to disk).

    ``dispatch_fn``, ``extract_fn``, and ``identity_fn`` are
    injected (defaulting to the real kaleidoscope dispatch, the real
    scfV reader, and the real ``datSkl.map`` reader) so the
    orchestration can be unit-tested with the toolchain seam mocked:
    end-to-end dispatch, the per-site ``scfV`` read, and the
    site-identity read all need a live Imago run (C74)."""

    manifest = load_manifest_v2(manifest_path)
    manifest_dir = os.path.dirname(manifest.manifest_path)
    guidance_root = os.path.join(data_root, "historicalGuidanceDB")
    dataspace = guidance_db.load(guidance_root)
    imago_commit = _git_sha()
    timestamp = _now_iso8601_utc()
    workspace = curation_workspace_root(pdb_root)

    # ----- Phase 1: build.  Refresh the isolated baselines, then
    # accumulate every solid's convergence units (and any loen
    # units) into ONE combined flight so the whole producer run
    # dispatches as a single flat parallel batch.
    elements = None if single_element is None else [single_element]
    databases = refresh_isolated_entries(
        pdb_root, manifest, imago_commit, timestamp, elements)

    all_units: list = []
    predictions: dict[str, Any] = {}
    struct_of: dict[str, str] = {}
    for ref in manifest.reference_solids:
        struct = materialize_structure(ref, manifest_dir, pdb_root)
        struct_of[ref.reference_id] = struct
        options = make_producer_options(ref, imago_commit)
        # One builder call per solid (DESIGN 6.2.9): a pinned
        # kpoint_spec.density is the curator override (predictor
        # bypassed); otherwise full predict-then-verify.
        flight_i, record_i = build_kpoint_convergence(
            struct, options, dataspace, ref.system_type,
            id=ref.reference_id,
            center=ref.kpoint_spec.get("density"))
        all_units.extend(flight_i.units)
        all_units.extend(build_loen_units(ref, struct))
        # Store the plain dict (metadata must be TOML-serializable).
        predictions[ref.reference_id] = asdict(record_i)

    flight = Flight(
        root=workspace, units=all_units,
        sweep=SweepRecord(varied_axes=("kpt-density",),
                          fixed_axes={}),
        metadata={"predictions": predictions})

    # ----- Phase 2: dispatch.  Kaleidoscope runs and tracks every
    # unit through the wingbeat seam and its run-reuse cache.
    dispatch_fn(flight, executor=curation_executor(force))

    # ----- Phase 3: harvest.  Per solid, pick the converged grid
    # point, extract the potential at each named site, and record
    # the run.  Non-converged solids are logged and skipped.
    per_run_log: list[dict[str, Any]] = []
    for ref in manifest.reference_solids:
        converged = pick_converged_unit(
            flight, ref.reference_id, workspace, ref.scf_threshold)
        if converged is None:
            per_run_log.append(make_nonconverged_log_entry(ref))
            continue
        unit, result_toml = converged
        per_run_log.append(
            make_run_log_entry(ref, unit, result_toml))
        scf_iterations = result_toml.get("scf_iterations")

        # The site-identity map (datSkl.map) is read at most once per
        # converged solid, and only if some entry needs a derived
        # label -- entries that pin an explicit label never touch it.
        site_identity: dict[int, tuple[str, int, int]] | None = None

        for spec in ref.entries:
            elem_key = spec.element.lower()
            if elem_key not in databases:
                continue        # filtered out by --element
            coefficients, alphas = extract_fn(
                result_toml, spec.atom_site)

            # The label is the curator's explicit override when given,
            # else derived from the run's site identity (DESIGN 5.2.1).
            if spec.label is not None:
                label = spec.label
            else:
                if site_identity is None:
                    site_identity = identity_fn(result_toml)
                _elem, species, type_number = (
                    site_identity[spec.atom_site])
                label = assemble_entry_label(
                    ref.reference_id, spec.element, species,
                    type_number, spec.atom_site)

            new_entry = ipdb.PotentialEntry(
                label=label,
                default=spec.default,
                description=spec.description,
                num_gaussians=len(coefficients),
                alpha_min=min(alphas),
                alpha_max=max(alphas),
                coefficients=coefficients,
                alphas=alphas,
                provenance=make_imago_provenance(
                    imago_commit, timestamp, ref, spec.atom_site,
                    scf_iterations),
                fingerprints=harvest_fingerprints(
                    flight, ref, spec, struct_of[ref.reference_id]))
            database = databases[elem_key]
            # Replace any prior entry with the same label (manifest
            # rule 6 makes the only possible prior the same entry
            # from a previous run).
            database.potentials = [
                entry for entry in database.potentials
                if entry.label != label]
            database.potentials.append(new_entry)

    # ----- Phase 3b: guidance contribution.  The same converged
    # grid points feed the historical-guidance dataspace staging,
    # so every solid the producer converges sharpens the predictor.
    guidance_harvest.harvest_flight(
        workspace, guidance_root, dataspace)

    # ----- Write outputs: every affected element file, plus the
    # run log the 5.8 validation harness reads.
    save_databases(databases, pdb_root)
    write_run_log(
        os.path.join(data_root, "curation", "run_log.toml"),
        imago_commit, timestamp, per_run_log)
    return per_run_log


# ============================================================
#  Command-line interface
# ============================================================

def _default_pdb_root() -> str:
    """``$IMAGO_DATA/atomicPDB`` (DESIGN 5.4 layout), or empty when
    $IMAGO_DATA is unset so the parser can demand ``--pdb-root``."""

    data_dir = os.environ.get("IMAGO_DATA", "")
    return os.path.join(data_dir, "atomicPDB") if data_dir else ""


def main(argv=None) -> int:
    """CLI entry point: run the producer over a curation manifest
    (DESIGN 5.7).  ``--element`` restricts the run to one element's
    database; ``--force`` bypasses kaleidoscope's run-reuse cache so
    every reference run re-executes."""

    parser = argparse.ArgumentParser(
        description="Build the augmented initial-potential database "
                    "from a curation manifest (DESIGN 5.7).")
    parser.add_argument(
        "--manifest", required=True,
        help="path to the curation manifest (schema v2)")
    parser.add_argument(
        "--pdb-root", default=_default_pdb_root(),
        help="the atomicPDB root (default: $IMAGO_DATA/atomicPDB)")
    parser.add_argument(
        "--element", default=None,
        help="restrict the build to this one element's database")
    parser.add_argument(
        "--force", action="store_true",
        help="bypass the run-reuse cache so every run re-executes")
    args = parser.parse_args(argv)

    if not args.pdb_root:
        parser.error("--pdb-root not given and $IMAGO_DATA is unset")
    data_root = os.path.dirname(args.pdb_root.rstrip("/"))

    per_run_log = build_initial_potentials(
        args.manifest, args.pdb_root, data_root,
        force=args.force, single_element=args.element)
    converged = sum(1 for row in per_run_log if row["converged"])
    print(f"producer: {converged}/{len(per_run_log)} reference "
          f"solids converged and harvested")
    return 0


if __name__ == "__main__":
    sys.exit(main())
