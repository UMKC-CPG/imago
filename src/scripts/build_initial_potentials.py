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

import os
import tomllib
from dataclasses import dataclass, field
from typing import Any

import initial_potential_db as ipdb


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
    """

    element: str
    atom_site: int
    label: str
    default: bool
    description: str
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
    """

    reference_id: str
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
       ``kpoint_spec``, and ``scf_threshold``.
    3. Every ``[[reference_solid.entry]]`` carries ``element``,
       ``atom_site``, ``label``, ``default``, ``description``.
    4. Exactly one of ``cod_id`` / ``structure_path`` per solid;
       ``cod_id`` must be a positive integer with a non-empty
       ``cod_revision``; ``structure_path`` must resolve to an
       existing file under the manifest's directory.
    5. ``reference_id`` is unique across the manifest.
    6. ``(element, label)`` is unique across the manifest.
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
    # elem -> count of default=true entries (rule 7, post-loop).
    default_per_element: dict[str, int] = {}
    seen_elements: set[str] = set()

    reference_solids: list[ReferenceSolid] = []

    for ref in raw.get("reference_solid", []):
        # ----- Rule 2: required per-solid fields.
        for field_name in ("reference_id", "kpoint_spec",
                            "scf_threshold"):
            _require(field_name in ref, path,
                     f"manifest rule 2: [[reference_solid]] "
                     f"missing field: {field_name}")

        rid = ref["reference_id"]

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
            # ----- Rule 3: required entry fields.
            for field_name in ("element", "atom_site", "label",
                               "default", "description"):
                _require(field_name in entry, path,
                         f"manifest rule 3: "
                         f"[[reference_solid.entry]] in {rid} "
                         f"missing field: {field_name}")

            elem = entry["element"]
            label = entry["label"]
            seen_elements.add(elem)

            # ----- Rule 6: (element, label) uniqueness across
            # the entire manifest -- two solids cannot both
            # produce the same database entry.
            key = (elem, label)
            _require(key not in seen_element_label, path,
                     f"manifest rule 6: duplicate "
                     f"(element, label): {key}")
            seen_element_label.add(key)

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
