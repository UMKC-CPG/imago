"""Curation-manifest schema library (DESIGN 5.7).

The curation manifest is the hand-authored file that drives the
initial-potential database build: it declares the database-wide
fingerprint recipe (the ``[characterization]`` block), names each
reference solid (its structure source and SCF run settings), and
carries optional per-environment customizations.  The build itself
auto-discovers one representative per distinct environment in every
converged calculation; an entry exists only to override what that
harvest would otherwise produce.  This module is the single, neutral
home for that schema -- the dataclasses a parsed manifest becomes and
the readers that parse and validate the file.

Keeping the schema here, rather than inside any one tool, lets every
side of the pipeline share one definition without depending on the
others.  The producer (``build_initial_potentials``) imports these to
run a manifest; the authoring tool (``expand_manifest``) imports them
to write one.  This library itself depends only on the lower-level
databases it validates against (``initial_potential_db`` for canonical
sub-spec equality, ``guidance_db`` for the valid system-type set), so
nothing in the tool layer is pulled in here.

Two readers are provided.  :func:`load_manifest_v2` is the strict,
complete read the producer uses: every run and harvest field must be
present and every validation rule must hold.  :func:`load_structure_
sources` is the relaxed read behind the ``--materialize-only``
pre-flight: it enforces just enough to materialize each structure, so
a half-written manifest can still be fetched and converted before its
SCF and harvest details are filled in.
"""

import os
import re
import tomllib
from dataclasses import dataclass, field
from typing import Any

import initial_potential_db as ipdb
from guidance_db import VALID_SYSTEM_TYPES


# ============================================================
#  Manifest dataclasses (DESIGN 5.7)
# ============================================================

@dataclass
class ManifestFingerprint:
    """One fingerprint declaration -- a ``method`` plus its
    ``sub_spec`` (DESIGN 5.7).

    The same shape serves both kinds of declaration a manifest
    carries, and the *kind* is fixed by which block holds it, not
    by an authored field:

    - In the database-wide ``[characterization]`` block it is the
      family's **preferred** record -- the one ``sub_spec`` per
      method the consumer queries for a file-dictated structure
      (DESIGN 5.6.5 step 2).  Lifting the recipe into one block
      makes manifest rule 11 (a uniform preferred ``sub_spec`` per
      family) hold *structurally*: a single declaration per method
      cannot diverge between elements.
    - In a ``[[reference_solid.entry.fingerprint]]`` block it is a
      rare **non-preferred** override: an extra ``sub_spec`` to
      harvest for one environment alongside the preferred recipe.

    ``method`` names a matcher (ARCHITECTURE 8.9) and ``sub_spec``
    is its method-specific parameter table.  Unlike a
    :class:`initial_potential_db.FingerprintRecord`, a declaration
    carries *no* payload: the payload is the SCF/loen output the
    producer will harvest, not something the curator writes.

    ``preferred`` is **derived by the reader from the containing
    block**, never read from the file: ``True`` for a
    ``[characterization]`` declaration, ``False`` for a per-entry
    override (manifest rule 10 forbids authoring ``preferred`` on a
    per-entry declaration).  The producer carries the flag onto
    each produced :class:`initial_potential_db.FingerprintRecord`,
    and the writer never serializes it -- a round-trip recovers it
    from the block the declaration lands in.
    """

    method: str
    sub_spec: dict[str, Any]
    preferred: bool = False


@dataclass
class ReferenceEntry:
    """One ``[[reference_solid.entry]]`` customization (DESIGN 5.7).

    An entry is an *optional customization* on an auto-discovered
    environment, not a harvest instruction: the harvest finds one
    representative per distinct environment on its own, and an
    environment no entry mentions is harvested all the same (5.2.2,
    5.2.3).  An entry exists only to override what the harvest would
    otherwise produce, so **every field is optional**:

    - ``atom_site`` -- the representative-atom handle that addresses
      the target environment (1-based site index).  The species
      number cannot be the handle because it is unknown until the
      run's grouping pass; a representative atom is known up front.
      When given it also *pins* that atom as the environment's
      representative.
    - ``element`` -- the element symbol the customization targets,
      cross-checked against the species at ``atom_site`` once Imago
      loads the structure.
    - ``default`` -- marks this environment the element's default
      (per-element-database rule 7).  When no customization
      designates one, the producer falls back to the element's
      ``"isolated"`` baseline, so each file still gets its single
      default entry.  Absent reads as ``False``.
    - ``description`` -- one-sentence prose for the database entry.
      When omitted the producer auto-composes one from the solid's
      ``source_description``, qualified by the species and site.
    - ``label`` -- the label this environment is written under
      (DESIGN 5.2.1).  When omitted the producer derives it at
      harvest from the run's site identity
      (``<reference_id>-<element><species>-t<type>-a<site>``); a
      present ``label`` overrides that derived default.
    """

    element: str | None = None
    atom_site: int | None = None
    default: bool = False
    description: str | None = None
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
    implicit default (VISION Principles 5 and 11).

    ``source_description`` is an optional one-line, structure-level
    description of the solid (the hint cod_fish reads from the CIF).
    The harvest composes each environment's auto-description from
    it, qualified by the species and site it discovers (5.2.1); a
    per-entry ``description`` customization can still replace any
    one by hand.  ``None`` when the curator authored no hint.
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
    source_description: str | None = None
    entries: list[ReferenceEntry] = field(default_factory=list)


@dataclass
class CurationManifest:
    """The whole curation manifest, parsed and validated.

    ``manifest_path`` is retained so later pipeline steps can
    resolve ``structure_path`` entries relative to the manifest's
    directory and report errors against the source file.

    ``characterization`` is the database-wide fingerprint recipe
    lifted out of the ``[characterization]`` block (DESIGN 5.7): at
    most one declaration per method, each marked ``preferred`` by
    the reader.  It is the one ``sub_spec`` per family the consumer
    queries (5.6.5 step 2) and the producer harvests for every
    environment.  :func:`load_manifest_v2` requires at least one
    declaration here (manifest rule 2), so a loaded manifest always
    carries a recipe; the field defaults to an empty list only for
    a manifest assembled in memory before its recipe is set.
    """

    schema_version: int
    manifest_path: str
    characterization: list[ManifestFingerprint] = field(
        default_factory=list)
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

    Implements ``load_manifest_v2`` of PSEUDOCODE 11.4 and the
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
    C54.  When supplied it is checked against every method, both in
    the ``[characterization]`` block and in per-entry overrides.
    Rule 8 (per-entry ``(method, sub_spec)`` uniqueness) is enforced
    regardless, using the same canonical sub-spec equality as the
    per-element-database reader.

    The validation rules (DESIGN 5.7):

    1. ``schema_version`` must equal 2.
    2. Every ``[[reference_solid]]`` carries ``reference_id``,
       ``system_type``, ``basis``, ``functional``,
       ``kpoint_integration``, ``kpoint_spec``, and
       ``scf_threshold``; ``system_type`` must be one of the four
       guidance system types (``crystalline`` / ``amorphous`` /
       ``nanostructure`` / ``molecular``).  A top-level
       ``[characterization]`` block declaring at least one
       fingerprint is also required: it sets the database-wide
       preferred recipe, so a manifest without one is refused.
    3. ``[[reference_solid.entry]]`` blocks are *optional*
       customizations (DESIGN 5.2.2): a solid may carry none, and
       every customization field (``element``, ``atom_site``,
       ``default``, ``description``, ``label``) is optional.  An
       entry only overrides what the harvest would otherwise
       produce; an unmentioned environment is harvested all the
       same.
    4. Exactly one of ``cod_id`` / ``structure_path`` per solid;
       ``cod_id`` must be a positive integer with a non-empty
       ``cod_revision``; ``structure_path`` must resolve to an
       existing file under the manifest's directory.
    5. ``reference_id`` is unique across the manifest and is
       label-safe (lowercase letters, digits, ``-``, ``_``),
       because it is embedded verbatim in every derived entry
       label and typed into ``-pot`` (DESIGN 5.2.1).
    6. For customizations with an explicit ``label`` (and
       ``element``), ``(element, label)`` is unique across the
       manifest.  Derived labels are unique by construction (each
       auto-discovered environment mints its own from the run
       identity, DESIGN 5.2.1), so they need no cross-manifest
       check.
    7. At most one ``default = true`` customization per element at
       load time (the "exactly one per harvested element" half is
       checked at harvest, once the run reveals which elements
       exist).  An element with no default customization falls back
       to its ``"isolated"`` baseline at harvest (DESIGN 5.7).
    8. Within one entry's fingerprint declarations,
       ``(method, sub_spec)`` is unique.
    9. Every fingerprint ``method`` is a registered matcher
       (checked only when ``known_methods`` is supplied).
    10. The ``[characterization]`` block declares *at most one*
       fingerprint per ``method``; that single declaration is the
       family's database-wide preferred record.  A method named
       twice there, or a per-entry declaration marked
       ``preferred = true``, is a hard error: the preferred recipe
       has exactly one home (DESIGN 5.6.5).
    11. Holds *structurally* (no runtime check): because the
       preferred recipe is the single ``[characterization]``
       declaration per method, the preferred ``sub_spec`` for a
       family cannot diverge between elements.  Non-preferred
       per-entry overrides may use any ``sub_spec``.
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

    # ----- Rules 9, 10: the [characterization] block is the
    # database-wide preferred recipe.  At most one fingerprint per
    # method (rule 10) and each method a registered matcher (rule
    # 9, gated on known_methods).  That single declaration per
    # method IS the family's preferred record, so the reader stamps
    # preferred=True on it; per-entry overrides (parsed below) may
    # never be preferred.  Rule 11 then holds by construction: one
    # declaration per method cannot diverge between elements.
    characterization: list[ManifestFingerprint] = []
    char_methods: set[str] = set()
    char_block = raw.get("characterization", {})
    for fp in char_block.get("fingerprint", []):
        _require("method" in fp and "sub_spec" in fp, path,
                 "manifest rule 10: [characterization] fingerprint "
                 "must carry both method and sub_spec")
        method = fp["method"]
        if known_methods is not None:
            _require(method in known_methods, path,
                     f"manifest rule 9: unknown matcher method "
                     f"{method!r} in [characterization]")
        _require(method not in char_methods, path,
                 f"manifest rule 10: method {method!r} declared "
                 f"twice in [characterization] (the preferred "
                 f"recipe has exactly one home)")
        char_methods.add(method)
        characterization.append(ManifestFingerprint(
            method=method, sub_spec=dict(fp["sub_spec"]),
            preferred=True))

    # ----- Rule 2: the recipe is required.  A manifest must declare
    # a [characterization] block with at least one fingerprint, so
    # the build cannot silently produce a database with no preferred
    # descriptors for the consumer to match against (DESIGN 5.7,
    # VISION Principles 5 and 11).  The relaxed --materialize-only reader
    # (load_structure_sources) does not apply this -- it only
    # materializes structures and never harvests.
    _require(len(characterization) > 0, path,
             "manifest rule 2: a [characterization] block declaring "
             "at least one fingerprint is required (it sets the "
             "database-wide preferred recipe the consumer matches "
             "against)")

    seen_ref_ids: set[str] = set()
    seen_element_label: set[tuple[str, str]] = set()
    # elem -> count of default=true customizations (rule 7 load
    # half: at most one per element; the exactly-one-per-harvested-
    # element half is checked at harvest, DESIGN 5.7).
    default_per_element: dict[str, int] = {}

    reference_solids: list[ReferenceSolid] = []

    for ref in raw.get("reference_solid", []):
        # ----- Rule 2: required per-solid fields.  The sub-model
        # triple (basis, functional, kpoint_integration) and
        # system_type are required alongside the run settings:
        # they select the guidance sub-model and land on every
        # produced entry, so nothing emitted rides on an implicit
        # default (VISION Principles 5 and 11).
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
            # ----- Rule 3: entries are OPTIONAL customizations and
            # every field is optional (DESIGN 5.2.2).  An absent
            # element / atom_site / description / label is filled by
            # the harvest; an absent default reads as False.
            elem = entry.get("element")
            label = entry.get("label")

            # ----- Rule 6: only an EXPLICIT label (with an element
            # to key it on) needs the cross-manifest (element,
            # label) uniqueness check -- two solids cannot both
            # produce, e.g., ("Au", "default_solid").  Derived
            # labels are unique by construction (each environment
            # mints its own from the run identity, DESIGN 5.2.1), so
            # they need no check here.
            if label is not None and elem is not None:
                key = (elem, label)
                _require(key not in seen_element_label, path,
                         f"manifest rule 6: duplicate "
                         f"(element, label): {key}")
                seen_element_label.add(key)

            # ----- Rule 7 (load half): count default=true
            # customizations per element; "at most one" is checked
            # after the full walk.  Only countable when the
            # customization names its element.
            if entry.get("default", False) and elem is not None:
                default_per_element[elem] = (
                    default_per_element.get(elem, 0) + 1)

            # ----- Per-entry fingerprint declarations are RARE
            # overrides: extra NON-preferred sub_specs harvested for
            # this environment alongside the [characterization]
            # recipe.  Rules 8 (uniqueness), 9 (known method), and
            # 10 (may-not-be-preferred) apply.
            fingerprints: list[ManifestFingerprint] = []
            seen_method_subspec: set = set()
            for fp in entry.get("fingerprint", []):
                _require("method" in fp and "sub_spec" in fp, path,
                         f"manifest rule 8: fingerprint "
                         f"declaration must carry both method and "
                         f"sub_spec ({rid}, label={label})")
                method = fp["method"]
                sub_spec = fp["sub_spec"]

                # Rule 10: a per-entry declaration may not be
                # preferred -- the preferred record of each family
                # is fixed database-wide by [characterization].
                _require(not fp.get("preferred", False), path,
                         f"manifest rule 10: per-entry fingerprint "
                         f"may not be preferred ({rid}, "
                         f"label={label}); set it in "
                         f"[characterization]")

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
                    method=method, sub_spec=dict(sub_spec),
                    preferred=False))

            entries.append(ReferenceEntry(
                element=elem,
                atom_site=entry.get("atom_site"),
                label=label,
                default=bool(entry.get("default", False)),
                description=entry.get("description"),
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
            source_description=ref.get("source_description"),
            entries=entries))

    # ----- Rule 7 post-loop (load half): at most one default=true
    # customization per element.  Zero is allowed -- an element may
    # take its default at harvest, or fall back to its isolated
    # baseline (DESIGN 5.7); the "exactly one per HARVESTED element"
    # half is checked once the run reveals which elements exist.
    # Sorting keeps the message deterministic.
    for elem in sorted(default_per_element):
        count = default_per_element[elem]
        _require(count <= 1, path,
                 f"manifest rule 7: element {elem} has {count} "
                 f"default=true customizations across the manifest "
                 f"(at most one)")

    return CurationManifest(
        schema_version=raw["schema_version"],
        manifest_path=os.path.abspath(path),
        characterization=characterization,
        reference_solids=reference_solids)


# ============================================================
#  Structure-source reader (DESIGN 5.7): the relaxed pre-flight
#  read behind --materialize-only
# ============================================================

def load_structure_sources(path: str) -> list[ReferenceSolid]:
    """Parse ONLY the structure sources from a curation manifest --
    the relaxed read behind ``--materialize-only`` (DESIGN 5.7).

    The full ``load_manifest_v2`` requires every run and harvest
    field (``kpoint_spec``, ``scf_threshold``, the basis / functional
    / integration sub-model, the per-site entries and their
    fingerprint declarations).  A curator who has just pinned a fresh
    structure set does not have those yet: the natural order of work
    is to get every structure fetching and converting cleanly FIRST,
    eyeball the skeletons, and only then fill in the SCF and harvest
    details.  This reader enforces just enough to materialize each
    solid -- schema version (rule 1), a label-safe and unique
    ``reference_id`` (rule 5), and exactly one valid structure source
    (rule 4).  ``system_type`` is captured when present (it is an
    intrinsic structure property a sketch already carries) but is not
    domain-checked here; the run and harvest fields are left at
    harmless placeholders, because this relaxed view is never
    dispatched, only materialized.  That lets a half-written manifest
    still be pre-flighted, and lets the authoring tool
    (``expand_manifest``) read the same sketch back to complete it."""

    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"curation manifest not found: {path} (a manifest "
            f"is required; create it or pass --manifest)")

    with open(path, "rb") as handle:
        raw = tomllib.load(handle)

    # ----- Rule 1: schema_version must equal 2 (same gate the full
    # loader applies, so a v1 file fails identically here).
    _require(raw.get("schema_version") == 2, path,
             f"manifest rule 1: schema_version must equal 2 "
             f"(found {raw.get('schema_version')!r})")

    seen_ref_ids: set[str] = set()
    sources: list[ReferenceSolid] = []

    for ref in raw.get("reference_solid", []):
        # reference_id is the one per-solid field the pre-flight
        #   genuinely needs: it names the cached CIF and skeleton.
        _require("reference_id" in ref, path,
                 "manifest rule 2: [[reference_solid]] missing "
                 "field: reference_id")
        rid = ref["reference_id"]

        # ----- Rule 5 (charset + uniqueness): the id is embedded in
        # the cached filenames and, later, in derived labels, so it
        # must be label-safe and appear at most once.
        _require(re.fullmatch(r"[a-z0-9_-]+", rid) is not None, path,
                 f"manifest rule 5: reference_id {rid!r} is not "
                 f"label-safe (use lowercase letters, digits, '-', "
                 f"'_'); it is embedded verbatim in derived labels")
        _require(rid not in seen_ref_ids, path,
                 f"manifest rule 5: duplicate reference_id: {rid}")
        seen_ref_ids.add(rid)

        # ----- Rule 4: exactly one structure source, validated the
        # same way the full loader validates it.
        has_cod = "cod_id" in ref
        has_path = "structure_path" in ref
        _require(has_cod != has_path, path,
                 f"manifest rule 4: [[reference_solid {rid}]] "
                 f"must set exactly one of cod_id or "
                 f"structure_path")
        cod_id = None
        cod_revision = None
        structure_path = None
        if has_cod:
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
            cod_revision = ref["cod_revision"]
        else:
            # File existence is verified at materialize time, which
            #   resolves it against the manifest directory and reports
            #   the resolved path -- no need to duplicate that here.
            structure_path = ref["structure_path"]

        sources.append(ReferenceSolid(
            reference_id=rid,
            system_type=ref.get("system_type", ""), basis="",
            functional="", kpoint_integration="", kpoint_spec={},
            scf_threshold=0.0, cod_id=cod_id,
            cod_revision=cod_revision,
            structure_path=structure_path,
            source_description=ref.get("source_description")))

    return sources


# ============================================================
#  Manifest writer (DESIGN 5.7): serialize a CurationManifest
# ============================================================
#
# The hand-rolled emitter mirrors initial_potential_db.save in
# spirit -- the project writes its own TOML so the exact layout
# stays under our control -- but its goals differ.  The database
# writer's payload floats use a 17-digit form to round-trip
# binary64 exactly; a manifest is human-authored configuration
# instead, so floats use Python's shortest round-trippable repr
# (0.05, 5.0, 1e-06) and inline-table keys keep the order they
# were given, both so the file reads naturally to a curator.  The
# output round-trips through load_manifest_v2: writing a manifest
# and reading it back yields the same reference solids.


## TOML 1.0 basic-string escape table.  Backslash and double
## quote get their literal escapes; the named control characters
## get their short sequences.  Any other control character (below
## 0x20 or 0x7F) is emitted as "\uXXXX" by _quote.
_TOML_ESCAPES = {
    "\\": "\\\\",
    "\"": "\\\"",
    "\b": "\\b",
    "\t": "\\t",
    "\n": "\\n",
    "\f": "\\f",
    "\r": "\\r",
}


def _quote(text: str) -> str:
    """Quote a string as a TOML 1.0 basic string.

    Backslash, double quote, and the named control characters get
    the standard escape sequences; any other control character
    (code point below 0x20 or equal to 0x7F) gets a ``\\uXXXX``
    escape.  Every other character passes through unchanged.  This
    matters mainly for a free-text ``description``; the other
    string fields (reference_id, element symbols, method names) are
    already restricted to safe characters.
    """

    out = ["\""]
    for ch in text:
        if ch in _TOML_ESCAPES:
            out.append(_TOML_ESCAPES[ch])
        elif ord(ch) < 0x20 or ord(ch) == 0x7F:
            out.append(f"\\u{ord(ch):04X}")
        else:
            out.append(ch)
    out.append("\"")
    return "".join(out)


def _scalar(value: Any) -> str:
    """Render a scalar in TOML form for a manifest.

    ``bool`` is tested before ``int`` because Python's ``bool`` is
    a subclass of ``int``.  Integers print as bare digits; floats
    use ``repr``, which gives the shortest decimal string that
    reads back as the same float (so ``0.05`` stays ``0.05`` and
    ``1e-06`` stays ``1e-06``) rather than a padded scientific
    form; strings are quoted as TOML basic strings.
    """

    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, str):
        return _quote(value)
    raise TypeError(
        f"curation_manifest writer: unsupported value type "
        f"{type(value).__name__} (value: {value!r})")


def _value(value: Any) -> str:
    """Render any manifest value: a nested dict becomes an inline
    table, a list becomes an inline array, and anything else a
    scalar (:func:`_scalar`).  Inline tables and arrays are used
    for ``kpoint_spec``, a fingerprint ``sub_spec``, and a
    ``shift`` triple -- the only non-scalar manifest values."""

    if isinstance(value, dict):
        return _inline_table(value)
    if isinstance(value, list):
        return "[" + ", ".join(_value(item) for item in value) + "]"
    return _scalar(value)


def _inline_table(table: dict[str, Any]) -> str:
    """Render a dict as a TOML inline table, keeping key order.

    Unlike the database writer's alphabetical-sort emitter, key
    order is preserved here so a curator-authored ``sub_spec`` such
    as ``{ level = 2, thick = 0.5, cutoff = 5.0, tolerance = 0.05 }``
    reads in its natural order rather than alphabetised.  An empty
    table renders as ``{}``.
    """

    if not table:
        return "{}"
    pairs = [f"{key} = {_value(val)}" for key, val in table.items()]
    return "{ " + ", ".join(pairs) + " }"


def _field(key: str, value: Any) -> str:
    """Render one top-level ``key = value`` manifest line, quoting
    or formatting the value per its type (:func:`_value`)."""

    return f"{key} = {_value(value)}"


def format_manifest(manifest: CurationManifest) -> str:
    """Serialize a :class:`CurationManifest` to schema-v2 TOML text.

    The emitted file round-trips through :func:`load_manifest_v2`:
    the manifest it yields equals the one written here.  The
    database-wide ``[characterization]`` recipe is emitted first
    (after ``schema_version``), then each ``[[reference_solid]]``
    emits its scalar fields and inline ``kpoint_spec``, then its
    ``[[reference_solid.entry]]`` sub-tables, each followed by its
    ``[[reference_solid.entry.fingerprint]]`` sub-tables -- the
    order TOML requires (a table's own keys before any of its
    sub-tables).  Blank lines separate blocks for readability.

    Customization fields are emitted only when set, since every one
    is optional (DESIGN 5.2.2): ``element`` / ``atom_site`` /
    ``description`` / ``label`` are emitted only when not ``None``,
    and ``default`` only when true (an absent flag reads as false).
    A per-solid ``source_description`` is emitted only when present.
    The ``preferred`` flag is never serialized: it is recovered
    from the block a declaration lands in (``True`` for a
    ``[characterization]`` record, ``False`` for a per-entry one),
    so the writer emits no ``preferred`` key anywhere.
    """

    lines = [_field("schema_version", manifest.schema_version)]

    # The database-wide preferred recipe (DESIGN 5.7).  Each method
    #   is its own [[characterization.fingerprint]] sub-table; the
    #   preferred flag is structural, so it is never written.
    if manifest.characterization:
        lines.append("")
        lines.append("[characterization]")
        for fingerprint in manifest.characterization:
            lines.append("[[characterization.fingerprint]]")
            lines.append(_field("method", fingerprint.method))
            lines.append(_field("sub_spec", fingerprint.sub_spec))

    for solid in manifest.reference_solids:
        lines.append("")
        lines.append("[[reference_solid]]")
        lines.append(_field("reference_id", solid.reference_id))
        lines.append(_field("system_type", solid.system_type))
        # Exactly one structure source is set (rule 4); emit
        #   whichever this solid carries.
        if solid.cod_id is not None:
            lines.append(_field("cod_id", solid.cod_id))
            lines.append(_field("cod_revision", solid.cod_revision))
        else:
            lines.append(
                _field("structure_path", solid.structure_path))
        lines.append(_field("basis", solid.basis))
        lines.append(_field("functional", solid.functional))
        lines.append(
            _field("kpoint_integration", solid.kpoint_integration))
        lines.append(_field("kpoint_spec", solid.kpoint_spec))
        lines.append(_field("scf_threshold", solid.scf_threshold))
        if solid.source_description is not None:
            lines.append(_field("source_description",
                                solid.source_description))

        for entry in solid.entries:
            lines.append("")
            lines.append("[[reference_solid.entry]]")
            # Every customization field is optional (DESIGN 5.2.2),
            #   so each is emitted only when the entry sets it.
            if entry.element is not None:
                lines.append(_field("element", entry.element))
            if entry.atom_site is not None:
                lines.append(_field("atom_site", entry.atom_site))
            if entry.default:
                lines.append(_field("default", entry.default))
            if entry.description is not None:
                lines.append(
                    _field("description", entry.description))
            if entry.label is not None:
                lines.append(_field("label", entry.label))

            for fingerprint in entry.fingerprints:
                lines.append("")
                lines.append(
                    "[[reference_solid.entry.fingerprint]]")
                lines.append(_field("method", fingerprint.method))
                lines.append(
                    _field("sub_spec", fingerprint.sub_spec))

    return "\n".join(lines) + "\n"


def write_manifest(manifest: CurationManifest, path: str) -> None:
    """Write :func:`format_manifest` output to ``path`` (overwriting
    any existing file), the file-producing counterpart of the
    pure-text formatter."""

    with open(path, "w") as handle:
        handle.write(format_manifest(manifest))
