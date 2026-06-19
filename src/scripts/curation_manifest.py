"""Curation-manifest schema library (DESIGN 5.7).

The curation manifest is the hand-authored file that drives the
initial-potential database build: it names each reference solid (its
structure source and SCF run settings) and the per-site potentials to
harvest from the converged calculation.  This module is the single,
neutral home for that schema -- the dataclasses a parsed manifest
becomes and the readers that parse and validate the file.

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
    """One ``[[reference_solid.entry.fingerprint]]`` declaration.

    Tells the producer to compute and harvest one fingerprint
    record alongside the numerical potential (the harvest itself
    is C60).  ``method`` names a matcher (ARCHITECTURE 8.9) and
    ``sub_spec`` is its method-specific parameter table.  Unlike a
    :class:`initial_potential_db.FingerprintRecord`, a declaration
    carries *no* payload: the payload is the SCF/loen output the
    producer will harvest, not something the curator writes.

    ``preferred`` marks the one record per matcher family the
    consumer matches against for a file-dictated (crystalline)
    structure (DESIGN 5.6.5 step 2).  Exactly one declaration per
    ``method`` is preferred for each element that carries any
    fingerprint of that family (manifest rule 10), and every
    preferred declaration of a family shares one ``sub_spec``
    across the whole manifest (manifest rule 11).  The flag rides
    through harvest onto the produced
    :class:`initial_potential_db.FingerprintRecord`.
    """

    method: str
    sub_spec: dict[str, Any]
    preferred: bool = False


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
    C54.  Rule 8 (per-entry ``(method, sub_spec)`` uniqueness) is
    enforced regardless, using the same canonical sub-spec equality
    as the per-element-database reader.

    The validation rules (DESIGN 5.7):

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
    10. For every ``(element, method)`` that appears with any
       fingerprint declaration across the manifest, exactly one
       declaration carries ``preferred = true``.  Zero or
       two-or-more is a hard error: the file-dictated consumer
       (DESIGN 5.6.5) needs one unambiguous record per family.
    11. All ``preferred = true`` declarations sharing a ``method``
       carry the identical ``sub_spec`` across the whole manifest,
       so the consumer queries one sub_spec per family (and a
       multi-element structure costs a single loen run).
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
    # Fingerprint preferred-flag bookkeeping (rules 10 and 11,
    # post-loop).  (element, method) -> count of declarations and
    # count of preferred declarations: every present (element,
    # method) needs exactly one preferred (rule 10).  method ->
    # the canonical sub_spec of the family's preferred record:
    # every preferred declaration of one method must agree (rule
    # 11), so the consumer queries one sub_spec per family.
    fingerprint_decls_seen: set[tuple[str, str]] = set()
    preferred_per_elem_method: dict[tuple[str, str], int] = {}
    preferred_subspec_per_method: dict[str, Any] = {}

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

            # ----- Fingerprint declarations (rules 8, 9, and the
            # tallies for rules 10 and 11, checked post-loop).
            fingerprints: list[ManifestFingerprint] = []
            seen_method_subspec: set = set()
            for fp in entry.get("fingerprint", []):
                _require("method" in fp and "sub_spec" in fp, path,
                         f"manifest rule 8: fingerprint "
                         f"declaration must carry both method and "
                         f"sub_spec ({rid}, label={label})")
                method = fp["method"]
                sub_spec = fp["sub_spec"]
                preferred = bool(fp.get("preferred", False))

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

                # Tally for rules 10 and 11 (checked post-loop).
                # Every (element, method) seen here is "present", so
                # it must end with exactly one preferred record; the
                # family's preferred sub_spec must be uniform.
                elem_method = (elem, method)
                fingerprint_decls_seen.add(elem_method)
                if preferred:
                    preferred_per_elem_method[elem_method] = (
                        preferred_per_elem_method.get(
                            elem_method, 0) + 1)
                    # Rule 11: a preferred record fixes the family's
                    # database-wide sub_spec; a later preferred
                    # record of the same method must match it.
                    prior = preferred_subspec_per_method.get(method)
                    if prior is None:
                        preferred_subspec_per_method[method] = canon
                    else:
                        _require(canon == prior, path,
                                 f"manifest rule 11: preferred "
                                 f"{method!r} sub_spec {sub_spec!r} "
                                 f"({rid}, label={label}) diverges "
                                 f"from the database-wide preferred "
                                 f"sub_spec for {method!r}; all "
                                 f"preferred records of a family "
                                 f"must share one sub_spec")

                fingerprints.append(ManifestFingerprint(
                    method=method, sub_spec=dict(sub_spec),
                    preferred=preferred))

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

    # ----- Rule 10 post-loop: every (element, method) present in
    # any fingerprint declaration carries exactly one preferred
    # record across the manifest.  Sorting keeps the message
    # deterministic.  (Rule 11's uniform-sub_spec check ran inline
    # above, as each preferred record was seen.)
    for elem, method in sorted(fingerprint_decls_seen):
        count = preferred_per_elem_method.get((elem, method), 0)
        _require(count == 1, path,
                 f"manifest rule 10: element {elem} method "
                 f"{method!r} has {count} preferred fingerprint "
                 f"declarations across the manifest (need exactly "
                 f"one so the file-dictated consumer has an "
                 f"unambiguous record to match)")

    return CurationManifest(
        schema_version=raw["schema_version"],
        manifest_path=os.path.abspath(path),
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
    (rule 4).  The run and harvest fields are left at harmless
    placeholders, because this relaxed view is never dispatched, only
    materialized; that lets a half-written manifest still be
    pre-flighted."""

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
            reference_id=rid, system_type="", basis="",
            functional="", kpoint_integration="", kpoint_spec={},
            scf_threshold=0.0, cod_id=cod_id,
            cod_revision=cod_revision,
            structure_path=structure_path))

    return sources
