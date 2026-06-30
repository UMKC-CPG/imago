"""initial_potential_db.py -- Per-element potential database
file-format library.

Role in the pipeline
--------------------
This module is the *library* half of the library / producer /
consumer split documented in DESIGN 5.4.  It is a small, passive
helper that knows exactly one thing: how to read, validate, look
up entries in, and write per-element ``s_gaussian_pot.toml``
files that live at ``share/atomicPDB/<elem>/``.

It contains no orchestration, no SCF runs, and no curation
logic.  Its only external dependency is ``tomllib`` (Python
standard library).

Other scripts in the chain import this module:

* ``build_initial_potentials.py`` (the *producer*) calls
  :func:`save` to write augmented per-element database files
  harvested from converged Imago SCF runs.
* ``makeinput.py`` (the *consumer*) calls :func:`load` and
  :func:`lookup` at input-generation time to pick the right
  entry for each element and emit it into the Imago input file.
* ``bench_initial_potential.py`` reaches the library indirectly
  through Imago runs that flow through ``makeinput.py``.

Keeping the file format isolated in one module means that
read-only callers never pull in SCF-runner or manifest code,
that a future schema-version bump touches one file (per
ARCHITECTURE 8.7), and that unit tests can cover the file
format with synthetic byte strings instead of real Imago runs.

Schema version
--------------
This module implements **schema version 2** (DESIGN 5.2).
Relative to the Phase-1 v1 schema it adds two things that let
``makeinput.py`` drive fingerprint-based manifest selection
(DESIGN 5.6) without disturbing the on-the-wire Imago input
format:

* a required per-entry ``default`` boolean -- exactly one entry
  per file carries ``default = true``, the entry chosen when no
  scheme matches and no ``-pot LABEL`` override is given; and
* an optional ``[[potential.fingerprint]]`` inner array on every
  entry, describing the local environment of the reference atom
  site at the moment that entry's potential was harvested.

There is no on-disk compatibility with v1 files: :func:`load`
rejects any ``schema_version != 2``.  No v1 files exist in
production; the producer (``build_initial_potentials.py``)
regenerates every file as v2.

Public surface
--------------
FingerprintRecord
    Dataclass holding one ``[[potential.fingerprint]]`` record:
    a matcher ``method`` name, a method-specific ``sub_spec``
    inline table that fully qualifies the fingerprint, and a
    ``payload`` dict carrying the remaining record fields
    (shape is matcher-defined -- e.g., ``values`` for the
    bispectrum matcher, ``shell_code`` for the reduce matcher).

PotentialEntry
    Dataclass holding one labeled potential entry: ``label``,
    ``default``, ``description``, ``num_gaussians``,
    ``alpha_min``, ``alpha_max``, ``coefficients`` (the
    representative atom's harvested potential), ``alphas``,
    ``provenance``, and the (possibly empty) list of
    :class:`FingerprintRecord` objects attached to the entry.

ElementDatabase
    Dataclass holding the top-level per-element record:
    ``schema_version``, ``element_symbol``, ``nuclear_z``,
    ``nuclear_alpha``, ``covalent_radius``, plus the list of
    :class:`PotentialEntry` objects parsed from the file.

load(path, known_methods=None) -> ElementDatabase
    Read and validate a per-element TOML file.  Enforces the
    set of validation rules from DESIGN 5.2.  Raises ``ValueError``
    on any rule violation, with a message naming the file path,
    the entry label when applicable, and the field at fault.
    ``known_methods`` is the optional set of registered matcher
    names used to enforce rule 9 (a fingerprint ``method`` must
    be a registered matcher); callers without a registry
    (isolated unit tests) pass ``None`` and rule 9 is skipped.
    Rule 10 requires exactly one preferred fingerprint record per
    method present in the file (DESIGN 5.6.5).

require_provenance(prov, path, label)
    Per-entry provenance validator used by :func:`load`.  Made
    public because it is part of the spec surface (PSEUDOCODE
    11.1) and is useful to producers that build provenance
    blocks programmatically before assembling an entry.

lookup(db, label) -> PotentialEntry
    Return the entry whose ``label`` matches.  Raises
    ``KeyError`` on miss.

baseline(db) -> PotentialEntry
    Convenience for ``lookup(db, "isolated")``.  Guaranteed to
    succeed on any database returned by :func:`load` because
    rule 6 enforces the ``"isolated"`` baseline at load time.

default_entry(db) -> PotentialEntry
    Return the entry whose ``default`` flag is true; guaranteed
    to succeed and be unique by rule 7.  Used by ``makeinput.py``
    (DESIGN 5.6) on the no-scheme fallback branch.  Distinct
    from :func:`baseline`: the curator may mark the isolated
    baseline as default (the two then return the same object)
    or mark a curated improved entry as default (they then
    return different objects).

find_fingerprint(entry, method, sub_spec) -> FingerprintRecord
    Return the fingerprint record on ``entry`` whose
    ``(method, sub_spec)`` matches, using the same canonical
    sub-spec equality as rule 8; raises ``KeyError`` on miss.
    Used by the matcher dispatch in ``makeinput.py`` (DESIGN
    5.6, 8.9) and by the producer when deciding whether a
    fingerprint already exists.

save(db, path) -> None
    Write ``db`` to ``path`` using the deterministic hand-
    formatted emitter from DESIGN 5.5 and PSEUDOCODE 11.2.
    Byte-identical output for byte-identical input -- the
    emitter level of the pipeline's layered-reproducibility
    contract.  See DESIGN 5.5 for what the emitter does and
    does not guarantee.
"""

import os
import tomllib
from dataclasses import dataclass, field
from typing import Any


# ============================================================
#  Dataclasses (DESIGN 5.4)
# ============================================================

@dataclass
class FingerprintRecord:
    """One ``[[potential.fingerprint]]`` record on an entry.

    A fingerprint describes the local environment of the
    reference atom site at the moment the parent entry's
    numerical potential was harvested.  The matcher named by
    ``method`` (ARCHITECTURE 8.9) interprets the record:

    * ``method`` -- matcher name, currently ``"reduce"`` or
      ``"bispectrum"``.  The schema validates only that the
      method is registered (rule 9); the matcher owns the
      meaning of the rest of the record.
    * ``sub_spec`` -- a method-specific inline table that fully
      qualifies the fingerprint (e.g., ``{twoj1 = 8, twoj2 =
      8}`` for bispectrum).  Two records on the same entry with
      the same ``method`` but different ``sub_spec`` are
      non-comparable and may coexist (rule 8); two with the
      same ``(method, sub_spec)`` are a hard error.
    * ``preferred`` -- the schema-v2 fingerprint-selection hint
      (DESIGN 5.6.5 step 2).  When a run lets the input file
      dictate the environment scheme (the crystalline /
      pre-assigned regime, with no ``-reduce`` or ``-bispec`` on
      the command line), the consumer must still decide which of
      the entry's many fingerprint records to match against.  It
      does so by matching the database's *preferred* record for
      each present method.  Rule 10 guarantees exactly one record
      per present method carries ``preferred = true``, so the
      choice is unambiguous.  Manifest rule 11 further requires
      the preferred sub-spec to be uniform database-wide so that
      every element is matched at the same fidelity.  Records that
      are not preferred still participate in user-scheme matching
      (when the command line names a method and sub-spec
      explicitly); ``preferred`` only governs the file-dictated
      regime.
    * ``payload`` -- every record field other than ``method``,
      ``sub_spec``, and ``preferred``.  Its shape is
      matcher-defined and is *not* validated by the library
      (e.g., bispectrum records carry a ``values`` real array;
      reduce records carry a ``shell_code`` inline table).  The
      matcher validates its own payload at lookup time.
    """

    method: str
    sub_spec: dict[str, Any]
    preferred: bool = False
    payload: dict[str, Any] = field(default_factory=dict)


@dataclass
class PotentialEntry:
    """One labeled potential entry within a per-element file.

    Mirrors a single ``[[potential]]`` block in the TOML file.
    ``provenance`` is kept as a free-form dict because its
    contents depend on the entry's source: ``atomSCF`` entries
    carry only the base provenance keys, while ``Imago`` entries
    additionally carry the reference-run metadata documented in
    DESIGN 5.2.

    The ``default`` flag is the schema-v2 selection hint: exactly
    one entry per file carries ``default = true`` (rule 7), and
    that entry is the one chosen when no environment scheme
    matches and no ``-pot LABEL`` override is given.  It is
    independent of the ``"isolated"`` baseline (rule 6): the
    curator may mark either the baseline or a curated improved
    entry as the default.

    ``fingerprints`` is the (possibly empty) list of
    :class:`FingerprintRecord` objects parsed from the entry's
    ``[[potential.fingerprint]]`` sub-blocks.  An entry with no
    fingerprints participates only in literal-label and
    default-tag selection, never in environment-scheme matching.

    The ``coefficients`` and ``alphas`` arrays hold the entry's
    potential: the harvested potential of the representative atom
    for this environment, stored verbatim (DESIGN 5.2.3).  The
    database stores one representative per distinct environment;
    a statistical mean across every atom that maps to the
    environment -- with its per-coefficient spread and the atom
    and model counts -- is a deferred refinement (DESIGN 5.2.3
    "Deferred", TODO C103).
    """

    label: str
    default: bool
    description: str
    num_gaussians: int
    alpha_min: float
    alpha_max: float
    coefficients: list[float]
    alphas: list[float]
    provenance: dict[str, Any]
    fingerprints: list[FingerprintRecord] = field(
        default_factory=list)


@dataclass
class ElementDatabase:
    """The full contents of one ``s_gaussian_pot.toml`` file.

    Top-level fields are per-element constants (``nuclear_z``,
    ``nuclear_alpha``, ``covalent_radius``); ``potentials`` is
    the list of labeled entries the file carries.  Rule 6
    guarantees that any database returned by :func:`load`
    includes at least one entry whose ``label`` is
    ``"isolated"``.

    ``nuclear_z`` is the atomic number Z.  It is nominally an
    integer, but is stored and emitted as a real because Imago
    consumes Z as a real number (and the legacy ``pot1`` file
    already records it as one).  :func:`load` coerces it to
    ``float`` so the in-memory value is a real regardless of
    whether the on-disk file spells it ``79`` or ``7.9e+01``.
    """

    schema_version: int
    element_symbol: str
    nuclear_z: float
    nuclear_alpha: float
    covalent_radius: float
    potentials: list[PotentialEntry] = field(
        default_factory=list)


# ============================================================
#  Reader (DESIGN 5.2; PSEUDOCODE 11.1)
# ============================================================

def load(path: str,
         known_methods: set[str] | None = None
         ) -> ElementDatabase:
    """Read and validate one per-element TOML database file.

    Implements PSEUDOCODE 11.1 and enforces the set of validation
    rules from DESIGN 5.2 (schema v2):

    1. ``schema_version`` must equal ``2``.
    2. ``element_symbol`` matches the parent directory name
       (case-insensitive).
    3. Every required field is present, at both the top level
       and inside each ``[[potential]]`` block -- including the
       new per-entry ``default`` flag.
    4. ``len(coefficients) == len(alphas) == num_gaussians``
       for every entry.
    5. Labels are unique within the file.
    6. At least one entry carries ``label == "isolated"`` --
       the legacy baseline must always be present so that the
       validation harness (DESIGN 5.8) can compare against it
       on the same code path.  Independent of rule 7.
    7. Exactly one entry per file carries ``default = true``.
       Zero or multiple defaults is a hard error.  The producer
       guarantees this by tagging the curator-designated entry
       when the manifest declares one and otherwise falling back
       to the ``"isolated"`` baseline (DESIGN 5.7 rule 7), so a
       hand-free element still loads with exactly one default.
    8. Within any one entry's fingerprint array, the pair
       ``(method, sub_spec)`` is unique.  Two records with the
       same method and the same sub-spec keys-and-values are a
       hard error; the same method with a differing sub-spec is
       allowed.
    9. Each fingerprint ``method`` must be a registered matcher
       name.  This is checked only when ``known_methods`` is
       supplied -- the registry lives in ``makeinput.py``
       (ARCHITECTURE 8.9), so to keep this library free of any
       import from ``makeinput.py`` the caller passes the set of
       known names in.  Callers without a registry (isolated
       unit tests) pass ``None`` and rule 9 is skipped, which
       decouples the library from the registry without weakening
       the rule for real consumer and producer runs.
    10. For every fingerprint ``method`` *present* anywhere in
       the file (across all entries), exactly one record carries
       ``preferred = true``.  Zero preferred records for a
       present method, or two or more, is a hard error: the
       file-dictated selection regime (DESIGN 5.6.5) needs an
       unambiguous representative record per method, and the
       curator must declare it explicitly.  A method that is
       wholly absent from the file is not subject to the rule.

    Any violation raises ``ValueError`` whose message names the
    file path, the entry label when applicable, and the field
    at fault.
    """

    with open(path, "rb") as handle:
        raw = tomllib.load(handle)

    # ----- Rule 3 (top-level half): required keys present.
    # Check before any value-level rule so that the error
    # message names the missing field rather than failing
    # later with a value mismatch on an absent key.
    top_required = (
        "schema_version", "element_symbol", "nuclear_z",
        "nuclear_alpha", "covalent_radius",
    )
    for field_name in top_required:
        if field_name not in raw:
            raise ValueError(
                f"{path}: missing top-level field: "
                f"{field_name}")

    # ----- Rule 1: schema_version must equal 2.
    # There is no on-disk compatibility with v1 files (DESIGN
    # 5.2); the producer regenerates every file as v2, so any
    # other version is a hard error rather than a migration.
    if raw["schema_version"] != 2:
        raise ValueError(
            f"{path}: unsupported schema_version "
            f"{raw['schema_version']!r} (expected 2)")

    # ----- Rule 2: element_symbol must match parent dir.
    expected_elem = os.path.basename(
        os.path.dirname(os.path.abspath(path)))
    if raw["element_symbol"].lower() != expected_elem.lower():
        raise ValueError(
            f"{path}: element_symbol "
            f"{raw['element_symbol']!r} does not match "
            f"parent directory name {expected_elem!r}")

    db = ElementDatabase(
        schema_version  = raw["schema_version"],
        element_symbol  = raw["element_symbol"],
        # Coerce Z to a real: nominally integral, but Imago uses
        # it as a real number, so the in-memory value is always
        # a float regardless of how the file spells it.
        nuclear_z       = float(raw["nuclear_z"]),
        nuclear_alpha   = raw["nuclear_alpha"],
        covalent_radius = raw["covalent_radius"],
    )

    # ----- Per-entry validation
    seen_labels: set[str] = set()
    # Count of entries flagged default = true; checked against
    # rule 7 (exactly one) after the loop so a "zero" and a
    # "multiple" condition report through the same message.
    default_count = 0
    # Rule 10 bookkeeping: how many fingerprint records carry
    # preferred = true for each method seen anywhere in the file.
    # A method appears as a key here as soon as any record names
    # it (preferred or not), so after the loop a present method
    # whose count is not exactly one is a hard error.  Methods
    # wholly absent from the file never become keys and are thus
    # exempt from the rule.
    preferred_count_by_method: dict[str, int] = {}
    entry_required = (
        "default", "description", "num_gaussians", "alpha_min",
        "alpha_max", "coefficients", "alphas", "provenance",
    )

    for entry_dict in raw.get("potential", []):
        # Rule 3 (per-entry half): check "label" first so that
        # all subsequent error messages can name the entry by
        # label.  The "label" check itself cites only the file
        # path and the fact that the [[potential]] is anonymous.
        if "label" not in entry_dict:
            raise ValueError(
                f"{path}: [[potential]] entry is missing "
                f"required field: label")

        label = entry_dict["label"]
        for field_name in entry_required:
            if field_name not in entry_dict:
                raise ValueError(
                    f"{path}: [[potential]] '{label}' is "
                    f"missing required field: {field_name}")

        # Rule 4: length consistency between the declared
        # num_gaussians and the per-coefficient arrays.
        n_gauss = entry_dict["num_gaussians"]
        coeffs = entry_dict["coefficients"]
        alphas = entry_dict["alphas"]
        if len(coeffs) != n_gauss or len(alphas) != n_gauss:
            raise ValueError(
                f"{path}: [[potential]] '{label}': coefficients/"
                f"alphas length "
                f"({len(coeffs)}/{len(alphas)}) does not match "
                f"num_gaussians ({n_gauss})")

        # Rule 5: labels must be unique within the file.
        if label in seen_labels:
            raise ValueError(
                f"{path}: duplicate [[potential]] label: "
                f"{label!r}")
        seen_labels.add(label)

        # Rule 7 (tally half): accumulate the default flag now;
        # the file-level "exactly one" check runs after the
        # loop.  The flag is a required field (checked above),
        # so the lookup is safe here.
        if entry_dict["default"]:
            default_count += 1

        # Per-entry provenance validation (DESIGN 5.2 lower
        # half; PSEUDOCODE 11.1 require_provenance).
        require_provenance(
            entry_dict["provenance"], path, label)

        # ----- Fingerprint sub-blocks (rules 8 and 9).
        # Each [[potential.fingerprint]] record becomes one
        # FingerprintRecord.  We enforce per-entry
        # (method, sub_spec) uniqueness (rule 8) and -- when a
        # matcher registry was supplied -- that the method is
        # registered (rule 9).
        fingerprints: list[FingerprintRecord] = []
        seen_method_subspec: set = set()
        for fp_dict in entry_dict.get("fingerprint", []):
            if "method" not in fp_dict:
                raise ValueError(
                    f"{path}: [[potential]] '{label}' "
                    f"fingerprint is missing required field: "
                    f"method")
            if "sub_spec" not in fp_dict:
                raise ValueError(
                    f"{path}: [[potential]] '{label}' "
                    f"fingerprint is missing required field: "
                    f"sub_spec")
            method = fp_dict["method"]
            sub_spec = fp_dict["sub_spec"]

            # Rule 8: per-entry (method, sub_spec) uniqueness.
            # The sub-spec is canonicalized so that equal
            # tables with different key order or int-vs-float
            # spelling compare equal.
            key = (method, canonicalize_sub_spec(sub_spec))
            if key in seen_method_subspec:
                raise ValueError(
                    f"{path}: [[potential]] '{label}': "
                    f"duplicate fingerprint (method={method!r}, "
                    f"sub_spec={sub_spec!r})")
            seen_method_subspec.add(key)

            # Rule 9: method must be registered.  Skipped when
            # the caller passed no registry (test contexts).
            if (known_methods is not None
                    and method not in known_methods):
                raise ValueError(
                    f"{path}: [[potential]] '{label}': "
                    f"fingerprint method {method!r} is not a "
                    f"registered matcher")

            # Rule 10 bookkeeping: register the method (so it is
            # counted as present) and tally its preferred records.
            # The flag defaults to false when the curator omits
            # it, so the common non-preferred record needs no key.
            preferred = bool(fp_dict.get("preferred", False))
            preferred_count_by_method.setdefault(method, 0)
            if preferred:
                preferred_count_by_method[method] += 1

            # The payload is every field other than method,
            # sub_spec, and preferred; the owning matcher
            # validates its shape.  preferred is a library-level
            # selection hint, not part of the matcher payload, so
            # it is excluded here and carried on its own field.
            payload = {k: v for k, v in fp_dict.items()
                       if k not in ("method", "sub_spec",
                                    "preferred")}
            fingerprints.append(FingerprintRecord(
                method    = method,
                sub_spec  = dict(sub_spec),
                preferred = preferred,
                payload   = payload,
            ))

        db.potentials.append(PotentialEntry(
            label           = label,
            default         = entry_dict["default"],
            description     = entry_dict["description"],
            num_gaussians   = n_gauss,
            alpha_min       = entry_dict["alpha_min"],
            alpha_max       = entry_dict["alpha_max"],
            coefficients    = list(coeffs),
            alphas          = list(alphas),
            provenance      = dict(entry_dict["provenance"]),
            fingerprints    = fingerprints,
        ))

    # ----- Rule 6: the 'isolated' baseline must be present.
    # This guarantees that downstream callers of baseline() --
    # and in particular the validation harness in DESIGN 5.8 --
    # can always reach the atomSCF-derived reference entry on
    # exactly the same code path used for any other label.
    if "isolated" not in seen_labels:
        raise ValueError(
            f"{path}: missing required 'isolated' baseline "
            f"entry (every database must carry the "
            f"atomSCF-derived baseline so the validation "
            f"harness can compare against it on the same "
            f"code path)")

    # ----- Rule 7: exactly one default entry per file.  This is
    # a structural invariant of every produced file; the producer
    # guarantees it by falling back to the 'isolated' baseline
    # when the manifest designates no default (DESIGN 5.7 rule 7).
    # Zero or multiple defaults in a loaded file is still a hard
    # error -- the loader trusts the producer to have tagged one.
    if default_count != 1:
        raise ValueError(
            f"{path}: expected exactly one [[potential]] with "
            f"default = true; found {default_count}")

    # ----- Rule 10: exactly one preferred fingerprint record per
    # present method.  The file-dictated selection regime (DESIGN
    # 5.6.5) matches against the preferred record, so an absent or
    # ambiguous preference would leave the consumer with no
    # well-defined representative to compare against.  Sorting the
    # methods keeps the error message deterministic.
    for method in sorted(preferred_count_by_method):
        count = preferred_count_by_method[method]
        if count != 1:
            raise ValueError(
                f"{path}: fingerprint method {method!r} is "
                f"present but has {count} record(s) flagged "
                f"preferred = true; exactly one is required so "
                f"the file-dictated selection regime has an "
                f"unambiguous representative")

    return db


def require_provenance(prov: dict[str, Any], path: str,
                       label: str) -> None:
    """Validate the per-entry ``[potential.provenance]`` block.

    Every provenance dict must carry ``source``, ``commit``, and
    ``generated_at``.  ``source`` must be either ``"atomSCF"``
    or ``"Imago"``.  Imago-source entries additionally require
    the reference-run fields (``reference_id``, ``atom_site``,
    ``kpoint_spec``, ``scf_threshold``,
    ``scf_iterations``) so the validation harness in DESIGN 5.8
    can identify and re-run the originating SCF on demand.

    Raises ``ValueError`` with a message that names the file
    path, the entry label, and the missing or invalid field.
    """

    base_required = ("source", "commit", "generated_at")
    for field_name in base_required:
        if field_name not in prov:
            raise ValueError(
                f"{path}: [[potential]] '{label}' provenance "
                f"is missing required field: {field_name}")

    if prov["source"] not in ("atomSCF", "Imago"):
        raise ValueError(
            f"{path}: [[potential]] '{label}' "
            f"provenance.source must be 'atomSCF' or "
            f"'Imago', got {prov['source']!r}")

    if prov["source"] == "Imago":
        imago_required = (
            "reference_id", "atom_site", "kpoint_spec",
            "scf_threshold", "scf_iterations",
        )
        for field_name in imago_required:
            if field_name not in prov:
                raise ValueError(
                    f"{path}: [[potential]] '{label}' Imago "
                    f"provenance is missing required field: "
                    f"{field_name}")


def lookup(db: ElementDatabase, label: str) -> PotentialEntry:
    """Return the entry whose ``label`` matches; KeyError else.

    Linear scan; the number of entries per file is small (a
    handful of labeled potentials), so a hash map would not
    pay for itself.
    """

    for entry in db.potentials:
        if entry.label == label:
            return entry
    raise KeyError(
        f"no [[potential]] with label {label!r} in database "
        f"for element {db.element_symbol!r}")


def baseline(db: ElementDatabase) -> PotentialEntry:
    """Return the ``'isolated'`` baseline entry.

    Guaranteed to succeed on any database returned by
    :func:`load` because rule 6 enforces the baseline at load
    time.  Producers that build a database in memory must
    ensure they add the baseline entry before calling
    :func:`save` (the emitter does not enforce rule 6 itself --
    enforcement happens on the read side so that round-trips
    through :func:`load` catch the omission).
    """

    return lookup(db, "isolated")


def default_entry(db: ElementDatabase) -> PotentialEntry:
    """Return the entry whose ``default`` flag is true.

    Guaranteed to succeed and to be unique on any database
    returned by :func:`load`, because rule 7 enforces exactly
    one default at load time.  Used by ``makeinput.py``
    (DESIGN 5.6) on the no-scheme fallback branch -- the entry
    picked when no environment scheme matches and the user gave
    no ``-pot LABEL`` override.

    Distinct from :func:`baseline`: a curator may flag the
    ``"isolated"`` baseline as the default (then both functions
    return the same object) or flag a curated improved entry as
    the default (then the two functions return different
    objects).  The ``RuntimeError`` below should be unreachable
    on a loaded database; it guards against callers that build
    a database in memory and forget to set a default before
    handing it to a consumer.
    """

    for entry in db.potentials:
        if entry.default:
            return entry
    raise RuntimeError(
        f"no default entry in database for element "
        f"{db.element_symbol!r}; load() enforces rule 7, so "
        f"this database was built in memory without a default")


def find_fingerprint(entry: PotentialEntry, method: str,
                     sub_spec: dict[str, Any]
                     ) -> FingerprintRecord:
    """Return the matching fingerprint record on ``entry``.

    A record matches when its ``method`` is equal and its
    ``sub_spec`` is canonically equal to ``sub_spec`` -- the
    same key-and-value equality (key order irrelevant,
    int-vs-float spelling irrelevant) that rule 8 uses to police
    uniqueness.  Raises ``KeyError`` on miss; the caller decides
    whether a miss is fatal (the consumer preflight in DESIGN
    5.6) or expected (the producer checking whether a
    fingerprint already exists before harvesting one).
    """

    target = canonicalize_sub_spec(sub_spec)
    for fp in entry.fingerprints:
        if (fp.method == method
                and canonicalize_sub_spec(fp.sub_spec)
                    == target):
            return fp
    raise KeyError(
        f"no fingerprint with method={method!r}, "
        f"sub_spec={sub_spec!r} on entry {entry.label!r}")


def find_preferred(db: ElementDatabase, method: str
                   ) -> FingerprintRecord | None:
    """Return this database's preferred record for ``method``.

    Implements PSEUDOCODE 11.3.d ``find_preferred``.  In the
    file-dictated selection regime (DESIGN 5.6.5 step 2) the
    consumer matches the input structure against the database's
    single representative record for each candidate method rather
    than against a command-line-named sub-spec.  This returns that
    representative: the one fingerprint record (across all entries
    in the file) whose ``method`` matches and whose ``preferred``
    flag is set.

    Returns ``None`` when the method is wholly absent from the
    file, which the caller reads as "this database offers no
    fingerprints of that family, fall through to the next method
    or to the default entry".  Because rule 10 guarantees exactly
    one preferred record per *present* method, a ``None`` result
    can only mean "family absent", never "present but
    unpreferred", and the first matching record is necessarily the
    unique preferred one -- so the scan can stop at the first hit.
    """

    for entry in db.potentials:
        for fp in entry.fingerprints:
            if fp.method == method and fp.preferred:
                return fp
    return None


def canonicalize_sub_spec(sub_spec: Any) -> Any:
    """Return a hashable canonical form of a sub-spec value.

    Used both by :func:`load` (rule 8 uniqueness) and by
    :func:`find_fingerprint` (lookup equality), so that two
    sub-specs compare equal exactly when they carry the same
    keys and values, regardless of:

    * **key order** -- dicts are frozen into a ``frozenset`` of
      (key, canonical-value) pairs, which is order-independent;
    * **int-versus-float spelling** -- numeric values are
      normalized to ``float`` so that ``twoj1 = 8`` and
      ``twoj1 = 8.0`` are treated as the same parameter.

    Each value is tagged with its kind (``"num"``, ``"str"``,
    ``"bool"``, ``"dict"``, ``"list"``) so that values of
    different types never collide (e.g., the number ``1.0`` and
    the string ``"1.0"``).  ``bool`` is handled before the
    numeric branch because Python's ``bool`` is a subclass of
    ``int``; a boolean ``sub_spec`` flag must stay distinct from
    the integers 0 and 1.  The returned object is fully
    hashable, so it can live in a ``set`` (rule-8 dedup) or be
    compared with ``==`` (find_fingerprint).
    """

    if isinstance(sub_spec, bool):
        return ("bool", sub_spec)
    if isinstance(sub_spec, (int, float)):
        return ("num", float(sub_spec))
    if isinstance(sub_spec, str):
        return ("str", sub_spec)
    if isinstance(sub_spec, dict):
        return ("dict", frozenset(
            (k, canonicalize_sub_spec(v))
            for k, v in sub_spec.items()))
    if isinstance(sub_spec, (list, tuple)):
        return ("list", tuple(
            canonicalize_sub_spec(v) for v in sub_spec))
    raise TypeError(
        f"canonicalize_sub_spec: unsupported value type "
        f"{type(sub_spec).__name__} ({sub_spec!r})")


# ============================================================
#  Emitter (DESIGN 5.5; PSEUDOCODE 11.2)
# ============================================================

def save(db: ElementDatabase, path: str) -> None:
    """Write ``db`` to ``path`` via the deterministic emitter.

    Implements PSEUDOCODE 11.2 and the layout rules from
    DESIGN 5.5:

    * Fixed key ordering inside each block.
    * Floats use ``"%.16e"`` (16 digits after the decimal
      point, 17 significant digits in total) so that IEEE-754
      binary64 values round-trip exactly through the
      text form.
    * Within each block, ``=`` signs are aligned to one space
      past the longest key in that block; alignment width is a
      pure function of the block's keys, so the layout is
      idempotent against re-emission.
    * ``coefficients`` and ``alphas`` are emitted as multi-line
      arrays, one value per line indented by three spaces,
      with a trailing comma on every element.  TOML 1.0 allows
      trailing commas, and using them keeps diffs clean when
      arrays grow or shrink at the tail.
    * Exactly one blank line separates the top-level block
      from the first ``[[potential]]`` block, separates each
      entry's body from its ``[potential.provenance]``
      sub-block, and separates consecutive entries.  The file
      ends with exactly one trailing newline.

    Determinism is a *bit-level* guarantee on the emitter
    output: given the same :class:`ElementDatabase`, the
    emitter writes the same bytes.  It does **not** promise
    file-level byte-identity across pipeline runs (provenance
    timestamps refresh on each regeneration, and SCF / fit
    numerical drift can perturb the numbers themselves).  See
    DESIGN 5.5 for the full layered-reproducibility contract.
    """

    lines: list[str] = []

    # ----- Top-level block (fixed key order).
    top_keys = [
        "schema_version", "element_symbol", "nuclear_z",
        "nuclear_alpha", "covalent_radius",
    ]
    top_view = {
        "schema_version":  db.schema_version,
        "element_symbol":  db.element_symbol,
        "nuclear_z":       db.nuclear_z,
        "nuclear_alpha":   db.nuclear_alpha,
        "covalent_radius": db.covalent_radius,
    }
    lines.extend(_format_block(top_view, top_keys))
    lines.append("")

    # ----- One block per labeled entry.
    # `default` slots between `label` and `description`, matching
    # the DESIGN 5.3 gold sketch.
    body_keys = [
        "label", "default", "description", "num_gaussians",
        "alpha_min", "alpha_max",
    ]
    # The "=" alignment width spans the scalar body keys and the
    # array keys, so every array opener aligns with the rest of
    # the entry's body.
    align_keys = body_keys + ["coefficients", "alphas"]
    body_width = max(len(k) for k in align_keys)

    for entry in db.potentials:
        lines.append("[[potential]]")
        entry_view = {
            "label":         entry.label,
            "default":       entry.default,
            "description":   entry.description,
            "num_gaussians": entry.num_gaussians,
            "alpha_min":     entry.alpha_min,
            "alpha_max":     entry.alpha_max,
        }
        for key in body_keys:
            lines.append(_format_kv(
                key, entry_view[key], body_width))

        lines.append(_format_array_open(
            "coefficients", body_width))
        for value in entry.coefficients:
            lines.append("   " + _fmt_float(value) + ",")
        lines.append("]")

        lines.append(_format_array_open(
            "alphas", body_width))
        for value in entry.alphas:
            lines.append("   " + _fmt_float(value) + ",")
        lines.append("]")
        lines.append("")

        # ----- Provenance sub-block (canonical key order).
        lines.append("[potential.provenance]")
        prov_keys = _ordered_provenance_keys(entry.provenance)
        lines.extend(
            _format_block(entry.provenance, prov_keys))
        lines.append("")

        # ----- Fingerprint sub-blocks (v2).  Emitted in
        # insertion order, which matches the reader's parse
        # order; this keeps round-trip load(save(db))
        # byte-deterministic without forcing a method-and-
        # sub_spec ordering on the producer.  An empty
        # fingerprint list produces no blocks and no extra
        # blank lines.
        for fingerprint in entry.fingerprints:
            _emit_fingerprint_block(lines, fingerprint)

    # Trim trailing blank lines and finish with exactly one
    # newline so the file ends as POSIX text files conventionally
    # end (and so re-running the emitter is idempotent).
    while lines and lines[-1] == "":
        lines.pop()
    text = "\n".join(lines) + "\n"

    with open(path, "w", encoding="utf-8", newline="\n") as h:
        h.write(text)


def _format_block(values: dict[str, Any],
                  keys: list[str]) -> list[str]:
    """Emit one aligned key/value block in ``keys`` order.

    Width is the longest key in *this* block, so each block
    computes its own alignment locally -- the top-level keys,
    the entry body, and the provenance sub-block each get
    their own width.
    """

    width = max(len(k) for k in keys)
    lines = []
    for key in keys:
        value = values[key]
        # A dict value (an Imago provenance ``kpoint_spec``) renders
        #   as a single-line inline table, like a fingerprint
        #   sub_spec; every other value is a plain aligned scalar.
        if isinstance(value, dict):
            lines.append(
                f"{key.ljust(width)} = {_format_inline_table(value)}")
        else:
            lines.append(_format_kv(key, value, width))
    return lines


def _format_kv(key: str, value: Any, width: int) -> str:
    """Render one ``key = value`` line, key padded to width."""

    return f"{key.ljust(width)} = {_format_scalar(value)}"


def _format_array_open(key: str, width: int) -> str:
    """Render an array-opener line: ``key  = [``.

    Caller appends the values (indented by three spaces) and
    a closing ``]`` line of its own.
    """

    return f"{key.ljust(width)} = ["


def _emit_fingerprint_block(lines: list[str],
                            fp: FingerprintRecord) -> None:
    """Append one ``[[potential.fingerprint]]`` block.

    Implements PSEUDOCODE 11.2 ``emit_fingerprint_block``.  The
    block opens with ``method`` and ``sub_spec`` (in that fixed
    order), optionally followed by ``preferred = true``, then the
    payload fields in the payload dict's iteration order.  All
    ``=`` signs align within the block, so the alignment width
    spans the fixed keys plus every payload key.

    The ``preferred`` line is emitted *only* when the record is
    preferred (DESIGN 5.6.5 step 2 / rule 10); a non-preferred
    record omits the line entirely, matching the gold sketch in
    which ``preferred`` defaults to false and is never spelled
    out.  When present it sits immediately after ``sub_spec`` and
    before the payload, and it joins the alignment-width key set.

    Payload values are rendered by kind:

    * a list of floats becomes a multi-line array, matching the
      ``coefficients``/``alphas`` layout (opener, one indented
      value per line with a trailing comma, closing ``]``);
    * a dict becomes a single-line inline table (e.g., the
      reduce matcher's ``shell_code``);
    * any other scalar is rendered like any body scalar.

    ``sub_spec`` is always emitted as an inline table; unlike
    the pseudocode's ``format_kv`` brace-detection, this routes
    inline tables directly so that a payload *string* value
    that happens to start with ``{`` is never mistaken for a
    pre-formatted table.  The emitted bytes are identical.
    """

    lines.append("[[potential.fingerprint]]")

    # The ``preferred`` key participates in alignment only when
    # the record is preferred, since a non-preferred record omits
    # the line and must not pad the block to its width.
    fixed_keys = ["method", "sub_spec"]
    if fp.preferred:
        fixed_keys.append("preferred")
    payload_keys = list(fp.payload.keys())
    width = max(len(k) for k in fixed_keys + payload_keys)

    lines.append(
        f"{'method'.ljust(width)} = "
        f"{_format_scalar(fp.method)}")
    lines.append(
        f"{'sub_spec'.ljust(width)} = "
        f"{_format_inline_table(fp.sub_spec)}")
    # Emit the selection hint only when set (gold sketch never
    # spells out the false default).
    if fp.preferred:
        lines.append(
            f"{'preferred'.ljust(width)} = true")

    for key in payload_keys:
        value = fp.payload[key]
        if _is_float_list(value):
            lines.append(_format_array_open(key, width))
            for element in value:
                lines.append(
                    "   " + _fmt_float(element) + ",")
            lines.append("]")
        elif isinstance(value, dict):
            lines.append(
                f"{key.ljust(width)} = "
                f"{_format_inline_table(value)}")
        else:
            lines.append(_format_kv(key, value, width))

    lines.append("")


def _format_inline_table(table: dict[str, Any]) -> str:
    """Render a dict as a deterministic TOML inline table.

    Keys are emitted in sorted (alphabetical) order with one
    space inside the braces, `` = `` between key and value, and
    ``, `` between pairs -- so two logically equal tables with
    different insertion order produce identical bytes.  Nested
    dicts recurse into nested inline tables; all other values
    use :func:`_format_scalar`.  Empty tables render as
    ``{}``.
    """

    if not table:
        return "{}"
    parts = []
    for key in sorted(table.keys()):
        value = table[key]
        if isinstance(value, dict):
            rhs = _format_inline_table(value)
        elif isinstance(value, list):
            # An inline array (e.g. an Imago provenance
            #   kpoint_spec's ``shift`` triple).  Elements are
            #   scalars; nested arrays are not expected here.
            rhs = "[" + ", ".join(
                _format_scalar(element) for element in value) + "]"
        else:
            rhs = _format_scalar(value)
        parts.append(f"{key} = {rhs}")
    return "{ " + ", ".join(parts) + " }"


def _is_float_list(value: Any) -> bool:
    """True for a non-empty list whose elements are all floats.

    Used by the fingerprint emitter to decide between the
    multi-line array layout (a real-valued payload vector such
    as the bispectrum ``values``) and a single-line scalar or
    inline-table rendering.  An empty list does not qualify --
    there is nothing to lay out across lines -- and is handled
    by the caller's scalar branch instead.
    """

    return (isinstance(value, list)
            and len(value) > 0
            and all(isinstance(x, float) for x in value))


def _format_scalar(value: Any) -> str:
    """Render a scalar in TOML form per DESIGN 5.5.

    ``bool`` is handled first because Python's ``bool`` is a
    subclass of ``int`` and would otherwise be picked up by the
    integer branch.  ``int`` is rendered as bare digits;
    ``float`` is rendered via :func:`_fmt_float`; ``str`` is
    rendered via :func:`_toml_quote` (TOML basic string with
    standard escapes).
    """

    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return _fmt_float(value)
    if isinstance(value, str):
        return _toml_quote(value)
    raise TypeError(
        f"initial_potential_db.save: unsupported scalar "
        f"type {type(value).__name__} for emitter "
        f"(value: {value!r})")


def _fmt_float(value: float) -> str:
    """Format a float as ``"%.16e"`` (17 significant digits).

    Seventeen significant decimal digits are the minimum
    needed for any IEEE-754 binary64 value to round-trip
    exactly through a decimal string.  The constant-width
    output (always one digit before the decimal point, 16
    after, plus the sign and exponent) keeps the multi-line
    array layout visually aligned without any further
    padding logic.
    """

    return f"{value:.16e}"


def _ordered_provenance_keys(
        prov: dict[str, Any]) -> list[str]:
    """Return provenance keys in the canonical emit order.

    Base keys (``source``, ``commit``, ``generated_at``)
    always come first; the Imago-source extras follow only
    when ``prov["source"] == "Imago"``.  The order is fixed
    regardless of the dict's iteration order, which keeps the
    output deterministic across Python versions and across
    in-memory constructions of the same logical record.
    """

    base = ["source", "commit", "generated_at"]
    extras = [
        "reference_id", "atom_site", "kpoint_spec",
        "scf_threshold", "scf_iterations",
    ]
    if prov.get("source") == "Imago":
        ordered = base + extras
        # system_type is an optional forensic extra the producer
        #   records (DESIGN 5.7 rule 2); emitted only when present
        #   so a hand-written Imago entry without it still saves.
        if "system_type" in prov:
            ordered = ordered + ["system_type"]
        return ordered
    return base


## TOML 1.0 basic-string escape table.  Backslash, double
## quote, and the named control-character escapes get their
## specific sequences; any other control character (U+0000
## through U+001F and U+007F per TOML 1.0) is emitted as
## "\uXXXX" by the caller.  Order matters: backslash must be
## first in the dict so that an iteration that wanted to
## re-escape an already-escaped string would not double-escape
## (the emitter never does this, but the convention guards
## against future misuse).
_TOML_ESCAPES = {
    "\\": "\\\\",
    "\"": "\\\"",
    "\b": "\\b",
    "\t": "\\t",
    "\n": "\\n",
    "\f": "\\f",
    "\r": "\\r",
}


def _toml_quote(value: str) -> str:
    """Quote a string per TOML 1.0 basic-string rules.

    Backslash, double quote, and the named control characters
    get the standard escape sequences from
    :data:`_TOML_ESCAPES`.  Other control characters (any
    code point below ``0x20`` or equal to ``0x7F``) get a
    ``\\uXXXX`` escape with uppercase hex digits for visual
    consistency with the TOML spec examples.  All other
    characters (printable ASCII and all non-control Unicode)
    pass through unchanged.
    """

    out = ["\""]
    for ch in value:
        if ch in _TOML_ESCAPES:
            out.append(_TOML_ESCAPES[ch])
        elif ord(ch) < 0x20 or ord(ch) == 0x7F:
            out.append(f"\\u{ord(ch):04X}")
        else:
            out.append(ch)
    out.append("\"")
    return "".join(out)
