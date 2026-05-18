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

Public surface
--------------
PotentialEntry
    Dataclass holding one labeled potential entry: ``label``,
    ``description``, ``num_gaussians``, ``alpha_min``,
    ``alpha_max``, ``coefficients``, ``alphas``, ``provenance``.

ElementDatabase
    Dataclass holding the top-level per-element record:
    ``schema_version``, ``element_symbol``, ``nuclear_z``,
    ``nuclear_alpha``, ``covalent_radius``, plus the list of
    :class:`PotentialEntry` objects parsed from the file.

load(path) -> ElementDatabase
    Read and validate a per-element TOML file.  Enforces the
    six validation rules from DESIGN 5.2.  Raises ``ValueError``
    on any rule violation, with a message naming the file path,
    the entry label when applicable, and the field at fault.

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
class PotentialEntry:
    """One labeled potential entry within a per-element file.

    Mirrors a single ``[[potential]]`` block in the TOML file.
    ``provenance`` is kept as a free-form dict because its
    contents depend on the entry's source: ``atomSCF`` entries
    carry only the base provenance keys, while ``Imago`` entries
    additionally carry the reference-run metadata documented in
    DESIGN 5.2.
    """

    label: str
    description: str
    num_gaussians: int
    alpha_min: float
    alpha_max: float
    coefficients: list[float]
    alphas: list[float]
    provenance: dict[str, Any]


@dataclass
class ElementDatabase:
    """The full contents of one ``s_gaussian_pot.toml`` file.

    Top-level fields are per-element constants (``nuclear_z``,
    ``nuclear_alpha``, ``covalent_radius``); ``potentials`` is
    the list of labeled entries the file carries.  Rule 6
    guarantees that any database returned by :func:`load`
    includes at least one entry whose ``label`` is
    ``"isolated"``.
    """

    schema_version: int
    element_symbol: str
    nuclear_z: int
    nuclear_alpha: float
    covalent_radius: float
    potentials: list[PotentialEntry] = field(
        default_factory=list)


# ============================================================
#  Reader (DESIGN 5.2; PSEUDOCODE 11.1)
# ============================================================

def load(path: str) -> ElementDatabase:
    """Read and validate one per-element TOML database file.

    Implements PSEUDOCODE 11.1 and enforces the six validation
    rules from DESIGN 5.2:

    1. ``schema_version`` must equal ``1``.
    2. ``element_symbol`` matches the parent directory name
       (case-insensitive).
    3. Every required field is present, at both the top level
       and inside each ``[[potential]]`` block.
    4. ``len(coefficients) == len(alphas) == num_gaussians``
       for every entry.
    5. Labels are unique within the file.
    6. At least one entry carries ``label == "isolated"`` --
       the legacy baseline must always be present so that the
       validation harness (DESIGN 5.8) can compare against it
       on the same code path.

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

    # ----- Rule 1: schema_version must equal 1.
    if raw["schema_version"] != 1:
        raise ValueError(
            f"{path}: unsupported schema_version "
            f"{raw['schema_version']!r} (expected 1)")

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
        nuclear_z       = raw["nuclear_z"],
        nuclear_alpha   = raw["nuclear_alpha"],
        covalent_radius = raw["covalent_radius"],
    )

    # ----- Per-entry validation
    seen_labels: set[str] = set()
    entry_required = (
        "description", "num_gaussians", "alpha_min",
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
        # num_gaussians and the two arrays.
        n_gauss = entry_dict["num_gaussians"]
        coeffs = entry_dict["coefficients"]
        alphas = entry_dict["alphas"]
        if len(coeffs) != n_gauss or len(alphas) != n_gauss:
            raise ValueError(
                f"{path}: [[potential]] '{label}': "
                f"coefficients/alphas length "
                f"({len(coeffs)}/{len(alphas)}) does not "
                f"match num_gaussians ({n_gauss})")

        # Rule 5: labels must be unique within the file.
        if label in seen_labels:
            raise ValueError(
                f"{path}: duplicate [[potential]] label: "
                f"{label!r}")
        seen_labels.add(label)

        # Per-entry provenance validation (DESIGN 5.2 lower
        # half; PSEUDOCODE 11.1 require_provenance).
        require_provenance(
            entry_dict["provenance"], path, label)

        db.potentials.append(PotentialEntry(
            label         = label,
            description   = entry_dict["description"],
            num_gaussians = n_gauss,
            alpha_min     = entry_dict["alpha_min"],
            alpha_max     = entry_dict["alpha_max"],
            coefficients  = list(coeffs),
            alphas        = list(alphas),
            provenance    = dict(entry_dict["provenance"]),
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

    return db


def require_provenance(prov: dict[str, Any], path: str,
                       label: str) -> None:
    """Validate the per-entry ``[potential.provenance]`` block.

    Every provenance dict must carry ``source``, ``commit``, and
    ``generated_at``.  ``source`` must be either ``"atomSCF"``
    or ``"Imago"``.  Imago-source entries additionally require
    the reference-run fields (``reference_id``, ``atom_site``,
    ``kpoint_spec``, ``convergence_threshold``,
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
            "convergence_threshold", "scf_iterations",
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
    body_keys = [
        "label", "description", "num_gaussians",
        "alpha_min", "alpha_max",
    ]
    # The "=" alignment width spans both the scalar body keys
    # and the two array keys, so that the array openers align
    # with the rest of the entry's body.
    align_keys = body_keys + ["coefficients", "alphas"]
    body_width = max(len(k) for k in align_keys)

    for entry in db.potentials:
        lines.append("[[potential]]")
        entry_view = {
            "label":         entry.label,
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
    return [_format_kv(k, values[k], width) for k in keys]


def _format_kv(key: str, value: Any, width: int) -> str:
    """Render one ``key = value`` line, key padded to width."""

    return f"{key.ljust(width)} = {_format_scalar(value)}"


def _format_array_open(key: str, width: int) -> str:
    """Render an array-opener line: ``key  = [``.

    Caller appends the values (indented by three spaces) and
    a closing ``]`` line of its own.
    """

    return f"{key.ljust(width)} = ["


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
        "convergence_threshold", "scf_iterations",
    ]
    if prov.get("source") == "Imago":
        return base + extras
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
