"""kaleidoscope.workspace -- run-directory layout, id rules,
and the status.toml lifecycle (DESIGN 6.2.4; PSEUDOCODE 13.3).

This module pins the workspace scheme that ARCHITECTURE 9.8 had
left open: the stable-id convention (a filesystem-safe slug),
the optional ``<calc>`` tag for one structure hosting several
calculations, and the ``status.toml`` schema.  It also provides
the small, deterministic TOML emitters the cache and the wingbeat
reuse (``emit_scalar`` / ``toml_line``), kept here so the whole
package writes TOML one way.
"""

import os
import re
import tomllib
from collections import Counter

from .model import CalcUnit, Flight, KaleidoscopeError, SweepRecord


# A run-directory id (and any <calc> tag) must be a slug:
#   lowercase letters, digits, underscore, and hyphen.  This
#   keeps ids safe as directory names and, by forbidding a
#   silent rewrite, keeps the cache hit-test honest (a rewritten
#   id would resolve to a different directory than the one a
#   prior run created).
_SLUG_RE = re.compile(r"^[a-z0-9_-]+$")


# ------------------------------------------------------------------
#  Minimal, deterministic TOML emitters (shared across the package)
# ------------------------------------------------------------------

def emit_scalar(value):
    """Render a Python scalar as a TOML value.  Supports the
    only types kaleidoscope writes: bool, int, float, and str
    (strings are quoted with backslash/quote escaped)."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    text = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{text}"'


def toml_line(key, value):
    """Render one ``key = value`` TOML line, or the empty string
    when value is None (so callers can omit absent fields).

    A list or tuple value renders as a TOML array
    (``key = [a, b, ...]``), each element passed through
    ``emit_scalar``; an empty sequence renders as ``key = []``.
    This is what lets a unit's ``calc`` tuple round-trip through
    both ``status.toml`` and ``flight.toml`` as an array."""
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        rendered = ", ".join(emit_scalar(item) for item in value)
        return f"{key} = [{rendered}]\n"
    return f"{key} = {emit_scalar(value)}\n"


# ------------------------------------------------------------------
#  Identifiers and run-directory paths
# ------------------------------------------------------------------

def require_slug(value, what):
    """Raise unless ``value`` is a slug; ``what`` names the field
    for the error message (DESIGN 6.2.4)."""
    if not isinstance(value, str) or not _SLUG_RE.match(value):
        raise KaleidoscopeError(
            f"{what} is not a slug ([a-z0-9_-]+): {value!r}"
        )


def unit_run_dir(flight, unit):
    """Absolute path of a unit's run directory:
    ``<root>/wingbeats/<id>[/<calc components...>]``.  Each element
    of the unit's ``calc`` tuple becomes one nested directory
    level (DESIGN 6.2.1/6.2.4), so a one-axis sweep adds a single
    ``<axis>-<value>`` level and an empty tuple adds none -- the
    common single-calculation case resolves to
    ``<root>/wingbeats/<id>``."""
    return os.path.join(
        flight.root, "wingbeats", unit.id, *unit.calc)


def derive_calc_tag(unit):
    """Derive a default ``<calc>`` tag for a unit that gave none
    but shares its id with other units (DESIGN 6.2.4).  For the
    Imago wingbeat this is ``"<job>-<scf_basis>"`` read from the
    makeinput options, falling back to ``scf``/``fb`` when those
    keys are absent.  Returned as a one-element tuple to match the
    ``CalcUnit.calc`` tuple shape (one component per axis); the
    derived tag is the single distinguishing directory level."""
    job = str(unit.options.get("job", "scf"))
    basis = str(unit.options.get("scf_basis", "fb"))
    return (f"{job}-{basis}".lower(),)


def validate_flight(flight):
    """Enforce the id/``<calc>`` scheme at flight-build time
    (DESIGN 6.2.4) and abort on any violation, naming the
    offenders -- a silent fix would break the cache hit-test.

    Side effect: a unit that gave no ``calc`` but shares its id
    with another unit has a derived tag assigned in place, so
    the two no longer collide on one run directory."""
    id_counts = Counter(unit.id for unit in flight.units)
    seen = {}                       # id -> set of calc tuples used
    for unit in flight.units:
        require_slug(unit.id, "id")
        # Every component of the calc tuple is its own directory
        #   level, so each must independently be a slug.
        for component in unit.calc:
            require_slug(component, "calc")

        tag = unit.calc
        # One id with several units needs a distinguishing tag;
        #   derive one (a one-element tuple) when the client gave
        #   no calc at all (the empty tuple).
        if tag == () and id_counts[unit.id] > 1:
            tag = derive_calc_tag(unit)
            for component in tag:
                require_slug(component, "derived calc")
            unit.calc = tag

        # The full calc tuple is the per-id collision key: two
        #   units sharing an id must differ in their calc tuple or
        #   they would resolve to the same run directory.
        used = seen.setdefault(unit.id, set())
        if tag in used:
            raise KaleidoscopeError(
                f"duplicate run directory for id={unit.id!r} "
                f"calc={tag!r}"
            )
        used.add(tag)


# ------------------------------------------------------------------
#  status.toml -- the per-run lifecycle record
# ------------------------------------------------------------------

def write_status(wingbeat_dir, **fields):
    """Update ``<wingbeat_dir>/status.toml`` with the given fields.

    The write *merges* into any existing file so that fields
    accumulated earlier in the lifecycle (``submitted_at`` at
    queue time, ``started_at`` when a worker picks the unit up)
    survive a later terminal write.  None-valued fields are
    skipped.  The five lifecycle ``status`` values are
    kaleidoscope-owned and generic (queued / running / done /
    failed / lost); convergence is never a status -- it rides in
    the wingbeat-supplied ``detail`` (DESIGN 6.2.4)."""
    os.makedirs(wingbeat_dir, exist_ok=True)
    merged = read_status(wingbeat_dir) or {}
    for key, value in fields.items():
        if value is not None:
            merged[key] = value
    path = os.path.join(wingbeat_dir, "status.toml")
    with open(path, "w") as status_file:
        for key, value in merged.items():
            status_file.write(toml_line(key, value))


def read_status(wingbeat_dir):
    """Return the parsed ``status.toml`` for a run directory, or
    None when the file does not exist."""
    path = os.path.join(wingbeat_dir, "status.toml")
    if not os.path.exists(path):
        return None
    with open(path, "rb") as status_file:
        return tomllib.load(status_file)


# ------------------------------------------------------------------
#  flight.toml -- the inspectable record of what was asked
# ------------------------------------------------------------------

def _serialize_table(flight_file, header, table):
    """Write a TOML table ``[header]`` from a plain dict of
    scalars, arrays, and (recursively) nested dicts.  All scalar
    and array lines are emitted first, then each nested-dict value
    as a ``[header.<subkey>]`` sub-table -- the order TOML
    requires, since within a table every bare ``key = value`` must
    precede any sub-table header.

    This emits the optional ``[flight.sweep]`` block and each
    opaque ``[flight.<key>]`` metadata table.  Metadata tables are
    written VERBATIM: the dispatch core never inspects their
    contents (Principle 9), it only round-trips them so a
    domain-aware client (the 6.2.8 helper's PredictionRecord) can
    read them back after the flight."""
    flight_file.write(f"\n[{header}]\n")
    nested = []
    for key, value in table.items():
        if isinstance(value, dict):
            nested.append((key, value))
        else:
            flight_file.write(toml_line(key, value))
    for key, value in nested:
        _serialize_table(flight_file, f"{header}.{key}", value)


def serialize_flight(flight):
    """Write ``<root>/flight.toml``: the authoritative record
    of *what was asked for*, kept separate from each run's
    ``status.toml`` record of *what happened* (DESIGN 6.2.1).

    This slice records the per-unit identity fields (id,
    structure, calc, wingbeat); the makeinput options and the
    cache key live with the run (cache_key.toml) and are not
    round-tripped here.  flight.toml is for inspection and as a
    resume record; a resume itself is just re-running the
    flight, whose cache hit-test skips the done units.

    When the flight carries a ``sweep`` (it was built by the
    predict-then-verify helper, DESIGN 6.2.8) a ``[flight.sweep]``
    block records the varied and fixed axes, and every
    ``metadata[key]`` is emitted as a verbatim ``[flight.<key>]``
    table (e.g. ``[flight.prediction]``) so the harvest step can
    recover it without path-parsing."""
    os.makedirs(flight.root, exist_ok=True)
    path = os.path.join(flight.root, "flight.toml")
    with open(path, "w") as flight_file:
        flight_file.write(toml_line("root", flight.root))
        flight_file.write(
            toml_line("default_wingbeat", flight.default_wingbeat)
        )
        for unit in flight.units:
            flight_file.write("\n[[unit]]\n")
            flight_file.write(toml_line("id", unit.id))
            flight_file.write(
                toml_line("structure", unit.structure)
            )
            # calc renders as a TOML array of directory components
            #   (``calc = []`` for a unit with no second level).
            flight_file.write(toml_line("calc", unit.calc))
            flight_file.write(toml_line("wingbeat", unit.wingbeat))
            # kind is the run-role label (DESIGN 6.2.9); written
            #   verbatim so each harvester can filter on it.
            flight_file.write(toml_line("kind", unit.kind))

        # The optional sweep description (DESIGN 6.2.1/6.2.8): its
        #   fixed_axes dict nests as a [flight.sweep.fixed_axes]
        #   sub-table.
        if flight.sweep is not None:
            _serialize_table(flight_file, "flight.sweep", {
                "varied_axes": flight.sweep.varied_axes,
                "fixed_axes": flight.sweep.fixed_axes,
            })

        # Opaque metadata tables, each emitted verbatim (the
        #   6.2.8 helper stashes its prediction under "prediction").
        for key, table in flight.metadata.items():
            _serialize_table(flight_file, f"flight.{key}", table)


def read_flight_toml(path):
    """Parse ``flight.toml`` back into a :class:`Flight` -- the
    disk-side inverse of :func:`serialize_flight` (DESIGN 6.2.1).

    The harvest step (DESIGN 7.8) is a separate process from the
    one that dispatched the flight, so it recovers the plan from
    disk rather than from an in-memory object.  Only the fields
    ``serialize_flight`` persists are restored: each unit's
    identity (id, structure, calc tuple, wingbeat, kind) -- but NOT
    its makeinput ``options`` or ``key_fields``, which live with
    the run, not the plan -- plus the optional ``[flight.sweep]``
    block and every opaque ``[flight.<key>]`` metadata table
    (e.g. ``[flight.prediction]``).  The sweep's ordered
    ``varied_axes`` is what lets the harvest read each swept value
    out of the calc tag without guessing the axis name.
    """
    with open(path, "rb") as flight_file:
        data = tomllib.load(flight_file)

    units = []
    for raw_unit in data.get("unit", []):
        units.append(CalcUnit(
            id=raw_unit["id"],
            structure=raw_unit["structure"],
            # calc is stored as a TOML array; restore the tuple
            #   shape so unit_run_dir can splat it (empty == no
            #   second directory level).
            calc=tuple(raw_unit.get("calc", [])),
            wingbeat=raw_unit.get("wingbeat"),
            # kind defaults to "convergence" for a pre-kind
            #   flight.toml (DESIGN 6.2.9).
            kind=raw_unit.get("kind", "convergence")))

    # Everything under the [flight] table: the optional sweep
    #   description plus the opaque metadata tables.  The sweep is
    #   pulled out into a SweepRecord; the rest stay as plain dicts
    #   in metadata, exactly as the dispatch core round-trips them.
    flight_block = data.get("flight", {})
    sweep = None
    if "sweep" in flight_block:
        raw_sweep = flight_block["sweep"]
        sweep = SweepRecord(
            varied_axes=tuple(raw_sweep.get("varied_axes", [])),
            fixed_axes=dict(raw_sweep.get("fixed_axes", {})))
    metadata = {key: value for key, value in flight_block.items()
                if key != "sweep"}

    return Flight(
        root=data["root"],
        units=units,
        default_wingbeat=data.get("default_wingbeat", "imago"),
        sweep=sweep,
        metadata=metadata)


def flight_id_of(workspace_root):
    """The flight's identifier for provenance (DESIGN 7.8): the
    basename of its workspace root directory.  A flight dispatched
    into ``.../flights/diamond-seed`` has flight_id
    ``"diamond-seed"``.  Kept as a one-liner here so the harvest
    and any future consumer agree on the convention."""
    return os.path.basename(os.path.normpath(workspace_root))
