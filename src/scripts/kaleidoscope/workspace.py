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

from .model import KaleidoscopeError


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
    when value is None (so callers can omit absent fields)."""
    if value is None:
        return ""
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
    ``<root>/wingbeats/<id>[/<calc>]``.  The ``<calc>`` level exists
    only when the unit carries a calc tag (DESIGN 6.2.4)."""
    base = os.path.join(flight.root, "wingbeats", unit.id)
    return os.path.join(base, unit.calc) if unit.calc else base


def derive_calc_tag(unit):
    """Derive a default ``<calc>`` tag for a unit that gave none
    but shares its id with other units (DESIGN 6.2.4).  For the
    Imago wingbeat this is ``"<job>-<scf_basis>"`` read from the
    makeinput options, falling back to ``scf``/``fb`` when those
    keys are absent."""
    job = str(unit.options.get("job", "scf"))
    basis = str(unit.options.get("scf_basis", "fb"))
    return f"{job}-{basis}".lower()


def validate_flight(flight):
    """Enforce the id/``<calc>`` scheme at flight-build time
    (DESIGN 6.2.4) and abort on any violation, naming the
    offenders -- a silent fix would break the cache hit-test.

    Side effect: a unit that gave no ``calc`` but shares its id
    with another unit has a derived tag assigned in place, so
    the two no longer collide on one run directory."""
    id_counts = Counter(unit.id for unit in flight.units)
    seen = {}                       # id -> set of calc tags used
    for unit in flight.units:
        require_slug(unit.id, "id")
        if unit.calc is not None:
            require_slug(unit.calc, "calc")

        tag = unit.calc
        # One id with several units needs a distinguishing tag;
        #   derive one when the client did not supply it.
        if tag is None and id_counts[unit.id] > 1:
            tag = derive_calc_tag(unit)
            require_slug(tag, "derived calc")
            unit.calc = tag

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

def serialize_flight(flight):
    """Write ``<root>/flight.toml``: the authoritative record
    of *what was asked for*, kept separate from each run's
    ``status.toml`` record of *what happened* (DESIGN 6.2.1).

    This slice records the per-unit identity fields (id,
    structure, calc, wingbeat); the makeinput options and the
    cache key live with the run (cache_key.toml) and are not
    round-tripped here.  flight.toml is for inspection and as a
    resume record; a resume itself is just re-running the
    flight, whose cache hit-test skips the done units."""
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
            flight_file.write(toml_line("calc", unit.calc))
            flight_file.write(toml_line("wingbeat", unit.wingbeat))
