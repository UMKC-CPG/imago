"""kaleidoscope.cache -- the run-reuse cache
(DESIGN 6.2.5; PSEUDOCODE 13.4).

The mechanism is owned by kaleidoscope; the key *fields* are
supplied by the client on each ``CalcUnit`` (DESIGN 6.2.5).
Before launching a unit the dispatcher asks ``is_cache_hit``;
on a hit it skips the launch and reports the existing outcome,
so resuming a campaign is nothing more than re-running it.

The key has two parts, mirroring the producer's original
``is_cached_v2`` (DESIGN 5.7) and generalizing it:

- **scalar fields** -- written verbatim into ``cache_key.toml``
  and compared field-by-field;
- **key files** -- compared by *byte-comparison* against the
  copy already staged in the run directory.  No hashing: a
  developer can diff the files to see exactly why a cache
  missed, which a hash would hide.
"""

import filecmp
import os
import tomllib

from .workspace import emit_scalar, toml_line


def write_cache_key(run_dir, unit):
    """Write ``<run_dir>/cache_key.toml`` -- the identity
    snapshot taken when a unit is launched (DESIGN 6.2.5).
    Scalars are emitted under a ``[scalars]`` table in sorted
    order (so the file is deterministic); the key-file names are
    recorded as a ``files`` array for inspection.  The actual
    file *contents* are not copied here -- the comparison reads
    the files staged in the run directory by the run itself."""
    os.makedirs(run_dir, exist_ok=True)
    path = os.path.join(run_dir, "cache_key.toml")
    names = [key_file.name for key_file in unit.key_fields.files]
    with open(path, "w") as key_file:
        rendered = ", ".join(emit_scalar(n) for n in names)
        key_file.write(f"files = [{rendered}]\n")
        key_file.write("\n[scalars]\n")
        for name in sorted(unit.key_fields.scalars):
            key_file.write(
                toml_line(name, unit.key_fields.scalars[name])
            )


def _read_cache_key(run_dir):
    """Return the parsed ``cache_key.toml``, or None if absent."""
    path = os.path.join(run_dir, "cache_key.toml")
    if not os.path.exists(path):
        return None
    with open(path, "rb") as key_file:
        return tomllib.load(key_file)


def cache_key_matches(unit, run_dir):
    """True iff the unit's *current* key matches the snapshot in
    the run directory (DESIGN 6.2.5).  The scalars must match
    field-by-field, and every declared key file must byte-equal
    its staged copy under the run directory."""
    saved = _read_cache_key(run_dir)
    if saved is None:
        return False

    # Scalar fields: verbatim, field-by-field.
    if saved.get("scalars", {}) != unit.key_fields.scalars:
        return False

    # Key files: byte-compare each current source against the
    #   copy staged in the run directory (no hashing).
    for key_file in unit.key_fields.files:
        staged = os.path.join(run_dir, key_file.name)
        if not os.path.exists(staged):
            return False
        if not filecmp.cmp(key_file.source, staged,
                           shallow=False):
            return False
    return True


def is_cache_hit(unit, run_dir):
    """True iff this unit can be skipped: the run directory
    exists, its recorded key still matches, and its status is
    ``done`` (DESIGN 6.2.5).  Anything else -- no directory, a
    key mismatch, or a non-done status -- is a miss, so the unit
    is (re)launched."""
    from .workspace import read_status
    if not os.path.isdir(run_dir):
        return False
    status = read_status(run_dir)
    if status is None or status.get("status") != "done":
        return False
    return cache_key_matches(unit, run_dir)
