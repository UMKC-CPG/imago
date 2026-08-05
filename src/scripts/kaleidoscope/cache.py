## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

"""kaleidoscope.cache -- the run-reuse cache
(DESIGN 6.2.5; PSEUDOCODE 13.4).

The mechanism is owned by kaleidoscope; the key *fields* are
supplied by the client on each ``CalcUnit`` (DESIGN 6.2.5).
Before launching a unit the dispatcher asks ``is_cache_hit``;
on a hit it skips the launch and reports the existing outcome,
so resuming a flight is nothing more than re-running it.

The key has two parts, mirroring the producer's original
``is_cached_v2`` (DESIGN 5.7) and generalizing it:

- **scalar fields** -- written verbatim into ``cache_key.toml``
  and compared field-by-field;
- **key files** -- compared by *byte-comparison* against the
  copy staged at the same path relative to the run directory.
  No hashing: a developer can diff the files to see exactly why
  a cache missed, which a hash would hide.

A key file is checked for two different things, and keeping
them apart is what lets the cache be both correct and usable.
*Identity* asks whether this is the same calculation, and is
settled on the staged copy at the declared path.  *Agreement*
asks whether the engine would read what the key describes, and
is settled on the flattened root copy -- but only when that
copy exists, because which names reach the root depends on
what the unit's job reads.
"""

import filecmp
import os
import tomllib

from .workspace import emit_scalar, toml_line


def write_cache_key(wingbeat_dir, unit):
    """Write ``<wingbeat_dir>/cache_key.toml`` -- the identity
    snapshot taken when a unit is launched (DESIGN 6.2.5).
    Scalars are emitted under a ``[scalars]`` table in sorted
    order (so the file is deterministic); the key-file paths are
    recorded as a ``files`` array for inspection.  The actual
    file *contents* are not copied here -- the comparison reads
    the files staged in the run directory by the run itself.

    Nothing compares that ``files`` array, so changing what the
    declared paths are called cannot invalidate a stored unit;
    only the scalars are compared as a table, which is why a new
    scalar field is the expensive kind of change and a new key
    file is not."""
    os.makedirs(wingbeat_dir, exist_ok=True)
    path = os.path.join(wingbeat_dir, "cache_key.toml")
    names = [key_file.path for key_file in unit.key_fields.files]
    with open(path, "w") as key_file:
        rendered = ", ".join(emit_scalar(n) for n in names)
        key_file.write(f"files = [{rendered}]\n")
        key_file.write("\n[scalars]\n")
        for name in sorted(unit.key_fields.scalars):
            key_file.write(
                toml_line(name, unit.key_fields.scalars[name])
            )


def _read_cache_key(wingbeat_dir):
    """Return the parsed ``cache_key.toml``, or None if absent."""
    path = os.path.join(wingbeat_dir, "cache_key.toml")
    if not os.path.exists(path):
        return None
    with open(path, "rb") as key_file:
        return tomllib.load(key_file)


def cache_key_matches(unit, wingbeat_dir):
    """True iff the unit's *current* key matches the snapshot in
    the run directory (DESIGN 6.2.5).  The scalars must match
    field-by-field, and every declared key file must pass two
    separate tests: it must byte-equal its staged copy at the
    declared path (identity), and any flattened copy of it at
    the run directory's root must byte-equal that staged copy
    too (agreement).

    A key file that cannot be read -- on EITHER side -- is a
    miss, never an error.  Either side may legitimately be
    absent: a prepare directory reclaimed as scratch, a
    structure cache that moved between runs, a run directory
    left half-written by a job that died.  All of those mean the
    same thing, that this unit's identity cannot be established,
    and the answer is always to re-run it rather than trust it.
    Raising instead would let one unreadable file abort a
    campaign that had already paid for hours of converged rungs.
    """

    saved = _read_cache_key(wingbeat_dir)
    if saved is None:
        return False

    # Scalar fields: verbatim, field-by-field.
    if saved.get("scalars", {}) != unit.key_fields.scalars:
        return False

    for key_file in unit.key_fields.files:
        # IDENTITY.  Byte-compare the current source against the
        #   copy staged at the SAME relative path (no hashing).
        #   Both sides are checked for existence first, because
        #   filecmp raises rather than returning False on a
        #   missing file.
        staged = os.path.join(wingbeat_dir, key_file.path)
        if not os.path.exists(staged):
            return False
        if not os.path.exists(key_file.source):
            return False
        if not filecmp.cmp(key_file.source, staged,
                           shallow=False):
            return False

        # AGREEMENT.  A run directory also carries some staged
        #   files flattened at its root, because that is where
        #   the engine reads them, and which names appear there
        #   depends on what the unit's job reads.  So an ABSENT
        #   root copy is not a fault -- it means this job does
        #   not read that file, and identity is already settled
        #   above.  A root copy that is PRESENT and disagrees
        #   means the engine would run inputs the key does not
        #   describe, which must miss.
        #
        # A client declaring a bare filename lands this on the
        #   staged file itself, so the comparison is self against
        #   self and the test costs nothing.
        root_copy = os.path.join(
            wingbeat_dir, os.path.basename(key_file.path))
        if os.path.exists(root_copy):
            if not filecmp.cmp(root_copy, staged, shallow=False):
                return False
    return True


def is_cache_hit(unit, wingbeat_dir):
    """True iff this unit can be skipped: the run directory
    exists, its recorded key still matches, and its status is
    ``done`` (DESIGN 6.2.5).  Anything else -- no directory, a
    key mismatch, or a non-done status -- is a miss, so the unit
    is (re)launched."""
    from .workspace import read_status
    if not os.path.isdir(wingbeat_dir):
        return False
    status = read_status(wingbeat_dir)
    if status is None or status.get("status") != "done":
        return False
    return cache_key_matches(unit, wingbeat_dir)
