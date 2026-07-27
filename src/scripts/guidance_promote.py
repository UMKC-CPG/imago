#!/usr/bin/env python3
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

"""guidance_promote.py -- the curator's promotion helper for the
historical-guidance dataspace (DESIGN 7.8 curator half;
PSEUDOCODE 15.7).

Role in the pipeline
--------------------
``guidance_harvest.py`` stages a rich ``GuidanceEntry`` per
converged sweep into ``staging/<system_type>/`` (DESIGN 7.8).
Staging is a deliberate checkpoint: an automated harvest is not
the same as scientific endorsement, and a buggy harvest, a sweep
that converged at a numerical artifact, or a mis-classified
structure could all stage an entry that should not feed future
predictions.  This script is where a human curator (or an
objective rule standing in for one) decides which staged entries
graduate into ``entries/<system_type>/``, where the predictor
reads them.

Promotion is a pure ``mv`` -- the file's bytes (and therefore its
provenance) never change as it crosses the staging boundary, so a
promoted entry is byte-identical to what the curator reviewed.

Four modes (DESIGN 7.8)
-----------------------
* **interactive** (default) -- print each staged entry's summary
  and ask the curator to PROMOTE / SKIP / DELETE.
* **--auto-promote** -- promote every staged entry that passes the
  objective acceptance test (:func:`auto_promote_ok`), leaving the
  rest in staging for review.  In practice this clears ~80% of a
  seed flight (TODO C75) so the curator reviews only the ~20%
  outliers.
* **--all** -- promote every staged entry without checking the
  rule (for when the curator has already eyeballed the directory).
* **--dry-run** -- report what each mode *would* do, moving
  nothing.

The acceptance rule reads the staged file *alone* -- this is why
the harvest records ``grid_energies`` (DESIGN 7.2 / 7.8): the
curator helper never needs the original flight workspace.

Re-run dedup (DESIGN 7.8)
-------------------------
Promotion is also the only stage that sees both an incoming
entry and the promoted corpus, so it is where a re-run of an
already-promoted solid is caught.  The harvest cannot: it writes
one entry per converged solid with no view of what a curator
accepted months ago, and a re-run mints a fresh ``entry_id`` by
construction (the slug hashes flight id, structure, and
timestamp), so no collision check ever fires.

This matters because the predictor treats every stored entry as
an independent observation.  With five nearest neighbours, five
copies of one calculation fill the whole neighbour set, the
weighted variance collapses to zero, and a single measurement is
delivered as near-certainty -- which then drives the flight
builder into its narrowest search.

Two entries make the same claim when :func:`dedup_key` matches,
and they agree when their converged meshes match.  Agreement
retires the staged copy to ``superseded/``; disagreement, or a
mesh that cannot be compared, is reported and left in staging
for the curator.  Every mode applies this, ``--all`` included:
refusing to store one claim twice is a correctness guard, not a
quality judgment, and ``--all`` waives only the latter.

Checking against ``entries/`` does widen what the helper reads
beyond the single staged file.  The intent of that rule survives
-- what it avoided was depending on the flight *workspace*,
which is large, remote, and reclaimable; ``entries/`` is small,
local, and already the thing being written into.

This script reuses ``guidance_db`` for the schema: it loads and
validates each staged file with :func:`guidance_db.load_entry`,
so a malformed staging file fails loudly here too.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from datetime import datetime

from guidance_db import (
    GuidanceEntry,
    VALID_SYSTEM_TYPES,
    load_entry,
)
from guidance_harvest import per_atom_ev


# ==============================================================
#  The objective acceptance test (DESIGN 7.8 / PSEUDOCODE 15.7)
# ==============================================================

def auto_promote_ok(entry: GuidanceEntry) -> bool:
    """Return True when ``entry`` clears the objective acceptance
    test (DESIGN 7.8), evaluated from the staged file alone.

    Three conditions, all of which a trustworthy converged sweep
    should satisfy:

    1. **Converged in the middle 60% of the grid.**  A sweep that
       converged at either endpoint is suspect -- the grid may
       have been too narrow to bracket the true converged point --
       so the converged k-density must sit in [0.2, 0.8] of the
       grid's span.
    2. **A convincingly flat top of grid.**  The SPREAD (max minus
       min) of the top-three grid points' energies, taken per atom
       and in eV (the basis ``metric_threshold`` uses), must be
       below ``metric_threshold * 10`` -- a converged *region*, not
       just one consecutive delta that dipped below threshold.  A
       spread is a like-for-like linear quantity; a variance would
       be an energy squared against a linear threshold, so it is
       deliberately not used (DESIGN 7.8).
    3. **gap_ev / gap_kind consistent.**  ``gap_kind == "none"``
       iff ``gap_ev == 0.0`` (a defense even though the schema's
       rule 6 also enforces it).

    An entry with no verification block, or one whose verification
    omits ``grid_energies`` (a hand-written manual entry), can
    never satisfy condition 2 and so is never auto-promoted -- it
    falls to interactive review.
    """

    verification = entry.verification
    if verification is None or verification.grid_energies is None:
        return False

    # 1. Converged density in the middle 60% of the grid.
    low = verification.grid_values[0]
    high = verification.grid_values[-1]
    if high == low:
        return False
    position = (verification.converged_at - low) / (high - low)
    if not (0.2 <= position <= 0.8):
        return False

    # 2. Top-three grid points convincingly flat.  grid_energies are
    #    RAW total-cell hartree (Option B), so normalize each to eV
    #    per atom (the basis metric_threshold uses) via the same
    #    per_atom_ev helper pick_converged uses, then require their
    #    SPREAD (max - min) below metric_threshold * 10.
    atom_count = entry.context.cell_atom_count
    top_three = [per_atom_ev(energy, atom_count)
                 for energy in verification.grid_energies[-3:]]
    if (max(top_three) - min(top_three)
            >= verification.metric_threshold * 10.0):
        return False

    # 3. gap_ev / gap_kind consistent.
    is_metal = (entry.measured.gap_ev == 0.0)
    if (entry.measured.gap_kind == "none") != is_metal:
        return False

    return True


# ==============================================================
#  Re-run dedup (DESIGN 7.8; PSEUDOCODE 15.7)
# ==============================================================

def dedup_key(entry: GuidanceEntry) -> tuple:
    """The identity of the *claim* ``entry`` makes: same system,
    same settings, same structure (DESIGN 7.8).

    Two entries sharing this key are re-runs of one another, and
    the dataspace should hold only one of them -- the predictor
    counts each stored entry as an independent observation, so a
    solid run ten times would fill the whole neighbor set with
    copies of itself and report near-perfect confidence from a
    single measurement.

    Each part earns its place:

    - The three context fields are already the predictor's
      sub-model partition, so a ``gaussian`` and a
      ``gaussian-0.1`` run of one solid are *not* re-runs of each
      other and must both survive.
    - The structure's BASENAME, not its full path: the path
      records only where the structure cache happened to sit, and
      that location has moved (ARCHITECTURE 8.1).  The basename
      is ``<reference_id>-<cell>.skl``, which is the identity
      actually wanted -- a ``full`` and a ``prim`` run of one COD
      entry are genuinely different structures.
    - ``imago_commit`` is deliberately absent.  It is what the
      comparison examines, not what the key partitions on.
    """

    return (
        entry.signature.system_type,
        entry.context.basis,
        entry.context.functional,
        entry.context.kpoint_integration,
        os.path.basename(entry.provenance.source_structure),
    )


def mesh_of(entry: GuidanceEntry):
    """The entry's converged mesh, or ``None`` when it has none.

    ``converged_mesh`` is optional in the schema and a manual
    entry carries no verification block at all (DESIGN 7.2 /
    7.9).  Both absences read the same way to the dedup: an
    answer that cannot be compared, which falls to the conflict
    branch rather than being waved through.
    """

    if entry.verification is None:
        return None
    return entry.verification.converged_mesh


def load_promoted_entries(db_root: str) -> dict:
    """Read every promoted entry under ``db_root/entries/`` and
    return them keyed by :func:`dedup_key`.

    This is the one judgment :func:`auto_promote_ok` cannot make
    from the staged file alone, so the promotion pass loads the
    promoted corpus once up front (DESIGN 7.8).  ``entries/`` is
    small and local, so this is a directory read -- not the
    flight-workspace re-read the single-file rule exists to
    avoid.  ``staging/`` and ``superseded/`` are deliberately NOT
    read: only a *promoted* entry can hold a claim already.

    A later file silently wins a key collision here.  Two
    promoted entries sharing a claim is a pre-existing state this
    pass cannot fix (it has no verb for removing a promoted
    entry), and refusing to run at all would leave the curator
    with no way to promote anything else.
    """

    promoted: dict = {}
    for system_type in VALID_SYSTEM_TYPES:
        subdir = os.path.join(db_root, "entries", system_type)
        for path in sorted(glob.glob(
                os.path.join(subdir, "*.toml"))):
            # Fresh seen_ids per file: this pass validates each
            #   entry on its own, and the entry_id uniqueness rule
            #   is enforced by the loader that builds a dataspace.
            entry = load_entry(path, system_type, {})
            promoted[dedup_key(entry)] = entry
    return promoted


def _by_generated_at(first: tuple, second: tuple) -> tuple:
    """Order two ``(path, entry)`` pairs sharing a claim into
    ``(keep, drop)``: the later ``generated_at`` is kept.

    Timestamps are ISO-8601 UTC, so a string comparison orders
    them correctly.  Ties fall back to the path so the outcome is
    deterministic rather than dependent on directory order.
    """

    first_stamp = (first[1].generated_at, first[0])
    second_stamp = (second[1].generated_at, second[0])
    if first_stamp >= second_stamp:
        return first, second
    return second, first


# ==============================================================
#  Moving a staged file across the promotion boundary
# ==============================================================

def move_to_superseded(path: str, db_root: str,
                       system_type: str) -> str:
    """Retire one staged file that a promoted entry already
    claims: move it into ``superseded/<system_type>/`` under
    ``db_root`` and return the new path.

    Mirrors :func:`move_to_entries` -- a pure rename that refuses
    a pre-existing destination rather than overwriting.  Retired
    and NOT deleted, so the record of what a re-run produced stays
    recoverable and ``staging/`` does not accrete files that every
    later promotion pass re-examines (ARCHITECTURE 10.1).  The
    predictor never reads this directory.
    """

    destination_dir = os.path.join(
        db_root, "superseded", system_type)
    os.makedirs(destination_dir, exist_ok=True)
    destination = os.path.join(
        destination_dir, os.path.basename(path))
    if os.path.exists(destination):
        raise ValueError(
            "move_to_superseded: " + destination + " already "
            "exists (superseded collision); resolve by hand")
    os.rename(path, destination)
    return destination


def move_to_entries(path: str, db_root: str,
                    system_type: str) -> str:
    """Promote one staged file: move it from
    ``staging/<system_type>/`` to ``entries/<system_type>/`` under
    ``db_root`` and return the new path.

    A pure rename -- the bytes never change, so the entry the
    predictor loads is byte-identical to the one the curator
    reviewed (DESIGN 7.8).  A pre-existing file at the destination
    is refused rather than overwritten: the slug is a hash over the
    provenance (DESIGN 7.5), so a collision means a genuine
    duplicate that the curator must resolve by hand."""

    destination_dir = os.path.join(db_root, "entries", system_type)
    os.makedirs(destination_dir, exist_ok=True)
    destination = os.path.join(
        destination_dir, os.path.basename(path))
    if os.path.exists(destination):
        raise ValueError(
            "move_to_entries: " + destination + " already exists "
            "(promoted-entry collision); resolve by hand")
    os.rename(path, destination)
    return destination


# ==============================================================
#  Human-readable summary for interactive review and --dry-run
# ==============================================================

def format_summary(entry: GuidanceEntry) -> str:
    """Render a compact multi-line summary of one staged entry for
    the curator: the signature, the measured quantities, the
    verification grid, and the provenance (DESIGN 7.8 interactive
    review).  Purely for display -- not parsed back."""

    sig = entry.signature
    measured = entry.measured
    verification = entry.verification
    provenance = entry.provenance

    lines = [
        "=" * 60,
        "entry_id : " + entry.entry_id,
        "source   : " + entry.source
        + "   generated_at: " + entry.generated_at,
        "signature: " + sig.system_type
        + ("  lattice=" + sig.lattice_family
           if sig.lattice_family else ""),
        "measured : gap=" + repr(measured.gap_ev) + " eV ("
        + measured.gap_kind + ")  mag="
        + repr(measured.total_magnetization)
        + "  kpd=" + repr(measured.kpoint_density),
        "context  : " + entry.context.basis + "/"
        + entry.context.functional + "/"
        + entry.context.kpoint_integration
        + "  atoms=" + repr(entry.context.cell_atom_count),
    ]
    if verification is not None:
        lines.append(
            "verify   : grid=" + repr(list(verification.grid_values))
            + "  converged_at=" + repr(verification.converged_at))
        lines.append(
            "           confidence="
            + repr(verification.predictor_confidence)
            + "  neighbors="
            + repr(list(verification.predictor_neighbor_ids)))
    lines.append(
        "provenance: flight=" + provenance.flight_id
        + "  structure=" + provenance.source_structure
        + "  commit=" + provenance.imago_commit)
    return "\n".join(lines)


def _format_conflict(staged: GuidanceEntry,
                     prior: GuidanceEntry) -> str:
    """Render the CONFLICT report: one claim, two converged
    meshes (DESIGN 7.8).

    Both sides are shown with their mesh and their Imago commit,
    because a mesh that changed across a commit is the expected
    shape of this -- the code's behaviour moved -- and which
    answer is right is then a physics judgment.  The report names
    what to compare rather than choosing, and says plainly that
    nothing was moved.
    """

    def describe(entry: GuidanceEntry) -> str:
        mesh = mesh_of(entry)
        return ("mesh=" + (repr(list(mesh)) if mesh is not None
                           else "<none recorded>")
                + "  commit=" + entry.provenance.imago_commit
                + "  generated=" + entry.generated_at)

    return "\n".join([
        "=" * 60,
        "CONFLICT: " + staged.entry_id + " re-runs a promoted "
        "claim with a different answer.",
        "  structure : "
        + os.path.basename(staged.provenance.source_structure),
        "  settings  : " + staged.context.basis + "/"
        + staged.context.functional + "/"
        + staged.context.kpoint_integration,
        "  promoted  : " + prior.entry_id + "  " + describe(prior),
        "  staged    : " + staged.entry_id + "  "
        + describe(staged),
        "  Nothing moved.  Compare the two runs and resolve by "
        "hand: keep",
        "  the promoted entry and delete the staged file, or "
        "remove the",
        "  promoted entry first and re-run this promotion.",
    ])


# ==============================================================
#  The promotion driver (DESIGN 7.8; PSEUDOCODE 15.7)
# ==============================================================

VALID_MODES = ("interactive", "auto-promote", "all", "dry-run")


def _ask_choice(ask) -> str:
    """Prompt the curator until they answer PROMOTE, SKIP, or
    DELETE (accepting the initials p / s / d, case-insensitively),
    and return the canonical word.  ``ask`` is injected so tests
    can feed canned responses; an empty answer defaults to SKIP
    (the safe, non-destructive choice)."""

    while True:
        answer = ask("PROMOTE / SKIP / DELETE [p/s/d]: ")
        normalized = answer.strip().lower()
        if normalized in ("p", "promote"):
            return "PROMOTE"
        if normalized in ("s", "skip", ""):
            return "SKIP"
        if normalized in ("d", "delete"):
            return "DELETE"
        # Anything else: re-prompt rather than guess.


def promote(db_root: str, mode: str = "interactive", *,
            ask=input, output=print):
    """Walk every ``staging/<system_type>/`` file under ``db_root``
    and apply ``mode`` (DESIGN 7.8; PSEUDOCODE 15.7).

    Returns a list of ``(entry_id, action)`` records -- where
    ``action`` is one of ``promoted`` / ``skipped`` / ``deleted``
    / ``superseded`` / ``conflicted`` (the acting modes) or
    ``would-promote`` / ``would-skip`` / ``would-supersede``
    (``dry-run``) -- so a caller or a test can see what happened
    without parsing printed text.  ``ask`` and ``output`` are
    injected so the interactive flow is testable; the default
    wiring is the real ``input`` / ``print``.

    Files are processed in sorted order within each system_type so
    the run is deterministic.  Each file is validated through
    ``guidance_db.load_entry`` before any decision, so a malformed
    staging file aborts the run loudly (naming the file).

    Before any mode runs, each staged entry is checked against the
    promoted corpus for a re-run of a claim already held (DESIGN
    7.8).  A re-run whose converged mesh agrees is retired to
    ``superseded/``; one whose mesh disagrees -- or cannot be
    compared -- is reported and left in staging for the curator,
    never promoted automatically.  Every mode applies this,
    ``--all`` included: refusing to store one claim twice is a
    correctness guard, and ``--all`` waives only the quality
    rule."""

    if mode not in VALID_MODES:
        raise ValueError(
            "unknown promote mode " + repr(mode)
            + " (one of " + repr(VALID_MODES) + ")")

    # The promoted corpus, keyed by claim, read once up front.
    promoted = load_promoted_entries(db_root)
    # dry-run evaluates every decision below and moves nothing, so
    #   a curator sees the whole outcome -- promotions,
    #   retirements, and conflicts -- before a file is touched.
    dry = (mode == "dry-run")

    results = []
    for system_type in VALID_SYSTEM_TYPES:
        subdir = os.path.join(db_root, "staging", system_type)
        staged = sorted(glob.glob(os.path.join(subdir, "*.toml")))

        # Resolve re-runs WITHIN this batch first, so at most one
        #   candidate per claim is compared against entries/.  This
        #   is where promotion stops judging each staged file in
        #   isolation: under the dedup rule staging IS a uniqueness
        #   namespace, where previously it was not.
        batch: dict = {}                 # dedup_key -> (path, entry)
        for path in staged:
            # Fresh seen_ids per file: entry_id uniqueness across
            #   the corpus is the dataspace loader's rule, and a
            #   cross-file entry_id clash here would be a false
            #   positive.  The CLAIM-level check is dedup_key.
            entry = load_entry(path, system_type, {})
            key = dedup_key(entry)
            if key not in batch:
                batch[key] = (path, entry)
                continue
            keep, drop = _by_generated_at(batch[key], (path, entry))
            batch[key] = keep
            drop_path, drop_entry = drop
            if not dry:
                move_to_superseded(drop_path, db_root, system_type)
            output(drop_entry.entry_id + ": "
                   + ("would be superseded" if dry else "superseded")
                   + " -- a later run in this batch makes the same "
                   "claim")
            results.append(
                (drop_entry.entry_id,
                 "would-supersede" if dry else "superseded"))

        # Sort on the path alone: GuidanceEntry is frozen but not
        #   ordered, so letting a tuple comparison fall through to
        #   it would raise the day two paths ever tie.
        for path, entry in sorted(batch.values(),
                                  key=lambda item: item[0]):
            prior = promoted.get(dedup_key(entry))
            if prior is not None:
                prior_mesh = mesh_of(prior)
                if (prior_mesh is not None
                        and prior_mesh == mesh_of(entry)):
                    # REDUNDANT.  The promoted entry stands
                    #   untouched -- promotion only ever ADDS to
                    #   entries/, so an entry the curator reviewed
                    #   stays byte-identical for as long as it
                    #   lives there.  The staged copy is retired.
                    if not dry:
                        move_to_superseded(
                            path, db_root, system_type)
                    output(entry.entry_id + ": "
                           + ("would be superseded" if dry
                              else "superseded")
                           + " -- already promoted as "
                           + prior.entry_id)
                    results.append(
                        (entry.entry_id,
                         "would-supersede" if dry else "superseded"))
                    continue
                # CONFLICT: one claim, two answers.  Never
                #   automatic in ANY mode -- deciding which is
                #   right is a physics judgment, not a timestamp
                #   comparison, and resolving it may mean removing
                #   an entry from entries/, which this tool has no
                #   verb for.  Report both and leave the staged
                #   file alone.
                output(_format_conflict(entry, prior))
                results.append((entry.entry_id, "conflicted"))
                continue

            if mode == "dry-run":
                ok = auto_promote_ok(entry)
                output(format_summary(entry))
                output("  WOULD " + ("PROMOTE" if ok
                       else "SKIP (fails auto-promote rule)"))
                results.append(
                    (entry.entry_id,
                     "would-promote" if ok else "would-skip"))

            elif mode == "all":
                new_path = move_to_entries(
                    path, db_root, system_type)
                output(entry.entry_id + ": promoted -> " + new_path)
                results.append((entry.entry_id, "promoted"))

            elif mode == "auto-promote":
                if auto_promote_ok(entry):
                    new_path = move_to_entries(
                        path, db_root, system_type)
                    output(entry.entry_id
                           + ": promoted -> " + new_path)
                    results.append((entry.entry_id, "promoted"))
                else:
                    output(entry.entry_id
                           + ": left in staging "
                           "(fails auto-promote rule)")
                    results.append((entry.entry_id, "skipped"))

            else:                                # interactive
                output(format_summary(entry))
                choice = _ask_choice(ask)
                if choice == "PROMOTE":
                    new_path = move_to_entries(
                        path, db_root, system_type)
                    output("  promoted -> " + new_path)
                    results.append((entry.entry_id, "promoted"))
                elif choice == "DELETE":
                    os.remove(path)
                    output("  deleted")
                    results.append((entry.entry_id, "deleted"))
                else:                            # SKIP
                    output("  skipped (left in staging)")
                    results.append((entry.entry_id, "skipped"))

    return results


# ==============================================================
#  Command-line interface
# ==============================================================

def _default_db_root() -> str:
    """The dataspace root under $IMAGO_DATA (DESIGN 7 layout:
    ``share/historicalGuidanceDB/``); empty when unset so the
    parser can demand ``--db-root``."""

    data_dir = os.environ.get("IMAGO_DATA", "")
    return (os.path.join(data_dir, "historicalGuidanceDB")
            if data_dir else "")


def main(argv=None):
    """CLI entry point: review and promote staged guidance entries
    (DESIGN 7.8)."""

    parser = argparse.ArgumentParser(
        description="Review and promote staged historical-guidance "
                    "entries from staging/ into entries/.  With no "
                    "mode flag the default is interactive review of "
                    "each staged entry.  In every mode, a staged "
                    "entry that re-runs a solid already promoted "
                    "under the same settings is retired to "
                    "superseded/ when its converged mesh agrees, "
                    "and reported as a conflict and left in "
                    "staging when it does not.")
    parser.add_argument(
        "--db-root", default=_default_db_root(),
        help="the historicalGuidanceDB root "
             "(default: $IMAGO_DATA/historicalGuidanceDB)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--auto-promote", action="store_const", dest="mode",
        const="auto-promote",
        help="promote only entries passing the objective rule; "
             "leave the rest in staging")
    group.add_argument(
        "--all", action="store_const", dest="mode", const="all",
        help="promote every staged entry without checking the rule")
    group.add_argument(
        "--dry-run", action="store_const", dest="mode",
        const="dry-run",
        help="report what would happen; move nothing")
    parser.set_defaults(mode="interactive")
    args = parser.parse_args(argv)

    if not args.db_root:
        parser.error(
            "--db-root not given and $IMAGO_DATA is unset")

    promote(args.db_root, args.mode)
    return 0


def record_command():
    """Append the issued command line to a file named "command" in
    the current directory, so the exact invocation can be recovered
    later.  This is a standing project convention: each run appends
    a dated block, so the file builds up a history of how the script
    was called."""

    with open("command", "a") as cmd:
        now = datetime.now()
        stamp = now.strftime("%b. %d, %Y: %H:%M:%S")
        cmd.write(f"Date: {stamp}\n")
        cmd.write("Cmnd:")
        for argument in sys.argv:
            cmd.write(f" {argument}")
        cmd.write("\n\n")


if __name__ == "__main__":
    record_command()
    sys.exit(main())
