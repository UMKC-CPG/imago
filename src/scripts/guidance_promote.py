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

This script reuses ``guidance_db`` for the schema: it loads and
validates each staged file with :func:`guidance_db.load_entry`,
so a malformed staging file fails loudly here too.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

from guidance_db import (
    GuidanceEntry,
    VALID_SYSTEM_TYPES,
    load_entry,
)


# ==============================================================
#  The objective acceptance test (DESIGN 7.8 / PSEUDOCODE 15.7)
# ==============================================================

def _population_variance(values) -> float:
    """The population variance (mean of squared deviations) of a
    short sequence.  Fewer than two values are perfectly "flat",
    so the variance is 0.0 -- which is what the auto-promote
    flatness test wants for a degenerate top-of-grid."""

    sample = list(values)
    count = len(sample)
    if count < 2:
        return 0.0
    mean = sum(sample) / count
    return sum((value - mean) ** 2 for value in sample) / count


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
    2. **A convincingly flat top of grid.**  The population
       variance of the top-three grid points' total energies must
       be below ``metric_threshold * 10`` -- a converged *region*,
       not just one consecutive delta that dipped below threshold.
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

    # 2. Top-three grid points convincingly flat.
    top_three = verification.grid_energies[-3:]
    if (_population_variance(top_three)
            >= verification.metric_threshold * 10.0):
        return False

    # 3. gap_ev / gap_kind consistent.
    is_metal = (entry.measured.gap_ev == 0.0)
    if (entry.measured.gap_kind == "none") != is_metal:
        return False

    return True


# ==============================================================
#  Moving a staged file across the promotion boundary
# ==============================================================

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
    (the acting modes) or ``would-promote`` / ``would-skip``
    (``dry-run``) -- so a caller or a test can see what happened
    without parsing printed text.  ``ask`` and ``output`` are
    injected so the interactive flow is testable; the default
    wiring is the real ``input`` / ``print``.

    Files are processed in sorted order within each system_type so
    the run is deterministic.  Each file is validated through
    ``guidance_db.load_entry`` before any decision, so a malformed
    staging file aborts the run loudly (naming the file)."""

    if mode not in VALID_MODES:
        raise ValueError(
            "unknown promote mode " + repr(mode)
            + " (one of " + repr(VALID_MODES) + ")")

    results = []
    for system_type in VALID_SYSTEM_TYPES:
        subdir = os.path.join(db_root, "staging", system_type)
        for path in sorted(glob.glob(
                os.path.join(subdir, "*.toml"))):
            # Fresh seen_ids per file: staging is not the unique-id
            #   namespace (entries/ is), and a cross-file dup here
            #   would be a false positive.
            entry = load_entry(path, system_type, {})

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
                    "entries from staging/ into entries/ "
                    "(DESIGN 7.8).")
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


if __name__ == "__main__":
    sys.exit(main())
