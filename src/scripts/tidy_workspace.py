#!/usr/bin/env python3
"""tidy_workspace.py -- reclaim the intermediate scratch a
kaleidoscope flight leaves behind (DESIGN 6.2.12; PSEUDOCODE 13.8).

A finished run directory spans two filesystems, and the difference
between them is almost entirely one of size.  The **kept** tier is
the run directory itself: the staged inputs, ``result.toml``,
``status.toml``, ``cache_key.toml``, the SCF potential, the log --
a couple of hundred kilobytes.  The **scratch** tier sits behind
the ``intermediate`` symlink that ``imago.py`` creates, pointing
into ``$IMAGO_TEMP``, and holds the engine's working files, above
all the HDF5 carrying the wavefunctions.  Measured on a seed-scale
producer run, scratch was 99.7% of the bytes -- 3.17 GB of 3.2 GB,
about 25 MB per calculation.

Reclaiming it is safe, and this is why.  The harvest reads only
the paths recorded in ``result.toml``'s ``outputs`` table, and
every one of them resolves inside the run directory; none points
through ``intermediate``.  The run-reuse cache decides a hit from
``status.toml`` and ``cache_key.toml`` alone.  So a reclaimed run
directory still answers every question the producer asks of it --
it still harvests, and it still counts as a cache hit.

Three refusals define the tool
------------------------------
It deletes, so what it will *not* do matters more than what it
will:

1. **Only scratch is removed, never the run directory.**  The kept
   tier is the record of the calculation and is two orders of
   magnitude smaller than the saving.  The ``intermediate`` link
   is deliberately left dangling, so a reclaimed run still shows
   where its scratch was.
2. **An unfinished unit is never touched.**  A ``running`` or
   ``queued`` status means the engine may still be writing, and a
   missing ``result.toml`` means the run produced nothing to keep
   -- usually the state a curator most wants to investigate.
3. **A link out of the scratch area is never followed.**  The
   target is resolved and checked to lie under the scratch root
   before anything is removed.  ``imago.py`` renames a stale link
   to ``intermediateFIXME``, and a hand-edited workspace can point
   anywhere; a cleanup tool a symlink can redirect is a hazard.

**Dry run is the default.**  The tool previews what it would
remove, with sizes, and removes nothing until ``--apply`` is
given.  An operation whose whole purpose is deletion should make
the destructive path the one the user typed on purpose.

Usage examples
--------------
Preview everything reclaimable in a workspace::

    tidy_workspace.py share/curation/workspace

Reclaim it for real::

    tidy_workspace.py share/curation/workspace --apply

Reclaim only one solid's convergence rungs, leaving its loen
descriptor runs alone::

    tidy_workspace.py <root> --id si_cmce_64_1999 \\
        --calc 'kpt-mesh-*' --apply

Reclaim only what has been sitting around for a week::

    tidy_workspace.py <root> --older-than 7 --apply
"""

import argparse
import fnmatch
import os
import shutil
import sys
import tomllib
from datetime import datetime


# The kaleidoscope workspace holds one run directory per
#   calculation under this subdirectory (DESIGN 6.2.4).
WINGBEATS = "wingbeats"

# The symlink imago.py plants in each run directory, pointing at
#   that run's scratch area (imago.py init_directories).
INTERMEDIATE = "intermediate"

SECONDS_PER_DAY = 86400.0


# ============================================================
#  Reclamation policy (the client's half; DESIGN 6.2.12)
# ============================================================

def default_reclaim_policy(run_dir):
    """Return ``(spent, reason)`` for one run directory under the
    conservative default policy.

    A unit is *spent* -- done with its scratch -- once it finished
    AND left a result.  ``done`` alone is deliberately not enough:
    a run that completed but wrote no ``result.toml`` is usually
    the state a curator most wants to look at, so its working
    files are preserved rather than reclaimed.

    A client whose flow needs the working files for longer (one
    that post-processes wavefunctions into a density of states,
    say) passes its own policy instead; only the client knows when
    a finished run is finished *with*.
    """

    status_path = os.path.join(run_dir, "status.toml")
    if not os.path.isfile(status_path):
        return False, "no status.toml"
    try:
        with open(status_path, "rb") as handle:
            status = tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        return False, f"unreadable status.toml ({exc})"

    state = status.get("status")
    if state != "done":
        return False, f"status is {state!r}, not done"
    if not os.path.isfile(os.path.join(run_dir, "result.toml")):
        return False, "no result.toml"
    return True, ""


# ============================================================
#  Scratch resolution (the refusal that matters most)
# ============================================================

def scratch_target(run_dir, scratch_root):
    """Resolve a run directory's scratch, or refuse it.

    Returns ``(target, reason)``: the resolved scratch path when it
    is safe to act on, else ``(None, reason)``.

    The check that earns its keep is the last one.  ``imago.py``
    renames a stale ``intermediate`` link to ``intermediateFIXME``
    and plants a fresh one, and a hand-edited workspace can point
    the link anywhere at all -- so the resolved target must be
    confirmed to lie under the scratch root before anything is
    removed.  A cleanup tool that a symlink can redirect is a
    hazard rather than a convenience.
    """

    link = os.path.join(run_dir, INTERMEDIATE)
    if not os.path.islink(link):
        return None, "no intermediate link"

    target = os.path.realpath(link)
    if not os.path.isdir(target):
        return None, "already reclaimed"

    root = os.path.realpath(scratch_root)
    # commonpath rather than startswith: a plain prefix test would
    #   accept "/scratch-other" for a root of "/scratch".
    try:
        shared = os.path.commonpath([target, root])
    except ValueError:              # different drives / relative
        return None, f"target not under scratch root: {target}"
    if shared != root:
        return None, f"target not under scratch root: {target}"
    return target, ""


def tree_size(path):
    """Total bytes under ``path``, following no symlinks.

    Symlinks are counted as their own (tiny) size rather than
    followed, so a stray link inside scratch cannot inflate the
    reported saving or, worse, lead the walk outside the tree.
    """

    total = 0
    for directory, _, files in os.walk(path, followlinks=False):
        for name in files:
            full = os.path.join(directory, name)
            try:
                total += os.lstat(full).st_size
            except OSError:
                pass                # vanished mid-walk; ignore
    return total


# ============================================================
#  The walk (the mechanism half; PSEUDOCODE 13.8)
# ============================================================

def find_run_dirs(root):
    """Yield every run directory under ``<root>/wingbeats``.

    A run directory is one that carries a ``status.toml``: the
    ``<calc>`` level is optional (DESIGN 6.2.4), so a unit may sit
    directly under its id or one or more levels below it, and
    keying on the status file finds either without assuming a
    depth.
    """

    wingbeats = os.path.join(root, WINGBEATS)
    if not os.path.isdir(wingbeats):
        return
    for directory, subdirs, files in os.walk(wingbeats):
        # Never descend through a symlink -- notably `intermediate`
        #   itself, which would walk us into the scratch area.
        subdirs[:] = [name for name in subdirs
                      if not os.path.islink(
                          os.path.join(directory, name))]
        if "status.toml" in files:
            yield directory


def plan_reclamation(root, scratch_root, policy=None,
                     ids=None, calc_pattern=None, older_than=None):
    """Build the reclamation plan for a workspace: one record per
    run directory, saying whether its scratch can go and why not
    when it cannot (PSEUDOCODE 13.8 ``reclaim_workspace``).

    Nothing is removed here.  The plan is the report, and applying
    it is a separate, explicit step -- which is what lets the same
    code path serve the preview, the standalone run, and the
    producer's ``--clean-after``.

    :param root: the flight workspace (the directory holding
        ``wingbeats/``).
    :param scratch_root: the area scratch must lie under, normally
        ``$IMAGO_TEMP``.
    :param policy: ``(run_dir) -> (spent, reason)``; defaults to
        :func:`default_reclaim_policy`.
    :param ids: restrict to these stable ids, or None for all.
    :param calc_pattern: a glob matched against the calc tag (the
        run directory's name), or None for all.
    :param older_than: only scratch untouched for at least this
        many days, or None for any age.
    """

    if policy is None:
        policy = default_reclaim_policy
    wingbeats = os.path.abspath(os.path.join(root, WINGBEATS))
    now = datetime.now().timestamp()
    plan = []

    for run_dir in sorted(find_run_dirs(root)):
        relative = os.path.relpath(os.path.abspath(run_dir),
                                   wingbeats)
        unit_id = relative.split(os.sep)[0]
        calc_tag = os.path.basename(run_dir)

        if ids is not None and unit_id not in ids:
            continue
        if (calc_pattern is not None
                and not fnmatch.fnmatch(calc_tag, calc_pattern)):
            continue

        spent, reason = policy(run_dir)
        if not spent:
            plan.append(dict(run_dir=run_dir, unit=relative,
                             ok=False, reason=reason, bytes=0))
            continue

        target, reason = scratch_target(run_dir, scratch_root)
        if target is None:
            plan.append(dict(run_dir=run_dir, unit=relative,
                             ok=False, reason=reason, bytes=0))
            continue

        if older_than is not None:
            age_days = (now - os.lstat(target).st_mtime) \
                / SECONDS_PER_DAY
            if age_days < older_than:
                plan.append(dict(
                    run_dir=run_dir, unit=relative, ok=False,
                    reason=f"only {age_days:.1f} days old", bytes=0))
                continue

        plan.append(dict(run_dir=run_dir, unit=relative, ok=True,
                         reason="", target=target,
                         bytes=tree_size(target)))
    return plan


def apply_reclamation(plan):
    """Remove the scratch of every reclaimable entry in ``plan``.

    Only the scratch tree is removed.  The run directory is left
    untouched, and its ``intermediate`` link is left in place and
    dangling on purpose, so the run still records where its
    scratch was (DESIGN 6.2.12).

    Returns ``(removed_count, removed_bytes, failures)``, where a
    failure is a ``(unit, message)`` pair: one unremovable tree
    (a permission problem, a busy filesystem) must not abandon the
    rest of the campaign.
    """

    removed, freed, failures = 0, 0, []
    for item in plan:
        if not item["ok"]:
            continue
        try:
            shutil.rmtree(item["target"])
        except OSError as exc:
            failures.append((item["unit"], str(exc)))
            continue
        removed += 1
        freed += item["bytes"]
    return removed, freed, failures


# ============================================================
#  Reporting
# ============================================================

def human_bytes(count):
    """Format a byte count for the report, in the largest unit
    that keeps the number readable."""

    size = float(count)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024.0 or unit == "TB":
            return f"{size:.1f} {unit}" if unit != "B" \
                else f"{int(size)} B"
        size /= 1024.0
    return f"{size:.1f} TB"


def print_report(plan, applied, removed=0, freed=0, failures=(),
                 verbose=False):
    """Print the plan, then the tally.

    Skips are summarized by reason rather than listed one per
    line: a workspace mid-flight has many units in the same state,
    and a hundred identical lines hide the one that matters.
    ``--verbose`` lists them individually.
    """

    reclaimable = [item for item in plan if item["ok"]]
    skipped = [item for item in plan if not item["ok"]]

    for item in reclaimable:
        mark = "removed" if applied else "would free"
        print(f"  {mark:11} {human_bytes(item['bytes']):>9}  "
              f"{item['unit']}")

    if skipped:
        if verbose:
            for item in skipped:
                print(f"  {'skipped':11} {'':>9}  {item['unit']} "
                      f"({item['reason']})")
        else:
            reasons = {}
            for item in skipped:
                reasons[item["reason"]] = \
                    reasons.get(item["reason"], 0) + 1
            for reason, count in sorted(reasons.items()):
                print(f"  {'skipped':11} {count:>9}  units: {reason}")

    total = sum(item["bytes"] for item in reclaimable)
    print()
    if applied:
        print(f"reclaimed {removed} of {len(reclaimable)} run "
              f"directories, freeing {human_bytes(freed)}")
        for unit, message in failures:
            print(f"  FAILED  {unit}: {message}")
    else:
        print(f"{len(reclaimable)} run directories reclaimable, "
              f"{human_bytes(total)} recoverable "
              f"({len(skipped)} skipped)")
        if reclaimable:
            print("re-run with --apply to remove them")


# ============================================================
#  CLI
# ============================================================

def build_parser():
    parser = argparse.ArgumentParser(
        description="Reclaim the intermediate scratch a "
                    "kaleidoscope flight leaves behind.  Removes "
                    "only scratch, never a run directory, never an "
                    "unfinished unit, and never a link pointing "
                    "outside the scratch area.  Previews by "
                    "default and removes nothing until --apply.")
    parser.add_argument(
        "root",
        help="the flight workspace to tidy: the directory holding "
             "wingbeats/ (for the producer, "
             "<data_root>/curation/workspace)")
    parser.add_argument(
        "--apply", action="store_true",
        help="actually remove the scratch (default: preview only, "
             "removing nothing)")
    parser.add_argument(
        "--scratch-root", default=os.environ.get("IMAGO_TEMP", ""),
        help="the area scratch must lie under; a resolved target "
             "outside it is refused (default: $IMAGO_TEMP)")
    parser.add_argument(
        "--id", action="append", dest="ids", metavar="ID",
        help="restrict to this stable id; repeatable (default: "
             "every id in the workspace)")
    parser.add_argument(
        "--calc", dest="calc_pattern", metavar="GLOB",
        help="restrict to run directories whose calc tag matches "
             "this glob, e.g. 'kpt-mesh-*' or 'loen-*' (default: "
             "every calc)")
    parser.add_argument(
        "--older-than", type=float, metavar="DAYS",
        help="restrict to scratch untouched for at least this many "
             "days (default: any age)")
    parser.add_argument(
        "--verbose", action="store_true",
        help="list every skipped unit individually instead of "
             "summarizing them by reason (default: summarize)")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    if not args.scratch_root:
        print("error: no scratch root: pass --scratch-root or set "
              "$IMAGO_TEMP.  Without it the safety check that "
              "scratch lies where it should cannot run.",
              file=sys.stderr)
        return 2
    if not os.path.isdir(os.path.join(args.root, WINGBEATS)):
        print(f"error: {args.root!r} holds no {WINGBEATS}/ "
              f"directory, so it is not a flight workspace.",
              file=sys.stderr)
        return 2

    plan = plan_reclamation(
        args.root, args.scratch_root, ids=(
            set(args.ids) if args.ids else None),
        calc_pattern=args.calc_pattern,
        older_than=args.older_than)

    if not args.apply:
        print_report(plan, applied=False, verbose=args.verbose)
        return 0

    removed, freed, failures = apply_reclamation(plan)
    print_report(plan, applied=True, removed=removed, freed=freed,
                 failures=failures, verbose=args.verbose)
    return 1 if failures else 0


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
