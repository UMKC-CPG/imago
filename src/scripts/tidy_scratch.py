#!/usr/bin/env python3
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

"""tidy_scratch.py -- reclaim the intermediate scratch an imago
run leaves behind (DESIGN 6.2.12; PSEUDOCODE 13.8).

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

Two kinds of root, one per call
-------------------------------
Scratch is not only a flight's problem.  Every ordinary
``imago.py`` run plants the same ``intermediate`` link and leaves
the same tens of megabytes behind it, so the tool recognizes two
kinds of root and decides which it has by looking:

* a **workspace** holds ``wingbeats/``; its units are kaleidoscope
  units, and each proves it is finished in ``status.toml``;
* a **job tree** holds ordinary run directories with no flight
  above them -- a student's job directory, a hand-driven
  convergence test.  These prove they are finished by holding no
  ``imagoLock`` and ending their ``runtime`` log with the
  completion marker.

One call handles exactly one kind.  Mixing them would gather two
different safety contracts under one set of totals, and they do
genuinely differ (see the refusals).

Why reclaiming a workspace unit is safe
---------------------------------------
The harvest reads only the paths recorded in ``result.toml``'s
``outputs`` table, and every one of them resolves inside the run
directory; none points through ``intermediate``.  The run-reuse
cache decides a hit from ``status.toml`` and ``cache_key.toml``
alone.  So a reclaimed run directory still answers every question
the producer asks of it -- it still harvests, and it still counts
as a cache hit.

Five refusals hold for every root
---------------------------------
It deletes, so what it will *not* do matters more than what it
will.  These five apply to a workspace and a job tree alike:

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
4. **A symlink is never descended while walking.**  The third
   refusal guards what is *removed*; this one guards what is even
   *considered*.  ``intermediate`` is itself a link into scratch,
   so a walk that followed links would leave the workspace and
   could plan a removal outside the tree it was pointed at.  The
   walk skips every symlinked subdirectory, keeping the candidate
   run directories inside the workspace by construction.
5. **A tree holding another run's scratch is deferred.**  Scratch
   mirrors the run directory's path, so a run nested inside
   another has its scratch nested too.  Removing the outer tree
   would take the inner one with it -- deleting working files the
   other refusals had just declined to touch.  It is skipped
   while any other run's scratch lies within it; once those are
   gone, a second pass reclaims it.

And two more for a job tree
---------------------------
A hand run writes no ``status.toml``, so the authority a
workspace unit carries is absent.  It is not silent, though: the
CLI is a thin wrapper over the same callable core, so an
``imago.py`` run leaves the same traces a flight unit does, and
two extra refusals rest on them:

6. **A run that has not declared it finished is never
   reclaimed.**  Its scratch must hold no ``imagoLock`` *and* its
   ``runtime`` log must end with the completion marker.  Absence
   of evidence is refusal: no log at all, or one ending
   mid-stream, is reported and left alone.
7. **A workspace is never descended into from a job tree.**
   Walking in would judge units that have a provable contract by
   the presumption-based one, so such directories are named in
   the report and pruned from the walk.

**Dry run is the default.**  The tool previews what it would
remove, with sizes, and removes nothing until ``--apply`` is
given.  An operation whose whole purpose is deletion should make
the destructive path the one the user typed on purpose.

Usage examples
--------------
Preview everything reclaimable in a workspace::

    tidy_scratch.py share/curation/workspace

Reclaim it for real::

    tidy_scratch.py share/curation/workspace --apply

Reclaim only one solid's convergence rungs, leaving its loen
descriptor runs alone::

    tidy_scratch.py <root> --id si_cmce_64_1999 \\
        --calc 'kpt-mesh-*' --apply

Sweep up after a season of hand runs (a job tree; selects on the
path relative to the root rather than on unit ids)::

    tidy_scratch.py ~/imago/jobs --match 'c/diamond/*' --apply

Reclaim only what has been sitting around for a week::

    tidy_scratch.py <root> --older-than 7 --apply
"""

import argparse
import fnmatch
import os
import shutil
import sys
import tomllib
from datetime import datetime

import imago


# The four names the job-tree contract shares with the engine are
#   sourced from imago.py rather than re-spelled here, so that a
#   change to how imago.py names or marks a run cannot silently
#   desync this recognizer (DESIGN 6.2.12).  The failure would be
#   in the safe direction -- "absence of evidence is refusal", so
#   an unrecognized marker refuses rather than deletes -- but it
#   would be silent, which is harder to diagnose than a loud one.
#
# INTERMEDIATE: the symlink imago.py plants in each run directory,
#   pointing at that run's scratch area (imago.py
#   init_directories).
INTERMEDIATE = imago.INTERMEDIATE_LINK

# IMAGO_LOCK: the per-run lock imago.py takes inside the SCRATCH
#   area before any work begins, and releases in the cleanup that
#   always runs (imago.py _run_core).  Its presence means the run
#   owns the directory now, or died without releasing it.
IMAGO_LOCK = imago.LOCK_FILE

# RUNTIME_LOG, COMPLETION_MARKER: the log imago.py keeps in the
#   RUN directory, and the exact line it appends as it closes that
#   log.  The log is opened in append mode, so a directory run
#   repeatedly accumulates one marker per run -- only the tail is
#   truthful (see hand_run_policy).
RUNTIME_LOG = imago.RUNTIME_FILE
COMPLETION_MARKER = imago.COMPLETION_MARKER

# The kaleidoscope workspace holds one run directory per
#   calculation under this subdirectory (DESIGN 6.2.4).  This one
#   is a kaleidoscope concept, not an imago filename, so it is
#   defined here.
WINGBEATS = "wingbeats"

# The two kinds of root a call may be pointed at (DESIGN 6.2.12).
WORKSPACE_ROOT = "workspace"
JOB_TREE_ROOT = "job-tree"

SECONDS_PER_DAY = 86400.0


# ============================================================
#  Reclamation policy (the client's half; DESIGN 6.2.12)
# ============================================================

def default_reclaim_policy(run_dir, target=None):
    """Return ``(spent, reason)`` for one WORKSPACE unit under the
    conservative default policy.

    A unit is *spent* -- done with its scratch -- once it finished
    AND left a result.  ``done`` alone is deliberately not enough:
    a run that completed but wrote no ``result.toml`` is usually
    the state a curator most wants to look at, so its working
    files are preserved rather than reclaimed.

    ``target`` is the run's already-resolved scratch.  It is
    unused here -- ``status.toml`` is authority enough -- but it
    is part of the policy contract because the job-tree policy
    below must look inside the scratch, and a client policy is no
    worse off for being handed the path it is deciding about.

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


def hand_run_policy(run_dir, target):
    """Return ``(spent, reason)`` for one JOB-TREE run directory
    -- an ordinary ``imago.py`` run with no flight above it
    (refusal 6).

    A hand run writes no ``status.toml``, so there is no status to
    read.  It is not silent, though: ``imago.py``'s CLI is a thin
    wrapper over the same callable core a flight uses, so a hand
    run leaves the same two traces, and BOTH are required here.

    The lock is the stronger of the two, and it lives in the
    scratch we are about to remove.  It is taken before any work
    begins and released in the cleanup that always runs, so its
    presence means the run owns this directory now or died without
    releasing it.  That is also what makes a reclamation racing a
    just-started run refuse rather than mis-time it: the new run
    takes the lock before it writes anything at all.

    The marker is the run's own word that its driver reached
    cleanup.  Note carefully that ``runtime`` is opened in APPEND
    mode, so a directory run four times holds four markers and
    only the LAST non-blank line describes the current state -- a
    run interrupted after three good ones ends in something else
    entirely.  Reading the tail is therefore correct where
    searching the file would be badly wrong.

    One asymmetry with the workspace contract is worth knowing:
    the marker is written from a ``finally``, so it records that
    the driver reached cleanup, not that the calculation
    succeeded.  A run that failed but exited tidily is reclaimable
    here, where a workspace unit would have been preserved for the
    curator by having no ``result.toml``.  A job tree carries no
    success signal with which to close that gap; ``--older-than``
    is the lever for anyone who wants recent failures kept.
    """

    if os.path.exists(os.path.join(target, IMAGO_LOCK)):
        return False, "imagoLock present: running or died"

    log_path = os.path.join(run_dir, RUNTIME_LOG)
    if not os.path.isfile(log_path):
        return False, "no runtime log"
    try:
        with open(log_path, errors="replace") as handle:
            lines = [line.strip() for line in handle
                     if line.strip()]
    except OSError as exc:
        return False, f"unreadable runtime log ({exc})"

    if not lines:
        return False, "empty runtime log"
    if lines[-1] != COMPLETION_MARKER:
        return False, "runtime log ends mid-run"
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

    if not _is_under(target, os.path.realpath(scratch_root)):
        return None, f"target not under scratch root: {target}"
    return target, ""


def _is_under(path, root):
    """Return whether ``path`` lies inside ``root``.

    ``commonpath`` rather than ``startswith``: a plain prefix test
    would accept "/scratch-other" for a root of "/scratch", which
    on a tool that deletes is the difference between a containment
    check and a near miss.
    """

    try:
        return os.path.commonpath([path, root]) == root
    except ValueError:              # different drives / relative
        return False


def tree_stats(path):
    """Return ``(total_bytes, newest_mtime)`` for the tree under
    ``path``, following no symlinks.

    Symlinks are counted as their own (tiny) size rather than
    followed, so a stray link inside scratch cannot inflate the
    reported saving or, worse, lead the walk outside the tree.

    The newest mtime is gathered here, on the walk that has to
    visit every file anyway, because taking it from the scratch
    *directory* would be wrong: a directory's mtime moves only
    when entries are added to or removed from it, so a job that
    has spent a week writing into an already-created HDF5 still
    presents a week-old directory.  Anything that ages scratch
    must age it by its newest file.
    """

    total = 0
    newest = 0.0
    for directory, _, files in os.walk(path, followlinks=False):
        for name in files:
            full = os.path.join(directory, name)
            try:
                info = os.lstat(full)
            except OSError:
                continue            # vanished mid-walk; ignore
            total += info.st_size
            newest = max(newest, info.st_mtime)
    return total, newest


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


def find_job_run_dirs(root):
    """Yield the run directories of a JOB TREE, and the
    workspaces it declines to enter.

    Each item is ``(kind, path)`` with ``kind`` either ``"run"``
    or ``"workspace"``.  A hand-run directory is one carrying an
    ``intermediate`` symlink: there is no ``status.toml`` to key
    on, and no fixed depth either, since job trees are organized
    however their author liked.

    Two refusals shape the walk.  As in :func:`find_run_dirs` it
    never descends a symlink (refusal 4) -- ``intermediate`` is
    one, and following it would leave for the scratch area.  And
    it never descends into a workspace (refusal 7): a job tree may
    hold one far below it, and walking in would judge units that
    have a provable contract by the presumption-based one.  Those
    directories are yielded rather than quietly dropped, so the
    bytes deliberately left behind stay visible in the report.
    """

    for directory, subdirs, files in os.walk(root):
        if os.path.isdir(os.path.join(directory, WINGBEATS)):
            yield "workspace", directory
            subdirs[:] = []             # refusal 7: stay out
            continue
        subdirs[:] = [name for name in subdirs
                      if not os.path.islink(
                          os.path.join(directory, name))]
        if os.path.islink(os.path.join(directory, INTERMEDIATE)):
            yield "run", directory


def detect_root_kind(root):
    """Decide what KIND of root ``root`` is, by looking at it
    rather than by being told: ``WORKSPACE_ROOT``,
    ``JOB_TREE_ROOT``, or None when it is neither.

    The decision is made once, up front, and fixes the contract
    for the whole call -- one call handles one kind, so that a
    single report never gathers two different safety contracts
    under one set of totals (DESIGN 6.2.12).
    """

    if os.path.isdir(os.path.join(root, WINGBEATS)):
        return WORKSPACE_ROOT
    for kind, _ in find_job_run_dirs(root):
        if kind == "run":
            return JOB_TREE_ROOT
    return None


def plan_reclamation(root, scratch_root, policy=None,
                     ids=None, calc_pattern=None, older_than=None,
                     kind=None, match=None):
    """Build the reclamation plan for a root: one record per run
    directory, saying whether its scratch can go and why not when
    it cannot (PSEUDOCODE 13.8).

    Nothing is removed here.  The plan is the report, and applying
    it is a separate, explicit step -- which is what lets the same
    code path serve the preview, the standalone run, and the
    producer's ``--clean-after``.

    The root's kind fixes the contract for the whole call, and
    with it the policy that fits when none is supplied: a
    workspace is judged by ``status.toml``, a job tree by its lock
    and its ``runtime`` log.  A job tree also yields the
    workspaces it declined to enter; those become records of their
    own so the bytes left behind stay visible.

    :param root: the directory to tidy -- a workspace (holding
        ``wingbeats/``) or a job tree of ordinary run directories.
    :param scratch_root: the area scratch must lie under, normally
        ``$IMAGO_TEMP``.
    :param policy: ``(run_dir, target) -> (spent, reason)``;
        defaults to the one matching the root's kind.
    :param ids: workspace only -- restrict to these stable ids, or
        None for all.
    :param calc_pattern: workspace only -- a glob matched against
        the calc tag (the run directory's name), or None for all.
    :param older_than: only scratch whose newest file is at least
        this many days old, or None for any age.
    :param kind: the root's kind, when the caller has already
        detected it; None to detect it here.
    :param match: job tree only -- a glob matched against the run
        directory's path relative to the root, or None for all.
    """

    if kind is None:
        kind = detect_root_kind(root)
    if kind is None:
        # Neither kind of root, so there is nothing here to
        #   reclaim -- and no contract under which to try.  Fall
        #   through to the job-tree walk and it would judge by the
        #   wrong rules; return empty instead.
        return []
    if policy is None:
        policy = (default_reclaim_policy if kind == WORKSPACE_ROOT
                  else hand_run_policy)

    if kind == WORKSPACE_ROOT:
        base = os.path.abspath(os.path.join(root, WINGBEATS))
        found = [("run", run_dir) for run_dir in find_run_dirs(root)]
    else:
        base = os.path.abspath(root)
        found = list(find_job_run_dirs(root))

    now = datetime.now().timestamp()
    plan = []

    # Resolve EVERY run directory's scratch before judging any of
    #   it.  Scratch mirrors the run directory's path, so a run
    #   nested inside another has its scratch nested too, and the
    #   containment refusal below needs the whole set to compare
    #   against -- including runs this call filtered out, since an
    #   excluded inner run is exactly the one an outer removal
    #   would take as collateral.
    #
    # The link is resolved twice: once here, to build that
    #   comparison set, and once inside plan_one_dir when it judges
    #   the directory.  Those are two cheap stats, and paying them
    #   keeps plan_one_dir usable on its own, with no caller
    #   obliged to hand it a pre-resolved target.
    resolved = []                    # (directory, relative, target)
    for item_kind, directory in sorted(found, key=lambda p: p[1]):
        relative = os.path.relpath(os.path.abspath(directory), base)

        # A workspace declined by a job-tree walk is not a
        #   candidate at all; it is recorded so the report can say
        #   where the untouched bytes went (refusal 7).
        if item_kind == "workspace":
            plan.append(dict(run_dir=directory, unit=relative,
                             ok=False, workspace=True, bytes=0,
                             reason="workspace: tidy it directly"))
            continue
        target, _ = scratch_target(directory, scratch_root)
        resolved.append((directory, relative, target))

    every_target = [target for _, _, target in resolved
                    if target is not None]

    for directory, relative, _ in resolved:
        if not _selected(relative, kind, ids, calc_pattern, match):
            continue

        # Every refusal from here down belongs to plan_one_dir, so
        #   this whole-tree sweep and the in-flight prune of layer
        #   (b) apply one set of rules rather than two.
        #
        # ``every_target`` is the comparison set for the nesting
        #   refusal: every run the WALK found, not the selected
        #   subset, so a filtered-out inner run still stops the
        #   outer removal -- precisely the case that refusal
        #   exists to stop.
        plan.append(plan_one_dir(directory, relative, scratch_root,
                                 policy, every_target, older_than,
                                 now=now))
    return plan


def plan_one_dir(run_dir, label, scratch_root, policy,
                 other_targets=(), older_than=None, now=None):
    """Judge ONE run directory and return its plan record
    (PSEUDOCODE 13.8).

    The whole per-directory decision -- resolve the link, refuse
    what must be refused, apply the policy, measure the tree --
    lives here so that :func:`plan_reclamation`, which sweeps a
    whole root, and :func:`reclaim_one_dir`, which prunes a single
    unit mid-flight, reach it by the same path and cannot drift
    apart.

    :param run_dir: the run directory to judge.
    :param label: how the report should name this run -- the path
        relative to the walk's base for a whole-tree plan, the
        unit's own directory name for a single prune.  Nothing is
        decided from it.
    :param scratch_root: the area scratch must lie under.
    :param policy: ``(run_dir, target) -> (spent, reason)``.
    :param other_targets: the comparison set for the nesting
        refusal, supplied by the CALLER because containment is the
        one refusal a single directory cannot judge on its own: a
        whole-tree sweep passes every run its walk found, and an
        in-flight prune passes the flight's other units.
    :param older_than: only scratch whose newest file is at least
        this many days old, or None for any age.
    :param now: the timestamp the age is measured against, so a
        sweep dates every entry from one instant; None to read the
        clock here.
    """

    if now is None:
        now = datetime.now().timestamp()

    target, reason = scratch_target(run_dir, scratch_root)
    if target is None:
        return _skipped(run_dir, label, reason)

    # REFUSAL 5: an outer tree holding another run's scratch is
    #   deferred, never removed.  Taking it would delete the inner
    #   run's working files as collateral, which would make the
    #   "never touch an unfinished run" refusal a formality.  A
    #   second pass reclaims it once the inner ones are gone, so
    #   nothing is lost -- only deferred.
    nested = [other for other in other_targets
              if other != target and _is_under(other, target)]
    if nested:
        return _skipped(
            run_dir, label,
            f"holds {len(nested)} nested run's scratch; "
            f"reclaim those first")

    spent, reason = policy(run_dir, target)
    if not spent:
        return _skipped(run_dir, label, reason)

    size, newest = tree_stats(target)
    if older_than is not None:
        age_days = (now - newest) / SECONDS_PER_DAY
        if age_days < older_than:
            return _skipped(run_dir, label,
                            f"only {age_days:.1f} days old")

    return dict(run_dir=run_dir, unit=label, ok=True, reason="",
                target=target, bytes=size)


def reclaim_one_dir(run_dir, scratch_root, policy=None,
                    other_targets=(), older_than=None, label=None):
    """Judge ONE finished run and remove its scratch when the
    policy calls it spent (PSEUDOCODE 13.8; DESIGN 6.2.12 layer
    (b)).

    This is the entry point a client uses to prune as a campaign
    advances, rather than sweeping the workspace after it ends.
    Nothing about flights appears here: the caller decides *when*
    to call it -- as a unit lands -- and *which* policy applies.
    This is only the mechanism.

    Returns the plan record :func:`plan_one_dir` built, with two
    fields added when a removal was attempted: ``removed`` says
    whether the tree actually went, and ``failure`` carries the
    message when it did not.  A caller can then report a single
    prune the way the standalone tool reports a whole sweep, and
    can tell a refusal (``ok`` false) from a failure (``ok`` true,
    ``removed`` false) -- a distinction that matters, because the
    first is the mechanism working and the second means a
    filesystem assumption has broken.

    The default policy is the workspace one, since a unit landing
    in a flight is what this layer exists for.
    """

    if policy is None:
        policy = default_reclaim_policy
    if label is None:
        label = os.path.basename(os.path.normpath(run_dir))

    record = plan_one_dir(run_dir, label, scratch_root, policy,
                          other_targets, older_than)
    if not record["ok"]:
        return record

    removed, _, failures = apply_reclamation([record])
    record["removed"] = removed == 1
    record["failure"] = failures[0][1] if failures else None
    return record


def _skipped(directory, relative, reason):
    """Build the plan record for a run directory that will not be
    reclaimed, carrying the reason so the report can explain the
    skip rather than passing over it silently."""

    return dict(run_dir=directory, unit=relative, ok=False,
                reason=reason, bytes=0)


def _selected(relative, kind, ids, calc_pattern, match):
    """Return whether one run directory passes the CLI filters.

    A workspace selects on the concepts its layout supplies -- the
    stable id (the first path component under ``wingbeats/``) and
    the calc tag (the directory's own name).  A job tree has
    neither, so it selects on the whole path relative to the root,
    which is the only handle its free-form layout offers.
    """

    if kind == WORKSPACE_ROOT:
        if ids is not None and relative.split(os.sep)[0] not in ids:
            return False
        if calc_pattern is not None and not fnmatch.fnmatch(
                os.path.basename(relative), calc_pattern):
            return False
        return True
    return match is None or fnmatch.fnmatch(relative, match)


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

    Workspaces a job-tree walk declined to enter are the one
    exception, and are always named in full.  They are not a unit
    the tool judged and passed over; they are a whole subtree of
    bytes it deliberately left for a separate call, and the user
    cannot act on that without knowing where they are.
    """

    reclaimable = [item for item in plan if item["ok"]]
    declined = [item for item in plan if item.get("workspace")]
    skipped = [item for item in plan
               if not item["ok"] and not item.get("workspace")]

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

    if declined:
        print(f"  {'skipped':11} {len(declined):>9}  workspaces, "
              f"not descended:")
        for item in declined:
            print(f"      {item['unit']}")

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
        if declined:
            print(f"{len(declined)} workspace(s) left alone; tidy "
                  f"each by pointing at it directly")
        if reclaimable:
            print("re-run with --apply to remove them")


# ============================================================
#  CLI
# ============================================================

def build_parser():
    parser = argparse.ArgumentParser(
        description="Reclaim the intermediate scratch an imago run "
                    "leaves behind, in either a kaleidoscope "
                    "workspace or a tree of ordinary run "
                    "directories.  Removes only scratch, never a "
                    "run directory, never an unfinished run, and "
                    "never a link pointing outside the scratch "
                    "area.  Previews by default and removes "
                    "nothing until --apply.")
    parser.add_argument(
        "root",
        help="the directory to tidy.  Its kind is detected: one "
             "holding wingbeats/ is a flight workspace (for the "
             "producer, <data_root>/curation/workspace), and one "
             "holding ordinary run directories is a job tree.  A "
             "single call handles one kind or the other")
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
        help="workspace only: restrict to this stable id; "
             "repeatable (default: every id in the workspace)")
    parser.add_argument(
        "--calc", dest="calc_pattern", metavar="GLOB",
        help="workspace only: restrict to run directories whose "
             "calc tag matches this glob, e.g. 'kpt-mesh-*' or "
             "'loen-*' (default: every calc)")
    parser.add_argument(
        "--match", metavar="GLOB",
        help="job tree only: restrict to run directories whose "
             "path relative to the root matches this glob, e.g. "
             "'c/diamond/*' (default: every run directory)")
    parser.add_argument(
        "--older-than", type=float, metavar="DAYS",
        help="restrict to scratch whose newest file is at least "
             "this many days old (default: any age)")
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
    if not os.path.isdir(args.root):
        print(f"error: {args.root!r} is not a directory.",
              file=sys.stderr)
        return 2

    # One call handles one kind of root, so the kind is settled
    #   here and passed down rather than re-detected per unit.
    kind = detect_root_kind(args.root)
    if kind is None:
        print(f"error: {args.root!r} holds neither a {WINGBEATS}/ "
              f"directory nor any run directory with an "
              f"'{INTERMEDIATE}' link, so there is no imago "
              f"scratch here to reclaim.", file=sys.stderr)
        return 2

    # Reject a filter belonging to the OTHER kind of root rather
    #   than ignoring it.  On a tool that deletes, a filter which
    #   silently does nothing would quietly widen the run to
    #   everything the user meant to narrow it down from.
    if kind == JOB_TREE_ROOT:
        misplaced = [name for name, value
                     in (("--id", args.ids),
                         ("--calc", args.calc_pattern)) if value]
        instead = "--match"
    else:
        misplaced = ["--match"] if args.match else []
        instead = "--id and --calc"
    if misplaced:
        print(f"error: {' and '.join(misplaced)} selects on the "
              f"other kind of root, but {args.root!r} is a "
              f"{kind}.  Use {instead} instead.", file=sys.stderr)
        return 2

    print(f"root: {args.root} ({kind})")
    plan = plan_reclamation(
        args.root, args.scratch_root, kind=kind, ids=(
            set(args.ids) if args.ids else None),
        calc_pattern=args.calc_pattern, match=args.match,
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
