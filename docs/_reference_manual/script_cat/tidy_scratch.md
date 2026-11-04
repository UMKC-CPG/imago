---
layout: default
title: tidy_scratch.py
script: tidy_scratch.py
description: 'Reclaims the intermediate scratch behind finished Imago runs, which is almost all of the bytes a run leaves, while never removing a run directory, an unfinished run, or anything a symlink points to outside $IMAGO_TEMP. Handles a kaleidoscope flight workspace or a tree of ordinary job directories. Previews by default; deletes only with --apply.'
parent: Script Catalog
---

# tidy_scratch

Reclaims the intermediate scratch that an Imago run leaves behind. A
finished run directory spans two filesystems, and the difference
between them is almost entirely one of size. The **kept** tier is
the run directory itself: the staged inputs, `result.toml`,
`status.toml`, `cache_key.toml`, the SCF potential and the log, a
couple of hundred kilobytes. The **scratch** tier sits behind the
`intermediate` symlink that `imago.py` creates, pointing into
`$IMAGO_TEMP`, and holds the engine's working files, above all the
HDF5 carrying the wavefunctions. Measured on a seed-scale producer
run, scratch was 99.7 percent of the bytes, about 25 MB per
calculation.

**Dry run is the default.** The tool previews what it would remove,
with sizes, and removes nothing until `--apply` is given.

## Two Kinds of Root

Every ordinary `imago.py` run plants the same `intermediate` link, so
the tool recognizes two kinds of root and decides which it has by
looking:

* a **workspace** holds `wingbeats/`; its units are kaleidoscope
  units, and each proves it is finished in `status.toml`;
* a **job tree** holds ordinary run directories with no flight above
  them, such as a student's job directory or a hand-driven
  convergence test. These prove they are finished by holding no
  `imagoLock` and ending their `runtime` log with the completion
  marker.

One call handles exactly one kind.

Reclaiming a workspace unit is safe because the harvest reads only
the paths recorded in `result.toml`, every one of which resolves
inside the run directory, and the run-reuse cache decides a hit from
`status.toml` and `cache_key.toml` alone. A reclaimed run still
harvests and still counts as a cache hit.

## Usage

**Commandline:**

```
tidy_scratch.py ROOT [--apply] [--scratch-root DIR]
                [--id ID]... [--calc GLOB] [--match GLOB]
                [--older-than DAYS] [--verbose] [--help]
```

**Examples:**

```
tidy_scratch.py share/curation/workspace                   # preview
tidy_scratch.py share/curation/workspace --apply           # reclaim
tidy_scratch.py <root> --id si_cmce_64_1999 --calc 'kpt-mesh-*' --apply
tidy_scratch.py ~/imago/jobs --match 'c/diamond/*' --apply # a job tree
tidy_scratch.py <root> --older-than 7 --apply
```

**Option Table:**

| Option          | Defaults                 | Description |
|:---------------:|:------------------------:|:------------|
|ROOT             |None; required            | The directory to tidy. One holding `wingbeats/` is a flight workspace; one holding ordinary run directories is a job tree. |
|--apply          |Preview only              | Actually remove the scratch. |
|--scratch-root   |$IMAGO_TEMP               | The area scratch must lie under. A resolved target outside it is refused. |
|--id             |Every id in the workspace | Workspace only: restrict to this stable id. Repeatable. |
|--calc           |Every calc                | Workspace only: restrict to run directories whose calc tag matches this glob, e.g. `kpt-mesh-*` or `loen-*`. |
|--match          |Every run directory       | Job tree only: restrict to run directories whose path relative to the root matches this glob, e.g. `c/diamond/*`. |
|--older-than     |Any age                   | Restrict to scratch whose newest file is at least this many days old. |
|--verbose        |Summarize by reason       | List every skipped unit individually. |
|--help           |N/A                       | Print usage information and exit. |

## What It Refuses to Do

The tool deletes, so what it will *not* do matters more than what it
will. Five refusals hold for every root:

1. **Only scratch is removed, never the run directory.** The kept
   tier is the record of the calculation. The `intermediate` link is
   deliberately left dangling, so a reclaimed run still shows where
   its scratch was.
2. **An unfinished unit is never touched.** A `running` or `queued`
   status means the engine may still be writing, and a missing
   `result.toml` means the run produced nothing to keep, usually the
   state a curator most wants to investigate.
3. **A link out of the scratch area is never followed.** The target
   is resolved and checked to lie under the scratch root before
   anything is removed. A cleanup tool a symlink can redirect is a
   hazard.
4. **A symlink is never descended while walking.** The walk skips
   every symlinked subdirectory, keeping the candidate run
   directories inside the root by construction.
5. **A tree holding another run's scratch is deferred.** Scratch
   mirrors the run directory's path, so a run nested inside another
   has its scratch nested too. The outer tree is skipped while any
   other run's scratch lies within it; a second pass reclaims it
   once those are gone.

Two more hold for a job tree, where no `status.toml` exists:

6. **A run that has not declared it finished is never reclaimed.**
   Its scratch must hold no `imagoLock` *and* its `runtime` log must
   end with the completion marker. No log at all, or one ending
   mid-stream, is reported and left alone.
7. **A workspace is never descended into from a job tree.** Such
   directories are named in the report and pruned from the walk.

## Known Blind Spot

Scratch whose run directory has already been deleted is invisible to
this tool, since it only ever starts from run directories. Such
orphans under `$IMAGO_TEMP` must be found and removed by hand.
