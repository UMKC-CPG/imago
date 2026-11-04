---
layout: default
title: guidance_promote.py
script: guidance_promote.py
description: 'The curator''s helper for the historical guidance database. Reviews the entries guidance_harvest.py staged and moves the accepted ones into entries/, where the k-point predictor reads them, either interactively, by an objective rule, or wholesale. Catches re-runs of an already-promoted solid and retires the duplicate.'
parent: Script Catalog
---

# guidance_promote

The curator's promotion helper for the historical guidance database.
[guidance_harvest.py](guidance_harvest.html) stages one entry per
converged sweep into `staging/<system_type>/`. Staging is a deliberate
checkpoint: an automated harvest is not the same as scientific
endorsement, and a buggy harvest, a sweep that converged at a
numerical artifact, or a mis-classified structure could all stage an
entry that should not feed future predictions. This script is where a
human curator, or an objective rule standing in for one, decides
which staged entries graduate into `entries/<system_type>/`, where
the predictor reads them.

**Promotion is a pure move.** The file's bytes, and therefore its
provenance, never change as it crosses the staging boundary, so a
promoted entry is byte-identical to what the curator reviewed.

## Usage

**Commandline:**

```
guidance_promote.py [--db-root DIR] [--auto-promote | --all | --dry-run]
                    [--help]
```

**Option Table:**

| Option           | Defaults                          | Description |
|:----------------:|:---------------------------------:|:------------|
|--db-root         |$IMAGO_DATA/historicalGuidanceDB   | The database root. Required explicitly if `$IMAGO_DATA` is unset. |
|(no mode flag)    |Interactive                        | Print each staged entry's summary and ask the curator to PROMOTE, SKIP, DELETE or, for a duplicate, REPLACE. |
|--auto-promote    |Off                                | Promote every staged entry that passes the objective acceptance test below, leaving the rest in staging for review. In practice this clears most of a seed flight so the curator reviews only the outliers. |
|--all             |Off                                | Promote every staged entry without checking the rule, for when the curator has already looked over the directory. |
|--dry-run         |Off                                | Report what the chosen mode would do, moving nothing. |
|--help            |N/A                                | Print usage information and exit. |

The three mode flags are mutually exclusive.

## The Objective Acceptance Test

> --auto-promote

An entry passes when all three conditions hold, evaluated from the
staged file alone:

1. **Converged in the middle of the grid.** A sweep that converged at
   either endpoint is suspect, since the grid may have been too
   narrow to bracket the true converged point, so the converged
   k-density must sit within the central 60 percent of the grid's
   span.
2. **A convincingly flat top of grid.** The spread (maximum minus
   minimum) of the top three grid points' energies, per atom and in
   eV, must be below ten times the metric threshold. This asks for a
   converged *region*, not one consecutive difference that happened
   to dip below threshold.
3. **Gap fields consistent.** The gap kind is `none` exactly when the
   gap energy is zero.

An entry with no verification block, or one that omits the grid
energies (a hand-written manual entry), cannot pass the rule and is
left for interactive review.

## Duplicate Detection

Promotion is the only stage that sees both an incoming entry and the
promoted corpus, so it is where a re-run of an already-promoted solid
is caught. The harvest cannot do this: it writes one entry per
converged solid with no view of what a curator accepted months ago,
and a re-run mints a fresh entry id by construction.

This matters because the predictor treats every stored entry as an
independent observation. With five nearest neighbours, five copies of
one calculation fill the whole neighbour set, the weighted variance
collapses to zero, and a single measurement is delivered as
near-certainty, which then drives the flight builder into its
narrowest search.

Two entries make the same claim when they describe the same structure
under the same settings, and they agree when their converged meshes
match. Agreement retires the staged copy to `superseded/`;
disagreement, or a mesh that cannot be compared, is reported and left
in staging for the curator. Every mode applies this check, `--all`
included: refusing to store one claim twice is a correctness guard,
not a quality judgement, and `--all` waives only the latter.
Interactive review additionally offers REPLACE, which is the only way
an entry ever leaves `entries/`.
