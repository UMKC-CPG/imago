---
layout: default
title: guidance_harvest.py
script: guidance_harvest.py
description: 'Reads a finished kaleidoscope k-point convergence flight off disk, finds the converged k-point density for each structure by a two-sided energy test, and stages one historical-guidance entry per structure under $IMAGO_DATA/historicalGuidanceDB/staging for a curator to promote with guidance_promote.py.'
parent: Script Catalog
---

# guidance_harvest

Turns a finished kaleidoscope flight into staged **historical
guidance** entries. Imago predicts the k-point density a new
structure needs from a database of past converged calculations, the
historical guidance database. This script is the producer half of
that loop: after kaleidoscope has run a structure across a small grid
of k-point densities, it reads the finished workspace, finds the
converged grid point for each structure, and writes a rich entry to
`staging/<system_type>/` so that a curator can later promote it with
[guidance_promote.py](guidance_promote.html). The promoted entries
are what the predictor learns from, so every converged sweep sharpens
the next prediction.

**Where each entry field comes from.** Three inputs, each with one
clear job:

* **`flight.toml`**, the plan: the unit list, the sweep block saying
  which axis varied, and the per-structure prediction blocks the
  flight builder stashed. Each block carries the structure's
  prediction and the basis, functional and k-point integration it
  ran under. The swept k-density of each grid point is read out of
  its calc tag.
* **Each run's `result.toml`**, the per-run facts: the final SCF
  total energy, the measured band gap and its kind, the total
  magnetization, and the SCF threshold the run used. `imago.py`
  writes all of these.
* **The structure `.skl`**, the structural facts: the harvest loads
  it to compute the signature (composition and lattice family) and
  the same load yields the atom count and cell volume.

**The convergence rule** is two-sided. A grid point counts as
converged only when *both* its neighbours' total energies are within
the metric threshold of it, per atom and in eV, so that a single
numerical fluke cannot masquerade as a flat region. A sweep whose
energy is still moving at the top of the grid earns no entry: it is
logged, the flight is tagged as a prediction mismatch, and the
structure is skipped. Widen the grid and re-run.

The metric threshold is the solid's k-point convergence threshold
from its manifest or the default, and is distinct from the run's own
SCF threshold. The grid energies are stored raw, in total-cell
hartree, and every comparison normalizes per atom at the point of
use.

## Usage

**Commandline:**

```
guidance_harvest.py WORKSPACE_ROOT [--db-root DIR] [--help]
```

**Variable Table:**

| Variable         | Accepted Values         | Description |
|:----------------:|:-----------------------:|:------------|
|WORKSPACE_ROOT    |String; directory path   | The flight workspace directory, holding `flight.toml` and the `wingbeats/` run tree. |
|DIR               |String; directory path   | The historical guidance database root. |

**Option Table:**

| Option      | Defaults                             | Description |
|:-----------:|:------------------------------------:|:------------|
|--db-root    |$IMAGO_DATA/historicalGuidanceDB      | The database root to stage into. Required explicitly if `$IMAGO_DATA` is unset. |
|--help       |N/A                                   | Print usage information and exit. |

## Output

One TOML entry file per converged structure under
`<db-root>/staging/<system_type>/`, plus a printed summary of each
structure's outcome. Nothing is written under `entries/`; that is the
curator's step. Entries are written by the shared `guidance_db`
library's byte-deterministic emitter, so a harvest re-run over the
same workspace produces identical files apart from the timestamp
and build identity in the provenance.
