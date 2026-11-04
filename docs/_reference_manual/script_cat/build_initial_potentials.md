---
layout: default
title: build_initial_potentials.py
script: build_initial_potentials.py
description: 'The producer of the augmented initial-potential database. Reads a curation manifest naming reference solids, fetches and converts each structure, runs or reuses a cached Imago SCF for it through kaleidoscope, harvests the converged potential at each distinct atomic environment, and writes the per-element s_gaussian_pot.toml files that makeinput.py reads to start an SCF from a solid-state potential instead of the free atom.'
parent: Script Catalog
---

# build_initial_potentials

Builds the **augmented initial-potential database**: the per-element
files `$IMAGO_DATA/atomicPDB/<element>/s_gaussian_pot.toml` that
`makeinput.py` reads when it starts a calculation. The default
potential database holds only free-atom potentials. This script adds
to each element's file one entry per distinct chemical environment
found in a curated set of **reference solids**, so that a new
calculation on, say, silicon in an oxide can start from a potential
already converged for silicon in an oxide and reach self-consistency
in fewer iterations.

This is the *producer* half of a library / producer / consumer
split. The consumer is `makeinput.py`; the shared library is
`initial_potential_db.py`.

**The curation manifest.** The pipeline's primary input is one
human-readable TOML file, the curation manifest. It declares which
solids, which atom sites, under which labels, and with which SCF
settings. It has three jobs:

1. **Declare the curation set.** The manifest *is* the curation
   strategy made explicit. Adding a reference solid means adding a
   manifest entry; reviewing the set means reading the manifest.
   Nothing is hidden inside the script.
2. **Tell the pipeline what to harvest.** For each reference solid,
   which atom sites' converged potentials enter the database and
   under what labels. The harvest auto-discovers one representative
   per distinct environment; a manifest entry exists only to
   override that.
3. **Record the SCF settings used.** k-points, convergence
   threshold, basis and functional are copied into each database
   entry's provenance, so every entry carries the conditions of its
   reference run.

The database must be regeneratable from the curated set, not a
hand-edited artifact, so regeneration is a deterministic function of
(manifest, structure files, Imago build). The manifest is the build
configuration, the structure files are the source, the Imago build is
the toolchain, this script is the build script, and the database is
the compiled output.

**Where the manifest comes from.** For a set of crystals from the
Crystallography Open Database, `cod_fish.py pin` writes a complete
manifest directly. For local structures, non-default settings or
per-structure customization, write a sketch and finish it with
[expand_manifest.py](expand_manifest.html). Non-crystalline
structures have their atoms grouped into species by
[makegroups.py](makegroups.html) as part of the build.

## Usage

**Commandline:**

```
build_initial_potentials.py [--manifest PATH] [--pdb-root DIR]
        [--element EL] [--force] [--verbose]
        [--dispatch {local,slurm-pooled,slurm-per-job}]
        [--partition NAME] [--nodes N] [--walltime HH:MM:SS]
        [--profile NAME] [--save-config]
        [--clean-after] [--tidy-run]
        [--materialize-only] [--materialize-dir DIR]
        [--submit] [--orchestrator-cores N]
        [--orchestrator-memory SIZE] [--orchestrator-walltime HH:MM:SS]
        [--help]
```

**Typical sequence on a cluster:**

```
cod_fish.py search --elements Si O
cod_fish.py pin 1 3 --save-cif cifs/ > manifest.toml
build_initial_potentials.py --materialize-only        # check the structures convert
build_initial_potentials.py --submit --verbose        # run as a batch job
```

**Option Table:**

| Option                  | Defaults                                     | Description |
|:-----------------------:|:--------------------------------------------:|:------------|
|--manifest               |manifest.toml in the current directory        | Path to the curation manifest, schema version 2. |
|--pdb-root               |$IMAGO_DATA/atomicPDB                         | The potential database root to write into. |
|--element                |Every element named in the manifest           | Restrict the build to one element's database. |
|--force                  |Reuse cached runs                             | Bypass the run-reuse cache so every reference run re-executes. |
|--verbose                |Print only failures, warnings and closing tallies | Narrate each step: every reference solid as it is fetched and converted, and every scratch tree as it is pruned or kept. |
|--dispatch               |The cluster settings file's default topology  | Where to run the flight. `local` runs every unit in process and needs no cluster settings file. `slurm-pooled` streams units through one shared allocation. `slurm-per-job` submits one scheduler job per unit. |
|--partition              |First partition in the cluster settings file  | Scheduler queue for a cluster dispatch. |
|--nodes                  |The cluster settings file's value             | Nodes per allocation. |
|--walltime               |The cluster settings file's value             | Time limit for a cluster dispatch, e.g. `02:00:00`. |
|--profile                |The settings file's base settings             | Named cluster profile to select from the settings file. |
|--save-config            |Off                                           | Write the resolved cluster dispatch choices beside the run for a reproducible record. |
|--clean-after            |Leave scratch in place                        | Once the harvest is written, reclaim the intermediate scratch of every finished unit, typically most of the bytes a run leaves behind. |
|--tidy-run               |Off                                           | Reclaim each unit's scratch as that unit finishes rather than at the end, for a campaign large enough to fill the scratch filesystem before it completes. Combines with `--clean-after`. |
|--materialize-only       |Run the full build                            | Fetch and convert every reference structure named in the manifest, then stop without running any SCF. A pre-flight to validate a freshly pinned structure set. |
|--materialize-dir        |The shared structure cache beside the databases | Directory the pre-flight writes fetched CIFs and converted skeletons into. |
|--submit                 |Run the build in this process                 | Materialize the structures on the login node, then submit the producer as its own batch job that runs the build under `--dispatch`. The driver job's resources come from the `orchestrator` block of `clusterrc.py`. |
|--orchestrator-cores     |The clusterrc orchestrator block              | Cores for the `--submit` driver job. |
|--orchestrator-memory    |The clusterrc orchestrator block              | Memory for the driver job, as a scheduler size such as `8G`. |
|--orchestrator-walltime  |The clusterrc orchestrator block, else `--walltime` | Time limit for the driver job. It must outlast the flight it supervises. |
|--help                   |N/A                                           | Print usage information and exit. |

## What a Build Does

1. **Refresh the isolated baseline.** For every element the manifest
   names, load or create its database file and rebuild its
   `"isolated"` entry from the current `pot1` and `coeff1` files.
   This alone forms a working producer-to-consumer loop with no
   SCF; it is what guarantees every element file has a default
   entry even when nothing else is curated for it.
2. **Materialize the structures.** Fetch each pinned COD structure
   with `cod_fish.py get`, verifying its revision, or read the local
   structure path, and convert it to a skeleton. Non-crystalline
   structures are grouped into species by `makegroups.py`.
3. **Run the reference calculations** as a kaleidoscope flight,
   reusing any run the content-keyed cache already holds. The
   k-point density is predicted from the historical guidance
   database and verified by a convergence sweep.
4. **Harvest.** For each converged solid, take the SCF potential at
   one representative atom of each distinct environment, attach the
   fingerprint that identifies that environment, and write the
   entry into the element's `s_gaussian_pot.toml` with its
   provenance.

## Reproducibility

Three layers, with different guarantees. The TOML emitter writes
byte-identical bytes for a fixed in-memory database. The numerical
outputs agree, for the same manifest, input files and Imago build, to
the precision the SCF and fit chain can reach; bit-identity is not
promised, since accumulation order, threading and library versions
can shift the last bits. Timestamps and commit hashes in the
provenance refresh every run and carry no reproducibility guarantee.
