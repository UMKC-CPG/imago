---
layout: default
title: graspElems.py
script: graspElems.py
description: 'Runs the Grasp2K relativistic atomic-structure pipeline (rnucleus, csl, rangular, rwfnestimate, rmcdhf, readrwf) for chosen elements and basis sets, either directly or by generating run scripts and SBATCH job files for a cluster. Produces the rwfn.out and isodata files that interFit.py fits.'
parent: Script Catalog
---

# graspElems

Runs the Grasp2K atomic-structure pipeline for specified elements
and basis sets. Grasp2K solves the multi-configuration Dirac-Hartree-
Fock equations for a free atom; its output is the set of relativistic
orbitals that [interFit.py](interFit.html) fits with Gaussians to
build Imago's relativistic basis. This script drives the six Grasp2K
programs in order, feeding each the answers it would otherwise
prompt for interactively:

1. `rnucleus`, the nuclear model;
2. `csl`, the configuration state list;
3. `rangular`, the angular coefficients;
4. `rwfnestimate`, the initial radial wavefunctions;
5. `rmcdhf`, the self-consistent Dirac-Hartree-Fock solve;
6. `readrwf`, which converts the binary orbitals to the readable
   `rwfn.out`.

By default the programs are executed directly. With `--batch` the
script instead writes `runC` and `runRSCF` shell scripts and SBATCH
job files, grouping elements across the requested number of cores,
for submission to a cluster.

**Requirements.** The Grasp2K binaries must be in the directory
named by `$GRASP_BIN`, or under `$GRASP_DIR/bin`; the script stops
if neither is set.

Author: Patrick Ryan Thomas.

**Defaults.** The batch flag and core count defaults are read from
`graspElemsrc.py` in the directory named by `$IMAGO_RC`.

## Usage

**Commandline:**

```
graspElems.py [-b] [-nc N] SPEC [SPEC ...]
```

**Examples:**

```
graspElems.py C MB                 # carbon, minimal basis
graspElems.py total FB             # every element, full basis
graspElems.py F MB H EB            # two element/basis pairs
graspElems.py C all                # carbon with MB, FB and EB
graspElems.py --batch C MB         # write scripts instead of running
graspElems.py --batch -nc 16 total # SBATCH files, 16 elements per job
```

**Variable Table:**

| SPEC form      | Meaning |
|:---------------|:--------|
|`ELEM BASIS`    | Run one element with one basis: MB, FB or EB. |
|`ELEM all`      | Run one element with all three bases. |
|`total`         | Run all elements with all bases. |
|`total BASIS`   | Run all elements with one basis. |

**Option Table:**

| Option              | Defaults | Description |
|:-------------------:|:--------:|:------------|
|-b / --batch         |Off       | Generate SBATCH scripts instead of running directly. |
|-nc / --num-cores    |1         | Number of cores for SBATCH job grouping. Batch mode only. |
|--help               |N/A       | Print usage information and exit. |

## Output

Each calculation runs in `<BASIS>/<element>/`, created under the
current directory, and leaves the Grasp2K files there, in
particular `isodata` and `rwfn.out`, which interFit reads. In batch
mode the current directory gains `SBATCHhead`, and numbered
`JobSBATCH-<n>` and `JobRSCF-<n>` job files, while each element
directory gains a `runC` script for the first four steps and a
`runRSCF` script for the self-consistent solve and conversion.
