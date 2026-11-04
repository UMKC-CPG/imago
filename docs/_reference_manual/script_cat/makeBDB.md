---
layout: default
title: makeBDB.py
script: makeBDB.py
description: 'Builds the default atomic basis set database in $IMAGO_DATA/atomicBDB. For every element it runs an atomic SCF calculation, fits the numerical free-atom potential with Gaussians, and writes the two contraction input files (core-orthogonalized and no-core) from which makeinput.py later contracts the basis.'
parent: Script Catalog
---

# makeBDB

Performs all the operations needed to construct the default atomic
basis set for Imago calculations, using the atomic data stored in
the `element_data.py` module. The result is the **atomic basis
database**, `$IMAGO_DATA/atomicBDB`, with one subdirectory per
element. There are three steps.

1. **SCF.** An atomic SCF input file is created and the `atomSCF`
   program is run on it to obtain a numerical free-atom potential
   function.
2. **Fit.** The numerical potential is fit with Gaussians as
   accurately as possible by the `gaussFit` program. This is a
   non-linear fit and may take considerable time, which is why the
   fitting can be run in parallel with `--fork`.
3. **Contraction files.** Two input files for contracting a Gaussian
   basis set are written for each element. The first sequesters
   specific orbitals to the core, where they can be orthogonalized
   out of the secular (eigenvalue) equation; the second includes
   those core orbitals in the valence section so that they are
   included in the secular equation.

**The contraction is not actually performed here.** The wave
functions are contracted at the time `makeinput.py` is invoked. This
lets you modify the basis set for any purpose, but the basis sets
are applicable to an exceptionally wide variety of systems and it is
recommended that they not be modified without cause to "improve" a
result for a specific system.

Each step records that it is done by the files it leaves, so a
re-run skips elements that are already complete unless `--force` is
given. The distributed database that `unpackImagoDB.py` unpacks is
the output of this script; you only need to run it to rebuild the
database from scratch or to add an element.

**Defaults.** Every option's default value is read from the resource
control file `makeBDBrc.py` in the directory named by `$IMAGO_RC`.

## Usage

**Commandline:**

```
makeBDB.py [-s {all,scf,fit,cont}] [-e Z] [-f NUMFORKS] [--force]
           [-v] [--help]
```

**Variable Table:**

| Variable   | Accepted Values     | Description |
|:----------:|:-------------------:|:------------|
|Z           |integer; atomic number | The single element to process. 0 means every element. |
|NUMFORKS    |integer              | Number of parallel processes for the Gaussian fitting. |

**Option Table:**

| Option          | Outputs                   | Defaults | Description |
|:---------------:|:-------------------------:|:--------:|:------------|
|-s / --step      |See below                  |all       | Which step to perform: `all`, `scf` only, `fit` only, or `cont` (contraction files) only. |
|-e / --element   |N/A                        |0 (all)   | Atomic Z number of a single element to process. |
|-f / --fork      |N/A                        |1         | Number of parallel processes for the Gaussian fitting. |
|--force          |N/A                        |Off       | Overwrite any previously obtained results for the requested operations. |
|-v / --verbose   |N/A                        |Off       | Print progress reports. |
|--help           |N/A                        |N/A       | Print usage information and exit. |

## Files Written

All files land in `$IMAGO_DATA/atomicBDB/<element>/`, where
`<element>` is the lowercase symbol.

| Step  | Files | Contents |
|:-----:|:------|:---------|
|scf    |`atomSCF.dat`, `atomSCF.out`, `numericalPot.dat` | The atomSCF input, its log, and the numerical free-atom potential on a radial grid. The presence of `atomSCF.out` marks the step done. |
|fit    |`gaussFit.out`, `gauss.fit`, `progress.fit` | The fit log, the Gaussian coefficients and exponents, and the fit's progress record. The presence of `gaussFit.out` marks the step done. |
|cont   |`contract1_<element>`, `nocore_contract1_<element>` | The two contraction input files: core orbitals orthogonalized away, and all orbitals in the valence set. |

`makePotDB.py` reads `numericalPot.dat` from this database, so run
makeBDB's SCF step before building the potential database.
