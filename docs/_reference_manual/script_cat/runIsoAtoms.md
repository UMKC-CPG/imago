---
layout: default
title: runIsoAtoms.py
script: runIsoAtoms.py
description: 'Runs a full Imago SCF calculation for every element as an isolated atom in a 20 Angstrom box and collects the converged potentials into an atomPDB directory, giving an isolated-atom potential database computed with the solid-state code itself rather than with atomSCF.'
parent: Script Catalog
---

# runIsoAtoms

Runs an isolated-atom calculation for every atom in the element
database, using Imago itself rather than the atomic SCF program, and
collects the converged potentials into an atomic potential function
database. The potentials in `$IMAGO_DATA/atomicPDB` come from a
spherical atomic SCF; the ones this script produces come from the
same periodic, Gaussian-expanded machinery that a solid-state run
uses, so they are the appropriate "isolated" reference when you want
to compare a potential in a solid against the free atom on equal
footing (as `makeFittedRhoV.py` does with its `coeff.isolated` files).

For each element, or for the single element given with `-e`, the
script:

1. Creates a directory named after the element.
2. Writes a minimal `imago.skl` with a single atom in a
   20 x 20 x 20 Angstrom cell.
3. Runs `makeinput.py` to generate the Imago input files.
4. Runs `imago.py` to perform the SCF calculation.

After all elements are computed it:

5. Creates an `atomPDB` directory.
6. For each element, copies the converged SCF potential
   coefficients from `gs_scfV-fb.dat` into `coeff1` and
   `coeff.isolated` files.
7. Copies the existing potential definition file `pot1` from
   `$IMAGO_DATA/atomicPDB/`.

The result has the same layout as `atomicPDB` and can be given to
`makeinput.py -potdb`.

**Defaults.** Every option's default value is read from the resource
control file `runIsoAtomsrc.py` in the directory named by
`$IMAGO_RC`.

## Usage

**Commandline:**

```
runIsoAtoms.py [-e Z] [--nocore] [--help]
```

Run it in an empty working directory; it creates one subdirectory
per element there, plus `atomPDB`.

**Variable Table:**

| Variable   | Accepted Values       | Description |
|:----------:|:---------------------:|:------------|
|Z           |integer; atomic number | The single element to process. 0 means every element. |

**Option Table:**

| Option          | Outputs                     | Defaults | Description |
|:---------------:|:---------------------------:|:--------:|:------------|
|-e / --element   |N/A                          |0 (all)   | Atomic Z number of a single element to process. |
|--nocore         |N/A                          |Off       | Pass `-nocore` to makeinput so that the all-electron electronic structure is computed. See below. |
|--help           |N/A                          |N/A       | Print usage information and exit. |

## Option Details

### All-electron Calculations

> --nocore

It is somewhat counter-intuitive that an option called `-nocore`
yields an all-electron result. The idea is that `-nocore` *will not*
include any core orbitals in the list of orbitals that are
orthogonalized away. Therefore all of the core orbitals *will be*
included in the list of orbitals to compute, and the SCF treats every
electron explicitly.
