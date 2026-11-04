---
layout: default
title: makePotDB.py
script: makePotDB.py
description: 'Builds the default atomic potential database in $IMAGO_DATA/atomicPDB. For every element it writes the potential definition file pot1 and fits the numerical free-atom potential from makeBDB with a linear Gaussian expansion to produce the initial SCF potential coefficients in coeff1.'
parent: Script Catalog
---

# makePotDB

Creates the default **atomic potential database**,
`$IMAGO_DATA/atomicPDB`, using information about each element stored
in the `element_data.py` module. It also creates an initial set of
potential coefficients for each element, based on the numerical
free-atom potential that `atomSCF` produced when the basis set
database was built with [makeBDB.py](makeBDB.html). Imago reads this
database when `makeinput.py` prepares a job: the `pot1` file fixes the
Gaussian exponents used to expand each element's potential, and
`coeff1` supplies the starting coefficients for the SCF cycle.

For each element the script does two things.

1. **Writes the potential definition file `pot1`** from the element
   data: the nuclear charge, a default alpha, a default covalent
   radius, and the number and range of Gaussian exponents.
2. **Fits the potential.** It reads `numericalPot.dat` from the
   element's `atomicBDB` directory, adjusts the header to this
   element's potential definition, runs `gaussFit` to perform a
   *linear* Gaussian fit onto the fixed exponents (quicker than the
   non-linear fit in makeBDB), appends three zero columns (total
   charge, valence charge and spin polarization) to each
   coefficient row, and saves the result as `coeff1`.

An element whose `pot1` or `coeff1` already exists is skipped. The
distributed database that `unpackImagoDB.py` unpacks is the output
of this script.

**Defaults.** Every option's default value is read from the resource
control file `makePotDBrc.py` in the directory named by `$IMAGO_RC`.

## Usage

**Commandline:**

```
makePotDB.py [-e Z] [-f NUMFORKS] [--help]
```

**Variable Table:**

| Variable   | Accepted Values       | Description |
|:----------:|:---------------------:|:------------|
|Z           |integer; atomic number | The single element to process. 0 means every element. |
|NUMFORKS    |integer                | Number of parallel processes for the Gaussian fitting. |

**Option Table:**

| Option          | Outputs                     | Defaults | Description |
|:---------------:|:---------------------------:|:--------:|:------------|
|-e / --element   |N/A                          |0 (all)   | Atomic Z number of a single element to process. |
|-f / --fork      |N/A                          |1         | Number of parallel processes for the Gaussian fitting. |
|--help           |N/A                          |N/A       | Print usage information and exit. |

## Files Written

Both land in `$IMAGO_DATA/atomicPDB/<element>/`.

**`pot1`**, the potential definition:

```
NUCLEAR_CHARGE__ALPHA
<Z> 20.000000
COVALENT_RADIUS
1.000000
NUM_ALPHAS
<number of Gaussian exponents>
ALPHAS
<minimum alpha> <maximum alpha>
```

**`coeff1`**, the fitted coefficients: the first line is the number
of terms, then one line per Gaussian with the coefficient followed
by the three zero columns. This is the same format `makeinput.py`
expects when you point it at an alternative database with `-potdb`.

Requires `$IMAGO_DATA` to be set and `makeBDB.py` to have completed
its SCF step for every element you ask for.
