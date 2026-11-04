---
layout: default
title: makeFittedRhoV.py
script: makeFittedRhoV.py
description: 'Evaluates the fitted SCF charge density and potential of an Imago calculation on a 3D real-space mesh, plus their differences from isolated-atom values, by building an input file for the imagoRhoV Fortran program and running it. Produces profile data along each lattice axis and openDX volume files.'
parent: Script Catalog
---

# makeFittedRhoV

Evaluates the charge density and potential functions from a completed
Imago SCF calculation on a 3D real-space mesh so that they can be
visualised. Data are produced for the potential, the total charge
density, the valence charge density, the potential difference and the
valence charge difference. The differences are taken with respect to
the potential and charge of a single isolated atom of each element,
as stored in the potential database, so they show how the solid has
rearranged charge relative to free atoms.

**It is important to understand that the charge density plotted here
is derived from the fitted charge.** Imago fits the density to a sum
of Gaussians during the SCF cycle, and this script evaluates that
fit, not the density of the wavefunctions directly.

The script itself does no numerical work. It builds an input file of
Gaussian coefficients and then runs the `imagoRhoV` Fortran program,
which evaluates the functions on the mesh and writes the outputs.

**Defaults.** Every option's default value is read from the resource
control file `makeFittedRhoVrc.py` in the directory named by
`$IMAGO_RC`. That file also sets the isolated-atom database path,
by default `$IMAGO_DATA/atomicPDB`, and the name of the temporary
input file, `tempPotRho.dat`.

## Usage

**Commandline:**

```
makeFittedRhoV.py -m A B C [-i INPUT_FILE] [-spindn] [-noDX]
                  [-nc NUM_CELLS] [--help]
```

Run it in a job directory containing `structure.dat` and the SCF
potential file.

**Variable Table:**

| Variable    | Accepted Values     | Description |
|:-----------:|:-------------------:|:------------|
|A, B, C      |integer, each > 0    | Number of mesh points along each of the a, b, c axes of the (possibly non-orthogonal) cell. Periodic boundary conditions are observed. |
|INPUT_FILE   |String; file path    | SCF charge and potential data file. |
|NUM_CELLS    |0, 1 or 2            | Levels of neighbouring replicated cells included when evaluating the functions. |

**Option Table:**

| Option              | Outputs      | Defaults          | Description |
|:-------------------:|:------------:|:-----------------:|:------------|
|-m / --mesh          |N/A           |None; required     | Number of mesh points along a, b, c. The script stops if this is not given. |
|-i / --input         |N/A           |gs_scfV-fb.dat     | Source of the SCF charge and potential data. |
|-spindn / --spindn   |N/A           |Off                | Use the spin-down section of the potential file. See below. |
|-noDX / --noDX       |N/A           |Off (openDX data is produced) | Do not create the openDX plottable data. |
|-nc / --num-cells    |N/A           |1                  | Levels of neighbouring cells: 0 is the original cell only, 1 is first neighbours (27 cells), 2 is two levels (125 cells). |
|--help               |N/A           |N/A                | Print usage information and exit. |

## Option Details

### Spin-polarized Calculations

> -spindn

The potential file has two sections. For a non-spin-polarized
calculation they are identical, so the second is ignored as
duplicate data. For a spin-polarized calculation they differ, and
`-spindn` selects the second section as the data source. To see both
spin up and spin down, run the script twice, once with and once
without this flag.

### openDX Data

> -noDX

By default the openDX volume data are produced along with the
profiles, because they cost little extra time and it saves running
the computation twice. Give `-noDX` to skip them when only the
profiles are wanted.

## What the Script Does

1. Reads the SCF potential file and `structure.dat`.
2. Reads the isolated-atom coefficient file, `coeff.isolated`, for
   each element from the potential database.
3. Writes the temporary input file for `imagoRhoV`, which has these
   columns for every Gaussian term:

   | Column | Contents |
   |:------:|:---------|
   | 1      | Exponential alpha of the Gaussian. |
   | 2      | Potential coefficient from the SCF calculation. |
   | 3      | Difference between the SCF and isolated-atom potential coefficients. |
   | 4      | Valence charge coefficient from the SCF calculation. |
   | 5      | Difference between the SCF and isolated-atom valence charge coefficients. |
   | 6      | Total charge coefficient from the SCF calculation. |

   There is no column for the difference in total charge at this
   time, although one could be added if the need arose.
4. Runs `imagoRhoV` with the mesh dimensions, the number of
   neighbouring cell levels, and the openDX flag as its arguments.

## Output

`imagoRhoV` writes the following files to the current directory.

| File                              | Contents |
|:----------------------------------|:---------|
|profile_a.dat, profile_b.dat, profile_c.dat | Each function profiled along one lattice axis: the potential is averaged over the perpendicular plane, and the charge densities are accumulated over it and corrected for the cross-sectional area. Always produced. |
|pot.dx                             | The SCF potential on the mesh. |
|potDelta.dx                        | Potential minus the isolated-atom potential. |
|rhoV.dx                            | The valence charge density. |
|rhoVDelta.dx                       | Valence charge density minus the isolated-atom valence density. |
|rhoT.dx                            | The total charge density. |
|lattice.dx                         | The cell outline for openDX. |

The `.dx` files are omitted with `-noDX`.
