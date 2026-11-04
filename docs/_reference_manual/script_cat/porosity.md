---
layout: default
title: porosity.py
script: porosity.py
description: 'Calculates the porosity of a solid from an imago.skl file by sampling the cell on a regular mesh and marking every point within a scaled covalent radius of an atom as solid. Porosity is the fraction of the cell volume left unmarked.'
parent: Script Catalog
---

# porosity

Calculates the porosity of a solid-state material. The unit cell is
sampled on a regular mesh, and each mesh point is classified as
*solid* if it lies within the (scaled) covalent radius of any atom,
or *void* otherwise. The porosity is the ratio of the void volume to
the total cell volume. This is a geometric measure that depends on
the radii you choose, so it is best used to compare models built the
same way rather than as an absolute number.

**Defaults.** Every option's default value is read from the resource
control file `porosityrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
porosity.py (-m NUMA NUMB NUMC | -r RESOLUTION) [-l LIMITDIST]
            [-i INFILE] [-sf SCALEFACTOR] [--help]
```

One of `-m` or `-r` must be given; they are mutually exclusive.

**Variable Table:**

| Variable         | Accepted Values     | Description |
|:----------------:|:-------------------:|:------------|
|NUMA, NUMB, NUMC  |integer              | Number of sampling points along the a, b, c lattice directions. |
|RESOLUTION        |float; Angstroms     | Target spacing between sampling points. |
|LIMITDIST         |float; Angstroms     | Radius within which an atom can contribute to the pore map. |
|INFILE            |String; file path    | Imago skeleton file. |
|SCALEFACTOR       |float                | Multiplier applied to the covalent radii from the element database. |

**Option Table:**

| Option   | Outputs      | Defaults       | Description |
|:--------:|:------------:|:--------------:|:------------|
|-m        |N/A           |None            | Explicit mesh counts along a, b, c. Mutually exclusive with `-r`. |
|-r        |N/A           |None            | Target resolution. The mesh counts are computed from it and the cell dimensions, then the resolution is refined so that an integer number of points fits each axis exactly. Mutually exclusive with `-m`. |
|-l        |N/A           |4.0             | Interaction distance limit. |
|-i        |N/A           |imago.skl       | Input skeleton file. |
|-sf       |N/A           |1.0             | Scale factor on the covalent radii. The scaled radius defines the pore boundary: mesh points within it of any atom centre are solid. |
|--help    |N/A           |N/A             | Print usage information and exit. |

## How It Works

1. Read the lattice vectors and atom positions from the skeleton.
2. Either compute the per-axis resolution from the mesh counts
   (`-m`), or compute the mesh counts from the target resolution
   (`-r`).
3. Scale the covalent radii by the scale factor and set the
   interaction distance limit.
4. Set up the 3D mesh of sampling points.
5. For every atom, mark every mesh point within its scaled covalent
   radius as solid.
6. Count the remaining unmarked points.
7. Report the porosity as pore volume divided by cell volume.

## Output

The script prints its progress and then two lines to standard
output:

```
total pore volume = <value> A^3
porosity = <fraction>
```

Increase the mesh density until the porosity stops changing; a
coarse mesh overestimates the solid fraction near sphere surfaces.
