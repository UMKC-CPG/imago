---
layout: default
title: xyz2skl.py
script: xyz2skl.py
description: 'Converts a plain XYZ Cartesian coordinate file into an imago.skl skeleton file with the atoms sorted by element, building an enclosing orthorhombic P1 cell, and writes a map file recording which XYZ atom became which skeleton atom.'
parent: Script Catalog
---

# xyz2skl

Converts a file of pure XYZ direct-space coordinates into an
`imago.skl` skeleton file. The atoms are rearranged so that they are
sorted according to their elements, and a **map file** is written that
records the correspondence between each atom in the original XYZ file
and its skeleton atom.

**Cell handling.** The XYZ format carries no lattice information, so
the script builds an orthorhombic P1 box that encloses all the atoms.
Imago always works with periodic boundary conditions, so a molecule
still needs a box. The box is assigned space group `1_a` (P1) so that
no symmetry operations are ever applied to it.

**Defaults.** Every option's default value is read from the resource
control file `xyz2sklrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
xyz2skl.py [-i INPUT] [-o OUTPUT] [-m MAPFILE] [--xyz-types] [--help]
```

**Variable Table:**

| Variable   | Accepted Values     | Description |
|:----------:|:-------------------:|:------------|
|INPUT       |String; file path    | Path to the XYZ input file. |
|OUTPUT      |String; file path    | Path to the imago skeleton output file. |
|MAPFILE     |String; file path    | Path to the XYZ-to-skeleton atom map file. |

**Option Table:**

| Option         | Outputs        | Defaults       | Description |
|:--------------:|:--------------:|:--------------:|:------------|
|-i / --input    |N/A             |model.xyz       | XYZ input file. |
|-o / --output   |Skeleton file   |imago.skl       | Imago skeleton output file. |
|-m / --map      |Map file        |sklXYZ.map      | File recording which XYZ atom became which skeleton atom. Same format as the map file described on the [pdb2skl.py](pdb2skl.html) page. |
|--xyz-types     |N/A             |Off             | Use the atom types as defined in the XYZ file as species tags in the skeleton. The default is to assign no types, so all atoms of the same element share a single species. |
|--help          |N/A             |N/A             | Print usage information and exit. |

## Input Format

The XYZ file has the usual layout:

```
<num_atoms>
<title or comment line>
<Element>  <x>  <y>  <z>
...
```

The first line is the atom count. The second line is a free-form
comment; its first word is stored as the skeleton's title. Each
following line is a one- or two-letter element symbol
(case-insensitive) followed by three Cartesian coordinates in
Angstroms.
