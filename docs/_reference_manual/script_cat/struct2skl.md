---
layout: default
title: struct2skl.py
script: struct2skl.py
description: 'Converts an Imago structure.dat file back into an imago.skl skeleton file with the atoms sorted by element, and writes a map file recording which structure.dat atom became which skeleton atom.'
parent: Script Catalog
---

# struct2skl

Converts a file in the format of an Imago `structure.dat` file back
into an `imago.skl` skeleton file. The atoms are rearranged so that
they are sorted according to their elements, and a **map file** is
written that records the correspondence between each atom in the
structure file and its skeleton atom.

A `structure.dat` file is what `makeinput.py` produces from a
skeleton: the fully expanded cell with every atom listed explicitly.
This script is the reverse direction. It is useful when a structure
has been modified at the expanded level, or when a `structure.dat`
is the only form of a structure you have, and you want to start a
fresh job from it.

**Defaults.** Every option's default value is read from the resource
control file `struct2sklrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
struct2skl.py [-i INPUT] [-o OUTPUT] [-m MAPFILE] [--structtypes]
              [--makemol] [--help]
```

**Variable Table:**

| Variable   | Accepted Values     | Description |
|:----------:|:-------------------:|:------------|
|INPUT       |String; file path    | Path to the structure.dat input file. |
|OUTPUT      |String; file path    | Path to the imago skeleton output file. |
|MAPFILE     |String; file path    | Path to the structure-to-skeleton atom map file. |

**Option Table:**

| Option         | Outputs        | Defaults         | Description |
|:--------------:|:--------------:|:----------------:|:------------|
|-i / --input    |N/A             |structure.dat     | Structure input file. |
|-o / --output   |Skeleton file   |imago.skl         | Imago skeleton output file. |
|-m / --map      |Map file        |sklStruct.map     | File recording which structure atom became which skeleton atom. Same format as the map file described on the [pdb2skl.py](pdb2skl.html) page. |
|--structtypes   |N/A             |Off               | Use the species types from the structure.dat file in the skeleton output. The default is to assign no types, so every atom has type number 1. |
|--makemol       |N/A             |Off               | Export the structure as a cluster (molecule) rather than inside a periodic cell. **NOT YET FUNCTIONAL**: the option is accepted for future use but has no effect. |
|--help          |N/A             |N/A               | Print usage information and exit. |
