---
layout: default
title: skl2pdb.py
script: skl2pdb.py
description: 'Converts an imago.skl skeleton file into a Protein Data Bank (.pdb) file with the atoms sorted by element, and writes a map file recording which skeleton atom became which PDB atom. Useful for viewing an Imago structure in molecular visualisation software.'
parent: Script Catalog
---

# skl2pdb

Converts an Imago skeleton `.skl` file into a Protein Data Bank
`.pdb` file, with the atoms rearranged so that they are sorted
according to their elements. A **map file** is also written recording
the correspondence between the skeleton atoms and the PDB atoms. This
is the reverse of [pdb2skl.py](pdb2skl.html), and the common use is
to open an Imago structure in a molecular viewer that reads PDB.

The PDB file always carries a `CRYST1` record with the cell
parameters, so periodic viewers can draw the cell.

**Defaults.** Every option's default value is read from the resource
control file `skl2pdbrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
skl2pdb.py [-i INPUT] [-o OUTPUT] [-m MAPFILE] [--skltypes] [--help]
```

**Variable Table:**

| Variable   | Accepted Values     | Description |
|:----------:|:-------------------:|:------------|
|INPUT       |String; file path    | Path to the imago skeleton input file. |
|OUTPUT      |String; file path    | Path to the PDB output file. |
|MAPFILE     |String; file path    | Path to the skeleton-to-PDB atom map file. |

**Option Table:**

| Option         | Outputs        | Defaults       | Description |
|:--------------:|:--------------:|:--------------:|:------------|
|-i / --input    |N/A             |imago.skl       | Imago skeleton input file. |
|-o / --output   |PDB file        |imago.pdb       | PDB output file. |
|-m / --map      |Map file        |sklPDB.map      | File recording which skeleton atom became which PDB atom. Same format as the map file described on the [pdb2skl.py](pdb2skl.html) page. |
|--skltypes      |N/A             |Off             | Use the species types defined in the skeleton file as atom types in the PDB output. The default is to assign no types. |
|--help          |N/A             |N/A             | Print usage information and exit. |
