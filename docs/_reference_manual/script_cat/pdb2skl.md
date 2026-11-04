---
layout: default
title: pdb2skl.py
script: pdb2skl.py
description: 'Converts a Protein Data Bank (.pdb) file into an imago.skl skeleton file with the atoms sorted by element, and writes a map file recording which PDB atom became which skeleton atom.'
parent: Script Catalog
---

# pdb2skl

Converts a Protein Data Bank `.pdb` file into an `imago.skl` skeleton
file. The atoms are rearranged so that they are sorted according to
their elements, which is the order Imago expects. Because that
reordering loses the original atom sequence, the script also writes a
**map file** that records the correspondence between each PDB atom and
its skeleton atom, so later analysis can be traced back to residues
and molecules in the source structure.

**Cell handling.** If the PDB file carries a `CRYST1` record its cell
parameters are used directly. If it does not, the script builds an
orthorhombic P1 box that encloses all the atoms, padded on every side
by the `-b` buffer. Imago always works with periodic boundary
conditions, so a molecule still needs a box; the buffer is what keeps
periodic images from interacting.

**Defaults.** Every option's default value is read from the resource
control file `pdb2sklrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
pdb2skl.py [-i INPUT] [-o OUTPUT] [-m MAPFILE] [-b BUFFER]
           [--pdbtypes] [--help]
```

**Variable Table:**

| Variable   | Accepted Values     | Description |
|:----------:|:-------------------:|:------------|
|INPUT       |String; file path    | Path to the PDB input file. |
|OUTPUT      |String; file path    | Path to the imago skeleton output file. |
|MAPFILE     |String; file path    | Path to the PDB-to-skeleton atom map file. |
|BUFFER      |float; Angstroms     | Padding added around the molecule on every side when the PDB file defines no cell. |

**Option Table:**

| Option         | Outputs        | Defaults       | Description |
|:--------------:|:--------------:|:--------------:|:------------|
|-i / --input    |N/A             |model.pdb       | PDB input file. |
|-o / --output   |Skeleton file   |imago.skl       | Imago skeleton output file. |
|-m / --map      |Map file        |skl2PDB.map     | File recording which PDB atom became which skeleton atom. See the map file section. |
|-b / --buffer   |N/A             |0.0             | Buffer size in Angstroms around the molecule when no cell is defined. |
|--pdbtypes      |N/A             |Off             | Use the atom types defined in the PDB file as species tags in the skeleton. The default is to assign no types, so all atoms of the same element share a single species. |
|--help          |N/A             |N/A             | Print usage information and exit. |

## Option Details

### Species Tags

> --pdbtypes

In a skeleton file, atoms of the same element can be split into
different *species* so that Imago can give them different basis sets
or potentials and report their properties separately. By default this
script assigns no species: every carbon is species `c1`, every oxygen
`o1`, and so on. With `--pdbtypes` the atom names from the PDB file
are used to distinguish species, so for example a carbonyl carbon and
a backbone alpha carbon become different species. Use this when you
intend to analyse the system per site rather than per element.

### The Map File

> -m / --map

The map file has one fixed header line followed by one line per atom
in the skeleton:

```
OATOM#   MOL#   RES#   RESNAME   ELEMENT   TYPE
     1       1       1       ALA          c      c1
     2       1       1       ALA          c      c1
```

| Column   | Meaning |
|:--------:|:--------|
|OATOM#    | Imago atom index, 1 through the number of atoms. |
|MOL#      | Molecule sequence number from the PDB input. |
|RES#      | Residue sequence number from the PDB input. |
|RESNAME   | Residue name from the PDB input. |
|ELEMENT   | Lowercase element symbol. |
|TYPE      | Imago species tag used in the skeleton file. |

Other Imago scripts read this file to translate between skeleton atom
numbers and the residues and molecules of the original structure.
