---
layout: default
title: vasp2skl.py
script: vasp2skl.py
description: 'Converts a VASP position file (CONTCAR or POSCAR) into an imago.skl skeleton file, taking the element names from the POTCAR file. Supports both VASP 4.6 and VASP 5.x formats and skips selective-dynamics lines.'
parent: Script Catalog
---

# vasp2skl

Reads a VASP atomic position file and produces an `imago.skl` skeleton
file. The usual case is importing a structure that VASP has relaxed,
which is why the default input is `CONTCAR` and the default output is
`imago.skl.relaxed`.

A VASP position file contains the lattice vectors in Cartesian x, y, z
form and the atomic positions in fractional (direct) coordinates. This
script extracts those data, computes the conventional cell parameters
(a, b, c, alpha, beta, gamma), and writes them in skeleton format.

**Element names.** VASP 4.6 position files do not name the elements,
so the names are obtained from the `POTCAR` file by scanning for the
`TITEL` lines, which carry the element symbol of each pseudopotential.
VASP 5.x adds an extra line of element symbols between the lattice
vectors and the atom counts; the script detects that line and handles
it automatically, but the atoms file is still read.

**Selective dynamics.** If the line after the atom counts is not
`direct` (case-insensitive) it is taken to be a selective-dynamics
flag line, skipped, and the next line is expected to be `direct`.

**Defaults.** Every option's default value is read from the resource
control file `vasp2sklrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
vasp2skl.py [-i INPUT] [-a ATOMFILE] [-o OUTPUT] [--help]
```

**Variable Table:**

| Variable   | Accepted Values     | Description |
|:----------:|:-------------------:|:------------|
|INPUT       |String; file path    | VASP position file, CONTCAR or POSCAR. |
|ATOMFILE    |String; file path    | File containing the element names in the order they appear in the position file, normally POTCAR. |
|OUTPUT      |String; file path    | Path to the imago skeleton output file. |

**Option Table:**

| Option         | Outputs        | Defaults            | Description |
|:--------------:|:--------------:|:-------------------:|:------------|
|-i / --input    |N/A             |CONTCAR              | Name of the VASP atomic position file after relaxation. If CONTCAR is not found the script falls back to POSCAR. |
|-a / --atoms    |N/A             |POTCAR               | Name of the file that contains the element names as they appear in the position file. |
|-o / --output   |Skeleton file   |imago.skl.relaxed    | Name of the output imago skeleton file. |
|--help          |N/A             |N/A                  | Print usage information and exit. |
