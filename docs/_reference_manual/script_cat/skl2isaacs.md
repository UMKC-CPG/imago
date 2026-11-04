---
layout: default
title: skl2isaacs.py
script: skl2isaacs.py
description: 'Converts an imago.skl skeleton file into the two input files of ISAACS, the Interactive Structure Analysis of Amorphous and Crystalline Structures tool: an XML project file (.ipf) and a Chem3D coordinate file (.chem3d).'
parent: Script Catalog
---

# skl2isaacs

Converts an Imago skeleton `.skl` file into a set of input files for
ISAACS (Interactive Structure Analysis of Amorphous and Crystalline
Structures), a tool for analysing the structural properties of both
crystalline and disordered materials. ISAACS computes quantities such
as pair distribution functions, coordination statistics and ring
statistics, which makes it a useful companion for the amorphous
models Imago is often applied to.

ISAACS uses a two-file scheme. The `.ipf` file is an XML project
file holding the chemistry, the box, the periodic boundary setting
and the cutoff information. The `.chem3d` file holds the Cartesian
atom coordinates. This script writes both from one skeleton.

**Defaults.** Every option's default value is read from the resource
control file `skl2isaacsrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
skl2isaacs.py [-i INPUT] [-o OUTPUTROOT] [--help]
```

**Variable Table:**

| Variable     | Accepted Values     | Description |
|:------------:|:-------------------:|:------------|
|INPUT         |String; file path    | Imago skeleton input file. |
|OUTPUTROOT    |String               | Root name for the two output files. |

**Option Table:**

| Option         | Outputs                              | Defaults    | Description |
|:--------------:|:------------------------------------:|:-----------:|:------------|
|-i / --input    |N/A                                   |imago.skl    | Input imago skeleton file. |
|-o / --output   |`<root>.ipf` and `<root>.chem3d`      |imago        | Root name for the output ISAACS files. |
|--help          |N/A                                   |N/A          | Print usage information and exit. |
