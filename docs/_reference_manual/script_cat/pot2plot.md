---
layout: default
title: pot2plot.py
script: pot2plot.py
description: 'Converts an Imago SCF potential file (scfV) and its structure.dat into per-element plottable files: one listing the Gaussian potential coefficients of every potential type of the element, and one evaluating each potential function on a radial grid.'
parent: Script Catalog
---

# pot2plot

Converts an Imago self-consistent potential file (`scfV`) and the
associated `structure.dat` into a set of plottable files, one pair per
element. The SCF potential in Imago is a sum of Gaussians for each
potential type, so this script lets you inspect how the fitted
potential differs from one type of an element to another, or watch
how it changes between SCF runs.

**Defaults.** Every option's default value is read from the resource
control file `pot2plotrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
pot2plot.py [-i SCFV_FILE] [-s STRUCT_FILE] [-d] [--help]
```

**Variable Table:**

| Variable     | Accepted Values     | Description |
|:------------:|:-------------------:|:------------|
|SCFV_FILE     |String; file path    | Imago SCF potential file. |
|STRUCT_FILE   |String; file path    | Imago structure file, used to learn which element each potential type belongs to. |

**Option Table:**

| Option              | Outputs      | Defaults          | Description |
|:-------------------:|:------------:|:-----------------:|:------------|
|-i / --infile        |N/A           |gs_scfV-fb.dat     | Input scfV file. |
|-s / --struct        |N/A           |structure.dat      | Input structure file. |
|-d / --down_spin     |N/A           |Off                | Use the second half of the potential file, which holds the spin-down potential of a spin-polarized calculation. |
|--help               |N/A           |N/A                | Print usage information and exit. |

## Output

For every element in the structure two files are written.

**`pot_coeffs_<element>.plot`** holds the raw Gaussian coefficients.
Each row is one potential type of that element, numbered from 1 in
the first column, and each following column is one potential function
coefficient, in the order of the exponential alphas in the scfV
file.

**`pot_eval_<element>.plot`** holds the potential functions evaluated
on a radial grid. The first column is the radius, 1000 points from 0
to 5 in the radial units of the potential file (atomic units), and
each following column is one potential type of that
element, giving the sum of its Gaussian terms at that radius. This
file is directly plottable with
`plotgraph.py -i pot_eval_<element>.plot`.
