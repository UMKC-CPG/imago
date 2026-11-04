---
layout: default
title: processPOPTC.py
script: processPOPTC.py
description: 'Partial optical properties post-processor. Splits the raw partial epsilon-2 file from an Imago POPTC calculation into its partial pairs, runs the Kramers-Kronig program imagoKKc on each, and gathers the resulting partial epsilon-1, ELF, refractive index, reflectivity and absorption into raw files that makePDOS.py can read.'
parent: Script Catalog
---

# processPOPTC

Reads the raw partial epsilon-2 (imaginary dielectric function) file
written by an Imago partial optical properties (POPTC) calculation,
takes each partial contribution in turn, runs the Kramers-Kronig
conversion program `imagoKKc` on it, and appends the results into new
raw files. The purpose is to allow calculation of the *partial*
energy-loss function (ELF), epsilon-1, imaginary epsilon-1 and
refractive index, which Imago itself only produces as totals.

The optical calculation writes one raw file containing every partial
epsilon-2 pair. For each pair this script creates a temporary input
file, runs `imagoKKc` on it to obtain epsilon-1, the ELF, the
refractive index, extinction coefficient, reflectivity and absorption
coefficient, and appends those curves into raw output files, one per
spectral quantity. When it finishes, each output file holds all
partial data for its quantity in the same raw format as a PDOS file,
so [makePDOS.py](makePDOS.html) can combine the partials into plottable
curves with a control file. Internally the script also uses makePDOS,
with a temporary control file, to extract each pair's epsilon-2 from
the raw input before handing it to `imagoKKc`.

**Defaults.** The spin channel default is read from the resource
control file `processPOPTCrc.py` in the directory named by
`$IMAGO_RC`.

## Usage

**Commandline:**

```
processPOPTC.py [-s {1,2}] [--help]
```

Run it in the job directory of a completed POPTC calculation. It
takes no file arguments; the files it reads and writes are fixed
Fortran unit numbers, listed below.

**Option Table:**

| Option      | Outputs      | Defaults  | Description |
|:-----------:|:------------:|:---------:|:------------|
|-s / --spin  |N/A           |1          | Spin channel to process: 1 is spin up, or the total for a non-spin-polarized calculation; 2 is spin down. |
|--help       |N/A           |N/A        | Print usage information and exit. |

## File Numbering Convention

**Input files**

| File       | Contents |
|:----------:|:---------|
|fort.250    | Raw epsilon-2 data, spin 1 (up, or total). |
|fort.251    | Raw epsilon-2 data, spin 2 (down). |
|fort.209    | KKC control data (scaling factors). |

**Intermediate files** (created and deleted for each partial pair)

| File       | Contents |
|:----------:|:---------|
|fort.700    | Control file for makePDOS. |
|fort.450    | Temporary epsilon-2 (input to imagoKKc). |
|fort.500    | Temporary epsilon-1, epsilon-2 and ELF (output from imagoKKc). |
|fort.510    | Temporary epsilon-1: total, x, y, z. |
|fort.520    | Temporary ELF: total, x, y, z. |
|fort.530    | Temporary refractive index: total, x, y, z. |
|fort.540    | Temporary extinction coefficient: total, x, y, z. |
|fort.550    | Temporary imaginary epsilon-1: total, x, y, z. |
|fort.560    | Temporary reflectivity: total, x, y, z. |
|fort.570    | Temporary absorption coefficient: total, x, y, z. |

**Output files** (spin 1 / spin 2)

| File             | Contents |
|:----------------:|:---------|
|fort.300 / 301    | Raw partial epsilon-1, epsilon-2 and ELF. |
|fort.310 / 311    | Raw partial epsilon-1: total, x, y, z. |
|fort.320 / 321    | Raw partial ELF: total, x, y, z. |
|fort.330 / 331    | Raw partial refractive index: total, x, y, z. |
|fort.340 / 341    | Raw partial extinction coefficient: total, x, y, z. |
|fort.350 / 351    | Raw partial imaginary epsilon-1: total, x, y, z. |
|fort.360 / 361    | Raw partial reflectivity: total, x, y, z. |
|fort.370 / 371    | Raw partial absorption coefficient: total, x, y, z. |

Each output file is a raw file in the makePDOS format. To plot, for
example, the partial ELF contributions of each element, run
`makePDOS.py -f fort.320 -i <control file> -o pELF.plot` with a
control file that groups the partial pairs as you wish, then plot the
result with [plotgraph.py](plotgraph.html).
