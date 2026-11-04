---
layout: default
title: makePDOS.py
script: makePDOS.py
description: 'Partial density of states post-processor. Reads the raw PDOS file from an Imago DOS calculation and combines its per-atom or per-type data units into the curves you ask for, either by a control file or by internal defaults, producing a PDOS.plot file for plotgraph.py.'
parent: Script Catalog
---

# makePDOS

Reads the raw partial density of states (PDOS) file written by an
Imago DOS calculation, collects its data components, and combines
them according to your requests into a plottable file. The raw file
holds one data unit per atom or per atomic type; this script is how
you turn hundreds of such units into a handful of meaningful curves,
such as "all silicon s states" or "the total DOS of oxygen species 2".
The output `PDOS.plot` is plotted with
[plotgraph.py](plotgraph.html).

**Defaults.** Every option's default value is read from the resource
control file `makePDOSrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
makePDOS.py [-i CONTROL_FILE] [-o OUTPUT_FILE] [-f RAW_DATA_FILE]
            [-xanes] [-negToZero] [--help]
```

**Variable Table:**

| Variable        | Accepted Values     | Description |
|:---------------:|:-------------------:|:------------|
|CONTROL_FILE     |String; file path    | A control file describing which data units and columns to combine into each output curve. See the control file section. |
|OUTPUT_FILE      |String; file path    | Plottable output file. |
|RAW_DATA_FILE    |String; file path    | Raw PDOS file from the DOS calculation. |

**Option Table:**

| Option      | Outputs      | Defaults            | Description |
|:-----------:|:------------:|:-------------------:|:------------|
|-i           |N/A           |None; internal defaults are used | Control file defining how to collect and organize the raw data. |
|-f           |N/A           |gs_dos-fb.p.raw      | Raw PDOS data file. |
|-o           |PDOS.plot     |PDOS.plot            | Output file. |
|-xanes       |N/A           |Off                  | Change the internally generated default control so the curves are decomposed into s+d and p states. Only meaningful when no `-i` control file is given. |
|-negToZero   |N/A           |Off                  | Replace any negative value in the final output with zero. Useful for cleaning up small numerical artifacts. |
|--help       |N/A           |N/A                  | Print usage information and exit. |

## The Raw PDOS File

To write a control file you first need to understand the raw input
file. It has a header followed by a series of data units.

1. **The header.** The first line is the number of data units in the
   file. The second line is the number of data points in the energy
   scale for each unit. The rest of the header is the energy scale
   itself, typically a few thousand points, one per line.

2. **The data units.** Each unit has a header and a block of data.

   * **The unit header** is a series of filter name / filter value
     pairs. The first filter name is always `SEQUENCE_NUM` and the
     last is always `COL_LABELS`. The first value is the atom number
     or the type number, depending on whether the calculation
     printed the PDOS per atom or per atomic type. The last value is
     a space-separated string naming the data columns printed for
     this unit. Per atom there is only one column, `TOTAL`. Per type
     there are many, for example `TOTAL 2s 3s 2p 3p 3d`. The filter
     names in between, such as `ELEMENT_NAME`, `SPECIES_ID` and
     `TYPE_ID`, identify the unit and are what a control file
     matches against.
   * **The unit data** is the columns named in `COL_LABELS`, one
     value per point of the energy scale. Per atom that is a single
     column holding the atom's total PDOS. Per type there is one
     column for the total and one for each valence orbital of the
     type.

## Default Behaviour (No Control File)

Without `-i` the script builds a control internally. Without
`-xanes` it produces one output column per data unit, each holding
the `TOTAL` PDOS of that unit, labelled from the unit's element,
species and type. With `-xanes` it produces two columns per unit, one
for the s+d orbitals and one for the p orbitals, a decomposition
commonly wanted for XANES analysis.

## The Control File

> -i

The control file is a plain text file of colon-separated fields.

1. **The filter definition line** is the first line. It is a
   colon-separated list of tags. The first tag is always `LABEL` and
   the last is always `COL_LABELS`. The tags in between can be any
   filter name present in the raw file's unit headers, such as
   `SEQUENCE_NUM`, `ELEMENT_NAME`, `SPECIES_ID` or `TYPE_ID`.

2. **The curve definition lines** follow, one per output column and
   therefore one per plotted curve. Each is a colon-separated list of
   values in the order the definition line set. A field may hold a
   *set* of values separated by spaces, or the word `all`, so that
   several data units can match one curve and be summed into it.

**Examples.** With the filter definition line

```
LABEL : ELEMENT_NAME : SPECIES_ID : TYPE_ID : COL_LABELS
```

each remaining line gives a label, which elements contribute, which
species of those elements, which types of those species, and finally
which data columns to include:

| Curve definition                        | Meaning |
|:----------------------------------------|:--------|
|`TOTAL : all : all : all : TOTAL`        | The total DOS. Every element, species and type matches, and only the total column of each is included. |
|`Si2_2 : si : 2 : 2 : TOTAL`             | The total DOS of silicon atoms of species 2 and type 2. |
|`Si_s : si : all : all : s`              | The combined s orbitals of all silicon. |
|`Si_2_3_4 : si : 2 3 4 : all : TOTAL`    | The combined DOS of silicon species 2, 3 and 4. |

Filter values are matched exactly against the raw file, which
writes element names in lowercase, so write `si`, not `Si`. The
label in the first field is free text and may use any case.

A complete control file for a silicate might therefore read:

```
LABEL : ELEMENT_NAME : SPECIES_ID : TYPE_ID : COL_LABELS
TOTAL : all : all : all : TOTAL
Si_s  : si  : all : all : s
Si_p  : si  : all : all : p
O_s   : o   : all : all : s
O_p   : o   : all : all : p
```

The output file has the energy scale in column one and one column
per curve definition, with the labels as column headers, ready for
`plotgraph.py -i PDOS.plot -ft dos`.
