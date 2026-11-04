---
layout: default
title: insert3cbo.py
script: insert3cbo.py
description: 'Merges three-centre bond order data into the two-centre bond order raw file and adds fake centroid atoms to structure.dat so that three-centre bonds can be visualised with makeBOND. NOT CURRENTLY FUNCTIONAL: it edits structure.dat while makeBOND now reads imago.skl, and must be re-merged with makeBOND before it works again.'
parent: Script Catalog
---

# insert3cbo

> **This script is currently broken.** Its own header says so: it
> worked by modifying `structure.dat`, but `makeBOND.py` now takes
> its atom and lattice data from `imago.skl`, so the modified
> structure is never read. The intended fix is for this script to
> write a separate file that makeBOND reads and applies on top of
> the skeleton data. Until that is done, running it produces files
> that nothing consumes. The page records what it is meant to do.

Creates new bond order raw and structure files that include
information describing **three-centre bonds**, with the two-centre
bond information removed except for real two-centre bonds, which are
identified by having no three-centre bond that contains them. Three-
centre bonds arise in electron-deficient systems such as boron
hydrides and carboranes, where one bonding pair is shared among three
atoms, and the ordinary pairwise bond order cannot represent them.

The merge represents each three-centre bond by a **fake atom at its
centroid**, added to the structure, so that the bond can be drawn as
two-centre bonds from each member atom to the centroid in the
visualisation tools [makeBOND.py](makeBOND.html) feeds.

**Defaults.** Every option's default value is read from the resource
control file `insert3cborc.py` in the directory named by `$IMAGO_RC`.

## Usage

**Commandline:**

```
insert3cbo.py [-3cdata FILE] [-2cdata FILE] [-sdata FILE]
              [-rawout FILE] [-sout FILE] [--help]
```

**Option Table:**

| Option    | Defaults               | Description |
|:---------:|:----------------------:|:------------|
|-3cdata    |gs_bond-mb.3c.raw       | Three-centre bond order raw file from an Imago bond calculation. |
|-2cdata    |gs_bond-mb.raw          | Two-centre bond order raw file. |
|-sdata     |structure.dat           | Structure file to read. |
|-rawout    |gs_bond-mb.23cbo.raw    | Output: the two-centre raw file manipulated to carry both two- and three-centre bond orders. |
|-sout      |structure.23cbo.dat     | Output: the structure file with the fake centroid atoms added. |
|--help     |N/A                     | Print usage information and exit. |
