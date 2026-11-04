---
layout: default
title: bond_analysis.py
script: bond_analysis.py
description: 'Geometric bond analysis of an atomic structure from its imago.skl, with atoms counted as bonded when closer than the sum of their covalent radii times a factor. One mode per run: ball-and-stick or VTK model files, bond angle or bond length lists, bond orientational order, openDX model, bond statistics with optional histogram, coordination numbers, Q^n distribution, ring distribution, or OpenSCAD files for 3D printing.'
parent: Script Catalog
---

# bond_analysis

Takes an atomic structure file and produces a new file describing the
nearest-neighbour relationships between its atoms. This can be
considered a bond analysis, but the data come purely from geometry,
not from any quantum mechanical calculation such as Imago's bond
order. Two atoms are "bonded" when their separation is at most the
sum of their covalent radii times a bonding factor, 1.1 by default.
Use it to characterise a model before or without a full calculation:
coordination numbers, bond-length and bond-angle statistics, the
Q^n distribution of a glass, ring statistics, or a quick picture.

**One analysis mode per run.** The mode flags are mutually
exclusive; the one given last on the command line is the one
performed. Each writes `bondAnalysis.<suffix>`, where the suffix
names the mode.

**Defaults.** Every option's default value is read from
`bond_analysisrc.py` in the directory named by `$IMAGO_RC`.

## Usage

**Commandline:**

```
bond_analysis.py [-o OUTFILE] [-i INFILE] [-s] [-dist LIMITDIST]
    [-bs | -ba | -bl | -boo [-Ylm_l L] | -dx | -st
     | -co | -qn [-b BRIDGES] [-c IONS] | -vtk
     | -rn MINRINGLEN MAXRINGLEN]
    [-3dfull  [-a ASCALE -b BSCALE -c CSCALE -r BONDRADIUSSCALE]]
    [-3dparts [-a ASCALE -b BSCALE -c CSCALE -f FSCALE
               -r BONDRADIUSSCALE -vol X Y Z]]
    [-bf BONDINGFACTOR] [-xyz | -abc | -fract]
    [-box MIN1 MAX1 MIN2 MAX2 MIN3 MAX3 ZONE]
    [-hist HISTDELTA] [--help]
```

**General options:**

| Option   | Defaults        | Description |
|:--------:|:---------------:|:------------|
|-o        |bondAnalysis     | Base name of the output file; the mode's suffix is appended. |
|-i        |imago.skl        | Input structure file. Other formats such as structure.dat and POSCAR are planned, but at present only the skeleton works. |
|-s        |Off              | Use Imago bond order results to support the analysis. The file expected is `S1__S1.dat`, made with `makeBOND.py -gc` and a group control file containing only `SYSTEM_NUM`, often with `-minbo 0.0`. |
|-dist     |4.0              | Cutoff distance in Angstroms beyond which atom pairs are not considered. Make sure the largest expected bond length does not exceed it. |
|-bf       |1.1              | Bonding factor: atoms are bonded if their distance is at most (sum of covalent radii) times this. |
|-xyz, -abc, -fract |abc fractional | Coordinate system for printed positions: Cartesian, direct a,b,c, or fractional. Also sets how `-box` values are interpreted. |
|-box      |Whole cell       | Use only atoms inside (ZONE 1) or outside (ZONE 2) the box bounded by MIN and MAX on each axis for the statistics. `max` may be given for any MAX to mean the full extent. |
|-hist     |Off              | Also write a histogram with the given bin spacing, to `<outfile>.hist`. |
|--help    |N/A              | Print usage information and exit. |

## Analysis Modes

| Mode      | Suffix | Output |
|:---------:|:------:|:-------|
|-bs        |.bs     | A ball-and-stick file for third-party drawing software: the lattice in a,b,c and x,y,z form, the atom count, then one line per atom with element, species, coordinates, and the numbers of the atoms it bonds to. |
|-vtk       |.vtk    | The same model in legacy VTK format, readable by Paraview. |
|-ba        |.ba     | Bond angles. One group per atom: the number of angles, then each angle as two triplets (element and species of the three atoms with the vertex in the middle, then their atom numbers) and the angle in degrees. |
|-bl        |.bl     | Bond lengths. One group per atom: the number of bonds and the common atom's element and species, then each bonded atom's element, species, bond length and number. |
|-boo       |.boo    | Bond orientational order number of each atom, from the Y_lm spherical harmonic of order L given by `-Ylm_l` (0, 1 or 6; default 6). |
|-dx        |.bx     | An openDX model with bonds from covalent radii and atom colours and sizes from the element data. |
|-st        |.st     | Statistics per element, species and atom: average bond length and its standard deviation, average bond angle and its standard deviation. |
|-co        |.co     | The coordination of each atom and a summary per element. |
|-qn        |.qn     | The Q^n distribution: for each network-forming cation (`-c`, default `si`), how many of its bridging anions (`-b`, default `o`) connect to another cation. Give several elements as a quoted, space-separated list. |
|-rn        |.rn     | The ring distribution for ring sizes from MINRINGLEN to MAXRINGLEN. |
|-3dfull    |.scad   | One OpenSCAD model of the whole system, atoms and bonds unioned, for a 3D printer. |
|-3dparts   |.scad files | A set of OpenSCAD part files (`atomTopSheet<n>`, `atomBotSheet<n>`, `bondSheet<n>`) that are printed separately and assembled into a large model. |

**3D printing sub-options**, for `-3dfull` and `-3dparts`:

| Sub-option | Default | Description |
|:----------:|:-------:|:------------|
|-a          |10.0     | Scale of the atom radius relative to the covalent radius. |
|-b          |1.0      | Scale of the bond stick diameter, by default 10 percent of the covalent radius of the largest atom. |
|-c          |12.0     | Overall scale of the lattice. |
|-r          |0.1      | Ratio between bond rod radius and covalent radius. |
|-f          |0.9      | `-3dparts` only: font size factor for the part labels, relative to the bond cylinder radius. |
|-vol        |150 150 150 | `-3dparts` only: the printable volume in mm. |

Bad choices for these values produce bad-looking models.

## Port Notes

This script is the Python port of the Perl `bondAnalysis`. Three
Perl bugs were fixed in the port: the `-Ylm_l` value was read off by
one, only three of the four `-3dfull` sub-options could be given in
one invocation, and the documented `-f` sub-option of `-3dparts` was
not implemented. The `-frag` mode of the Perl script was disabled
there and is not carried over.
