---
layout: default
title: makeSGDB.py
script: makeSGDB.py
description: 'Builds the space group database consumed by applySpaceGroup. Reads a plain-text table of all 230 space groups in x,y,z operation notation, decomposes every operation into a rotation matrix and fractional shift, expands non-primitive lattices, writes one file per space-group setting and creates the numeric symbolic links by which skeletons name a setting.'
parent: Script Catalog
---

# makeSGDB

Creates the database of space groups and their symmetry operations
that Imago uses when it expands a skeleton's asymmetric unit into a
full cell. The database is consumed by the `applySpaceGroup` Fortran
program, and the `space` line of every `imago.skl` names one of its
entries (for example `227_b`).

The script reads a specially formatted plain-text data file listing
every crystallographic space group with its symmetry operations:
rotations, reflections, inversions and fractional translations. Each
space group may have one or more sub-groups, which are its alternate
settings and origin choices. For every sub-group the script:

1. Parses the operations written in compact `x,y,z` form. Here `x`,
   `y` and `z` refer to the possibly non-orthogonal `a`, `b`, `c`
   lattice axes, the convention used by *A Hypertext Book of
   Crystallographic Space Group Diagrams and Tables* from Birkbeck
   College, University of London.
2. Decomposes each operation into a 3x3 rotation or reflection
   matrix and a three-component fractional shift vector.
3. Handles non-primitive lattices by applying shifted repetitions of
   the base operations (body-centred, face-centred and so on).
4. Writes the processed operations to one file per sub-group, named
   after the space group symbol.
5. Creates symbolic links so that each setting can be referenced by
   its numeric ID with a letter suffix: `227_a`, `227_b`, and so on.

The files and links are written in the **current working
directory**, so run the script inside the `spaceDB` directory you
intend to build. The distributed `spaceDB.tgz` that
`unpackImagoDB.py` unpacks already contains the result, and its
`remake` script re-runs this conversion, so you need this script
only when the source table is changed.

## Usage

**Commandline:**

```
makeSGDB.py -i SPACEGROUPFILE [--help]
```

**Option Table:**

| Option         | Outputs                          | Defaults          | Description |
|:--------------:|:--------------------------------:|:-----------------:|:------------|
|-i / --input    |One file per setting, plus numeric links | None; required | The data file of space group definitions to convert. |
|--help          |N/A                               |N/A                | Print usage information and exit. |

## Naming Convention

The end of the script carries a table of the uniform ASCII notation
used for point groups, space groups and crystal systems, since
Hermann-Mauguin symbols contain overbars and subscripts that plain
text cannot carry. In that notation a `~` after an axis marks an
inversion axis, the neutral screw axes (2₁, 4₂, 6₃) use `'`, the
right-handed screw axes use `]` (or `}` for 6₂), the left-handed ones
use `[` (or `{` for 6₄), and the point-group column uses Schoenflies
symbols. For example, space group 227 (diamond) is `oh..7  Fd3m`.
Read that table when you need to match a symbol in the source file
to an International Tables entry.
