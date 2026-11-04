---
layout: default
title: mod_struct.py
script: mod_struct.py
description: 'Modifies an imago.skl structure by applying an ordered sequence of operations given on the command line: supercell replication with optional mirroring, translation, rotation about an axis or plane normal, reorientation of a Miller direction to z for surfaces, vacuum insertion, sphere and block cutting, orthorhombic conversion and random perturbation.'
parent: Script Catalog
---

# mod_struct

Reads an Imago skeleton file (`imago.skl` by default) and applies a
user-specified sequence of modifications to the structure, writing
the result to a new skeleton. The available operations are supercell
generation, translation, rotation, surface preparation, filtering,
vacuum insertion, sphere and block cutting, orthorhombic conversion
and random perturbation.

**Operations are modular and sequential.** Each of the operation
flags (`-sc`, `-trans`, `-rotA`, `-rotP`, `-rothkl`, `-filter`,
`-addvac`, `-cutsphere`, `-cutblock`, `-ortho`, `-perturb`) can be
repeated any number of times, and they are applied **in the order
given on the command line**. So a surface slab is typically built as
"supercell, then translate to choose the surface plane, then add
vacuum", and the command line reads in that order.

**Defaults.** The input and output file names and the default
axis-assignment order are read from the resource control file
`mod_structrc.py` in the directory named by `$IMAGO_RC`.

## Usage

**Commandline:**

```
mod_struct.py [-i INFILE] [-o OUTFILE] [-t TITLE]
              [-abcxyz ABC1 XYZ1 ABC2 XYZ2]
              [-sc SC1 SC2 SC3 MIRROR1 MIRROR2 MIRROR3]
              [-trans TX TY TZ]
              [-rotA RX RY RZ ANGLE [-orig OX OY OZ]]
              [-rotP P1X P1Y P1Z P2X P2Y P2Z P3X P3Y P3Z
                     [-angle ANGLE] [-orig OX OY OZ]]
              [-rothkl H K L]
              [-filter MINDIST]
              [-addvac AXIS AMOUNT]
              [-cutsphere [-atom ATOM | -atxyz X Y Z | -atabc A B C]
                          [-radius R] [-zone in|out]]
              [-cutblock [-zone in|out]
                         <-xyz FROMX TOX FROMY TOY FROMZ TOZ |
                          -abc FROMA TOA FROMB TOB FROMC TOC>]
              [-ortho]
              [-perturb MAGNITUDE]
              [--help]
```

**Examples:**

```
# 2x2x2 supercell, no mirroring
mod_struct.py -sc 2 2 2 0 0 0 -o imago.skl.222

# Slab: double along c, shift so the wanted plane is at the cell
#   boundary, then open a 10 A vacuum gap along c
mod_struct.py -sc 1 1 2 0 0 0 -trans 0 0 1.5 -addvac c 10

# Cluster: keep only atoms within 6 A of atom 17
mod_struct.py -cutsphere -atom 17 -radius 6 -zone out

# Shake every atom by up to 0.05 A for a disorder test
mod_struct.py -perturb 0.05
```

**Variable Table:**

| Variable                | Accepted Values           | Description |
|:-----------------------:|:-------------------------:|:------------|
|INFILE, OUTFILE          |String; file path          | Input and output skeleton files. |
|TITLE                    |String; quote if it has spaces | Title line of the output skeleton. |
|ABC1, ABC2               |1, 2 or 3                  | Lattice vector indices (a=1, b=2, c=3). |
|XYZ1, XYZ2               |1, 2 or 3                  | Cartesian axis indices (x=1, y=2, z=3). |
|SC1, SC2, SC3            |integer                    | Replication counts along a, b, c. |
|MIRROR1, MIRROR2, MIRROR3|0 or 1                     | Whether every other replica along that axis is mirrored. |
|TX, TY, TZ               |float; Angstroms           | Translation vector in direct Cartesian coordinates. |
|RX, RY, RZ               |float                      | Rotation axis direction in Cartesian coordinates. |
|ANGLE                    |float; degrees             | Rotation angle. |
|OX, OY, OZ               |float; Angstroms           | Origin through which the rotation axis passes. |
|P1, P2, P3               |three floats each          | Three points defining a plane. |
|H, K, L                  |float                      | Miller-index direction, in units of the lattice-vector magnitudes. |
|MINDIST                  |float; Angstroms           | Minimum allowed interatomic distance. |
|AXIS                     |a, b, or c                 | Lattice axis to lengthen. |
|AMOUNT                   |float; Angstroms           | Vacuum thickness to insert. |
|ATOM                     |integer; atom number       | Atom at the centre of a sphere. |
|X Y Z / A B C            |floats                     | Sphere centre in Cartesian or fractional coordinates. |
|R                        |float; Angstroms           | Sphere radius. |
|FROM.., TO..             |floats                     | Block bounds in Cartesian (`-xyz`) or fractional (`-abc`) coordinates. |
|MAGNITUDE                |float; Angstroms           | Maximum random displacement. |

**Option Table:**

| Option        | Repeatable | Defaults            | Description |
|:-------------:|:----------:|:-------------------:|:------------|
|-i             |No          |imago.skl            | Input skeleton file. |
|-o             |No          |imago.skl.new        | Output skeleton file. |
|-t             |No          |hostname and working directory | Title for the output file. |
|-abcxyz        |No; once only |1 1 2 2            | Orientation of the lattice vectors with respect to the Cartesian axes. See below. |
|-sc            |Yes         |No supercell         | Replicate the cell, optionally mirroring alternate replicas. Overrides any supercell line in the input file. |
|-trans         |Yes         |None                 | Translate every atom. `0 0 0` requests centring the model in the cell. |
|-rotA          |Yes         |Origin 0 0 0         | Rotate every atom about an axis. |
|-rotP          |Yes         |Origin 0 0 0; angle from the points | Rotate every atom about the normal of a plane. |
|-rothkl        |Yes         |None                 | Rotate a Miller direction to z and enlarge the cell for a surface. **Presently only works for simple cases.** |
|-filter        |Yes         |None                 | Remove atoms closer than MINDIST to another atom. **Does not really work yet.** |
|-addvac        |Yes         |None                 | Lengthen one lattice axis to create a vacuum region. |
|-cutsphere     |Yes         |Radius 4 A; zone in; centre by sub-option | Remove atoms inside or outside a sphere. |
|-cutblock      |Yes         |Zone in              | Remove atoms inside or outside a rectangular block. |
|-ortho         |Yes         |None                 | Convert the lattice to orthorhombic. **Hexagonal cells only, after suitable supercell doubling.** |
|-perturb       |Yes         |None                 | Random shift of every atom. |
|--help         |N/A         |N/A                  | Print usage information and exit. |

## Option Details

### Axis Orientation

> -abcxyz ABC1 XYZ1 ABC2 XYZ2

Defines how the a, b, c lattice vectors sit in Cartesian space. ABC1
is aligned with XYZ1, and ABC2 lies in the XYZ1-XYZ2 plane, so that
the third lattice vector can point in an arbitrary direction. The
default `1 1 2 2` aligns a with x and puts b in the xy plane. This
option can be applied only once and governs every Cartesian-space
operation that follows (translation, rotation, cutting by `-xyz`).

### Supercell

> -sc SC1 SC2 SC3 MIRROR1 MIRROR2 MIRROR3

Replicates the cell SC1, SC2 and SC3 times along a, b and c. The
three MIRROR values (0 or 1) control whether every other replica
along that axis is a mirror image, which is useful for building
twin boundaries or for making a cell orthorhombic. This overrides any
`supercell` line already in the input skeleton.

### Translation

> -trans TX TY TZ

Translates all atomic positions by (TX, TY, TZ) in direct Cartesian
coordinates, wrapping under periodic boundary conditions. The special
value `0 0 0` instead centres the model within the cell, which is
the usual want for a molecule in a box.

### Rotation

> -rotA RX RY RZ ANGLE [-orig OX OY OZ]

Rotates all atoms by ANGLE degrees about the axis (RX, RY, RZ), a
direction in Cartesian coordinates. The axis passes through the
origin unless `-orig` moves it. Periodic boundary conditions are
retained: atoms rotated out of the cell are wrapped back in as if
they came from a neighbouring replica.

> -rotP P1X P1Y P1Z P2X P2Y P2Z P3X P3Y P3Z [-angle ANGLE] [-orig OX OY OZ]

Rotates all atoms about the normal of the plane through the three
points P1, P2, P3. The angle is either given with `-angle` or, when
omitted, computed as the angle between the second and third points.
Periodic boundary conditions are retained as for `-rotA`.

### Surface Preparation

> -rothkl H K L

Reorients and resizes the crystal lattice so that the Miller
direction (H, K, L), expressed in units of the lattice-vector
magnitudes, coincides with the Cartesian z axis. The simulation cell
is then enlarged to retain periodicity. Follow it with `-addvac c`
to open the surface. **Presently this only works for simple cases**;
check the output cell carefully.

### Vacuum

> -addvac AXIS AMOUNT

Increases the cell length along AXIS (a, b or c) by AMOUNT Angstroms
to create a vacuum region. The atoms keep their direct-space
coordinates, so the vacuum opens at the far end of the axis. A
translation is typically performed beforehand to choose which plane
becomes the surface.

### Cutting

> -cutsphere [-atom ATOM | -atxyz X Y Z | -atabc A B C] [-radius R] [-zone in|out]

Removes atoms inside (`-zone in`, the default) or outside
(`-zone out`) a sphere of radius R (default 4 A) centred on an atom
number, a Cartesian point, or a fractional point. `-zone out` with a
large radius carves a cluster out of a bulk model.

> -cutblock [-zone in|out] <-xyz FROMX TOX FROMY TOY FROMZ TOZ | -abc FROMA TOA FROMB TOB FROMC TOC>

Removes atoms inside (default) or outside a rectangular block given
as coordinate ranges in either Cartesian (`-xyz`) or fractional
(`-abc`) space.

### Orthorhombic Conversion

> -ortho

Converts the lattice to orthorhombic type. This should only be used
on cells that have been carefully analysed. At present only hexagonal
cells can be converted, and only after they have been doubled with
an appropriate `-sc` so that an orthorhombic cell is commensurate.

### Perturbation

> -perturb MAGNITUDE

Applies a random shift of up to MAGNITUDE Angstroms to every atom.
Useful for breaking symmetry before a relaxation or for testing the
sensitivity of a result to small disorder.
