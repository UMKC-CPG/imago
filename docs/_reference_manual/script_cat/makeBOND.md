---
layout: default
title: makeBOND.py
script: makeBOND.py
description: 'Parses the raw bond order and Q* (effective charge) output of an Imago bond calculation and writes files for visualisation: scatter-plot data grouped by element and species, openDX or VTK structural models with colour-coded bonds and atoms, a bond profile along a lattice axis, or a 3D openDX mesh of bond order.'
parent: Script Catalog
---

# makeBOND

Parses the bond order and Q\* output file of an Imago bond
calculation and creates a number of other files that can be plotted
by openDX, Paraview, VisIt, Origin, or similar tools. Bond order
measures the strength of each atom pair's bond and Q\* is the
effective charge on each atom; together they are the main chemical
bonding picture Imago offers, and this script is how you get them out
of the raw file and into a figure.

By default the bond order and Q\* data are read from `gs_bond-mb.raw`
and the atomic positions from `structure.dat`. The position file must
be in the format of the current Imago output, which has the positions
sorted by type.

**Defaults.** Every option's default value is read from the resource
control file `makeBONDrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
makeBOND.py [-data DATA_FILE] [-pos POS_FILE]
            [-bc BOND_CONTROL_FILE] [-gc GROUP_CONTROL_FILE]
            [-maxbo MAX_BO] [-minbo MIN_BO]
            [-maxbl MAX_BL] [-minbl MIN_BL]
            [-scatter]
            [-dipole]
            [-profile PROFILE_AXIS]
            [-model [-atom | -bond] [-rf RADIUS_FACTOR]
                    [-grey] [-bo3c]]
            [-vtk]
            [-mesh [-comp_atom COMP_ATOM] [-steps A B C]
                   [-fwhm FWHM] [-limit LIMIT_DIST]]
            [--help]
```

Only one output mode is active at a time. If none is requested the
default is `-scatter`.

**Variable Table:**

| Variable             | Accepted Values       | Description |
|:--------------------:|:---------------------:|:------------|
|DATA_FILE             |String; file path      | Raw bond order and Q\* file from the bond calculation. |
|POS_FILE              |String; file path      | Imago structure file with lattice vectors and atom positions. |
|BOND_CONTROL_FILE     |String; file path      | File of per-element-pair bond order and bond length limits. See below. |
|GROUP_CONTROL_FILE    |String; file path      | File naming the filters used to tag and group atoms and bonds. See below. |
|MAX_BO, MIN_BO        |float                  | Bond order limits applied to every bond regardless of element. |
|MAX_BL, MIN_BL        |float                  | Bond length limits applied to every bond regardless of element. |
|PROFILE_AXIS          |a, b, or c             | Lattice axis along which to compute the bond profile. |
|RADIUS_FACTOR         |float                  | Multiplier on the covalent radius used to draw atom spheres. |
|COMP_ATOM             |integer; atom index    | Atom whose properties are subtracted from every mesh point to form a difference mesh. 0 means no subtraction. |
|A, B, C               |integer                | Number of mesh points along each lattice direction. |
|FWHM                  |float; Angstroms       | Full width at half maximum of the Gaussian used to weight nearby bond values onto a mesh point. |
|LIMIT_DIST            |float; Angstroms       | Maximum separation between a bond point and a mesh point for the bond to contribute. |

**Option Table:**

| Option      | Outputs                                    | Defaults          | Description |
|:-----------:|:------------------------------------------:|:-----------------:|:------------|
|-data        |N/A                                         |gs_bond-mb.raw     | Raw bond order data file. |
|-pos         |N/A                                         |structure.dat      | Structure file with lattice and positions. |
|-bc          |N/A                                         |None               | Bond control file for element-pair bond order and length filtering. Applied in addition to the command-line limits. |
|-gc          |N/A                                         |ELEMENT_NAME and SPECIES_ID | Group control file defining the output grouping tags. |
|-maxbo       |N/A                                         |1000.0             | Exclude bonds with bond order above this. |
|-minbo       |N/A                                         |-1000.0            | Exclude bonds with bond order below this. |
|-maxbl       |N/A                                         |1000.0             | Exclude bonds longer than this. |
|-minbl       |N/A                                         |0.0                | Exclude bonds shorter than this. |
|-scatter     |`ATOM_<tag>.dat` and `<bond_tag>.dat`       |On if no other mode| Scatter-plot data files grouped by atom tag and bond tag. |
|-dipole      |N/A                                         |Off                | Compute and print dipole information. Accepted for compatibility; at present it has no effect. |
|-profile     |bondProfile.plot                            |Off                | Bond profile along the given lattice axis. |
|-model       |`atomPos_<tag>.dx`, `<tag>.dx`, lattice.dx  |Off                | openDX structural model files with colour-coded atoms and bonds. |
|-atom        |N/A                                         |On                 | Make the openDX model, or the mesh, atom-based. |
|-bond        |N/A                                         |Off                | Make the openDX model bond-based. |
|-rf          |N/A                                         |1.0                | Multiplicative factor for atom sphere radii in the model. |
|-grey        |N/A                                         |Off                | Draw atom spheres in grey-scale instead of colour. |
|-bo3c        |N/A                                         |Off                | Show three-centre bonds in the model where they occur. |
|-vtk         |`<tag>.vtk`                                 |Off                | Structural model in VTK legacy format. |
|-mesh        |`<tag>.dx`, lattice.dx                      |Off                | 3D openDX mesh of bond order values. |
|-comp_atom   |N/A                                         |0                  | Atom index for a difference mesh. |
|-steps       |N/A                                         |10 10 10           | Mesh points along a, b, c. |
|-fwhm        |N/A                                         |2.0                | Gaussian FWHM for mesh weighting. |
|-limit       |N/A                                         |4.0                | Maximum bond-to-mesh-point distance in Angstroms. |
|--help       |N/A                                         |N/A                | Print usage information and exit. |

## Atoms: Element, Species and Type

Four distinctions matter for how this script groups its output.
There is the total number of atoms in the system, and those atoms are
divided three ways.

**Element** is the elemental identity. An atom is a carbon, a
silicon, a nitrogen, and so on. The set of elements is the smallest
set. It may hold only two members (Si and N for g-Si3N4). If a carbon
is substituted for a silicon in the tetrahedral site the set grows by
one (Si, N and C); a second carbon substituted into the octahedral
site leaves the set unchanged. The element set deals only with which
elements are present, not with their relations within the system.

**Species** groups atoms that share the same global symmetry or a
similar local configuration before any modification of the system.
Crystalline g-Si3N4 has two species of Si, one tetrahedral and one
octahedral, and still one species of N. For complex systems with
little or no symmetry, species are defined by local configurational
similarity: atoms with similar neighbours at similar distances are
treated as the same species.

**Type** groups atoms that the calculation treats as equivalent, i.e.
that share one set of potential coefficients. Often the types are
exactly the species. Sometimes atoms of one species are deliberately
given different types so that they have different potentials, as in
XANES calculations, or when all atoms in a region are assigned
distinct types because they cannot be grouped any other way.

The raw bond order file records all three for every atom, and the
group control file picks which of them build the tags that split the
output into files.

## Option Details

### Bond Control File

> -bc

Filters bonds by bond length and bond order for each pair of
elements. The file is a series of lines, each with five
space-separated columns:

```
BO  si  n   0.1  1.5
BL  si  n   1.5  2.0
```

The first column is a tag, `BO` or `BL`, saying whether the line
restricts bond order or bond length. The second and third are the two
element names. The fourth and fifth are the minimum and maximum of
the allowed range for that element pair. Give as many lines as you
wish; they apply *in addition to* the command-line limits.

### Group Control File

> -gc

Controls how the output is grouped, which often corresponds to the
definitions in a HIN file. The file is a series of space-separated
filter names available in the raw data file. Each name contributes to
the tag built for every atom and bond. For example, with the filter
names `ELEMENT_NAME SPECIES_ID` each atom is tagged with its element
name and species number, and those tags are what separate the atoms
into files. The default grouping is exactly that pair.

### Global Limits

> -maxbo, -minbo, -maxbl, -minbl

Each takes one value and excludes from every printout any bond whose
bond order or bond length is beyond it. These are often used in place
of a bond control file, but they apply to all atoms regardless of
element, so they are only useful when you want one condition for the
whole system.

### Scatter Plot

> -scatter

The default output. Creates multiple files that an associated Origin
script, or any plotting tool, can turn into a scatter plot. Each
`ATOM_<tag>.dat` holds the per-atom Q\*, bond order and bond length
statistics and the a, b, c position of every atom carrying that tag.
Each `<bond_tag>.dat` holds the bond order, bond length, atom indices
and a, b, c midpoint of every bond carrying that tag. The tags come
from the group control file or, by default, from element and species.

### Bond Profile

> -profile a|b|c

Produces a simple data file, `bondProfile.plot`, that can be plotted
as a curve indicating the bond profile along one lattice axis. The
bond profile at a point along the axis is the sum of the bond orders
of all bonds crossing that point divided by the number of such bonds.
At present one file is produced for the whole system; it is not
separated by any criterion.

### openDX Model

> -model [-atom | -bond] [-rf RADIUS_FACTOR] [-grey] [-bo3c]

Produces multiple openDX files describing a model of the system that
can be drawn with colour coding to show relative bond order and
charge transfer. For each unique tag it writes `atomPos_<tag>.dx`
with atom positions, colours and radii, and `<tag>.dx` with the Q\*
data and the bond order data (positions, connections and values), plus
`lattice.dx` for the cell outline. The atoms are drawn as spheres
coloured statically by element from the database, with radius equal
to the covalent radius by default; `-rf` multiplies that radius and
`-grey` switches to grey-scale. `-bo3c` alters the structure to show
three-centre bonds where they occur.

The most important sub-options are `-atom` and `-bond`. With `-atom`
(the default) each model file contains a specific set of atoms, and a
bond is included only if at least one of its atoms is in that set.
With `-bond` each file contains a specific set of bonds, and an atom
is included only as the end of an included bond.

### VTK Model

> -vtk

The same as the openDX model, but produces a single file per tag in
the VTK legacy format for Paraview or VisIt. Colour coding has yet to
be implemented.

### openDX Mesh

> -mesh [-comp_atom COMP_ATOM] [-steps A B C] [-fwhm FWHM] [-limit LIMIT_DIST]

Produces a 3D mesh of data points for display in openDX. The value at
each mesh point is a weighted average of nearby bond data, such as
the bond order at each bond midpoint. `-steps` defines the mesh grid
along a, b and c. `-fwhm` is the full width at half maximum of the
Gaussian used to weight nearby values. `-limit` is the maximum
distance between a bond point and a mesh point for the bond to
contribute. `-comp_atom` subtracts the properties of one specific
atom from every mesh point to create a difference mesh.
