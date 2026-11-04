---
layout: default
title: make_reactions.py
script: make_reactions.py
description: 'Creates the pre-reaction and post-reaction molecule template files and the atom map that LAMMPS fix bond/react needs for one pair of precursor molecules, by identifying leaving (surface) atoms on each, orienting the molecules for bonding, merging them, and pruning the result to the atoms within a chain length of the new bond.'
parent: Script Catalog
---

# make_reactions

Creates a set of pre- and post-reaction molecule files, along with a
mapping file, for a `fix bond/react` based LAMMPS calculation. These
are the **reaction templates** that tell LAMMPS what a reacting pair
of fragments looks like before and after a bond forms. The script is
run by [condense.py](condense.html) for any reacting pair the
precursor database does not already hold; you run it by hand to add
a reaction to the database or to tune one.

## Algorithm

The script takes two directories, each identifying a precursor
molecule, and creates a new directory named from their combination
(`DIR1__DIR2`). Inside it the skeleton files are copied from the
original directories and processed.

**Surface atoms and the delete list.** Each skeleton is parsed to
find every **surface atom**, called an S atom (not to be confused
with sulfur), whose element is set by `-s`. For each S atom the
script identifies the set of attached atoms to be deleted when the
molecules bond, the **delete list**. By default only the S atoms
themselves are deleted. The atoms bonded to the S atoms are the
**trigger atoms**, the points at which the new bond forms. The
physical picture is that the S atoms break off to form S2 and the
bare trigger atoms bond; in the simulation the S2 is never created,
the two S atoms are simply deleted. With `-d1` or `-d2` a deeper
tree of atoms is removed instead, and the binding point moves inward.

**Orientation and merging.** Once all atoms are catalogued, the
second molecule is systematically translated and oriented so that
each of its S atoms sits in a proper position for reacting with each
S atom of the first molecule. For each such pairing the two skeletons
are merged into one file.

**Templates.** From each merged file two LAMMPS molecule files are
created with the same atoms and coordinates but different bonding
topology. The **pre-reaction** template holds the relevant segments
of both molecules, where relevant is governed by the **chain length**,
4 by default: the atom bonded to the S atom, its nearest neighbours,
and its next-nearest neighbours. The **post-reaction** template holds
the same atoms with the new bond and the deletions applied. Four is
the default because the angles at the S-bonded atom are changed, and
an angle can reach two bonds deeper into the molecule, so for a
chain S-B-C-X the X atom must be present for a possible B-C-X angle.

**Defaults.** Every option's default value is read from
`make_reactionsrc.py` in the directory named by `$IMAGO_RC`.

## Usage

**Commandline:**

```
make_reactions.py -m1 DIR1 -m2 DIR2 [-c1 N] [-c2 N]
                  [-d1 DEPTH TARGET] [-d2 DEPTH TARGET]
                  [-s ELEMENT DIST DEL_BOTH] [-uptypes STEP]
                  [-consttypes] [--help]
```

**Variable Table:**

| Variable   | Accepted Values          | Description |
|:----------:|:------------------------:|:------------|
|DIR1, DIR2  |String; directory path    | Directories containing the two molecules' skeleton files. |
|N           |integer                   | Chain length: bond-hops from the S atom to include in the template. |
|DEPTH       |integer                   | Depth of the delete tree. |
|TARGET      |element_species, e.g. `si_1` | The atom that becomes the binding point when the tree is deleted. |
|ELEMENT     |Element symbol            | The surface-atom element. |
|DIST        |float; Angstroms          | Separation between the two S atoms in the templates. |
|DEL_BOTH    |1 or 0                    | Delete both S atoms, or only the one from molecule 2. |
|STEP        |integer                   | Increment in type number for bonded atoms. |

**Option Table:**

| Option        | Defaults                 | Description |
|:-------------:|:------------------------:|:------------|
|-m1, -m2       |None; required            | Directories of the first and second molecule. |
|-c1, -c2       |4                         | Chain length for molecule 1 and 2. For CH4 with hydrogen as S, a chain length of 3 includes the S atom, the central C, and the other three H. |
|-d1, -d2       |Depth 0; only the S atom deleted | Delete tree for molecule 1 and 2. See below. |
|-s             |h 0.5 1                   | Surface-atom element, S-S distance, and whether to delete both S atoms. |
|-uptypes       |0                         | Atoms that bind, and their neighbours, are assigned a type number higher by STEP than the highest existing type for that element, so bonded atoms become distinguishable. 0 keeps types constant. |
|-consttypes    |Off                       | Keep the type numbers of each molecule unchanged when merging, so a type-1 carbon from each molecule stays type 1 in the merged system. |
|--help         |N/A                       | Print usage information and exit. |

## Option Details

### Delete Trees

> -d1 DEPTH TARGET, -d2 DEPTH TARGET

Consider hexamethyldisilazane, a chain of CH3 groups around
Si-N-Si. With `-d1 1 si_1` the far-left CH3 group is added to the
delete tree and the left Si becomes the binding point. With
`-d1 2 n_1` every atom to the left of the N is deleted and the N
becomes the binding point. The target names the element and species
of the atom that is left bare.

## Output

The directory `DIR1__DIR2/` with working subdirectories `centered/`,
`merged/` and `pruned/`, and `rxnTemplates/` holding, for each pair
of sites, `preRxn.<SITE1>_<SITE2>.data`,
`postRxn.<SITE1>_<SITE2>.data` and `<SITE1>_<SITE2>.map`. Those three
are what condense.py copies into a run and names in `lammps.in`.
