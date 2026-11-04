---
layout: default
title: condense.py
script: condense.py
description: 'Creates every input file for a LAMMPS condensation of precursor molecules with bond/react reaction chemistry. Reads a condense.in keyword file naming the molecules, cell, compression stages and reactions, gathers reaction templates from the precursor database, packs the molecules with packmol, and writes the LAMMPS data, input and SLURM files.'
parent: Script Catalog
---

# condense

Creates all the input files needed to run a LAMMPS **condensation**:
a molecular dynamics run in which a box of precursor molecules is
compressed while `fix bond/react` lets them react and link into an
amorphous network. This is how Imago's amorphous models of
precursor-derived ceramics (a-SiCN:H, a-BC:H and the like) are
built. The script is driven by one keyword file, `condense.in`, and
is itself driven by [deadmd.py](deadmd.html) when the condensation
parameters are being searched.

**What a run does:**

1. **Prepare molecules and reactions.** Reads the input file and, if
   the precursor database does not already hold them, runs
   [make_reactions.py](make_reactions.html) for each reacting pair of
   molecules to produce the pre- and post-reaction templates that
   drive the bonding.
2. **Run packmol.** Grabs the skeleton file of each molecule,
   converts it to PDB if needed, and packs the requested number of
   each into the cell with the external `packmol` program.
3. **Prepare the LAMMPS input.** Converts the packed PDB back to a
   skeleton and then to a LAMMPS data file, gathers the reaction
   templates and reconciles every bond and angle type across all
   files, and writes the LAMMPS input file and a SLURM submission
   file.

**Defaults.** The input file name, collision mode, cell size, speed
and template directory defaults are read from `condenserc.py` in the
directory named by `$IMAGO_RC`.

## Usage

**Commandline:**

```
condense.py [-i INPUTFILE] [-f] [-s SEED] [--help]
```

**Option Table:**

| Option          | Defaults                     | Description |
|:---------------:|:----------------------------:|:------------|
|-i / --input     |condense.in                   | Name of the condense input file. |
|-f / --force     |Off                           | Force-collision mode: arrange the molecules on a collision course near the cell centre instead of packing them. Only works for three or fewer molecules. |
|-s / --seed      |Drawn at random               | Random seed, a positive integer, for every random draw: the packmol seed, the bond/react seeds and the molecule velocities. The same seed and input give the same LAMMPS input. Overrides a `seed` line in the input file. The seed used is recorded in the file `seeds` either way. |
|--help           |N/A                           | Print usage information and exit. |

## The Input File

The file is a sequence of keyword sections. A keyword may be
followed by a value giving the number of lines in its section.

### composition

```
composition 3
b10c2h12_1 Family2 30
c2h6_1     Family1 3
ch4_1      Family1 8
```

The value is the number of different molecule types. One line
follows per type: **molecule name**, **family name**, **count**. The
lines **must be alphabetized** by molecule name.

A molecule name has two parts joined by `_`: the chemical formula and
an integer distinguishing a particular atom-type assignment or
conformation of that formula. The precursor database must hold a
directory of that name containing at least `<name>.skl`. For
`b10c2h12`, say, some boron atoms have carbon neighbours and some do
not, and the two carbons can sit para, meta or ortho; each variant
gets its own integer.

All molecules with the same **family name** share element-species
types: a `c1` from ch4 and a `c1` from c2h6 in the same family are
one LAMMPS type, while in different families they are distinct. To
keep a molecule out of any family, give it its own name as family.

### cell_size, max_speed, seed

`cell_size` is the cell edge in Angstroms, `max_speed` the maximum
initial speed of each molecule, and `seed` an optional positive
integer that seeds every random draw (the `-s` option overrides it).

### stage

```
stage 0.33 500 nvt 800 800 500 100000
stage 1.0  500 nvt 800 300 500 50000
```

Seven values: the **squish factor** (final fractional cell size
after compression; 0 means no deformation), the **squish step size**
(apply the deformation every this many timesteps), the **ensemble**
(`nve`, Newtonian, or `nvt`, Nose-Hoover thermostat), the **start
and end temperatures** in K, the **thermostat damping** in fs (the
Tdamp of `fix nvt`, ignored for nve), and the **run steps**. The
keyword may be repeated for successive stages.

### target_density (optional)

A density in g/cc. Every compressing stage then stops squishing as
soon as the live density reaches the target, checked every 10 steps,
and runs the rest of its steps at the fixed box. The squish factor
becomes the most the box may shrink in that stage. The box changes
only every squish step, so a smaller step gives a smaller overshoot.

### reactions

```
reactions 4
c6h19nsi2_1 si-1 c6h19nsi2_1 c-1 0.99 inter
c6h19nsi2_1 c-1  c6h19nsi2_1 c-1 0.01 inter
c6h19nsi2_1 si-1 c6h19nsi2_1 n-1 0.99 inter
c6h19nsi2_1 si-1 c6h19nsi2_1 h-1 0.0005 5.0 intra rmin 2.0
```

The value is the number of reaction kinds, one line each:

```
<MOL1> <SITE1> <MOL2> <SITE2> <PROB> [RMAX] inter|intra [rmin R] [max_rxn N]
```

The molecule names must match the composition section and the sites
(`element-type`) must exist in those molecules. PROB is the
probability, 0 to 1, that the reaction fires once the fragments are
close enough. The optional RMAX is the largest initiator-atom
distance in Angstrom at which `fix bond/react` may fire (default 5.0).
`inter` means the bond forms between two molecules, `intra` within
one. `rmin R` sets the smallest initiator distance (default 0.0).
`max_rxn N` caps the reaction at N events; it is a cap, never a
quota.

The working expectation is that every reaction discharges an H atom
from each molecule before the underlying atoms bond. The templates
come from the precursor database, from the directory
`<MOL1>__<MOL2>/rxnTemplates/`. Every angle implied by a post
template's bonds but missing from it is added automatically, and any
new bond or angle type is added to the data file.

### bond_override (optional)

```
bond_override 2
B B 1.7800 uff
B C 1.7100 uff
```

One line per overridden element pair: the two elements, the
equilibrium bond length R0 in Angstrom replacing the UFF value, and
the force constant K as either a number used verbatim (kcal/mol/A^2,
with E = K(r - r0)^2, and *not* scaled by `bond_parameter_scale`),
`uff` to rescale the UFF constant to the new length (K scales as
r0^-3, so a longer bond becomes softer), or `keep` to leave K at the
UFF value for the original length. The keyword is purely additive:
an input without it produces the same output as before it existed.

This exists because the UFF r0 is a sum of element natural radii
with no input from the molecule being modelled, and can be far from
the real bond length of a specific compound. Measured examples:
silica Si-O 1.716 against 1.608 in the glass; o-carborane B-B 1.676
against about 1.78 experimental.

## Output

A `packmol/` directory with the packmol input and packed PDB, and a
`lammps/` directory with `lammps.dat`, `lammps.in`, the reaction
template and map files, and a SLURM submission script. The generated
input uses `special_bonds lj/coul 0 0 0.5`, a 0.1 fs timestep, logs
every reaction's event count every 500 steps with the final counts
appended to `rxn_counts_evolve`, and writes `final.data` after the
final minimisation.
