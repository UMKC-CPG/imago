---
layout: default
title: deadmd_validity.py
script: deadmd_validity.py
description: 'Pass/fail gate for a condensed LAMMPS model before deadmd.py scores it. Checks that every atom has an allowed bond count for its element, that the final atom count matches the reactions that fired, and that no bond is strained beyond a tolerance. Writes validity.txt and exits 0 for valid, 1 for invalid.'
parent: Script Catalog
---

# deadmd_validity

A pass/fail check of a condensed model before it is scored. A model
that fails is not a *worse* model, it is an *invalid* one, so
[deadmd.py](deadmd.html) gives it the crash sentinel instead of a
fitness. The principle is "gate first, score second": a structure
with an impossible bonding pattern must never be allowed to win a
generation on energy.

The script reads, from a member's `lammps/` directory, all of which
[condense.py](condense.html) wrote:

* `final.data`, the structure after the final minimisation;
* `lammps.dat`, for the starting atom count;
* `lammps.in`, to learn which map file belongs to which reaction;
* the reaction map files, which list the atoms each event deletes;
* `rxn_counts_evolve`, the number of events per reaction.

## Usage

**Commandline:**

```
deadmd_validity.py [-dir DIR] [-valence EL N [N ...]]... [-max_strain S]
```

Run it inside a member's `lammps/` directory, or point `-dir` at one.

**Option Table:**

| Option        | Defaults  | Description |
|:-------------:|:---------:|:------------|
|-dir           |.          | The lammps directory to check. |
|-valence       |See below  | Allowed bond counts for one element, e.g. `-valence C 4 6`. Repeatable, one element per use. |
|-max_strain    |0.30       | Maximum allowed fractional deviation of a bond length from its `Bond Coeffs` r0. |

## Checks

**Hard checks.** Any failure makes the model invalid.

| Check         | Rule |
|:--------------|:-----|
|coordination   | Every atom's bond count is in the allowed set for its element. Defaults: H 1, C 4, N 3, O 2, F 1, Si 4, P 3 or 4, S 2, Cl 1, Ge 4, B 3 to 6. This catches reactions that delete their leaving atoms but never write their new bond, leaving an atom short. |
|bookkeeping    | Final atoms = starting atoms minus the sum over reactions of (atoms deleted per event times events). |
|strain         | Every bond length is within plus or minus `max_strain` of its r0. |

**Reported only.** These never fail the gate: the number of connected
fragments, angles missing or extra relative to the bond topology
(linkage angles are deliberately floppy in some runs), and non-bonded
close contacts under 1.0 Angstrom.

## Output

`validity.txt` in the checked directory, holding the result of each
check, and an exit code of 0 for valid or 1 for invalid. In a deadmd
run the gate is switched on with `validity_gate yes` in `deadmd.in`,
and `valence` and `max_strain` lines there pass through to these
options.
