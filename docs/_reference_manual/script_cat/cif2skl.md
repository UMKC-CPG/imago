---
layout: default
title: cif2skl.py
script: cif2skl.py
description: 'Converts a CIF crystal structure file into an imago.skl skeleton file while preserving the space group. The skeleton stores the asymmetric unit plus a space-group token rather than a flattened P1 cell.'
parent: Script Catalog
---

# cif2skl

Converts a CIF (Crystallographic Information File) into an `imago.skl`
skeleton file while preserving the crystal's space group. CIF is the
standard interchange format for crystal structures: it is what the
Crystallography Open Database serves and what journals publish. The
skeleton is Imago's own structure format, so this script is the bridge
between the two.

**Why the space group matters.** An `imago.skl` is *not* a flat list
of every atom in the cell. It stores the **asymmetric unit** together
with a **space-group token**, and Imago regenerates the full cell from
those two things when it reads the skeleton. Imago's Brillouin-zone
integration samples the irreducible wedge *using that space group*, so
flattening a crystal to P1 would silently corrupt every
k-point-dependent result. This converter therefore recovers the CIF's
space group rather than discarding it, and it refuses to fall back to
P1 for a crystal.

**Suggested usage:** run it once on the CIF you downloaded, then copy
the resulting skeleton into its own job directory as `imago.skl`
before running `makeinput.py`. To obtain CIFs directly from the
Crystallography Open Database see
[cod_fish.py](cod_fish.html).

## Usage

**Commandline:**

```
cif2skl.py INPUT.cif [OUTPUT.skl] [-s TOKEN] [-t STR] [--help]
```

**Examples:**

```
cif2skl.py au_fcc.cif                 # writes au_fcc.skl beside it
cif2skl.py si.cif si_diamond.skl      # writes the named skeleton
cif2skl.py si.cif --space 227_b       # force a spaceDB setting
cif2skl.py si.cif -t "Si diamond"     # set the skeleton title
```

**Variable Table:**

| Variable   | Accepted Values               | Description |
|:----------:|:-----------------------------:|:------------|
|INPUT.cif   |String; file path              | Path to the input CIF file. Required. |
|OUTPUT.skl  |String; file path              | Path for the output skeleton. Optional. |
|TOKEN       |spaceDB setting token          | A space-group setting token of the form `<IT#>_<letter>`, e.g. `227_b`. The integer is the International Tables space-group number; the letter picks one of the settings stored in Imago's `spaceDB`. |
|STR         |String                         | Free text recorded in the skeleton's title block. |

**Option Table:**

| Option        | Outputs      | Defaults                                   | Description |
|:-------------:|:------------:|:------------------------------------------:|:------------|
|INPUT.cif      |N/A           |None; required                              | The CIF to convert. |
|OUTPUT.skl     |Skeleton file |Input path with its extension changed to `.skl` | Where to write the skeleton. |
|-s / --space   |N/A           |Automatic resolution                        | Force a specific `spaceDB` setting token instead of resolving it automatically. The forced setting is still verified against the CIF's expansion. |
|-t / --title   |N/A           |A note naming the source CIF                | Title recorded in the skeleton's title block. |
|--help         |N/A           |N/A                                         | Print usage information and exit. |

## How the Space Group is Resolved

The converter works in three steps.

1. **ASE parses the CIF.** This is the only step that touches the
   Atomic Simulation Environment, and it spares Imago a hand-written
   CIF parser: ASE handles CIF's flexible `loop_` blocks and its
   symmetry expansion. From it the script reads the *authored
   asymmetric unit* (the raw `_atom_site_*` list), the cell, the
   International Tables space-group number and origin setting, and
   the fully symmetry-expanded cell. That expanded cell is kept as
   the verification reference for the next step.

2. **The space group is resolved by verification, not by parsing
   operation tables.** The `spaceDB` settings for the CIF's
   International Tables number are the candidate tokens
   (`<IT#>_a`, `<IT#>_b`, ...). For each candidate, the authored
   asymmetric unit is written with that token and run through
   Imago's own space-group expansion. The candidate whose regenerated
   cell reproduces ASE's full expansion (same atom count and the same
   positions, up to lattice translation and permutation, within a
   tolerance) is the answer. Because the symmetry translations are
   origin-dependent, the origin choice falls out of the match
   automatically. Imago's own symmetry engine is the judge, so no
   `spaceDB` operation set is ever parsed by this script.

3. **The skeleton is written** as the asymmetric unit plus the matched
   token. If no candidate verifies, that is a hard error: a P1
   fallback is deliberately *not* offered for a crystal. If you know
   the correct setting you can force it with `--space`, and the
   forced setting is still checked the same way.

No external symmetry library (for example spglib) is used. The CIF's
symmetry is taken as authored and checked against Imago's expansion,
never re-derived.

## Limitations

* **Partial or mixed site occupancy is refused.** A CIF whose sites
  carry fractional occupancies cannot be represented as a single
  skeleton, and the script stops with an error rather than guessing.
* **A crystal that fails verification is an error, not a P1 cell.**
  See step 3 above. Check the CIF's `_symmetry_space_group_name_H-M`
  and `_space_group_IT_number` entries, and try `--space` with each
  of the settings listed for that number in `spaceDB`.
* **Requires ASE.** The `ase` Python package must be importable from
  the interpreter running the script.

## Importable Usage

The conversion is also available as a function for use from other
Python code:

```python
from cif2skl import convert
convert("au_fcc.cif", "au_fcc.skl")
```
