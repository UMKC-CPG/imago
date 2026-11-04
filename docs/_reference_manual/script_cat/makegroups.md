---
layout: default
title: makegroups.py
script: makegroups.py
description: 'Groups the atoms of a non-crystalline (P1) skeleton into species by their bispectrum local-environment fingerprints. Runs makeinput and an imago.py -loen pass to compute the descriptors inside the engine, clusters atoms of each element by fingerprint distance, and rewrites the skeleton with explicit per-element species tags.'
parent: Script Catalog
---

# makegroups

Groups the atoms of a non-crystalline structure into **species** by
the similarity of their local environments. In a crystal, species
follow from symmetry: every atom on the same Wyckoff site is the same
species. In an amorphous model there is no symmetry, yet atoms still
fall into a handful of environment types (a silicon bonded to four
oxygens, a silicon with one dangling bond, and so on), and Imago can
give each type its own potential only if the skeleton tags them as
distinct species. This script produces those tags.

The descriptor is the **bispectrum**, which measures the angular
arrangement of an atom's neighbours. It is computed inside the Imago
engine's local-environment path (`-loen`), not in Python, so grouping
by it is a short, explicit sequence of ordinary program runs driven
from outside `makeinput.py`:

1. **First makeinput.** Write a provisional `imago.dat` from the
   skeleton with no environment grouping, so every atom keeps the
   species the skeleton already declares. The `LOEN_INPUT_DATA`
   block of that file carries the chosen bispectrum parameters.
2. **Run loen.** Invoke `imago.py -loen -scf no` against that input.
   The `-scf no` skips the self-consistent cycle the bispectrum does
   not need. The run writes a self-describing `fort.21` with one row
   per potential site: site number, element, species, type, then
   the descriptor components.
3. **Bucket and rewrite.** Read `fort.21`, group the atoms of each
   element by how close their fingerprints are (leader clustering
   against a similarity floor), and rewrite the skeleton with
   explicit per-element species tags: `Si1, Si2, ...` then
   `O1, O2, ...`, each element restarting at 1. `makeinput.py` then
   reads those tags on a later run like any other explicit species
   assignment.

**Crystals are refused.** Grouping applies only to P1 systems. A
symmetry-bearing skeleton is rejected up front, because rewriting its
species in P1 would force it out of its space group, and the
Brillouin-zone k-point folding samples the irreducible wedge *using
that space group*. Dropping it would silently corrupt every
k-point-dependent result.

The script is also importable: `build_initial_potentials.py` calls
`group_by_bispectrum()` for every non-crystalline reference solid.

## Usage

**Commandline:**

```
makegroups.py SKELETON [-twoj1 N] [-twoj2 N] [-maxneigh N]
              [-cutoff R] [-anglesqueeze W] [-floor D]
              [-workdir DIR] [-keepwork] [--help]
```

The skeleton is rewritten **in place**; the original is kept beside
it as `<skeleton>.orig`.

**Variable Table:**

| Variable   | Accepted Values     | Description |
|:----------:|:-------------------:|:------------|
|SKELETON    |String; file path    | Path to the P1 imago skeleton to regroup. |
|N           |integer              | Angular-momentum parameter or neighbour cap. |
|R           |float; Bohr          | Radial neighbour cutoff. |
|W           |float                | Angular weighting. |
|D           |float                | Largest fingerprint distance at which two atoms still share a species. |
|DIR         |String; directory path | Scratch directory for the makeinput and loen runs. |

**Option Table:**

| Option          | Defaults                                  | Description |
|:---------------:|:-----------------------------------------:|:------------|
|-twoj1           |4                                          | Larger angular-momentum parameter 2j1 of the bispectrum. |
|-twoj2           |4                                          | Smaller angular-momentum parameter 2j2. The descriptor has 2j2 + 1 components. |
|-maxneigh        |50                                         | Per-site neighbour-list cap. |
|-cutoff          |9.0                                        | Radial neighbour cutoff in Bohr. |
|-anglesqueeze    |0.85                                       | Angular weighting (angleSqueeze). |
|-floor           |The bispectrum matcher's built-in 0.10     | Similarity floor: the largest fingerprint distance at which two atoms still share a species. Lower it to split environments more finely. |
|-workdir         |`makegroups_work` beside the skeleton      | Scratch directory for the makeinput and loen runs. |
|-keepwork        |Remove it                                  | Keep the scratch directory for inspection. |
|--help           |N/A                                        | Print usage information and exit. |

## Output

The rewritten skeleton, with the original saved as
`<skeleton>.orig`, and a report of how many species each element was
split into. Inspect the result by comparing the two files; if an
element has split into more species than its chemistry justifies,
raise `-floor`, and if distinct environments have merged, lower it.
