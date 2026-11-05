---
title: Introduction
layout: default
nav_order: 1
---

# Introduction

Imago is the successor to the long-running OLCAO codebase. The
method keeps its established name; only the code has been rebranded.
It applies to a broad range of material systems: crystals, amorphous
solids, nanoparticles, molecules, interfaces, and grain boundaries.

Key characteristics:

- **All-electron.** There are no pseudopotentials; core electrons are
  treated explicitly. This is what makes core-level spectroscopy
  (XANES and ELNES) a native capability rather than an add-on.
- **Periodic boundary conditions.** Used throughout, even for
  non-periodic systems, which are placed in a box large enough that
  their images do not interact.
- **LCAO basis.** Wavefunctions are expanded in atomic orbitals,
  themselves built from Gaussians, so that every integral is
  analytic and every result can be decomposed atom by atom and
  orbital by orbital.
- **Orthogonalization.** Valence orbitals are orthogonalized to the
  core, shrinking the secular (eigenvalue) equation to the valence
  problem while keeping the all-electron description.

The method is described in full in the [Theory](../reference_manual/theory/index.html)
section and in the book by its authors, W.-Y. Ching and P. Rulis,
*Electronic Structure Methods for Complex Materials: The
Orthogonalized Linear Combination of Atomic Orbitals* (Oxford
University Press, 2012).

Source and development happen on
[GitHub](https://github.com/umkc-cpg/imago).

## What a calculation looks like

Every Imago calculation has the same shape. You describe a structure
in a short text file, the **skeleton**; `makeinput.py` expands it
into the engine's inputs; `imago.py` runs a self-consistent field
calculation and then whichever property you asked for; and a
post-processing script turns the result into a figure. The
[Imago Workflow](workflow.html) page walks through those stages, and
the [first calculation tutorial](../tutorials/first_calculation.html)
performs them on silicon with every command shown.

## Getting started

Imago is built from source. Follow the
[Prerequisites](../tutorials/prerequisites.html) and
[Installing Imago](../tutorials/installation.html) tutorials, or the
condensed version on the [Setup](setup.html) page. Then run the
[first calculation](../tutorials/first_calculation.html).

For the options of any individual script, see its page in the
[Script Catalog](../reference_manual/script_cat/index.html). Every
script also prints its options with `--help`.
