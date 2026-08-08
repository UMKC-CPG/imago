---
title: Introduction
layout: default
nav_order: 1
---

# Introduction

Welcome to the Imago wiki! Imago is the successor to the long-running OLCAO codebase. The method
keeps its established name; only the code has been rebranded. It applies
to a broad range of material systems -- crystals, amorphous solids,
nanoparticles, molecules, interfaces, and grain boundaries.

Key characteristics:

- **All-electron** -- no pseudopotentials; core electrons are treated
  explicitly.
- **Periodic boundary conditions** -- used throughout, even for
  non-periodic systems.
- **LCAO basis** -- wavefunctions expressed as linear combinations of
  atomic orbitals.
- **Orthogonalization** -- valence orbitals are orthogonalized to the
  core, shrinking the secular (eigenvalue) equation.

Source and development happen on
[GitHub](https://github.com/umkc-cpg/imago). 

## Getting Started

It is recommended that you build and run Imago from source at the moment. Thus, you should clone the github
repository and follow the installation instructions [here](setup.html).

After getting Imago installed, familiarize yourself with the contents of [Imago Workflow](workflow.html) for
and overview on how to use Imago. For script-specific help see the individual script documentations under the
[Script Catalog](../reference_manual/script_cat/index.html). 

See the Reference Manual section for other documentation about concepts and objects universal to the Imago program.

