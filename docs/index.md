---
# Front matter drives how Just the Docs treats this page.  `title`
# is the label the theme would show in the left-hand sidebar, but
# `nav_exclude: true` keeps this page out of the sidebar altogether:
# the site title in the header already links here, so a "Home" entry
# would be a duplicate.  `nav_order: 1` is kept so the page sorts
# first if the exclusion is ever lifted.  The presence of front
# matter is also what tells Jekyll to render the page through the
# active theme's layout.
title: Home
nav_order: 1
nav_exclude: true
---

# Imago

Imago is an all-electron electronic structure code implementing the
OLCAO (Orthogonalized Linear Combination of Atomic Orbitals) method.
It computes the electronic properties of crystals, amorphous solids,
nanoparticles, molecules, interfaces, and grain boundaries from a
single framework: no pseudopotentials, periodic boundary conditions
throughout, and a basis of atomic orbitals orthogonalized to the core
so the secular equation stays small.  Imago is the successor to the
long-running OLCAO codebase.  The method keeps its established name;
only the code has been rebranded.

Imago is an academic code, developed by the Computational Physics
Group at the University of Missouri--Kansas City and used for both
research and student training.  Source and development live on
[GitHub](https://github.com/umkc-cpg/imago).

## Where to start

The pages below are ordered the way a first-time user meets them.

1. **[Introduction](overview/introduction.html)** -- what the OLCAO
   method is and what makes Imago different from a pseudopotential
   code.
2. **[Setup](overview/setup.html)** -- prerequisites, cloning the
   repository, building with CMake, installing, and validating the
   install.  Imago is built from source; there are no binary
   releases yet.
3. **[Imago Workflow](overview/workflow.html)** -- the shape of a
   calculation from an atomic structure file to plotted results:
   convert the structure to a skeleton file, generate inputs with
   `makeinput.py`, run the job with `imago.py`, then post-process.

## Reference

Once you know the shape of a calculation, the Reference Manual holds
the details.

- **[Script Catalog](reference_manual/script_cat/index.html)** -- a
  searchable list of the scripts that ship with Imago.  The two every
  calculation uses are
  [makeinput.py](reference_manual/script_cat/makeinput.html), which
  turns a skeleton file into a complete job directory, and
  [imago.py](reference_manual/script_cat/imago.html), the driver
  that runs the self-consistent field and post-SCF calculations.
- **[Post Processing](reference_manual/post_processing.html)** -- the
  scripts that turn raw output into plots and into files other tools
  can read.  For most results `plotgraph.py` is the right first stop.

## Status of this documentation

This site is **incomplete** and grows over time.  Some links point at
pages that do not exist yet; those are left in place deliberately, as
reminders to the maintainers of what still needs writing.  If a page
you need is missing or a link is broken, you have three options:

- **Read the script's own help and source.**  Every script answers
  `--help` with its options and defaults, and most scripts open with
  a long comment block describing what they do and why.  This is the
  preferred route, because the source is always current.
- **Contact Dr. Paul Rulis** with questions, difficulties, or
  solutions.  Good questions become documentation.
- **Fix it yourself.**  The [Contributing](overview/contributing.html)
  page explains how the site is built and how to submit changes.

## Citation and license

If you use Imago in published work, please cite it.  Machine-readable
metadata is in the repository's `CITATION.cff`, which GitHub renders
as a "Cite this repository" panel.  A DOI will be minted at the first
tagged release; until then, cite the repository directly.

Imago is licensed under the Educational Community License, Version
2.0 (ECL-2.0).  Copyright (c) 2026 Paul Rulis.
