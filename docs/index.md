---
# Front matter drives where Just the Docs places this page.  `title`
# is the label shown in the left-hand navigation sidebar, and
# `nav_order: 1` pins this landing page to the very top of that list.
# Its presence is also what tells Jekyll to render the page through the
# active theme's layout.
title: Home
nav_order: 1
nav_exclude: true
---

<!-- The Imago logo sits to the right of the heading and opening
     paragraphs.  Width is fixed so the SVG, which declares a
     1024-pixel natural size, does not dominate the page. -->
<img src="assets/images/imago_logo.svg" alt="Imago logo" width="180"
     style="float: right; margin: 0 0 1rem 1.5rem;">

# Imago

Imago is an all-electron electronic structure code implementing the
OLCAO method, the Orthogonalized Linear Combination of Atomic
Orbitals. It computes the electronic ground state of crystals,
glasses, molecules, nanoparticles, interfaces and grain boundaries,
and from it the properties that experiments measure: densities of
states, band structures, bond orders and charges, optical spectra,
and core-level spectra. It is developed by the Computational Physics
Group at the University of Missouri-Kansas City and released under
the Educational Community License 2.0.

This site is the user documentation. The source, the issue tracker
and the design documents are on
[GitHub](https://github.com/umkc-cpg/imago).

<div style="clear: both;"></div>

## Where to start

| If you want to...                          | Go to |
|:-------------------------------------------|:------|
|Understand what Imago does and how a calculation flows | [Introduction](overview/introduction.html) and [Imago Workflow](overview/workflow.html) |
|Install it and run your first calculation   | The [Tutorials](tutorials/index.html), which take you from a bare machine to a plotted density of states and band structure |
|Look up an option of a script               | The [Script Catalog](reference_manual/script_cat/index.html), one page per script, searchable from the box above |
|Understand the method behind the numbers    | The [Theory](reference_manual/theory/index.html) section, which follows the OLCAO book chapter by chapter and names the Imago option that controls each piece |
|Turn output files into figures              | The [post-processing section](reference_manual/script_cat/index.html#post-processing) of the Script Catalog |
|Fix or extend the documentation             | [Contributing](overview/contributing.html) |

## Getting help

If you cannot find what you need, the scripts themselves are the
next place to look: every one carries a full description in its
docstring and prints its options with `--help`. Beyond that, contact
Dr. Paul Rulis at the University of Missouri-Kansas City. When
reporting a problem, include the output of `env | grep IMAGO`, the
error text, and the `runtime` file from the job directory; the
installation tutorial explains why those three are what is needed.

## Developed by

<!-- Affiliation marks on one row.  The sizes are set with inline
     styles rather than height attributes because the theme's
     stylesheet gives every image `height: auto`, which overrides the
     attribute and lets an SVG render at its natural size. -->
<p style="display: flex; align-items: center; gap: 2.5rem;
          flex-wrap: wrap; margin: 1rem 0;">
  <img src="assets/images/cpg_logo.svg"
       alt="Computational Physics Group"
       style="height: 64px; width: auto;">
  <img src="assets/images/umkc_sse_wordmark.svg"
       alt="UMKC School of Science and Engineering"
       style="height: 48px; width: auto;">
</p>

The Computational Physics Group in the School of Science and
Engineering at the University of Missouri-Kansas City.
