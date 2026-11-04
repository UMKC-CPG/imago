---
# This page is the parent of the Theory section.  `has_children`
# makes Just the Docs nest every page that names `parent: Theory`
# beneath it in the sidebar, and `nav_order` places the section
# ahead of the Script Catalog so that a reader meets the method
# before the tools.
layout: default
title: Theory
nav_order: 1
has_children: true
---

# Theory

These pages describe the OLCAO (Orthogonalized Linear Combination of
Atomic Orbitals) method as Imago implements it.  They are written for
a reader who knows basic quantum mechanics and has met density
functional theory, but who has never seen an atomic-orbital solid
state code from the inside.  Each page explains one piece of the
method and, where it matters, names the Imago option that controls
it.

The primary reference for the method is the book by its authors:

> W.-Y. Ching and P. Rulis, *Electronic Structure Methods for Complex
> Materials: The Orthogonalized Linear Combination of Atomic
> Orbitals*, Oxford University Press, 2012.
> [doi:10.1093/acprof:oso/9780199575800.001.0001](http://dx.doi.org/10.1093/acprof:oso/9780199575800.001.0001)

Chapter 3 of that book is the source for most of what follows.  The
book describes the OLCAO program as it stood in 2012.  Imago has since
changed some of what it does, most visibly in Brillouin zone
integration and in how the initial potential is chosen, and where the
two differ these pages describe Imago and say so.

Read the pages in order the first time through.  Each one builds on
the previous one.

1. **[The OLCAO method](olcao_method.html)** -- what the method is,
   where it sits among electronic structure methods, and why it is
   built the way it is.
2. **[Atomic basis functions](basis_functions.html)** -- the orbitals
   the wave function is expanded in, how they are built from
   Gaussians, and what the minimal, full, and extended basis sets
   are.
3. **[Bloch functions and the Kohn-Sham equation](kohn_sham.html)** --
   the one-electron equation Imago solves and the secular equation it
   becomes in the atomic-orbital basis.
4. **[Potential and charge density](potential_and_density.html)** --
   the site-decomposed Gaussian representation of the potential and
   density, the Gaussian transformation that makes the integrals
   analytic, and the self-consistent cycle.
5. **[Core orthogonalization](core_orthogonalization.html)** -- the
   step that removes the core states from the secular equation and
   gives the method its name.
6. **[Brillouin zone integration](brillouin_zone.html)** -- the
   k-point mesh, its reduction by symmetry, and the choice between
   Gaussian broadening and the linear analytic tetrahedron method.
7. **[Calculated properties](properties.html)** -- what each
   post-self-consistent property is as physics, and the `imago.py`
   option that computes it.
8. **[Core level spectroscopy](core_spectroscopy.html)** -- the
   supercell method behind the photo-absorption cross section, with
   an explicit core hole.
