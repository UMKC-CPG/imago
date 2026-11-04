---
layout: default
title: Brillouin zone integration
parent: Theory
nav_order: 6
---

# Brillouin zone integration

*Source: Ching and Rulis (2012), section 3.6, and Imago's design
documents for the k-point mesh and the tetrahedron method, which
postdate the book.*

## Why k-points are needed

For a periodic solid the secular equation gives band energies and
eigenvectors at one wave vector k at a time.  Every quantity that
sums over occupied states, from the electron density that drives the
self-consistent cycle to the density of states and the optical
spectra, is in principle an integral over the whole Brillouin zone.
In practice the integral is approximated by solving the secular
equation on a finite set of k-points and combining the results with
weights.  The choice of that set, and of how the results are combined,
is the subject of this page.

Two facts set the scale of the problem.  The Brillouin zone shrinks
as the real space cell grows, so a large supercell needs few
k-points and a small primitive cell needs many.  And the integrand is
smooth for an insulator but has a step at the Fermi surface for a
metal, so metals need denser sampling than insulators to converge the
same quantity.

## The single k-point limit

Most applications of OLCAO are to large cells: amorphous models,
defect supercells, interfaces, biomolecules.  For these the Brillouin
zone is so small that a single k-point at its center, the Gamma
point, is sufficient, and the whole calculation collapses to one
diagonalization per iteration.  This is also the cheapest possible
k-point: at Gamma the Bloch phase factors are all unity, the overlap
and Hamiltonian matrices are real, and a real diagonalization is
several times faster than a complex one.  The same is true at the
corners of the reduced zone.  Where more than one k-point is wanted
for a large cell, those corners are the next points to add.

When Imago is given a single k-point it places it at Gamma, and the
real-arithmetic path is taken automatically.

## The mesh

For crystals with small cells, Imago samples the Brillouin zone on a
uniform Monkhorst-Pack mesh.  The mesh is described either by its
axial counts, with the `-kp` option of `makeinput.py` and its SCF and
post-SCF variants, or by a target k-point density, with `-kpd`.  The
density form is the one to prefer for systematic work: it is a single
number, in k-points per unit volume of reciprocal space, that means
the same thing for every cell, so a density that converges one
material is a sensible starting guess for another.  Given a density,
Imago chooses the axial counts itself, as uniformly as the cell shape
and its symmetry allow.

The mesh may be shifted off the origin with `-kpshift`.  A shifted
mesh avoids placing points on the high-symmetry planes where bands
are degenerate and samples the zone more evenly for a given count.
Left unspecified, Imago chooses a shift from the mesh counts and the
point group.

**Reduction to the irreducible zone.**  The point group of the
crystal maps the mesh onto itself, and the band energies at two
k-points related by a symmetry operation are equal.  Imago therefore
folds the mesh to its irreducible wedge and gives each surviving
k-point a weight equal to the number of points it stands for.  For a
cubic crystal this reduces the work by a factor approaching
forty-eight.  The space group database shipped with Imago supplies
the point operations, and `makeinput.py` writes them into the k-point
file.

The eigenvalues are invariant under the symmetry operations, but the
eigenvectors are not: a rotation that carries one atom onto another
carries that atom's orbital coefficients with it.  A quantity that is
resolved per atom, such as the bond order of a particular pair or the
partial density of states of a particular site, cannot simply be
multiplied by the weight of the irreducible k-point.  It must be
unfolded by permuting the atom labels under each operation of the
star.  Imago applies that permutation for the effective charges, bond
orders, atom-resolved partial densities of states, and atom-resolved
partial optical properties, so that a site-resolved result from the
reduced mesh equals the result from the full mesh.  Quantities
resolved only by atom type need no correction, because a type is
closed under the symmetry operations.  Partial densities of states
resolved to individual m components are the one case not handled:
those mix under rotation and are rarely wanted.

## Two integration methods

Once the eigenvalues are in hand at the mesh points, there are two
ways to turn the sum into an integral, and Imago offers both.  The
method is chosen with the `-kpint` option of `makeinput.py` and its
SCF and post-SCF variants: 0 selects Gaussian broadening, 1 selects
the linear analytic tetrahedron method.  The choice is written into
the k-point file and read by every stage.

**Gaussian broadening** (the default) replaces each eigenvalue by a
Gaussian of fixed width and sums the Gaussians with the k-point
weights.  For the occupation numbers inside the self-consistent
cycle the smearing is thermal, controlled by a width in electron
volts set by `-thermsmear`, with a default of zero.  For a density of
states or a spectrum the width is a separate parameter, with
defaults in the `makeinput.py` resource-control file.  The method is
simple and robust, and for a large cell with few k-points it is the
only sensible choice.  Its weakness is that the width is arbitrary.
It does two jobs at once, turning a discrete sum into a smooth curve
and modelling the physical broadening a measured spectrum has, and
the two cannot be told apart in the result.  A converged Gaussian
calculation is therefore a surface over width and mesh density
rather than a single number, and metals in particular need dense
meshes before the total energy settles.

**The linear analytic tetrahedron method**, after Bloechl, Jepsen,
and Andersen (1994), divides the Brillouin zone into tetrahedra with
corners on the mesh, interpolates each band linearly inside every
tetrahedron, and integrates the interpolated band analytically.  The
density of states, the occupation numbers, and the Fermi level all
follow from closed-form expressions in the four corner energies.  No
width enters.  The van Hove singularities of a crystalline density of
states come out sharp, the same mesh gives better accuracy than
broadening does, and the converged answer depends on mesh density
alone.  Imago uses the full mesh rather than the irreducible wedge
for the tetrahedra, and includes the Bloechl curvature correction
that compensates for the linear interpolation.

In Imago the tetrahedron method is available in the self-consistent
cycle, where it determines the Fermi level by requiring the
integrated electron count to equal the number of electrons, and in
the post-self-consistent properties: total and partial densities of
states, effective charges and bond orders, and the linear optical
properties.  For a tetrahedron run the optical broadening width is
still honored, but its meaning changes: it is applied afterwards as a
physical lifetime broadening and may be set to zero to see the raw
joint density of states.

**Both methods are kept on purpose.**  The two fail in different
ways, so a disagreement between them on the same mesh is a diagnostic
rather than a nuisance, and existing work done with broadening stays
reproducible.  Where a quantity has no tetrahedron path yet, that is
unfinished work rather than a decision against it.

## Converging the mesh

The quantity to watch is the total energy against mesh density.  For
an insulator it settles quickly and monotonically.  For a metal under
Gaussian broadening it oscillates, because a uniform quadrature over
a discontinuous occupation converges slowly and non-monotonically,
and the amplitude of the oscillation sets a floor below which no
convergence threshold is meaningful.  The tetrahedron method is the
standard remedy.  Imago's calculation flight tools can climb a ladder
of meshes automatically and record the density at which a material
converged, so that later calculations on similar materials can start
from it.

## Band structures

A band structure plot is not a Brillouin zone integral but a sampling
of the eigenvalues along a path of straight segments between
high-symmetry points of the zone.  For a symmetric band structure
calculation Imago ignores the mesh and instead reads the path from a
database of standard paths, one per Bravais lattice, following the
conventions of Setyawan and Curtarolo (2010).  The path is chosen
with the `-sybdpath` option of `makeinput.py`, and the calculation is
run as a post-self-consistent job with the `-sybd` option of
`imago.py`.
