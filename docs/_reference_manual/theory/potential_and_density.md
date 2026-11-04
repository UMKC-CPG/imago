---
layout: default
title: Potential and charge density
parent: Theory
nav_order: 4
---

# Potential and charge density

*Source: Ching and Rulis (2012), sections 3.3, 3.4, and 12.1.2 and
Appendix B.*

## Site-decomposed functions

The most distinctive feature of the OLCAO method is how it stores the
crystal potential and the electron density.  Most electronic structure
codes tabulate them on a three-dimensional grid of points.  Imago
instead writes each of them as a sum of atom-centered functions, and
each atom-centered function as a sum of spherical Gaussians:

```
rho(r) = sum over atoms A of  rho_A(r - t_A),
rho_A(r) = sum over j of  B_j * exp(-beta_j r^2)
```

The Coulomb and exchange-correlation parts of the potential are
written the same way:

```
V_Coul(r) = sum over atoms A of  V_C(r - t_A),
V_C(r)  = -(Z_A / r) * exp(-zeta r^2)  -  sum over j of  D_j * exp(-beta_j r^2)

V_xc(r)  = sum over atoms A of  V_x(r - t_A),
V_x(r)  = sum over j of  F_j * exp(-beta_j r^2)
```

The first term of V_C is the bare nuclear attraction, screened by a
rapidly decaying Gaussian with a fixed exponent zeta (20 in Imago, the
`NUCLEAR_ALPHA` of the database entries) so that its integrals reduce
to error functions.  Everything else is a plain Gaussian.  The total
crystal potential is the sum of V_C and V_x on every site.

Two things about this representation are deliberate and carry most of
the method's efficiency.

**The same exponent set serves the density and both parts of the
potential.**  For a given element the exponents beta_j are fixed in
advance, as a geometric series from beta_min to beta_max, and only
the coefficients B_j, D_j, F_j change from one iteration to the next.
Every multi-center integral between two basis functions and a
potential Gaussian can therefore be computed once, before the
self-consistent cycle starts, and reused at every iteration by
multiplying it by the current coefficient.  Sharing the set between
density and potential also means the electron-electron interaction
needs only three-center integrals, not four-center ones.

**Spherical Gaussians, non-spherical result.**  Each atom-centered
function is spherically symmetric, but the superposition of functions
on neighboring atoms is not, and it reproduces the directional
bonding of real solids with no shape approximation.  The book's
figure 3.4 sketches this for a two-atom cell.  The known weakness is
that purely spherical site functions cannot describe the lobed
density of p or d electrons as well as they might, which slows SCF
convergence a little and is one reason forces are hard; adding
higher-angular-momentum Gaussians to the potential is a possible
future refinement that the current code does not make.

The exponent set for each element is chosen the same way as the basis
set: by experience, and by comparing calculated lattice parameters,
bulk moduli, and spectra with experiment.  Heavier atoms need more
terms and a larger beta_max; beta_min lies between about 0.08 and
0.15, smaller for alkali atoms whose potential reaches far, larger for
nearly closed-shell atoms.  For silicon the current database uses
16 Gaussians between 0.15 and 10^6.

## Why Gaussians: the Gaussian transformation

The reason for building everything from Gaussians is that integrals
of products of Gaussians on different centers have closed forms.  The
product of two Gaussians centered at A and B is a single Gaussian
centered at a point on the line between them, with an exponent equal
to the sum of the two.  This observation, the Gaussian transformation,
dates to Boys in 1950 and is the foundation of every quantum
chemistry code built on Gaussian orbitals.  OLCAO applies it to a
periodic solid.

For s-type Gaussians, every matrix element in the secular equation
reduces to one of five integrals:

1. the **overlap** of two Gaussians on sites A and B,
2. the **kinetic energy** integral between them,
3. the **three-center integral** of two basis Gaussians with a
   potential Gaussian on a third site C,
4. the **three-center integral over 1/r**, where the potential
   Gaussian is the screened nuclear term, and
5. the **momentum** integral, needed for optical matrix elements.

The first three and the fifth have closed forms.  The fourth reduces
to a one-dimensional error-function integral.  The kinetic and
momentum integrals are obtained from the overlap by differentiating
with respect to the exponent or the center.  Integrals for p, d, and
f orbitals are obtained from the s-type ones by successive
differentiation with respect to the Cartesian components of the
centers, which is where the real spherical harmonics in Cartesian
form come in.  Formulae up to f orbitals are built into Imago.

The three-center integrals are far more numerous than the two-center
ones and dominate the integral evaluation time.  For a large cell,
however, the overwhelming majority of them are negligible, because
the Gaussians on well-separated sites do not overlap, and they are
skipped.  Two cutoff exponents, one for basis function interactions
and one for electrostatic ones, set how small an interaction must be
to be dropped.  They are read from the `makeinput.py` resource-control
file and can be overridden on its command line.  This is what makes
the method scale nearly linearly with system size.

## The self-consistent cycle

Putting the pieces together, one self-consistent iteration in Imago
runs as follows.

1. The current potential coefficients D_j and F_j are combined with
   the stored integrals to assemble the Hamiltonian at each k-point.
2. The secular equation is solved and the states are occupied up to
   the Fermi level.
3. The new electron density is built from the occupied eigenvectors
   and fitted onto the atom-centered Gaussian set, giving new B_j.
   The fit is checked by comparing the integrated charge of the
   fitted density with the number of electrons in the cell.
4. The Coulomb potential follows from the density analytically.  The
   exchange-correlation potential, which is nonlinear in the density,
   is evaluated numerically on atom-centered spherical meshes and
   fitted onto the same Gaussian set, giving new D_j and F_j.
5. The new potential is mixed with the old one before the next
   iteration, to damp the oscillations that a straight replacement
   would cause.  Imago uses a fixed relaxation factor of 0.2 by
   default: the next potential is one fifth new and four fifths old.
6. Convergence is tested.  Imago measures the change in the potential
   itself: on every site's radial mesh it forms the weighted
   root-mean-square difference between the new and old potential in
   real space, takes the largest value over all sites, and declares
   convergence when that falls below the threshold.  The threshold is
   set by the `-converg` option of `makeinput.py`, with a default of
   0.0001, and the iteration count is capped at 50 by default.  The
   book describes convergence in terms of the total energy change;
   Imago reports the total energy every iteration but tests on the
   potential.

Insulators converge quickly, because a gap separates the occupied
from the unoccupied states unambiguously.  Metals converge more
slowly and need more k-points, because the Fermi level must be
located accurately for the occupation to be stable from one iteration
to the next.

## Sites, types, and the initial potential

The potential functions are attached to sites, and sites are grouped
into types.  In a high-symmetry crystal, every site that symmetry
makes equivalent shares one set of coefficients.  In a glass, every
site can be its own type, which is the most accurate representation
of a disordered potential.  In between, when many sites are similar
but not equivalent, grouping them trades a little accuracy for a
smaller fit.  Imago's structure tools assign the types.

The exponents come from the potential database shipped with Imago,
`atomicPDB`, which holds one directory per element.  The coefficients
in the database are the starting point only.  The traditional entry
for each element is the potential of the isolated atom, computed by
Imago's `atomSCF` program and fitted to the Gaussian set.  Because the
coefficients are iterated to self-consistency in every calculation,
the starting values affect how many iterations are needed, not where
they end up.  This is the sharpest contrast with pseudopotential
methods, where the pseudopotential is chosen in advance and fixes the
physics.

Imago is extending the database beyond isolated atoms.  A converged
calculation on a simple solid yields site potentials that are already
close to those of a more complex system containing the same element
in a similar environment, and starting from them can save many
iterations.  The database format now carries such entries alongside
the isolated-atom ones, and `makeinput.py` selects among them; the
`-pot`, `-potdb`, `-subpot`, and `-modpot` options control the choice
and allow the exponent set to be overridden for testing.

## Outputs worth knowing about

Because the density and potential are analytic, Imago can evaluate
either of them at any point in space without a separate calculation.
The site-decomposed density and potential are written at the end of a
self-consistent run and are the inputs to the charge density and
potential visualization tools described on the
[Post Processing](../post_processing.html) page.
