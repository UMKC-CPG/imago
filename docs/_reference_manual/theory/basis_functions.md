---
layout: default
title: Atomic basis functions
parent: Theory
nav_order: 2
---

# Atomic basis functions

*Source: Ching and Rulis (2012), sections 3.1 and 12.1.1 and
Appendix A.*

## The expansion

Imago describes a solid strictly as a band structure problem, so every
electronic state carries a band index n and a wave vector k, even when
the "solid" is a single molecule in a large cell.  The wave function
of a state is expanded in Bloch sums:

```
Psi_nk(r) = sum over (i, gamma) of  C_{i gamma}^n(k) * b_{i gamma}(k, r)
```

Here gamma runs over the atoms in the cell and i over the orbitals on
each atom.  A Bloch sum is the same atomic orbital placed on every
periodic image of one atom, each copy carrying the phase factor
Bloch's theorem demands:

```
b_{i gamma}(k, r) = (1 / sqrt N) * sum over lattice vectors R of
                    exp(i k . R) * u_i(r - R - t_gamma)
```

where t_gamma is the position of atom gamma in the cell.  The
coefficients C are what the self-consistent calculation determines.
Everything else on this page concerns the atomic orbitals u_i, which
are fixed before the calculation begins.

## One orbital

Each atomic orbital is a radial function times a real spherical
harmonic:

```
u_i(r) = [ sum over j of  A_j * r^l * exp(-alpha_j r^2) ] * Y_lm(theta, phi)
```

The index i stands for the full set of quantum numbers n, l, m.  The
radial part is a linear combination of Gaussian-type orbitals (GTOs),
each with a decay exponent alpha_j.  The angular part uses real
spherical harmonics up to l = 3, so s, p, d, and f orbitals are
available.  Real rather than complex harmonics are used so that the
Hamiltonian is real at k-points where symmetry allows it.

The reason for Gaussians, rather than the exponential (Slater-type)
functions that actually solve the hydrogen atom, is that products of
Gaussians on different centers are themselves Gaussians.  That single
identity is what lets every integral in the Hamiltonian be evaluated
in closed form; see the
[potential and charge density](potential_and_density.html) page.  The
price is that many Gaussians are needed to reproduce the cusp of a
true orbital at the nucleus and its exponential tail far away.

## The exponent set

The exponents alpha_j are not optimized orbital by orbital.  Instead
each element is given a single set of N exponents, spaced as a
geometric series from alpha_min to alpha_max, and **every orbital of
that element uses the same set**.  Orbitals differ only in their
coefficients A_j, in how many of the N Gaussians they use, and in their
spherical harmonic.

This is the first of the method's deliberate economies.  With one
exponent set per element, an integral between two Gaussians depends
only on the pair of elements and the geometry, not on which orbitals
the Gaussians belong to.  The same integral serves every orbital pair
built from that exponent pair, and the total number of distinct
integrals for a thousand-atom cell stays manageable.

The size of the set and its range are chosen per element from long
experience.  The number of terms is typically between 16 and 30.  The
smallest exponent is near 0.12 for nearly every element: a smaller
value gives over-extended, nearly linearly dependent orbitals and a
longer interaction range, while a larger one localizes the orbitals
too tightly for the solid.  The largest exponent grows with atomic
number, from around 10^5 for light atoms to 10^9 for the heaviest,
because the 1s orbital of a heavy atom is very compact.  For silicon,
the current database entry uses 21 Gaussians between 0.12 and
2 x 10^5.

The coefficients A_j come from a self-consistent calculation on the
isolated atom, done in the same density functional theory with the
same Gaussian set.  Imago's atomic program, `atomSCF`, performs that
calculation; its numerical orbitals are then fitted to the Gaussian
set, and the results are stored in the basis database.  The isolated
atom's orbitals are more extended than they would be in a solid, but
this does not matter much: the solid state wave function readjusts
itself through the coefficients C during the self-consistent cycle,
and the basis only has to span the right space.

## Core, valence, and the three basis sets

The orbitals of an element are divided into two groups.  **Core**
orbitals are the deep, chemically inert ones, typically everything
more than about 30 eV below the highest occupied state.  **Valence**
orbitals, in Imago's usage, are all the rest: the occupied valence
shell and any unoccupied orbitals included for variational freedom.
The core orbitals are present in the calculation but are removed from
the eigenvalue problem by
[orthogonalization](core_orthogonalization.html).

Three basis sets of increasing size are defined for every element:

- **Minimal basis (MB)** -- the core plus the orbitals of the
  valence shell, occupied or not.  For iron that is 1s, 2s, 2p, 3s,
  3p as core and 4s, 4p, 3d as valence.  The minimal basis is the
  right choice for Mulliken effective charges and bond orders, which
  are most consistent when no fully unoccupied orbitals are present
  to spread charge into.
- **Full basis (FB)** -- the minimal basis plus one more shell of
  unoccupied orbitals.  For iron that adds 5s, 5p, and 5d.  The full
  basis is sufficient for most purposes and is Imago's default for
  the self-consistent calculation.
- **Extended basis (EB)** -- the full basis plus a further shell.  It
  is needed when high-lying unoccupied states matter, as in optical
  spectra and core level spectroscopy, where transitions reach far
  above the Fermi level.

Imago's `imago.py` driver chooses the basis for the self-consistent
and post-self-consistent stages by calculation type, and the choice
can be overridden with its `-scf` and `-pscf` options; see the
[imago.py](../script_cat/imago.html) page for the defaults per
property.

Larger is not always better.  Adding orbitals increases variational
freedom and reaches higher energies, but an atomic orbital basis has
no single parameter that converges it systematically the way a plane
wave cutoff does, and past a point the added orbitals become nearly
linearly dependent and destabilize the calculation.  The default
definitions have been tested for decades on the first four rows of
the periodic table and rarely need changing.

## The basis database

The basis functions for 103 elements ship with Imago in the
`atomicBDB` database, one directory per element named by its
lowercase symbol.  Each directory holds the atomic calculation that
produced the orbitals and one or more contraction files that define
the basis sets.  A contraction file records, per angular momentum,
how many orbitals are core and how many are valence in each of the
minimal, full, and extended sets; the number of Gaussians and the
exponent range; and which of the Gaussians each angular momentum
actually uses.  Silicon's entry, for example, places 1s, 2s, and 2p
in the core, 3s and 3p in the minimal valence set, adds a d orbital
for the full basis, and uses the same orbitals again for the extended
basis.

The database can be redirected with the `-basisdb` option of
`makeinput.py`, and an alternative contraction file can be chosen for
one element or species with `-subbasis`.  The `-basisVis` option
writes the radial functions out for plotting.  Including the core
orbitals in the eigenvalue problem instead of orthogonalizing them
away is a separate choice, made with `-nocore`, and is described on
the [core orthogonalization](core_orthogonalization.html) page.
