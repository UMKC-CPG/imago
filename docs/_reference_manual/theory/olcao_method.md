---
layout: default
title: The OLCAO method
parent: Theory
nav_order: 1
---

# The OLCAO method

*Source: Ching and Rulis (2012), sections 1.4, 2.4, 2.5, and 3.7.*

## What it is

OLCAO stands for Orthogonalized Linear Combination of Atomic Orbitals.
It is a method for solving the Kohn-Sham equation of density
functional theory for a periodic arrangement of atoms, and Imago is
the program that implements it.  The method is a direct extension of
the classical LCAO (Linear Combination of Atomic Orbitals) method of
band theory.  Its one distinguishing step, the orthogonalization of
the valence orbitals against the core orbitals so that the core can
be removed from the eigenvalue problem, was borrowed from the
orthogonalized plane wave (OPW) method and gave OLCAO its name.

Although OLCAO is classified as a solid state method, it is better
thought of as a hybrid of quantum chemistry and solid state physics,
sharing strengths and limitations with both.  Like a quantum
chemistry code it expands the wave function in atomic orbitals built
from Gaussians, evaluates every integral analytically in real space,
and treats all electrons explicitly.  Like a solid state code it uses
Bloch's theorem, so the system is always an infinitely repeated cell
and the results are naturally expressed in reciprocal space, as bands
and densities of states.

## What it is built for

The method was designed in 1975 to study amorphous solids represented
by large periodic models.  That origin explains most of its design
choices.  If the periodic model is large compared with the mean free
path of an electron, the artificial periodicity does no harm, and a
disordered solid can be treated with the same machinery as a crystal.
Everything about the method is arranged to keep that machinery
affordable when the cell contains hundreds or thousands of atoms with
no symmetry at all:

- **Element-specific, predetermined basis functions.**  Every atom of
  a given element uses the same basis, read from a database, so the
  number of distinct integrals to evaluate does not grow with the
  number of chemically distinct sites.
- **One set of Gaussian exponents per element.**  All the orbitals of
  an element share the same exponents and differ only in their
  coefficients.  The potential and the charge density of that element
  share a second set of exponents between them.  Integrals computed
  once for an exponent pair are reused across every orbital and every
  self-consistent iteration.
- **Analytic integrals.**  Because basis, potential, and density are
  all Gaussians, every multi-center integral in the Hamiltonian has a
  closed form or reduces to the error function.  There is no real
  space grid for the potential and no numerical quadrature for the
  matrix elements.
- **Core orthogonalization.**  Core states are removed from the
  secular equation before it is diagonalized, which matters most for
  heavy elements, where the core is large and chemically inert.

Together these make OLCAO effectively an order-N method for large
systems: beyond a certain separation, nearly all multi-center
integrals are negligible and are skipped.

The same machinery serves crystals, defect supercells, surfaces,
interfaces, grain boundaries, glasses, liquids, nanoparticles, and
molecules.  A molecule or an isolated cluster is handled by placing it
in a cell large enough that it does not interact with its images.

## Where it sits among methods

Two families of methods dominate electronic structure work on solids.
Plane wave pseudopotential codes expand the wave function in plane
waves, replace the core with a pseudopotential, and converge
systematically by raising a single energy cutoff.  All-electron
methods such as augmented plane wave (APW) and linearized
muffin-tin orbital (LMTO) codes keep the core but partition space
around each atom.  OLCAO is an all-electron method with neither a
pseudopotential nor a space partition.  Its basis is a fixed set of
atomic orbitals, so there is no single parameter that converges it
systematically; accuracy is adjusted instead by the choice of basis
set, as described on the [basis functions](basis_functions.html)
page.

What that trade buys is interpretability.  Because the wave function
is written in atomic orbitals, the chemistry is visible directly in
the eigenvectors.  Partial densities of states, Mulliken effective
charges, bond orders, and the localization of each state all fall out
of the solution without further approximation.  The matrix to be
diagonalized stays small enough, even with an extended basis, that
the full spectrum of eigenvalues and eigenvectors at a k-point can be
obtained by direct diagonalization, which is what makes the
high-energy unoccupied states needed for optical spectra accessible.

The two families are complementary rather than competing.  A common
pattern is to relax a structure with a plane wave code, where forces
are cheap and accurate, and then compute its electronic properties
and spectra with OLCAO.  Imago's structure converters exist for
exactly this hand-off; see the [Workflow](../../overview/workflow.html)
page.

## A short history

The OLCAO method was created by W.-Y. Ching and C. C. Lin in 1975,
tested on crystalline silicon, and immediately applied to models of
amorphous silicon, then amorphous silica and metallic glasses.  Those
early calculations were not self-consistent.  Full self-consistency
was added in the 1980s, after which the method was applied over the
following decades, almost entirely by the electronic structure group
at the University of Missouri-Kansas City, to permanent magnets,
oxide superconductors, fullerene crystals, nonlinear optical and
laser host crystals, spinel nitrides, biomolecules such as vitamin
B12, grain boundaries, and core level spectroscopy.

The book that documents the method closed with a stated need: a
complete, distributable, parallelized package with a manual, so that
the method could be used outside the handful of groups that built it.
Imago is that package.  It is the successor to the OLCAO codebase; the
method keeps its name and only the program has been rebranded.

## What the method does not do well

The book is candid about the limits, and they still hold.  Density
functional theory in the local density or generalized gradient
approximation does not describe strongly correlated systems or heavy
elements with f electrons in the valence shell adequately, and OLCAO
inherits that limit.  Forces are harder than in plane wave codes
because the basis functions move with the atoms; Imago's force
calculation is experimental.  And because the basis is fixed per
element, a system that puts an element in an unusual electronic
configuration can occasionally need a basis tuned for it, which is
why the basis database can be overridden.
