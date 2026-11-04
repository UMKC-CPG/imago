---
layout: default
title: Bloch functions and the Kohn-Sham equation
parent: Theory
nav_order: 3
---

# Bloch functions and the Kohn-Sham equation

*Source: Ching and Rulis (2012), sections 3.2 and 3.4 (first part).*

## The one-electron equation

Imago works within density functional theory.  The many-electron
problem is replaced by a set of one-electron Kohn-Sham equations, one
for each state, in which every electron moves in an effective
potential built from the total electron density.  In atomic units:

<!-- Plain-text form of the equation below, for readers of the
     Markdown source:
       [ -grad^2 + V_en(r) + V_ee(r) + V_xc[rho](r) ] Psi_nk(r)
           = E_n(k) Psi_nk(r)                                      -->

$$
\Bigl[ -\nabla^2 + V_{e\text{-}n}(\mathbf r) + V_{e\text{-}e}(\mathbf r)
       + V_{xc}[\rho](\mathbf r) \Bigr] \Psi_{n\mathbf k}(\mathbf r)
  = E_n(\mathbf k)\, \Psi_{n\mathbf k}(\mathbf r)
$$

The first term is the kinetic energy.  The three potentials are the
attraction of the electron to the nuclei, the classical Coulomb
repulsion of the electron from the whole electron density, and the
exchange-correlation potential, which carries everything about the
many-electron problem that the first two terms leave out.  The
electron density that enters the potentials is built from the
occupied solutions:

<!-- Plain-text form of the equation below, for readers of the
     Markdown source:
       rho(r) = sum over occupied (n, k) of  |Psi_nk(r)|^2          -->

$$
\rho(\mathbf r) = \sum_{n,\mathbf k}^{\mathrm{occ}}
                  \bigl| \Psi_{n\mathbf k}(\mathbf r) \bigr|^2
$$

The potential depends on the density, and the density depends on the
solutions, so the equation is solved by iteration: guess a potential,
solve, build the density, rebuild the potential, and repeat until
nothing changes.  That cycle is described on the
[potential and charge density](potential_and_density.html) page.
This page is about what one pass of the cycle solves.

## Exchange and correlation

Imago evaluates the exchange-correlation potential in the local
density approximation (LDA) by default, with the Wigner interpolation
formula for correlation.  Several other functionals are available,
including Ceperley-Alder and Hedin-Lundqvist parameterizations of the
LDA, their spin-polarized counterparts, the PBE generalized gradient
approximation, and scalar relativistic variants of each.  The
functional is selected by number with the `-xccode` option of
`makeinput.py`; the full table of codes is on the
[makeinput.py](../script_cat/makeinput.html) page and in the
`xc_code.dat` data file.

Unlike the Coulomb terms, the exchange-correlation potential is a
nonlinear function of the density and has no closed form in the
Gaussian representation Imago uses for everything else.  It is
therefore evaluated numerically on a set of atom-centered spherical
meshes and then fitted back onto the Gaussian representation.  The
`-xcmesh` options of `makeinput.py` control the construction of those
meshes.

## The secular equation

Substituting the Bloch sum expansion of the wave function (see the
[basis functions](basis_functions.html) page) into the Kohn-Sham
equation and projecting onto each basis function turns the
differential equation into a matrix equation at each wave vector k:

<!-- Plain-text form of the equation below, for readers of the
     Markdown source:
       | H(k) - E(k) S(k) | = 0                                     -->

$$
\bigl| H(\mathbf k) - E(\mathbf k)\, S(\mathbf k) \bigr| = 0
$$

This is the secular equation, a generalized eigenvalue problem.  S is
the overlap matrix, whose elements are the overlaps of pairs of Bloch
sums, and H is the Hamiltonian matrix, whose elements are the
expectation values of the Kohn-Sham operator between pairs of Bloch
sums.  Both are square matrices of dimension equal to the number of
basis functions in the cell.  The overlap matrix is not the identity
because atomic orbitals on different atoms are not orthogonal to one
another, which is the ordinary situation in any LCAO method.

Each matrix element is a lattice sum: the orbital on atom gamma in
the home cell against the orbital on atom delta in every periodic
image, each image weighted by its Bloch phase factor.  Because the
orbitals are localized, the sum converges after a few neighboring
cells, and for a large supercell only the nearest images matter.  For
a large supercell the Brillouin zone is correspondingly small and a
single k-point at the zone center usually suffices, which keeps the
total number of integrals manageable even for thousands of atoms.

Imago builds the overlap matrix, the kinetic energy and nuclear
attraction parts of the Hamiltonian, and the overlap of each basis
function pair with each Gaussian term of the potential once, before
the self-consistent cycle begins, because none of these depend on the
potential coefficients.  Each iteration then assembles the Hamiltonian
by combining the stored integrals with the current potential
coefficients, which is cheap, and diagonalizes it.

## Real and complex arithmetic

At a general k-point the Bloch phase factors are complex and so are H
and S.  At the zone center, and at the corners of the reduced
Brillouin zone where symmetry makes the phase factors real, the
matrices are real.  This matters for cost: a complex diagonalization
takes several times as long as a real one of the same dimension.  The
use of real spherical harmonics in the basis is what makes the real
case available at all.

## Solving it

The secular equation is solved by direct diagonalization, giving the
full set of eigenvalues E_n(k) and eigenvectors C^n(k) at each
k-point.  Before that happens the core orbitals are removed from H
and S by [orthogonalization](core_orthogonalization.html), which
reduces the dimension to the number of valence basis functions.

The eigenvalues are the band energies.  Occupying them from the
bottom up with the right number of electrons fixes the Fermi level,
and integrating over the Brillouin zone gives the density needed for
the next iteration.  How that integration is done, by Gaussian
broadening or by the linear analytic tetrahedron method, is a user
choice in Imago and is described on its own page.

The eigenvectors are the quantities every post-self-consistent
property is built from.  Because the basis functions are atomic
orbitals, the coefficient of a state on a given orbital has a direct
chemical meaning, and partial densities of states, effective charges,
bond orders, and optical matrix elements follow from the eigenvectors
without further approximation.

## Total energy

The total energy of the cell follows from the band energies with the
usual density functional corrections for double counting:

<!-- Plain-text form of the equation below, for readers of the
     Markdown source:
       E_total = sum over occupied (n, k) of E_n(k)
               + integral of rho(r) [ eps_xc - V_xc - V_ee / 2 ] dr
               + (1/2) sum over atom pairs of  Z_a Z_b / |R_a - R_b| -->

$$
E_{\mathrm{total}}
  = \sum_{n,\mathbf k}^{\mathrm{occ}} E_n(\mathbf k)
  + \int \rho(\mathbf r)
      \Bigl[ \varepsilon_{xc} - V_{xc} - \tfrac12 V_{e\text{-}e} \Bigr]
      d\mathbf r
  + \frac12 \sum_{a \ne b} \frac{Z_a Z_b}{|\mathbf R_a - \mathbf R_b|}
$$

The first term is the sum of occupied one-electron energies, the
second removes the electron-electron interaction and
exchange-correlation energy that the band energies count twice, and
the third is the nuclear repulsion.  Imago reports the kinetic,
electrostatic, exchange-correlation, and total energies at every
iteration.  The total energy is used to compare the stability of
structures with the same composition and as the basis of the
experimental force calculation.
