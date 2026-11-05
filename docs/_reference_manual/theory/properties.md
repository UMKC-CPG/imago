---
layout: default
title: Calculated properties
parent: Theory
nav_order: 7
---

# Calculated properties

*Source: Ching and Rulis (2012), chapter 4.*

Once the self-consistent potential is in hand, every property Imago
computes is built from the same three ingredients: the band energies,
the eigenvectors, and matrix elements between pairs of states.  This
page describes each property as physics and names the `imago.py`
option that computes it.  The option page lists output files and
formats; the [Post Processing](../post_processing.html) page covers
turning them into plots.

Every property can be run in two ways.  The plain option, such as
`-dos`, runs the self-consistent cycle first and then the property.
The `-scf` prefixed form, such as `-scfdos`, reuses an existing
self-consistent potential and runs only the property.  The second is
the usual choice after the first calculation on a structure, and it
is where a denser k-point mesh or a larger basis is typically asked
for.

## Band structure and band gap

The band structure E_n(k) is the most fundamental property of a
crystal.  Plotted along high-symmetry directions of the Brillouin
zone, it shows at a glance whether a material is a metal, a
semimetal, a semiconductor, or an insulator, and whether its gap is
direct or indirect.  Imago computes it with the `-sybd` option along
the standard path for the lattice, as described on the
[Brillouin zone](brillouin_zone.html) page.

For an amorphous solid, k is not a meaningful quantum number and the
band structure loses its meaning; the density of states carries the
information instead.  For a large cell of any kind the bands fold
into the small zone and become a dense set of levels at Gamma, which
is why band structures are drawn for primitive cells.

Two cautions that the book states plainly.  Density functional
theory is a ground state theory, and the local density and gradient
approximations underestimate the band gaps of semiconductors and
insulators by thirty to fifty percent.  Imago inherits this.  And for
a metal there is no gap, only a Fermi level, and the quantity of
interest becomes the density of states there.

## Density of states

The total density of states G(E) counts states per unit energy per
cell.  Because the eigenvectors are written in atomic orbitals, Imago
can resolve it into partial densities of states (PDOS) by atom, by
element, by orbital angular momentum, by individual m component, and
by spin, with no further approximation.  The projection uses
Mulliken's population analysis: the fraction of a normalized state m
that belongs to orbital i on atom alpha is

<!-- Plain-text form of the equation below, for readers of the
     Markdown source:
       rho_{i alpha}^m = sum over (j, beta) of
                         C_{i alpha}^m * C_{j beta}^m * S_{i alpha, j beta} -->

$$
\rho^{\,m}_{i\alpha}
  = \sum_{j,\beta} C^{\,m}_{i\alpha}\, C^{\,m}_{j\beta}\,
    S_{i\alpha,\, j\beta}
$$

the state's coefficient on the orbital times its overlap with every
other orbital, which sums to one over all orbitals.  Peaks that line
up in the partial densities of two atoms show that their orbitals mix
in the same states, which is the direct signature of a bond.

The `-dos` option computes the total density of states, the
localization index described below, and a raw partial density file
that the `makePDOS.py` script resolves into whichever decomposition
is wanted.  The energy grid and broadening width defaults are in the
`makeinput.py` resource-control file, and the choice between Gaussian
broadening and the tetrahedron method is made as described on the
[Brillouin zone](brillouin_zone.html) page.

## Effective charges, bond order, and localization

Summing the Mulliken fractions over all occupied states and all
orbitals of one atom gives that atom's **effective charge** Q*.
Summing instead over the cross terms between the orbitals of two
atoms gives their **bond order**, also called the bond overlap
population, a measure of the strength of the bond between them.
Covalent bonds have larger bond orders than ionic ones.  Bond order
tends to track bond length, but since it is built from the wave
function the bond angles matter too.

The charge transfer an atom undergoes is the difference between its
effective charge and the number of valence electrons of the neutral
atom.  Effective charges are not formal charges: sodium chloride does
not come out as Na with zero electrons and Cl with eight.  Both Q*
and bond order are basis-dependent, so values are comparable only
between calculations with the same method and basis.  For this reason
they are always computed in the minimal basis, which is the most
localized and the most faithful to Mulliken's scheme.  The `-bond`
option computes both and the `makeBOND.py` script post-processes the
result.

The **localization index** of a state is the sum of the squares of
its Mulliken fractions.  It runs from one, for a state confined to a
single orbital, down to the reciprocal of the number of states, for
one spread evenly over the whole cell.  It is the inverse of the
participation ratio, and in a large amorphous model it gives a
realistic measure of how localized each one-electron state is across
the whole spectrum.  The `-dos` option writes it alongside the total
density of states.

## Spin polarization

Magnetic materials need the local spin density approximation, in
which the exchange-correlation potential depends on the spin-up and
spin-down densities separately.  The spin-polarized functionals in
Imago, von Barth-Hedin and Ceperley-Alder among them, are selected by
their codes with the `-xccode` option of `makeinput.py`.  The
self-consistent cycle then solves the Kohn-Sham equation once for
each spin, starting from an initial rigid splitting of the potential
whose size is set in the resource-control file, and merges the two
densities.  Convergence is slower than for the unpolarized case,
because the Fermi level must be located for two coupled spectra.

In a spin-polarized run every property splits into up and down
components, and the difference of the spin-resolved effective
charges on an atom is its spin magnetic moment.  Because that moment
is a Mulliken quantity, it involves no atomic sphere radius or other
adjustable parameter.  The calculation is collinear: every moment
points along one axis.  Non-collinear magnetism and spin-orbit
coupling are not implemented.

## Relativistic corrections

Core electrons of heavy elements move fast enough that the Dirac
equation, not the Schroedinger equation, describes them.  Reducing
the Dirac equation to its large component adds three terms to the
Hamiltonian: the mass-velocity and Darwin terms, which are scalars
and can be folded into the one-electron potential, and the spin-orbit
term, which couples the two spins and is generally an order of
magnitude smaller.  The book describes both the scalar relativistic
correction and a parameterized spin-orbit treatment that doubles the
matrix dimension.

Imago carries scalar relativistic variants of every functional in its
exchange-correlation table, and the atomic basis program can produce
scalar relativistic orbitals.  The relativistic input path of
`makeinput.py`, the `-rel` option, is not yet functional, and
spin-orbit coupling is not implemented.  Treat relativistic
calculations as unfinished work.

## Linear optical properties

The optical response of a solid in the random phase approximation
is an interband absorption: an electron in an occupied state l is
promoted to an unoccupied state n by a photon of energy equal to the
difference, with a probability given by the momentum matrix element
between the two.  The real part of the optical conductivity is

<!-- Plain-text form of the equation below, for readers of the
     Markdown source:
       sigma_1(E) ~ (1 / E) sum over (l occupied, n unoccupied, k) of
                    |<n| p |l>|^2  delta( E_n(k) - E_l(k) - E )     -->

$$
\sigma_1(E) \;\propto\; \frac{1}{E}
  \sum_{l\,\mathrm{occ}}\ \sum_{n\,\mathrm{unocc}}\ \sum_{\mathbf k}
  \bigl| \langle n | \mathbf p | l \rangle \bigr|^2\,
  \delta\bigl( E_n(\mathbf k) - E_l(\mathbf k) - E \bigr)
$$

and the imaginary dielectric function follows from it,
epsilon_2 = 4 pi sigma_1 / omega.  The momentum matrix elements are
two-center integrals between Gaussians, evaluated exactly as
described on the
[potential and charge density](potential_and_density.html) page.
The energy-conserving delta function is handled by whichever
Brillouin zone integration method the run uses.

The real dielectric function epsilon_1 is obtained from epsilon_2 by
the Kramers-Kronig transform over a finite energy range, which is one
reason the extended basis is wanted: the transform needs epsilon_2
well above the energies of interest.  From epsilon_1 and epsilon_2
every other linear optical constant follows: the refractive index and
extinction coefficient, the absorption coefficient, the reflectivity,
and the energy loss function, whose peak is the plasmon energy.

Imago resolves each quantity into x, y, and z components plus their
average, which is the isotropic polycrystalline value.  The `-optc`
option computes the full set, and a partial optical calculation
decomposes epsilon_2 by atom or type pairs for the `processPOPTC.py`
script.  The book's remark on intraband transitions applies: in a
large supercell at Gamma the bands fold so densely that the
low-frequency Drude response of a metal is captured by the same
interband formula.

## Conductivity in metals

For metals the quantity of interest is the energy-dependent
conductivity function sigma(E) of the Kubo-Greenwood formula, a
generalization of the interband conductivity in which both states
lie at the energy E.  Only states near the Fermi level matter for
transport.  The zero-temperature d.c. conductivity is sigma at the
Fermi level, its reciprocal the residual resistivity, and these
account for elastic scattering only, not phonons.  In a finite
supercell the spectrum is discrete and each level is given a
Gaussian of unit area; the larger the cell, the narrower the Gaussian
can be.  The `-sige` option computes this curve.

## Nonlinear optical properties

The polarization of a medium is a power series in the applied field,
and the second and third order terms are the nonlinear
susceptibilities responsible for second and third harmonic
generation.  The second order term vanishes in centrosymmetric
crystals.  The calculation resembles the linear one but also needs
momentum matrix elements between pairs of unoccupied states, because
the virtual transitions through the conduction band dominate.  The
`-nlop` option computes the imaginary part of the susceptibility and
recovers the real part by Kramers-Kronig transform, reusing the
linear optical machinery.

## Total energy, forces, and geometry

The total energy is computed every iteration and is the quantity used
to compare structures of the same composition and to converge the
k-point mesh.  Forces are harder.  Because the basis functions ride
on the nuclei, the Hellmann-Feynman force carries a large Pulay
correction that is difficult to evaluate in a local orbital method.
The book's route was finite differences of the total energy, driven
by an external optimizer, which is practical for simple crystals with
few parameters and not for complex ones.  Imago has an experimental
direct force calculation, selected with `-force`, whose results and
output format are preliminary.  The practical recommendation remains
the one the book gives: relax the structure in a plane wave code and
bring it to Imago for the electronic properties.

## Other analyses

Three further options fall outside the book's chapter.  The `-dimo`
option computes the electric dipole moment and polarization of the
cell from the electron and nuclear charge distributions.  The
`-field` option evaluates the wave function, charge density, and
potential on a real space grid, in absolute form and as the
difference from superposed neutral atoms, for visualization.  The
`-loen` option computes, for every atom, a set of rotationally
invariant bispectrum components describing the geometry of its
neighbor shell; this is not an electronic property but a fingerprint
of the local environment used to compare and classify sites.
