---
layout: default
title: Core level spectroscopy
parent: Theory
nav_order: 8
---

# Core level spectroscopy

*Source: Ching and Rulis (2012), section 11.1.*

## What is measured

X-ray absorption near edge structure (XANES) and electron energy loss
near edge structure (ELNES) both probe the same thing: the
excitation of an electron from a deep core level of one atom into
the unoccupied states of the solid.  The spectrum is a function of
the energy transferred, starting at the edge onset where the core
electron just reaches the lowest empty state and extending some tens
of electron volts above it.  Because the core level is atomic and
sharp, the shape of the spectrum is set by the unoccupied states as
seen from that one atom, which makes the technique a site-specific
and element-specific probe of local bonding.

For small momentum transfer and energies well above the plasma
frequency the cross section reduces, by Fermi's golden rule in the
dipole approximation, to a sum over final states f of the squared
dipole matrix element between the core state g and f, times a delta
function at the transition energy.  The dipole matrix element is
equivalent to the momentum matrix element of the optical calculation
and is computed the same way.  Dipole selection rules follow
automatically from the symmetry of the wave functions: a 1s core
state reaches only p-like final states, a 2p state only s and d.

## Why the simple approximations fail

Early calculations replaced the matrix elements by the orbital
resolved partial density of states of the conduction band, on the
grounds that the core level is localized and orthogonal to the final
states.  This loses the modulation of peak intensities by the matrix
elements, and it presumes a decomposition of delocalized conduction
states into orbital components that has no precise meaning.

The second difficulty is the core hole.  Removing an electron from
the core leaves a positive charge that the excited electron feels,
much as an exciton does in a semiconductor.  In metals the conduction
electrons screen it; in insulators it reshapes the spectrum.  Two
standard shortcuts exist.  The Z+1 approximation replaces the
excited atom by the next element in the periodic table, which cannot
distinguish holes in different core levels.  The Slater transition
state removes half an electron from the core and places it in the
lowest empty state, which does not describe the hole and electron
interacting.

## The supercell OLCAO method

The method Imago uses was built on OLCAO specifically for this
problem, around 2000.  It rests on the core orthogonalization step
that gives OLCAO its name, and is described on the
[core orthogonalization](core_orthogonalization.html) page: the core
orbitals of the one target atom are retained in the secular equation
while every other atom's core is orthogonalized away.  The target
atom's core states, and transitions out of them, are then explicit in
the solution.

The core hole is treated by two separate self-consistent
calculations rather than by a shortcut.

1. **The ground state.**  The supercell is solved with the target
   atom's core retained.  This gives the initial core state and the
   ground state total energy.
2. **The excited state.**  One electron is removed from the chosen
   core orbital of the target atom and placed in the lowest
   conduction band state, keeping the cell neutral, and the
   Kohn-Sham equation is solved again to self-consistency.  The
   excited electron and the core hole relax together, so their
   interaction and all multiple scattering are included in the final
   states.  This gives the excited state total energy and the final
   state wave functions.

The spectrum is then the golden rule sum with the ground state core
orbital as initial state and the excited state conduction states as
final states.  The transition energy, which fixes where the edge
onset falls on the absolute energy axis, is the difference of the two
total energies, an N-electron system against an N-electron system
with a core hole and a conduction electron.  Being able to place the
onset absolutely matters because the measured spectrum of an element
is the sum over all its inequivalent sites, each with a slightly
different onset, and small shifts between sites change the summed
shape substantially.  A Gaussian broadening of about one electron
volt is applied to represent the finite lifetime of the core hole and
the instrumental resolution.

**Why a supercell.**  The core hole is a localized perturbation, and
periodic boundary conditions replicate it in every image cell.  The
cell must be large enough that the holes do not interact, which
means roughly nine to ten angstroms between images of the target
atom.  For a compact cubic cell that is a modest supercell; for an
elongated cell, a model of an interface, or a defect structure it
can mean several hundred atoms.  Since the Brillouin zone of such a
cell is small, a single k-point suffices for the self-consistent
parts, which keeps the cost in hand.

## What the method offers

The book lists the advantages that have made this approach widely
used, and they follow from the design above.  The ground state and
excited state are computed separately, so the electron-hole
interaction is one step beyond a single-particle description.  The
core states are physically present, not represented by an effective
potential.  The dipole matrix elements come from the ab initio wave
functions, impose the selection rules by themselves, and resolve
into Cartesian components so that the anisotropy of the spectrum
can be studied.  The onset is absolute.  And because no atomic
radius enters anywhere, the method applies to any element, any
edge, and any kind of material, from hydrogen in a biomolecule to a
rare earth in an oxide, with the same procedure.

## In Imago

The calculation is prepared and run in two stages.

- `makeinput.py -xanes` chooses the target atoms and writes one
  complete input set per target.  By default one atom of each species
  is excited; the `-atom` suboption names specific atoms instead.
  In each input set the target atom's core orbitals are included in
  the basis, and every atom within a radius of the target, 3.5
  angstroms by default and adjustable with `-sphere`, is given its
  own potential type so that the screening response to the core hole
  can be represented site by site.
- `imago.py -pacs EDGE` runs the calculation for one input set.  The
  edge, such as `1s` or `2p`, is required: it names the core orbital
  the electron is removed from and selects the excited state
  potential.  The driver runs the ground state and excited state
  self-consistent cycles, evaluates the transitions, places the onset
  at the total energy difference, and writes a single spectrum with
  x, y, z, and averaged columns.  The extended basis is used for the
  final states, since the spectrum reaches forty to fifty electron
  volts above the edge.

Spectra for an element with several inequivalent sites are obtained
by running each target's input set and summing the results with the
site multiplicities.
