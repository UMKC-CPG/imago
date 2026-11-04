---
layout: default
title: Core orthogonalization
parent: Theory
nav_order: 5
---

# Core orthogonalization

*Source: Ching and Rulis (2012), section 3.5.*

## The problem it solves

An all-electron basis carries every orbital of every atom, including
the deep core orbitals.  For a large cell, or for heavy elements,
those core orbitals make up a substantial fraction of the basis, and
the secular equation grows with them.  Diagonalization costs roughly
the cube of the matrix dimension, so the core states, which have no
chemical interest in nearly all calculations, can dominate the cost.

Pseudopotential methods solve this by removing the core electrons
altogether and replacing their effect on the valence electrons with a
fitted potential.  OLCAO keeps the core electrons in the density and
the potential, where they belong physically, and removes the core
orbitals only from the eigenvalue problem.  This is the step the
method is named for, and it was borrowed in 1975 from the
orthogonalized plane wave method, which had used the same trick for
decades.

## The idea

Partition the basis into core and valence Bloch sums and arrange the
overlap and Hamiltonian matrices in blocks:

```
            core        valence
core      [ cc    |    cv    ]
valence   [ vc    |    vv    ]
```

Two approximations are made.  First, the core orbitals of one atom are
taken to have negligible overlap with the core orbitals of any other
atom, so the core-core block is diagonal in the atom index.  This is
very accurate for genuine core states, which are compact.  Second,
each valence Bloch sum is replaced by an orthogonalized version:
itself minus its projection onto every core Bloch sum.

```
b'_valence = b_valence  -  sum over core states c of  <b_c | b_valence> b_c
```

The orthogonalized valence functions are, by construction, orthogonal
to the core.  In this new basis the core-valence blocks vanish, the
core block decouples, and the eigenvalue problem for the valence
states closes on the valence-valence block alone:

```
| H'_vv(k) - E(k) S'_vv(k) | = 0
```

The orthogonalized matrices H' and S' are the original
valence-valence blocks plus correction terms built from the
core-valence overlaps and Hamiltonian elements.  The corrections are
cheap: they involve only matrix products of blocks already computed.
The dimension of the problem drops from the total number of basis
functions to the number of valence basis functions.

## What it costs in accuracy

The only approximation is the neglect of core-core overlap between
different atoms.  The original test, on silicon with the 1s, 2s, and
2p orbitals in the core, found that the eigenvalues near the band gap
changed by less than 0.0001 eV compared with the full problem.  This
is negligible for electronic structure work.

The accuracy can be tuned by moving the boundary between core and
valence.  Any orbital that is shallow enough to overlap its
neighbors, or that participates in bonding, should be in the valence
set, and the basis database defines the boundary per element with
that in mind; see the [basis functions](basis_functions.html) page.
A semicore state that is treated as core when it should not be will
show up as a small but systematic error in states near the gap.

## When the core is kept

There are two reasons to keep core orbitals in the eigenvalue
problem.

The first is to see the core eigenvalues themselves.  Imago's
`-nocore` option of `makeinput.py` moves every orbital the database
designates as core into the valence set, so that no orthogonalization
happens and the core states appear in the band structure and density
of states.  The calculation is correspondingly larger.

The second is core level spectroscopy.  A photo-absorption cross
section or an electron energy loss near-edge spectrum involves a
transition from a specific core orbital of a specific atom to the
unoccupied states.  For that, the core orbital of the one excited atom
is kept in the problem, together with a core hole in its occupation,
while every other atom's core is orthogonalized away as usual.  The
`-xanes` option of `makeinput.py` chooses the excited atoms and
prepares one input set per atom; the `-pacs` option of `imago.py`
then runs each one.  The method is described in chapter 11 of the
book.

## Where it fits in the program

Orthogonalization happens after the overlap and Hamiltonian integrals
are assembled at a k-point and before the secular equation is
diagonalized.  Everything upstream of it, the basis, the potential,
and the integrals, is unaware of the core-valence distinction.
Everything downstream, the eigenvalues, the eigenvectors, and every
property built from them, lives in the reduced valence space.  One
consequence is worth remembering when reading Imago's output: the
dimension reported for the secular equation is the valence dimension,
and the number of states a calculation produces is a multiple of the
number of valence electrons, set by the `-statefactor` option of
`makeinput.py`, not of the total electron count.
