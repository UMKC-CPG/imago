---
layout: default
title: org_radwavefn.py
script: org_radwavefn.py
description: 'Parser for Grasp2K Dirac four-component relativistic wavefunctions. Reads rwfn.out and isodata from the current directory and writes a CSV of orbital eigenvalues plus one CSV per orbital of the large and small components on the radial grid. Also the library that interFit.py and make_veusz_graph.py import.'
parent: Script Catalog
---

# org_radwavefn

Reads the Grasp2K output file `rwfn.out`, which holds the numerical
radial wavefunctions of a relativistic atom, the large component P(r)
and the small component Q(r) on a log-linear grid, and organizes
them for Gaussian fitting. It also parses the companion `isodata`
file for the nuclear charge-distribution parameters. The module is
primarily a library: [interFit.py](interFit.html) and
[make_veusz_graph.py](make_veusz_graph.html) import it. Run as a
script it dumps everything it parsed to plain CSV files for
inspection.

**Relativistic bookkeeping.** In the Dirac equation each
non-relativistic orbital with angular momentum l splits into two
states with total angular momentum j = l plus or minus one half.
Grasp2K labels them by the kappa quantum number, which encodes both:
kappa = -(l+1) for j = l + 1/2 and kappa = l for j = l - 1/2. The
small component has angular momentum differing from the large by one
unit. The module carries the conversion functions between kappa, l
and j, and three fitting modes that trade relativistic detail for
simplicity: single-component (j-split pairs averaged into one
orbital per n and l), two-component (each j kept) and
four-component (large and small fitted per j). A Gaussian shrinking
function can be applied before fitting so the orbitals vanish at a
finite cutoff radius.

## Usage

**Commandline:**

```
org_radwavefn.py
```

There are no options. Run it in a directory holding `isodata` and
`rwfn.out`, such as a `<BASIS>/<element>/` directory produced by
[graspElems.py](graspElems.html).

## Output

**`eigenvalues`**, a CSV with one row per orbital:

```
n_quantum_number,k_quantum_number,large_l_quantum_number,small_l_quantum_number,energy_eigenvalue
```

**`<Element>_n_<n>_k_<kappa>.dat`**, one CSV per orbital, with the
radial grid, the large component and the large component divided by
r, and the small component and the small component divided by r, in
columns. Plot these with any CSV-aware tool, or use
make_veusz_graph.py for a ready-made Veusz document.
