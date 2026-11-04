---
layout: default
title: make_veusz_graph.py
script: make_veusz_graph.py
description: 'Writes Veusz plot documents of the Grasp2K relativistic orbitals in the current directory: one page per orbital showing the large and small components divided by r against radius, plus a second document of the fitted orbital equations when a fit is present.'
parent: Script Catalog
---

# make_veusz_graph

Creates Veusz plotting documents for the relativistic atomic orbitals
in a Grasp2K output directory, so that the numerical orbitals, and
the Gaussian fits made to them by [interFit.py](interFit.html), can
be inspected visually. Veusz is a scientific plotting package whose
documents are plain text scripts; this tool writes such a script,
with the data embedded, that Veusz opens directly.

The data are read with [org_radwavefn.py](org_radwavefn.html), and
two documents are attempted.

* **Numerical orbitals.** One page per orbital, named from the
  element, n and kappa, plotting the large component P(r)/r and the
  small component Q(r)/r against radius in Angstroms, with a legend
  giving the orbital as n, l and j (for example `2p(3/2)`).
* **Fitted equations.** The same pages drawn from the Gaussian
  expansions, when an equations file from a fit is present. If none
  is, the script says so and writes only the first document.

## Usage

**Commandline:**

```
make_veusz_graph.py
```

There are no options. Run it in a directory holding `isodata` and
`rwfn.out`, such as a `<BASIS>/<element>/` directory produced by
[graspElems.py](graspElems.html).

## Output

`<Element>_grasp_numerical_orbitals.vsz` in the current directory,
and the fitted-equation document when available. Open them with
`veusz <file>.vsz`.
