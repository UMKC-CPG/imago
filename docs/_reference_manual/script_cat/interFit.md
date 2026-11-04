---
layout: default
title: interFit.py
script: interFit.py
description: 'Gaussian fitting of Grasp2K relativistic atomic orbitals. Reads the rwfn.out and isodata files of a Grasp2K run, sweeps the maximum Gaussian exponent and the number of terms per angular momentum channel, keeps the expansion with the lowest weighted error, and writes the atom-type and potential-type block that Imago reads, in one-, two- or four-component form.'
parent: Script Catalog
---

# interFit

Fits the relativistic atomic orbitals computed by Grasp2K with
Gaussian expansions, producing the basis-function block that Imago
reads for an element. This is the relativistic counterpart of the
non-relativistic basis construction in [makeBDB.py](makeBDB.html):
where makeBDB fits a potential from the atomic SCF program, interFit
fits orbitals from a Dirac-equation atomic code, so that a
relativistic Imago calculation has a basis built on the same
footing.

The script reads the Grasp2K output through `org_radwavefn.py`,
performs a parameter sweep over the maximum Gaussian exponent and the
number of terms in each angular-momentum channel, and keeps the
expansion with the lowest weighted root-mean-square error. The
exponents form a geometric series, and for each trial set the
coefficients come from the minimum-norm least-squares system solved
by a pivotless LU decomposition.

**Three fitting modes** trade relativistic detail for simplicity:

| `-comp` | Mode             | What is fitted |
|:-------:|:-----------------|:---------------|
|1        | single-component | The two j-split partners of each (n, l) are averaged into one orbital. |
|2        | two-component    | Each j component is kept and fitted separately. |
|4        | four-component   | Both the large and the small Dirac component of each j are fitted. |

A Gaussian damping, or shrinking, function can be applied before
fitting so that the orbitals vanish at a finite cutoff radius.

**Defaults.** Every option's default value is read from the resource
control file `interFitrc.py` in the directory named by `$IMAGO_RC`.

## Usage

**Commandline:**

```
interFit.py [-comp {1,2,4}] [-shrink RC SIG] [-tol TOL] [--help]
```

Run it in the directory of a completed Grasp2K calculation holding
`isodata` and `rwfn.out`, such as one produced by
[graspElems.py](graspElems.html).

**Examples:**

```
interFit.py -comp 1              # single-component (averaged j)
interFit.py -comp 4              # four-component (large + small)
interFit.py -comp 1 -shrink 5 1  # shrinking function, rc=5, sigma=1
interFit.py -tol 1e-4            # tighter fitting tolerance
```

**Option Table:**

| Option    | Defaults  | Description |
|:---------:|:---------:|:------------|
|-comp      |1          | Number of components to generate: 1, 2 or 4. |
|-shrink    |0 0 (off)  | Apply the shrinking function with critical radius RC and width SIG. |
|-tol       |1e-3       | RMSE tolerance for accepting a fit. |
|--help     |N/A        | Print usage information and exit. |

## Output

`waveFn.dat`, holding the combined atom-type and potential-type
sections in the format read by Imago's atomic-type and
potential-type readers. The sweep's progress and the best fit found
are printed as it runs.

## Notes

The script's own docstring records a list of performance
improvements that have been identified but not yet made (vectorising
the design-matrix construction, compiling the LU solve, and
parallelising the outer exponent loop), and a planned decoupling of
the s and p term counts, which currently share one value. The sweep
can take a long time for heavy elements; the `-tol` setting is the
main lever on run time.
