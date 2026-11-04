---
layout: default
title: pdf_neutron.py
script: pdf_neutron.py
description: 'Computes the normalised, optionally neutron-weighted total pair distribution function g(r) or reduced G(r) of a periodic model from an imago.skl file, including periodic images, for comparison with diffraction experiments. Used by deadmd.py for its PDF fitness metric.'
parent: Script Catalog
---

# pdf_neutron

Computes the normalised, optionally neutron-weighted pair
distribution function of a periodic model. It reads an Imago skeleton
file (the `cell` line and the `fractional N` atom list, as written by
`dump2skl.py`), counts every atom pair out to `rmax` including
periodic images, and writes two columns: r and the requested
function.

**Relation to rpdf.** This is *not* a replacement for the `rpdf`
geometric pair distribution tool in `structure_control.py`. That tool
counts every pair as 1 and offers bounding-box and element-pair
filters, but its histogram is currently divided by N r^2 only, so its
large-r level follows the model's number density instead of settling
at 1. This script computes the *total* function to compare with a
diffraction experiment: normalised so that g(r) tends to 1 (or G(r)
to 0), and weighted by neutron scattering lengths by default, which
makes H pairs negative as in a measured neutron PDF. With
`-weighting none` it gives the normalised geometric g(r).
[deadmd.py](deadmd.html) uses it for its PDF fitness metric.

## Definitions

Following Keen, *J. Appl. Cryst.* **34**, 172 (2001):

```
partials   g_ij(r) = n_ij(r) / (N_i 4 pi r^2 dr rho_j)
total      g(r)    = 1 + sum_ij c_i c_j b_i b_j (g_ij(r) - 1) / <b^2>
                     with <b^2> = sum_i c_i b_i^2,  so g -> 1 at large r
reduced    G(r)    = 4 pi rho r (g(r) - 1),  oscillating about 0
```

The b_i are the bound coherent neutron scattering lengths in fm from
NIST, or 1 for `-weighting none`. H is negative, so H-X peaks are
negative in a neutron-weighted curve.

The total is normalised by the mean square length, not the textbook
square of the mean. In H-rich models the negative H cancels the rest
and the mean length is near zero (for HMDS its square is about a
hundredth of the mean square), so dividing by it would amplify every
small deviation a hundredfold and the curve would no longer settle at
1. The two differ only by a constant factor on g minus 1, which the
scale fit in deadmd.py absorbs anyway.

## Usage

**Commandline:**

```
pdf_neutron.py [-i SKL] [-o OUT] [-rmax R] [-dr DR] [-sigma S]
               [-weighting {neutron,none}] [-function {G,g}] [-full]
```

**Example, a calculated reference for deadmd:**

```
pdf_neutron.py -i ref.skl -o ref_G.dat -rmax 12 -full
```

**Option Table:**

| Option       | Defaults            | Description |
|:------------:|:-------------------:|:------------|
|-i            |imago.skl            | Input skeleton file. |
|-o            |pdf_neutron.plot     | Output file: two columns, r and the function. |
|-rmax         |10.0                 | Largest r in Angstrom. |
|-dr           |0.01                 | Histogram bin width in Angstrom. |
|-sigma        |0.05                 | Gaussian broadening in Angstrom applied to the histogram, standing in for thermal motion and finite-Q resolution. |
|-weighting    |neutron              | `neutron` for scattering-length weighting, `none` for the plain number-weighted function. |
|-function     |G                    | Write the reduced G(r) or the total g(r). |
|-full         |Off                  | Write out to rmax regardless of the cell size. See below. |

## Notes

Use the **same** sigma, weighting and function for the reference
curve and for the models being compared.

By default the output stops at half the shortest perpendicular cell
width, or at rmax if that is smaller. Beyond that distance the
periodic images repeat the same neighbours and the curve is an
artefact of the box. `-full` writes out to rmax regardless, which is
what a deadmd reference file needs so that it reaches 10 Angstrom.
