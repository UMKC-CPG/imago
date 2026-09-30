#!/usr/bin/env python3
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Mohammed Belhadj Larbi

"""
pdf_neutron.py -- normalised, optionally neutron-weighted pair distribution
function of a periodic model.

USAGE:
    pdf_neutron.py [-i imago.skl] [-o pdf_neutron.plot] [-rmax 10.0]
                  [-dr 0.01] [-sigma 0.05] [-weighting neutron|none]
                  [-function G|g] [-full]

Reads an imago skeleton file (the `cell` a b c alpha beta gamma line and a
`fractional N` atom list, as written by dump2skl.py), counts every atom pair
out to rmax including periodic images, and writes two columns: r and the
requested function.

Relation to rpdf (structure_control.py compute_rpdf): NOT a replacement.
rpdf is the geometric pair distribution tool (every pair counts 1) with
bounding-box and element-pair filters; its histogram is currently divided
by N r^2 only, not by 4 pi rho r^2 dr, so its large-r level follows the
model's number density instead of 1 (a fix is planned).  pdf_neutron.py
computes the TOTAL function to compare with a diffraction experiment:
normalised so g(r) -> 1 (or G(r) -> 0), and weighted by neutron scattering
lengths by default, which puts H pairs negative as in measured neutron
PDFs.  With "-weighting none" it gives the normalised geometric g(r).
deadmd.py uses it for its PDF fitness metric.  Once rpdf is normalised,
the neutron weighting could move into rpdf and the two tools merge.

Definitions (Keen, J. Appl. Cryst. 34, 172 (2001)):
    partials   g_ij(r) = n_ij(r) / (N_i 4 pi r^2 dr rho_j)
    total      g(r)    = 1 + sum_ij c_i c_j b_i b_j (g_ij(r) - 1) / <b^2>
               with <b^2> = sum_i c_i b_i^2  -> 1 at large r.
               b_i = 1 for "-weighting none" (plain number-weighted g(r)).
    reduced    G(r)    = 4 pi rho r (g(r) - 1)  -> oscillates about 0.
b_i are the bound coherent neutron scattering lengths (fm, NIST).  H is
negative, so H-X peaks are negative in a neutron-weighted g(r) and G(r).
Normalised by <b^2>, not the textbook <b>^2: in H-rich models the negative
H cancels the rest and <b> is near 0 (HMDS: <b>^2 is ~1/100 of <b^2>), so
dividing by <b>^2 amplifies every small deviation ~100x and the curve no
longer settles at 1.  The two differ only by a constant factor on g - 1,
which the scale fit in deadmd.py absorbs anyway.

-sigma is a Gaussian broadening (Angstrom) applied to the histogram, standing
in for thermal motion and finite-Q resolution.  Use the SAME sigma, weighting
and function for the reference curve and for the models being compared.

By default the output stops at half the shortest perpendicular cell width
(or rmax if smaller): beyond it the periodic images repeat the same
neighbours and the curve is an artefact of the box.  -full writes out to
rmax regardless.
"""

import argparse
import math
import sys

import numpy as np

# Bound coherent neutron scattering lengths, fm (NIST, natural abundance).
SCATTERING_LENGTH = {
    "h": -3.7390, "d": 6.671, "he": 3.26, "li": -1.90, "be": 7.79,
    "b": 5.30, "c": 6.6460, "n": 9.36, "o": 5.803, "f": 5.654,
    "na": 3.63, "mg": 5.375, "al": 3.449, "si": 4.1491, "p": 5.13,
    "s": 2.847, "cl": 9.5770, "k": 3.67, "ca": 4.70, "ti": -3.438,
    "fe": 9.45, "ga": 7.288, "ge": 8.185, "as": 6.58, "zr": 7.16,
}


def read_skl(path):
    """Return (lattice 3x3 rows a,b,c in Angstrom, fractional Nx3, elements)."""
    with open(path) as f:
        lines = [ln.split() for ln in f]
    for k, words in enumerate(lines):
        if words and words[0].lower() == "cell":
            a, b, c, al, be, ga = map(float, lines[k + 1][:6])
        if words and words[0].lower() == "fractional":
            n = int(words[1])
            atoms = lines[k + 1:k + 1 + n]
            break
    else:
        sys.exit(f"pdf_neutron.py: no 'fractional' atom list in {path}")
    al, be, ga = (math.radians(x) for x in (al, be, ga))
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    lattice = np.array([[a, 0.0, 0.0],
                        [b * math.cos(ga), b * math.sin(ga), 0.0],
                        [cx, cy, math.sqrt(c * c - cx * cx - cy * cy)]])
    elements = [w[0].lower() for w in atoms]
    frac = np.array([[float(x) for x in w[1:4]] for w in atoms])
    return lattice, frac, elements


def pair_histograms(lattice, frac, types, ntypes, rmax, dr):
    """Ordered-pair counts n_ij per bin, including periodic images."""
    volume = abs(np.linalg.det(lattice))
    # Perpendicular widths decide how many image cells rmax reaches.
    widths = [volume / np.linalg.norm(np.cross(lattice[(k + 1) % 3],
                                               lattice[(k + 2) % 3]))
              for k in range(3)]
    reach = [int(math.ceil(rmax / w)) for w in widths]
    shifts = np.array([[i, j, k]
                       for i in range(-reach[0], reach[0] + 1)
                       for j in range(-reach[1], reach[1] + 1)
                       for k in range(-reach[2], reach[2] + 1)], float)
    cart = (frac % 1.0) @ lattice
    images = (shifts @ lattice)                       # (S, 3)
    # Bin k is centred on r = k * dr, so the output grid is 0.01, 0.02, ...
    # (the grid deadmd.py aligns reference curves on).
    nbins = int(round(rmax / dr)) + 1
    hist = np.zeros((ntypes, ntypes, nbins))
    for i in range(len(cart)):
        # Vectors from atom i to every atom in every image cell.
        d = (cart[None, :, :] + images[:, None, :] - cart[i]).reshape(-1, 3)
        r = np.sqrt(np.einsum("ij,ij->i", d, d))
        tj = np.tile(types, len(images))
        bins = np.rint(r / dr).astype(int)
        keep = (r > 1e-8) & (bins < nbins)
        bins = bins[keep]
        np.add.at(hist[types[i]], (tj[keep], bins), 1.0)
    return hist, volume, widths


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("-i", default="imago.skl")
    ap.add_argument("-o", default="pdf_neutron.plot")
    ap.add_argument("-rmax", type=float, default=10.0)
    ap.add_argument("-dr", type=float, default=0.01)
    ap.add_argument("-sigma", type=float, default=0.05)
    ap.add_argument("-weighting", choices=("neutron", "none"),
                    default="neutron")
    ap.add_argument("-function", choices=("G", "g"), default="G")
    ap.add_argument("-full", action="store_true")
    args = ap.parse_args()

    lattice, frac, elements = read_skl(args.i)
    species = sorted(set(elements))
    types = np.array([species.index(e) for e in elements])
    counts = np.array([np.sum(types == t) for t in range(len(species))], float)
    natoms = len(elements)
    conc = counts / natoms
    if args.weighting == "neutron":
        missing = [e for e in species if e not in SCATTERING_LENGTH]
        if missing:
            sys.exit(f"pdf_neutron.py: no scattering length for {missing}; "
                     "add it to SCATTERING_LENGTH")
        b = np.array([SCATTERING_LENGTH[e] for e in species])
    else:
        b = np.ones(len(species))

    hist, volume, widths = pair_histograms(lattice, frac, types,
                                           len(species), args.rmax, args.dr)
    rho = natoms / volume
    nbins = hist.shape[2]
    r = np.arange(nbins) * args.dr

    # Area-preserving Gaussian broadening of every partial histogram.
    if args.sigma > 0:
        half = int(math.ceil(5 * args.sigma / args.dr))
        x = np.arange(-half, half + 1) * args.dr
        kernel = np.exp(-0.5 * (x / args.sigma) ** 2)
        kernel /= kernel.sum()
        hist = np.apply_along_axis(
            lambda h: np.convolve(h, kernel, mode="same"), 2, hist)

    shell = 4.0 * math.pi * r ** 2 * args.dr
    shell[0] = np.inf                     # r = 0: no pairs, avoid 0/0
    g_minus_1 = np.zeros(nbins)
    for i in range(len(species)):
        for j in range(len(species)):
            g_ij = hist[i, j] / (counts[i] * shell * rho * conc[j])
            g_minus_1 += conc[i] * conc[j] * b[i] * b[j] * (g_ij - 1.0)
    g_total = 1.0 + g_minus_1 / np.dot(conc, b * b)

    out = (g_total if args.function == "g"
           else 4.0 * math.pi * rho * r * (g_total - 1.0))
    rlimit = args.rmax if args.full else min(args.rmax, 0.5 * min(widths))
    keep = (r > 0) & (r <= rlimit + 1e-9)
    with open(args.o, "w") as f:
        for ri, vi in zip(r[keep], out[keep]):
            f.write(f"{ri:.4f} {vi:.6f}\n")
    print(f"pdf_neutron.py: {natoms} atoms {species}, rho {rho:.5f} /A^3, "
          f"weighting {args.weighting}, {args.function}(r), sigma "
          f"{args.sigma}, written to r = {r[keep][-1]:.2f} A -> {args.o}")


if __name__ == "__main__":
    main()
