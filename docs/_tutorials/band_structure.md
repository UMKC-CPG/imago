---
layout: default
title: A band structure
nav_order: 4
---

# A band structure: silicon along the high-symmetry path

This tutorial draws the electronic band structure of the silicon
crystal from [Your first calculation](first_calculation.html),
reusing that job's converged SCF. It takes two commands. Most of the
page is about reading the result.

## 1. What a band structure calculation is

The density of states counted the electronic states at each energy
summed over the whole Brillouin zone. A band structure instead
follows each state's energy as the k-vector moves along chosen
straight lines through the zone, the lines that connect its
high-symmetry points. Because the points are chosen by symmetry, the
path is fixed by the lattice type, and the plot answers questions the
DOS cannot: where in the zone the valence band maximum and
conduction band minimum sit, whether the gap is direct, and how
heavy the bands are.

The path was already recorded when the job was prepared, by the
`-sybdpath fcc` option to `makeinput.py`. For a face-centred-cubic
lattice it is the standard tour

```
GAMMA - X - W - K - GAMMA - L - U - W - L - K | U - X
```

in two segments: the first runs through the zone and the second,
after a break, covers the remaining edge from U to X. The definition
lives in `$IMAGO_DATA/sybdDB/fcc`; the other lattice types have
their own files there. If you prepared the job without `-sybdpath`,
re-run `makeinput.py` with it now; the converged potential is not
affected.

## 2. Run it

In the silicon job directory:

```
imago.py -sybd
```

The driver finds `gs_scfV-fb.dat`, skips the SCF pass, and runs the
band-structure pass against it. Band structures are evaluated at
general k-points, so this pass uses the general `imago` engine even
when the SCF ran on `imagoG`; the driver chooses. When the engine
finishes, the driver itself calls `makeSYBD.py` to merge the raw
eigenvalues with the k-point labels into the plot file, so the
directory now holds:

| File                 | Contents |
|:---------------------|:---------|
|`gs_pscf+sybd-fb.out` | Full log of the band-structure pass. |
|`gs_sybd-fb.plot`     | The plottable bands: one line per k-point along the path, the distance along the path first, then every band energy in eV relative to the Fermi level, with the high-symmetry points' positions and labels appended to the first lines. |
|`gs_sybd-fb.raw.31`   | The raw eigenvalues the engine wrote. |
|`gs_sybd-fb.raw.33`   | The valence-dimension metadata: which atom and orbital each basis function belongs to, for later decomposition. |

## 3. Plot it

```
plotgraph.py -i gs_sybd-fb.plot -ft sybd
```

The `sybd` figure type reads the plot file's trailing columns to
draw a vertical line at each high-symmetry point and to label the
horizontal axis with the point names, with `G` standing for GAMMA.
Only the matplotlib backend supports this figure type.

What to look for in silicon:

* **The valence band maximum at GAMMA**, at 0 eV, where three bands
  meet.
* **The conduction band minimum along GAMMA to X**, at roughly 80
  percent of the way to X, not at a symmetry point. This is why
  silicon's gap is **indirect**: the lowest empty state and the
  highest filled state sit at different k. The `Gap_Kind` column of
  `gs_iter-fb.dat` reported the same thing from the mesh.
* **The gap** between those two, around 0.5 to 1 eV in the
  local-density approximation. The experimental value is 1.1 eV;
  the underestimate is a property of the approximation, not an
  error in the calculation.
* **Band degeneracies** along high-symmetry lines that split away
  from them, which is symmetry made visible.

Use `-yi -13 -ya 6` to frame the valence and lower conduction bands
and hide the deep bands.

## 4. Vary it

* **Path resolution.** The number of k-points along the path is set
  in the sybdDB file (300 for fcc). Edit a copy of the file in your
  own directory and point `-sybdpath` at it by path to change the
  sampling or the route.
* **A different lattice.** For a non-cubic crystal, choose the
  sybdDB file that matches the lattice type: `hex` for hexagonal,
  `tet` for simple tetragonal, and so on. `makeinput.py` warns if
  the file's lattice type does not match the cell.
* **Spin-polarized systems** produce `gs_sybd-fb.up.plot` and
  `gs_sybd-fb.dn.plot` instead of one file; plot each.

The band structure is documented on the
[imago.py](../reference_manual/script_cat/imago.html) page under
Symmetric Band Structure, and the merge step on the
[makeSYBD.py](../reference_manual/script_cat/makeSYBD.html) page.
