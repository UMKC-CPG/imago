---
layout: default
title: Your first calculation
nav_order: 3
---

# Your first calculation: the density of states of silicon

This tutorial takes crystalline silicon from the skeleton file
shipped with Imago to a plot of its density of states. Along the way
it explains every file that appears, because those same files appear
in every Imago calculation you will ever run. Silicon is chosen
because its two-atom primitive cell runs in seconds, and because its
result, a semiconductor with a gap of about a volt, is one you can
check by eye.

It assumes Imago is [installed](installation.html) and your shell
has sourced `imagorc`.

## 1. Make a job directory

Every calculation lives in its own directory, and every file the
calculation produces lands there. Keep them under one `jobs/`
directory so that they are easy to find and easy to clean up:

```
mkdir -p $IMAGO_DIR/jobs/si
cd $IMAGO_DIR/jobs/si
cp $IMAGO_DIR/skl/silicon imago.skl
```

The copy must be named `imago.skl`; that is the name `makeinput.py`
looks for.

## 2. Read the skeleton file

The skeleton is the only input you write by hand. Here is silicon:

```
title
Si crystal
end
cell
   5.43070 5.43070 5.43070 90.000 90.000 90.000
fract 1
Si 0.0000000 0.00000000 0.00000000
space 227_a
supercell 1 1 1
prim
```

Line by line:

* **`title` ... `end`** is free text carried into the output files.
* **`cell`** gives the lattice parameters a, b, c in Angstroms and
  the angles alpha, beta, gamma in degrees. Silicon is cubic with a
  5.43 Angstrom edge.
* **`fract 1`** says that one atom follows, in fractional
  coordinates. This is not the whole cell. It is the **asymmetric
  unit**: the smallest set of atoms from which the space group
  regenerates every other.
* **`space 227_a`** names the space group, Fd-3m, number 227 in the
  International Tables, in its first origin setting `a`. Imago
  applies its symmetry operations to the one silicon atom and
  produces the eight atoms of the conventional cell. This is why
  a crystal is entered as one atom: the symmetry is the data, and
  Imago uses it again later to reduce the k-point mesh.
* **`supercell 1 1 1`** asks for no replication.
* **`prim`** asks for the primitive cell, two atoms, rather than the
  eight-atom conventional cell (`full`). The primitive cell is four
  times cheaper and gives the same physics.

## 3. Generate the inputs

`makeinput.py` turns the skeleton into everything the engine needs:

```
makeinput.py -kpd 200 -sybdpath fcc -slurm -p PARTITION -a ACCOUNT
```

Three choices are being made here.

* **`-kpd 200`** sets the k-point density: the number of k-points
  per unit volume of reciprocal space, in Bohr^-3. Imago turns it
  into a mesh as uniform as the cell allows and then reduces it by
  symmetry. The default without this option is the Gamma point
  alone, which is right for a large amorphous cell but far too
  coarse for a two-atom crystal. 200 is a comfortable value for a
  first look; converge it later by raising it until the total
  energy stops moving.
* **`-sybdpath fcc`** records the high-symmetry k-point path through
  the Brillouin zone of a face-centred-cubic lattice. It costs
  nothing now and lets the [next tutorial](band_structure.html)
  reuse this job for a band structure.
* **`-slurm -p PARTITION -a ACCOUNT`** names the queue and the
  account for the batch script. The defaults in `makeinputrc.py`
  are the authors' own and will not exist on your cluster, so give
  yours. Add `-t HH:MM:SS` for a wall time other than one hour and
  `-m 4G` for memory other than ten gigabytes.

The directory now holds:

| File or directory | Contents |
|:------------------|:---------|
|`inputs/`          | The engine's input set: `imago.dat` (every setting the engine reads), `structure.dat` (the expanded cell, every atom explicit), `scfV.dat` (the starting potential, assembled from the atomic potential database), `kp-scf.dat` and `kp-pscf.dat` (the k-point meshes), and the expanded skeletons `imago.fract-mi` and `imago.cart-mi` with the atom-number map `datSkl.map`. |
|`summary`          | A one-screen digest: the lattice vectors in atomic units, the basis dimensions, the number of atoms, electrons and states, the k-point settings, and the convergence limit. Read it every time; it is where a wrong skeleton shows first. |
|`slurm`            | A ready-to-submit batch script. |
|`command`          | A dated log of every Imago command run in this directory, appended to by every script. |

Open `summary`. For silicon you should see two atoms, 28 electrons,
and `SCF KPoints = density=200.0`.

## 4. Run it

Open `slurm`. Its last line runs the engine:

```
"$IMAGO_BIN/imago.py" -bond
```

The generated script requests a bond-order calculation by default.
Change that line to request the density of states instead:

```
"$IMAGO_BIN/imago.py" -dos
```

then submit it:

```
sbatch slurm
```

On a workstation, or in an interactive session on a compute node,
run the same command directly:

```
imago.py -dos
```

`-dos` is a **two-pass job**. The first pass is the self-consistent
field (SCF) calculation: starting from the atomic potentials, the
engine solves the eigenvalue problem, builds a charge density,
derives a new potential, and repeats until the potential stops
changing. The second pass takes that converged potential and
computes the density of states on the post-SCF k-point mesh. Each
pass can use a different basis; both default to the full basis,
`fb`, here. The single-pass variant `-scfdos` computes the DOS from
the SCF pass itself and is a little faster, but the two-pass form is
the one to learn, because its SCF is **reusable**: once
`gs_scfV-fb.dat` exists, any later `-dos`, `-bond` or `-sybd` job in
this directory skips straight to its own second pass.

Silicon at this density takes well under a minute.

## 5. Watch it converge

Two files tell you what happened.

**`gs_iter-fb.dat`** has one line per SCF iteration:

```
 Iter. Occ._Energy   El._Error   Convergence   Total_Energy(au)      Mag._Mom.          Gap(au)  Gap_Kind
    1   ...
```

`Convergence` is the change in the potential between iterations and
must fall below the convergence limit in `summary` (0.0001 by
default). `Total_Energy(au)` is in hartree and should settle to a
steady value. `Gap(au)` is the gap between the highest occupied and
lowest unoccupied state on the mesh, and `Gap_Kind` whether it is
direct or indirect. A metal shows a gap of zero.

**`runtime`** is the driver's log. Each pass writes a header, its
timings, and the line `Program Sequence Complete.` when it finished.
For a `-dos` job you should see two such blocks, `gs_scf-fb` and
`gs_pscf+dos-fb`. If the second is missing, the SCF did not
complete; its own log is `gs_scf-fb.out`.

## 6. Know the output files

The naming convention is `<edge>_<job>-<basis>.<kind>`. `gs` is the
ground state, as opposed to a core-hole edge such as `1s`; the basis
is the one that pass used. For this job:

| File                   | Contents |
|:-----------------------|:---------|
|`gs_scf-fb.out`         | Full log of the SCF pass. |
|`gs_scfV-fb.dat`        | The converged potential: the reusable result of the SCF. |
|`gs_iter-fb.dat`        | The iteration history above. |
|`gs_enrg-fb.dat`        | The converged total energy broken into kinetic, electrostatic and exchange-correlation parts. |
|`gs_pscf+dos-fb.out`    | Full log of the DOS pass. |
|`gs_dos-fb.t.plot`      | The **total** DOS: two columns, energy in eV relative to the Fermi level and the DOS. Directly plottable. |
|`gs_dos-fb.p.raw`       | The **partial** DOS, one data unit per atom type with one column per orbital. Not plottable until combined by `makePDOS.py`. |
|`gs_dos-fb.loci.plot`   | The localization index of each state, a measure of how spatially confined it is. |
|`intermediate`          | A symbolic link into `$IMAGO_TEMP` where the engine's large working files, the HDF5 matrices above all, were kept. The link is what `tidy_scratch.py` follows to reclaim that space later. |

## 7. Plot the total DOS

```
plotgraph.py -i gs_dos-fb.t.plot -ft dos
```

A matplotlib window opens. You should see the valence band ending at
0 eV, a gap of roughly 0.5 to 1 eV (the local-density approximation
underestimates silicon's 1.1 eV), and the conduction band beyond it.
Hover over the curve for its values. The window was drawn by
`mplscript.py`, which plotgraph wrote into the directory and ran;
edit and re-run that file for any change the command line does not
offer, or add `-xi -15 -xa 10` to set the energy range directly.

## 8. Plot the partial DOS

The raw partial file holds, for the one silicon type, a column for
the total and a column for each valence orbital. To see the s and p
contributions separately, write a control file `pdos.ctl`:

```
LABEL : ELEMENT_NAME : SPECIES_ID : TYPE_ID : COL_LABELS
Total : all : all : all : TOTAL
Si_s  : si  : all : all : s
Si_p  : si  : all : all : p
```

Each line after the first defines one curve: a label, which atoms
contribute, and which orbital columns to sum. Then:

```
makePDOS.py -i pdos.ctl
plotgraph.py -i PDOS.plot -ft dos -cs 3 -pl
```

`makePDOS.py` writes `PDOS.plot` with the energy and the three
curves; `-cs 3` puts all three in one panel and `-pl` adds a legend.
The s states dominate the bottom of the valence band and the p
states its top, which is the textbook picture of the sp3 solid.

## Where to go from here

* The same directory can now answer other questions without
  repeating the SCF: `imago.py -bond` for bond orders and effective
  charges, then [makeBOND.py](../reference_manual/script_cat/makeBOND.html);
  `imago.py -optc` for optical properties.
* The [next tutorial](band_structure.html) draws the band structure.
* To converge the result, repeat with `-kpd 400` in a fresh
  directory and compare the total energies in `gs_enrg-fb.dat`.
* Every option used here is documented in full on the
  [makeinput.py](../reference_manual/script_cat/makeinput.html),
  [imago.py](../reference_manual/script_cat/imago.html),
  [makePDOS.py](../reference_manual/script_cat/makePDOS.html) and
  [plotgraph.py](../reference_manual/script_cat/plotgraph.html)
  pages.
