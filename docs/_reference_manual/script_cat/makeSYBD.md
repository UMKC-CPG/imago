---
layout: default
title: makeSYBD.py
script: makeSYBD.py
description: 'Merges the three files of an Imago symmetric band structure (SYBD) calculation into one plot-ready file: each k-point line carries its path coordinate and band energies, and the high-symmetry k-points are annotated with their position and label. imago.py runs it automatically at the end of a -sybd job.'
parent: Script Catalog
---

# makeSYBD

Writes the data file that plotting software uses to draw a symmetric
band structure (SYBD) diagram. An Imago SYBD calculation produces its
band energies, its k-point path, and its high-symmetry k-point names
in three separate files. This script merges them into a single
`.plot` file in which every k-point line holds the path coordinate
followed by the band energies, and the lines belonging to
high-symmetry k-points end with the point's position and label. That
file is what [plotgraph.py](plotgraph.html) reads with `-ft sybd`,
and what Origin, Veusz or gnuplot can read directly.

**You normally do not run this by hand.** When `imago.py` runs a
`-sybd` or `-scfsybd` job it calls this script itself, immediately
after the Fortran executable finishes, and renames the result to
`gs_sybd-fb.plot` (and `gs_sybd-fb.up.plot`, `gs_sybd-fb.dn.plot` for
a spin-polarized run). Run it yourself only to regenerate the plot
file, for example after editing the k-point labels.

**Defaults.** Every option's default value is read from the resource
control file `makeSYBDrc.py` in the directory named by `$IMAGO_RC`.
The file-name defaults are built from the `basis` and `edge` tags
there. Note that those defaults predate the current `imago.py`
naming, so a manual run must name the files explicitly (see below).

## Usage

**Commandline:**

```
makeSYBD.py [-dat DATFILE] [-out OUTFILE] [-raw RAWFILE]
            [-plot PLOTFILE] [-basis BASIS] [-edge EDGE] [--help]
```

**To regenerate the plot file in a finished job directory:**

```
makeSYBD.py -dat inputs/imago.dat -out gs_sybd-fb.out \
            -raw gs_sybd-fb.raw.31 -plot gs_sybd-fb.plot
```

For the spin-down bands of a spin-polarized run use
`-raw gs_sybd-fb.raw.32 -plot gs_sybd-fb.dn.plot`.

**Variable Table:**

| Variable   | Accepted Values     | Description |
|:----------:|:-------------------:|:------------|
|DATFILE     |String; file path    | The Imago input file containing the `SYBD_INPUT_DATA` section, which names the high-symmetry k-points. In a prepared job this is `inputs/imago.dat`. |
|OUTFILE     |String; file path    | The main output of the SYBD run, containing the `HIGH SYMMETRY K POINT` index lines and the k-point table. |
|RAWFILE     |String; file path    | The raw eigenvalue file written by the Fortran program. |
|PLOTFILE    |String; file path    | The merged output file to write. |
|BASIS       |fb, mb, eb           | Basis tag used to build default file names. |
|EDGE        |gs, or an excitation edge | Edge tag used to build default file names. |

**Option Table:**

| Option   | Outputs      | Defaults                  | Description |
|:--------:|:------------:|:-------------------------:|:------------|
|-dat      |N/A           |setup-\<basis\>.dat        | Input file holding the `SYBD_INPUT_DATA` section. |
|-out      |N/A           |\<edge\>_sybd-\<basis\>.out| SYBD output file with the k-point index numbers and positions. |
|-raw      |N/A           |\<edge\>_sybd-\<basis\>.raw| SYBD raw file with the eigenvalues. |
|-plot     |Plot file     |\<edge\>_sybd-\<basis\>.plot| Output plot file. |
|-basis    |N/A           |fb                         | Basis tag for default names. |
|-edge     |N/A           |gs                         | Edge tag for default names. |
|--help    |N/A           |N/A                        | Print usage information and exit. |

## What the Script Does

1. **Reads the k-point names** from the `SYBD_INPUT_DATA` section of
   the input file. That section gives the number of discontinuous
   paths, the number of vertices on each path, and then one line per
   high-symmetry k-point whose last token is its name. The name
   `GAMMA` is shortened to `G` for display.

2. **Reads the k-point positions** from the SYBD output file. It
   first collects the index numbers on the lines containing
   `HIGH SYMMETRY K POINT`, then reads the k-point table that follows
   to find the distance along the path of each of those indices.

3. **Detects path breaks.** When two consecutive high-symmetry
   k-points share the same position, the band path has jumped from
   one segment to another. The two names are merged into one
   comma-separated label, such as `X,U`, and the duplicate position
   is dropped, so the break is drawn as one labelled vertical line.

4. **Reads the raw eigenvalue file and writes the plot file.** The
   raw file begins with a header line holding the number of k-points
   on the path and the number of states. Each k-point then occupies
   one line with its distance along the path followed by its
   energies, written 11 per line, in eV relative to the Fermi
   level. The script concatenates each k-point's lines into one
   output line and, for the high-symmetry k-points in order, appends
   the position value and the label.

## Output Format

The plot file has no header. Each line is one k-point:

```
<path coordinate> <E1> <E2> ... <EN> [<hskp position> <hskp label>]
```

The first column is the distance along the k-path, the next N
columns are the N band energies, and the final two columns, present
only on the first N<sub>hskp</sub> lines, are the position and label
of each high-symmetry k-point in turn. `plotgraph.py -ft sybd`
reads those two trailing columns to draw the vertical lines and
relabel the x axis.
