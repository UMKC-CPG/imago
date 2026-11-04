---
layout: default
title: viewCell.py
script: viewCell.py
description: 'Interactive 3D viewer for the real-space cell, the reciprocal cell, a chosen Brillouin zone, the k-point mesh, the folded k-points scaled by weight, and the band-structure path, read from a BZ.<n> file produced by makeinput.py -printbz. Built on vedo.'
parent: Script Catalog
---

# viewCell

Opens an interactive 3D view of the cell and the k-points used by an
Imago calculation. Six objects are drawn together and each has a
slider to control its transparency, so you can fade the real-space
cell out to study the Brillouin zone, or fade the full mesh out to
see only the folded k-points:

| Object               | Drawn as |
|:---------------------|:---------|
|Real-space cell       | The a, b, c lattice vectors as a box. |
|Reciprocal cell       | The reciprocal lattice vectors as a box. |
|Brillouin zone        | The polyhedron of the zone selected by `-bz`. |
|Mesh k-points         | Red spheres, one per point of the full k-point mesh. |
|Folded k-points       | White glyphs at the symmetry-reduced k-points, sized by weight. |
|Path k-points         | The high-symmetry band-structure path. |

**Input.** The script reads a file named `BZ.<n>` from the current
directory, where `<n>` is the zone number given by `-bz`. That file
is written when you prepare the job with `makeinput.py -printbz <n>`
(see the [makeinput.py](makeinput.html) page), so run this script in
a job directory prepared that way. The file is read as data: it is
parsed into literals and never executed.

**Requirements.** The `vedo` Python package and a working display
(for example an SSH session with X forwarding).

**Defaults.** Every option's default value is read from the resource
control file `viewCellrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
viewCell.py [-bz N] [-rla ALPHA] [-rca ALPHA] [-bza ALPHA]
            [-mkpa ALPHA] [-fkpa ALPHA] [-pkpa ALPHA] [--help]
```

**Variable Table:**

| Variable   | Accepted Values     | Description |
|:----------:|:-------------------:|:------------|
|N           |integer              | Which Brillouin zone to show; selects the `BZ.N` file. |
|ALPHA       |float; 0 to 1        | Initial opacity of an object: 0 is invisible, 1 is solid. Each can be changed with its slider once the viewer is open. |

**Option Table:**

| Option                  | Defaults | Description |
|:-----------------------:|:--------:|:------------|
|-bz / --brillouinzone    |1         | Which Brillouin zone to show. |
|-rla / --real_alpha      |0.1       | Opacity of the real-space cell. |
|-rca / --recip_alpha     |0.4       | Opacity of the reciprocal cell. |
|-bza / --bz_alpha        |0.4       | Opacity of the Brillouin zone. |
|-mkpa / --mesh_kp_alpha  |0.4       | Opacity of the mesh k-points. |
|-fkpa / --fold_kp_alpha  |0.4       | Opacity of the folded k-points. |
|-pkpa / --path_kp_alpha  |0.4       | Opacity of the path k-points. |
|--help                   |N/A       | Print usage information and exit. |

## Notes on the Display

The mesh k-point spheres are sized from the mesh spacing. The folded
k-points are drawn as glyphs whose size is proportional to each
point's weight, scaled so that the most heavily weighted point
matches the mesh sphere size; the lightest points can therefore be
small, which is a faithful picture of how little they contribute.
The view uses an orthographic projection so that parallel edges of
the cells stay parallel on screen.
