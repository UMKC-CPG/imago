---
layout: default
title: hdf5Comp.py
script: hdf5Comp.py
description: 'Debugging tool that compares one matrix dataset in an Imago HDF5 file against another, in the same or a different file, prints both and their difference, and shows the difference as a colour map. Can unpack packed-triangle real, imaginary or magnitude storage.'
parent: Script Catalog
---

# hdf5Comp

A debugging tool for direct comparison of the matrices stored in
Imago's HDF5 files. Imago keeps its integral and Hamiltonian
matrices in HDF5 groups keyed by k-point and interaction type, and
when two builds, two compilers, or two code paths are suspected of
disagreeing, this is the quickest way to see where. It opens two
datasets, which may live in the same file or in two files in two
directories, prints both matrices and their difference to the
terminal, and opens a colour-map image of the difference.

The default group and dataset names point at the atom-overlap
matrices of the first k-point in a post-SCF file, comparing dataset
`0000001` with `0000002`; change them to any group and dataset path
present in the file. Use `h5ls -r` on the file to see what is
available.

**Defaults.** Every option's default value is read from the resource
control file `hdf5Comprc.py` in the directory named by `$IMAGO_RC`.

## Usage

**Commandline:**

```
hdf5Comp.py [-f1 FILE] [-f2 FILE] [-d1 DIR] [-d2 DIR]
            [-g1 GROUP] [-g2 GROUP] [-s1 DATASET] [-s2 DATASET]
            [-u] [-i] [-m] [--help]
```

**Option Table:**

| Option              | Defaults                                     | Description |
|:-------------------:|:--------------------------------------------:|:------------|
|-f1 / --file1        |gs_pscf-fb.hdf5                               | First HDF5 file. |
|-f2 / --file2        |gs_pscf-fb.hdf5                               | Second HDF5 file. |
|-d1 / --dir1         |.                                             | Directory holding the first file. |
|-d2 / --dir2         |.                                             | Directory holding the second file. |
|-g1 / --group1       |00001_00001_00001/atomIntgGroup/atomOverlap   | Group path of the first dataset. |
|-g2 / --group2       |00001_00001_00001/atomIntgGroup/atomOverlap   | Group path of the second dataset. |
|-s1 / --dataset1     |0000001                                       | Name of the first dataset. |
|-s2 / --dataset2     |0000002                                       | Name of the second dataset. |
|-u / --unpack        |Off                                           | The datasets are packed upper triangles of symmetric or Hermitian matrices; unpack them into full square matrices before comparing. |
|-i / --imaginary     |Off                                           | With `-u`, take the imaginary part of the packed complex matrix (antisymmetric on unpacking). |
|-m / --magnitude     |Off                                           | Take the absolute magnitude of the complex matrix, i.e. the square root of real squared plus imaginary squared. With `-u` this reads the packed pair; without it, it reads separate `real<dataset>` and `imag<dataset>` datasets. |
|--help               |N/A                                           | Print usage information and exit. |

## Output

The two matrices and their difference (second minus first) are
printed to standard output, then a matplotlib window shows the
difference as an image with a colour bar. A uniformly dark image
means agreement; a structured pattern locates the disagreeing
blocks, which for an atom-indexed matrix points at the atom pairs
involved. Requires `h5py`, `numpy` and `matplotlib`.
