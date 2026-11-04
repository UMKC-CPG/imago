---
layout: default
title: Prerequisites
nav_order: 1
---

# Prerequisites

Imago is built from source. Before cloning it, make sure the machine
has the three things the build needs: a Fortran compiler, an HDF5
library built with that same compiler, and a Python 3 interpreter
with a handful of packages. This page explains what each is for and
how to check it. On a shared cluster most of this is already
provided and your job is only to find it.

## A Fortran and C compiler

The engine is Fortran 90; the build also compiles a little C. Any of
`gfortran` or `ifort` works, with the matching C compiler. Check:

```
gfortran --version
gcc --version
```

Imago is built with CMake, version 3.1 or later, and GNU make. On
some machines the CMake executable is called `cmake3`.

```
cmake --version    # or: cmake3 --version
make --version
```

## HDF5, built with your compiler

Imago stores its matrices in HDF5 files and reads them back through
HDF5's Fortran interface. That interface produces Fortran module
files, and **module files from one compiler cannot be read by
another**. So the HDF5 library must have been built with the same
Fortran compiler you will build Imago with. This is the single most
common cause of a failed build.

HDF5 ships a compiler wrapper, `h5fc`, that knows where its headers
and libraries are. Imago's build is told to use that wrapper as the
Fortran compiler, which is how the two stay matched. Check what it
wraps:

```
h5fc -showconfig
```

In the output, look for three things:

* the **Fortran compiler** it was built with, which must be the same
  kind as yours (gfortran with gfortran, ifort with ifort);
* the **Fortran interface**, which must be enabled;
* the **deflate filter**, which must be on. Imago compresses its
  matrices with it.

If there is no `h5fc`, HDF5 is either absent or installed without
its Fortran interface. Ask your system administrator, install it
with conda (`conda install -c conda-forge hdf5`, which provides a
gfortran-built `h5fc`), or build it yourself in your home
directory. If you later use the parallel build of Imago, the
MPI-enabled HDF5 provides `h5pfc` instead; the two cannot coexist in
one conda environment, which is why the group keeps a separate
environment for the parallel toolchain.

## Linear algebra

The engine links against BLAS and LAPACK. Any implementation works:
the reference libraries, OpenBLAS, MKL or ATLAS. CMake finds them
automatically; if it cannot, point it at them when you configure.

## Python 3 and its packages

Every user-facing tool in Imago is a Python 3 script. Beyond the
standard library, the scripts import these packages:

| Package       | Used by |
|:--------------|:--------|
|numpy, scipy   | Structure handling, fitting, the k-point and integral tools. |
|pandas         | plotgraph and the relativistic orbital tools. |
|matplotlib     | plotgraph's default backend, plot_deadmd, hdf5Comp. |
|mplcursors     | Hover labels in plotgraph's matplotlib figures. |
|h5py           | hdf5Comp, the HDF5 comparison tool. |
|sympy          | The integral code generator osrecurintg. |
|vedo           | viewCell, the 3D cell and k-point viewer. |
|plotly         | Optional: plotgraph's plotly backend. |
|veusz          | Optional: plotgraph's veusz backend and make_veusz_graph. |
|ase            | cif2skl, which converts CIF files to skeletons. |

Install them into a virtual environment rather than the system
Python so that nothing else on the machine is disturbed:

```
python3 -m venv ~/imago-venv
source ~/imago-venv/bin/activate
pip install numpy scipy pandas matplotlib mplcursors h5py sympy vedo ase
```

Python 3.11 or later is required: some tools read TOML files with the
standard library's `tomllib`, which arrived in 3.11.

## Optional external programs

A few scripts drive other codes and need them on the path only when
you use those scripts: LAMMPS (`lmp`) and packmol for the
condensation and DEAD-MD tools, Grasp2K for the relativistic basis
tools, and Veusz, openDX, Paraview or OpenSCAD to view some of the
visualisation outputs.

## What to do next

With a compiler, an `h5fc` that wraps it, and the Python packages
installed, go on to [Installing Imago](installation.html).
