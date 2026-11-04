---
layout: default
title: Installing Imago
nav_order: 2
---

# Installing Imago

This tutorial installs Imago under your own home directory, which is
the normal arrangement for a student or a single researcher. A group
install differs only in where the clone lives: one person clones
and builds it in a shared location, and everyone else sources that
install's `imagorc` file instead of their own. The steps below are
the same in both cases.

By the end you will have the engine compiled, every script on your
path, the databases unpacked, and a shell that knows where all of it
is. It assumes the [Prerequisites](prerequisites.html) are met.

## 1. Clone the repository

```
cd ~
git clone git@github.com:umkc-cpg/imago.git
```

or, over HTTPS, `git clone https://github.com/umkc-cpg/imago.git`.
This creates `~/imago`. The environment file assumes that name and
location by default, so keeping it saves an edit later. If you clone
elsewhere, note the path; you will put it in the environment file.

## 2. Set up the environment file

Everything Imago does is located through a handful of environment
variables, and all of them are set by one file that you source,
`~/imago/.imago/imagorc`. Open it and check the first few
definitions:

| Variable      | Default                       | Meaning |
|:--------------|:------------------------------|:--------|
|`IMAGO_BRANCH` | `imago`                       | The name of the clone directory. |
|`IMAGO_DIR`    | `$HOME/$IMAGO_BRANCH`         | The root of the install: the clone. Change this if you cloned elsewhere. |
|`IMAGO_RC`     | `$IMAGO_DIR/.imago`           | Where the scripts look for their resource-control files, the `*rc.py` files holding each script's defaults. |
|`IMAGO_BIN`    | `$IMAGO_DIR/bin`              | Where the engine and the scripts are installed. Added to your `PATH`. |
|`IMAGO_DATA`   | `$IMAGO_DIR/share`            | The databases: basis sets, potentials, space groups, band paths. |
|`IMAGO_TEMP`   | `$HOME/data/scratch/imago`    | Scratch space for the engine's large working files. Point this at a filesystem with room; a run can leave tens of megabytes here and a campaign gigabytes. |
|`IMAGO_QUEUE`  | `3`                           | The batch system: 0 for plain bash scripts, 1 PBS, 2 LSF, 3 SLURM. Controls what kind of submission script `makeinput.py` writes. |

The `GRASP_DIR` and `VASPPOT_*` variables point at external packages
and can be left alone until you need them.

Create the scratch directory the file names:

```
mkdir -p ~/data/scratch/imago
```

Then make your shell source the file at login by adding one line to
`~/.bashrc`:

```
source "$HOME/imago/.imago/imagorc"
```

The file ends by sourcing `envs.sh` from the same directory, which
defines two helper functions you will meet below. Open a new shell,
or source `~/.bashrc`, and confirm:

```
echo $IMAGO_DIR $IMAGO_BIN $IMAGO_DATA
```

## 3. Tell the build which compiler to use

The build reads the `FC` environment variable for the Fortran
compiler, and it must name the HDF5 wrapper so that the engine links
against an HDF5 built with the same compiler:

```
export FC=h5fc
```

Put this line in `~/.bashrc` as well if `h5fc` is always the right
choice on your machine. On a machine with both a serial and a
parallel HDF5, the wrapper differs between the two environments, so
set `FC` in the shell where you build rather than globally.

## 4. Configure and build

Imago builds out of source, in a directory under `build/`. The
production tree is `build/base`:

```
mkdir -p $IMAGO_DIR/build/base
cd $IMAGO_DIR/build/base
cmake ../..
make install
```

`cmake ../..` configures an optimised build by default and reads
`IMAGO_DIR` to learn where to install. `make install` compiles the
engine and copies the executables, every script, and the
resource-control files into `$IMAGO_DIR/bin` and `$IMAGO_DIR/.imago`.
Two engine executables are produced: `imago`, the general build for
any k-point mesh, and `imagoG`, the Gamma-point-only build that is
faster and uses about half the memory. `imago.py` picks between them
for you.

If you only ever change scripts afterwards, there is no need to
recompile. The `imago_scripts` shell function that `envs.sh` defined
reinstalls just the Python side:

```
imago_scripts
```

## 5. Unpack the databases

The basis sets, potentials, space groups and band paths ship as
compressed archives in `share/`. Unpack them once:

```
unpackImagoDB.py
```

This takes a few seconds and prints nothing on success. Afterwards
`$IMAGO_DATA` holds `atomicBDB/`, `atomicPDB/`, `precursorDB/`,
`spaceDB/` and `sybdDB/`.

## 6. Check the install

Three quick checks, each confirming one piece:

```
which imago.py makeinput.py      # both under $IMAGO_BIN
ls $IMAGO_BIN/imago $IMAGO_BIN/imagoG   # the two engines exist
ls $IMAGO_DATA/atomicBDB | head   # element directories
```

The Python test suite exercises the scripts without running the
engine. With the virtual environment active and `IMAGO_DATA` set:

```
cd $IMAGO_DIR
python src/tests/test_o.py -m "not fortran"
```

The first real test of the engine is the next tutorial,
[Your first calculation](first_calculation.html), which runs a
small silicon crystal end to end.

## Keeping it up to date

Imago is developed in the open, and pulling the latest version is a
three-step routine:

```
cd $IMAGO_DIR
git pull
cd build/base && make install
```

`make install` re-runs CMake if the source tree's build files
changed. Re-run `unpackImagoDB.py` only if the `share/*.tgz`
archives changed in the pull; `git pull` lists the files it updated.

## If something is wrong

Nearly every problem on first use is an environment variable.
Check, in this order:

1. `echo $IMAGO_DIR` prints your clone. If it is empty, the
   `imagorc` line is not in `~/.bashrc` or the shell has not been
   restarted.
2. `echo $FC` prints `h5fc` and `which h5fc` finds it. If the build
   failed with errors about `.mod` files, HDF5 was built with a
   different compiler.
3. `echo $IMAGO_TEMP` names a directory that exists and is writable.
   `imago.py` creates its scratch there.
4. `python3 -c "import numpy, pandas, matplotlib"` succeeds in the
   shell you run scripts from. The virtual environment must be
   activated **before** sourcing `imagorc`, because activation
   resets `PYTHONPATH` and `imagorc` needs to add `$IMAGO_BIN` to
   it afterwards.

When asking for help, send the output of `env | grep IMAGO`, the
errors that appeared on screen, and the `runtime` file from the job
directory if a run got that far.
