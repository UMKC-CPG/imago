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

What you have just built is the *production* build, the one every
shell uses unless told otherwise. Built as above, with `FC=h5fc` and
no extra options, it is an optimised serial engine: full compiler
optimisation, the serial HDF5, and no MPI. The word "production" names
the place the build is installed, `$IMAGO_DIR/bin`, and not a
particular way of compiling it. A second engine compiled differently
(in parallel, or with debugging checks) is kept beside it as a
*flavor* rather than built over it; see
[Other builds of the engine: flavors](#other-builds-of-the-engine-flavors)
below once the production build is working.

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

## Other builds of the engine: flavors

The steps above give you one engine, the production build in
`$IMAGO_DIR/bin`. Sometimes you want a second engine compiled
differently without giving up the first: a parallel build for a large
cell, a build full of runtime checks to locate a crash, a build with
profiling symbols to see where the time goes. Imago calls each of
these a *flavor*. You do not need one to follow the tutorials; come
back to this section when you do.

A flavor lives in its own directory, `$IMAGO_DIR/envs/<name>/bin`,
and never touches the production install. Only the compiled Fortran
programs differ from one flavor to the next. The Python scripts and
the databases do not depend on how the engine was compiled, so a
flavor's `bin` is a set of symbolic links to the production `bin/`
with the engine programs replaced by the flavor's own. Nothing is
copied twice, and a script you update in production is updated in
every flavor at the same moment.

### The flavors available

| Name      | What is compiled | Use it when |
|:----------|:-----------------|:------------|
|`serial`   | Optimised code, serial HDF5, no MPI. This is the same engine as a production build made by the steps above. | Production was built some other way (in parallel, say) and you want a serial engine beside it. On a standard install it adds nothing. |
|`mpi`      | Optimised code with MPI and the ELPA parallel eigensolver, plus profiling symbols. | The cell is large enough that one core is too slow. |
|`debug`    | Low optimisation, debugging symbols, every runtime check gfortran offers, and a backtrace on failure. | A run crashes and you need to know where. |
|`audit`    | The debug build plus floating-point exception traps, reals initialised to signaling NaN, and extra compiler warnings. | You suspect an uninitialised variable or a numerical error such as a division by zero. |
|`asan`     | The debug build plus AddressSanitizer, which includes leak detection. | You suspect memory is being overrun or leaked. |
|`profile`  | Optimised code plus symbols and stack frames; the generated code is the production code. | You want timings, or a profile from a tool such as callgrind. |
|`gprof`    | The profile build plus gprof instrumentation. | You want call counts or a call graph. Its timings are disturbed by the instrumentation, so take timings from `profile`. |

Each name corresponds to a CMake *preset* called `gfortran-<name>`,
defined in `CMakePresets.json` at the top of the repository. That
file is the place to look for the exact compiler options behind each
row.

### What a flavor needs

* **CMake 3.21 or later.** Presets are a newer CMake feature than the
  plain build in step 4 requires. Check with `cmake --version`.
* **The production scripts installed.** A flavor borrows its scripts
  from `$IMAGO_DIR/bin`, so complete step 4 first (or at least run
  `imago_scripts`). The engine half of a flavor does not depend on the
  production engine at all.
* **The compiler wrapper the preset names.** A preset chooses its own
  compiler and does not consult your `FC`: `h5fc` for every flavor
  except `mpi`, which uses `h5pfc`, the wrapper that an MPI-enabled
  HDF5 provides.

The `mpi` flavor needs more than the others: an MPI library, an HDF5
built against it, ScaLAPACK, and ELPA, with ELPA visible to
`pkg-config`. A serial HDF5 and an MPI-enabled HDF5 generally cannot
share one environment (see [Prerequisites](prerequisites.html)), so
the parallel libraries usually live in a second environment that you
activate in place of the serial one. With that environment active,
two checks tell you it is complete:

```
h5pfc -showconfig | grep -i "parallel hdf5"   # should say yes
pkg-config --modversion elpa                  # should print a version
```

### Building a flavor

The `imago_env` shell function, defined by `envs.sh` like
`imago_scripts`, does the whole job. To build the `mpi` flavor:

```
imago_env --add mpi
```

This configures and compiles the `gfortran-mpi` preset in its own
tree, `build/gfortran-mpi`, and then assembles `envs/mpi/bin`. It
finishes by reporting how many programs the flavor holds. Any other
name from the table works the same way, for example
`imago_env --add asan`. Run the command in the environment whose
libraries the flavor links against: the parallel environment for
`mpi`, the ordinary one for the rest.

Two messages are worth recognising. If `imago_env` notes that
`$IMAGO_DIR/bin` is empty, the flavor has an engine but no scripts;
run `imago_scripts` and add the flavor again. If it reports that the
flavor has no engine, the compile failed, and the cause is in the
CMake output just above the message.

### Using a flavor

Activating a flavor points `IMAGO_BIN`, `PATH` and `PYTHONPATH` at
its `bin` for the current shell only:

```
imago_env mpi       # this shell now runs the mpi engine
imago_env --list    # show production and every flavor that is built
imago_env           # back to the production build
```

`imago.py` starts whichever engine `$IMAGO_BIN` holds, so nothing
else changes in how you work. Because activation belongs to one
shell, two terminals can run two different flavors at once. A new
shell always starts on production, since `imagorc` sets `IMAGO_BIN`
afresh.

An engine finds its libraries when it starts, so the environment a
flavor was built in must also be active wherever it runs. For the
`mpi` flavor that means the parallel environment, both in your shell
and in a batch job.

Batch jobs deserve one more remark, because a job starts in a new
shell and would otherwise fall back to production. `makeinput.py`
handles this for you: when a flavor is active in the shell where you
run it, the submission file it writes selects the same flavor before
starting the engine, and stops the job at once if that flavor is not
built on the machine the job lands on. So activate the flavor
*before* running `makeinput.py`, not afterwards. When you also ask
`makeinput.py` for more than one core, the submission file starts
the engine through `mpirun` with that many processes; that request
only makes sense with the `mpi` flavor active.

### After changing the source

A flavor's engine programs are links into its build tree, so
recompiling that tree is enough to bring the flavor up to date:

```
cmake --build --preset gfortran-mpi
```

Running `imago_env --add mpi` again does the same compile and also
reassembles the flavor's `bin`, which is what you want when the set
of scripts or programs has changed.

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

`make install` in `build/base` rebuilds the production engine only.
If you keep flavors, rebuild each one you use after a pull as well,
or it will go on running the engine from before the pull:

```
imago_env --add mpi      # once per flavor in use
```

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
5. `echo $IMAGO_BIN` names the build you mean to run. It ends in
   `/bin` directly under `$IMAGO_DIR` for production and in
   `envs/<name>/bin` for a flavor. If `imago_env <name>` says the
   flavor is not built, add it with `imago_env --add <name>`.

When asking for help, send the output of `env | grep IMAGO`, the
errors that appeared on screen, and the `runtime` file from the job
directory if a run got that far.
