---
layout: default
title: Setup
nav_order: 2
---
# Setting up Imago

Imago is intended to be cloned from GitHub, built on the local
machine, and run from source. There are no downloadable release
builds. This page is the condensed version of the setup process; the
step-by-step guides are the tutorials:

- [Prerequisites](../tutorials/prerequisites.html): what the machine
  needs before Imago can be built, and how to check each item.
- [Installing Imago](../tutorials/installation.html): a complete
  install under one user's home directory. A group install follows
  the same steps in a shared location, with each member sourcing
  that install's `imagorc`.

The install process in outline:

- Setting up the compilers and HDF5
- Cloning from GitHub
- Setting up the user environment
- Compiling and installing Imago
- Optionally, building other flavors of the engine
- Validating the install
- Updating Imago

## Prerequisites

##### Compilers and build environment

* A C compiler.
* A Fortran compiler capable of Fortran 90/95 (gfortran or ifort).
* CMake (version 3.1 or later) and GNU make. The CMake executable
  may be `cmake` or `cmake3` depending on your machine.

##### Python

* Python 3.11 or later.
* numpy, scipy, pandas, matplotlib, mplcursors, h5py, sympy, vedo
  and ase.
* Optional: plotly and veusz, for plotgraph's alternative backends.

Perl is not required. Earlier versions of the code used it; every
tool is now Python or bash.

##### Other requirements

* Linear algebra libraries: BLAS and LAPACK in any implementation
  (reference, OpenBLAS, MKL, ATLAS).
* HDF5, including the Fortran interface and the deflate filter.

***IMPORTANT:*** the Fortran compiler used to build HDF5 must be the
same one you build Imago with. The two compilers' module files are
not interchangeable.

## Setting up Compilers

We assume you have a Fortran and C compiler, whether a system-wide
gfortran and gcc or a compiler installed with conda. Either way you
need HDF5 correctly wrapping that compiler. Run `h5fc -showconfig`
and look for three things: that a "deflate" filter is on, that the
Fortran interface is included, and which compilers built HDF5. It is
vital that the compiler used to build HDF5 be the same kind you have
on your system. If HDF5 was built with gfortran but your system has
ifort, or the reverse, you will have a problem, because the two
Fortran compilers cannot read each other's module files.

If you do not have HDF5, ask your system administrators to install
it, install it with conda, or build it yourself in your `$HOME`. In
any case make sure the deflate filter is on and the Fortran
interface is included. The Imago build is told to use `h5fc` as its
Fortran compiler through the `FC` environment variable, which is how
the match is enforced:

```
export FC=h5fc
```

For Python, use a virtual environment and install the packages
listed above into it with `pip`.

## Cloning from GitHub

Once you have the dependencies, clone the repository. If you already
know how, use your preferred method and move on; otherwise:

- Clone to your desired location, ideally `$HOME`, over SSH or
  HTTPS:
   - `git clone git@github.com:umkc-cpg/imago.git`
   - `git clone https://github.com/umkc-cpg/imago.git`
- This creates an `imago/` directory. For a different name, add it
  to the command: `git clone <url> <directoryName>`. The environment
  file assumes `$HOME/imago` by default, so another location means
  one edit in the next step.

## Setting up the User Environment

Regardless of your scenario, each user should set up an "Imago
friendly" environment. It makes setup and use much easier.

- Go to the hidden `.imago` directory in the repository.
- Edit the `imagorc` file (not any `*rc.py` file) so that its
  environment variables point to the right directories:
   - `IMAGO_BRANCH` is the name of the clone directory, `imago` by
     default.
   - `IMAGO_DIR` is the full path to the clone, `$HOME/$IMAGO_BRANCH`
     by default.
   - `IMAGO_RC` is where the scripts look for their `*rc.py` default
     files; it may differ from `IMAGO_DIR/.imago` if you keep
     personal defaults elsewhere.
   - `IMAGO_TEMP` is scratch space for the engine's large working
     files; point it at a filesystem with room and create it.
   - `IMAGO_QUEUE` is the batch system: 0 bash, 1 PBS, 2 LSF, 3
     SLURM.
- Edit your `.bashrc` or bash profile to source that `imagorc` file.
  It sources `envs.sh` in turn, which defines the `imago_env` and
  `imago_scripts` helper functions.
- Optionally add aliases such as these:
     + `alias cdbase='cd $IMAGO_DIR/build/base/'`
     + `alias cddebug='cd $IMAGO_DIR/build/debug/'`
     + `alias cdsrc='cd $IMAGO_DIR/src/imago/'`
     + `alias cdjobs='cd $IMAGO_DIR/jobs/'`
     + `alias interactive='srun -p interactive --mem=8G --pty /bin/bash'`
     + `alias interactive_login='srun -p interactive --mem=20G --pty /bin/bash --login'`
     + `alias interactive_x11='srun -p interactive --x11=first --mem=20G --pty /bin/bash --login'`
     + `alias la='ls -latr'`
     + `alias shist='sacct --starttime $(date -d yesterday +%D-%R) --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss,MaxVMSize,nnodes,ncpus,nodelist'`

If you use a Python virtual environment, activate it **before**
sourcing `imagorc`: activation resets `PYTHONPATH`, and `imagorc`
needs to add `$IMAGO_BIN` to it afterwards.

After editing, source your profile or restart the terminal.

## Compiling and Installing Imago

Imago builds out of source. The production tree is `build/base`; a
debug build goes in `build/debug`. The first time:

```
mkdir -p $IMAGO_DIR/build/base
cd $IMAGO_DIR/build/base
cmake ../..          # or: cmake -DCMAKE_BUILD_TYPE=Debug ../..
make install
unpackImagoDB.py
```

`make install` builds the engine and installs it, every script and
the resource-control files into `$IMAGO_DIR/bin` and
`$IMAGO_DIR/.imago`. `unpackImagoDB.py` unpacks the basis, potential,
space-group and band-path databases in `$IMAGO_DATA` and takes a few
seconds.

Later, after editing source files, go to the build directory and run
`make install` again; it re-runs CMake if needed. If you only
changed scripts, `imago_scripts` reinstalls the Python side with no
compile at all.

## Building Other Flavors of the Engine

The build above is the *production* build. Made as shown, with
`FC=h5fc` and no extra options, it is an optimised serial engine: full
compiler optimisation, the serial HDF5, and no MPI. "Production" names
where the build is installed, `$IMAGO_DIR/bin`, and not a particular
way of compiling it. A second engine compiled differently is kept
beside it as a *flavor* under `$IMAGO_DIR/envs/<name>/bin`, never
installed over production. Only the Fortran programs differ between
flavors, so a flavor's `bin` is symbolic links to the production
`bin/` with the engine programs replaced by its own. Flavors are
optional; the production build is all the tutorials need.

| Name      | What is compiled |
|:----------|:-----------------|
|`serial`   | Optimised, serial HDF5, no MPI. The same engine as a production build made as above; useful only when production was built some other way. |
|`mpi`      | Optimised, with MPI and the ELPA parallel eigensolver. |
|`debug`    | Low optimisation with gfortran's runtime checks and a backtrace. |
|`audit`    | Debug plus floating-point traps, signaling-NaN initialisation and extra warnings. |
|`asan`     | Debug plus AddressSanitizer, with leak detection. |
|`profile`  | Optimised code plus symbols, for timings and profilers. |
|`gprof`    | Profile plus gprof instrumentation, for call counts. |

Each name is the CMake preset `gfortran-<name>` in
`CMakePresets.json`. The `imago_env` function builds, updates, lists
and selects them:

```
imago_env --add mpi     # once: compile the preset, assemble envs/mpi/bin
imago_env --update mpi  # later: recompile after a source change or pull
imago_env --update      # the same, for every flavor that is built
imago_env mpi           # use it in this shell (repoints IMAGO_BIN)
imago_env --list        # production and every flavor that is built
imago_env               # back to production
```

Points to know before building one:

* Presets need CMake 3.21 or later, newer than the plain build does.
* A preset picks its own compiler wrapper and ignores `FC`: `h5fc`
  for every flavor except `mpi`, which uses `h5pfc`.
* The `mpi` flavor additionally needs an MPI library, an MPI-enabled
  HDF5, ScaLAPACK and ELPA (found through `pkg-config`). These usually
  sit in a separate environment from the serial HDF5; build and run
  the flavor with that environment active.
* A flavor borrows its scripts from `$IMAGO_DIR/bin`, so install the
  production build, or at least run `imago_scripts`, first.
* `--update` must be run in the environment the flavor was built in.
  With no name it rebuilds what it can, leaves the others untouched,
  and lists the flavors that failed so they can be updated by name.
* Activation lasts for one shell. A batch job starts a new shell, so
  activate the flavor *before* running `makeinput.py`; the submission
  file it writes then selects the same flavor inside the job.

The [installation tutorial](../tutorials/installation.html) walks
through each of these. `BUILD.md` at the repository root adds the
developer detail: the individual instrumentation options behind the
presets and the group's own compiler environments.

## Validating the Imago Installation

* `unpackImagoDB.py` should run for a few seconds and return with no
  output. Errors here mean `$IMAGO_DATA` is wrong or the scripts are
  not on your path.
* `which imago.py` should find it under `$IMAGO_BIN`, and
  `ls $IMAGO_BIN/imago $IMAGO_BIN/imagoG` should list both engines.
* The Python test suite exercises the scripts without the engine.
  With your virtual environment active, run from `$IMAGO_DIR`:
  `python src/tests/test_o.py -m "not fortran"`.
* The real test is the
  [first calculation tutorial](../tutorials/first_calculation.html),
  which runs silicon end to end in under a minute.

If things are not working at this point, contact your system
administrator or Dr. Rulis. Provide your environment
(`env | grep IMAGO`), a copy of any errors dumped to the screen, the
`runtime` file if present, the `intermediate/fort.20` file if
present, and the directory where you are doing the work.

## Updating Imago

If you cloned from this repository and are not doing development
yourself, keeping up to date is simple:

* `cd $IMAGO_DIR`.
* Check that your remote is this repository with `git remote -v`; it
  should show `github.com:umkc-cpg/imago.git` or the HTTPS
  equivalent.
* You are on the `main` branch (`git checkout main` if not).
* `git fetch origin` shows what has changed; `git pull` merges it.
* Then `cd build/base` and `make install`. Re-run `unpackImagoDB.py`
  only if the pull updated a `share/*.tgz` archive.
* If you keep flavors, update each one you use with
  `imago_env --update <name>`; `make install` in `build/base`
  refreshes the production engine only.

You are now up to date.
