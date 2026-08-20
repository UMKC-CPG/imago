# Building Imago

This note covers how to build the Imago engine: the normal
production build, and the opt-in debugging builds used by the
bug-squashing campaign in `dev/DEBUG.md`. It is meant to be read
cold -- you should be able to pick "I want a leak-checking build"
without reading any CMake.

## Prerequisites (environment)

Imago is built with CMake and a Fortran compiler that **must match
the compiler that built the HDF5 library** -- the two cannot read
each other's `.mod` files. On this system the compiler is reached
through `FC=h5fc`, an HDF5 wrapper around gfortran. The build also
reads `IMAGO_DIR` to know where to install.

The normal session setup (the shell aliases expand to this):

- `cpg` -- activate the conda environment and set
  `LD_LIBRARY_PATH`.
- `simago` -- activate the project virtualenv and `source
  .imago/imagorc` (which sets `IMAGO_DIR`, `IMAGO_BIN`, etc.).

With that done, `h5fc`, `cmake`, and the databases are all on the
path.

## Normal (production) build

The hand-driven production trees live in `build/release` and
`build/debug`. From a configured session:

```
cdrelease            # cd $IMAGO_DIR/build/release
cmake ../..          # configure (defaults to a RELEASE build)
make install         # build and install into $IMAGO_DIR/bin
```

For a debug build use `build/debug` and
`cmake -DCMAKE_BUILD_TYPE=Debug ../..`. These two trees, and the
installed `bin/`, are the known-good production artifacts.

The two executables produced are `imago` (the general multi-k,
complex build) and `imagoG` (the gamma-point-only, real build,
which is faster and uses about half the memory). `imago.py`
selects between them automatically from the k-point files.

## Debugging builds (the bug-squashing campaign)

The top-level `CMakeLists.txt` adds opt-in instrumentation options,
all **default OFF**. With every option off, the build is
byte-for-byte identical to the production build above -- the
instrumentation is only ever *added*, never substituted, so the
known-good default is never disturbed.

| Option             | Adds                                       |
| ------------------ | ------------------------------------------ |
| `IMAGO_CHECKS`     | runtime checks: array bounds, etc.         |
| `IMAGO_FPE_TRAP`   | trap FP exceptions (invalid / divide-by-   |
|                    | zero / overflow) at the point of origin    |
| `IMAGO_INIT_SNAN`  | initialize reals to signaling NaN, so any  |
|                    | use of an uninitialized value traps        |
| `IMAGO_WARN_EXTRA` | the extra, audit-grade compiler warnings,  |
|                    | CORRECTNESS diagnostics only               |
| `IMAGO_WARN_PERF`  | performance warnings: hidden array copies. |
|                    | Kept apart from the audit set on purpose;  |
|                    | see the note below                         |
| `IMAGO_SANITIZE`   | a sanitizer: `address`, `undefined`, or    |
|                    | `leak` (default `none`)                     |
| `IMAGO_PROFILE`    | `-g -fno-omit-frame-pointer`: symbols and  |
|                    | stack frames for profilers, code unchanged |
| `IMAGO_GPROF`      | `-pg` on compile and link: gprof hooks     |
|                    | (call counts / call graph; perturbs time)  |
| `IMAGO_MPI`        | `-DIMAGO_MPI`; FC must be an MPI wrapper   |
|                    | (`h5pfc`); see "Parallel toolchain"        |
| `IMAGO_ELPA`       | `-DIMAGO_ELPA` + ELPA/ScaLAPACK include and|
|                    | link lines via pkg-config; needs IMAGO_MPI |

**Why the two warning options are separate.** `IMAGO_WARN_PERF`
turns on `-Warray-temporaries`, which reports hidden array copies.
Those are a COST, not a defect. One sweep produced 187 of them
against 4 genuine correctness warnings, so bundling the two makes
the correctness warnings impossible to find -- which is the whole
purpose of the audit build. Turn `IMAGO_WARN_PERF` on when you are
hunting for speed, and leave it off when you are hunting for bugs.

You can compose these by hand onto any build, e.g.:

```
cmake -DCMAKE_BUILD_TYPE=Debug -DIMAGO_SANITIZE=address ../..
```

### Presets (the easy way)

`CMakePresets.json` bundles the common combinations. Each preset
builds into its **own** tree `build/<preset-name>/` and installs to
a throwaway prefix inside that tree, so a preset can never overwrite
the production install or the `build/release` / `build/debug` trees.

```
cmake --preset gfortran-asan     # configure
cmake --build --preset gfortran-asan   # build (or: make -C build/gfortran-asan)
```

| Preset             | What you get                               |
| ------------------ | ------------------------------------------ |
| `gfortran-release` | optimized, no instrumentation (reference)  |
| `gfortran-debug`   | the historical debug flags                 |
| `gfortran-audit`   | debug + checks + FP traps + signaling NaN  |
|                    | + extra warnings (static/runtime hunting)  |
| `gfortran-asan`    | debug + AddressSanitizer (+ leak) + checks |
| `gfortran-profile` | RELEASE code + symbols + frame pointers    |
|                    | (timings and callgrind; dev/PERFORMANCE.md)|
| `gfortran-gprof`   | profile preset + `-pg` (gprof call counts  |
|                    | and call graph; its timings are perturbed) |

So: "I want a leak-checking build" -> `gfortran-asan`; "I want
every warning the compiler can give" -> `gfortran-audit`; "I want
to know where the time goes" -> `gfortran-profile` under
`valgrind --tool=callgrind`, and `gfortran-gprof` only when a
call count or the call graph is the question.  The profile
preset's binaries run at release speed, so a wall-clock baseline
taken with them is a baseline for the production build.

Intel (`ifort`) presets are deferred: they need an `h5fc` that
wraps `ifort` (a matching HDF5 build), tracked as step 0c in
`dev/DEBUG.md`.

## Parallel toolchain (the `cpgp` environment)

The MPI + ELPA build (ARCHITECTURE 6.5-6.6, `dev/PERFORMANCE.md`)
compiles against a second conda environment, **`cpgp`** -- the
group's general `cpg` env plus a trailing `p` for *parallel*.
It exists because conda cannot hold the serial (`nompi`) and the
MPI builds of HDF5 in one environment: `cpg` keeps the serial
HDF5 that the serial binaries and the performance baselines were
built with, and `cpgp` carries the MPI-enabled HDF5 (`h5pfc`, a
wrapper around `mpifort`), openmpi, ScaLAPACK and ELPA, on the
same gfortran major and the same OpenBLAS family, so a binary
from `cpgp` differs from a `cpg` build only in the parallel
libraries. The recipe is `dev/env/cpgp.yml`:

```
mamba env create -f dev/env/cpgp.yml     # first time
conda activate cpgp
h5pfc -showconfig | grep -i "parallel hdf5"   # -> yes
```

The Python side does not move: imago's venv is self-contained
and links neither HDF5 nor MPI, so one venv serves both envs.

Build with the preset, which sets `FC=h5pfc` and turns on
`IMAGO_MPI` and `IMAGO_ELPA` (plus `IMAGO_PROFILE`, so parallel
runs profile like serial ones):

```
conda activate cpgp
cmake --preset gfortran-mpi
cmake --build --preset gfortran-mpi
```

`IMAGO_MPI` refuses to configure unless `FC -show` expands to an
MPI command line -- the "compiler must match HDF5" rule of the
prerequisites section, applied to the parallel HDF5. `IMAGO_ELPA`
finds ELPA through pkg-config (`elpa.pc` also names ScaLAPACK,
LAPACK and BLAS on its link line). Neither option adds compiler
flags of its own beyond `-DIMAGO_MPI` / `-DIMAGO_ELPA` and the
ELPA include path; the MPI flags come from the wrapper.

Launch (verified on two nodes, 2026-08-19): `mpirun -np
"$SLURM_NTASKS" …` inside the batch script -- the group's LAMMPS
precedent -- spans nodes correctly, and the cross-node ELPA
handshake (`elpa_init` / `elpa_allocate`) succeeds under it.
`srun --mpi=pmix` behaves identically and is the working
alternative. Do NOT use `srun --mpi=pmi2`: OpenMPI 5 dropped
PMI-2 support, so it silently launches N independent one-rank
worlds AND returns exit code 0 -- only a rank/size printout
betrays the fragmentation. (A benign PMIx stderr warning about
a missing `munge` component accompanies every launch; PMIx
falls back to its native security plugin and the runs are
unaffected.) Note also that `mpirun` inside a SLURM job takes
its slot count from the allocation: request `-n` at least as
large as the widest `-np` you will launch, or mpirun refuses
with "not enough slots" before any rank starts.

Logs under many ranks: rank 0 writes `fort.20` exactly as the
serial program does, and it is the only log a normal parallel
run leaves -- the other ranks' copies of the same output are
discarded. Set `IMAGO_RANK_LOGS=1` in the environment to give
every worker rank its own `fort.20.rNNNN` file instead, for
debugging sessions where one rank's view of events is the
question. Error messages never depend on this switch: aborts
are rank-stamped and written to standard error, which the
launcher aggregates into the job's error file.

## Build flavors: running an instrumented binary on a real job

The instrumented builds are never installed over the production
`bin/`. Instead they live as **flavors** under
`$IMAGO_DIR/envs/<flavor>/bin`, and you switch the active toolchain
to one with a single command. A flavor's bin is just symlinks to
the production `bin/` with the engine executables (`imago`,
`imagoG`, ...) overlaid by the flavor build -- so every helper
script and the `share/` database are reused unchanged, nothing is
duplicated, and the production install is never touched.

The switcher is `envs.sh`, installed to `bin/`. Source it once in
your shell startup, right after the Imago rc file:

```
source "$IMAGO_DIR/.imago/imagorc"
source "$IMAGO_BIN/envs.sh"
```

Then:

```
imago_env --build asan    # build the gfortran-asan preset and
                          #   assemble envs/asan/bin (one-time / on
                          #   rebuild)
imago_env asan            # activate it: repoints IMAGO_BIN at the
                          #   asan flavor for this shell
imago.py                  # now runs the instrumented engine
imago_env                 # back to the production install
imago_env --list          # show which flavors are built
```

`imago.py` launches the engine as `$IMAGO_BIN/imagoG`, so activating
a flavor selects which build runs with no code changes. Activation
is per-shell, so different terminals can run different flavors at
once. For an `asan` flavor, `envs.sh` sets a sensible default
`ASAN_OPTIONS` (leak detection on) that you can override. `imago.py`
also has a built-in `-valgrind` flag for running the stock binary
under valgrind without any rebuild.

A worked example (a full Gamma SCF under AddressSanitizer) is
recorded in `dev/DEBUG.md` (Phase 0d).
