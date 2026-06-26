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
| `IMAGO_WARN_EXTRA` | the extra, audit-grade compiler warnings   |
| `IMAGO_SANITIZE`   | a sanitizer: `address`, `undefined`, or    |
|                    | `leak` (default `none`)                     |

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

So: "I want a leak-checking build" -> `gfortran-asan`; "I want
every warning the compiler can give" -> `gfortran-audit`.

Intel (`ifort`) presets are deferred: they need an `h5fc` that
wraps `ifort` (a matching HDF5 build), tracked as step 0c in
`dev/DEBUG.md`.

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
