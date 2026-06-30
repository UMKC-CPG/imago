# Architecture

> **Document hierarchy:** VISION -> **ARCHITECTURE** -> DESIGN
> -> PSEUDOCODE -> Code. For goals and principles, see
> `VISION.md`.

---

## 1. Repository Layout

```
src/
  imago/              Primary production program
    imago.F90          Top-level dispatcher
    kpoints.f90        O_KPoints module (mesh, weights, phases)
    dos.F90            O_DOS module (TDOS, PDOS)
    bond.F90           Bond order (Mulliken overlap population)
    populate.F90       Electron population / occupancy
  makeKPoints/
    makekpoints.F90    Standalone k-point mesh generator
  scripts/             Python tools: flat .py files are CLI
                       entry points; subdirectories are Python
                       packages
    imago.py           Imago driver: CLI and callable API (9.2)
    makeinput.py       Input file orchestrator
    makegroups.py      Bispectrum type-grouping helper (8.9,
                       DESIGN 5.10): CLI + importable
    matchers.py        Matcher protocol library: Reduce/Bispec
                       matchers + MATCHERS registry (8.9)
    ase_imago.py       ASE Calculator, ImagoCalculator (9.3)
    cod_fish.py        COD acquisition front-end -> CIF (9.5)
    cif2skl.py         CIF -> imago.skl converter (9.5)
    kaleidoscope/      Parsl flight dispatcher, "kaleidoscope" (9.4)
  kinds.f90            Shared precision kinds
  constants.f90        Shared physical/mathematical constants
dev/
  VISION.md            Goals and principles
  ARCHITECTURE.md      This document
  DESIGN.md            Algorithmic design
  PSEUDOCODE.md        Algorithm specifications
  TODO.md              Task list by level
```

---

## 2. Module Map

Modules directly affected by the current development:

- **O_KPoints** (`kpoints.f90`): Stores the k-point mesh,
  weights, and phase factors. Declares `kPointIntgCode`
  (0=Gaussian, 1=LAT). LAT data already implemented:
  `numTetrahedra`, `tetraVol`, `tetrahedra(:,:)`,
  `fullKPToIBZKPMap(:)`, `generateTetrahedra`,
  `computeTetraVol`, `initializeKPointMesh` (with IBZ
  folding). Phase F additions for correct IBZ unfolding
  of eigenvector-dependent quantities:
  `fullKPToIBZOpMap(:)` (which point group operation
  mapped each full-mesh k-point to its IBZ
  representative), `atomPerm(:,:)` (atom permutation
  table under each point group operation), and
  `buildAtomPerm` (builds the table from point ops and
  fractional atom positions).
- **O_DOS** (`dos.F90`): Two DOS paths: `computeIterationTDOS`
  (in-SCF convergence monitoring) and `computeDOS` (full
  TDOS/PDOS post-processing). Will gain a LAT branch in
  `computeDOS` dispatched on `kPointIntgCode`.
- **bond.F90**: Computes bond order and effective charge
  (Q*). Already uses `electronPopulation_LAT` when
  `kPointIntgCode == 1`. Phase F additions: the
  accumulation loop will buffer per-IBZ-kpoint
  projections and distribute them across the star of
  each IBZ k-point using `atomPerm` and
  `fullKPToIBZOpMap` from O_KPoints, giving correct
  per-atom Q* and per-pair bond order under IBZ
  reduction.
- **populate.F90** (`O_Populate`): Computes
  `electronPopulation` (Gaussian/Fermi-filling path)
  which folds in `kPointWeight`.
  `computeElectronPopulation_LAT` (already implemented)
  produces `electronPopulation_LAT` using tetrahedra
  data from O_KPoints. Both occupation arrays serve the
  same downstream consumers (bond order, effective
  charge); the integration method (`kPointIntgCode`)
  selects which one is used. Upstream of bond order and
  effective charge.
- **makeKPoints** (`makekpoints.F90`): Legacy standalone
  program that generates explicit k-point lists with IBZ
  reduction (`foldMesh`). Separate executable, not linked
  into imago. Being retired -- its mesh-building and IBZ
  reduction capabilities have been ported into imago's
  `initializeKPointMesh`. Retained only for backward
  compatibility; new workflows should not depend on it.
- **makeinput.py** (`src/scripts/makeinput.py`): Top-level
  orchestrator that prepares all input files for Imago. For
  k-points it supports three pipelines:
  1. *Mesh mode* (`-kp`, `-scfkp`, `-pscfkp`): writes a
     style-code-1 k-point file with axial counts and shift.
     Imago builds the full mesh internally and reduces to
     the IBZ, giving it full access to the mesh topology
     and symmetry maps needed for correct decomposition
     properties.
  2. *Density mode* (`-kpd`, `-scfkpd`, `-pscfkpd`):
     writes a style-code-2 k-point file with density and
     shift. Imago computes axial counts from the density
     and reciprocal cell geometry, then builds and reduces
     the mesh identically to mesh mode.
  3. *Explicit mode* (style code 0 in the k-point file):
     a pre-built list of k-points with weights, read
     directly by Imago with no internal mesh construction.
     Supported for special cases (e.g., hand-crafted
     k-point sets), but Imago emits a prominent warning
     that decomposition properties (effective charge,
     bond order, PDOS) will not be correct unless the
     user has taken extreme care to provide a symmetric
     mesh. Not produced by makeinput.
- **makegroups.py** (`src/scripts/makegroups.py`): the
  bispectrum type-grouping helper (8.9, DESIGN 5.10).  A
  Fortran-side descriptor can only come from a completed
  Imago run, so grouping atoms by bispectrum is a *sequence*
  -- makeinput (no grouping) -> `imago -loen -scf no` ->
  read the enriched `fort.21` -> bucket with `BispecMatcher`
  -> rewrite the skeleton with explicit per-element species
  tags -- orchestrated outside `makeinput.py`.  This script
  owns that sequence and is **dual-mode**: an importable
  `group_by_bispectrum(...)` the producer
  (`build_initial_potentials.py`) calls, plus a `__main__`
  CLI for manual use.  `makeinput.py` stays a plain
  input-writer that reads whatever explicit types the
  skeleton already carries; the bispectrum reasoning never
  enters it.  Reduce grouping remains in `makeinput.py` for
  now (the two environment schemes are deliberately
  asymmetric); the long-term direction is for `makegroups.py`
  to absorb *all* type assignment (reduce, target, block),
  leaving `makeinput.py` purely an input-writer.

---

## 3. Dependency Graph

```
makeinput.py (top-level orchestrator)
  +-- (mesh mode: -kp)
  |     writes kp-scf.dat / kp-pscf.dat (style 1)
  +-- (density mode: -kpd)
  |     writes kp-scf.dat / kp-pscf.dat (style 2)
  +-- makeKPoints (legacy, no longer called by makeinput)

imago.F90 (top-level dispatcher)
  +-- O_KPoints (kpoints.f90)
  |     +-- O_Lattice (recipCellVolume, invRealVectors)
  |     +-- O_AtomicSites (atom positions for buildAtomPerm)
  +-- O_DOS (dos.F90)
  |     +-- O_KPoints (eigenvalues, tetrahedra, weights)
  |     +-- O_PSCFIntg / O_SCFIntg (eigenvectors via HDF5)
  +-- populate.F90 (O_Populate)
  |     +-- O_KPoints (kPointWeight, tetrahedra,
  |     |     eigenvalues, fullKPToIBZKPMap for LAT)
  +-- bond.F90
  |     +-- O_Populate (electronPopulation or
  |     |     electronPopulation_LAT)
  |     +-- O_KPoints (fullKPToIBZKPMap,
  |     |     fullKPToIBZOpMap, atomPerm,
  |     |     numFullMeshKP for star distribution)
  |     +-- O_PSCFIntg (eigenvectors via HDF5)
```

---

## 4. Build System

- Fortran 90; CMake out-of-source build in `build/`
- Compiler set via `$FC` (must match HDF5-compiled compiler)
- Install prefix via `$IMAGO_DIR`
- HDF5 required for primary I/O in imago
- Supported compilers: gfortran, ifort

### 4.1 Debug instrumentation and build flavors

The top-level `CMakeLists.txt` carries opt-in instrumentation
options -- `IMAGO_CHECKS`, `IMAGO_FPE_TRAP`, `IMAGO_INIT_SNAN`,
`IMAGO_WARN_EXTRA`, and `IMAGO_SANITIZE` -- all default OFF.  With
every option off the compile and link flags are byte-for-byte
identical to the historical build, so the production default is
never disturbed; the instrumentation is only ever appended.
`CMakePresets.json` bundles the common combinations as named
presets (`gfortran-release` / `-debug` / `-audit` / `-asan`), each
building into its own tree and installing to a per-preset sandbox,
so a preset can never overwrite the production install.

A *build flavor* is an alternative build of the engine kept under
`$IMAGO_DIR/envs/<flavor>/bin` -- e.g. `asan` for leak hunting,
`audit` for checks + traps + warnings, and later profiling and
parallel-debug builds.  Because only the compiled Fortran
executables differ between flavors (the Python helper scripts and
the `share/` database are compiler-agnostic and identical), a
flavor bin is just symlinks to the production bin with the engine
executables overlaid.  The `envs.sh` switcher (installed to `bin/`,
sourced into the shell) provides `imago_env <flavor>`, which
activates a flavor by repointing only `IMAGO_BIN` (and
`PATH` / `PYTHONPATH`) at it -- `share/` and `.imago/` are reused,
so no data is duplicated and the production install is untouched.
This is the durable replacement for ad-hoc overlay directories.
Activation is per-shell, so different terminals can run different
flavors at once without interference.  Operational details live in
`BUILD.md`; the bug-squashing campaign that motivates the
instrumented flavors is in `dev/DEBUG.md`.

---

## 5. Key Existing Infrastructure

These elements exist in the code and are leveraged by
the current design:

- `kPointStyleCode`: integer controlling how k-points
  are specified (kpoints.f90)
  - 0 = explicit list (legacy; not produced by
    makeinput; Imago emits a warning that
    decomposition properties may be incorrect)
  - 1 = axial counts + shift (mesh built internally;
    primary mode for `-kp` / `-scfkp` / `-pscfkp`)
  - 2 = minimum density + shift (mesh built internally
    via `computeAxialKPoints`; used by `-kpd` etc.)
- `kPointIntgCode`: integer read from input
  (0=Gaussian, 1=LAT) (kpoints.f90)
- `readKPoints`: reads the k-point input file; branches
  on `kPointStyleCode`. For style codes 1 and 2, reads
  point group operations into `convAbcPointOps` and
  `convAbcFracTrans` (conventional-cell-abc fractional
  form, written straight from `share/spaceDB/<sg>`),
  plus the `CONV_LATTICE` block (the conventional cell
  in Bohr) and the `CELL_MODE` flag (`full` or `prim`).
  See DESIGN 2.7 for the basis-invariance design.
  (kpoints.f90)
- `computeAxialKPoints`: converts `minKPointDensity`
  into `numAxialKPoints(:)` using `recipMag` and
  `recipCellVolume` (kpoints.f90)
- `initializeKPointMesh`: builds the uniform mesh from
  `numAxialKPoints` and `kPointShift`; when
  `applySymmetry=1`, folds the mesh to the IBZ using
  `abcRecipPointOps` and saves `fullKPToIBZKPMap`
  (kpoints.f90)
- `computeRecipPointOps`: conjugates the conv-abc point
  group operations in `convAbcPointOps` into the basis
  of whichever reciprocal lattice O_Lattice currently
  holds, using the composed change-of-basis matrix
  `C = invRealVectors^T * M_conv` (with `M_conv` read
  from the kp file's `CONV_LATTICE` block).  In `full`
  mode the `CELL_MODE` flag triggers an identity
  shortcut so the loop becomes a copy.  Output is
  `abcRecipPointOps` (kpoints.f90).  See DESIGN 2.7.
- `computeRealPointOps`: real-space sibling of
  `computeRecipPointOps`. Conjugates `convAbcPointOps`
  and `convAbcFracTrans` into the basis of the loaded
  real lattice using the same `C`, producing
  `abcRealPointOps` and `abcRealFracTrans` for
  `buildAtomPerm` (kpoints.f90, DESIGN 2.7).  Same
  full-mode identity shortcut as the reciprocal-space
  routine.
- `convAbcPointOps`, `convAbcFracTrans`: on-disk
  conventional-cell-abc fractional form of the symmetry
  operations read from the kp file (kpoints.f90).
  Values match `share/spaceDB/<sg>` entry-for-entry;
  the consumer-side conjugation rebases them into the
  loaded cell at runtime.
- `convLattice`, `cellMode`: conventional-cell matrix
  (Bohr) and `full`/`prim` flag read from the kp file's
  `CONV_LATTICE` and `CELL_MODE` blocks.  Together they
  give the consumer the information needed to form the
  change-of-basis matrix `C` and to decide between the
  identity shortcut and full conjugation (kpoints.f90,
  DESIGN 2.7).
- `abcRealPointOps`, `abcRealFracTrans`: symmetry
  operations expressed in the basis of the real lattice
  currently in O_Lattice (full conventional cell or
  primitive reduction).  Consumed by `buildAtomPerm`
  (kpoints.f90 / atomicSites.f90, DESIGN 2.7)
- `sc.full_cell_real_lattice` (structure_control.py):
  snapshot of the conventional lattice captured before
  `apply_space_group()` may overwrite `sc.real_lattice`
  with a primitive reduction; sibling of
  `full_cell_mag` / `full_cell_angle`; same 1-indexed
  layout as `real_lattice`; written into the kp file's
  `CONV_LATTICE` block by the kp-file writer in
  `makeinput.py` (DESIGN 2.7).
- Note: an earlier iteration of this design used a
  producer-side `_to_cartesian_ops` helper to convert
  spaceDB operations to a Cartesian xyz intermediate
  on disk.  Under the current design the producer is a
  near-passthrough: spaceDB operations are written
  unchanged; the only added responsibilities are
  emitting `CONV_LATTICE` from `sc.full_cell_real_-
  lattice` and `CELL_MODE` from the skeleton flag.
  See DESIGN 2.7 for the diagnostic history.
- `generateTetrahedra`: tiles the full uniform mesh
  with tetrahedra (6 per sub-cube) using
  `getIndexFromIndices` for periodic wrapping.
  Called from `initializeKPoints` when
  `kPointIntgCode == 1` (kpoints.f90)
- `getIndexFromIndices(a,b,c)`: converts mesh indices
  to the linear k-point index (kpoints.f90)
- `energyEigenValues(n, i, spin)`: fully in memory
  after SCF
- `recipCellVolume`, `invRealVectors`: available from
  O_Lattice after `initializeLattice`
- `atomSites(:)%cartPos`, `atomSites(:)%atomTypeAssn`:
  atom positions and type assignments from
  O_AtomicSites, available after input parsing

---

## 6. Compute Architecture Direction

### 6.1 Configurable Precision

The program will support compilation in both double
precision (64-bit) and single precision (32-bit) via a
compile-time kind parameter (e.g., `wp` for "working
precision") defined in `kinds.f90`. All floating-point
declarations, literal constants, and MPI/HDF5 type tags
must use this parameter so that a single preprocessor
flag or CMake option switches the entire build.

**Motivation.** Single precision doubles SIMD throughput
on CPU (8 floats vs. 4 doubles in AVX-256) and is
essential for GPU performance, where double-precision
throughput is 2x slower on data-center GPUs and up to
32x slower on consumer hardware.

**Numerical stability.** Single precision provides only
~7 significant digits, which raises concerns in
accumulation loops where small contributions are added
to large sums. The lattice loop in the integral engine
(`gaussOverlapOL` and siblings in `integrals.F90`) sums
contributions from nearest to farthest replicated cells;
far-away terms are individually small but collectively
significant. At double precision this is safe (~15
digits of headroom). At single precision, compensated
(Kahan) summation or smallest-first accumulation order
may be required for specific accumulation sites.
Stability analysis must accompany the precision switch,
targeting at minimum the lattice-sum accumulators and
the eigenvalue solver interface.

### 6.2 Inner-Loop Vectorization

The alpha-pair iteration inside the lattice loop of
`gaussOverlapOL` (and the analogous nuclear and
three-center routines) is currently control-flow-
dominated: a `do while (.true.)` loop with mode
switching, early exits, and per-pair negligibility
tests. This structure prevents SIMD vectorization.

**Target restructure.** Separate the selection phase
(which alpha pairs survive the `alphaDist` threshold for
a given atom pair and lattice distance) from the compute
phase (evaluate the Gaussian overlap integral). The
selection phase produces a packed list of surviving alpha
pair indices; the compute phase processes that list in a
tight, branchless loop amenable to SIMD.

This "gather, filter, compute in bulk" pattern is the
same restructure needed for GPU offload, so the two
goals reinforce each other.

### 6.3 GPU Offload

Long-term, the integral engine and eigenvalue solver are
candidates for GPU acceleration. The restructured inner
loops from 6.2 translate almost directly to GPU kernels:
the packed alpha-pair list becomes a kernel launch over
surviving pairs. Single precision from 6.1 is required
to achieve full GPU throughput.

**Practical path:**
1. Introduce `wp` kind parameter (6.1) and validate
   numerical stability at single precision.
2. Restructure alpha-pair loops into gather/compute
   phases (6.2); verify SIMD vectorization on CPU.
3. Offload the compute phase to GPU (OpenACC, CUDA
   Fortran, or OpenMP target) as a later step.

Each phase is independently useful: phase 1 halves
memory footprint and improves cache behavior; phase 2
speeds up CPU execution; phase 3 adds GPU capability.

### 6.4 Reference: Prior Vectorized Integrals

An earlier vectorized implementation exists outside this
repository at:

```
/home/rulisp/lewis/CPG/cpg-repo/v34/src/olcao/
  integralsSCF.vec.F90    (SCF integrals)
  integralsPSCF.vec.F90   (post-SCF integrals)
  gaussIntegrals.vec.f90  (vectorized primitives)
```

The vectorized `gaussIntegrals.vec.f90` primitives are
also tracked in this repo under `src/imago/`
`src/imago/`.

These files demonstrate the gather/compute separation
described in 6.2. The alpha-pair selection loop
(branchy, scalar) collects surviving pairs into
`lmAlphaPairs` organized by angular momentum type. A
packing step linearizes them into `orderedAlphaPairs`
with `segIndices` marking contiguous same-type regions.
A single batched call (`overlap2CIntgVec`) then
processes all pairs, enabling SIMD over the contiguous
segments. This pattern is the starting point for the
restructure in A4.

### 6.5 Distributed-Memory Parallelism (MPI)

The directions above (6.1-6.4) attack *intra-node*
performance. This section addresses *intra-problem*
scaling across nodes: distributing one large secular
problem over many MPI ranks so a system too large for a
single node's memory still runs. This is the first of the
two parallel axes named in VISION Goal 7; the second,
inter-problem throughput, is the flight layer of section
9.

A sibling branch of the common-ancestor OLCAO code (see
6.7) already located where MPI must reach into the
algorithm. That seam map carries directly to imago:

- **Real-space grid work** -- the site and grid loops in
  `elecStat` (electrostatic potential) and `exchCorr`
  (exchange-correlation). These are embarrassingly
  parallel: a one-dimensional load balance hands each
  rank a disjoint range of sites, and an `MPI_REDUCE`
  accumulates the partials. This is the cheap, proven
  win and the right first target.
- **Integral assembly** -- the atom-pair loops in
  `integrals.F90`. The interim, correct-but-non-scaling
  pattern is replicate-and-broadcast: rank 0 computes,
  then `MPI_BCAST` makes every rank consistent. This
  provides *no* speedup and must not be mistaken for the
  finished state. Genuine distribution assigns atom pairs
  to ranks under the block-cyclic layout (DESIGN section
  9).
- **The secular solve** -- the eigenvalue problem
  `H c = e S c`. This is the actual scaling wall and the
  subject of 6.6.

The MPI lifecycle, a one-dimensional load balancer, and
parallel-HDF5 helpers belong in a dedicated module
(`mpi.f90` in the sibling branch). imago does not yet
carry this module; introducing it is the first build-
level step, paired with switching the gamma and k-point
targets from the serial data reader (`readDataSerial`) to
the parallel one.

### 6.6 Eigensolver Backend Abstraction

The secular solve must sit behind one interface with
swappable backends, chosen by problem size and available
hardware:

- **Serial LAPACK** (`ZHEGV` / `DSYGV`) -- today's path,
  correct for systems that fit one node.
- **Distributed (ScaLAPACK / ELPA)** -- for problems too
  large for one node. Both consume the same block-cyclic
  BLACS data layout, so the distribution machinery is
  shared. ELPA is the leading candidate: its two-stage
  solver outperforms ScaLAPACK on the dense generalized
  problem typical of an LCAO method, and a single ELPA
  API spans CPU and GPU.
- **GPU** -- via ELPA's GPU backend, or a vendor solver
  (cuSOLVER / MAGMA) for the single-node case.

This is the single most important architectural boundary
for parallelization. The sibling branch *declared* a
ScaLAPACK generalized-eigensolver interface (`PZHEGVX`)
but never called it; imago should define the boundary
cleanly and can then complete the distributed arm
(likely with ELPA) rather than inherit an unfinished
stub. Because the abstraction names CPU and GPU backends
behind one call site, it is also where VISION's "device
placement is per-kernel" principle (Principle 14) is
realized for the most expensive kernel in the program.

### 6.7 Reference: upolcao MPI Exploration

A sibling branch of the common-ancestor OLCAO code lives
outside this repository at:

```
~/olcao/src/upolcao/
```

It is *not* a parent of imago and *not* a merge source.
upolcao branched from an older serial OLCAO to add MPI,
while imago rebranded and further developed that same
ancestor on a separate path. The two have since diverged
file by file, so upolcao is mined for design, not
cherry-picked.

What is real and reusable:
- `mpi.f90` -- MPI lifecycle (`initMPI` / `closeMPI` /
  `stopMPI`), a one-dimensional load balancer
  (`loadBalMPI`), a most-square process-grid helper, and
  a parallel-HDF5 attribute check (`checkAttributeHDF5`).
- Genuine grid parallelism in `elecStat` and `exchCorr`
  (the `loadBalMPI` + `MPI_REDUCE` pattern of 6.5).
- ScaLAPACK interface declarations (`PZHEGVX`, `pzgemm`,
  `pdgemm`, `numroc`, a BLACS descriptor type) that are
  tedious to write correctly and reusable as-is.
- The block-cyclic distribution reasoning in `mpi.f90`'s
  header, including the deliberate decision to recompute
  a few atom pairs redundantly rather than communicate
  partial sub-matrices (transcribed in DESIGN section 9).

What is design-only and must not be mistaken for working
code:
- `createProcessGrid` and `createBlockCyclicDistribution`
  are empty stubs.
- The ScaLAPACK eigensolver interface has *zero* call
  sites; the SCF solve still runs serial LAPACK on a
  broadcast copy of the full matrix, on every rank.
- The integral and secular paths use replicate-and-
  broadcast, which is correct but does not scale.

The honest summary: upolcao finished the grid-parallel
half and left the distributed-linear-algebra half as an
elaborate design with stubs. imago can take the finished
half directly and treat the unfinished half as a vetted
design study rather than a starting implementation.

---

## 7. Python Scripts Refactoring Direction

The Python side of the toolchain is currently anchored
by `src/scripts/structure_control.py`, a single module
holding the `StructureControl` class plus a large set
of related utility routines (9,335 lines, 171 tests).
It was produced by the Perl-to-Python port of
`StructureControl.pm` and inherits that file's
monolithic shape: everything related to an atomic
structure -- coordinate manipulation, fractional-to-
Cartesian conversion, neighbor finding, bond analysis,
PDB/skeleton I/O, symmetry operations, supercell
construction, reaction-template surgery, and more --
lives in one class and one file.

This shape has carried the codebase well enough to
reach feature completeness on the port, but it is
approaching the limits of what a single file should
be asked to hold:

- **Navigation cost.**  A file of 9,000+ lines is
  hard to search, hard to diff meaningfully in code
  review, and hard for students to study as an
  example.  The "Prefer concise, self-documenting
  names" rule in `CLAUDE.md` aims at readability;
  file size is the same problem at a higher scope.
- **Domain bleeding.**  New force-field concerns
  (UFF bond parameters, geometry-derived angle
  clustering, LAMMPS type-ID unification) that
  arise in `condense.py` and `make_reactions.py`
  want a natural home.  Adding them to
  `structure_control.py` would couple its
  atomic-structure responsibility to force-field
  data prep -- a domain mismatch.  Forcing
  everything that *touches* a structure into
  `structure_control.py` turns it into a
  miscellaneous bin.
- **Test concentration.**  171 tests in one suite
  already means slow runs and noisy coverage
  reports; adding more erodes the value of running
  the full suite on small changes.
- **Change risk.**  Every edit to
  `structure_control.py` now touches a file
  imported by almost every script.  Small
  localized changes carry outsized blast radius
  on import-time errors and accidental namespace
  collisions.

**Proposed direction (future work; no timeline).**
Split `structure_control.py` along natural domain
seams, with the `StructureControl` class kept as
the main public entry point but its implementation
distributed across focused modules.  Candidate
seams, in rough order of lowest churn first:

1. **`element_data.py`** already exists; the
   element-lookup code that currently lives in
   `structure_control.py` should migrate there.
2. **`geometry.py`** -- coordinate math,
   fractional/Cartesian conversion, distance and
   angle kernels, minimum-image wrapping.  Pure
   functions, heavily reused, no state.
3. **`symmetry.py`** -- space group operations,
   supercell construction, atom permutation
   tables.
4. **`neighbors.py`** -- neighbor search, bond
   analysis, coordination shells.
5. **`io_skeleton.py`**, **`io_pdb.py`**,
   **`io_lammps.py`** -- format-specific I/O,
   each file small and focused on one external
   schema.
6. **`structure_control.py`** (remaining) -- the
   `StructureControl` class itself, delegating
   concrete operations to the modules above.
   Imports stay stable from external callers'
   perspective; internal implementation becomes
   discoverable.

Force-field concerns (UFF bond parameters, angle
clustering) live in their own separate modules
(`bond_utils.py`, `angle_utils.py`, or a shared
`ff_utils.py`) and never mix into
`structure_control.py`'s domain.

**Criteria for starting the split.**  This is a
refactor, not a bug fix, so it should be scheduled
deliberately rather than bolted onto feature work.
Reasonable triggers:
- The file crosses 12,000 lines or the test suite
  crosses 250 tests.
- Two independent pieces of feature work stall on
  each other due to merge conflicts inside
  `structure_control.py`.
- A new concern (e.g., machine-learning potential
  support, a new I/O format) wants to live near
  the structure code but clearly does not belong
  in the monolith.

Until one of those triggers hits, the conservative
choice is to keep `structure_control.py` as-is and
route new force-field-adjacent logic into small
dedicated modules (the path taken for
`angle_utils.py` in the angle-handling rework).

---

## 8. Initial SCF Potential Database

The third active prong (VISION Goal 3) augments today's
isolated-atom potential database with potentials extracted
from converged Imago runs on curated reference solids. This
section covers on-disk layout, file format, the lookup path
through `makeinput.py`, the build pipeline, and the
validation harness. Algorithmic details (interpolation across
nearby labels, environment-descriptor computation) are
deferred to DESIGN.

### 8.1 Layout

The augmented database lives inside the existing per-element
directories under `share/atomicPDB/`. Each element keeps its
own unique store, into which additional files may be added
over time without disturbing the rest of the tree.

```
share/
  atomicPDB/         Per-element database directory
    au/
      pot1               Existing: nuclear/electronic metadata
      coeff1             Existing: Gaussian coefficients + alphas
      s_gaussian_pot.toml  New: augmented database, multiple
                           labeled potentials per element
    o/
      pot1
      coeff1
      s_gaussian_pot.toml
    ti/
      pot1
      coeff1
      s_gaussian_pot.toml
    ...
```

`pot1` and `coeff1` remain the canonical `atomSCF` output
for the isolated-atom potential. The new per-element file is
a *derived* artifact regenerated from `atomSCF` output plus
converged Imago runs on curated reference solids (8.5). It
contains the same isolated-atom data as one labeled entry,
plus one or more additional labeled entries from the
solid-state runs.

### 8.2 File Format

The augmented database uses TOML.  Files are small (well
under 1 MB per element), human inspection of checked-in
files is a hard requirement, and TOML reads via the Python
stdlib (`tomllib`).  Writes are hand-formatted from inside
the helper module (`initial_potential_db.py`) -- the subset
of TOML emitted is small enough that a focused emitter is
simpler than pulling in a third-party writer dependency.

Architectural invariants the format enforces:

- **`schema_version` on every file.**  Lets the reader
  refuse unknown versions and lets future schema changes
  be additive without breaking older files.  Bumped from
  v1 (Phase 1) to v2 (Phase 2) to carry the per-entry
  default tag and the per-entry fingerprint sub-blocks.
- **Provenance is required, not optional.**  Every
  numerical-potential entry must carry enough information
  to retrace its origin -- source structure, atom site,
  Imago commit, convergence parameters, SCF iteration
  count.  Non-negotiable per VISION Principle 5 and feeds
  the validation harness (8.6).
- **Geometric-alpha layout preserved but explicit.**
  `alpha_min`, `alpha_max`, and `num_gaussians` define an
  implicit geometric series; `alphas[]` is written
  explicitly to allow future non-geometric layouts
  without a format change.
- **Numerical content is the same as the legacy
  `pot1` + `coeff1` pair** (radial Gaussian fit of the
  nuclear and electronic potentials); the TOML file wraps
  it in a metadata envelope.

Full schema (top-level keys, per-entry keys, fingerprint
sub-blocks, provenance, validation rules), the
deterministic emitter contract, and a worked sketch live
in DESIGN §5.2, §5.3, and §5.5.  This section names the
architectural choice; DESIGN owns the bytes.

A future arrays-grow-long-enough-to-hurt-readability
contingency (hybrid TOML metadata + sidecar columnar
block referenced by relative path) is the format-level
escape hatch.  The reader/writer abstraction in §8.7 is
designed so that swap is local.

### 8.3 Data Flow

```
share/atomicPDB/<elem>/{pot1,coeff1}   (atomSCF output;
                                        canonical
                                        isolated-atom data)
            |
            | ingested at regeneration time (8.5)
            v
share/atomicPDB/<elem>/s_gaussian_pot.toml
                                       (augmented database;
                                        "isolated" entry
                                        mirrors atomSCF;
                                        extra labeled
                                        entries from
                                        curated Imago runs)
            |
            | makeinput.py reads, picks a label per atom
            | site, emits the chosen coefficients into the
            | Imago input file in the existing on-the-wire
            | format
            v
Imago input file (unchanged format)
            |
            v
Imago SCF (Fortran consumes input identically regardless of
           database source)
```

The Fortran side does not change. All format awareness lives
in `makeinput.py`. The augmentation is a Python-side feature
end to end.

### 8.4 Labeling and Lookup

A "label" is an arbitrary string key under which one
numerical potential entry lives within an element's
database file.  The format imposes no naming convention
on labels beyond one reserved value: `"isolated"`, which
is the atomSCF-derived baseline and must be present per
DESIGN §5.2 rule 6.  Every other label is curator-chosen
(`"default_solid"`, `"fcc_bulk_metal"`, `"tetrahedral_O"`,
etc.) and conveys meaning only to humans reading the file.

Selection -- the question of *which labeled entry to use
for each atom* -- happens through two cooperating
mechanisms:

- **The `default` tag (DESIGN §5.2 rule 7).**  Exactly
  one entry per file carries `default = true`.  This is
  the entry picked for any atom that no scheme grouped
  by local environment.
- **Fingerprint records (DESIGN §5.2 fingerprint
  sub-blocks).**  Each entry optionally carries one or
  more `[[potential.fingerprint]]` records keyed by
  `(method, sub_spec)`.  When the user requests an
  environment-based grouping (`-reduce`, `-bispec`), the
  matcher protocol (§8.9) computes a per-atom fingerprint,
  buckets atoms into species by descriptor similarity,
  and picks the manifest entry whose recorded fingerprint
  best matches each species.

**Phase 1** ships only the literal-label override
(`-pot LABEL` applied uniformly across the structure)
plus the `default` tag fallback.  No fingerprint records
are required in Phase-1 files; entries that lack them
simply don't participate in fingerprint matching.

**Phase 2** layers the matcher-driven path on top of
Phase 1 without changing the legacy override: `-pot LABEL`
still wins when given.  The matcher protocol is the
extension point for new descriptor families (element-
aware bispectrum, SOAP, future schemes) -- adding one is
a new matcher class in `makeinput.py` plus a new
`[[potential.fingerprint]]` shape, with no schema
rewrite.

Full algorithm (CLI surface, mutual-exclusion rules,
spatial scoping, per-species pick, type inheritance) is
in DESIGN §5.6.  The Phase-3 interpolation question --
what to do when the best fingerprint match exceeds the
matcher's similarity floor -- is parked in DESIGN §5.9.

### 8.5 Curation and Regeneration

The augmented database is a build product, not a
hand-edited artifact (VISION Principle 5).

```
src/scripts/
  build_initial_potentials.py    Curates and regenerates the
                                 augmented per-element files
                                 under share/atomicPDB/<elem>/
```

**The producer is a kaleidoscope client.**  It does not
run reference SCFs itself.  For each curated reference
solid it builds a small *predict-then-verify* flight,
hands that flight to kaleidoscope to dispatch and track,
and then harvests the converged potential from the run
directories.  This is the relationship set out in 9.7:
the producer is a *what-to-compute* script, while
kaleidoscope owns running, caching, and tracking the
batch.

This replaces the earlier shape in which the producer
ran SCFs through its own driver and kept its own
per-solid SCF cache.  That bespoke cache is gone:
kaleidoscope's run-reuse cache (9.6, keyed on the
structure file plus the makeinput options and build
identity) subsumes it, so editing an entry
customization no longer re-triggers SCF -- only the cheap harvest step
re-runs.

Inputs:
- A curation manifest in TOML, schema v2 (full spec in DESIGN
  5.7).  It declares the reference solids, each solid's
  `system_type`, a database-wide `[characterization]` fingerprint
  recipe, and optional per-entry customizations.  The harvest is
  automatic -- one representative per distinct environment (DESIGN
  5.2.2/5.2.3) -- so the manifest does not enumerate atom sites;
  a customization, when present, only annotates an auto-discovered
  environment (its `default` tag, description, label, or a pinning
  `atom_site`).
- The historical guidance dataspace (section 10), consulted
  per reference solid to predict the converged k-point density
  and the width of the verification grid around it.  When the
  dataspace is too sparse to predict, the producer falls back
  to a wide default grid.
- The existing `pot1` / `coeff1` files (for `"isolated"`
  entries).
- An Imago build, used *through kaleidoscope* to run the
  verification flights and, for Fortran-side fingerprint
  matchers, the follow-on `imago.py -loen -scf no` runs per
  Fortran-side characterization fingerprint.

**Structure materialization.**  A reference solid is named
either by a local `structure_path` or by a `cod_id` pinned to
a `cod_revision`.  A thin materialization step turns either
form into one local structure file: a `structure_path` is read
from disk; a `cod_id` is fetched once from the Crystallography
Open Database at its pinned revision and written to a plain
local location.  This is the producer's only network access,
and it is deliberately decoupled from any cache -- its sole
job is to hand kaleidoscope a local structure file for its
run-reuse cache to key on.  Fetch failures are strict (named
error, no silent fallback to another revision).

**Manifest authoring.**  A complete manifest is not
hand-written from scratch.  `cod_fish.py` (9.5) discovers and
pins structures and prints a complete sketch -- a `schema_version`
header plus `[[reference_solid]]` stubs whose `reference_id` is
auto-derived from each CIF's metadata -- so `cod_fish.py pin <ids>
> sketch.toml` writes it in one step; `expand_manifest.py`
reads the sketch and fills in the database-wide
`[characterization]` recipe and the shared method defaults, then
writes the finished manifest (per-structure customizations are
optional -- added interactively or by hand).  The producer reads
only that finished manifest, and `cod_fish.py` never writes one.  The manifest schema -- the
dataclasses, the strict and relaxed readers, and the writer --
lives in the shared leaf library `curation_manifest.py`, imported
by both `expand_manifest.py` (to write a manifest) and the
producer (to read one), so the two cannot drift.  The relaxed
reader (`load_structure_sources`) also backs the producer's
`--materialize-only` pre-flight, which fetches and converts every
reference structure before the run and harvest fields are filled
in.

Outputs:
- Updated augmented database file in each affected
  `share/atomicPDB/<element>/` directory (grown in place;
  see the Incremental property below).
- A guidance contribution.  The same converged grid point that
  feeds a potential entry is also harvested into the historical
  guidance dataspace's staging area (section 10), so every
  reference solid the producer converges sharpens the predictor
  for the next one.
- A run log capturing SCF iteration counts and convergence
  metrics per reference run (input to 8.6).

Properties:
- **Reproducible at supported precision.** The TOML
  emitter is deterministic at the bit level: given a
  fixed in-memory database, it always writes the same
  bytes.  The pipeline's numerical outputs are
  reproducible only to within SCF / Gaussian-fit
  precision -- not bit-identical -- because
  floating-point accumulation order, threading, and
  external library versions can shift the last few
  bits, and development changes to the SCF or fit
  code that improve converged potentials are expected
  to perturb the numbers.  Provenance metadata (e.g.,
  `generated_at` timestamps) refreshes on each run
  and is exempt from any reproducibility guarantee.
- **Incremental, on two levels.** The *database* grows
  in place: the producer loads each existing per-element
  file, appends the environments new to it, and skips
  those already present (DESIGN 5.2.3), so several
  manifests accrete into one database and re-running an
  unchanged manifest changes nothing.  The *compute* is
  incremental too: reusing kaleidoscope's run-reuse cache
  (9.6), a run re-executes SCF only for reference solids
  whose structure, options, or build identity changed,
  and adding harvest declarations to an unchanged solid
  re-runs only the harvest.
- **Scriptable.** Invokable from CI or a developer
  workstation without manual intervention.

### 8.6 Validation Harness

VISION Principle 7 sets the headline metric: improved-
initial runs must average at least 20% fewer iterations
than isolated-initial runs across a representative
benchmark.

```
src/scripts/
  bench_initial_potential.py    Runs the benchmark set
                                under both "isolated" and
                                "default_solid" labels,
                                records iteration counts,
                                emits a comparison report.
```

The benchmark set should partially overlap with the
curation manifest (validating that the database performs
on the solids that produced it) and partially extend
beyond it (guarding against overfitting). Both sets are
checked into the repo as part of the deliverable. The
harness report is the artifact that decides whether the
20% target has been met; it is regenerated whenever the
augmented database is regenerated.

### 8.7 Module and Script Impact

Python:
- `makeinput.py`: gains a database reader (the format
  helper from below), a CLI argument that selects the
  global default label, and the augmented-vs-legacy
  lookup path. Legacy `pot1`/`coeff1` lookup remains for
  backward compatibility but is no longer the default
  once the augmented database is populated.
- `build_initial_potentials.py`: new (8.5).
- `bench_initial_potential.py`: new (8.6).

The database reader/writer is encapsulated in a small
helper module (proposed `src/scripts/initial_potential_db.py`)
so that any future format swap or schema-version bump is a
one-file change rather than a spread of edits across the
three scripts above.

The curation-manifest schema -- its dataclasses, the strict
`load_manifest_v2` and relaxed `load_structure_sources` readers,
and the `format_manifest` writer -- is likewise encapsulated in a
small leaf library, `src/scripts/curation_manifest.py`, imported
by both the producer and the authoring tool; it depends only on
the lower libraries it validates against (`initial_potential_db`,
`guidance_db`).  New scripts: `curation_manifest.py` (the schema
library) and `expand_manifest.py` (the sketch-to-manifest
authoring tool, 8.5).

Fortran: no changes for the Phase-1 chain or the Phase-2
base chain.  Imago consumes the same input file format;
the choice of initial potential is opaque to it.  The
Phase-2 follow-up (element-aware bispectrum, TODO C62)
requires extensions to `O_Input` (a new `bispecByElement`
parameter) and `loen.f90` (per-neighbor-element
accumulation in `computeBispectrumComponent` plus an
extended `fort.21` output format).  Those extensions are
scoped under DESIGN §5.9 and TODO D10 / C62 -- they are
not part of the Phase-2 base deliverable.

Shared modules: no changes. The new files live inside the
existing `share/atomicPDB/<elem>/` directories alongside
`pot1` and `coeff1`; no new top-level directory is
introduced.

### 8.8 Open Architectural Questions

Flagged here so DESIGN can pick them up rather than
litigating them inline:

- **Curation manifest format.** *Resolved (2026-05-11).* TOML,
  schema v1 (Phase 1), bumped to v2 in Phase 2 to carry the
  default tag, the database-wide `[characterization]` fingerprint
  recipe, and optional per-entry customizations.  Full schema,
  validation rules, cache layout, and COD-fetch contract are
  specified in DESIGN 5.7.  TOML was chosen so the manifest
  reader uses the same `tomllib` stdlib machinery as the
  per-element database file.
- **Phase 2 selection / interpolation method.** *Selection
  resolved (2026-05-19; decoupled from grouping 2026-06-17,
  C93); interpolation parked for Phase 3.*  The per-species
  database pick (DESIGN 5.6.5) runs for *every* species,
  independent of how it was grouped -- crystallographic,
  position-based (`-target`/`-block`), or environment-based
  (`-reduce`/`-bispec`) -- and is default-on, with
  `-nofingerprint` the opt-out.  The pick forms one
  order-independent representative per species and matches it
  against the element database: when the user chose an
  environment scheme it matches that family at the user's
  `sub_spec`; otherwise it reads the database's `preferred`
  record per family (5.2.2, the `[characterization]`
  convention).  The interpolation question (what to do when
  the best match exceeds the similarity floor) is parked for
  Phase 3; Phase 2 falls back to the default-tagged entry with
  a warning.  Full algorithm in DESIGN 5.6, parameter mapping
  in 5.10.

### 8.9 Matcher Protocol

The Phase-2 selection algorithm dispatches on a small
abstraction called a **matcher**.  Each matcher knows
exactly one descriptor family (e.g., reduce shell-codes,
bispectrum components) and exposes a uniform interface
so the species pass and the producer can call it without
caring which family it is.

**Location.**  The matcher protocol and its concrete
implementations live in `src/scripts/matchers.py`, a
neutral *library* module (not a CLI script, so it does not
add to script proliferation).  The protocol first lived
inside `makeinput.py`, but moved here once `makegroups.py`
needed `BispecMatcher`: keeping it in `makeinput` would
have forced `makeinput` to import upward from `makegroups`
(a cycle).  Every caller now imports it downward --
`makeinput.py` (`ReduceMatcher` / `ReduceStructureView`
for the reduce grouping), `makegroups.py` (`BispecMatcher`
for the sequential loen flow), and
`build_initial_potentials.py` (the `MATCHERS` registry).
This also stages the eventual split in which `makegroups`
owns all environment-based grouping.

**Protocol surface.**  Each matcher class exposes:

  Member                          Purpose
  -----------------------------------------------------
  name                            Matcher identifier
                                  (the string written
                                  into manifest
                                  `method` fields).
                                  E.g., `"reduce"`,
                                  `"bispectrum"`.
  needs_loen_run (bool)           True for matchers
                                  whose descriptor is
                                  computed by Imago
                                  (loen path), false
                                  for matchers that
                                  compute in Python.
                                  True means the
                                  descriptor is obtained
                                  by the sequential loen
                                  flow that `makegroups`
                                  orchestrates (DESIGN
                                  5.10), not in-process.
  default_similarity_floor (float) Per-matcher default
                                  distance threshold
                                  used by the
                                  fingerprint-match
                                  step (DESIGN 5.6.5).
                                  Atoms whose nearest
                                  manifest fingerprint
                                  is farther than this
                                  fall back to the
                                  default tag.  Users
                                  can override per
                                  scheme on the CLI.
  to_loen_input(sub_spec)         Translates a
                                  `sub_spec` inline
                                  table into the
                                  parameter dict the
                                  LOEN block of
                                  `imago.dat`
                                  expects.  Only
                                  meaningful when
                                  `needs_loen_run` is
                                  true.
  parse_loen_output(path,         Reads the enriched,
    sub_spec)                     self-describing
                                  `fort.21` of a loen run
                                  (DESIGN 5.10.3) into
                                  per-site records --
                                  each carrying the row's
                                  identity (element,
                                  species, type) plus its
                                  fingerprint vector.
                                  Only meaningful when
                                  `needs_loen_run` is
                                  true.
  compute_query(structure,        Python-side matchers
    sub_spec)                     only: computes
                                  per-atom fingerprint
                                  vectors from
                                  `StructureControl`
                                  in-process.  Loen-side
                                  matchers do not
                                  implement it -- their
                                  vectors come from the
                                  sequential loen flow
                                  (DESIGN 5.10), read off
                                  `fort.21` by
                                  `parse_loen_output`.
  distance(vec_a, vec_b)          Symmetric scalar
                                  distance in the
                                  matcher's descriptor
                                  space, used for both
                                  species bucketing
                                  and manifest-entry
                                  selection.
  representative(members)         Reduces a list of
                                  member-atom
                                  fingerprints into one
                                  representative
                                  fingerprint per
                                  species (DESIGN
                                  5.6.5 step 2).
                                  `BispecMatcher`
                                  returns the
                                  element-wise mean;
                                  `ReduceMatcher`
                                  returns the first
                                  member's shell-code;
                                  future matchers may
                                  use a medoid or any
                                  scheme appropriate
                                  to their descriptor
                                  space.  The protocol
                                  does not pin the
                                  choice, only the
                                  shape (members in,
                                  one fingerprint
                                  out).
  build_payload(vector)           Wraps one fingerprint
                                  into the dict stored
                                  as a
                                  `FingerprintRecord`
                                  payload (DESIGN 5.2,
                                  5.4): `BispecMatcher`
                                  emits `{values: ...}`,
                                  `ReduceMatcher` emits
                                  `{shell_code: ...}`
                                  (element-only, so the
                                  descriptor transfers
                                  across structures).
                                  The producer (the
                                  fingerprint harvest)
                                  calls it.
  extract_query_vector(payload)   The inverse: reads a
                                  stored payload back
                                  into the form this
                                  matcher's `distance`
                                  expects, so producer
                                  and consumer agree on
                                  field naming.

**Registry.**  `matchers.py` defines a module-level dict
mapping matcher names to matcher classes (imported by every
caller -- `makeinput.py`, `makegroups.py`, and
`build_initial_potentials.py`):

```python
MATCHERS = {
    "reduce":     ReduceMatcher,
    "bispectrum": BispecMatcher,
}
```

Adding a new descriptor family (e.g., element-aware
bispectrum, SOAP) is a new class plus a new dict
entry; no other code needs to change.  The
`initial_potential_db.load()` validator consults
`MATCHERS.keys()` to enforce per-element-database
rule 9 (unknown `method` is a hard error).

**Concrete Phase-2 matchers.**

- `ReduceMatcher` (`needs_loen_run = false`).  Wraps
  the existing reduce algorithm (DESIGN 5.6.4); the
  current in-place implementation in `group_reduce`
  is refactored to live behind this class so producer
  and consumer reach it through the same surface.
- `BispecMatcher` (`needs_loen_run = true`).  Maps
  `sub_spec = {twoj1, twoj2}` to the LOEN input
  block; parses `fort.21` rows into per-site records
  whose vector has `twoj2 + 1` components (the coupling
  channels `j` in `|j1 - j2| <= j <= j1 + j2`, with
  `twoj1 >= twoj2`).  Element-aware mode is gated by an
  optional `by_element` key in `sub_spec`, currently
  ignored (Phase-2 follow-up; DESIGN 5.9 and
  TODO).

**The pick is a per-species step, decoupled from
grouping.**  The fingerprint-based potential pick is not a
side effect of an environment grouping flag: it runs for
*every* species, however that species was grouped --
crystallographic, position-based, or environment-based
(the C93 decoupling; DESIGN 5.6.4/5.6.5).  It is on by
default and disabled per run by `-nofingerprint`.  Which
matcher and parameters the pick uses split by regime:

- *User grouped with an environment scheme* (reduce in
  makeinput; bispectrum is grouped upstream by
  `makegroups`).  That family and the user's `sub_spec`
  drive the match, reusing the per-atom descriptors the
  grouping pass already computed.  The database never
  overrules the user: a database lacking that `sub_spec` is
  a silent best-effort miss to the default entry, never an
  error.

- *Species are file-dictated* (crystalline or pre-assigned).
  The matcher carries no user `sub_spec`, so the database
  supplies it: each element's database marks one
  `preferred` fingerprint record per family, and the
  preferred `sub_spec` for a family is uniform across the
  whole database (the curation convention, DESIGN 5.7).  The
  pick reads the preferred record -- bispectrum if present,
  else reduce -- computes exactly one query, and matches.
  No search across families or sub_specs, and at most one
  loen run per structure.

This adds one protocol obligation beyond grouping: a
`needs_loen_run` matcher must be able to produce a query
for the pick's file-dictated branch.  It does so through
the same loen seam used for grouping (`to_loen_input` ->
loen run -> `parse_loen_output`), reusing a `makegroups`
run when one exists and otherwise paying one fast loen run.
The bispectrum compute is cheap enough that the consumer
simply pays it; optimizing it away is not warranted unless
it ever becomes a measured problem.

**Why this is architectural, not just design.**  The
matcher abstraction is what isolates the Imago Fortran
side from the manifest schema's growth.  Adding a
descriptor family does not touch `imago.f90` unless
the family needs new loen capability; conversely,
adding loen capability (element-aware bispectrum, a
SOAP path) does not touch the manifest schema or the
selection algorithm.  Listing the protocol here pins
that contract down before code lands.

---

## 9. High-Throughput Calculation Flights

This section specifies the infrastructure of VISION
Goal 4: shared machinery for submitting, tracking, and
harvesting batches of Imago calculations.  Its purpose
is that the scripts which decide *what* to compute --
the initial-potential database build (Section 8.5),
convergence sweeps, validation harnesses, and future
ab-initio molecular dynamics and high-throughput
screening -- do not each reinvent *how* to submit and
watch jobs on a cluster.  The flight dispatcher is named
**kaleidoscope**.

The layers and the four scripts/packages that realize
them:

- `imago.py` -- runs one calculation per invocation;
  both a CLI and a callable Python API (9.2).
- `ase_imago.py` -- an ASE Calculator,
  `ImagoCalculator`, adapting Imago to the
  materials-simulation community (9.3).
- `kaleidoscope/` -- a Parsl-based flight dispatcher that
  drives many calculations in parallel and/or series
  batches (9.4).
- `cod_fish.py` + `cif2skl.py` -- the structure-
  acquisition front-end: pull CIFs from the
  Crystallography Open Database, convert them to
  `imago.skl` (9.5).

Three VISION principles are the load-bearing
constraints.  Principle 8 (decouple adapter from
orchestrator) keeps 9.3 and 9.4 independent.  Principle
9 (domain-specific machinery lives at the adapter
layer; the flight layer is ordinary scientific
Python) keeps kaleidoscope free of materials-specific
coupling.  Principle 10 (complete-and-report at the
flight level) means one failed calculation never
fails the whole flight by default.

### 9.1 Layering and dependency direction

The dependency arrows point in one direction only, so
each layer can be replaced without disturbing the
others:

```
cod_fish.py --> cif2skl.py --> imago.skl
                                   |
   kaleidoscope/  --drives-->  makeinput.py --> inputs
        |                                          |
        +----------------- calls ---------> imago.py (API)
                                                   ^
                                                   |
                                   ase_imago.py (ImagoCalculator)
```

`imago.py` is the foundation: it has no dependency on
ASE, on Parsl, or on the flight layer.  The ASE
adapter and kaleidoscope both sit *above* it and call
*down* into it.  This is what lets the cluster jobs and
kaleidoscope run without ASE installed, and lets a
future orchestrator (Snakemake on top, a different
dispatch backend) slot in without touching `imago.py`
or the adapter.

### 9.2 imago.py: CLI and callable API

Today `imago.py` is a command-line driver: it stages
prepared input files, selects the gamma vs. non-gamma
executable, runs the Fortran binaries, collects and
renames outputs, checkpoints completed work, and
manages a lock file.  `imago.py` on its own performs a
plain SCF.

The refactor exposes that same orchestration as a
callable Python API, with the existing CLI reduced to a
thin wrapper over it.  The name stays `imago.py` (not a
separate `run.py`): the module is new, so there is no
namespace collision to avoid, and one obvious name for
"run an Imago calculation" is better than two.

The API offers two entry granularities, so a caller can
join at whichever level it already has inputs for:

- *Prepared-directory mode* -- given a run directory
  that already holds the Imago inputs (`imago.dat`,
  `structure.dat`, `scfV.dat`, kp files), run it as-is.
  No makeinput call; the caller (or a prior step)
  produced the inputs.
- *Structure-and-options mode* -- given a structure and
  a set of makeinput options, drive `makeinput.py` to
  build the run directory first, then run it.

Both return the same small result object
(success/failure, SCF iteration count, paths to the
output files such as the converged
`gs_scfV-<basis>.dat`).  Input *preparation* still lives
in `makeinput.py`; structure-and-options mode simply
calls it on the caller's behalf.

This API is the single seam every higher layer reaches
through, so its contract is deliberately Imago-native
and dependency-free.

### 9.3 ase_imago.py: the ASE Calculator adapter

`ImagoCalculator` subclasses ASE's `Calculator`.  It
translates between ASE conventions (an `Atoms` object,
a results dict keyed by ASE property names, eV/Angstrom
units) and Imago's native world (an `imago.skl`
structure, makeinput options, Hartree/Bohr), and runs
the calculation by calling the 9.2 API.

It is a *separate* module, not folded into `imago.py`,
for three reasons.  (1) Dependency isolation: `ase` is
a third-party package; keeping it out of `imago.py`
means the callable API, the cluster jobs, and
kaleidoscope do not require ASE to be installed.  (2)
VISION Principles 8 and 9 treat ASE as a swappable
adapter, with domain-specific machinery confined to
this layer.  (3) Contract mismatch: the native API is
Imago-shaped while a `Calculator` must be ASE-shaped;
separating them lets each stay clean.

The module is named `ase_imago.py`, *not* `ase.py`: a
flat module named `ase.py` on the scripts path would
shadow the real `ase` package on import.  It stays a
single flat module: `ImagoCalculator` is one class that
exposes many ASE properties (energy, forces, stress,
charges, dipole, and Imago specialties such as DOS,
bands, and bond order as custom result keys) and plugs
into ASE's optimizer, MD, and analysis ecosystem -- all
of which is *one* calculator, not many entry points.
It would be promoted to an `imago_ase/` package only if
auxiliary adapter modules later accrete (a dedicated
DOS/bands adapter, an `ase.db` writer), the same
trigger-based split policy as Section 7.

Committing to ASE buys the materials-community
interoperability VISION Goal 4 describes -- LAMMPS,
ASE's MD integrators for future AIMD, and acceptance by
ASE-consuming workflow tools.

The `Atoms` -> Imago structure translation splits along
the ASE boundary, so that no ASE dependency leaks into
the core:

- `structure_control.py` gains an **ASE-free factory**
  that builds a `StructureControl` from plain arrays
  (lattice vectors, fractional coordinates, element
  symbols).  This is the genuinely shared, broadly
  reused primitive; it has no knowledge of ASE.
- The small **ASE-specific step** -- reading those
  arrays off an `Atoms` object -- stays here at the
  adapter layer.

This matters because `structure_control.py` is imported
by nearly every script (Section 7); putting ASE-domain
machinery there would risk making the whole toolchain
ASE-dependent the first time someone added an `import
ase` for an `isinstance` check or to *construct* an
`Atoms`.  Keeping the factory ASE-free and the
`Atoms`-reading glue at the adapter layer honors
Principle 9 and keeps `import ase` confined to the
modules that genuinely need it.  `cif2skl.py` (9.5)
reuses the same factory.

### 9.4 kaleidoscope: the flight dispatcher

`kaleidoscope/` is a Parsl-based package that drives a
*set* of calculations.  Given a flight specification
(which structures, with which makeinput options), it
dispatches the per-structure work -- makeinput to build
inputs, then the 9.2 `imago.py` API to run -- across
SLURM via Parsl, handling both embarrassingly parallel
sweeps (thousands of independent SCFs) and tightly
iterative inner loops (adaptive convergence, future
AIMD) under one model.

Kaleidoscope dispatches a *unit of work* through a
pluggable wingbeat seam.  The default wingbeat is the
`imago.py` API directly -- the flight layer is
ordinary scientific Python (Principle 9), and most
flights (the database build, convergence sweeps) need
nothing more.  But a unit may also be dispatched
*through the ASE adapter* (for ASE-MD or ASE relaxation
semantics), or through a future adapter, and a single
flight may blend wingbeats -- some units plain Imago
SCF, others adapter-wrapped -- so that Imago
calculations can be mixed with other ASE-compatible
calculations and dispatch activities under one flight.
Keeping the wingbeat pluggable is what lets new adapters
and new blends slot in without changing kaleidoscope's
dispatch core (Principle 8).

Per Principle 10, kaleidoscope records and surfaces
per-job outcomes (converged, non-converged,
cluster-side loss, post-processing error) but does not
abort the flight on a single failure; the client
script decides whether the aggregate result is
acceptable for its scientific purpose.  It is a package
(not a flat module) because it carries real substance:
dispatch, status tracking, harvest hooks, and workspace
management (9.6).

**Cluster dispatch is configuration, not code.**  The
mapping onto SLURM lives entirely in the flight's Parsl
`Config` (DESIGN 6.2.3): the same dispatch core serves a
laptop, an interactive node, and a batch allocation, and
only the `Config` changes.  Two cluster topologies are
both supported because they suit opposite regimes (VISION
Goals 4 and 7):

- **One shared allocation (pooled).**  Parsl requests one
  (optionally auto-scaled) block of nodes and streams many
  units through its workers.  Best for many small, similar
  units -- the convergence sweeps and the database seed.
- **One scheduler job per unit.**  Each unit is its own
  SLURM submission.  Best for large or heterogeneous units
  -- future MPI/GPU solves of differing size -- where a
  single uniform worker slice would waste resources.

No client hand-writes a SLURM script.  A `Config` is
assembled from three layers: a per-site resource-control
file (queues, account, per-node cores and accelerators,
and the commands a worker runs to bring up the imago
environment), a few per-run choices (topology, partition,
nodes, walltime), and the per-unit resource request that
the resource-and-cost dataspace (section 11) predicts.  A
generator fills the defaults and surfaces only the
decisions that are genuinely the user's.  Right-sizing
*heterogeneous* parallel units -- routing each to a block
sized for it -- needs per-size executors keyed on a cost
hint and is staged with section 11, not built up front.
How those layers are settled -- the rc-file home, the CLI choices,
the right-sizing deferral, and the generator's home -- is recorded
in DESIGN 6.2.11, resolving the questions once gathered in 9.8.

### 9.5 Structure acquisition: cod_fish.py + cif2skl.py

Structure acquisition is a front-end, separate from the
flight dispatcher, in two small CLI tools.  Both lean on
existing services rather than reinventing them: COD's own
HTTP query endpoint for search, ASE for CIF parsing, and
Imago's own `apply_space_group` for symmetry expansion.
Neither tool adds a dependency beyond `urllib` (standard
library) and `ase` (already required to read CIF).

**`cod_fish.py` -- COD acquisition front-end.**  Four
verbs, all thin wrappers over COD endpoints; it keeps no
local database and writes no SQL of its own:

- `get <cod_id> [--revision REV]` fetches one structure's
  CIF via `urllib` and is **strict**: with `--revision`
  it verifies the CIF's `$Revision:$` and refuses on a
  mismatch; a network/COD failure or a missing pinned
  revision errors rather than silently substituting a
  different revision (the reproducible path the producer
  and a manifest `cod_id` use; matches DESIGN 5.7).
- `search --elements E [E ...] [--text AUTHOR]
  [--max-extra K]` queries COD's `result.php` and prints
  a **numbered** candidate table (cod_id, formula /
  mineral / common name, spacegroup + number, the six
  cell parameters, volume).  Exact composition is the
  default: the query sets the distinct-element count
  bounds to the length of the element list
  (`strictmin = strictmax = len`), so `--elements Si`
  returns single-element Si phases, not every
  Si-containing compound; `--max-extra K` loosens the
  upper bound to admit up to K further elements.  Search
  also writes the result set to a small session file in
  the working directory so later verbs can refer to rows
  by their **index** rather than their eight-digit id.
- `pin <index|id> [...]` resolves the chosen rows'
  revisions (fetching only those few) and prints a complete,
  ready-to-read sketch manifest: a `schema_version` header
  plus one `[[reference_solid]]` stub per structure, each
  with `cod_id`/`cod_revision` filled in and a `reference_id`
  derived from that CIF's own metadata --
  `<formula>_<H-M symbol>_<IT number>_<year>`, e.g.
  `si_fd-3m_227_2010`.  The year both dates the entry and
  separates phases that share a space group (4H and hcp
  silicon are both P 6_3/m m c); a residual name clash gets
  a trailing counter so every id stays unique (manifest rule
  5).  Each stub also carries two discovery hints read from
  the same CIF -- the composition (`elements`) and a
  `source_description` (chemical name + space group + year) --
  which the authoring tool uses to auto-fill each
  customization's element and description so the curator
  invents neither.  The two part ways in the finished
  manifest: `elements` is a transient, non-schema hint the
  producer ignores and the finished manifest omits, while
  `source_description` is *persisted* as a `reference_solid`
  field (DESIGN 5.7), the structure-level description the
  harvest qualifies per environment.  So `cod_fish.py pin
  <ids> > sketch.toml` writes a sketch the authoring tool reads
  directly.  This is the bridge from browsing to a reproducible
  pull; indices come from the saved search, so a student never
  retypes long ids, and the auto-named id makes a wrong pick
  obvious.
- `rank` (or `search --rank`) is **advisory triage**, not
  validation: when a composition has many COD entries it
  annotates and orders them by interpretable signals
  drawn only from the returned metadata -- ambient vs
  high-P/T conditions, presence of a mineral/common name,
  spacegroup consensus across the candidates, and
  cell-volume-per-formula-unit outliers -- and prints the
  reasons so the curator decides.  It narrows "which of
  these is the real phase," it does not answer it.

Discovery (`search`/`rank`/the index session) is kept
strictly separate from acquisition (`get`/`pin` by pinned
id+revision): a manifest only ever pins an id and a
revision, so reproducibility never depends on a search.

**`cif2skl.py` -- CIF to `imago.skl`, preserving the
space group.**  An `imago.skl` is *not* a flat P1 list of
atoms: it stores the asymmetric unit plus a space-group
token, and `read_imago_skl`'s `apply_space_group`
regenerates the full cell.  Preserving the space group
matters for correctness, not just compactness -- the
Brillouin-zone integration samples the irreducible wedge
*using that space group*, so flattening a crystal to P1
would silently corrupt every k-point-dependent result
(the same hazard `makegroups.py` guards against).  So the
converter recovers the CIF's space group rather than
discarding it:

1. ASE parses the CIF (`ase.io.cif`), which yields the
   authored asymmetric unit (the raw `_atom_site_*`
   list), the cell, the IT number and origin setting, and
   -- read the ordinary way -- the fully symmetry-expanded
   cell.  This is the only step that touches `ase`, and it
   spares us a hand-written CIF parser (CIF's `loop_`
   blocks and symmetry expansion are exactly what ASE
   already handles).  Genuine partial/mixed occupancy is
   refused for now; a refined occupancy within a small
   tolerance of 1 (e.g. 0.9999, a rounding on a physically
   full site) counts as full.
2. The space group is resolved by **verification, not by
   parsing operation tables**.  The `spaceDB` entries for
   the CIF's IT number are the candidate settings (the
   `<IT#>_a`, `<IT#>_b`, ... variants, ordered by ASE's
   reported origin setting so the likely one is tried
   first).  For each candidate the authored asymmetric
   unit is written with that `space` token and run through
   `apply_space_group`; the variant whose regenerated cell
   matches ASE's full expansion (atom count and positions,
   up to lattice translation and permutation, within a
   tolerance) is the answer.  This reuses Imago's own
   symmetry engine as the judge, so no `spaceDB` operation
   set is ever parsed or canonicalized, and the origin
   choice falls out of the match (translations are
   origin-dependent).  It is the automated form of the
   established manual practice -- default to `_a`, check
   the resulting structure, try the others if it is wrong.
3. On a unique match the skeleton is written as the
   asymmetric unit plus that token.  No candidate
   verifying is a **hard error** that names the IT number
   and the tried variants; a `--space <token>` override
   forces a specific setting (still verified).  P1
   fallback is deliberately *not* offered for a crystal,
   since it would reintroduce the corruption above.

No symmetry library (e.g. `spglib`) is added: the CIF's
symmetry is taken as authored (ASE reports it) and
checked against Imago's expansion, rather than
re-derived.  The ASE-free `StructureControl` factory
(9.3) remains the shared primitive for genuinely P1
inputs and the ASE adapter; `cif2skl`'s symmetric path
writes the asymmetric unit and token directly and
verifies through `apply_space_group`.

Splitting acquisition from kaleidoscope keeps the flight
layer agnostic about *where* a structure came from:
kaleidoscope consumes `imago.skl` files, whether they
came from COD via this front-end, from a hand-authored
`structure_path`, or from any other source.

### 9.6 Organizational layout (flight workspace)

A flight that touches hundreds or thousands of
structures needs a directory and naming scheme so the
inputs and outputs do not become an unnavigable mess.
A flight owns a workspace rooted at a single
directory, keyed throughout by a stable per-structure
id (e.g., a COD id or a curation `reference_id`).
The layout is pinned in DESIGN 6.2.4 (the id
charset/uniqueness rule, the `<calc>` tag format, and
the `status.toml` schema); the shape is:

```
<flight_root>/
  flight.toml            What to run + global options.
  structures/<id>/       Acquired inputs.
      <id>.cif
      <id>.skl
  wingbeats/<id>[/<calc>]/     One working dir per calculation:
      <makeinput inputs: imago.dat, structure.dat,
       scfV.dat, kp-*>
      <run outputs: gs_scfV-<basis>.dat, imago.out>
      cache_key.toml      Cache identity snapshot.
      result.toml         Wingbeat-persisted native result.
      status.toml         queued / running / done / failed /
                          lost, plus detail and timings.
  results/               Harvested / aggregated outputs.
  logs/
```

A single structure may host more than one calculation
(e.g., different bases or property runs), hence the
optional `<calc>` level under `wingbeats/<id>/`.

This layout also subsumes the producer's content-keyed
SCF cache (DESIGN 5.7's
`share/atomicBDB/cache/scf/<reference_id>/`): a
kaleidoscope `wingbeats/<id>/` directory *is* a cached run,
so cache-hit/miss logic becomes "is there a completed
run directory for this id whose inputs still match?"

The cache is a **general kaleidoscope feature**, split
into mechanism and policy so that generality does not
cost correctness:

- *Mechanism (kaleidoscope).*  Write a key snapshot into
  the run directory, compare it on the hit-test, skip a
  run whose snapshot still matches, and resume a
  flight over already-completed runs.
- *Policy (client).*  The client supplies the *key
  fields* -- only it knows which inputs define identity
  for its calculations.  The database producer, for
  instance, declares its key as `kpoint_spec` +
  `scf_threshold` + `imago_commit` + the
  structure bytes (DESIGN 5.7).

This keeps kaleidoscope from guessing input identity (a
too-broad key risks false hits and wrong results; a
too-narrow key risks needless re-runs) while still
giving every flight one cache implementation.  The
boundary with `imago.py`'s existing checkpointing is
clean: `imago.py` resumes *within* a run directory;
kaleidoscope decides whether to *launch* the run
directory at all.

### 9.7 Clients, and the producer relationship

Kaleidoscope's clients are the *what-to-compute*
scripts: `build_initial_potentials.py` (the
initial-potential producer, C48),
`bench_initial_potential.py` (the validation harness,
C50), and the future AIMD and screening flights.

This reshapes the producer.  As currently written,
DESIGN 5.7 / PSEUDOCODE 11.4 / ARCHITECTURE 8.5 have
`build_initial_potentials.py` run reference SCFs itself,
with its own COD fetch and its own per-solid cache.
Under this section the producer instead becomes a
kaleidoscope *client*: it hands kaleidoscope the curated
structures and the makeinput options, lets kaleidoscope
run and track the batch, and then harvests the converged
potentials from the run directories.  The converged `scfV`
output carries every potential type in the material (a
`NUM_TYPES` header + per-type Gaussian blocks under a
`TOTAL__OR__SPIN_UP` channel); the harvest selects each named
site's type block and takes its coefficients and alphas
together (the site's type number comes from the `datSkl.map`
described below; the producer runs non-spin, so the
`TOTAL__OR__SPIN_UP` channel is the total potential -- spin
handling is deferred; see DESIGN 5.7 / PSEUDOCODE 11.4).
Those three
producer sections will be revised to delegate to
kaleidoscope once this section stabilizes; that
revision is tracked as follow-up work, not performed
here.

The harvest contract reads one companion output beyond the
converged `scfV`.  Because the producer-derived entry label
encodes the OLCAO `(species, type)` of each harvested site
(DESIGN 5.2.1), and those numbers are assigned by makeinput
during input preparation rather than known at manifest
time, the harvest stage needs them back at storage time.
Rather than add a new output, the producer co-opts the
`datSkl.map` file makeinput already writes (the sorted-dat
to skeleton atom-number mapping): it gains three columns
carrying each site's element, `atom_species_id`, and
`atom_type_id`.  The harvester looks up the manifest
entry's `atom_site` (a skeleton-numbering index) and reads
its `(species, type)` straight off the map.  The assigner
thus records its own verdict where the run lives, and the
harvest stage reads it back to mint the storage label
without re-parsing the run's Imago input.  This keeps the
makeinput/harvest interface explicit: the run produces the
numbers, the producer only consumes them.  The
producer-side change is tracked as TODO C87.

One seam this section once left open -- how the producer
*supplies* the dispatch `Config` (9.4) -- is now closed by
C100.  The producer no longer forces a local executor: it
turns its dispatch choice into a flight `Config` through
kaleidoscope's shared generator and lets the driver
auto-select Local versus Parsl (9.4; DESIGN 6.2.11), so the
seed and database builds reach the scheduler rather than
running one calculation at a time on a login node.
Scheduler dispatch is the default -- a run with no site rc
file present is a configuration error, not a quiet local
fall-back -- and an in-process local run is the deliberate
opt-out that tests, a laptop, and the materialize pre-flight
request explicitly (DESIGN 6.2.11, decision 2).  The wiring is
settled in DESIGN 6.2.11 -- the tiered site rc file, the
per-run CLI choices, the `Config` generator's home in
kaleidoscope, and the uniform-slice deferral of right-sizing --
and the code is C100 (done).

### 9.8 Open architectural questions

Most of the early questions are resolved and recorded
in the subsections above: the callable API offers both
a prepared-directory and a structure-and-options entry
(9.2); the run-reuse cache is a general kaleidoscope
*mechanism* with client-supplied *key fields* (9.6);
the ASE adapter stays a flat `ase_imago.py` (9.3); the
`Atoms`-to-`StructureControl` translation splits into
an ASE-free factory in `structure_control.py` plus
adapter-layer glue (9.3, 9.5); and kaleidoscope
dispatches through a pluggable wingbeat seam, defaulting
to `imago.py` but able to use the ASE adapter or future
adapters (9.4).

The **workspace scheme (9.6)** was deferred here and is
now resolved: the stable-id convention, the `<calc>`
tag format, and the `status.toml` schema are pinned in
DESIGN 6.2.4 (with each `wingbeats/<id>[/<calc>]/` cache
directory carrying `cache_key.toml`, `result.toml`, and
`status.toml`).

What remains open:

- **Producer-section revisions (9.7) -- RESOLVED (C69).**
  DESIGN 5.7, PSEUDOCODE 11.4, and ARCHITECTURE 8.5 have
  been rewritten so the producer delegates SCF to
  kaleidoscope: it materializes each reference structure
  (a thin local read, or a one-shot fetch of the pinned
  COD revision), builds a predict-then-verify flight,
  dispatches one flat batch, harvests each solid's
  converged grid point, and contributes that point back to
  the guidance dataspace (section 10).  The bespoke
  per-solid SCF cache is gone; kaleidoscope's run-reuse
  cache (9.6) subsumes it.  The matching code is TODO C74.

- **Cluster dispatch configuration (9.4) -- RESOLVED
  (DESIGN 6.2.11).**  The four decisions are settled:
  - **Site-config home.**  A dedicated, tiered `*rc.py`
    resource-control file (pure data, like every other
    `*rc.py`) -- a tiny required core (queues, `worker_init`,
    and account where the cluster demands one) with every
    performance and advanced knob optional and defaulted --
    rather than a section of `imagorc`.  A separate tool,
    `cluster_probe.py`, reads the discoverable tiers off the
    scheduler and node hardware (`sinfo`, `scontrol`,
    `lscpu`) and emits a starter file, leaving the
    non-discoverable required core as blanks to complete.
  - **Per-run choices.**  CLI flags (`--dispatch`,
    `--partition`, `--nodes`, `--walltime`) defaulting from
    the rc file, with an optional resolved-config file
    written beside the run for a reproducible record.  The
    command-line default is `slurm-per-job`: on a cluster the
    producer and seed reach the scheduler with no flags, and
    a run with no settings file is a configuration error
    rather than a quiet local fall-back.  `local` is the
    deliberate opt-out (no settings file, no `Config`), which
    the test suite and laptop sessions request explicitly.
  - **Per-unit right-sizing.**  Deferred: both shapes give
    every unit a uniform slice now; right-sizing waits for
    a parallel imago and the section-11 cost predictor that
    would feed it.
  - **Generator sharing.**  The `Config` generator lives in
    `kaleidoscope`, the dispatcher every flight already
    imports, so the producer and future clients share one
    copy.

  Recorded at VISION Goals 4/6/7 and DESIGN 6.2.11; the
  producer change-over and the generator are TODO C100.

## 10. Historical Guidance Dataspace

The fifth prong (VISION Goal 5) accumulates a curated
**dataspace** of converged calculations and a small
**predictor** trained on it.  Each entry is a datapoint:
the chemistry-and-structure signature of a converged
system, the resulting electronic-structure character (band
gap, spin polarization), and the convergence settings
(initially k-point density) that worked.  The predictor is
a two-stage k-nearest-neighbor regression: first it
predicts the electronic character from chemistry, then it
predicts the k-density from the electronic character.  New
calculations query the predictor, get a predicted operating
point plus an uncertainty, and run a verification grid
whose width is set by the uncertainty.  Each successful
flight appends back into the dataspace, so the predictor
gets better as it accumulates evidence.

This section covers on-disk layout, file format choice,
the data-flow through kaleidoscope's flight builder, the
harvest pipeline, and the module boundaries.
Algorithmic details (feature-vector definition, k-NN
metric, predictor stages, confidence measure, grid
widening) live in DESIGN section 7.

### 10.1 Layout

The dataspace lives under a new top-level directory in
`share/`, parallel to `share/atomicPDB/` (Goal 3) and
`share/atomicBDB/` (curation manifest cache):

```
share/
  historicalGuidanceDB/
    entries/                Canonical entries.  One TOML
                            file per entry.  Two-level
                            partition by system_type:
                            crystalline/, amorphous/,
                            nanostructure/, molecular/.
                            File slug:
                            `<system_type>-<short_sha>.toml`.
    staging/                Auto-harvested entries
                            awaiting curator review.
                            Mirrors `entries/`s system_type
                            partition.  Once promoted,
                            files move from staging/ into
                            entries/; the staging area is
                            not consumed by the predictor.
    SCHEMA_VERSION          Single line containing a bare
                            decimal integer (e.g. `1\n`).
                            Bumped on schema change.
                            Readers refuse files whose
                            `schema_version` mismatches.
    elemental_groups.toml   Element-group classification
                            table (alkali / alkali-earth /
                            halide / chalcogen / ... / H).
                            Loaded by the library at init
                            to compute composition
                            vectors.  Checked-in data,
                            not code, per Principle 11.
```

**Why partition only by system_type, not by gap or
chemistry.**  The predictor operates on a continuous
feature space (composition vector + lattice family +
measured gap and spin-polarization); the only *categorical*
axis is system_type, because the predictor's stage-2
relationship (gap+spin → k-density) is qualitatively
different across system_types (a crystalline insulator's
k-density depends on chemistry; an amorphous insulator's
is set by cell size via the density convention; a
molecular cluster is always Γ-only).  Partitioning by
system_type lets the predictor switch sub-models cleanly;
partitioning further by gap or chemistry would force a
continuous variable into bins prematurely.

**One file per entry, not one big file.**  Easy append,
easy git diff, no atomic-write contortions during a
500-entry parallel seed run, and a directory of 500-5000
small files stays tractable.  The predictor loads all
entries into an in-memory dataspace at init time -- TOML
parsing of a few thousand small files takes well under a
second on commodity hardware.

### 10.2 File Format

TOML, same rationale as ARCHITECTURE 8.2: small files,
human inspection a hard requirement, Python stdlib
`tomllib` for reads, hand-formatted emitter for writes.
The full schema, deterministic emitter contract, and
worked sketch live in DESIGN 7.2 / 7.3 / 7.5.
Architectural invariants:

- **`schema_version` on every file.**  Lets the reader
  refuse unknown versions; lets future schema additions
  remain additive.  Day-1 ships v1.
- **Signature block first.**  Every entry carries an
  `[entry.signature]` block declaring `system_type`,
  the 13-dimensional `composition_vector` keyed by
  element-group, and (for crystalline only) the 6-axis
  one-hot `lattice_family`.  This is the predictor's
  feature input.
- **Measured-quantities block.**  `[entry.measured]`
  carries the quantities harvested from the converged
  calculation: `gap_ev`, `gap_kind`
  (`direct`/`indirect`/`none`), `spin_polarization`,
  `total_magnetization`, and the converged
  `kpoint_density` (the target the predictor learns to
  produce).
- **Context block.**  `[entry.context]` carries
  parameters that influenced the measurement: `basis`
  (`mb`/`fb`/`eb`), `functional`, and
  `scf_threshold`.  Future predictors may use
  these as additional regression features; v1's
  predictor conditions on `basis` and `functional`
  (separate sub-models) and ignores the rest.
- **Provenance is required.**  Every entry carries the
  flight id, the source structure id, the imago
  commit, and a UTC timestamp.  Non-negotiable per
  Principle 11.

### 10.3 Data Flow

```
new calculation (structure, options, system_type)
    |
    | The kaleidoscope flight-builder helper
    | (DESIGN 6.2.8) asks the predictor:
    |    "given this structure, predict the converged
    |     k-density and tell me how confident you are."
    v
guidance_db.predict(dataspace, query_signature,
                    basis, functional, kpoint_integration)
    |   (free function over the loaded Dataspace, DESIGN
    |    7.4; query_signature carries system_type, and
    |    (basis, functional, kpoint_integration) selects
    |    the sub-model.)
    |   1. Switch on system_type:
    |      - amorphous, nanostructure, molecular:
    |        consult the small dedicated sub-model
    |        (typically returns "use the density
    |        convention floor", i.e. ~Gamma-only for
    |        large cells), and return.
    |      - crystalline: run the two-stage regression
    |        below.
    |   2. Stage 1 (chemistry -> electronic character):
    |      k-NN regression in composition+lattice
    |      feature space, predicting (gap, intensive
    |      magnetization = |M|/atom).
    |   3. Stage 2 (electronic character -> k-density):
    |      k-NN regression in (gap, magnetization) space,
    |      predicting kpoint_density.
    |   4. Confidence: variance over the k nearest
    |      neighbors at each stage, combined into one
    |      uncertainty score in [0.0, 1.0].
    |   5. Return PredictionResult(
    |        predicted_kpoint_density, confidence,
    |        neighbor_entry_ids, predicted_gap,
    |        predicted_magnetization,
    |        is_under_trained,
    |      ).
    v
The helper builds the verification grid using the
variance-aware widening function (DESIGN 7.7): narrow
when confidence is high, wide when low, wide-grid
fallback when `is_under_trained` is set.
    |
    | Kaleidoscope dispatches the grid (DESIGN 6.2).
    v
flight harvest hook (10.5)
    |
    | 1. Pick the converged grid point per structure
    |    (smallest density where consecutive grid
    |    points' energy delta < threshold).
    | 2. Read measured gap and total magnetization from
    |    the converged calc's result.toml (the kpd comes
    |    from the calc tag; cell facts from the structure).
    | 3. Build a richly-populated GuidanceEntry and
    |    emit it to staging/<system_type>/.
    v
Curator promotion (manual one-at-a-time, batch-promote
for trusted automated flights, or dry-run preview)
moves staging entries into entries/<system_type>/.
    |
    v
The promoted entry joins the in-memory dataspace on
the next library load and is visible to all future
predict() calls.
```

The predictor runs at flight-construction time -- well
before any Imago run.  The harvest runs at flight-
completion time.  Imago itself is unaware of the
dataspace.  All format awareness lives in the new helper
module (10.6).

### 10.4 Feature Space and the Predictor

The predictor operates on a **feature space** that
combines a chemistry vector with a lattice-family vector
(crystalline only) and a coarse system-type partition:

- **Composition vector** (13-dim): atom-fraction weight
  in each of 13 element-group buckets -- alkali,
  alkali-earth, halide, chalcogen, pnictogen, group-IV,
  group-III, transition-metal (lumped 3d/4d/5d),
  lanthanide, actinide, metalloid, noble-gas, hydrogen.
  The group classification table (`elemental_groups.toml`,
  10.1) is a checked-in data file; Principle 11
  requires the chemistry knowledge to be auditable,
  not buried in code.
- **Lattice family** (6-axis one-hot, crystalline
  only): cubic / hex / tet / ortho / mono / tri.
  Cheap-to-extract from the structure file via
  StructureControl.  Mitigates the gap-vs-polymorph
  smoothness risk (DESIGN 7.6).
- **System type partition** (4-way): crystalline,
  amorphous, nanostructure, molecular.  Hard switch
  for which sub-model runs.

The predictor is **k-nearest-neighbor with inverse-
distance weighting**, in two stages for crystalline:
stage 1 maps composition+lattice → (gap, intensive
magnetization); stage 2 maps (gap, magnetization) →
k-density.  The split
exploits transferability: stage 1 is chemistry-heavy and
needs broad chemistry coverage; stage 2 is physics-heavy
and is roughly material-independent.  Each stage carries
its own confidence (variance of the k neighbors); the
two combine into one uncertainty for the verification
grid.

For non-crystalline system types (amorphous,
nanostructure, molecular), the predictor returns a
canonical result driven by the density convention plus
any system_type-specific corrections; chemistry plays
little role.  The implementation is correspondingly
simpler.

Full algorithm (k value, distance metric, normalization,
confidence formula, weighting under sparse data) is in
DESIGN 7.6.  The exact bootstrap behavior when the
dataspace is too thin for stage 1 to be trustable is in
DESIGN 7.9.

### 10.5 Curation, Regeneration, and Harvest

The dataspace is a build product like the initial-
potential database (8.5), but with one key difference:
every successful flight is a *potential contributor*,
not just a hand-curated reference set.  This makes the
harvest hook central machinery:

```
src/scripts/
  guidance_db.py          Library: read,
                                     compute composition
                                     vectors, run the
                                     predictor, emit
                                     entries.  No
                                     orchestration.
  guidance_harvest.py                Producer-side
                                     helper: given a
                                     finished flight,
                                     examine each
                                     structure's
                                     verification grid,
                                     read measured
                                     quantities from
                                     each result.toml,
                                     write staged
                                     entries.
  guidance_promote.py                Curator helper:
                                     review and promote
                                     staging entries
                                     into the canonical
                                     entries directory.
```

- `guidance_harvest.py` is invoked at flight-completion
  time (either as a post-step the flight driver calls,
  or as a standalone CLI run after the fact).  It reads
  the flight's workspace, identifies each verification
  grid, picks the converged point per structure, reads
  the measured electronic-structure quantities from the
  converged calc's result.toml, and emits one staged
  TOML per structure into the appropriate
  `share/historicalGuidanceDB/staging/<system_type>/`
  subdirectory.
- `guidance_promote.py` is the curator's tool.  Four
  modes: interactive one-at-a-time review (default);
  batch `--auto-promote` for trusted automated flights
  meeting an objective acceptance rule (the converged
  k-density landed in the middle 60% of the verification
  grid AND the top three grid points' total-energy
  variance -- read from the staged entry's `grid_energies`
  array -- is below threshold); `--all` to promote the
  whole staging directory after manual review; `--dry-run`
  preview.
  The auto-promotion rule lets a 500-entry seed flight
  promote ~80% of entries unattended, with the curator
  reviewing only the ~20% outliers.
- Both databases grow incrementally and in place (the
  initial-potential DB does too, DESIGN 5.2.3); the
  difference is what gates the growth.  The
  initial-potential DB grows only from a curated manifest
  of reference solids, whereas the dataspace grows
  monotonically from *any* successful flight, with no
  curation manifest in the loop.  Old entries are not
  deleted on schema bumps; a migration tool
  (`guidance_migrate.py`) rewrites them in place.

The **seed flight** (TODO C75) is what populates the
initial dataspace.  It is a one-time stratified sweep
across element-group pairs and common stoichiometry
patterns (~150-250 calculations covering the chemistry
surface representatively rather than at random) feeding
the auto-promotion rule above.  After the seed lands,
ongoing flights (the C48.3 producer, future
characterization runs, etc.) contribute additional
entries as they finish.

### 10.6 Module and Script Impact

Python:

- `src/scripts/guidance_db.py`: new library.
  TOML reader + validation; element-group classifier
  (loads `elemental_groups.toml`); composition-vector and
  lattice-family computation from a StructureControl;
  the in-memory dataspace; the two-stage k-NN
  predictor; the deterministic hand-formatted emitter.
  Imports only `tomllib`, `math`, and the existing
  `structure_control.py`.  Module-level docstring
  describes its role as the **library** half of the
  library/producer/consumer split (10.5).
- `src/scripts/kaleidoscope/` flight-builder helper
  (DESIGN 6.2.8): consumes a `PredictionResult` and
  builds a `Flight` of `CalcUnit`s laid out per the
  tag convention of DESIGN 6.2.4.  This is the
  option-axis half of the builder split per VISION
  Principle 12; the structure-axis half remains
  domain-specific and lives in `structure_control` /
  acquisition.
- `src/scripts/guidance_harvest.py`: new producer
  helper.
- `src/scripts/guidance_promote.py`: new curator
  helper.
- `src/scripts/guidance_migrate.py`: future schema-
  migration tool.  Not in day-1 scope.

Fortran:

- A small extension to imago.py's harvest path so that
  `gap_ev`, `gap_kind`, and `total_magnetization` are
  surfaced (alongside the total energy already present).
  These quantities are computable from existing SCF
  output (the eigenvalue spectrum); the change exposes
  them in the iteration data so any plain SCF run yields
  them (DESIGN 6.1).  Tracked as TODO C76 (under Phase K).

The library / producer / consumer split mirrors DESIGN 5:
the library knows the format and runs the predictor; the
producer (harvest + promote) writes entries; the
consumers (kaleidoscope builder and any future client)
read entries.

### 10.7 Relationship to Other Prongs

- **DESIGN 5 (initial potential database):** the
  guidance dataspace stores convergence *settings* +
  electronic-structure character; the potential DB
  stores converged *potentials*.  They share the
  library/producer/consumer discipline (Principle 11)
  and the same `share/` shape (per-element or
  per-signature TOML), but they are independent
  artifacts with independent lifetimes and no
  cross-references.  A guidance entry never names a
  potential-DB entry, and vice versa; the two share
  only the curation discipline, not their contents.
  Considered and rejected in DESIGN 7.10 ("Closed by
  decision"): the two serve different audiences and
  update cadences, and entangling their schemas would
  couple their lifetimes unnecessarily.
- **DESIGN 6 (kaleidoscope):** kaleidoscope is the
  dispatch layer that runs the verification grid the
  predictor produces.  The dependency goes one way:
  the kaleidoscope flight-builder helper (6.2.8)
  reads the dataspace; the dataspace does not depend
  on kaleidoscope.  Kaleidoscope still works without
  the dataspace -- callers can construct `CalcUnit`s
  directly without going through the helper.
- **C48.3 producer (initial-potential database
  build):** the **first major consumer**.  Its
  workflow under the new prong: for each curated
  reference solid in the potential-DB manifest, the
  producer queries the guidance predictor (system_type
  = crystalline; structure file already known),
  receives a predicted k-density and confidence, lets
  the kaleidoscope helper build a verification sub-
  grid around the prediction, dispatches the sub-grid,
  and harvests the converged potential from the run
  dir whose grid point converged.  Without the
  dataspace seeded, this would run as a wide-grid
  sweep per solid; with the seed in place, the sub-
  grid shrinks to 3-5 points per solid (predict-then-
  verify acceleration).  The seed flight (C75)
  therefore directly accelerates Goal 3.

### 10.8 Open Architectural Questions

- **File naming uniqueness under parallel harvest.**
  Two flights finishing nearly simultaneously could
  collide on a slug.  Plan: include a short SHA over
  (flight_id, source_structure, generated_at) in the
  filename; the flight+structure pair is unique by
  construction, so collisions only occur if two threads
  in the same flight harvest the same structure at
  the same instant.  Detailed slug algorithm in DESIGN
  7.5.
- **Polytype confusion within the predictor.**  Two
  polymorphs of the same compound at the same composition
  vector but different lattice family currently produce
  two entries with the same chemistry features but
  different lattice features.  The k-NN distance metric
  weights these features; whether the default weights
  separate polytypes cleanly enough is an empirical
  question that won't be settled until the seed flight
  lands.  Open knob: the relative weight of composition
  vs lattice in the stage-1 distance metric.  Reasonable
  default in DESIGN 7.6; calibration after seed.
- **Spin polarization unit and interpretation.**
  Recorded as `total_magnetization` (Bohr magnetons per
  formula unit); but for antiferromagnets the net
  magnetization can be zero while the magnetic ordering
  is non-trivial.  A future schema bump may add a
  separate `local_moment_per_atom` field.  Day-1 we
  document the limitation and proceed.
- **Functional / basis as sub-model dimensions vs
  features.**  The predictor conditions on basis and
  functional by running separate sub-models per
  (basis, functional) combination.  Alternative: treat
  them as additional regression features.  The split
  approach is cleaner (no spurious cross-functional
  interpolation) but proliferates sub-models.  Open;
  default approach in DESIGN 7.6 is sub-models per
  combination, with fallback to the most-similar
  combination when the queried one is empty.

## 11. Resource & Cost Guidance Dataspace

### 11.1 Overview and Relationship to Section 10

VISION Goal 6 adds a second curated dataspace that records,
for every imago run, what it *cost* to compute -- peak
memory, disk footprint, and walltime -- as a function of the
problem size and the parallel configuration the run used.
It is a deliberate sibling of the historical-guidance
dataspace (section 10), not an extension of it: the two
share the library / producer / consumer discipline
(Principle 11), the staging-then-promote curation flow, the
schema-versioning-plus-migration pattern, and the
registry-validated key discipline -- but they are
independent artifacts with independent schemas, predictors,
and lifetimes.

The reason they must stay separate is portability.  A
section-10 entry ("MgO converges at k-density 60") is
hardware-independent: it is equally true on a laptop and a
national supercomputer, so that dataspace accumulates
globally and never goes stale when hardware changes.  A
resource entry ("this run used 42 GB and 3.1 h") is the
opposite -- it is meaningful *only* on the machine that
produced it.  Folding cost fields into section-10 entries
would contaminate the one artifact whose value depends on
being portable.  The resource dataspace is therefore
partitioned by a **hardware fingerprint** (11.2), the way
section 10 partitions by `system_type`.

### 11.2 Layout

The artifact lives under its own root in `share/` (proposed
`share/resourceGuidance/`, final name in DESIGN 8),
mirroring section 10's `entries/` + `staging/` split and
its `SCHEMA_VERSION` marker.  Entries are partitioned one
directory level down by hardware fingerprint rather than by
system_type:

```
<root>/
  SCHEMA_VERSION                bare-integer marker.
  hardware_registry.toml        known fingerprints + the
                                probed attributes that
                                define each (analogous to
                                section 10's elemental_groups).
  entries/<fingerprint>/*.toml  promoted observations.
  staging/<fingerprint>/*.toml  harvested, awaiting a
                                curator.
```

The **hardware fingerprint** is a short stable key (e.g.
`<cpu-model>-<cores-per-node>-<mem-per-node>`, exact recipe
in DESIGN 8) naming the node type a run executed on.  It is
the coarse partition within which cost is comparable.  The
finer parallel-configuration knobs (rank/thread counts,
pinning) and the build configuration (11.3) instead live
inside each observation and feed the predictor as features,
not partitions: hardware is partitioned because cost is never
comparable across machines, but parallel-config and build
choices are kept comparable on purpose -- the whole point of
recording them is to learn which configuration, and which
build, is cheapest.

The **atomic unit is one execution observation** -- a single
run under a single configuration -- never collapsed to a
per-system or per-fingerprint summary (VISION Goal 6).  Two
runs of the same structure under different rank counts are
two observations.  This granularity is what lets the same
artifact answer "how much will *this* configuration need"
(provisioning, near-term) and, later, "which configuration
is most efficient" (optimization) and "how does cost scale
with size" (scaling studies), with no schema change.

### 11.3 The Four Blocks of an Observation

Each observation TOML carries four content blocks plus
provenance (field-level schema in DESIGN 8):

- **Size signature** -- the cost-driving dimensions of the
  problem: atom count, electron count (with the core/valence
  split, since the orthogonalized secular dimension is the
  real driver), basis-function count, wave-function
  representation (1-component Schrodinger vs 4-component
  Dirac -- the 4-spinor structure multiplies the secular
  dimension), k-point count, and spin channels.  Derived
  from the makeinput inputs and the structure; known before
  the run.
- **Execution configuration** -- how the run was launched:
  node count, cores per node, total cores, MPI rank count,
  OpenMP threads per rank, and process/thread binding (core-
  or socket-pinning).  Captured as an **extensible,
  registry-validated key-value table**: a checked-in
  `EXECUTION_KNOB_REGISTRY` enumerates the recognized knobs,
  and a new knob (a GPU count, a NUMA policy) is added by
  extending the registry and bumping the schema version --
  never by silently redefining fields.  Known at dispatch
  time from the flight's Parsl provider / launch spec.
- **Build configuration** -- the toolchain the binary was
  compiled with, recorded in *two layers* so it serves both
  the accumulating artifact and a targeted experiment.  The
  **coarse layer** is a registry-validated set of normalized
  knobs that act as predictor features and the axis the
  artifact groups by: compiler family and major version,
  optimization level and a coarse arch/SIMD tag, and for each
  key library (HDF5, ScaLAPACK, BLAS/LAPACK, MPI) its
  implementation, major version, and cost-relevant variant
  (HDF5 parallel vs serial, BLAS threaded vs sequential).
  Like the execution configuration this is an extensible
  `BUILD_KNOB_REGISTRY`, but deliberately bucketed (an
  optimization *level*, not a flag string; a *major* version,
  not a patch) so a build is a comparable feature, not a
  fragmenting one.  The **fidelity layer** is the complete
  compile string and full library build details, stored
  verbatim as provenance (never a predictor feature, just a
  string -- cheap to keep), so nothing is lost.  A later
  study of one specific flag then has two routes: mine its
  state out of the verbatim string post-hoc, or *promote* the
  flag to a named knob in the registry so it becomes a
  first-class comparable feature for that study -- without
  coarsening the global artifact.  Such a study is itself
  just a flight whose varied axis is the flag (on/off) across
  a limited grid of systems (the DESIGN 6 sweep).  Known at
  compile time (11.4).
- **Measured resources** -- what the run actually used: peak
  resident memory, disk footprint (output plus scratch
  high-water mark), and walltime, with optional per-phase
  timings (setup, SCF iteration, eigensolve,
  post-processing).  Like section 10's `METRIC_REGISTRY`,
  the metric set is itself a registry so new measurements
  append cleanly.  Captured after the run (11.4).

### 11.4 Data Flow and Capture

The producer is again the flight layer; the consumer is the
provisioner.  Capture draws from four sources, joined at
harvest:

1. **Dispatch-time (known before the run):** the size
   signature (from makeinput + structure) and the execution
   configuration (from the kaleidoscope flight's Parsl
   provider and launch spec).  The wingbeat records these
   into the run directory alongside its inputs.
2. **Build manifest (compile time):** the CMake build emits a
   `build_info.toml` -- compiler and full flag string, plus
   the detected versions and variants of HDF5, ScaLAPACK,
   BLAS/LAPACK, and MPI -- installed beside the binary, and
   imago can echo its compiled-in configuration.  The harvest
   reads this for both layers of the build block (11.3): the
   coarse knobs and the verbatim compile string.
3. **Scheduler accounting (after the run):** SLURM `sacct`
   supplies `MaxRSS` (peak memory), disk read/write
   high-water marks, and `Elapsed` (walltime); a
   `/usr/bin/time -v` wrapper is the fallback off-scheduler.
4. **Imago self-report (optional, after the run):** the
   Fortran side may write its own peak allocation and
   per-phase wall timings into `result.toml`, which is more
   precise than scheduler granularity for the eigensolve
   phase specifically.

A dedicated harvest helper (`resource_harvest.py`, parallel
to section 10's harvest) walks the finished flight's
workspace, assembles one observation per run dir from these
sources, and writes it to `staging/<fingerprint>/`.  A
curator promotes with the same staging-then-promote
discipline as section 10.  The two harvests are separate
tools reading the *same* workspace; a single completed run
can feed both a section-10 convergence observation and a
section-11 cost observation.

Unlike section 10, **failed runs are not always discarded**:
an out-of-memory or walltime-exceeded run is positive
evidence that a configuration is insufficient at that size,
which is exactly the signal a provisioner needs.  How much
of that negative data to retain -- and how to mark a
censored "needed at least X" so the predictor treats it
differently from a measured "used exactly X" -- is an open
question (11.8).

### 11.5 Feature Space and the Predictor

Within a hardware fingerprint, cost is a smooth,
physics-grounded function of the size signature and the
execution configuration.  This argues for a
**physics-informed regressor** rather than the pure k-NN of
section 10: memory grows roughly as the square of the
secular dimension and the eigensolve as its cube, so the
model can fit scaling exponents (a power law in the secular
dimension, scaled by the parallel configuration) instead of
interpolating raw neighbors.  A k-NN fallback is reasonable
when a fingerprint is thinly populated.  The predictor
choice is deferred to DESIGN 8 and is explicitly *not*
required for the artifact to begin accumulating data.

The **near-term consumer is provisioning** (VISION Goal 6):
given a chosen parallel configuration for a new run, the
predictor estimates the memory, disk, and walltime the run
will need, and the flight layer turns those into SLURM
resource requests with a safety margin.  This directly
serves Principle 6 (cost discipline at SCF setup) and
replaces today's hand-guessed requests, whose failure modes
are wasted allocation (over-request) or killed jobs
(under-request).

Because every observation stores its full configuration, the
same artifact later supports three forward-looking consumers
with no schema change: **configuration optimization**
(compare predicted cost across candidate configurations to
recommend the most efficient), **build comparison** (compare
cost across compiler flags or library builds on one
fingerprint -- the targeted flag study of 11.3), and
**scaling studies** (fit and report cost-vs-size curves).
These are not near-term, but the granular-observation design
(11.2) is chosen precisely so they need no migration.

### 11.6 Module and Script Impact

- `src/scripts/resource_db.py` (new): the library --
  load / validate / hand-formatted emit, the hardware
  fingerprint plus registries, and the predictor.  Built on
  the same patterns as `guidance_db.py` (section 10 / C70).
- `src/scripts/resource_harvest.py` (new): the harvest
  producer (11.4).
- A provisioning consumer in the flight layer (the
  kaleidoscope flight-builder helper, DESIGN 6.2.8, or a thin
  sibling) that reads the predictor and annotates the Parsl
  provider's resource request.  As with section 10, the
  dependency is one-way: the flight layer reads the
  dataspace; the dataspace does not depend on the flight
  layer.
- `src/scripts/resource_migrate.py` (future): the schema
  migration tool, mirroring the planned `guidance_migrate`.
  Not day-1 scope.
- CMake build-system hook: emit `build_info.toml` (compiler
  and full flag string, detected library versions/variants)
  at configure/install time, so the build block is captured
  without hand-entry (11.4).
- Wingbeat / `imago.py` capture hooks (11.4) to record the
  dispatch-time config and scrape scheduler accounting into
  the run directory.

### 11.7 Relationship to Other Prongs

- **Section 10 (historical guidance):** orthogonal axes of
  the same run.  Section 10 answers "what operating point is
  *accurate*" (portable, chemistry-keyed); section 11
  answers "what will it *cost*" (hardware-keyed).  One
  completed run can produce one observation in each.  They
  share library scaffolding and curation discipline only,
  never contents or cross-references -- the same boundary
  section 10 draws against DESIGN 5 (10.7).
- **DESIGN 6 (kaleidoscope):** the flight layer is the
  producer (it captures config + accounting) and the home of
  the provisioning consumer.  The dataspace does not depend
  on kaleidoscope; flights run fine without it (callers
  request resources manually).
- **DESIGN 5 (initial-potential build) / C48.3:** a heavy
  batch consumer.  Once seeded, the potential-DB build can
  size its many per-solid SCF jobs from predicted cost
  rather than a one-size-fits-all request, cutting both
  queue rejection and wasted allocation across a large
  flight.

### 11.8 Open Architectural Questions

- **Hardware fingerprint granularity.**  Too coarse (just
  CPU model) lumps differently-provisioned partitions
  together; too fine (every BIOS/microcode revision)
  fragments the data so no fingerprint accumulates enough
  observations to predict from.  The right granularity is
  empirical; DESIGN 8 proposes a default recipe and a
  fallback when a fingerprint is under-populated.
- **Reliable peak-memory capture.**  `sacct MaxRSS` is
  sampled and can miss short-lived spikes; `/usr/bin/time`
  is per-process, not per-job; the Fortran self-report is
  precise but only for what imago itself allocates (not MPI
  buffers or libraries).  Which source is authoritative, and
  how to reconcile them, is open.
- **Node sharing and contention noise.**  A run that shares a
  node with other jobs sees inflated, noisy walltime and
  memory pressure.  Whether to record an exclusivity flag and
  weight or filter shared-node observations is open.
- **Censored (failed-run) data.**  As in 11.4, an OOM or
  timeout is a bound, not a point measurement.  Whether and
  how the schema and predictor represent censored
  observations is open.
- **Any portable normalization.**  Walltime is not portable,
  but core-hours or a hardware-normalized cost index might
  transfer partially across fingerprints.  Whether a
  normalized cost is worth recording (to bootstrap a new
  fingerprint from related ones) is open.
- **Which build knobs deserve coarse-layer promotion.**  The
  two-layer build record (11.3) bounds this -- the full
  compile string is always kept, so nothing is lost -- but
  *which* knobs are stable and cost-relevant enough to live
  in the coarse `BUILD_KNOB_REGISTRY` (and thus drive the
  predictor by default) is a judgment call DESIGN 8 must
  seed and the curator revisits as evidence accumulates.
- **Build effects on numerics, not just cost.**  A build
  block primarily conditions cost here, but the same choices
  (fast-math, a different BLAS summation order) can perturb
  low-order digits of the result.  That is a reproducibility
  concern straddling section 10, where today only
  `imago_commit` keys the binary's identity.  Whether the
  build block should also be referenced from the convergence
  side -- against the no-cross-reference boundary (11.7) --
  is an open question to settle in DESIGN, not assume now.
