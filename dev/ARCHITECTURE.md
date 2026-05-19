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
    imago.py          Imago driver (also a callable API; 9.4)
    makeinput.py       Input file orchestrator
    imago/             The `imago` Python package (Section 9)
      run.py           Callable API for one Imago calculation
      ase.py           ASE Calculator (ImagoCalculator)
      campaign/        Parsl-based campaign helpers
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
through `makeinput.py`, the regeneration pipeline, and the
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

Inputs:
- A curation manifest in TOML, schema v2 (full spec in DESIGN
  5.7), listing reference solids, atom sites to harvest from
  each, the labels to assign, which entry per element carries
  the `default` tag, and which `(method, sub_spec)` fingerprints
  to harvest alongside each numerical potential.  Reference
  structures are fetched from the Crystallography Open Database
  at regeneration time using a pinned revision, with a
  `structure_path` escape hatch for materials not in COD.
- The existing `pot1` / `coeff1` files (for `"isolated"`
  entries).
- An Imago build (for running reference SCF calculations and,
  for Fortran-side fingerprint matchers, follow-on
  `imago.py -loen -scf no` runs per declared fingerprint).

Outputs:
- Regenerated augmented database file in each affected
  `share/atomicPDB/<element>/` directory.
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
- **Incremental.** Adding manifest entries updates only
  the affected element files; existing entries are not
  touched unless the manifest requests a rebuild.
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
  default tag and the per-entry fingerprint declarations.
  Full schema, validation rules, cache layout, and COD-fetch
  contract are specified in DESIGN 5.7.  TOML was
  chosen so the manifest reader uses the same `tomllib` stdlib
  machinery as the per-element database file.
- **Phase 2 selection / interpolation method.** *Selection
  resolved (2026-05-19); interpolation parked for Phase 3.*
  Per-site label selection is now driven by the matcher
  protocol (8.9): environment-based grouping flags
  (`-reduce`, `-bispec`) compute per-atom fingerprints,
  bucket atoms into species by fingerprint similarity, and
  pick the manifest entry whose recorded fingerprint best
  matches the species centroid.  Atoms outside any
  environment scheme fall through to the file's
  default-tagged entry.  The interpolation question (what
  to do when the best fingerprint match exceeds the
  similarity floor) is parked for Phase 3; Phase 2 falls
  back to the default tag with a warning.  Full algorithm
  in DESIGN 5.6, parameter mapping in 5.10.

### 8.9 Matcher Protocol

The Phase-2 selection algorithm dispatches on a small
abstraction called a **matcher**.  Each matcher knows
exactly one descriptor family (e.g., reduce shell-codes,
bispectrum components) and exposes a uniform interface
so the species pass and the producer can call it without
caring which family it is.

**Location.**  The matcher protocol and its concrete
implementations live inside `src/scripts/makeinput.py`,
not in a new script.  Per the parked memory note's
decision to avoid script proliferation, the matcher
classes are co-located with the species-pass machinery
that drives them.  Producer-side use of matchers
(`build_initial_potentials.py`) imports the protocol
from `makeinput.py`.

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
                                  Drives the
                                  nested-makeinput
                                  bootstrap of DESIGN
                                  5.10.
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
  parse_loen_output(path,         Reads `fort.21` from
    sub_spec)                     a loen run and
                                  returns per-site
                                  fingerprint vectors.
                                  Only meaningful when
                                  `needs_loen_run` is
                                  true.
  compute_query(structure,        For Python-side
    sub_spec)                     matchers, computes
                                  per-atom fingerprint
                                  vectors from
                                  `StructureControl`.
                                  For loen-side
                                  matchers, this is
                                  the outer entry
                                  point that triggers
                                  the bootstrap of
                                  DESIGN 5.10 and
                                  returns the parsed
                                  vectors.
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

**Registry.**  `makeinput.py` maintains a module-level
dict mapping matcher names to matcher classes:

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
  block; parses `fort.21` rows into vectors of length
  `2 * twoj2 + 1`.  Element-aware mode is gated by an
  optional `by_element` key in `sub_spec`, currently
  ignored (Phase-2 follow-up; DESIGN 5.9 and
  TODO).

**Why this is architectural, not just design.**  The
matcher abstraction is what isolates the Imago Fortran
side from the manifest schema's growth.  Adding a
descriptor family does not touch `imago.f90` unless
the family needs new loen capability; conversely,
adding loen capability (element-aware bispectrum, a
SOAP path) does not touch the manifest schema or the
selection algorithm.  Listing the protocol here pins
that contract down before code lands.
