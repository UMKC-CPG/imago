# Debug

## Purpose

This document is the bug-squashing campaign ledger for the `imago`
electronic-structure engine -- the Fortran application under
`src/imago/` together with the shared modules it links (`kinds`,
`constants`, `elementData`, `readDataSerial`, `writeData`). It
holds two things: the *methodology* for hunting bugs and memory
leaks in that code, and the *findings ledger* of what we discover,
ordered by severity so that we can work through fixes later.

The motivating goal is to squash the bugs in the serial code
*before* a parallelized version is developed. Bugs are far harder
to find and reproduce in a parallel environment, so the serial
code should be made as clean as we can manage first.

This is a *tracking artifact*, not a sixth level of the design
chain (VISION -> ARCHITECTURE -> DESIGN -> PSEUDOCODE -> source).
It sits alongside those documents in `dev/` and references them
where a finding touches a design decision, but it does not itself
specify behavior. When a bug here reveals a flaw at the DESIGN or
ARCHITECTURE level, that flaw should be propagated up the chain in
the normal way, with the DEBUG entry left as a pointer.

## Status

- Date opened: 2026-06-25
- Current phase: Phase 0 (build harness). **Step 0a COMPLETE**
  (2026-06-26): the opt-in gfortran instrumentation options are
  wired into the top-level `CMakeLists.txt` and verified -- the
  default build is byte-for-byte unchanged, and the options inject
  correctly when enabled. Next: 0b (presets), 0d (harness
  validation), 0e (docs). 0c (ifort) deferrable.
- Phase 2 audit mechanism: multi-agent workflow (decided)
- Findings recorded so far: 0

### Phase 0a verification record (2026-06-26)

Confirmed in throwaway build trees (since removed), environment
replicated as `cpg` + venv + `imagorc`, `FC=h5fc`:

- **Default RELEASE**, no options: real/complex both
  `-O3 -fimplicit-none -Wall` (+`-DGAMMA` real) -- matches baseline.
- **Default DEBUG**, no options:
  `-Og -g -fcheck=all -fimplicit-none -Wall -fbacktrace -DDEBUG`
  -- matches baseline.
- **Debug + all options on** (`IMAGO_CHECKS`, `IMAGO_FPE_TRAP`,
  `IMAGO_INIT_SNAN`, `IMAGO_WARN_EXTRA`, `IMAGO_SANITIZE=address`):
  every flag is appended to both variants and `-fsanitize=address`
  reaches the linker. (Two harmless cosmetic effects appear only
  in instrumented builds: a leading space on the flag string, and
  a duplicated `-fcheck=all` / `-fbacktrace` already present in the
  Debug build type. Neither affects the default build or compiler
  behavior.)

## Phase 0a baseline (ground truth, captured 2026-06-26)

The exact flag strings the working build emits today, read from
the already-configured `build/release` (RELEASE) and `build/debug`
(DEBUG) trees with `FC=h5fc` (gfortran 15.2). These are the
reference that step 0a.4 must reproduce byte-for-byte after the
new options are added. The only per-variant difference is the
`-DGAMMA` define contributed by `src/imago/real/CMakeLists.txt`.

```
RELEASE  real (imagoG):    -O3 -fimplicit-none -Wall      + (-DGAMMA)
RELEASE  complex (imago):  -O3 -fimplicit-none -Wall
DEBUG    real (imagoG):    -Og -g -fcheck=all -fimplicit-none -Wall
                           -fbacktrace -DDEBUG             + (-DGAMMA)
DEBUG    complex (imago):  -Og -g -fcheck=all -fimplicit-none -Wall
                           -fbacktrace -DDEBUG
```

## The two code bases (real vs complex)

The single most important structural fact for this campaign is
that one source tree compiles into two distinct executables, and a
bug can live in one while the other is perfectly fine.

The split is driven by a single preprocessor definition, `-DGAMMA`,
added only in `src/imago/real/CMakeLists.txt`:

- `imagoG` -- the **real** (gamma-point) build. `GAMMA` is
  *defined*, so every `#ifndef GAMMA` block is **excluded**.
- `imago` -- the **complex** (multi-k) build. `GAMMA` is
  *undefined*, so every `#ifndef GAMMA` block is **included**.

The files that diverge most heavily between the two builds (by
count of `GAMMA`-related preprocessor lines) are `integrals.F90`,
`secularEqn.F90`, `field.F90`, `optc.F90`, and
`integrals3Terms.F90`. These are where a fix applied to one branch
but not the other, or a real-vs-complex type mismatch, is most
likely to hide.

**Methodological rule:** the manual audit (Phase 2) reads the
*preprocessed* source for each variant (`gfortran -E -DGAMMA` for
real, `gfortran -E` for complex), never the raw source. This
guarantees we only ever reason about code paths the compiler
actually sees, and that any divergence between the two variants is
made explicit rather than missed.

## Severity scale

Each finding is ranked so the ledger can be worked top-down. The
scale weighs *consequence* (corruption / wrong answer / crash /
nuisance) against *reach* (every run / common path / rare path).

- **S1 -- Critical.** Memory corruption, out-of-bounds writes,
  silently wrong scientific results on a common path, or crashes
  on common input. Fix before anything else.
- **S2 -- High.** Leaks that grow with SCF iteration or k-point
  count (these dominate long production runs), wrong results on an
  uncommon-but-real path, or use of uninitialized memory whose
  effect is data-dependent.
- **S3 -- Medium.** Bounded leaks (a fixed amount once per run),
  latent bugs reachable only on rare paths, missing-`else`
  fallthroughs with currently-benign effect.
- **S4 -- Low.** Robustness and hygiene: missing `stat=` on
  allocation, unchecked I/O status, dead code, fragile interfaces
  -- no present misbehavior but a hazard for future change.

A separate, orthogonal tag flags **parallelization hazards**:
mutable `SAVE`/module state, non-reentrant routines, and
race-prone I/O. These may be perfectly correct in the serial code
yet block or endanger the planned parallel version, so we surface
them now even when their serial severity is low.

## Variant tags

Every finding records which executable(s) it affects:

- **[BOTH]** -- present in code compiled into both builds.
- **[GAMMA]** -- only in the real (`imagoG`) build (`GAMMA`
  defined).
- **[COMPLEX]** -- only in the multi-k (`imago`) build (`GAMMA`
  undefined; inside a `#ifndef GAMMA` block).

## Bug categories (audit taxonomy)

A fixed taxonomy keeps the fan-out subagents consistent and the
ledger searchable:

- **LEAK** -- allocatable or pointer never freed. Must be
  scope-classified (see note below), not counted blindly.
- **ALLOC** -- allocation-status faults: double `allocate`,
  `deallocate` of an unallocated entity, use-after-deallocate,
  missing `stat=` on a failure-prone allocation.
- **UNINIT** -- use of a variable before it is set, including
  partially-initialized arrays and derived-type components.
- **BOUNDS** -- array index out of range, shape/rank mismatch,
  off-by-one in loop bounds.
- **PTR** -- dangling pointer, undefined association status, use
  of a pointer after its target is freed, missing `nullify`.
- **NUM** -- numerical hazards: division by zero, `sqrt`/`log` of
  a non-positive argument, NaN/Inf propagation, `kind`/precision
  mismatches in mixed expressions.
- **HDF5** -- HDF5 handle leaks: a dataset, dataspace, group,
  property list, or file opened and not closed. These leak library
  resources even when Fortran memory is clean; the HDF5-heavy
  `hdf5*.F90` files get a dedicated balance check.
- **IFACE** -- interface and `intent` mismatches: implicit vs
  explicit interface disagreement, wrong `intent`, argument
  shape/type mismatch, especially across the `#ifdef` boundary.
- **LOGIC** -- control-flow errors: missing `else`, wrong branch,
  incorrect loop nesting, fallthrough.
- **DIVERGE** -- real-vs-complex divergence: a fix or guard
  present in one `GAMMA` branch but absent in the other.

**Leak scope note (false-positive guard).** In Fortran 95 and
later, a local `allocatable` array is automatically deallocated
when its procedure returns, so an `allocate` with no matching
`deallocate` is *not* automatically a leak. A leak claim must
identify why the entity persists: it is `SAVE` or module-scope, it
is a `pointer`, it is an allocatable function result or
`intent(out)` dummy handed back to a caller, or it is
re-`allocate`d inside a loop without an intervening `deallocate`.
The raw repo-wide count of 162 `allocate` against 73 `deallocate`
is therefore a starting map, not a leak tally.

## Entry schema

Each finding in the ledger uses this fixed shape so that commits,
`TODO.md`, and future sessions can reference a stable ID:

```
### BUG-NNN -- <short title>
- File:     <path>:<line>  (and the procedure/module name)
- Variant:  [BOTH] | [GAMMA] | [COMPLEX]
- Category: LEAK | ALLOC | UNINIT | BOUNDS | PTR | NUM | HDF5 |
            IFACE | LOGIC | DIVERGE
- Severity: S1 | S2 | S3 | S4   (+ PARALLEL-HAZARD if applicable)
- Status:   open | confirmed | fixed | wontfix | duplicate
- Evidence: how we know -- compiler warning, sanitizer/valgrind
            trace, or the reasoning from the manual read.
- Fix:      the suggested remedy (and any design-chain
            propagation it implies).
```

## Campaign plan

The phases run cheapest-and-broadest first, so mechanical findings
are captured before the labor-intensive reading begins.

### Phase 0 -- Build harness

Rework the CMake build so that debugging features toggle
independently and the compiler choice is clean, *without* changing
the current, working default build. This is the prerequisite for
the runtime phase: it is what produces the instrumented binaries
we will run the test decks under.

The current state of the build (for reference):

- Compiler is selected through the `$FC` environment variable,
  which on this machine is `h5fc` -- an HDF5 wrapper that itself
  wraps gfortran 15.2 (conda-forge). The top-level `CMakeLists.txt`
  detects the underlying compiler ID and sets one bundled flag
  string per build type (`RELEASE`, `DEBUG`) for GNU and Intel.
- A key constraint, already documented in `CMakeLists.txt`: imago
  must be built with the *same* compiler that built the HDF5
  library, because the two compilers cannot read each other's
  `.mod` files. So "use ifort" really means "have an `h5fc` that
  wraps ifort," which is a toolchain precondition, not just a flag.

Proposed shape -- granular, default-OFF cache options that compose
onto the existing flags, plus convenience presets:

- `IMAGO_CHECKS`     -> `-fcheck=all` / `-check all`
- `IMAGO_FPE_TRAP`   -> `-ffpe-trap=invalid,zero,overflow
                         -fbacktrace` / `-fpe0 -traceback`
- `IMAGO_INIT_SNAN`  -> `-finit-real=snan -finit-integer=-99999` /
                         `-init=snan,arrays` (flushes use of
                         uninitialized values into an immediate
                         trap rather than a silent wrong answer)
- `IMAGO_SANITIZE`   -> gfortran `-fsanitize=address|undefined|leak`
- `IMAGO_WARN_EXTRA` -> `-Wextra -Wuninitialized -Wrealloc-lhs-all
                         -Warray-temporaries -Wimplicit-interface`

plus a `CMakePresets.json` (cmake 3.26 supports it) giving named
presets for {gfortran, ifort} x {release, debug, audit, asan},
each building both `imagoG` (real) and `imago` (complex).

**Phase 0 stabilization doctrine (the off-ramp).** The conda /
HDF5 / ifort matching problem can become a swamp. To keep it from
corrupting the working build, Phase 0 obeys these invariants and
is tiered so we can stop early at the safe subset:

- *The default build is sacrosanct.* Running `cmake` with no
  options and `FC=h5fc` must reproduce today's flags exactly. All
  new machinery is additive and defaults to OFF. Doing nothing
  must behave identically to today.
- *Separate build trees.* Instrumented builds live in their own
  directories (the existing `debug/`, `release/` pattern), never
  overwriting the known-good tree.
- *Each increment is its own commit*, so any step can be undone
  with a single `git revert` if it destabilizes the build.

The work is tiered with explicit stop points. The full task
breakdown follows; each item names what it does, why, its risk,
and the condition under which it is "done."

**0a -- Baseline capture + zero-risk gfortran toggles**
*(the stabilization target / lock point).* This is the safe
core; needs no new toolchain at all, since gfortran is already
live. If the rest of Phase 0 turns swampy, we lock the campaign
here and run Phases 1-3 gfortran-only.

- *0a.1 Capture the baseline.* Record the exact flag strings
  today's `RELEASE` and `DEBUG` builds emit for gfortran, for both
  `imagoG` and `imago`, so we can later prove "no change."
- *0a.2 Add granular, default-OFF cache options* that *append*
  onto the existing flags -- the `IMAGO_CHECKS`, `IMAGO_FPE_TRAP`,
  `IMAGO_INIT_SNAN`, `IMAGO_SANITIZE`, and `IMAGO_WARN_EXTRA`
  options listed under "Proposed shape" above.
- *0a.3 Wire both variants* so `imagoG` and `imago` inherit the
  options identically.
- *0a.4 Prove the invariant.* `cmake` with no options and
  `FC=h5fc` reproduces the 0a.1 baseline byte-for-byte.
  Done-when: the default build is unchanged, and each option
  flips its flags on and off as expected.

**0b -- gfortran convenience presets** *(pure ergonomics, no
change to compile semantics).*

- `CMakePresets.json` (cmake 3.26 supports it) with named presets
  `gfortran-release`, `gfortran-debug`, `gfortran-audit` (checks +
  SNaN + FPE + extra warnings), and `gfortran-asan` (sanitize =
  address + leak).
- Each preset binds a dedicated build directory (the existing
  `debug/` / `release/` pattern) so trees never overwrite the
  known-good build.

**0c -- ifort enablement** *(deferrable -- the swamp).* Explicitly
optional and gated; can be closed off without blocking 0a, 0b, or
Phases 1-3.

- *0c.1 Locate the toolchain.* Find ifort/ifx (likely a module
  load) *and* an ifort-wrapped `h5fc` -- HDF5 `.mod` files are not
  cross-compiler readable, so the HDF5 build must match. This is
  the precondition, not merely a flag.
- *0c.2 Mirror the toggles in Intel syntax* (`-check all`,
  `-fpe0 -traceback`, `-init=snan,arrays`, `-warn all`).
- *0c.3 Add ifort presets* to `CMakePresets.json` (stubbed until
  0c.1 succeeds), completing the {gfortran, ifort} x {release,
  debug, audit, asan} matrix.

**0d -- Harness validation** *(the acceptance test that proves
Phase 3 is unblocked).*

- Build one instrumented binary (the `audit` or `asan` preset) and
  confirm it *links and runs one small test deck to completion*
  without the harness itself erroring -- the sanitizer runtime
  resolves and the traps do not fire spuriously. Full runtime
  hunting stays in Phase 3; this only proves the harness works.

**0e -- Build-harness documentation.**

- A short note (in `dev/` or a build README) listing each option,
  what it does, and the preset names, so a student can select
  "I want a leak-checking build" without reading the CMake. Per
  the project's documentation policy.

**Critical-path subset.** 0a -> 0d yields instrumented gfortran
binaries -- enough to run all of Phases 1-3. So the must-haves are
**0a, 0d, 0e**; **0b** (convenience) and **0c** (the second
compiler) are nice-to-haves we can defer without blocking the
campaign.

### Phase 1 -- Compiler sweep

Build both variants (and, once 0c lands, both compilers) with
`IMAGO_WARN_EXTRA`, and cross-diff the warning logs. Each compiler
catches diagnostics the other misses (uninitialized use, unused
entities, implicit-interface calls, suspect conversions). The
deduplicated warnings become the first tranche of ledger entries,
essentially for free.

### Phase 2 -- Manual audit (multi-agent workflow)

The human-judgment read that compilers and sanitizers cannot do:
scope-classified leaks, logic errors, real-vs-complex divergence,
`intent`/interface mismatches, HDF5 handle balance, and numeric
hazards. Executed as a multi-agent workflow that fans the file
groups out to parallel subagents against the fixed taxonomy and
schema above, with a synthesis pass to dedup and severity-rank.
Each subagent reads the *preprocessed* source for the variant it
is assigned, per the methodological rule above.

Provisional file grouping (independent subsystems, each one
subagent): input/parsing; lattice/k-points; basis/atomic data;
the integrals family (`integrals*.F90`, `intg*.F90`); the secular
solve (`secularEqn`, `matrixSubs`); SCF charge/potential
(`coreCharge`, `valeCharge`, `potential*`, `elecStat`, `exchCorr`,
`populate`); the HDF5 I/O family; post-processing (`dos`, `bond*`,
`optc*`, `field`, `forces`, `mtop`, `dimo`, `loen`); and the
driver (`imago.F90`, `commandLine`, `interfaces`).

### Phase 3 -- Runtime instrumentation

The highest-yield phase, now unblocked because representative
input decks are available. Run the decks through the Phase-0
instrumented binaries:

- valgrind `--leak-check=full` for real heap leaks and invalid
  reads/writes,
- gfortran `-fsanitize=address` (with leak detection) as a faster
  cross-check,
- `IMAGO_INIT_SNAN` + `IMAGO_FPE_TRAP` to catch use of
  uninitialized reals and NaN/Inf propagation at the point of
  origin.

Run both variants. Findings here outrank static suspicions because
they are observed, not inferred.

### Phase 4 -- Synthesize and rank

Merge the findings from all phases into the ledger below, dedup,
assign final severities and IDs, and order by severity. From there
the ledger drives the fix work and can be cross-linked into
`TODO.md`.

## Decisions log

- **Phase 2 mechanism:** multi-agent workflow (parallel subagent
  fan-out), chosen 2026-06-25.
- **Generated integral files:** `gaussIntegrals.f90` (135K lines)
  and `gaussIntegrals.vec.f90` (39K lines) are machine-generated
  and treated as *trusted*. They receive only a structural
  spot-check (allocate/deallocate balance, interface sanity), not
  a full read. Their *generator* programs are a separate, later
  effort; a hand-bug in the output is really a generator bug.
- **Phase 0 sequencing:** proceed slowly and carefully through the
  build-system work, with the stabilization doctrine above as the
  guaranteed off-ramp so the working compile process is never put
  at risk.

## Findings ledger

No findings recorded yet. Entries will be added here, ordered by
severity (S1 first), as the phases proceed.

<!-- Template -- copy for each new finding:

### BUG-001 -- <short title>
- File:     <path>:<line>  (procedure/module)
- Variant:  [BOTH] | [GAMMA] | [COMPLEX]
- Category: LEAK | ALLOC | UNINIT | BOUNDS | PTR | NUM | HDF5 |
            IFACE | LOGIC | DIVERGE
- Severity: S1 | S2 | S3 | S4   (+ PARALLEL-HAZARD if applicable)
- Status:   open
- Evidence: ...
- Fix:      ...

-->
