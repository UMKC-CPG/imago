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
- **Phase 0 is COMPLETE except 0c (ifort), which stays deferred.**
  Re-verified against the repository 2026-08-09, because the status
  below had gone stale: it still named 0b and 0e as "next" when both
  had in fact been done.
  - *0a* -- the opt-in gfortran instrumentation options are wired
    into the top-level `CMakeLists.txt` and verified: the default
    build is byte-for-byte unchanged and each option injects
    correctly. `IMAGO_CHECKS`, `IMAGO_FPE_TRAP`, `IMAGO_INIT_SNAN`,
    `IMAGO_WARN_EXTRA` and `IMAGO_SANITIZE` are all present.
  - *0b* -- `CMakePresets.json` exists and carries
    `gfortran-release`, `gfortran-debug`, `gfortran-audit` and
    `gfortran-asan`, each binding its own build tree.
  - *0d* -- an instrumented `imagoG` ran a full all-electron Gamma
    SCF to completion with zero AddressSanitizer errors.
  - *0e* -- `BUILD.md` (140 lines, at the repository root, NOT in
    `dev/`) documents every option and preset well enough to pick a
    build cold.
- **Phase 1 (compiler sweep) is SUBSTANTIALLY DONE, 2026-08-09.**
  676 warnings in hand-written source at the start; **116** now,
  with zero errors in either variant. An 83 percent reduction.

  ```
    676  at the start of Phase 1
   -187  array temporaries, MOVED to dev/PERFORMANCE.md as HOT-001
          rather than fixed -- they are a cost, not a defect
   -373  unused entities removed (imports, locals, whole dead
          "use" statements) across 26 files
    116  remaining, per the counts below
  ```

  Remaining, by class:

  ```
    37  -Wunused-variable        all single-variant; need reading
    33  -Wunused-dummy-argument  an argument passed but never used
                                   can mean the routine FORGOT it
    20  -Wcompare-reals          genuine hazard class
    19  -Wconversion             precision / kind, bears on A3
     4  -Wrealloc-lhs            silent shape change on assignment
     2  -Wimplicit-interface
  ```

  Note `-Wimplicit-interface` fell from 22 to 2 not through any
  fix but because the other 20 are in AUXILIARY programs, which
  the `imagoG` and `imago` targets do not build. Those are outside
  this document's stated scope and would need their own pass.

- **Phase 1 correction: build each variant into its OWN log.**
  The first sweep classified each warning's variant by the nearest
  preceding "Building Fortran object" line in a single `-j8` log.
  That is guesswork: under a parallel build the real and complex
  target lines INTERLEAVE. Three entities were removed on the
  strength of it that are genuinely used -- `secularEqn`'s
  `atomKOverlap_did`, `atomKOverlapPlusG_did` and `fullVVDims`
  with their PSCF twins, and `valeCharge`'s `h`, `l` and `skipKP`.
  The compiler caught all three, but that was the build catching
  the method rather than the method working.

  The fix is `make imagoG` and `make imago` into separate logs, so
  a warning's variant is which FILE it is in. Rebuilt that way
  2026-08-09.

- **This also makes the DIVERGE list real** rather than inferred:
  29 sites appear only in the real build and 23 only in the
  complex one. `field.F90` dominates both, and its pairing of
  `valeValeOLGamma` (real only) against `valeValeOL` (complex
  only) is the textbook case of each build using its own type.
  More interesting are six `-Wconversion` sites in `field.F90`
  present ONLY in the complex build: a real-to-complex conversion
  that exists solely on the multi-k path.

- **Phases 2 and 3 have not started. Findings ledger: 1 entry.**
  `build/gfortran-asan` was configured 2026-06-26 and predates
  months of source change; reconfigure before trusting it.
- Phase 2 audit mechanism: multi-agent workflow (decided)

### Tool inventory on this machine (measured 2026-08-09)

Checked rather than assumed, because the campaign plan names tools
without recording whether they are actually reachable here.

```
  valgrind             /usr/bin/valgrind   3.22.0
  gprof                /usr/bin/gprof      2.30
  callgrind_annotate   /usr/bin/callgrind_annotate
  perf                 NOT on PATH
  massif-visualizer    NOT on PATH
  hpcrun               NOT on PATH
  nvidia-smi           NOT on PATH  (no GPU visible from the head node)
```

`module` is available, so `perf` and others may be loadable; that
was not chased. The absence of `perf` matters for the performance
campaign rather than for this one -- see `dev/PERFORMANCE.md`.

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

### Phase 0d validation record (2026-06-26)

The instrumented binary was exercised end-to-end without touching
the installed `bin/` (which stays the stock build):

- **Build:** `cmake -DCMAKE_BUILD_TYPE=Debug -DIMAGO_CHECKS=ON
  -DIMAGO_SANITIZE=address` in a throwaway tree (`build/_0d_asan`);
  `make imagoG` and `make imago` both linked cleanly. The
  instrumented `imagoG` carries 638 `asan` symbols (vs 0 in the
  stock install) and is ~13 MB (vs 3.5 MB).
- **Run:** a plain Gamma SCF on `jobs/knbo3/cubic/debug` (true
  Gamma, `-kp 0 0 0`). imago.py was pointed at the instrumented
  binary through a scratch *overlay bin* -- symlinks to the real
  `bin/` with only `imagoG` replaced -- so the install was never
  modified. `IMAGO_BIN` set to the overlay for the run only.
- **Result:** ran the full 50-iteration SCF to its natural cap
  (`not_converged`, expected for a single Gamma point on this
  perovskite) with **zero AddressSanitizer errors** and no spurious
  traps. Leak detection was disabled for this pass
  (`ASAN_OPTIONS=detect_leaks=0`); the leak hunt itself is Phase 3.

This confirms the Phase-0 harness produces a usable instrumented
binary, so the runtime phase (Phase 3) is unblocked.

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

## Doctrine: drive the warning count to ZERO

Adopted 2026-08-09, during Phase 1, and it changes what the sweep
is for.

A warning that has been investigated once and found harmless does
not stay investigated. It sits in the log, and the next person to
read that log -- or the same person six weeks later -- has to
re-derive the same conclusion from the same code. There are far
too many to hold in a human head, so in practice nobody
re-derives anything: they learn to skim past warnings, and the
one that matters goes past with the rest.

So the target is not "understand the warnings." It is **no
warnings**, by one of three routes, in this order of preference:

1. **Fix the code**, when the warning names something real.
2. **Make the warning structurally impossible**, when the code is
   correct but only by a convention a future edit could break.
   BUG-001 below is the model: the pointers were nullified by
   hand at four sites and it worked; default initialization makes
   it a property of the type instead of a habit, and the warning
   disappears because the situation it warns about can no longer
   arise.
3. **Move the class out of the bug audit**, when it is a real
   diagnostic about something other than correctness.
   `-Warray-temporaries` is the case in point: 187 of them, none
   a bug, every one a hidden array copy that the PERFORMANCE
   campaign wants and this one does not.

What is NOT acceptable is a standing list of known-benign
warnings, because that list is indistinguishable from a standing
list of unread ones.

**A caution that Phase 1 produced immediately.** Silencing a
warning can introduce a defect. The naive fix for BUG-001 -- add
`=> null()` everywhere a pointer is declared -- would have given
two LOCAL pointers the implicit SAVE attribute, making the
routine non-reentrant: a PARALLEL-HAZARD, manufactured by a
campaign that exists partly to find them. Type components take
default initialization; local pointers take an executable
`nullify`. Check what a silencing edit actually means before
making it.

## Findings ledger

Ordered by severity (S1 first).

### BUG-001 -- list-node pointers had undefined association status
- File:     `src/imago/bond3C.F90:16-52` (module O_Bond3C types),
            uses at `:309`, `:371`
- Variant:  [BOTH]
- Category: PTR
- Severity: S4  (+ the silencing edit had a PARALLEL-HAZARD trap)
- Status:   fixed 2026-08-09
- Evidence: `-Wmaybe-uninitialized` on `currentNode2` and
            `currentNode3` in the audit build, twice each (once
            per variant).
- Analysis: **The warning was a false positive as a bug report,
            and a true positive as a fragility report.** The three
            linked-list types declared five pointer components
            with no default initialization, so a freshly allocated
            node had UNDEFINED association status -- and
            `associated()` on an undefined pointer is undefined
            behaviour, not a question with an answer. The code
            compensated by nullifying every one by hand: all list
            heads at `:259-261` before the main loop, and each new
            node's components at `:331`, `:334`, `:412`. Verified
            by reading that this covers every path, so the `else`
            branches at `:309` and `:371` can never be reached
            with an unset cursor. No defect in behaviour.
- Fix:      `=> null()` on the four type COMPONENTS, making the
            invariant structural. An executable
            `nullify (currentNode2, currentNode3)` for the two
            LOCAL pointers -- deliberately not `=> null()`, which
            would confer implicit SAVE and break reentrancy.
            Verified: both variants recompile with the warning
            gone and no new diagnostic.

### BUG-002 -- an unsound test recovers a boolean the caller knew
- File:     `src/imago/integrals3Terms.F90:1057` and `:1377`
            (gaussKOverlap and its sibling)
- Variant:  [COMPLEX]
- Category: NUM  (+ LOGIC)
- Severity: S4 -- latent only; see the consequence note
- Status:   open
- Evidence: `-Wcompare-reals` on `sum(PlusG(:,:)) == 0.0_double`.
- Analysis: `plusG` is a 3x3 matrix and the two callers pass
            either `zeroVectors` or `recipVectors`. The routine
            recovers WHICH call it is by summing all nine
            components and testing against exactly zero. That is
            unsound in principle: a reciprocal lattice whose
            components cancel would sum to zero while being
            nothing like the zero matrix, and the routine would
            take the wrong branch.

            **It is harmless today, and only because the two
            branches are behaviourally IDENTICAL.** Both read the
            same status attribute, both test it the same way, both
            return the same way. They differ in the text of three
            log and error strings and in nothing else. So taking
            the wrong branch currently produces a misleading log
            line and no wrong science.

            That is exactly what makes it worth recording rather
            than dismissing: the moment anyone makes the branches
            differ, an unsound test becomes a real bug, and the
            person making that change has no reason to suspect the
            condition guarding it.
- Fix:      Pass the distinction instead of deducing it. The
            callers already know which case they are; a logical
            argument, or two entry points, is correct by
            construction and needs no floating-point comparison.

### BUG-003 -- zgetrf is called without an explicit interface
- File:     `src/imago/mtop.F90:847`
- Variant:  [BOTH]
- Category: IFACE
- Severity: S3  (+ PARALLEL-HAZARD: none)
- Status:   open
- Evidence: `-Wimplicit-interface`. It is one of only two such
            warnings left in the engine.
- Analysis: Without an explicit interface the compiler cannot
            check argument types, kinds, or shapes at the call
            site, so a mismatch is silent and corrupts memory
            rather than failing to compile. The call looks correct
            as written; the hazard is that nothing enforces it.

            **The codebase already has the convention this
            violates.** `interfaces.F90` wraps seven LAPACK and
            BLAS routines in explicit-interface modules --
            `O_LAPACKZHEGV`, `O_LAPACKDSYGV`, `O_LAPACKDPOSVX`,
            `O_BLASZHER`, `O_BLASDSYR`, `O_BLASZGERC` -- and
            callers `use` them. `zgetrf` is the one that skips it.
- Fix:      Add an `O_LAPACKZGETRF` module to `interfaces.F90`
            following the existing pattern, and `use` it in
            `mtop.F90`. Mechanical, and it brings the last
            unwrapped call into line.

### BUG-004 -- an explicit allocate defeated by auto-reallocation
- File:     `src/imago/intgSaving.F90:779`; same class at
            `src/imago/mtop.F90:370` and `:846`
- Variant:  [BOTH]
- Category: ALLOC
- Severity: S4
- Status:   open
- Evidence: `-Wrealloc-lhs`.
- Analysis: `currentPairGammaTranspose` is allocated to
            `(maxNumStates,maxNumStates)` and then assigned
            `transpose(currentPairGamma)`. In Fortran 2003 and
            later, assigning to an allocatable whose shape differs
            SILENTLY reallocates it to match the right-hand side.
            So the explicit allocate does not constrain anything:
            if the two shapes ever disagree the array quietly
            becomes whatever the expression produced, and the
            allocate reads as an intent the language then ignores.
- Fix:      Either drop the allocate and let the assignment size
            it, which is honest about what happens, or assign into
            an explicit section so the shape is enforced. Choosing
            between them needs a reader who knows whether the
            declared shape is a requirement or a guess.

### BUG-005 -- a unitarity check summed the wrong quantity
- File:     `src/imago/mtop.F90:374` (computeMTOP)
- Variant:  [BOTH]
- Category: NUM
- Severity: S3 -- a diagnostic reports wrong, not the science
- Status:   fixed 2026-08-09
- Evidence: `-Wconversion`, COMPLEX(8) to REAL(8).
- Analysis: The line read
            `idenDiff = idenDiff + sqrt((1.0_double - unitary(m,m))**2)`
            with `unitary` complex, under a comment saying "compute
            the deviation from the identity matrix".

            For a COMPLEX z, `sqrt(z**2)` is not `abs(z)`. It
            returns plus or minus z on the principal branch, and
            assigning that to the real `idenDiff` then discarded
            the imaginary part. So the check accumulated
            `Re(1-u)` -- a SIGNED quantity that can cancel across
            the loop -- where it meant to accumulate a magnitude.
            A badly non-unitary matrix could therefore report a
            small deviation, which is the one thing a check like
            this must never do.
- Fix:      `abs(1.0_double - unitary(m,m))`, the magnitude the
            check was always after.

### BUG-006 -- cmplx() without a kind truncated to single precision
- File:     `src/imago/field.F90:2237` (the neutral-atom
            coefficients); the same idiom at `field.F90` x9 and
            `intgSaving.F90` x2
- Variant:  [COMPLEX]
- Category: NUM
- Severity: S2 at the one live site, S4 at the others
- Status:   fixed 2026-08-09
- Evidence: `-Wconversion`, "REAL(8) to default-kind COMPLEX(4)".
- Analysis: `cmplx(a,b)` returns DEFAULT kind -- single precision
            -- regardless of the kind of its arguments. The third
            argument is not decoration.

            At `field.F90:2237` the arguments were the
            double-precision `neutralCoeffs`, so every digit past
            the seventh was discarded before the square root and
            only then widened back to double. That is a real loss
            in a computed quantity.

            The other eleven uses pass 0 and -1, which are exact
            in single precision, so their VALUE never suffered.
            They were fixed anyway: the habit is what produced the
            live defect, and leaving eleven examples of it in the
            tree invites the twelfth.
- Fix:      `cmplx(..., double)` throughout.

### BUG-007 -- a complex square root whose result is thrown away
- File:     `src/imago/field.F90:2243`
- Variant:  [COMPLEX]
- Category: NUM (+ LOGIC)
- Severity: S3 if reachable, S4 if not -- UNRESOLVED
- Status:   open; the conversion made explicit, the question left
- Evidence: `-Wconversion`, COMPLEX(8) to REAL(8), surviving the
            BUG-006 fix.
- Analysis: `accumWaveFnCoeffsNeut` is REAL, and it is assigned
            `sqrt(cmplx(neutralCoeffs(:),0,double)/2 * weight/spin)`
            -- a complex square root whose imaginary part is then
            silently dropped.

            The `cmplx()` detour only earns its place if
            `neutralCoeffs` can be NEGATIVE, since that is the case
            a real `sqrt` cannot take. But if it ever IS negative,
            the complex square root is purely imaginary and this
            line stores ZERO, losing the value rather than
            reporting anything.

            So the code is either doing something unnecessary or
            hiding a hole, and which one depends on whether a
            fitted neutral-atom charge coefficient can go below
            zero. That is a physics question and is deliberately
            not answered here.
- Fix:      Pending the answer. The assignment now carries an
            explicit `real(...)`, which changes nothing but states
            what was already happening.

### The rest of -Wconversion: implicit, justified, now explicit
The remaining `-Wconversion` sites were all implicit narrowing
that is correct but undocumented, and each now says so:

  - `secularEqn.F90` x3 and `forces.F90` take the real part of a
    DIAGONAL element of a Hermitian matrix, or of a force. Both
    are real by construction, so the discard is exact. Written
    with an explicit `real(..., double)`.
  - `potentialUpdate.F90:801` narrows an 8-byte HDF5 integer to a
    4-byte one; now an explicit `int()`.
  - `mathSubs.f90:353` routed an integer through a real and back
    (`halfk = real(twok)/2.0_double`) purely to truncate. Now
    plain integer division, which truncates identically. Note it
    still truncates when `twok` is odd; whether that can happen
    depends on whether the m values are half-integral, and that
    question was left where it was rather than settled by an edit
    made for a warning.

One incidental find: the `#ifndef GAMMA` around the `plusUJ`
assignment in `secularEqn.F90` had two arms that were
character-for-character IDENTICAL. The branch said nothing and is
gone.

### The rest of -Wcompare-reals: no defects, and now no warnings
All eighteen of the remaining `-Wcompare-reals` sites were read.
None was a defect. They share one pattern: a test against exactly
zero where zero is a SENTINEL some earlier line assigned
literally, not a computed value. `bondLength(i,j) /= 0` means "a
bond was recorded here" in an array zeroed at allocation;
`thermalSigma /= 0` tests an input field the user either set or
left at zero; `multFactor == 0` tests a coefficient the caller
passes as a literal. Comparing a float to a literal it was
literally assigned is exact and safe.

**They were corrected anyway, under the drive-to-zero doctrine.**
A warning that has to be re-adjudicated by hand every sweep is
the thing that doctrine exists to prevent, and "we checked, it is
fine" does not survive into the next reader's head.

The correction changes no behaviour, because gfortran draws this
warning on `==` and `/=` but NOT on the ordered comparisons, and
for these quantities an ordered form asks exactly the same
question:

```
  non-negative quantity   x /= 0   ->  x > 0
    (bond lengths, a smearing width, a Fermi occupation in
     [0,1], a strictly positive nuclear alpha)
  possibly-signed         x == 0   ->  abs(x) <= 0
    (a bond projection, a matrix coefficient, a complex LU
     pivot).  abs is never negative, so it is at most zero only
     when it IS zero -- the same test, exactly.
```

Eight sites in five files: `bond.F90` (four), `matrixSubs.F90`
(two), `gaussRelations.f90`, `populate.F90` (two). Each carries a
comment saying which form was used and why, so the next reader
does not undo it. Verified: both variants report zero
`-Wcompare-reals`, and a full SCF through optical output is byte
identical.

The only two left are BUG-002's, which need the signature fix
described there rather than a rewritten comparison.

Two sites are still worth a look by someone interested in
conditioning rather than correctness, and both are commented in
place. `matrixSubs.F90:49` returns a zero determinant on an
exactly-zero LU pivot, which is an early exit rather than a
singularity test -- a pivot of 1e-300 multiplies through to the
same answer anyway. `populate.F90:620` exits a loop when the
Fermi function reaches exactly zero, which works because IEEE
underflow really does reach exactly zero.

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
