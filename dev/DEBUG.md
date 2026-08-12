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
- **Phase 1 (compiler sweep) is DONE, 2026-08-09.**
  676 warnings in hand-written source at the start; **75** at the
  Phase 1 close, with zero errors in either variant. An 89 percent
  reduction. Counts taken from a FULL recompile of both variants
  into separate logs, not from an incremental build.

  ```
    676  at the start of Phase 1
   -187  array temporaries, MOVED to dev/PERFORMANCE.md as HOT-001
          rather than fixed -- they are a cost, not a defect
   -373  unused entities removed (imports, locals, whole dead
          "use" statements) across 26 files
    -41  the hazard classes: compare-reals, conversion and
          realloc-lhs driven to zero, yielding BUG-004 through
          BUG-007
     75  remaining
  ```

  Remaining, recounted 2026-08-10: after the BUG-007 and BUG-008
  resolutions a full clean recompile of both variants into
  separate logs gave **62 unique warning sites**, 73 warning lines
  (a site flagged by both variants counts once as a site, twice as
  a line). The BUG-002 and BUG-003 fixes later the same day
  removed three more sites (both compare-reals, the one
  implicit-interface), leaving **59 sites, 69 lines** -- that last
  step by arithmetic over the incremental per-variant logs, whose
  recompiled files were checked for new warnings and had none. By
  class, on the unique-site basis -- and none of it is mechanical:

  ```
    39  -Wunused-variable        mostly single-variant; each needs
                                   a read against the other build
    19  -Wunused-dummy-argument  an argument passed but never used
                                   can mean the routine FORGOT it
     1  gfortran "Extension" note: unary minus after `*` without
                                   parentheses (mtop.F90:1076)
  ```

  The earlier tally in this section (75 remaining; 37/33/2/2/1 by
  class) was taken from logs that no longer exist, and its counting
  basis was not recorded, so the two breakdowns cannot be compared
  class by class. The basis is now stated so the next recount is
  comparable. Of the deliberate change between the two, exactly
  three warnings were removed: BUG-008's dummy-argument trio, gone
  with the #if 0 guard.

  So what is left is exactly the part that needs a human rather
  than a rule. Every class that could be settled by inspection has
  been.

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

- **This also makes the DIVERGE list real** rather than inferred.
  On the 2026-08-10 recount: 41 sites appear only in the real
  build, 10 only in the complex one, and 11 in both. `imago.F90`
  and `field.F90` dominate the real-only list, largely imports of
  complex-path HDF5 handles and k-point machinery that the Gamma
  build never touches -- the expected face of the #ifdef
  divergence, but each one still needs its read against the other
  build before it can be called expected.

- **RESUME HERE (2026-08-12).** The Phase 1 warning residue is now
  READ, not just counted: every remaining `-Wunused` site was
  adjudicated in classes this session, each against its callers
  and against the other build. Both variants compile clean with
  **six unique warning sites**, every one explained in place:

  - **The generated-code pair is CLEARED (2026-08-12).** The
    current generator (`src/scripts/osrecurintgana.py`, its
    KOverlap printer) already emits a single, used `preFactorKO`;
    only the checked-in artifact carried the two dead names,
    inherited from whatever older generator state produced it
    before the initial Imago commit. A statement-level
    comparison (continuations joined, spaces stripped, case
    folded) of the artifact's KOverlap subroutine against a
    fresh `-p -ko` generation shows the two texts identical in
    all 2197 statements EXCEPT that declaration -- so the
    artifact's single declaration line was synced to what the
    generator emits, and nothing else. The regeneration work
    also settled two things worth keeping (ledger below):
    **BUG-013** records that `-a` crashes on a missing
    `nuclearbb` function -- known to the programmer, deferred --
    and that the working invocation is the explicit integral
    list, which reproduced the artifact statement-for-statement
    (58,876 statements) except for one loop bound; that bound is
    **BUG-014**, a Boys-order issue the programmer has deferred
    until g-type orbitals join the method, and it must NOT be
    cleared by syncing the artifact to the generator.
  - **Four accepted and documented at the site**: readDataSCF's
    `h` and `numStates` (its argument list deliberately mirrors
    readDataPSCF, which needs both in both builds), ortho3Terms'
    `did2`, and mtop's `inSCF`. Each is a shared routine whose
    argument only the multi-k build consumes; each declaration
    now carries a comment saying which build uses it, why it
    stays, and that the warning is accepted. The gamma log will
    always show these four; that is the recorded end state.

  *What the reading changed (2026-08-12, coded, ledger below):*
  - The unused-variable classes went first: true dead locals
    deleted, variant-divergent locals (`h`, `skipKP`, `xyzP`, the
    readMatrix staging buffers) moved under the guard of the
    build that uses them, with the eigenvector re-read comments
    corrected to say what the variables are actually for.
  - The four k-point allocators were then variant-split into
    honest-signature pairs (`allocateIntegralsSCF` /
    `allocateIntegralsSCFGamma`, likewise PSCF, reallocatePSCF,
    and allocateIntegrals3Terms), matching the readMatrix pair
    convention: the gamma form no longer accepts a `numKPoints`
    it cannot use. The five call sites in `imago.F90` select
    under the same guard, and the gamma-side `numKPoints`
    imports those calls had justified are guarded too.
  - **BUG-012**: `allocateIntegralsForce` had NO callers at all
    (computeForceIntg allocates the same arrays inline) and is
    removed.
  - True dead arguments removed through their call chains:
    `atom1`/`atom2` from the two Full-saving routines in
    `intgSaving.F90` (the Full forms write both triangle blocks
    explicitly, so the diagonal test that needed atom identity
    vanished with the packing), `weight` from `gaussCalc` (only
    its square is consumed), and the `attribInt*` handles from
    the three HDF5 routines that never create attributes -- the
    two access-side EigVec routines and the unused dims of the
    ElecStat init. Dropping the uniform attribInt argument
    bundle across the hdf5 files was an explicit programmer
    decision (2026-08-12), not an oversight.
  - The lone gfortran Extension note (unary minus directly after
    `*` in computeForceGamma) is rewritten as a plain negation,
    which is algebraically identical.

  **Verification status: the paired A/B PASSED (2026-08-12).**
  Twin copies of the reduced KNbO3 deck
  (`jobs/knbo3/o3/ab_unusedsweep_{old,new}`, left in place as
  evidence) ran `imago.py -optc` TWICE each -- the second pass
  deliberately exercises the restart path this batch touched
  (the "already exists, skipping" branch, the access-side HDF5
  routines, the readDataSCF re-reads). Old side used the
  installed Aug-10 binaries; new side used a fresh release build
  of the working tree via `IMAGO_BIN` pointed at a staging dir
  (`jobs/knbo3/o3/ab_stage_bin`). Every physics output
  (all `.plot` and `.dat` files) is byte-identical; the text
  logs differ only in dates and wall-clock lines. The two HDF5
  intermediates differ at the raw-byte level (embedded object
  metadata), so they were adjudicated at CONTENT level instead:
  `h5diff` reports zero differences, and a structural walk of
  all 695 datasets in the two file pairs confirms identical
  trees, shapes, dtypes, and per-dataset allocated storage --
  including that h5diff's 348 "not comparable" objects are
  datasets empty on BOTH sides (declared, never written).

  Coverage caveat, same as every prior A/B on this deck: the
  4 k-point job runs the COMPLEX binary only, so the gamma
  build's changes (its allocator pair members, the negation
  rewrite in computeForceGamma) are verified by clean
  compilation, not by execution; and no force run exists to
  exercise forces.F90 at runtime. The negation rewrite is
  IEEE-exact regardless (sign-bit change vs multiply by -1.0
  are the same operation).

  *How to VERIFY a change against the KNbO3 reference:* do not
  compare a fresh run against the accumulated outputs in
  `jobs/knbo3/o3/reduced` -- those grew through many restarts, and
  a from-scratch rerun differs from them by up to 8.6 in pointwise
  eps2 (converged-potential drift moving sharp peaks; protocol,
  not defect). Instead make TWO fresh copies of the job, run one
  with the old binary and one with the new, and byte-compare the
  pair. Copies must live inside `jobs/`; imago.py sets the copied
  `intermediate` link aside as `intermediateFIXME` and gives the
  copy its own scratch, which is correct behavior. On 2026-08-10
  this A/B came back byte-identical on every data file for the
  BUG-007 + BUG-008 changes; only timestamps differed.

  *How to regenerate the evidence:*
  ```
    cmake --preset gfortran-audit
    cd build/gfortran-audit
    make -j8 imagoG > real.log 2>&1
    make -j8 imago  > complex.log 2>&1
  ```
  Separate logs, never one combined log. About 20 minutes,
  dominated by compiling the 135k-line generated integral file
  twice.

- **Phases 2 and 3 have not started. Findings ledger: 14 entries,
  four real defects fixed** (BUG-001, BUG-005, BUG-006, and
  BUG-009's near miss), **two open** (BUG-010's physics question,
  BUG-011's silent mtop death), **and two known-deferred by the
  programmer** (BUG-013's -a path, BUG-014's eighth Boys order,
  waiting on the g-orbital work). `build/gfortran-asan` was
  configured 2026-06-26 and predates months of source change;
  reconfigure before trusting it.
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
- **BUG-007 physics answer (2026-08-10):** a fitted neutral-atom
  charge coefficient can never be negative (programmer verified),
  so the complex-sqrt detour in `field.F90` was unnecessary and is
  now a plain real square root.
- **BUG-008 disposition (2026-08-10):** the four unreachable
  SCF-blending routines are RETAINED behind an `#if 0` guard, not
  deleted -- they are an Anderson-mixing convergence effort the
  programmer intends to return to. Guard chosen over per-line
  commenting so the bodies stay byte-for-byte intact.
- **Verification protocol (2026-08-10):** byte-identity is judged
  between two fresh copies of the KNbO3 job run with the old and
  new binaries -- never between a fresh run and the accumulated
  reference outputs, which drift with every restart.

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
            (both inside gaussKOverlap)
- Variant:  [COMPLEX]
- Category: NUM  (+ LOGIC)
- Severity: S4 -- latent only; see the consequence note
- Status:   fixed 2026-08-10
- Evidence: `-Wcompare-reals` on `sum(PlusG(:,:)) == 0.0_double`.
- Analysis: `plusG` is a 3x3 matrix and the two callers pass
            either `zeroVectors` or `recipVectors`. The routine
            recovers WHICH call it is by summing all nine
            components and testing against exactly zero. That is
            unsound in principle: a reciprocal lattice whose
            components cancel would sum to zero while being
            nothing like the zero matrix, and the routine would
            take the wrong branch.

            CORRECTED 2026-08-10: the original entry claimed both
            branch pairs were behaviourally identical. That is true
            only at the first site (log strings). At the second the
            branches call `ortho3Terms(8,...)` versus
            `ortho3Terms(9,...)`, and that code selects WHICH HDF5
            dataset family the results are written under -- plain
            KOverlap or KOverlapPlusG. The unsound test was already
            load-bearing: a misfire would silently file results
            under the wrong name.
- Fix:      Pass the distinction instead of deducing it. Done: a
            `plusGVariant` logical argument, set `.false.`/`.true.`
            at the four call sites in `imago.F90`, replaces both
            tests -- correct by construction, no floating-point
            comparison. In the same change the `plusG` matrix
            argument itself was commented out rather than deleted,
            at the programmer's direction: the flag was its last
            live use, but it belongs to the disabled +G pathway
            whose physics is an open question (BUG-010), so the
            full infrastructure -- argument, declarations, callers'
            zeroVectors/recipVectors -- stays in the source,
            commented and documented, ready to re-instate. Verified
            by paired A/B: byte-identical data outputs.

### BUG-003 -- zgetrf is called without an explicit interface
- File:     `src/imago/mtop.F90:847`
- Variant:  [BOTH]
- Category: IFACE
- Severity: S3  (+ PARALLEL-HAZARD: none)
- Status:   fixed 2026-08-10
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

            Done: interface-only module (no solver wrapper --
            zgetrf takes no work arrays, so unlike the seven
            solver modules there is no ceremony worth hiding), the
            `external :: zgetrf` declaration removed from
            matrixDet, and the call now compiler-checked. The
            engine's `-Wimplicit-interface` count is zero.

### BUG-004 -- an explicit allocate defeated by auto-reallocation
- File:     `src/imago/intgSaving.F90:793`; same class at
            `src/imago/mtop.F90:370` and `:853`
- Variant:  [BOTH]
- Category: ALLOC
- Severity: S4
- Status:   fixed 2026-08-09
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
- Fix:      A section assignment, `x(:,:) = ...`, at all three
            sites. A section cannot resize, so the allocated shape
            is enforced and any future mismatch becomes an error
            rather than a silent resize -- which is what each
            explicit `allocate` was trying to express in the first
            place.

            Choosing this over "drop the allocate" was possible
            because all three shapes are provably equal today, and
            that was checked rather than assumed:
            `currentPairGamma` is DECLARED
            `(maxNumStates,maxNumStates)` so its transpose is too;
            `stateStateMat` is allocated
            `(maxOccupiedState,maxOccupiedState,spin)` against a
            `unitary` of `(maxOccupiedState,maxOccupiedState)`; and
            `Ac = A` sits behind a guard that returns unless
            `size(A,1) == size(A,2)`, with `n` taken from A. So the
            edit preserves behaviour today and tightens what
            happens tomorrow.

            Verified: zero `-Wrealloc-lhs` in both variants, and a
            full SCF through optical output byte identical.

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
- Severity: S4 -- the silent-zero hole was unreachable
- Status:   closed 2026-08-10; physics question answered
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
            zero. That is a physics question and was deliberately
            left to the programmer.

            ANSWERED 2026-08-10: the programmer verified the
            coefficient can never be negative -- which the site's
            own construction agrees with, since each entry is an
            electron count spread evenly over the QN_m orbitals of
            one QN_l. The complex detour was unnecessary, not a
            hole.
- Fix:      The assignment is now a plain real `sqrt` with the
            non-negativity invariant recorded in a comment at the
            site. Verified numerically inert by the paired A/B
            protocol (RESUME HERE above): byte-identical data
            outputs against the pre-change binary.

### BUG-008 -- two potential-blending routines are unreachable
- File:     `src/imago/potentialUpdate.F90:2033-2313`
            (blendPotentialsSCF) and `:2318-2557`
            (blendPotentialsTE)
- Variant:  [BOTH]
- Category: LOGIC (dead code)
- Severity: S4 -- nothing misbehaves; roughly 520 lines mislead
- Status:   resolved 2026-08-10; guarded out, kept for a return
- Evidence: three `-Wunused-dummy-argument` warnings that made no
            sense -- `firstTerm` unused while `numTerms` is used,
            and `totalEnergyRecord` unused by a routine whose whole
            purpose is blending by total energy.
- Analysis: The explanation is that neither routine is ever
            called. Every call site is commented out
            (`potentialUpdate.F90:1540, 1544, 1550, 1559, 1573`)
            and a search of the whole tree finds nothing but the
            definitions and their `end subroutine` lines.

            That also explains `firstTerm` specifically: the dead
            call sites all pass a literal 1, so a body that ignores
            it and starts from 1 was never wrong.

            The hazard is not misbehaviour, it is 520 lines of
            plausible physics that a reader will assume runs.
- Fix:      A decision rather than an edit, and the programmer's to
            make: delete them, or mark them clearly as retained on
            purpose. Deleting 520 lines of potential-blending logic
            is not something to do on the strength of a warning
            sweep.

            DECIDED 2026-08-10: the routines are part of an effort
            to improve SCF convergence (Anderson mixing) that never
            reached working form and WILL be returned to, so they
            are retained. The resolution also swept in
            `shiftPotentials` and `blendJointPotentials`, whose
            only call sites are commented out too -- the same
            effort, found unreachable by the same search. All four
            now sit inside a single `#if 0` / `#endif` guard with a
            banner comment naming the effort, the reactivation
            steps (delete the two guard lines, restore the call
            sites) and this entry. The guard was chosen over
            per-line commenting so the bodies stay byte-for-byte
            intact; gfortran was verified to emit nothing for the
            skipped block. The three nonsensical dummy-argument
            warnings are gone, and the paired A/B run (RESUME HERE
            above) confirmed the build is numerically unchanged.

### BUG-009 -- an unused argument that must NOT be "fixed"
- File:     `src/imago/bond3C.F90:717` (compute3CBO)
- Variant:  [BOTH]
- Category: IFACE
- Severity: S4 as found; S2 if someone had corrected it wrongly
- Status:   fixed 2026-08-09
- Evidence: `-Wunused-dummy-argument` on `kPointWeight` and `spin`
            in a bond-order routine -- which reads like the routine
            forgot to weight its result.
- Analysis: **It did not.** The weighting is already inside
            `chargeScaleFactor`, which the caller sets from
            `electronPopulation`, and that array carries the factor
            `kPointWeight/spin` within it. The optical code shows
            the same fact from the other side: `computePairs`
            DIVIDES by `kPointWeight/spin` to recover a plain
            zero-to-one occupancy from the same array.

            This is the reason the bucket was read rather than
            cleared mechanically. The obvious "fix" -- multiply the
            accumulation by `kPointWeight` since it is right there
            unused -- would count the weight TWICE and silently
            corrupt every bond order on a multi-k run. An unused
            argument can be an invitation to a bug rather than
            evidence of one.
- Fix:      Both arguments removed from the signature and from the
            three call sites, with the reasoning recorded at the
            routine so the question cannot be re-asked from a
            signature that still offers them.

### BUG-010 -- the KOverlapPlusG datasets duplicate plain KOverlap
- File:     `src/imago/integrals3Terms.F90` (gaussKOverlap, the
            disabled KOverlap2CIntg form); consumer at
            `src/imago/mtop.F90:286-294` via matrix codes 6-8 in
            `secularEqn.F90`
- Variant:  [COMPLEX]
- Category: LOGIC (possibly NUM if the answer is "wrong")
- Severity: S2 if the wrap step needs the shift; S4 if redundant
            -- UNRESOLVED, a physics question for the programmer
- Status:   open; found 2026-08-10 while fixing BUG-002
- Evidence: Reading, not a warning: the only line that ever used
            the `plusG` argument physically -- adding it to the
            k-displacement passed to KOverlap2CIntg -- is commented
            out, and has been since before the OLCAO import (the
            initial Imago commit already carries it disabled).
- Analysis: With the shift disabled, both gaussKOverlap
            invocations compute identical integrals, so the
            KOverlapPlusG datasets are byte-duplicates of plain
            KOverlap under a different name. Yet mtop reads the
            PlusG datasets deliberately, at exactly one place: the
            final link of each k-point string, where the walk
            wraps across the Brillouin zone boundary and the true
            displacement differs from an interior step by a full
            reciprocal lattice vector. So either (1) the wrap step
            uses an unshifted overlap where a G-shifted one
            belongs -- a per-string error of just the sort the
            BUG-005 unitarity check monitors -- or (2) the LCAO
            phase convention makes the two matrices genuinely
            equal, and the whole PlusG apparatus (second dataset
            family, matrix codes 6-8, the duplicate invocation) is
            redundant scaffolding.
- Fix:      Blocked on the physics answer. The revival recipe is
            preserved IN THE SOURCE, commented and documented, at
            the programmer's direction: the alternative signature
            and plusG declaration in gaussKOverlap, the plusG form
            of the four calls in imago.F90 with their zeroVectors
            and recipVectors support, and the shifted
            KOverlap2CIntg call itself, each annotated with a
            pointer here. A physics note above the KOverlap2CIntg
            call carries the full explanation. If the answer is
            (2), delete the apparatus instead. Note the question
            cannot be settled by running mtop until BUG-011 is
            fixed.

### BUG-011 -- mtop dies silently after tetrahedron construction
- File:     unknown; last output from `src/imago/kpoints.f90`
            (the RESOLVED_KP_CLASSES emitter near `:1312`)
- Variant:  [COMPLEX] (mtop is #ifndef GAMMA)
- Category: LOGIC (+ a fail-loudly violation)
- Severity: S2 -- the mtop feature is unusable on this deck
- Status:   open; found 2026-08-10 by the first live `-mtop` run
            ever recorded (no `command` file in the tree had one)
- Evidence: `imago.py -mtop` on the KNbO3 reduced deck: the
            Fortran run dies with no error message ("Fortran
            success file missing" is all the user sees), and its
            fort.20 ends immediately after the tetrahedron
            construction check. The RESOLVED_KP_CLASSES line
            prints `20 -1927434624 22014` on one run and
            `20 -1615442304 22047` on the next -- values that
            change between runs on the same binary are
            uninitialized memory being printed.
- Analysis: Not yet performed; recording the reproduction. The
            failure is bit-identical in shape on the pre- and
            post-BUG-002/003 binaries, so it is pre-existing and
            unrelated to those fixes. Two independent defects are
            visible already: whatever kills the run, and the
            uninitialized values reaching a structured output line
            that downstream tooling parses.
- Fix:      Pending investigation. Start at the
            RESOLVED_KP_CLASSES emitter in kpoints.f90 and at
            whatever runs next after the construction check; an
            IMAGO_CHECKS or asan build of the mtop path would
            likely name the death site directly.

### BUG-012 -- a public allocator no caller has ever used
- File:     `src/imago/forces.F90:39` (allocateIntegralsForce)
- Variant:  [BOTH]
- Category: IFACE
- Severity: S4
- Status:   fixed 2026-08-12 (removed)
- Evidence: `-Wunused-dummy-argument` on `numKPoints` led to the
            caller read the dummy-argument class requires, and
            the read found NO callers anywhere in the tree,
            case-insensitive.
- Analysis: computeForceIntg performs the identical three force
            matrix allocations inline (`forces.F90:212`-`218` at
            the time of removal), so the routine was superseded
            plumbing that survived because nothing referenced it.
            Read as a warning site alone it would have passed for
            one more variant-divergence dummy -- `numKPoints` IS
            used by its multi-k arm -- and the accepted-warning
            comment it would have received would have documented
            a routine that does not exist in any execution path.
            Only the caller search showed the truth.
- Fix:      Routine deleted. The inline allocations in
            computeForceIntg are unchanged and remain the single
            allocation site for the force matrices.

### BUG-013 -- the generator's -a path is broken (known, deferred)
- File:     `src/scripts/osrecurintg.py:1385` (main) against
            `src/scripts/osrecurintgana.py`; artifact is
            `src/imago/gaussIntegrals.f90`
- Variant:  [BOTH] (build tooling, not the engine)
- Category: IFACE
- Severity: S4 -- the working invocation is the explicit list
- Status:   recorded 2026-08-12; the programmer confirmed the
            breakage is known and will be resolved later. Not an
            open engine defect.
- Evidence: `osrecurintg.py -p -a` (production output, all
            integral types) dies with AttributeError: module
            'osrecurintgana' has no attribute 'nuclearbb'.
            Per the programmer, `-a` is not how the artifact is
            produced: the program is called with an explicit
            list of the required integrals. The invocation that
            matches the artifact's twelve subroutines, verified
            2026-08-12 to run to completion and reproduce them
            in the same order, is:
            `osrecurintg.py -p -o -k -e -n -m -mv -d -dk -dncb
            -decb -ko`
- Analysis: with that explicit regeneration in hand, the whole
            135k-line artifact was compared against it at
            normalized statement level (continuations joined,
            spaces stripped, case folded): all 58,876 statements
            are IDENTICAL except two. One was the KOverlap
            prefactor declaration -- the Class 4 warning pair,
            synced to the generator's form. The other is a loop
            bound in dnuclear3CIntgCB where the artifact and the
            generator are BOTH defective in different ways; that
            is BUG-014, and it must NOT be "fixed" by syncing
            the artifact to the generator. The case and
            line-wrap style of the KOverlap section also differs
            from the current generator's output (the artifact
            predates the initial commit), which any future
            regeneration diff will show as bulk churn.
- Fix:      Deferred by the programmer. Until then, regenerate
            with the explicit list above, never `-a`.

### BUG-014 -- dnuclear3CIntgCB's eighth Boys order (deferred to g)
- File:     `src/imago/gaussIntegrals.f90:88127` (the m loop in
            dnuclear3CIntgCB) and the `boys` routine above it;
            generator side in `src/scripts/osrecurintgana.py`
- Variant:  [BOTH]
- Category: UNINIT (artifact) / BOUNDS (generator output)
- Severity: S4 today by the programmer's determination; S2 the
            day g-type orbitals join the basis
- Status:   recorded 2026-08-12; known to the programmer,
            deferred until g-type orbitals are added to the
            Imago method (not for a while)
- Evidence: the single semantic difference the BUG-013
            statement-level comparison found in 58,876
            statements: the artifact initializes the Boys-order
            prefactors with `do m = 1, 7` where the current
            generator emits `do m = 1, 8`. In the artifact,
            `preFactorN` is dimensioned (8), filled to 7, and
            READ at index 8 in three statements of the final
            l1l2switch branch (switch 136), so those three
            products draw an uninitialized value. The
            generator's form is no fix: its loop reads `F(8)`
            while `F` is dimensioned (7) and `boys` exports only
            seven orders -- an out-of-bounds read in place of an
            uninitialized one. Inside `boys`, the small-T series
            path computes an eighth order internally (local `S`
            is dimension (8)) but discards it, and the analytic
            large-T path has closed-form expressions only
            through `F(7)`, so exporting order 8 requires
            deriving that expression, not just widening arrays.
            Corroboration: the optimized release build's
            `-Wmaybe-uninitialized` independently flags the
            artifact's read (`gaussIntegrals.f90:104198`) in
            both variants' compiles.
- Analysis: the programmer knows this issue: the top-order terms
            exist to support the derivative recursion's raised
            angular momentum, and the defect only becomes live
            once g-type orbitals are included in the method,
            which is future work. Until then the affected
            entries do not reach the returned integrals. Note
            the force path that calls dnuclear3CIntgCB is also
            currently unexercised by any recorded run.
- Fix:      Deferred with the g-orbital work, as one coordinated
            change: derive the analytic F(8), export eight
            orders from `boys` with `F` dimensioned (8), keep
            the generator's `do m = 1, 8`, and regenerate. Do
            NOT sync the artifact's loop bound to the generator
            before `boys` is extended -- that converts the
            uninitialized read into an out-of-bounds read.

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
