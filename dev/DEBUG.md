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

- **RESUME HERE (2026-08-16, mid-session).** Where things
  stand and what is next, in order:

  1. **Phase 1 is CLOSED** (commits `8c68021`, `478d1dd`,
     `3e170d6`): audit-build warnings at their recorded end
     state of four accepted-and-documented sites; A/B passed;
     binaries installed.
  2. **Phase 2 was resequenced** (decisions log): tool-seeded
     evidence before open-ended reading. Step 1 -- the
     release-build `-Wuninitialized` tranche -- is DONE
     (section above): 34 groups adjudicated, 30 benign or
     already known, three accepted into the ledger as
     **BUG-015 (ETA), BUG-016 (gamma force dump), BUG-017
     (spin clamp)**.
  3. **BUG-015/016/017 are FIXED and A/B-VERIFIED 2026-08-12.**
     015 became `GZ = 0` with the reference form and the ETA
     declaration retained as commented-out code (programmer's
     choice); 016 became the `(:,l,j,i)` index fix, compile-
     verified only since no recorded run reaches the gamma force
     dump; 017 gained the `spin == 2` guard. The post-fix release
     rebuild is clean of all three warnings. The
     paired `-optc` x2 A/B (`jobs/knbo3/o3/ab_uninit_{old,new}`,
     staged bin `jobs/ab_stage_bin`) was byte-identical on every
     data file with only wall-clock timestamps differing, and the
     HDF5 files content-identical: h5diff clean, structural walks
     matched at 143 and 676 objects, and every "not comparable"
     object is an empty dataset present identically in both
     files. Verified binaries installed to `$IMAGO_DIR/bin`,
     cmp-identical to the A/B builds. The evidence copies
     (`ab_uninit_{old,new}`, `ab_stage_bin`, and their scratch
     dirs) were deleted at the programmer's direction
     2026-08-12; the from-clean rebuild of `build/release`
     (2026-08-13) replaced the post-fix logs with the current
     `uninit_{real,complex}.log` pair.
  4. **The standing release-warning list is machine-guarded
     (2026-08-13).** Decision recorded in the log below: the
     release build's adjudicated-benign `-Wmaybe-uninitialized`
     groups plus the four accepted unused-dummy sites are
     encoded in `dev/tools/release_warning_manifest.tsv` (44
     per-variant rows: 29 distinct uninit groups union'd across
     variants, 4 dummy sites, gamma build only), and
     `dev/tools/check_release_warnings.py` diffs a fresh
     per-variant from-clean build log against it -- any new or
     vanished warning prints loudly and exits nonzero. The
     2026-08-13 logs match the manifest, and both failure
     directions were test-fired with doctored logs.
  5. **Phase 3 runtime instrumentation is IN PROGRESS; the
     2026-08-13 candidates are all adjudicated and FIXED
     2026-08-14 as BUG-018 (PSCF radial-fn leak), BUG-019
     (PACS uninitialized conversion), and BUG-020 (stale
     `check_gamma_kp` offsets)** -- ledger entries carry the
     fixes and verification. Fixing them widened the harvest
     into a second wave (harvest section, 2026-08-14 wave), of
     which the SNaN blocker is now adjudicated and FIXED
     2026-08-16 as **BUG-021** (the integral reader loops read
     `oneAlphaPair` columns above the current alpha's lm
     coverage; fixed by a per-state skip at all six reader
     sites plus block-restricted hand-backs in `nuclearPE` and
     `electronicPE` -- strictly less work than before). **The
     whole SCF is now SNaN-clean on both variants**, and the
     two SNaN runs that showed it double as the batch A/B for
     BUG-018/019/021: line-for-line identical trajectories
     against unfixed release controls on both variants
     (`jobs/knbo3/o3/phase3_b21_ctrl`, `phase3_asan_g`).
     Nothing is committed or installed yet. **The imagoG
     "divergence" candidate is RESOLVED: not a defect in either
     variant.** Its bisection found the first differing quantity
     was the k-point itself and led to **BUG-022 (FIXED
     2026-08-16 through the chain)**: the mesh formula's `-1/2`
     offset inverted the shift semantics for odd counts, so the
     canonical Gamma deck sampled the zone corner in the complex
     binary. Given a true Gamma point the complex binary
     reproduces imagoG line for line -- imagoG's numerics are
     validated for the first time -- and the Gamma-only SCF
     oscillation is real behaviour of that sampling. The last
     candidate of the wave, the stale SYBD gamma demotion, is
     **BUG-023 (FIXED 2026-08-16, verified live: a gamma-deck
     band structure now completes on the general binary)**.
     BUG-018..022 are committed (`6377ba9`) and BUG-023
     (`0ba89ff`); both binaries and `imago.py` were INSTALLED
     and the commits PUSHED by the programmer 2026-08-16 (the
     installed complex binary verified Gamma-centred on a
     3x3x3 mesh). **The coverage sweep then ran 2026-08-16:
     all six runnable families x both decks x both instrumented
     variants (24 cells; matrix and per-cell notes in the
     harvest section). `-dos` and `-bond` are clean everywhere;
     `-loen`, `-sybd`, `-force`, `-mtop` each fail with a
     located, deterministic signature -- FIVE CANDIDATES (A-E)
     await programmer review in the "2026-08-16 wave"
     subsection, no BUG numbers. E is BUG-011's mechanism (the
     MTOP path never sets up the point ops that `buildAtomPerm`
     reads).** No sanitizer report appeared in the sweep. The
     harvest section also records three harness traps
     (`$IMAGO_BIN` must be exported for every overlay run --
     three void "verification" runs happened without it; the
     overlay bins hold COPIES that must be refreshed after each
     rebuild; a re-created deck of the same name silently
     inherits the dead deck's scratch), and the discovery that
     `imago.py` has a built-in valgrind mechanism that may
     settle the still-undecided valgrind approach. **NEXT:
     adjudicate candidates A-E one at a time (suggested order
     E, A, B, C, D), fix, re-run the failing cells; then the
     `-optc` SNaN/gamma cells still marked todo; then decide
     the valgrind approach.** Working state at session end:
     tree clean except this DEBUG.md update; overlay bins are
     current with the committed source (rebuilt 2026-08-16
     16:49); the four sweep decks hold their run logs.
  6. **Then reassess the shrunken fan-out** (step 3 of the
     resequencing). The preprocessed variant texts prepared for
     it lived in session scratch and are gone; regenerate with
     `gfortran -cpp -E [-DGAMMA] <file>` per variant when
     needed.

  Evidence locations: the current release warning logs
  (post-BUG-015/016/017 state, regenerated from clean
  2026-08-13) in `build/release/uninit_{real,complex}.log` --
  the pre-fix tranche content survives as the grouped table in
  the tranche section below; the Phase 1 class logs in
  `build/gfortran-audit/class*_{real,complex}.log`.

  *The Phase 1 record (kept for reference):* every remaining
  `-Wunused` site was adjudicated in classes, each against its
  callers and against the other build. Both variants compile
  clean with six unique warning sites at the time of closure,
  every one explained in place:

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
  (`jobs/knbo3/o3/ab_unusedsweep_{old,new}`, since inspected and
  deleted along with their scratch) ran `imago.py -optc` TWICE
  each -- the second pass
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

- **Phase 2 step 1 (release-uninit tranche) is DONE and its
  three findings are FIXED. Findings ledger: 17 entries, seven
  real defects fixed** (BUG-001, BUG-005, BUG-006, BUG-009's
  near miss, and the tranche's three: BUG-015 ETA, BUG-016
  gamma force dump, BUG-017 spin-clamp guard -- fixed and
  A/B-verified 2026-08-12, binaries installed), **two open**
  (BUG-010's physics question, BUG-011's silent mtop death),
  **and two known-deferred by the programmer** (BUG-013's -a
  path, BUG-014's eighth Boys order, waiting on the g-orbital
  work). Next: Phase 3 runtime instrumentation.
  `build/gfortran-asan` was deleted in the 2026-08-12 cleanup;
  reconfigure from the preset when Phase 3 needs it.
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

## Release-build uninitialized tranche (harvested 2026-08-12)

The resequenced Phase 2's step 1. Fresh from-scratch release
(-O3) builds of both variants, separate logs
(`build/release/uninit_{real,complex}.log`), yield 45 raw
`-Wuninitialized`/`-Wmaybe-uninitialized` sites that collapse to
**34 (file, variable) groups**: 1 definite, 33 maybes. Group
counts by file: field.F90 21, potentialUpdate.F90 5, forces.F90
2, optcSpectra.F90 2, and one each in optc.F90, bond.F90,
dos.F90, gaussIntegrals.f90 (that last is BUG-014's
`preFactorN`, already adjudicated). The full grouped table is
reproducible by re-parsing the two logs.

Adjudication complete (2026-08-12). The programmer reviewed the
three candidates and accepted all three into the ledger:

- `ETA` uninitialized in the PBE correlation copy -> **BUG-015**
- gamma force dump's stale-index OOB read -> **BUG-016**
- SXS/SYS/SZS clamp missing its spin guard -> **BUG-017**
- **28 remaining maybes adjudicated BENIGN, plus one already
  known (gaussIntegrals `preFactorN` = BUG-014).** Every benign
  group is the same flow-analysis blind spot: setup and use
  guarded by the same module-variable condition, which the
  optimizer must assume could change between reads. The
  clusters: all 21 field.F90 groups (allocation, initialization,
  accumulation, and deallocation each under matching
  `doPsiFIELD`/`doWavFIELD`/`doRhoFIELD`/`doPotFIELD` tests, or
  scalars set and used under the OR of those flags);
  potentialUpdate `currentPot` (`rel == 1` on allocate, use, and
  deallocate); optc `momentPairTemp` (allocated for
  `numStoredCompPOPTC > 1`, written in the `else` of `== 1`,
  read under `> 1`, deallocated behind `allocated()`); forces
  `zFactor` (assigned and consumed inside the same
  nuclear-potential branch of the term loop); optcSpectra
  `sigma`/`energyDelta` (assigned for doOPTC 1 and 2; the
  routine's comment documents that Sigma(E) runs never call it);
  dos `pdosAccum` (Gaussian path allocates and uses; the LAT
  path allocates its own inside computeProjections_LAT, as
  commented at the allocation); bond `bondCompleteAtom`
  (`excitedAtomPACS /= 0` guards allocation, zeroing, and
  output).

Post-fix standing list (2026-08-13): with BUG-015/016/017
fixed, the union across variants is **29 (file, variable)
groups** -- the 34 above minus ETA, forces `q`, and
SXS/SYS/SZS. That list, plus the four Phase 1 unused-dummy
sites, is encoded per variant in
`dev/tools/release_warning_manifest.tsv` and checked by
`dev/tools/check_release_warnings.py` (decision below).

## Phase 3 runtime harvest (started 2026-08-13, in progress)

Harness: asan binaries rebuilt from the `gfortran-asan` preset
(both variants link clean, 676/677 asan symbols), run through an
overlay bin (`jobs/phase3_stage_bin` -- symlinks to the install
with only `imagoG`/`imago` replaced, selected via `IMAGO_BIN`) on
fresh deck copies `jobs/knbo3/cubic/phase3_asan` (from
`cubic/debug`) and `jobs/knbo3/o3/phase3_asan` (from
`o3/reduced`). Leak detection ON -- Phase 0d had it explicitly
off. Two operational notes from getting the decks running: an
`imago.dat` generated before the O3/O11 optical-direction fields
(2026-08-08) no longer parses -- old decks need `makeinput.py`
rerun before reuse -- and a copied deck restarts from its
converged `gs_scfV` file, so scratch-coverage runs must set that
file aside or they exercise a single iteration.

AddressSanitizer + LeakSanitizer results (KNbO3 decks). NOTE:
all four runs exercised the COMPLEX binary -- including the two
on the cubic deck that were briefly recorded here as gamma runs
-- because `imago.py` never selects `imagoG` anymore (the
executable-selection candidate below). `imagoG` therefore has
had NO runtime coverage yet in this phase.

- cubic deck SCF, fresh 13-iteration run to convergence: CLEAN.
- o3 deck SCF, fresh 15-iteration run to convergence: CLEAN.
- o3 deck `-optc` (pscf+optc, 4m02s): **one leak** (below).
- cubic deck `-optc`: the SAME leak, byte-for-byte (6480 bytes
  in 6 allocations), so the defect is path-dependent (SCF clean,
  PSCF leaks). Whether it is also variant-independent is
  UNVERIFIED until `imagoG` can actually be run, though the
  implicated code (`atomicTypes.f90`, `cleanUpPSCF`) is shared
  between the variants.

**ACCEPTED as BUG-018 and FIXED 2026-08-14 (see the ledger entry
for the fix): every PSCF-family run leaked the atomic-type
radial-function arrays.**
`readAtomicTypes` allocates the pointer components
`coreRadialFns` and `valeRadialFns` per atomic type
(`atomicTypes.f90:240` and `:265`; 3 types x 2 arrays = the 6
leaked allocations). On the SCF path they are freed mid-setup by
`cleanUpRadialFns`, called only from `setupSCF`
(`imago.F90:414`). The PSCF path never calls it, and
`cleanUpPSCF` then calls `cleanUpAtomTypes`, which deallocates
`atomTypes` itself (`atomicTypes.f90:530`) without touching the
two radial-fn components -- the pointers are orphaned and
LeakSanitizer reported a direct leak. Exit-time and small on
this deck, but structural: it scales with atom-type count and
would bite any future path that tears down and re-reads types
inside one process. Fixed by making `cleanUpAtomTypes` a
complete guarded teardown (the SCF's early shed stays); the
post-fix asan rerun on this deck is sanitizer-silent.

SNaN+FPE results (`build/gfortran-snan`, configured Debug +
`IMAGO_CHECKS` + `IMAGO_FPE_TRAP` + `IMAGO_INIT_SNAN`; decks
copied to `phase3_snan` with converged potentials set aside):
every run on both decks trapped at PARSE TIME, at the same
site, so the pass could not see past input reading until that
candidate was adjudicated and fixed (it since was, as BUG-019;
the pass now reaches integral setup and traps at the
`gaussOverlapNP` candidate in the 2026-08-14 wave below). The
trap was the finding:

**ACCEPTED as BUG-019 and FIXED 2026-08-14 (see the ledger
entry): `readPACSControl` converted an uninitialized real when
PACS was off.** The unit
conversion `totalEnergyDiffPACS = totalEnergyDiffPACS / hartree`
(`input.f90:785`) runs unconditionally, but the variable is only
assigned inside the core-state loop above it, and only for the
entry matching the excited QN pair. With PACS off (these decks:
`excitedAtomPACS = 0`, `numCorePACS = 0`) the loop never runs,
so the conversion reads an uninitialized real -- a SIGFPE under
SNaN init, silently-stored garbage in production builds. Benign
today only if nothing reads the value when PACS is off; the same
shape as BUG-015's silently-right-by-luck ETA. Fix shape:
initialize `totalEnergyDiffPACS` (and review its three sibling
`*PACS` state variables) before the loop, or guard the
conversions on a matched entry.

**ACCEPTED as BUG-020 and FIXED 2026-08-14 (see the ledger
entry): `imago.py` could no longer select `imagoG` -- every run
since the July k-point format change silently used the complex
binary.** The k-point
file grammar gained `NUM_TETRA_DIAGONALS` and
`SYMMETRIZE_LAT_PARTIALS` for ALL style codes (`readKPoints`,
`kpoints.f90:274-276`, read before the style branch), shifting
every later field down four lines. `check_gamma_kp` in
`imago.py` still reads the mesh counts and shift at their old
fixed offsets (`lines[5]`, `lines[7]`), so a style-1 gamma file
(`1 1 1`, zero shift) is misread as `[4]` vs `[1,1,1]` and the
function returns False; style-0 files misread the same way. The
docstring says the layouts are mirrored from `readKPoints` --
they no longer are. Consequences: gamma decks silently run the
slower complex binary (results correct, speed lost), and the
gamma executable has had zero live coverage since the change.
Worse, with `-scftetradiag 1` a style-0 file would make
`check_gamma_kp` read the SYMMETRIZE value line as a k-point
line and crash on a missing token. Fix shape: parse the file by
LABEL rather than fixed line offsets, so the next grammar change
cannot silently break the mirror again. Found because the
SNaN crash line printed the invoked command
(`.../imago 2 0 0 0 0 0`) for a deck regenerated with
`-kp 0 0 0`.

Operational corollary, recorded for honesty: Phase 0d's June
validation ran `imagoG` only because the old grammar predated
the drift; and this phase's cubic-deck runs are complex-binary
runs (labels corrected above).

### 2026-08-14 wave: adjudication, and what the fixes uncovered

All three 2026-08-13 candidates were reviewed and accepted as
BUG-018, BUG-019, and BUG-020 (mechanisms, fixes, and
verification in their ledger entries). Restoring `imagoG`
selection and unblocking the SNaN pass immediately produced a
second harvest wave: the three NEW candidates below.

Two harness traps burned during verification, recorded so they
are never re-learned:

- **`bin_dir` comes from `$IMAGO_BIN`, never from which copy of
  `imago.py` was invoked.** Every overlay run must export
  `IMAGO_BIN` to the overlay directory. Invoking the overlay's
  `imago.py` without it silently runs the INSTALLED release
  binaries, which trap nothing and leak-report nothing, so a
  "verification" run passes vacuously. Three such void runs
  happened before the trap was caught -- by a stale trap line
  number and by a 36-second "asan" run that should have taken
  four minutes. Judge an instrumented run only by markers an
  instrumented binary alone can produce: sanitizer reports,
  backtraces, characteristic run time, or a `pgrep` of the
  running process.
- **The overlay bins hold COPIES of the two binaries, not
  symlinks.** A rebuild does not reach them; refresh the copies
  after every rebuild. (The overlay `imago.py` IS a symlink,
  and now points at the repo source copy, so script fixes take
  effect without an install.)
- **`rm -rf <deck>` does not remove the deck's scratch, and a
  new deck of the same name silently inherits it.** The
  `intermediate` symlink points under `$IMAGO_TEMP`; deleting the
  deck leaves that directory and its `gs_scf-fb.hdf5` in place,
  and a fresh deck created under the same name finds the old
  integral sets marked complete and REUSES them. Caught
  2026-08-16 when a re-created deck reproduced a different
  k-point's run to eight digits in one second. Before re-creating
  a deck name, `rm -f $(readlink -f <deck>/intermediate)/*`
  first (or use a new name); judge freshness by run time and by
  the integral timestamps in the log.

One tooling note for the pending valgrind pass: `imago.py`
already carries a built-in valgrind mechanism (`settings.
valgrind` wraps execution in `time valgrind --leak-check=yes`),
which may make the rejected wrapper-script overlay unnecessary;
the approach still awaits the programmer's decision.

**Coverage matrix -- the plan of record for this phase.** The
job menu is twelve post-SCF families, each in `-X` (pscf) and
`-scfX` forms, plus the plain SCF; instrumented coverage so far
touches two families. Legend: ok = run clean; B18 = ran and
found BUG-018, clean after the fix; B21 = trapped on BUG-021,
clean after the fix (2026-08-16: the whole SCF is SNaN-clean on
both variants); A..E = failed on the 2026-08-16 sweep with the
signature of that letter in the wave below (candidates awaiting
review); gen = a gamma deck routes this family to the general
binary (BUG-023), so the cell IS the general binary; deck =
needs deck work first. Valgrind is a whole pending column
(approach undecided). The `-scfX` forms run their analysis from
different call sites than the `-X` forms; whether they need
their own rows is an open scope question.

    family    imago-asan  imago-snan  imagoG-asan  imagoG-snan
    SCF       ok          B21         ok           B21
    -optc     B18         todo        todo         todo
    -dos      ok          ok          ok           ok
    -bond     ok          ok          ok           ok
    -loen     A           A           A            A
    -sybd     B           B           C (gen)      C (gen)
    -force    D           D           D            D
    -mtop     E           E           E            E'
    -pacs     deck        deck        deck         deck
    -nlop     deck        deck        deck         deck
    -sige     deck        deck        deck         deck
    -dimo     deck        deck        deck         deck
    -field    deck        deck        deck         deck

Sweep mechanics (2026-08-16, all 24 cells): decks
`jobs/knbo3/o3/phase3_{asan,snan}` and
`jobs/knbo3/cubic/phase3_{asan,snan}_g`, overlays refreshed from
the committed source, one log per run (`matrix_<family>.log`
in each deck) plus `matrix_summary.txt`; driver script kept in
session scratch only (trivial: loop the six families with
`$IMAGO_BIN` exported). Each `-X` job first re-runs the SCF from
the stored potential, so every gamma cell also re-ran the full
50-iteration Gamma SCF under asan and SNaN -- clean. The gamma
`-dos`/`-bond` cells report `not_converged` from the script;
that is the inherited Gamma-SCF verdict, not a failure (outputs
written, no marker). No AddressSanitizer or LeakSanitizer
report appeared anywhere in the sweep; every failure is a
deterministic logic defect with a located site.

Notes tying rows to the ledger: `-mtop` is BUG-011's
silent-death path -- candidate E below is its mechanism;
`-force` reaches BUG-016's gamma force dump, which was fixed by
inspection only (and candidate D shows the PSCF force job never
gets there); `-pacs` is where BUG-019's variable becomes live,
and needs an excited-atom deck the KNbO3 skeletons do not
provide.

### 2026-08-16 wave: five candidates from the coverage sweep

All five await programmer review; NO BUG numbers. Suggested
order: E first (it closes BUG-011), then A, B, C, D.

**CANDIDATE A -- `-loen` bare invocation builds a job the Fortran
cannot run (all four cells).**
- `imago.py` `JOB_DEFS["loen"]` gives no default bases, so bare
  `-loen` means `scf=fb, pscf=no`. The Fortran then runs an SCF
  pass and calls `loen(0)`, whose `parseInput(0)` indexes
  `tempIntArray(basisCode_PSCF)` with basisCode 0
  (`atomicTypes.f90:238`; bounds error under -fcheck, garbage
  otherwise).
- The pipeline's form `-loen -scf no` WORKS (both codes zero, so
  `commandLine.f90:178` forces PSCF to 1); confirmed live on the
  o3 asan deck: bispectrum computed, exit 0.
- `-loen -pscf fb` (any SCF basis) dies differently: "Attempting
  to allocate already allocated variable angsamplevectors" --
  `parseInput` runs twice (SCF, then loen) and the loen reader
  does not guard its allocation.
- Fix shape: default `JOB_DEFS["loen"]` to `("no", ...)`
  script-side (mirroring `build_initial_potentials.py:1069`),
  and either refuse SCF+loen or guard the re-parse Fortran-side.

**CANDIDATE B -- `-sybd` with LAT integration segfaults in
`latElectronCount` (o3, both variants).**
- The o3 decks were made with `-pscfkpint 1`. `printSYBD ->
  populateStates` takes the LAT branch (`populate.F90:198`) on
  the PATH k-points, where no mesh, tetrahedra or IBZ maps
  exist; `latElectronCount` (`populate.F90:946`) walks
  unallocated tetrahedra -> SIGSEGV.
- A band structure has no zone integral to perform; SYBD must
  not enter LAT populate (Fortran guard on the SYBD path), or
  makeinput must not write intg code 1 into the pscf file for a
  sybd job -- probably both.

**CANDIDATE C -- `-sybd` on a gamma deck whose scratch holds an
imagoG SCF file: `STOP Failed to create atom overlapCV did
SCF` (both gamma cells, general binary).**
- BUG-023 routes the job to the general binary; its SCF pass
  opens the existing `gs_scf-fb.hdf5` that imagoG wrote (real
  layout, `atomOverlapCV/real` only) and tries to create the
  complex layout's datasets -> the create fails -> STOP.
- The BUG-023 verification probe deck had NO SCF HDF5 (only the
  potential was copied in), so it created a fresh file and
  passed; the ledger's "gamma-SCF-then-general-SYBD works"
  holds only in that case. Real usage (imagoG SCF, then -sybd)
  fails until the scratch is cleared. Cross-variant HDF5
  incompatibility is the underlying issue; options: the SYBD
  job's SCF pass should not need the SCF HDF5 at all (`-scf no`
  semantics), or the file must be detected and re-created.

**CANDIDATE D -- `-force` (PSCF, job 209) is a stub: the run
"succeeds" and computes no forces (all four cells).**
- The Fortran completes (`Program Sequence Complete`, fort.2
  written) but never computes forces: `makeValenceRho`, the only
  caller of `computeForce`, is called only from `mainSCF`
  (`imago.F90:572`); the PSCF call at `imago.F90:1187` is
  commented out. `-scfforce` (109) would reach `computeForce` on
  the converged SCF iteration; 209 cannot.
- What the script harvests as force output (`fort.98`) is unit
  `97+spin` -- the debug dump in `forces.F90:821-838`, BUG-016's
  neighbourhood -- so the script's contract is with a debug
  file.
- Needs the programmer's classification: feature-incomplete
  (forces are work in progress) or dead job to remove from the
  menu until finished. Not a crash; misleading success.

**CANDIDATE E -- `-mtop` = BUG-011, no longer silent: SIGSEGV in
`buildAtomPerm` (`atomicSites.f90:378`) from `intgPSCF`
(`imago.F90:724`), o3 both variants and gamma asan.**
- The MTOP branch of `initializeKPoints` (`kpoints.f90:1221-
  1232`) builds its mesh from `numAxialMTOP_KP` and never calls
  `computeRealPointOps`, so `abcRealPointOps`/`abcRealFracTrans`
  are unallocated when `buildAtomPerm` reads them. The comment at
  `imago.F90:712-723` documents exactly this hazard for SYBD and
  guards it (`if (doSYBD_PSCF /= 1)`); MTOP needs the same
  treatment or the point-op setup on its path. This also explains
  BUG-011's "RESOLVED_KP_CLASSES varies between runs":
  `axisClass` is never computed on that path either.
- **E' (gamma, SNaN):** the trap fires earlier, `1/
  numAxialKPoints` at `kpoints.f90:965` -- the gamma deck's
  `MTOP_INPUT_DATA` carries mesh counts of 0, a division by zero
  the release build sails through into E. Whether the counts
  should ever be 0 (makeinput's default for a gamma deck) is a
  second question.

**ACCEPTED as BUG-023 and FIXED 2026-08-16 (ledger entry; the
decision moved into `init_exes`): the SYBD gamma demotion tests
the legacy executable name.**
`imago.py:1533`: before a band-structure job (job_id 108/208)
the script demotes the gamma executable to the general one via
`if exe.startswith('g'): exe = exe[1:]` -- the old OLCAO
convention of a `g`-PREFIXED gamma binary. The current name is
`imagoG`, suffix-marked, so the test never fires, and a
gamma-deck `-sybd` run would keep `imagoG` for a k-path walk
the gamma binary cannot perform. Unreachable while BUG-020
stood (nothing selected imagoG); live again now that it is
fixed. Fix shape: test the current naming (the `G` suffix) and
strip it, or better, derive both names from one place so the
next rename cannot split them again.

**RESOLVED 2026-08-16 -- NOT A DEFECT IN EITHER VARIANT; it
exposed a k-point sentinel inversion (next candidate) and a
harness trap.** The bisection ran both binaries fresh on
byte-identical inputs and compared the HDF5 intermediates
dataset by dataset. The first differing quantity was not an
integral but the k-point itself: the canonical Gamma file
(`1 1 1`, shift `0 0 0`) is resolved by `kpoints.f90`'s mesh
formula `k = (i - 1 + s)/n - 1/2` to k = (-1/2, -1/2, -1/2) --
the zone corner R -- and the complex binary sampled R with
phases of -1 on the odd cells while `imagoG` computes at Gamma
by construction. Given an explicit style-0 k-point at (0,0,0),
the complex binary reproduces `imagoG` EXACTLY: raw pair
accumulations, core-valence and core-core matrices, the
orthogonalized overlap (diagonal 1.64434995 in both), and the
entire 50-iteration SCF trajectory line for line (fresh decks
`jobs/knbo3/cubic/bisect2_g` and `bisect2_c`, current-source
release binaries; the installed binary agrees byte for byte).
So `imagoG`'s numerics are validated against the complex
arithmetic at Gamma -- the first such comparison ever -- and
the oscillation to the 50-cap is what this 5-atom cubic cell
does under Gamma-only sampling in BOTH variants (June's Phase
0d record was that behaviour, correctly computed). A stray
mid-bisection result that made the complex binary look
k-blind (Gamma == R to eight digits) was a deck-name reuse
trap, recorded in the harness list above. The original record
follows for reference. Fresh decks, byte-identical inputs
(`diff -r` verified): the complex binary's first iteration
gives total energy -102.06884276 and it converges in 13
iterations to -103.79932729 with a 0.198 au gap -- a proper
insulator. `imagoG`'s first iteration gives -100.73536664 --
already different where the variants should be numerically
equivalent -- and it then oscillates to the 50-iteration cap,
never below convergence 0.2, with the gap bouncing near zero
and the run intermittently classifying the system as metallic.
Both trajectories are deterministic: today's `imagoG` run
reproduces June's Phase 0d trajectory to every printed digit
(that June record read as a threshold artifact at the time; it
was this), which also proves the compile-verified-only gamma
changes since June did not alter imagoG's SCF numerics -- a
free A/B of that work. No prior A/B ever compared imagoG
numerics; every one ran complex-only. Evidence decks:
`jobs/knbo3/cubic/phase3_asan_g` (fresh imagoG, release
binary) and `jobs/knbo3/cubic/phase3_fresh_complex_ctrl`
(fresh complex control, run via the installed pre-fix
`imago.py` whose stale checker still selects the complex
binary on a gamma deck). NOT connected to BUG-021 (the
`gaussOverlapNP` read below): with that fix in, the SNaN
`imagoG` run reproduces the same divergent trajectory digit for
digit and completes all 50 iterations without a trap, so no
uninitialized-real read anywhere in the SCF path is the cause.

**ACCEPTED as BUG-022 and FIXED 2026-08-16 by option (a) below,
through the chain (DESIGN 3.6/3.8/3.9 -> PSEUDOCODE 4c.4 ->
`kpoints.f90`; ledger entry has the verification): the Gamma
sentinel and the mesh formula are inverted with respect to each
other.** Found by the bisection above. The record as written for
review:

*What the chain says (three statements that cannot all hold):*

- DESIGN 3.8 and the code (`kpoints.f90:980-983`): mesh points
  are `k = (i - 1 + s)/n - 1/2`, i = 1..n.
- DESIGN 3.9: "`s = 0` places a sample on the origin (Gamma-
  centred); `s = 1/2` centres the samples between nodes, so
  Gamma is absent."
- DESIGN 3.6: "on an axis with a single point the shift becomes
  that lone point's absolute coordinate", and the canonical
  Gamma file is therefore `1 1 1` with shift `0 0 0`.

3.6 and 3.9 both assume the un-offset formula `(i - 1 + s)/n`.
Under the coded formula the truth per case is:

- `n` even, `s = 0`: Gamma present (the point `(i-1)/n = 1/2`).
- `n` even, `s = 1/2`: Gamma absent.
- `n` odd, `s = 0`: Gamma ABSENT (`(i-1)/n` never equals 1/2).
- `n` odd, `s = 1/2`: Gamma PRESENT (standard Monkhorst-Pack).
- `n = 1`, `s = 0`: the lone point is at -1/2 (zone corner).
- `n = 1`, `s = 1/2`: the lone point is at 0 (Gamma).

*The single-point rule compounds it:* `resolveShift`
(`kpoints.f90:1699-1701`, implementing 3.6) zeroes the shift on
any single-point axis, so under the coded formula a style-1
request can NEVER put a lone point at the origin.

*Measured on the cubic deck (complex binary):*

- Canonical Gamma file (`1 1 1`, shift `0 0 0`) resolves to
  (-1/2, -1/2, -1/2) = R; the log prints
  `Kpoints ... -0.41313600` on all three axes (= -1/2 x 2pi/a).
- Explicit `-kpshift 0.5 0.5 0.5` on `1 1 1`: the shift is
  dropped by the single-point rule; also resolves to R.
- Hand-written style-0 list with one point at (0,0,0): the only
  way to reach Gamma; reproduces `imagoG` exactly.

*Consequences:*

- Since the July mesh rework, every complex-binary run of a
  "Gamma" deck has sampled R -- and thanks to BUG-020 that was
  EVERY Gamma deck.
- Every odd-count `s = 0` mesh has been a Gamma-free mesh while
  labelled Gamma-centred.

*Fix options (a DESIGN decision, not a code one):*

- (a) Drop the `- 1/2` from the formula (standard Gamma-centred
  convention). 3.6, 3.9, the sentinel and `check_gamma_kp` all
  become true as written. Cost: every odd-count mesh moves,
  which touches stored convergence baselines and the guidance
  database keyed on meshes.
- (b) Keep the formula; rewrite 3.6/3.9 and the sentinel around
  it (Gamma = `1 1 1` with shift 1/2), and change the single-
  point rule to FORCE the half shift instead of zeroing it.

The programmer chose (a) (baselines and the guidance database
are recomputable). Why the `-1/2` was there, so nobody restores
it: it is the classic Monkhorst-Pack prescription
`(2m - n + 1)/(2n)` = `(m + 1/2)/n - 1/2` with the built-in half
offset promoted to a parameter -- inherited from the initial
commit, and still built literally by the legacy `makeKPoints`
program (which then SUBTRACTS its shift, a third convention, and
whose single-point case places a non-Gamma point at
(0.125, 0.25, 1/3) -- the ancestor of DESIGN 3.6's old "shifted
mean-value sample" wording). No live consumer needed the offset:
tetrahedra and the MTOP map are index-based and the fold
compares modulo 1.

**ACCEPTED as BUG-021 and FIXED 2026-08-16 (see the ledger
entry; the mechanism turned out to be the reader loops, not the
kernels): `gaussOverlapNP` reads an uninitialized real in its
overlap accumulation (`integrals.F90:1909`), BOTH variants;
BLOCKED the SNaN pass.** As first recorded:
with BUG-019 fixed, fresh SNaN runs of both binaries trapped
at the same statement during `setupSCF` integral setup
(`imago.F90:337`): the accumulation `pairXBasisFn2(...) =
pairXBasisFn2(...) + oneAlphaPair(:currentlmAlphaIndex(...),
currentlmIndex(m,2)) * currentBasisFns(...)`. The accumulator
`pairXBasisFn2` is zeroed over its used range at line 1804
(every sibling routine carries the same zeroing), so the
suspect operand is `oneAlphaPair`: a fixed 16x16 local filled
per alpha pair by the integral kernels, read here at row range
`:currentlmAlphaIndex(alphaIndex(1),1)` and a column the
kernel may never have written for low-angular-momentum pairs
-- the same generated-kernel-coverage family as BUG-014's
`preFactor` gaps. Adjudication needs the kernel read: which
elements does the generated code actually write, per l1/l2
switch? Until fixed this blocks SNaN for the plain SCF and for
any fresh-deck pipeline that passes through SCF setup, on both
variants.

## Decisions log

- **Phase 2 mechanism:** multi-agent workflow (parallel subagent
  fan-out), chosen 2026-06-25.
- **Phase 2 resequenced (2026-08-12):** the open-ended nine-group
  fan-out is deferred behind cheaper evidence. Rationale: every
  ledger finding to date was SEEDED by a tool (a warning, a live
  run, a generator comparison) and settled by targeted reading
  around the seed; none came from open-ended reading. New order:
  (1) harvest the RELEASE-build `-Wuninitialized` and
  `-Wmaybe-uninitialized` tranche -- the optimizer's flow
  analysis fires only at -O2 and above, so the audit-build logs
  Phase 1 was scored on never contained this family -- and
  adjudicate those sites seed-first, definite hits before
  maybes; (2) run Phase 3 runtime instrumentation (asan,
  valgrind, SNaN+FPE) ahead of any fan-out, since it hunts the
  leak/uninit classes better than static reading and its
  findings are observed rather than inferred; (3) only then
  reassess a SHRUNKEN Phase 2: a mechanical diff of the
  preprocessed real-vs-complex texts feeding a small
  verification pass, an HDF5 handle-balance check, and the
  parallel-hazard inventory. Candidate findings are reviewed by
  the programmer BEFORE receiving BUG numbers.
- **Standing release-build warnings: manifest, not source
  separation (2026-08-13).** The `-Wmaybe-uninitialized` family
  fires only at -O2 and above, so the release build permanently
  shows the adjudicated-benign tranche groups that the audit
  build -- the tree the zero-warning doctrine was scored on --
  can never show. Leaving that standing list unguarded is
  exactly what the doctrine forbids, so route 3 (move the class
  out of the audit) is made checkable:
  `dev/tools/release_warning_manifest.tsv` records the expected
  (variant, file, class, variable) groups, and
  `dev/tools/check_release_warnings.py` diffs fresh per-variant
  from-clean build logs against it, printing loudly and exiting
  nonzero on any NEW or VANISHED warning. The manifest is
  refreshed only via `--write-manifest`, after the change to the
  expected set has been reviewed and recorded here. The tool
  lives under `dev/tools/` as campaign harness, beside the
  document that governs it, not under `src/`. A full
  gamma/non-gamma source separation was considered for the same
  problem and REJECTED: measurement showed the benign groups
  are guard-condition optimizer blind spots orthogonal to the
  variant split (a potentialUpdate.F90 home to several of them
  never branches on GAMMA, and both separated trees would keep
  the whole list), while the split would duplicate roughly 22k
  lines of shared physics across 22 GAMMA-branching files to
  dissolve only the four unused-dummy sites.
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

### BUG-015 -- ETA is never set in the PBE correlation potential
- File:     `src/imago/potentialUpdate.F90:3913` (the PBE
            correlation section; declared at `:3781`)
- Variant:  [BOTH]
- Category: UNINIT
- Severity: S2 -- silently right by luck, NaN-poisons the SCF on
            the wrong stack bits
- Status:   FIXED and A/B-VERIFIED 2026-08-12 (paired `-optc` x2
            runs byte-identical, timestamps aside); found
            2026-08-12 by the release-build `-Wuninitialized`
            harvest (the tranche above), accepted by the
            programmer 2026-08-12
- Evidence: gfortran's definite (not "maybe") flow verdict: ETA
            is read at :3913 and assigned nowhere. Both
            variants' compiles agree.
- Analysis: the section is a copy of the Burke reference CORPBE
            routine, where `ETA=1.D-12` is a regularization
            constant inside `((1+/-ZET)**2+ETA)**(-1/6)` guarding
            the spin-polarization derivative at |ZET| = 1. This
            copy is the spin-restricted (ZET = 0) specialization:
            both terms collapsed to `(1D0+ETA)`, and the ETA
            assignment was dropped. So
            `GZ = ((1+ETA)**(-1/6) - (1+ETA)**(-1/6))/3`
            subtracts two identical terms (the source spells the
            exponents `-1D0/6D0` and `-1/6D0`, but both evaluate
            to -1/6): for any garbage with 1+ETA > 0 they cancel
            EXACTLY, GZ = 0,
            which is the correct spin-restricted value -- results
            are right, by luck, almost always. If the stack frame
            ever holds bits with 1+ETA < 0 (or a NaN pattern),
            the fractional power is NaN, GZ = NaN, and it
            propagates through HZ into the correlation potential
            and the whole SCF cycle.
- Fix:      APPLIED 2026-08-12, shape chosen by the programmer:
            `GZ = 0D0` stated directly, with a comment giving the
            symmetry argument (phi(zeta) is even in zeta, so its
            derivative vanishes at zeta = 0). At the programmer's
            direction the reference form is RETAINED as
            commented-out code, and ETA's declaration is likewise
            retained commented-out beside the other correlation
            locals with a pointer to that block, so restoring the
            general formula restores both together.

### BUG-016 -- gamma force dump reads stale loop indices
- File:     `src/imago/forces.F90:850` (end-of-routine force
            matrix print, gamma branch)
- Variant:  [GAMMA]
- Category: BOUNDS
- Severity: S2 -- memory-unsafe reads and garbage output
            whenever the gamma force path runs (currently
            unexercised by any recorded run)
- Status:   FIXED 2026-08-12, compile-verified only -- no
            recorded run reaches the gamma force dump, though the
            fixed binary passed the SCF-path A/B; found
            2026-08-12 adjudicating the `-Wmaybe-uninitialized`
            tranche ('q' may be used uninitialized), accepted by
            the programmer 2026-08-12
- Evidence: the complex branch of the same print block writes
            `valeValeF(m,l,k,j,i)` using its own loop indices.
            The gamma branch loops over i (xyz), j (spin), and l
            but prints `valeValeFGamma(:,l,q,r)` -- and q, r are
            not loop variables there. They hold whatever the
            save loops above left behind: a completed
            `do q = 1, spin` / `do r = 1, 3` exits with
            q = spin+1, r = 4.
- Analysis: every print statement therefore reads
            `valeValeFGamma(:,l,spin+1,4)` -- out of bounds on
            BOTH trailing dimensions of the
            (valeDim,valeDim,spin,3) array -- and emits the same
            wrong slice for every (i,j) combination instead of
            walking the spin and direction axes. A `-fcheck=all`
            build aborts here the first time the gamma force
            dump executes. The compiler surfaced it as
            "maybe-uninitialized" only because the zero-trip
            possibility of the atom loops is the one path where
            q is never assigned at all.
- Fix:      APPLIED 2026-08-12: the print now reads
            `valeValeFGamma(:,l,j,i)`, mirroring the complex
            branch's index use.

### BUG-017 -- the GGA spin-gradient clamp misses its spin guard
- File:     `src/imago/potentialUpdate.F90:958-966` (GGA
            ray-point loop)
- Variant:  [BOTH]
- Category: UNINIT
- Severity: S4 -- a real uninitialized read, consequence-free
            today
- Status:   FIXED and A/B-VERIFIED 2026-08-12 (paired `-optc` x2
            runs byte-identical, timestamps aside); found
            2026-08-12 adjudicating the `-Wmaybe-uninitialized`
            tranche, accepted by the programmer 2026-08-12
- Evidence: the spin-difference gradient sums SXS/SYS/SZS are
            assigned under `if (spin == 2)` and consumed
            (written into exchCorrRhoSpin) under
            `if (spin == 2)`, but the smallThresh clamp between
            those sits under `if (GGA == 1)` alone. Every
            spin-restricted GGA ray point therefore compares
            three garbage values against smallThresh.
- Analysis: consequence-free as written: the clamp writes back
            only to the locals, and nothing reads them when
            spin == 1 (a NaN bit pattern simply skips the
            clamp). But it is undefined behavior standing one
            refactor away from mattering, and it is exactly the
            read the compiler flagged.
- Fix:      APPLIED 2026-08-12: the three clamps now sit under
            the same `spin == 2` guard their assignments and
            their consumers carry, with a comment saying why.

### BUG-018 -- PSCF teardown orphans the radial-function arrays
- File:     `src/imago/atomicTypes.f90` (type declaration,
            `cleanUpRadialFns`, `cleanUpAtomTypes`)
- Variant:  [BOTH] -- shared module, though observed only on the
            complex binary (the only variant with Phase 3
            runtime coverage so far)
- Category: LEAK
- Severity: S4 -- exit-time leak today, structural tomorrow
- Status:   FIXED 2026-08-14, leak-verified by an asan rerun;
            found 2026-08-13 by the Phase 3 LeakSanitizer pass,
            accepted by the programmer 2026-08-14; paired A/B
            still to run with this batch's other fixes
- Evidence: every PSCF-family run leaked 6480 bytes in 6
            allocations (3 atom types x 2 arrays on KNbO3);
            plain SCF runs were clean. `readAtomicTypes`
            allocates the pointer components `coreRadialFns` and
            `valeRadialFns` per type; the only routine freeing
            them, `cleanUpRadialFns`, is called solely from
            `setupSCF`, so the SCF path sheds them mid-setup.
            The PSCF path instead reached `cleanUpAtomTypes`,
            which deallocated the parent `atomTypes` array
            without touching the two components -- orphaning
            them.
- Analysis: exit-time and small on this deck, but it scales with
            atom-type count and becomes a live in-process leak
            for any future path that tears down and re-reads
            atom types inside one process (driver loops, library
            embedding, the parallelization work). The components
            also had UNDEFINED association status before first
            allocation -- the BUG-001 class -- which the fix
            settles as a prerequisite.
- Fix:      APPLIED 2026-08-14, shape 1 by the programmer's
            choice: `cleanUpAtomTypes` is now a complete
            teardown. Three parts: the two components are
            default-initialized to `null()` in the type
            declaration (defined association status from the
            moment `atomTypes` exists); `cleanUpRadialFns` now
            nullifies after each deallocation so the early SCF
            shed leaves a defined "gone" state; and
            `cleanUpAtomTypes` frees both components behind
            `associated()` guards before deallocating the
            parent. The rejected alternative (calling
            `cleanUpRadialFns` from `cleanUpPSCF`) would have
            left the invariant "free the components before the
            parent" enforced nowhere, for each future caller to
            rediscover. Verified: both asan variants rebuilt
            clean (`bug018_asan_{real,complex}.log`), the
            overlay copies refreshed, and the o3 `-optc` rerun
            that previously reported the leak -- run with
            `IMAGO_BIN` pointing at the asan overlay, 4m11s, so
            demonstrably the instrumented binary -- now runs
            sanitizer-silent to Program Sequence Complete.  (A
            first verification attempt silently ran the
            installed release binary and proved nothing; the
            harness trap is recorded in the Phase 3 harvest
            section.)

### BUG-019 -- readPACSControl converts an uninitialized real
- File:     `src/imago/input.f90` (`readPACSControl`)
- Variant:  [BOTH]
- Category: UNINIT
- Severity: S4 -- consequence-free when PACS is off, but a PACS
            run whose core-state list has no entry matching the
            excited QN pair would carry garbage into live state
- Status:   FIXED 2026-08-14; found 2026-08-13 by the Phase 3
            SNaN+FPE pass (which it BLOCKED at parse time on
            every deck), accepted by the programmer 2026-08-14
- Evidence: the PACS block is parsed for every run.  The trio
            `firstInitStatePACS`/`lastInitStatePACS`/
            `totalEnergyDiffPACS` is filled by the core-state
            loop only for the entry matching the excited QN pair
            from the command line; the two integers were
            zero-initialized before the loop, the real was not
            -- an asymmetry that reads as an oversight.  The
            unit conversion `totalEnergyDiffPACS =
            totalEnergyDiffPACS / hartree` then runs
            unconditionally, so with PACS off (or no matching
            entry) it divides an undefined real: SIGFPE under
            SNaN init at parse, silently-stored garbage in
            production.  The four sibling PACS reals are read
            unconditionally from the file and are always
            defined; this was the only conditionally-assigned
            one.  Same silently-right-by-luck shape as BUG-015.
- Fix:      APPLIED 2026-08-14, the programmer's choice: one
            line, `totalEnergyDiffPACS = 0.0_double` beside the
            existing zero-inits, with the comment extended to
            cover the trio and say why zero must be defined
            before the unconditional conversion.  Verified: both
            SNaN variants rebuilt clean
            (`bug019_snan_{real,complex}.log`), and the fresh o3
            run under the rebuilt binaries passes the former
            trap site -- the SNaN pass now reaches integral
            setup, where it trapped at the NEXT candidate
            (`gaussOverlapNP`, see the Phase 3 harvest section).

### BUG-020 -- check_gamma_kp parsed stale offsets; imagoG was
### never selected after the July k-point grammar change
- File:     `src/scripts/imago.py` (`check_gamma_kp`; the new
            label helper `_kp_value_tokens`)
- Variant:  selection layer (script) -- gates which variant runs
- Category: STALE-MIRROR
- Severity: S3 -- no wrong numbers, but every gamma deck since
            July silently ran the slower complex binary, the
            gamma executable lost all live coverage, and a
            style-0 file with `-scftetradiag 1` would crash the
            script on a missing token
- Status:   FIXED 2026-08-14; found 2026-08-13 (the Phase 3 SNaN
            crash line printed the invoked command for a deck
            regenerated `-kp 0 0 0`), accepted by the programmer
            2026-08-14
- Evidence: July's LAT work made `readKPoints` read
            `NUM_TETRA_DIAGONALS` and `SYMMETRIZE_LAT_PARTIALS`
            for ALL style codes before the style branch
            (`kpoints.f90:274-276`), shifting every later field
            down four lines.  `check_gamma_kp` still read
            `lines[5]`/`lines[7]`, so a style-1 gamma file was
            misread as `[4]` vs `[1,1,1]` and every style
            returned False.  Its docstring claimed the layouts
            were mirrored from `readKPoints`; they no longer
            were.
- Fix:      APPLIED 2026-08-14: values are located by LABEL via
            `_kp_value_tokens`, mirroring how the Fortran side
            reads the file (`readAndCheckLabel` then the value
            line), so an inserted or reordered field cannot
            silently shift what the script reads, and a missing
            label exits loudly naming the file and field --
            never a silent False, which is the failure mode that
            hid this defect.  Verified three ways: the new
            regression suite `src/tests/test_imago_gamma_kp.py`
            (9 tests: all three styles, the tetra-settings axis
            that broke, and the loud-failure contract; the
            companion makeinput gamma tests stay green);
            unit-level checks on the real deck files (cubic
            gamma True, o3 4x4x4 False); and direct observation
            of a live run invoking `imagoG 2 0 0 0 0 0` on the
            cubic deck (the first imagoG selection since July)
            while o3 runs still select `imago`.
- Note:     the coverage this restored immediately exposed the
            imagoG SCF divergence candidate recorded in the
            Phase 3 harvest section, and reading around the fix
            surfaced the stale SYBD gamma-demotion candidate
            (`exe.startswith('g')`) -- both awaiting review.

### BUG-021 -- the integral accumulation loops read oneAlphaPair
### columns the current alpha never defined
- File:     `src/imago/integrals.F90` (reader loops in
            `gaussOverlapOL`, `KE`, `MV`, `NP`, `EP` and
            `gaussOverlapHamPSCF`; hand-back copies in `nuclearPE`
            and `electronicPE`)
- Variant:  both
- Category: UNINIT (hot loop; the BUG-015/019 class of a
            garbage-times-zero read, but executed per alpha pair
            per lattice cell rather than once)
- Severity: S4 -- numerically inert with finite garbage (the
            multiplier is exactly zero by construction), and the
            paired A/B confirms no printed digit changes; but an
            Inf/NaN in the undefined slot would poison the sum,
            SNaN traps it, and the loops spent a row of multiplies
            per skipped state adding nothing
- Status:   FIXED 2026-08-16; found 2026-08-14 (first SNaN run to
            reach integral setup once BUG-019 was fixed), accepted
            by the programmer 2026-08-16
- Evidence: Each reader loop runs over EVERY state m of atom 2
            and reads column `currentlmIndex(m,2)` of the 16x16
            `oneAlphaPair` for the current atom 2 alpha.  But an
            alpha serves only the states up to its lm coverage
            (`currentlmAlphaIndex`, 1/4/9/16, non-increasing
            along the alpha list because the alphas that also
            serve higher l come first), and the kernels write only
            that block.  A state whose lm slot exceeds the coverage
            therefore reads a column nothing wrote for this alpha.
            The two-centre routines survive this by accident:
            `oneAlphaPair` persists across pairs and the first
            pair processed has the widest coverage, so later
            narrower pairs read STALE values from earlier pairs --
            wrong pair's numbers, but defined.  `nuclearPE` and
            `electronicPE` break even that: they accumulate into a
            FRESH 16x16 local zeroed only over the coverage block,
            then hand it back with a full-width copy
            (`oneAlphaPair(:,:) = nucPotAlphaOverlap(:,:)`), which
            transplants never-written elements over the caller's
            defined leftovers.  A temporary print of the loop
            indices at the trap caught it exactly: the SIGFPE fired
            at atom 2 alpha 17 (coverage 4, sp-only) meeting state
            9 (lm slot 5, a d state), coefficient
            `currentBasisFns` = 0.0000E+00; the 896 loop passes
            before it all had slot <= coverage.  The zero
            multiplier is guaranteed, not lucky: `renormalizeBasis`
            zeroes the whole basis array and fills each state only
            through its own orbital's alpha count.  Whether the
            states an alpha serves form a contiguous range was
            checked (they do NOT -- states run valence then core,
            each in radial-function order, so core s/p states
            follow valence d states), which is why the fix is a
            per-state test rather than a loop bound.
- Fix:      APPLIED 2026-08-16, the programmer's choice after
            weighing a loop-bound form (impossible without state
            reordering) and precomputed per-tier state lists (an
            indirect load in place of a compare, plus new arrays
            threaded through every routine -- no gain).  Two parts:
            (1) every reader loop skips a state whose lm slot lies
            above the current atom 2 alpha's coverage (`if
            (currentlmIndex(m,2) > currentlmAlphaIndex(alphaIndex
            (2),2)) cycle`) -- one well-predicted integer compare
            replacing a row of multiply-adds whose result was zero
            by construction, so the fixed loops do strictly LESS
            work; (2) `nuclearPE` and `electronicPE` hand back only
            the coverage block, so a half-defined accumulator can
            no longer leak into the caller's array.  Verified: all
            four instrumented binaries rebuilt warning-free; the
            SNaN imagoG run on the cubic gamma deck now clears the
            integrals and runs the entire 50-iteration SCF with no
            trap, matching the release imagoG trajectory
            (`phase3_asan_g`) digit-for-digit; the SNaN complex
            run on the o3 deck converges in 15 iterations, no trap,
            with the iteration record and the energy decomposition
            file line-for-line identical to a fresh release control
            (`jobs/knbo3/o3/phase3_b21_ctrl`, byte-identical
            inputs, installed unfixed binary).  Because the SNaN
            binaries carry BUG-018/019/021 together and the
            controls carry none, those two runs are the batch A/B
            for all three: numerically inert on both variants.
- Note:     the whole SCF is now SNaN-clean on both variants, so
            the SNaN column of the coverage matrix is unblocked.
            The imagoG divergence candidate is NOT this defect:
            with the fix in, imagoG reproduces its divergent
            trajectory unchanged, and a full SNaN-clean SCF rules
            out any uninitialized-real read in that path.  A
            reader loop bound over states cannot replace the test
            unless the basis states are ever regrouped by l; if
            these loops become a profiling target, per-radial-
            function iteration is the sharper lever, and belongs to
            the performance campaign, not this ledger.

### BUG-022 -- the mesh formula's -1/2 offset inverted the shift
### semantics for odd counts; the Gamma sentinel sampled R
- File:     `src/imago/kpoints.f90` (`initializeKPointMesh`, one
            statement); `dev/DESIGN.md` 3.6, 3.8, 3.9;
            `dev/PSEUDOCODE.md` 4c.4
- Variant:  complex (the gamma binary computes at Gamma by
            construction and never read the mesh); the defect is
            in the shared k-point layer
- Category: DESIGN-INCONSISTENCY (spec statements that could not
            all hold, with the code on the losing side)
- Severity: S2 -- wrong k-point sampling with no warning: every
            complex-binary run of the canonical Gamma deck since
            the July mesh rework sampled the zone corner
            (-1/2,-1/2,-1/2), and every odd-count zero-shift mesh
            was Gamma-free while documented as Gamma-centred.
            Even-count meshes were unaffected.
- Status:   FIXED 2026-08-16 (found the same day by the imagoG
            "divergence" bisection; the programmer chose the fix
            after a pros/cons review)
- Evidence: The bisection compared imagoG and the complex binary
            on byte-identical Gamma decks: the FIRST differing
            quantity was the resolved k-point the log printed
            (`-0.41313600` x3 = -1/2 x 2pi/a). DESIGN 3.8 and the
            code placed mesh points at `(m + s)/n - 1/2`; DESIGN
            3.6 and 3.9 reasoned as if they were at `(m + s)/n`
            (`s = 0` on the origin; a lone point's shift is its
            coordinate). Under the offset form the two shift
            values swap roles for odd `n`, and `resolveShift`'s
            zeroing of a single-point axis then pins that axis at
            -1/2, so no style-1 request could reach Gamma at all
            (an explicit half shift on `1 1 1` was dropped too).
            Given a hand-written style-0 point at (0,0,0) the
            complex binary reproduced imagoG exactly, integral by
            integral and over the full 50-iteration SCF -- which
            both validated imagoG's numerics and pinned the
            defect on the mesh layer. Heritage: the offset form is
            the classic Monkhorst-Pack `(2m-n+1)/(2n)` with its
            built-in half shift made a parameter, from the initial
            commit; the legacy `makeKPoints` program still builds
            that grid.
- Fix:      APPLIED 2026-08-16 through the chain. DESIGN 3.9 now
            states the convention explicitly (`k = (m + s)/n`,
            parity-independent meaning of the shift, why there is
            no offset, and the makeKPoints heritage so it is not
            read back in); 3.8's expression corrected; 3.6's
            "shifted mean-value sample" wording for `1 1 1`
            replaced by what the single-point rule actually does.
            PSEUDOCODE 4c.4 carries the same statement and the
            formula without the offset. `kpoints.f90` drops the
            `-0.5_double` (one statement) with the convention in
            its comment. Verified: (1) the canonical Gamma deck on
            the complex binary now resolves to (0,0,0) and its
            50-iteration trajectory is line-for-line identical to
            imagoG's fresh run (`jobs/knbo3/cubic/meshfix_c_gamma`
            vs `bisect2_g`); (2) the even 4x4x4 half-shift o3 deck
            folds to the same 4 IBZ points and reproduces its
            pre-change control to the 8th-9th digit, energy file
            identical at print precision -- the representatives
            moved by a reciprocal lattice vector, the physics did
            not (`jobs/knbo3/o3/meshfix_444` vs `phase3_b21_ctrl`);
            (3) an odd 3x3x3 auto-shift cubic mesh is now
            Gamma-centred (Gamma is IBZ point 1 of 4) and converges
            (`meshfix_333`). The 91 k-point/gamma tests pass
            unchanged (their expectations were already the
            un-offset semantics).
- Note:     stored convergence baselines and guidance-database
            entries keyed on ODD-count meshes describe a
            different point set than the same key now produces;
            the programmer has ruled they are recomputable (TODO
            entry). The Gamma-only SCF oscillation on the cubic
            KNbO3 cell is real behaviour of that sampling in both
            variants, not a defect.

### BUG-023 -- the SYBD gamma demotion tested the legacy binary
### name; a gamma-deck band structure ran on imagoG and died
- File:     `src/scripts/imago.py` (`init_exes`,
            `execute_program`)
- Variant:  selection layer (script) -- gates which variant runs
- Category: STALE-MIRROR (a rule written against the old
            `g`-prefixed OLCAO binary name)
- Severity: S3 -- a `-sybd`/`-scfsybd` job on any gamma deck
            failed outright (`DSYGV error code = 80` in the
            secular equation, then `makeSYBD.py` crashing on the
            missing `fort.31`); no wrong numbers, no silent
            path. Unreachable while BUG-020 stood, live once
            imagoG was selectable again.
- Status:   FIXED 2026-08-16; found 2026-08-14 reading around
            the BUG-020 fix, accepted 2026-08-16
- Evidence: `execute_program` demoted the gamma executable for
            job 108/208 with `if exe.startswith('g'): exe =
            exe[1:]`; the current name is `imagoG`, so the test
            never fired. Live probe
            (`jobs/knbo3/cubic/sybd_gamma_probe`): `-sybd` on the
            canonical Gamma deck invoked the gamma binary (the
            log shows the real-symmetric solver `DSYGV`; the
            general binary calls `ZHEGV`), which announced 297
            path k-points and died in the first diagonalization.
- Fix:      APPLIED 2026-08-16, the programmer's choice of two
            shapes: the executable is now decided ONCE, in
            `init_exes`, where the rule is stated in words -- a
            band-structure job (108/208) always runs on the
            general executable, whatever the k-point files say,
            because the path k-points are generated inside imago
            and only the complex arithmetic can evaluate them;
            the k-point files decide every other job as before.
            The name-prefix demotion in `execute_program` is
            deleted, so no second place re-decides and no
            spelling of the binary name is tested anywhere.
            Verified live: the same probe deck rerun through the
            fixed script completes on the general binary -- exit
            0, `gs_sybd-fb.plot` with all 297 path points, no
            solver error -- with the SCF potential from a Gamma
            SCF (the potential file is the same text file either
            executable writes). The nine gamma-detector tests
            stay green.
- Note:     the probe deck held NO SCF HDF5, so the general
            binary's SCF pass created a fresh one. The coverage
            sweep the same day showed the real sequence -- imagoG
            SCF leaving its real-layout `gs_scf-fb.hdf5` in
            scratch, then `-sybd` -- STOPs on
            `Failed to create atom overlapCV` (harvest candidate
            C). So this fix makes the routing right; the
            cross-variant HDF5 collision on the SYBD job's SCF
            pass is a separate, still-open defect.

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
