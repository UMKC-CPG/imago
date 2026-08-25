# Performance

## Purpose

This document is the campaign ledger for making the `imago`
engine FAST, and for getting it into a state where
parallelization can begin. It is the sibling of `dev/DEBUG.md`,
which hunts correctness bugs; this one hunts time and memory.

Like that document it is a *tracking artifact*, not a sixth level
of the design chain (VISION -> ARCHITECTURE -> DESIGN ->
PSEUDOCODE -> source). The parallelization DESIGN already exists
and lives in the chain proper -- ARCHITECTURE section 6 -- and
this document does not restate it. What lives here is the
measurement apparatus, the numbers it produces, and the ledger of
hotspots those numbers expose.

**Why it is separate from DEBUG.md.** The two campaigns share a
build harness and nothing else. A bug ledger ranks findings by
consequence and reach; a performance ledger ranks them by how much
time they take and how well they will parallelize. Merging them
would force one severity scale onto two unrelated questions.

## Status

- Date opened: 2026-08-09
- **ACTIVE since 2026-08-18.** The bug campaign reached its Phase
  3 exit for the parallelization targets (`dev/DEBUG.md`,
  decisions log 2026-08-18): the SCF core and the `-optc` re-entry
  are clean under asan, SNaN and valgrind memcheck on both
  variants; the remaining candidates were punted behind this
  work by the programmer.
- **PF1 and PF2 COMPLETE; PA1 CLOSED; PA2 at the chain gate
  (2026-08-20).** Presets, decks, all six baselines, the full
  time map, the PA1 cost measurement, the ARCHITECTURE 6 and
  DESIGN 9 revisions and PSEUDOCODE 24 (the O_MPI module) are
  committed. The LARGE pair was harvested 2026-08-19/20 (jobs
  16632708/16632709; Baseline table, time map and ARCHITECTURE
  6.8 updated): at 1296 atoms the solve is 62 % of the complex
  run at ONE k-point and 37 % of the real run, and the
  once-per-run electrostatic setup (21 %/11 %) is the emerging
  Amdahl bound. A /refine pass over the parallelization chain
  (2026-08-20, `762bb9d`) propagated the PA1 decisions upward
  (VISION Goal 7, ARCH 6.5) and flagged the write-only share of
  the three-centre stage as UNMEASURED -- PA3 stamps it. The
  two-node launch test is harvested (see the toolchain bullet
  below): the launcher is `mpirun -np $SLURM_NTASKS`, `srun
  --mpi=pmix` the working alternative, `--mpi=pmi2` forbidden.
- **PA2 (the O_MPI module) CODED and ACCEPTED 2026-08-20.**
  PSEUDOCODE 24.4 was revised first (programmer's ruling):
  workers' unit 20 goes to `/dev/null` by default so a large
  run sheds ONE log file, `IMAGO_RANK_LOGS=1` switches to
  per-rank files, errors go rank-stamped to stderr. The code:
  `src/imago/mpi.F90` plus the four seam edits (imagoWrap,
  parseCommandLine's open, the fort.2 certificate behind
  `barrierMPI`, both variant CMakeLists). All four 24.7 checks
  passed (jobs 16670777/16670874; harness and throwaway decks
  under `jobs/bench/pa2*`): serial and one-rank-parallel are
  digit-identical to the bn_small baselines, the 4-rank
  replicate proof ran to the accepted scratch-HDF5 collision
  (workers die at the file lock, no fort.2, imago.py reports
  failure -- PA3's motivation, recorded), and stopMPI takes a
  2-rank job down promptly. One harness lesson worth keeping:
  `mpirun` inside a SLURM job takes its slot count from the
  ALLOCATION, so an acceptance job must request `-n` >= the
  widest `-np` it launches or mpirun refuses before any rank
  starts.
- **PA3 CODED and ACCEPTED 2026-08-21** (chain: DESIGN 9.5 as
  amended + PSEUDOCODE 25, both reviewed 2026-08-20). The
  three-centre term stage is distributed by the snake deal
  with HDF5's file lock as the write mutex and the
  root-serial remainder. All five acceptance checks pass; the
  scaling and balance numbers are in "PA3 term-stage scaling"
  below, the closure record in TODO PA3. Two latent serial
  defects were exposed and fixed on the way (a missing
  deallocate in closeSCFIntegralHDF5; packedVVDims set only
  by root-run routines), and the acceptance work itself
  produced four recorded traps (concurrent mpirun core
  binding, /usr/bin/time swallowing signals, cross-build and
  cross-node BLAS last-bit noise, relative tolerances on
  near-zero elements -- see the scaling section).
- **PA4 STAGE A (the solve boundary + k-point deal) CODED and
  ACCEPTED 2026-08-22** (chain: DESIGN 9.6 as amended +
  PSEUDOCODE 26, reviewed 2026-08-21; numbers and findings in
  "PA4-A solve deal" below, closure in TODO PA4). Root ships
  each dealt k-point's packed H and S to its owner in ROUNDS
  of one task per rank; owners solve with the serial backend;
  root writes. The first acceptance run deadlocked and bought
  the round discipline (blocking sends meet MPI's rendezvous
  path on large matrices); the thread-lever measurement
  showed ZHEGV nearly thread-proof, making ELPA the only path
  for the one-k-point solve.
- **PA4 STAGE B (the collective ELPA solve) CODED and ACCEPTED
  2026-08-22** (chain: DESIGN 9.6 stage B + PSEUDOCODE 27,
  reviewed 2026-08-22; numbers in "PA4-B collective ELPA
  solve" below, closure in TODO PA4). The one-k-point solve is
  distributed over a block-cyclic grid behind the same
  boundary, entered through one more server control code. This
  is the increment the campaign was built for: the 1296-atom
  glass now runs in 4h33m against 15h55m serial (3.5x), its
  solve 11.3x faster. Three ELPA protocol facts and the
  full-triangle unpack were measured, not read from the
  documentation, and are recorded in 27.
- **PA5a (the electrostatic setup deal) CODED and ACCEPTED
  2026-08-22** (chain: DESIGN 9.2/9.8 as re-ordered +
  PSEUDOCODE 28, reviewed 2026-08-22; numbers in "PA5a
  electrostatic setup deal" below, closure in TODO PA5a). The
  once-per-run setup -- 44 % of the post-PA4 headline run and
  the standing serial ceiling -- is dealt by site range and
  reduced onto root: ideal scaling, serial bit-identical, the
  smallest increment of the campaign. NEXT: PA4's
  valence-density tail (the third cost, 10.8 % of the headline
  run), then the atom-pair distribution of the remaining
  integral stages, per DESIGN 9.8.
- **Install no profiling tools yet** (TODO PF6). `gprof` and
  `callgrind` are already present and cover Phases P0 through P2.
  The investigation into `perf` and `gperftools` is recorded below
  so it does not have to be repeated.
- **Parallel toolchain STANDING as of 2026-08-18** (programmer's
  ruling: go straight for MPI + ELPA, the dominant use case). The
  `cpgp` conda environment (`dev/env/cpgp.yml`; `BUILD.md`
  "Parallel toolchain") carries openmpi 5.0.10, MPI-enabled HDF5
  1.14.6 (`h5pfc`), ScaLAPACK 2.2 and ELPA 2025.06 on the same
  gfortran 15.2 / OpenBLAS family as `cpg`. Verified: an `mpi_f08`
  hello runs; ELPA links and `elpa_init(20241105)` handshakes on 2
  ranks; the `gfortran-mpi` preset (`IMAGO_MPI`, `IMAGO_ELPA`)
  builds the UNCHANGED serial source against it, and those
  MPI-linked binaries reproduce the `bn_small_{c,g}` baseline
  energies to every printed digit. The two-node SLURM launch
  test (job 16634074, c[083,085], harvested 2026-08-20; outputs
  `~/data/scratch/imago/mpitest/`) settled the launcher: `mpirun
  -np $SLURM_NTASKS` spans nodes correctly and the cross-node
  ELPA handshake succeeds under it; `srun --mpi=pmix` behaves
  identically; `srun --mpi=pmi2` is a TRAP -- OpenMPI 5 dropped
  PMI-2, so it silently launches N one-rank worlds with exit
  code 0, and only a rank/size printout betrays it. The recipe
  is recorded in BUILD.md ("Parallel toolchain"). Environment
  naming convention: a trailing `p` means parallel (`cpg` ->
  `cpgp`); the Python venv is unchanged.

## Benchmark deck set (designated 2026-08-18)

Six decks under the gitignored `jobs/bench/`, one per (size,
variant) cell, all generated 2026-08-18 by the installed
`makeinput.py` (source `9f21e96`) from the historical OLCAO
skeletons with `-nofingerprint` (each species takes its database
default potential entry -- deterministic, and for the glasses it
is type-per-element, 2 types, potential dimension 32, which is
what the historical runs used; the programmer ruled out
bispectral typing for the baseline). Every baseline is a plain
SCF from the initial potential (`imago.py` with no job option).

    cell           deck               atoms  k-points          history
    small/imago    bn_small_c         8      4 4 4 (auto shift) BN cubic,
    small/imagoG   bn_small_g         8      Gamma              29.5 s SCF at
                                                                7x7x7 (2025)
    medium/imago   knbo3_333_med_c    135    2 2 2, no shift    cubic KNbO3
                                             (= 4 TRIM orbits)  3x3x3 supercell
    medium/imagoG  sio2_243_med_g     243    Gamma              19m42s SCF
                                                                (2025)
    large/imago    sio2_1296_large_c  1296   1 1 1 shift 0.5    see below
    large/imagoG   sio2_1296_large_g  1296   Gamma              2h32m for one
                                                                restart iter
                                                                + PSCF (2024)

Sources: `~/olcao/jobs/bn/cubic/olcao.skl` (sphalerite BN, space
group 216, conventional cell), `~/olcao/jobs/glass/sio2/{243,
1296}/olcao.skl` (amorphous SiO2 models), and the KNbO3
skeleton behind the `o3` decks used throughout `dev/DEBUG.md`
(which is CUBIC KNbO3, space group 221, `supercell 1 1 1`).

How medium/imago was sized (2026-08-18). The first choice, that
5-atom KNbO3 cell with 4 reduced k-points, ran its whole SCF in
4 s -- smaller than the "small" BN deck -- so the supercell line
in the skeleton was used as the size dial, keeping the cubic
symmetry and the same 2 2 2 unshifted mesh (whose four TRIM
orbits under O_h, Gamma/X/M/R, give exactly 4 k-points):
2x2x2 (40 atoms) took 1m24s, 3x3x3 (135 atoms) 25 min. The
135-atom cell is the designated medium; the 5- and 40-atom
decks (`knbo3_med_c`, `knbo3_222_med_c`) stay in `jobs/bench/`
as auxiliary size points for the scaling picture, not as cells.

Why these. Small must finish in seconds so callgrind (20-50x
slowdown) can run on it. Medium is the honest profile target:
minutes, representative of a real run. Large is the one that
hurts, which is what parallelization is for. The large-complex
cell is the SAME 1296-atom glass at ONE shifted k-point
(1 1 1, shift 0.5): a 1296-atom cell does not want a real mesh,
Gamma is its physically right sampling, but a complex-binary
large baseline is still needed, and a single non-Gamma point
gives the same problem size and the same integrals in the
complex code path -- so the pair isolates the cost of complex
arithmetic against real on the largest case. (Programmer's
ruling 2026-08-18.) `imago.py` routes the shifted point to
`imago` because it is not Gamma.

Two facts about the environment that the baseline record must
carry: (1) the head node is cgroup-limited to ONE visible CPU
and runs at load ~34, so timings taken there are noise -- every
baseline runs as a SLURM job on an exclusive node, one core;
(2) the link line says `-llapack -lblas`, but `ldd` on the built
binary -- checked on the head node AND printed by every baseline
job on its compute node -- resolves both to `libopenblas.so.0`
from the `cpg` mamba environment (`$CPG_SHARE/.../envs/cpg/lib`).
So the baselines are OpenBLAS, pinned to one thread
(`OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`), NOT reference
BLAS; a different environment could silently link a different
library, which is why the job header records what resolved. If
the secular solve dominates, the first question is still the
library (threads, or a vendor library), before any code change.

Harness: `jobs/bench/baseline.sbatch CELL...` (`REPEATS=n`,
`BENCH_BIN=` overlay, default `jobs/bench_profile_bin` holding
COPIES of the `build/gfortran-profile` binaries). It runs the
named cells sequentially on an exclusive node, resets each deck
and its scratch before every repeat so every run starts from
the initial potential, wraps `imago.py` in `/usr/bin/time -v`,
and prints host, CPU model, commit and the resolved BLAS. Two
traps it had to learn: `imagorc` takes an OPTIONAL positional
argument as the install root (the nanoHUB hook) and a sourced
file sees the caller's `$1`, so the cell list must be set aside
around the `source` or the first cell name becomes `$IMAGO_DIR`
(and `IMAGO_TEMP` becomes node-local `/tmp`); and a DANGLING
`intermediate` link (left by that very failure) makes
`imago.py`'s `os.path.exists()` test say "no link", after which
its `os.symlink` fails with EEXIST -- recorded as a candidate in
`dev/DEBUG.md`; the harness removes dead links itself.

## Measurement layers (agreed 2026-08-18)

Cheapest first; every layer runs on the designated decks and is
recorded here with commit hash, compiler, BLAS and node type.

1. **Wall clock + peak memory.** `gfortran-profile` binaries
   (release code, symbols only -- open question 4 answered: the
   generated arithmetic is the release binary's, so no
   perturbation), each deck x variant, plain SCF from the initial
   potential, under `/usr/bin/time -v` on an exclusive SLURM node
   (wall, user, max RSS). Small and medium get three repeats to
   establish the noise floor before any change is called an
   improvement. This is PF1's "single command -> comparable
   number".
2. **Coarse time map (PF2).** Free: imago already stamps `Date
   is: ... Time is: ...` at the start and end of each of the 34
   `O_TimeStamps` operations in the `.out` file. A small
   `dev/tools/timemap.py` diffs consecutive stamps into a
   per-operation table (integrals / elecStat / exchCorr / secular
   solve / potential update ...). It costs nothing, so it runs on
   the large deck too, and it directly tests ARCHITECTURE 6.5's
   inherited claims.
3. **Fine map.** `valgrind --tool=callgrind` on the small decks
   with the profile binaries: exact instruction and call counts
   per routine and per line, `callgrind_annotate` for the report;
   `gfortran-gprof` only as a cross-check for call counts and the
   call graph.
4. **Memory shape (P3, later).** Max RSS from layer 1 first;
   massif only if allocation churn becomes the question.

## Baseline (layer 1), captured 2026-08-18

Commit `9d7cdad`, preset `gfortran-profile` (`-O3 -g
-fno-omit-frame-pointer`, gfortran via `h5fc`), OpenBLAS
single-threaded (see above), node c085 of `rulisp-lab` (AMD EPYC
7713, exclusive), one core, plain SCF from the initial potential,
`imago.py` wall clock and max RSS from `/usr/bin/time -v`. Three
repeats where shown; the spread IS the noise floor.

    cell            deck                 wall (repeats)          RSS   iters
    small/imago     bn_small_c           7.7 / 7.1 / 7.3 s       47 MB  10
    small/imagoG    bn_small_g           5.8 / 5.8 / 5.8 s       32 MB   9
    medium/imago    knbo3_333_med_c      24m55 / 25m07 / 25m07  413 MB  13
    medium/imagoG   sio2_243_med_g       11m42 / 11m41 / 11m39  292 MB  11
    large/imago     sio2_1296_large_c    30h44m (1 repeat)   11.2 GB  11
    large/imagoG    sio2_1296_large_g    15h56m (1 repeat)    6.3 GB  11
    aux             knbo3_med_c (5 at.)  4.1 / 4.0 / 4.0 s       47 MB  15
    aux             knbo3_222_med_c (40) 1m24 / 1m25 / 1m24      73 MB  14

Every repeat reproduced its final total energy to all printed
digits, and `sio2_243_med_g` reproduced the 2025 OLCAO value
(-2925.62720444 Ha) exactly. So the noise floor is under 1 % on
the medium cells and a few percent on the seconds-long small
ones; a claimed improvement smaller than that is not one.

The large pair (harvested 2026-08-19: jobs 16632708/16632709,
nodes c083/c084, one repeat each -- at 16 and 31 hours a repeat
is bought only when a claim needs it) is its own cross-check
instead: the real/Gamma and complex/one-shifted-k runs of the
SAME glass agree on the final total energy to the eighth decimal
(-15614.45846254 vs -15614.45846252) with identical iteration
counts and convergence traces -- two different binaries and
arithmetic paths landing on the same physics.

## Coarse time map (layer 2), first reading 2026-08-18

`dev/tools/timemap.py` over the repeat-1 outputs. Share of the
stamped span; the stamps cover >99.8 % of every run (unstamped
time between operations is ~0.1 %), so the map is complete.

    deck              3-centre   secular   elecStat  valence  SCF
                      potential  equation  matrices  density  potential
    bn_small_c   (c)    66 %       9 %      0.1 %      5 %      8 %
    bn_small_g   (g)    77 %     0.6 %      0.2 %    0.3 %      8 %
    knbo3 5-atom (c)    45 %      14 %      0.2 %     10 %     21 %
    knbo3 40-at. (c)    37 %      37 %      0.3 %     14 %      7 %
    knbo3 135-at.(c)    40 %      46 %      0.5 %      9 %      1 %
    sio2_243     (g)    45 %      27 %     14.5 %      5 %      3 %
    sio2_1296    (g)    35 %      37 %       21 %      3 %    0.2 %
    sio2_1296    (c)    22 %      62 %       11 %      3 %    0.1 %

("3-centre potential" = the Electronic Potential Integrals stamp,
computed once per run; "secular equation" is per iteration and
per k-point; "elecStat matrices" is the once-per-run
Electrostatic Matrices setup.)

What it says, re-read with the large pair (2026-08-19):
- The two costs that matter are the ONE-TIME three-centre
  electronic-potential integrals and the PER-ITERATION secular
  solve; together they are 72-86 % of every run above the toy
  size, and the large pair sits inside that band (72 % real,
  84 % complex). Everything ARCHITECTURE 6.5 once ranked first
  -- the real-space grid work in `Make SCF Potential` (exchCorr)
  -- is 1-3 % on the medium cells and 0.2 % on the large ones:
  the ranking inversion the medium cells suggested, the large
  cells confirm.
- The solve overtakes the integrals as the cell grows, and at
  1296 atoms it HAS: 9 % -> 14 % -> 37 % -> 46 % -> 62 % along
  the complex column (the 62 % at ONE k-point), and 27 % -> 37 %
  real at Gamma. The scaling wall is no longer a projection;
  the largest deck spends more time in the solve than in
  everything else combined (complex), and the solve edges out
  the integrals even in the real binary.
- The complex-vs-real cost ratio, measured clean (same 1296-atom
  glass, one shifted k-point vs Gamma): whole run 1.93x, secular
  solve 3.27x (68708 vs 21022 s -- consistent with complex
  arithmetic's ~4x flops in ZHEGV vs DSYGV), three-centre
  integrals only 1.17x.
- A third cost emerges with size: the once-per-run
  `Electrostatic Matrices` setup is ~12050 s on BOTH variants
  (real-space work, independent of k) -- 21 % of the real run,
  11 % of the complex. Once the integrals and the solve are
  distributed, this serial 3.3 h becomes the Amdahl bound
  (about 5x max speedup on the real run if untouched). It is
  exactly the `elecStat` work upolcao already parallelized
  (ARCHITECTURE 6.7), so PA5 is last but not optional at scale.

## PA1 cost distributions (measured 2026-08-18)

The decomposition and load-balance decision for the three-centre
electronic-potential integral stage (ARCHITECTURE 6.5 first
bullet; TODO PA1) needed the per-term and per-pair cost
distributions rather than a guess. A throwaway build
(`build/pa1cost`, release+symbols) stamped every term
(`PA1TERM`: wall time of each `gaussOverlapEP` call) and every
outer-atom row of the pair loop (`PA1ROW`: pairs (i, i..N) within
a term); the stamps were removed from the source after the runs.
Decks: fresh copies `pa1_sio2_243_g` and `pa1_knbo3_333_c` under
`jobs/bench/` (outputs kept there), one exclusive-node run each;
both reproduced their baseline energies to every printed digit,
and the summed term stamps match the baseline stage stamps
(315.6 vs 316.2 s; 606.4 vs 599.4 s), so the instrumentation
cost is ~1 %.

    deck            terms  mean s  max s  max/mean   rows  row max/mean
    sio2_243 (g)      32     9.9   15.4     1.57      243      1.91
    knbo3_135 (c)     62     9.8   14.3     1.46      135      2.08

- **Both margins are mild.** Term cost varies smoothly and
  MONOTONELY with the alpha exponent (diffuse = costlier; the
  most diffuse alpha of each type tops its family), max/mean
  under 1.6. Row cost, summed over terms or within one term,
  has max/mean about 2: the negligibility cutoff makes an
  atom's cost proportional to its neighbours in range, so the
  triangular row length never shows. Nothing here needs a
  dynamic work queue: a static deal of terms, largest
  (most-diffuse) first, is "good enough" (the programmer's
  stated bar), and a cost model is available for free if wanted
  since cost tracks the exponent.
- **The surprise, and it matters more than the balance: the
  pair loop is NOT the whole stage.** Row sums account for
  73 % of the stage on the real/Gamma glass (230.6 of 315.6 s)
  but only 16 % on the complex 4-k-point cell (97.5 of
  606.4 s). The rest is per-term work OUTSIDE the pair loop --
  dominated by the core-orthogonalization and dataset write
  (`ortho(4,...)`), whose matrix products scale with k-points --
  which is why the complex share is so much larger. A
  decomposition BY TERM distributes that cost automatically; a
  decomposition by pair within a term would leave 84 % of the
  complex stage serial unless the ortho were also distributed.
  This settles the decomposition question in favour of BY TERM
  for the first implementation, with pair-level splitting (or
  in-rank threading) as the later refinement for when ranks
  exceed `potDim` -- exactly the widths involved: 32-62 terms
  here, more for richer chemistries.

## PA3 term-stage scaling (measured 2026-08-21)

The acceptance record for the distributed three-centre stage
(PSEUDOCODE 25.7; TODO PA3). Correctness ran as one concurrent
suite plus follow-ups (jobs 16677143, 16688102/450/488,
16688165); the numbers below come from the SEQUENTIAL series on
an exclusive rulisp-lab node (job 16688480, c084, one mpirun at
a time), which is the only honest way to time this: concurrent
mpirun launches each bind their ranks starting at core 0, so
simultaneous runs stack onto the same cores and corrupt both
the absolute times and the apparent balance.

    deck (terms)              np1      np2      np4      np8
    sio2_243_med_g (32)     313.5 s  158.4 s  81.6 s  43.4 s
    knbo3_333_med_c (62)    604.4 s     --   163.4 s     --

Stage-stamp wall times; ratios 1.98 / 1.94 / 1.88 per doubling
on the glass (7.2x at 8 ranks) and 3.70x at 4 ranks on the
ortho-heavy multi-k cell. The one-rank parallel stage matches
the serial baseline (313.5 vs 316-317 s; 604.4 vs 599-602 s):
the distribution machinery costs nothing when there is nothing
to distribute. Snake-deal balance: compute max/mean <= 1.05 at
every width (e.g. np4 glass: 78.0/78.5/74.5/74.6 s). The
lock-discipline overhead DESIGN 9.5 flagged as unmeasured:
lock waits 0.0-3.6 s and dataset writes ~4 % of compute at
every width -- the write-in-turn mutex is benign at these rank
counts, and the 9.7 collective form stays a later calibration.

Correctness, per the amended 25.7 criteria: serial and the
one-rank forms are BIT-exact (all four decks h5diff-identical
to the pre-PA3 serial reference; bare singleton == mpirun -np 1
bit-identical); np 2/4/8 are CLEAN at 1e-10 relative against
the same-build same-node one-rank file with eigenvectors
excluded; every run's iteration and energy traces are
digit-identical to the recorded baselines; a 4-rank run killed
mid-stage (20 of 32 terms done) restarted with exactly the 12
undone terms redealt 3/3/3/3 and finished identical to
baseline; the default run sheds no worker logs and
IMAGO_RANK_LOGS=1 sheds exactly fort.20.r0001-3.

Traps this campaign paid for, so they are paid for once:

- **Concurrent mpirun runs collide on cores.** Each mpirun
  binds its ranks starting at core 0. Nineteen concurrent runs
  put every rank 0 on the same core: a one-rank run under
  mpirun ran 5x slower than the same binary bare, and rank 0
  looked slowest in every balance table. Ratios survived;
  absolutes did not. Timing runs must be sequential (or bind
  explicitly).
- **/usr/bin/time does not forward signals.** imago.py wraps
  the engine in `time`, so TERM-ing (or KILL-ing) mpirun
  orphans the engine ranks, which keep the HDF5 file lock and
  poison any restart. Killing a parallel run from a script
  needs a process-SESSION kill (setsid + pkill -s) or must
  target the engine pids; scancel's cgroup teardown is what
  real crashes look like and needs none of this.
- **Bit-exactness has a scope.** Different OpenBLAS builds
  (cpg vs cpgp) differ in the last bit; so do different CPU
  models dispatching different kernels (~1e-19 absolute on
  near-zero elements, cross-node). Bit-level comparisons are
  meaningful only same-build AND same-node; across those
  boundaries the digit-identical printed outputs are the
  check.
- **Relative tolerances lie near zero.** h5diff -p flags
  1e-19 absolute differences on 1e-9-magnitude elements as
  1e-10 "relative disagreement". Any relative criterion needs
  the element magnitudes in view before it is believed.

## PA4-A solve deal (measured 2026-08-22)

The acceptance record for the k-point deal (PSEUDOCODE 26.6;
TODO PA4 stage A; job 16700619, sequential on exclusive c084,
in-job references from the accepted PA3 binaries so every
bit-comparison is same-build-lineage and same-node).

    knbo3_333_med_c, Secular Equation stage (14 calls):
                        np1        np2        np4
    stage wall        695.0 s    438.8 s    308.3 s
    per call           49.6 s     31.3 s     22.0 s
      = assembly       ~12 s      ~12 s      ~12 s   (root-serial)
      + solve chain    ~37.6 s    ~18.8 s    ~9.4 s  (ideal 2x
                                                      per doubling)

The deal parallelizes the SOLVES perfectly -- the workers'
reported per-solve times (9.6-10.4 s) match root's own -- and
the stage is bounded by root's serial per-k-point ASSEMBLY
(reading the ~60 packed matrices per k-point from the file:
~3 s each, ~12 s per call), exactly the designed root-side
cost: predicted bound 2.3x for this deck, measured 2.25x.
That bound is a known input to stage B: ELPA distributes the
solve itself, and the assembly share is what remains until it
is distributed or absorbed into the redistribution.

Correctness: serial bit-exact vs the PA3 binaries (bn pair and
the 25-minute medium, outputs and h5diff); bare singleton ==
mpirun -np 1 bit-identical; np 2/4 criterion-clean vs the
same-build same-node np1 file on both the 8-k-point small and
the 4-k-point medium with digit-identical energies; width-one
solve stamps unchanged from PA3 (188.0 s on the glass, PA3-era
~190 s); a mid-iteration session-kill restarted from the
checkpointed potential and converged in ONE further iteration
to the SCF tolerance (5e-7 Ha from the uninterrupted run --
checkpoint-restart of an iterative solver reconverges, it does
not replay); clean worker lifecycle, no worker files.

Two findings, both feeding stage B:

- **Blocking sends deadlock past the eager threshold.** The
  first acceptance run hung 6.7 h on np2_medc:
  dispatch-all-then-collect gave a worker a second task while
  its first reply was unsent, and once the packed matrices
  (~64 MB) crossed into MPI's rendezvous path, root's send and
  the worker's reply blocked head-on with neither receive
  posted. The small deck could never expose it (eager-path
  messages complete immediately) -- which is why the medium
  deck is in the acceptance set. The fix is the round
  discipline of PSEUDOCODE 26.3: at most one outstanding task
  per worker, so every blocking send faces a receiver already
  waiting.
- **The OpenBLAS thread lever is nearly worthless for
  ZHEGV.** 16 threads on the one-k-point glass: solve 188.0 ->
  180.1 s (~4 %). The driver is dominated by the serial QL
  iteration and BLAS-2 tridiagonalization that threads cannot
  touch. ELPA's two-stage algorithm exists precisely because
  of this; it is now measured to be the ONLY path to a faster
  single-k-point solve.

## PA4-B collective ELPA solve (measured 2026-08-22)

The acceptance record for the distributed single solve
(PSEUDOCODE 27.5; TODO PA4 stage B; jobs 16701920 and 16701921
on exclusive c084, 16702210 on c049; in-job references from the
accepted PA3 and PA4-A binaries).

    sio2_243_med_g, Secular Equation stage (12 calls):
                        np1        np2        np4        np8
    stage wall        188.0 s     74.3 s     62.5 s     57.2 s
    per call           15.7 s      6.2 s      5.2 s      4.8 s
                      LAPACK      ELPA       ELPA       ELPA

At valence dimension 2349 the collective arm beats the serial
one by 2.5x at two ranks and then flattens: the matrix is too
small for eight ranks to help, which is the expected shape and
the reason the headline deck exists. The 135-atom KNbO3 3x3x3
supercell at one k-point (the COMPLEX arm, job 16702210) gives
the same picture with more room: whole run 527.3 -> 197.4 s and
the solve 197.9 -> 65.6 s at np8 (3.0x), energies identical.

The headline -- `sio2_1296_large_g`, valence dimension 12528
(numStates 5184, potDim 32), at np8 against the serial baseline:

    stage                     serial       np8    speedup
    whole run               15h55m      4h33m       3.5x
    Secular Equation        21022 s     1863 s     11.3x
    3-centre term stage     20316 s     3342 s      6.1x
    Electrostatic Matrices  12045 s     7237 s      (serial;
                                                  PA5a's target)

Energies IDENTICAL to the baseline. The solve's 11.3x is the
number to read: ELPA's two-stage algorithm scales far better at
12528 than at 2349, so the solve fell from 37 % of the run to
11 % and the once-per-run electrostatic setup became the new
leader at 44 %. The term stage's per-rank table is the best
balance the campaign has recorded: compute 3258-3328 s across
all eight ranks (max/mean 1.007), lock waits <= 14.5 s, writes
~13 s.

Correctness, per the amended 27.5 criteria: serial bit-exact
against the PA3 binaries; the 4-k-point medium np4 deal stamp
unchanged at 308.6 s (the deal path is untouched); the glass at
np 2/4/8 clean against the same-build same-node np1 file at
1e-9 ABSOLUTE whole-file and 1e-12 absolute on the eigenvalues,
eigenvectors excluded, with energies digit-identical. The
absolute form is required because this comparison crosses
ALGORITHMS -- the 1e-10 relative criterion of PA3 mis-fires on
near-zero density coefficients whose absolute differences sit
at the eigensolver floor, the trap PA3 already recorded. On the
complex supercell the energies are identical while the
ITERATION traces differ: floor-level differences move the
intermediate iterations, and the converged answer does not
care.

Three ELPA protocol facts the coding measured, each after a
failed run, all folded into PSEUDOCODE 27:

- **The handle setup is collective.** `elpaHandle%setup()`
  splits the row and column communicators, so every worker must
  be woken into its own `ensureELPA` BEFORE root enters its
  own; the first acceptance run hung with root inside the
  collective and the workers waiting for a control message.
- **The generalized path requires a BLACS context.** The
  standard eigenproblem needs none, but
  `generalized_eigenvectors` runs ScaLAPACK operations for its
  Cholesky transformation and refuses without one ("BLACS
  context has not been set beforehand"). The grid is created
  column-major so the BLACS view and the scattered data agree.
- **A fixed block edge can starve a process column.** With
  nblk = 64 a matrix smaller than 64 x npcols left one column
  of the grid owning nothing, and the generalized solve
  returned garbage on the degenerate layout. The edge is now
  adaptive and the arm policy additionally requires
  valeDim >= mpiSize. Production sizes (valeDim 1620-12528 on
  the benchmark decks) were never affected; a 60-row acceptance
  deck exposed it.

And one data-shape fact: ELPA consumes the FULL Hermitian
matrix, while LAPACK's solvers read only the upper triangle --
so the serial unpack, which fills only that half, fed ELPA
zeros for the rest and produced garbage eigenvectors and a
POSITIVE total energy. The collective arm unpacks both
triangles.

## PA5a electrostatic setup deal (measured 2026-08-22)

The acceptance record for the dealt electrostatic setup
(PSEUDOCODE 28.3; TODO PA5a; job 16707849, sequential on
exclusive c084, in-job references from the accepted PA4
binaries).

    sio2_243_med_g, Electrostatic Matrices stage (1 call):
                        np1        np2        np4        np8
    stage wall         58.6 s     29.4 s     14.8 s      7.8 s

Ideal halving at every doubling (7.5x at eight ranks), with the
per-rank loop times inside 3 % of one another at every width --
the near-uniform site costs DESIGN 9.2 assumed, confirmed
rather than asserted. Applying that ratio to the headline run's
7237-second setup projects the 1296-atom np8 run at about 2.9 h
-- 5.5x over serial -- which is the point of the re-ordering
recorded in DESIGN 9.8.

Correctness: the serial build is BIT-IDENTICAL to the pre-PA5a
binaries on both the small pair and the 243-atom glass, so the
two range-safety corrections (refresh on TYPE CHANGE rather
than on the `firstPotType` flag; alpha indices from
`cumulAlphaSum` rather than a sequential carry) are exactly
equivalent to what they replaced. At np 2/4/8 the files are
clean against the same-build same-node np1 file at 1e-6
ABSOLUTE with eigenvectors excluded and the energies are
identical. That criterion is measured, not chosen for comfort:
the reduce reorders the accumulator sums, and the potential
fit's linear solve carries the floor into `potCoeffs` as about
2e-9 absolute on order-40 values -- 5e-11 relative, physically
nothing, but past any near-zero-safe relative form. Every
physical magnitude in the file is below 1e3, so 1e-6 absolute
is at worst 1e-9 relative where it matters. A run killed
mid-SCF restarted with the finished setup skipped (37 skip
lines, no dispatch, no orphans) and reconverged to the baseline
total energy.

Two observations worth keeping:

- **The cost is the reciprocal-space sub-stage.** On the glass,
  `neutralAndNuclearQPot` takes 0.099 s and `residualQ` takes
  99.97 s. Both are dealt, but only the second matters, and it
  is the one that scales with the reciprocal-cell count.
- **The stage is about 1.7x faster in the MPI build.** Same
  node, same job: 100.1 s under the serial `cpg` build against
  58.6 s under the `cpgp` MPI build at one rank; the same ratio
  separates the 1296-atom baseline's 12045 s from the headline
  run's 7237 s. ATTRIBUTED 2026-08-23 -- see "The two
  toolchains are not equally fast" below. PA5a's own np-series
  is same-build throughout and is unaffected by it.

## The two toolchains are not equally fast (attributed 2026-08-23)

The `cpgp` (parallel) build emits calls to glibc's VECTOR math
routines and the `cpg` (serial) build does not. `imagoG` from
`jobs/pa5a_mpi_bin` carries undefined symbols `_ZGVbN2v_exp`,
`_ZGVbN2v_sin`, `_ZGVbN2v_cos` and `_ZGVbN2vv_pow` and a NEEDED
entry for `libmvec.so.1`; the same program from
`jobs/pa5a_release_bin` carries none of them and links only
scalar `libm`. A probe isolating one exponential loop
(`vecmath_probe.f90`, the shape `residualQ` runs) compiled at
`-O3` by each environment's own gfortran, on one node:

    environment   vector symbols        exp loop     checksum
    cpg           none                   0.257 s   999999.49793899385
    cpgp          _ZGVbN2v_exp           0.083 s   999999.49793899385

3.1x on the loop, with the summed result identical to every
printed digit. The cause is not a flag imago sets: both builds
compile at `-O3 -fimplicit-none -Wall` (the MPI preset adds
only `-g -fno-omit-frame-pointer` from IMAGO_PROFILE, plus the
ELPA includes), and neither HDF5 wrapper injects optimization
flags. It is the COMPILER PACKAGE. Both environments carry
gfortran 15.2.0, but from different conda-forge builds --
`15.2.0-7` in `cpg` against `15.2.0-20` in `cpgp` -- configured
against different sysroots, and only the newer one's target
glibc lets GCC's vectorizer emit libmvec calls. Fortran gets
this without `-ffast-math` because gfortran already implies
`-fno-math-errno`.

What follows, in order of consequence:

1. **The serial build is leaving a free speedup on the floor.**
   Every exp/erf/pow-heavy stage -- the electrostatic setup
   above all, and plausibly the exchange-correlation mesh and
   the Gaussian evaluation inside the integrals -- runs slower
   in `cpg` for no reason but a package pin. Upgrading `cpg`'s
   compiler (or pinning both environments to one build) is a
   no-code change with a measured payoff.
2. **Cross-build ratios in this file are contaminated, and the
   headline is one of them.** The 3.5x whole-run figure
   compares a `cpg` serial baseline against a `cpgp` np8 run.
   The electrostatic setup alone accounts for 4808 s of that
   gap from the toolchain rather than from parallelism, which
   puts the honest comparison at about 14h35m vs 4h33m = 3.2x.
   Other stages may add to that correction. The stage-level
   speedups measured WITHIN one build -- the solve's 11.3x, the
   term stage's 6.1x, PA5a's 7.5x -- are untouched, because
   both sides of each of those ratios came from the same
   binary.
3. **Cross-build bit-comparison now has a second named cause.**
   PA3 recorded OpenBLAS build differences; vector versus
   scalar `exp`/`sin`/`cos`/`pow` is the other. The rule is
   unchanged and now better founded: bit-level comparisons are
   meaningful only same-build AND same-node.

What is NOT yet done, and is TODO PF7: rebuild the serial
preset against a toolchain that emits libmvec, confirm the
243-atom glass setup lands near 58 s rather than 100 s, and
re-take the serial baselines that every parallel ratio in this
file is quoted against.

## Integral storage layout and read cost (measured 2026-08-23)

Read from the files themselves with `h5ls -v`, and from the
dataset-creation code in `hdf5SCFIntg.F90`, while settling
whether parallel reads could serve the valence density step
(DESIGN 9.7). Two facts about how the integrals are stored
decide more than any timing does.

**Every packed matrix is ONE compressed chunk.** The datasets
are created chunked with DEFLATE level 1, and the chunk is the
whole dataset unless `packedVVDims(1) * packedVVDims(2)` passes
250 million elements, which engages only above valence
dimension 22,360 real / 15,811 complex. The largest benchmark
cell is 12,528. Confirmed on disk rather than inferred from the
code -- every potential-term dataset of the 1296-atom run
reports `Chunks: {78481656, 1} 627853248 bytes`, one chunk
holding the entire matrix.

**Compression is strong and very uneven.** Per dataset on the
1296-atom real run, all of them 627,853,248 logical bytes:

    dataset                    on disk    ratio
    atomPotOverlap term 1     191.6 MB     3.3x
    atomPotOverlap term 16     11.6 MB    54.1x
    atomPotOverlap term 32     13.6 MB    46.3x
    atomOverlap               209.4 MB     3.0x
    atomNPOverlap             220.1 MB     2.9x
    atomKEOverlap             242.8 MB     2.6x
    atomMVOverlap                   0      (unwritten; rel=0)

The spread is the negligibility mask showing up in the file: a
diffuse term reaches most pairs and barely compresses, a tight
one is mostly exact zeros and compresses fifty-fold. Totals for
the 35 datasets the valence density reads per k-point: 21.97 GB
logical, 3.85 GB on disk, 5.7x overall. For the 64 datasets the
secular assembly reads per k-point on the 4-k-point medium
cell: 1.345 GB logical, 0.185 GB on disk, 7.3x.

**Per-iteration read volume.** The assembly and the valence
density read essentially the same set, so the SCF loop inflates
about 44 GB per iteration on the large cell -- data that has
not changed since the integral stage produced it.

**The read rate, and what it does NOT establish.** PA4-A
measured root's assembly at about 12 s per call on the medium
cell; that call covers all FOUR k-points, so it is about 3 s
per k-point -- roughly 448 MB/s of uncompressed output, or
62 MB/s off the device. Both figures sit inside the plausible
band for their own mechanism (DEFLATE level 1 inflates at a few
hundred MB/s per core; 62 MB/s is slow for this filesystem but
not impossible across 64 separate reads), so this measurement
does NOT decide whether the cost is inflation or storage. An
earlier note in DESIGN 9.7 claimed inflation confidently on an
arithmetic slip -- it read the 12 s as one k-point's worth and
derived 26 MB/s, which nothing here would explain. Corrected
here and there. TODO PF8 settles it directly.

## Valence density split and thread sweep (measured 2026-08-23)

PSEUDOCODE 30's instrument, run to decide the order of DESIGN
9.6's three one-k-point candidates. Four regions are stamped
inside `makeValenceRho`: R reads eigenvectors, A accumulates the
density by rank-1 update, I reads the integral datasets, M is
the contract arithmetic. Binary `jobs/p30_check/bin` at commit
e353f74; one exclusive node (c085, AMD EPYC 7713, 128 cores);
the SAME binary at every thread setting, so the toolchain
difference recorded above cancels out of every ratio here.

**The instrument checks out.** Check 2: the four regions sum to
the stage's own start-to-end stamp within 9 ms on every one of
the ten readings, always slightly under and never over, which
is the `date_and_time` resolution floor rather than unaccounted
work. Check 4: the 1296-atom cell reports 3456 rank-1 updates,
and 432 SiO2 units at 16 valence electrons each over 2 is
exactly 3456 -- written down as a prediction before the job
ran. Check 5: A moves by 2.5x across the sweep while I and M
sit still, so the thread setting is reaching the library and
the sweep is not measuring its own plumbing.

**Where the cost sits depends entirely on cell size.** One
k-point, gamma-point real build, seconds per call:

    deck              valeDim  updates  triangle       A        I    A share
    bn_small_g             64       16     16 KB   0.0000   0.0025      1.5%
    sio2_243_probe       2349      648     22 MB   1.029    3.247      23.1%
    sio2_1296_probe     12528     3456    628 MB 127.96    30.51       79.4%

The change-over sits between the middle and the large cell, and
it sits exactly where the accumulator's stored triangle stops
fitting in a last-level cache. Below that size the integral
read dominates and the accumulation is nearly free; above it
the accumulation is the stage. A design chosen on either of the
two smaller decks would have parallelized the wrong half, which
is what PSEUDOCODE 29.7's deck constraint was written to
prevent. (The 243-atom probe's thread setting was not recorded,
so it belongs to this ladder but not to the sweep below.)

**The thread sweep on the large cell.** Two calls per setting:

    threads  A call 1  A call 2  A mean  speedup  eff    GB/s       I       M
    1         128.00    127.91   127.96    1.00  100%    33.9   30.51   2.782
    2         103.61    103.84   103.72    1.23   62%    41.8   30.67   2.782
    4          69.05     75.77    72.41    1.77   44%    59.9   30.50   2.781
    8          46.01     55.63    50.82    2.52   31%    85.4   30.60   2.781

**DESIGN 9.6's traffic arithmetic held.** The stored triangle is
12528 * 12529 / 2 = 78,481,656 doubles = 627,853,248 bytes,
which is the same number `h5ls` reports for a packed matrix
chunk in the section above. Each `dsyr` reads and writes all of
it, so 3456 updates move 4.34 TB against 542 GFLOP of
arithmetic: 0.125 flop per byte, the eighth that 9.6 predicted
from the array shapes. The serial rate follows -- 4.24 GFLOP/s
on a core with roughly fifty available. The arithmetic unit is
idle while the memory system works flat out.

**Threads are not a substitute for fixing that.** Eight of them
return 2.52x at 31 % efficiency, and the run-to-run spread
grows from 0.07 % at one thread to 19 % at eight. The sweep
sets no `OMP_PROC_BIND` or `OMP_PLACES`, so on a 128-core part
the four or eight threads may land inside one core complex on
one call and spread across several on the next, changing how
many memory controllers serve the traffic. Treat the t4 and t8
points as carrying +/- 10 to 20 %; t1 and t2 are tight. Note
also that the curve has NOT flattened at eight threads -- each
doubling still returns 1.3 to 1.6x -- so this sweep did not
find the bandwidth ceiling, only established that it is above
85 GB/s.

**The multi-k-point decks behave oppositely.** Seconds per
call, complex build:

    deck              nk  valeDim  updates      R       A       I    I share
    bn_small_c         8       64      128  0.0031  0.0007  0.0478     89.9%
    knbo3_333_probe    4     1620     1944  0.509   1.040   8.439      79.0%

Region R is exactly zero on every one-k-point reading and 4.8 %
of the stage at four k-points: the eigenvector re-read that
DESIGN 9.6 worried about does not exist in the case 9.6 is
about.

**What this settles for DESIGN 9.6, and what it changes.** The
ordering is confirmed, but the second candidate is worth more
than 9.6 gave it.

1. **Candidate (i), the rank-k recast, goes first** -- and by a
   wider margin than 9.6 argued. It is the only move that
   attacks the 0.125 flop per byte itself. Blocking the update
   over columns of width nb cuts the traffic to roughly
   (3456 / nb) * 1.26 GB, which at nb = 64 is 68 GB instead of
   4.34 TB. The measured alternative -- eight threads on the
   existing rank-1 loop -- buys 2.52x with 19 % noise and costs
   eight cores that MPI ranks would otherwise use.
2. **Candidate (ii), the dataset deal, is promoted to a near
   sibling of (i) rather than a distant second.** I is
   invariant: 30.5 s at every thread setting, because nothing
   in it is threaded. It is 19 % of the stage today, but once
   (i) lands and A falls toward the tens of seconds, the
   integral read becomes the largest region in the stage. 9.6
   treats (ii) as a follow-on; the measurement says it is the
   next bottleneck, immediately.
3. **Candidate (iii), the state axis, recedes further.** After
   (i) there is little left for it to divide, and it remains
   the only one of the three that would break
   parallel-equals-serial.

**What is NOT established here.** The bandwidth ceiling, per
above. Whether I's 30.5 s is inflation or storage -- that is
PF8, unchanged. And the recast's actual speed: 68 GB of traffic
predicts a compute-bound `dsyrk`, but no `dsyrk` has been timed
on this shape, and that measurement belongs in the section that
specifies candidate (i).

## The situation, as of 2026-08-09

Three things are true at once, and the campaign has to be planned
around all three.

**1. The debugging environment exists and is unused.** Phase 0 of
`dev/DEBUG.md` is complete: `CMakePresets.json` gives
`gfortran-release`, `gfortran-debug`, `gfortran-audit` and
`gfortran-asan`; `BUILD.md` documents them; an instrumented binary
has been shown to run a real deck clean. But Phases 1-3 never ran
and the findings ledger holds ZERO entries. The instrumented build
trees date from 2026-06-26 and predate months of source change.

**2. The performance environment does not exist at all.** There is
no profiling build, no benchmark deck set, no timing baseline, and
no ledger of where the time goes. The only instrumentation in the
engine is `O_TimeStamps` (`src/imago/timeStamps.f90`, 147 lines),
which writes wall-clock start and end STRINGS for 34 named
operations into the log. That tells a reader when a phase began.
It accumulates nothing, counts no calls, and cannot say which
routine inside a phase costs anything. It is a progress indicator,
not a profile.

**3. The parallelization design is largely written already.**
ARCHITECTURE section 6 covers it in seven subsections, and the
work is already decomposed in TODO as A3 through A6. This campaign
does not need to invent that plan; it needs to produce the
measurements that say which part of it to do first.

## What ARCHITECTURE 6 already decides

Read it before proposing anything -- the design questions are
settled and only the ordering is open.

```
  6.1  configurable precision (a `wp` kind)      -> TODO A3
  6.2  inner-loop vectorization of alpha pairs   -> TODO A4
  6.3  GPU offload of the compute phase          -> TODO A6
  6.4  reference: a prior vectorized integrals
         implementation OUTSIDE this repo, at
         /home/rulisp/lewis/CPG/cpg-repo/v34/src/olcao/
  6.5  distributed memory (MPI), with a seam map
         inherited from a sibling OLCAO branch
  6.6  eigensolver backend abstraction
  6.7  reference: the upolcao MPI exploration
```

Two claims in 6.5 are worth carrying forward because they shape
the first move:

- The real-space grid loops in `elecStat` and `exchCorr` are
  called out as *embarrassingly parallel* and "the cheap, proven
  win and the right first target."
- The secular solve `H c = e S c` is named as "the actual scaling
  wall."

Both are assertions inherited from a sibling code base. **Neither
has been measured on imago.** Establishing whether they are true
here is exactly what this campaign is for: if the grid loops turn
out to be five percent of runtime, "the right first target" is the
wrong first target.

## Tool inventory on this machine (measured 2026-08-09)

```
  valgrind             /usr/bin/valgrind   3.22.0
  callgrind_annotate   /usr/bin/callgrind_annotate
  gprof                /usr/bin/gprof      2.30
  perf                 NOT on PATH
  massif-visualizer    NOT on PATH
  hpcrun               NOT on PATH
  nvidia-smi           NOT on PATH
```

### Can `perf` be obtained? Investigated 2026-08-09; answer: not
### without root, and two plausible routes are traps

**Policy is NOT the obstacle.** `/proc/sys/kernel/perf_event_paranoid`
is **2** on this machine, which permits unprivileged USER-SPACE
profiling and blocks only kernel profiling. For a compute-bound
Fortran code user space is all anyone wants, so perf would work
here the moment the binary existed.

**Route 1, conda or pip: a trap.** There IS a package named `perf`
on conda-forge and on PyPI. It is not this tool. Versions 1.5-1.6
with `py27` / `py36` / `py37` build strings give it away: it is
Victor Stinner's PYTHON microbenchmark library, since renamed
`pyperf`. It cannot profile Fortran and installing it would waste
an afternoon before anyone noticed.

**Route 2, the real thing: needs root.** Linux perf ships in the
`perf` RPM, built against kernel headers and versioned to the
running kernel (here `4.18.0-553.141.2.el8_10`). `rpm -q perf`
reports it is not installed. Being kernel-tied is exactly why
conda-forge carries no such package -- there is no
kernel-agnostic build to ship. So obtaining perf means asking a
system administrator for a stock RHEL package. That is a
reasonable request, and `paranoid = 2` means it would work on
arrival; it simply is not self-service.

### The gperftools module: sampler yes, reader no

`module load gperftools/v2.17.2` resolves to a SPACK prefix that
contains

```
  lib/libprofiler.so      the CPU sampler -- present
  bin/                    EMPTY. there is no pprof.
```

Modern gperftools dropped its bundled Perl `pprof` and expects
Google's Go `pprof` instead; the install even ships
`share/doc/gperftools/pprof_integration.adoc` saying so. So it
will sample happily and write a binary profile that nothing here
can read. Conda-forge has `pprofile` (an unrelated PYTHON line
profiler) but not Go pprof, and there is no `go` on the system,
though `go-cgo` is installable from conda-forge if someone wanted
to build it.

**Keep gperftools in mind anyway**, for one specific reason: it
samples through `SIGPROF` in user space and needs no kernel perf
events at all. On a machine where `perf_event_paranoid` were 3 it
would work where perf could not. That is insurance against a
restriction this machine does not currently impose.

### What this means for the plan

**Install nothing yet.** `gprof` and `callgrind` are already
present and cover Phases P0 through P2 completely -- the coarse
map and the exact per-routine costs on a small deck -- without a
single new package. perf's advantage is low overhead on large,
long runs; gperftools' advantage is surviving stricter
permissions than this machine imposes. Both answer problems not
yet encountered.

The moment to ask an administrator for the perf RPM is when a
callgrind run on a production deck proves too slow to be
practical. Then the request carries a measurement behind it
instead of a preference.

**The absence of `perf` is the one that constrains the plan.**
Without it there is no low-overhead sampling profiler, which
leaves:

- *callgrind* (a valgrind tool) -- exact call counts and
  instruction costs, but runs the program roughly 20-50x slower.
  Excellent for "where does the time go" on a SMALL deck; useless
  on a production one.
- *gprof* (`-pg`) -- cheap, but its sampling is coarse, it
  attributes badly through library calls, and it needs its own
  build.
- *self-instrumentation* -- extending `O_TimeStamps` into real
  accumulating timers. Costs almost nothing at runtime, survives
  into production runs, and is the only option that will still
  work once the code is parallel and running on many ranks.

That last point is why self-instrumentation is likely the
backbone rather than an afterthought: callgrind and gprof both
become awkward the moment MPI enters, and this campaign exists to
lead into MPI.

## Open questions to settle before building anything

These are decisions, not tasks. They are written here so the next
session can put them to the programmer rather than guess.

1. ~~**Bugs first, or performance first?**~~ **SETTLED
   2026-08-09: bugs first.** The programmer confirmed the
   `dev/DEBUG.md` doctrine holds -- squash the bugs in the serial
   code *before* a parallelized version is developed, because bugs
   are far harder to find and reproduce in a parallel environment.

   So `dev/DEBUG.md` Phases 1, 2 and 3 are the active work and
   this campaign waits behind them. Two consequences worth being
   explicit about, because they are what the decision actually
   costs and buys:

   - **Nothing here is blocked from starting**, only from acting.
     Profiling finds no bugs and changes no code, so PF1 (baseline
     and benchmark decks) and PF4 (parallelization-hazard sweep)
     can proceed alongside the bug hunt if there is appetite. What
     must wait is RESTRUCTURING -- A3 through A6 -- since that
     rewrites the very code the bug campaign is reading.
   - **The benchmark decks are shared infrastructure.** Phase 3 of
     the bug campaign needs representative inputs to run under
     valgrind, and PF1 needs representative inputs to profile.
     They are largely the same decks, so whichever campaign
     designates them first should do it for both.

2. ~~**What is the benchmark deck set?**~~ **SETTLED 2026-08-18;
   see "Benchmark deck set" above.** Optimization without a
   fixed, versioned set of representative inputs measures noise.
   Needs at least: a small fast deck (callgrind-able, seconds), a
   medium one (the honest profile target), and a large one that
   hurts. The gitignored `jobs/` tree holds many candidates but
   none is designated. Note also that both variants must be
   covered -- `imagoG` (real, gamma) and `imago` (complex,
   multi-k) are different programs and will have different
   profiles.

3. **What is "faster" measured against?** A baseline has to be
   captured and stored before any change, or improvement claims
   are unfalsifiable. This project's own recent history is the
   argument: several optical changes looked obviously right and
   were wrong, and only a measurement that could fail caught them.

4. **Does the profiling build perturb what it measures?** `-pg`
   and heavy instrumentation change inlining and layout. The
   baseline must be captured from the SAME build configuration
   that later comparisons use.

## Proposed phases (draft, not agreed)

Written cheapest-first, mirroring the DEBUG.md structure. Nothing
here is committed to.

**Phase P0 -- Measurement harness.** A profiling build (a
`gfortran-profile` preset carrying `-pg` and/or a callgrind-clean
`-O2 -g`), a designated benchmark deck set, and a recorded
baseline for both variants. Done when a single command produces a
comparable number.

**Phase P1 -- Coarse map.** Where does the time actually go, at
the level of the 34 operations `O_TimeStamps` already names?
Cheapest possible answer, and enough to test ARCHITECTURE 6.5's
inherited claims about `elecStat`, `exchCorr` and the secular
solve.

**Phase P2 -- Fine map.** callgrind on the small deck for exact
call counts and per-routine cost; identify the actual hot loops
rather than the assumed ones.

**Phase P3 -- Memory behaviour.** Peak footprint and allocation
churn, which decide what can be distributed and what must be. This
is where the campaign touches `dev/DEBUG.md`'s leak hunt, and the
two should share findings.

**Phase P4 -- Parallelization readiness.** The PARALLEL-HAZARD tag
that `dev/DEBUG.md` already defines -- mutable `SAVE`/module
state, non-reentrant routines, race-prone I/O. These block
threading regardless of where the time goes, and finding them is
static work that can proceed in parallel with everything above.

## Findings ledger

### HOT-001 -- 121 hidden array copies, from the compiler, for free
- Site:      121 distinct sites, engine-wide. Worst offenders:
             ```
               48  imago/integrals.F90
               43  imago/auxiliary/imagoKKc.f90
               37  imago/potentialUpdate.F90
               22  imago/integrals3Terms.F90
               20  imago/secularEqn.F90
               12  imago/forces.F90
                9  imago/field.F90
                8  imago/coreCharge.f90
             ```
             (187 warnings over 121 sites, since each site is
             reported once per variant that compiles it.)
- Variant:   [BOTH]
- Deck:      none -- this is a COMPILE-time diagnostic, obtained
             from the Phase 1 warning sweep of `dev/DEBUG.md` at
             no runtime cost whatever
- Cost:      UNKNOWN and deliberately not guessed. An array
             temporary is a full copy of an array the compiler had
             to materialize because it could not prove an operation
             was safe in place. Whether that matters depends
             entirely on the array's size and on how often the line
             executes -- a temporary in a setup routine costs
             nothing, one inside the integral loops could be
             enormous. Ranking these needs PF1 and PF2.
- Scaling:   unknown, per above
- Parallel:  unknown
- Status:    open
- Evidence:  `-Warray-temporaries` under the `gfortran-audit`
             build, 2026-08-09. Moved out of the correctness audit
             and behind the new `IMAGO_WARN_PERF` option in the
             same change, because 187 of these against 4 genuine
             correctness warnings made the latter unfindable.

**Why this is not yet actionable, and why it is still worth
having.** It says WHERE copies happen but nothing about what they
cost, and the distribution above is suggestive rather than
conclusive: `integrals.F90` leading the count is consistent with
ARCHITECTURE 6.2's plan to restructure exactly those alpha-pair
loops, but a count is not a cost. Cross this list against the PF2
time map and the intersection -- a hot routine that also copies --
is where the cheap wins are.

Note `imagoKKc.f90` at 43. That is post-processing, not the SCF,
so it will not show up in a time map of the engine at all; it is
its own small optimization target, separate from the parallel
work.

Entry schema, deliberately different from the bug ledger's --
these rank by cost and by how well the site will parallelize:

```
### HOT-NNN -- <short title>
- Site:      <path>:<line>  (procedure/module)
- Variant:   [BOTH] | [GAMMA] | [COMPLEX]
- Deck:      which benchmark input, which build
- Cost:      percent of run, absolute time, call count
- Scaling:   how the cost grows (atoms, k-points, basis, SCF
             iterations)
- Parallel:  embarrassingly | reducible | serialized | unknown
- Status:    open | confirmed | addressed | wontfix
- Evidence:  the measurement, not the expectation
```
