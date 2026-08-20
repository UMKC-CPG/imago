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
