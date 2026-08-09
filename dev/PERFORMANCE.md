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
- **Nothing built yet.** This document is an orientation written
  before the work starts, so that the situation does not have to
  be re-derived.
- **This campaign is WAITING.** Sequencing settled 2026-08-09:
  the bug campaign in `dev/DEBUG.md` goes first (see open question
  1 below, and TODO PF5). Its Phases 1, 2 and 3 are the active
  work.
- **Install no profiling tools yet** (TODO PF6). `gprof` and
  `callgrind` are already present and cover Phases P0 through P2.
  The investigation into `perf` and `gperftools` is recorded below
  so it does not have to be repeated.

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

2. **What is the benchmark deck set?** Optimization without a
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

No findings recorded yet.

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
