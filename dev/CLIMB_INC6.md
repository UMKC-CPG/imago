# Adaptive Mesh Climb -- Increment 6 working note

Scratch working doc for C118 increment 6 (the final wiring
increment).  It captures the flow-design analysis and the OPEN
design questions so we can resume without re-deriving them.  Once
these questions resolve into DESIGN 5.7 / PSEUDOCODE 11.4, this
note is obsolete and should be deleted (like `DEBUG.md` after its
campaign).  Nothing here is a spec -- the chain edits come after
the questions below are answered.

## Where C118 stands

Increments 0-5 are DONE and committed:

- Inc 0-3b -- `5f9f472` (primitives, decision helpers, climb
  driver, mesh dispatch).
- Inc 4 -- `9659bee` (`converged_mesh` through the guidance
  schema).
- Inc 5 -- `9e621b5` (climb knobs from the manifest
  `[harvest.kpoint_climb]` sub-table).

Increment 6 = wire the producer main to the climb.  Split, with
the programmer's decisions recorded:

- **6a** (design/pseudocode, the gate): rewrite DESIGN 5.7 +
  PSEUDOCODE 11.4 producer flow for the climb; retire the
  density-grid builder and reconcile 15.6.
- **6b** (Python code): rewire `build_initial_potentials` main to
  the climb; retire `build_kpoint_convergence`; orchestration
  tests with the toolchain seam mocked.
- **6c** (Fortran + Python -- IN SCOPE per the programmer): the
  imago `RESOLVED_KP_CLASSES` emit (PSEUDOCODE 4d.5, already
  written) + the producer axis-class self-test (4c.7) + a deep
  fcc-prim class check.  Needs a rebuilt imago binary.

Programmer's session decisions: **talk through the flow design
before rewriting the settled sections**, and **include 6c in
Inc 6**.

## Chain-gate finding (why 6a must come first)

PSEUDOCODE 11.4 ("Build Pipeline") and DESIGN 5.7 still describe
the **density-grid one-shot** producer that C69 wrote:
`build_kpoint_convergence` (a density verify grid) -> one flat
flight -> single dispatch -> `pick_converged_unit`.  Neither
describes the iterative climb.  So the producer main may not be
rewired until 11.4 / 5.7 are rewritten for the climb.

## Current density-grid flow (what 6b replaces)

In `build_initial_potentials.py::build_initial_potentials`:

- Phase 1 build: refresh isolated entries; per solid
  `build_kpoint_convergence` (density grid) + `build_loen_units`
  -> ONE flat flight.
- Phase 1b: `prepare_units` (driver-side makeinput).
- Phase 2: one `dispatch_fn` call.
- Phase 3: per solid `pick_converged_unit` ->
  `discover_environments` -> `extract_potential` per environment
  -> `insert_or_skip`.
- Phase 3b guidance: `guidance_harvest.harvest_flight(workspace)`
  -- RE-READS the whole flight from disk and stages guidance
  entries.

Already built and unit-tested (inc 0-5), waiting to be wired:
`predict_kpoint_density`, `build_mesh_unit`, `make_dispatch_round`,
`climb_action`, `converge_by_climb`, `record_converged`
(producer); `axis_classes_for_cell`, `climb_policy_from_manifest`,
`resolve_climb_policy`, `initial_meshes`, `at_ceiling` (mesh_climb);
`ClimbConfig` / `Rung` types.

## Open design questions (resume HERE)

**Q1 -- Guidance contribution: in-memory vs workspace re-read.**
The climb holds the converged rung + full ladder in memory, and
4e.6 has `record_converged` *feed* `build_entry`.  So the climb
path is: per converged material, `record_converged(rung, rungs,
config)` -> grid arrays; read the converged rung's `result.toml`
for gap/mag/scf; `build_entry(...)` -> `save_entry`.  That
REPLACES the `harvest_flight(workspace)` re-read for the producer.
- PROPOSAL: climb producer contributes guidance in-memory; the
  standalone `harvest_flight` CLI (guidance_harvest.py:519, its
  `__main__` at :682) stays as-is for density-grid flights --
  untouched, not retired.
- REFACTOR WORRY: `harvest_flight` and the climb path both need
  the per-structure "pick -> build_entry -> save_entry" core.
  PROPOSAL: factor that core into a shared helper both call, so
  they cannot diverge.  << programmer to confirm >>

**Q2 -- `build_entry`'s interface** (guidance_harvest.py:378).
Today it takes the density-grid arrays (`kpoint_densities`,
`energies`, `result_tomls`, `idx`, collapsed arrays).  The climb
supplies `(converged_kpoint_density, converged_mesh, grid_values,
grid_energies)` from `record_converged` plus the chosen run's
`result.toml`.  PROPOSAL: refactor `build_entry` to take the
already-chosen `(grid_values, grid_energies, converged_density,
converged_mesh, chosen_result)` explicitly -- one shape both the
density harvest (via `collapse_by_mesh`) and the climb (via
`record_converged`) produce, so there is ONE entry builder.
<< programmer to confirm; this is the other real architecture call
 alongside Q1 >>

**Q3 -- `-loen` fingerprint units (SETTLED-looking).**
Geometry-only, mesh-independent (`build_loen_units`).  They do not
belong in a climb round (the round adapter builds `kpt-mesh`
units).  PROPOSAL: dispatch all `-loen` units once in a small
separate pre-flight batch before the climb; their run dirs persist
for the Phase 3 fingerprint harvest exactly as now.

**Q4 -- Potential extraction (SETTLED-looking, mechanical).**
The converged rung's mesh -> its `kpt-mesh-<a>-<b>-<c>` run dir;
Phase 3 `extract_potential` / `discover_environments` / the
fingerprint harvest read that dir.  The producer maps
`outcomes[m].mesh` -> the run dir to locate the potential.

## After the questions resolve

1. 6a: rewrite DESIGN 5.7 + PSEUDOCODE 11.4 for the climb (build
   ClimbConfig + seed per solid; pre-flight `-loen` batch;
   iterative `converge_by_climb`; in-memory guidance via
   `record_converged` -> shared `build_entry`); retire
   `build_kpoint_convergence` and reconcile PSEUDOCODE 15.6.
2. 6b: rewire the producer main; retire the density-grid builder;
   orchestration tests with the seam mocked.
3. 6c: imago `RESOLVED_KP_CLASSES` emit (4d.5) + producer
   axis-class self-test (4c.7) + fcc-prim check; needs a rebuilt
   binary.
4. Follow-on (cluster): live seed re-run -> closes C116, unblocks
   C117.
