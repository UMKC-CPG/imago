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

## Open design questions

**Q1 + Q2 -- RESOLVED (2026-07-13): chosen-facts `build_entry`.**
Q1 (in-memory vs workspace re-read) and Q2 (`build_entry`'s
interface) collapse into one refactor.  The two guidance paths
differ ONLY in how they *pick* the converged point; once picked,
both hand the identical already-chosen facts to one entry builder.

Refactor `build_entry` (guidance_harvest.py:378) to the
chosen-facts shape:

    build_entry(workspace_root, source_structure, prediction,
                dataspace, structure, kpoint_threshold,
                grid_values, grid_energies,        # the ladder
                converged_density, converged_mesh, # the chosen rung
                chosen_result)                     # result.toml

- The DENSITY harvest assembles those in `harvest_flight`'s loop
  from `collapse_by_mesh` + `pick_converged` (the pick stays
  where it already lives).
- The CLIMB assembles them from `record_converged` (grid arrays +
  converged density/mesh) plus the converged rung's `result.toml`.
- Shared stage core: `save_entry(build_entry(...))` -- a thin pair
  both callers go through, so the two paths CANNOT drift on a
  schema change.  Drops the density-only args (`grid`,
  `kpoint_densities`, `energies`, `result_tomls`, `idx`,
  `collapsed_*`) and the internal `result_tomls[idx]` re-pick.
- The standalone `harvest_flight` CLI (guidance_harvest.py:519,
  `__main__` at :682) stays intact for density-grid flights --
  untouched, not retired.

This is the last architecture call; Q3/Q4 below were already
settled.  6a may now proceed.

**Refinement while writing 6a (2026-07-14):** `build_entry` takes
**10** args, not 11 -- `converged_mesh` is NOT passed.  Both paths
already hand `build_entry` the chosen rung's `result.toml`, and the
exact mesh lives in it (`kpoint_mesh`, DESIGN 6.1.2), so build_entry
reads the mesh from `chosen_result` in BOTH paths.  `record_converged`
supplies only the density and the flatness ladder.  This drops one
arg from the Q1-Q2 preview and removes the need for the climb to
thread a mesh through.  Landed in PSEUDOCODE 15.7 (definition), 11.4
(climb call), and the 4e.6 note.

**Scoping move (2026-07-14):** retiring the pseudocode for
`build_kpoint_convergence` (15.6) and `pick_converged_unit` is moved
to **6b**, done together with removing their code.  Deleting the
pseudocode while its code still exists would create
code-without-pseudocode -- the wrong direction through the gate.  So
after 6a the state is deliberately: 11.4/5.7 describe the climb, and
15.6 + the retired-helper code still stand until 6b removes both at
once.  DESIGN 6.2.8 now carries a forward-pointer (added in 6a) so a
reader lands on the climb; 15.6 reconciliation itself waits for 6b.

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

1. 6a (DONE 2026-07-14): rewrote DESIGN 5.7 + PSEUDOCODE 11.4 for
   the climb (build ClimbConfig + seed per solid; pre-flight `-loen`
   batch; iterative `converge_by_climb`; in-memory guidance via
   `record_converged` -> shared 10-arg `build_entry`); added the
   DESIGN 6.2.8 forward-pointer.  15.6 / `pick_converged_unit`
   retirement deferred to 6b (see the scoping move above).
2. 6b: rewire the producer main; retire the density-grid builder
   (`build_kpoint_convergence`) AND its pseudocode 15.6, plus
   `pick_converged_unit` and its pseudocode, together with the code;
   orchestration tests with the seam mocked.
3. 6c: imago `RESOLVED_KP_CLASSES` emit (4d.5) + producer
   axis-class self-test (4c.7) + fcc-prim check; needs a rebuilt
   binary.
4. Follow-on (cluster): live seed re-run -> closes C116, unblocks
   C117.
