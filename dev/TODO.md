# Task List

> **Document hierarchy:** Tasks are organized by the level of
> the design chain they affect. Each item cites the relevant
> document section.

---

## VISION

- [ ] V1. Confirm design principles against additional test
  cases beyond KNbO3 (VISION Principles.1)
- [x] V2. Add Goal 5 (historical-guidance dataspace),
  Principle 11 (experience as a curated artifact, not tribal
  knowledge), and Principle 12 (the flight layer stays
  dumb; flight description lives in Python) to VISION.
  Done 2026-05-28 as part of the original DESIGN 7 / ARCH 10
  chain landing (categorical signature shape); revised
  2026-05-29 (Path B) when Goal 5 was rewritten to frame
  the artifact as a regression-trained dataspace rather than
  a categorical lookup -- a small k-NN predictor learns
  composition -> electronic-character -> k-density from a
  curated set of converged calculations, and verifies its
  prediction with a small grid whose width tracks the
  predictor's uncertainty.  Principle 11 generalises this
  curation discipline back to the initial-potential DB
  (Goal 3) so both Goals share one philosophy.  Principle 12
  promotes the no-DSL strategic decision out of memory and
  governs DESIGN 6.2 and 7.
- [x] V3. Add Goal 6 (resource & cost guidance dataspace) to
  VISION.  Done 2026-05-29.  A hardware-aware sibling of Goal
  5: records per-run cost (peak memory, disk, walltime) vs
  problem size, parallel configuration, and build/toolchain so
  the flight layer provisions SLURM requests instead of
  guessing.  Kept separate from Goal 5 because cost is
  hardware-specific where convergence guidance is portable
  (ARCHITECTURE 11, DESIGN 8).

---

## ARCHITECTURE

- [x] A1. Add numTetrahedra, tetraVol, tetrahedra(:,:),
  numFullMeshKP, fullKPToIBZKPMap(:) to O_KPoints module
  (ARCHITECTURE 2, DESIGN 2.5)
- [x] A2. electronPopulation_LAT lives in populate.F90
  (O_Populate), alongside electronPopulation. Tetrahedra
  data passed from O_KPoints as arguments.
  (ARCHITECTURE 2)
- [ ] A3. Introduce compile-time working precision kind
  (`wp`) in kinds.f90; propagate to all floating-point
  declarations, constants, and I/O type tags so a single
  build flag switches between double and single precision
  (ARCHITECTURE 6.1)
- [ ] A4. Restructure alpha-pair inner loops in
  integrals.F90 (gaussOverlapOL and siblings): separate
  selection phase (which pairs survive alphaDist test)
  from compute phase (evaluate integral), producing a
  packed list for branchless SIMD execution
  (ARCHITECTURE 6.2)
- [ ] A5. Identify accumulation sites that require
  compensated (Kahan) summation or reordering for
  numerical stability at single precision; implement
  guards and validate (ARCHITECTURE 6.1)
- [ ] A6. GPU offload of restructured integral compute
  phase via OpenACC, CUDA Fortran, or OpenMP target
  (ARCHITECTURE 6.3)
- [x] A8. Write ARCHITECTURE §10 (historical guidance
  dataspace, VISION Goal 5).  Done 2026-05-28 (categorical
  signature shape); rewritten 2026-05-29 (Path B) when
  signature changed from (elements + stoichiometry +
  structural_class) to (system_type + composition_vector +
  lattice_family) with k-NN regression replacing the
  categorical lookup.  Final shape (eight subsections):
  10.1 layout under share/historicalGuidanceDB/ with
  entries/<system_type>/, staging/<system_type>/,
  SCHEMA_VERSION marker, and elemental_groups.toml table; 10.2
  TOML format with signature-first + measured + context +
  verification + provenance invariants; 10.3 data flow
  diagram showing the two-stage k-NN predictor and the
  variance-aware widening; 10.4 feature space (13-d
  composition + 6-axis lattice family + 4-way system_type)
  and the k-NN predictor; 10.5 curation/regeneration/harvest
  with seed-flight auto-promote rule; 10.6 module impact
  (guidance_db.py library, kaleidoscope helper,
  guidance_harvest.py, guidance_promote.py, future
  guidance_migrate.py, plus a small imago.py extension to
  expose gap/spin/dos in result.toml, tracked as C76); 10.7
  relationship to other prongs (no-link with DESIGN 5
  closed by decision); 10.8 open architectural questions
  including polytype confusion + spin-pol-vs-AFM
  interpretation + functional/basis as sub-models vs
  features.
- [x] A7. Write ARCHITECTURE §9 (high-throughput
  calculation flights / "kaleidoscope", VISION Goal
  4).  Layers + one-directional dependency graph (9.1);
  imago.py CLI+callable API with two entry modes (9.2);
  ase_imago.py ASE Calculator kept separate, with the
  ASE-free StructureControl factory / adapter-glue split
  (9.3); kaleidoscope/ Parsl dispatcher with a pluggable
  wingbeat seam (9.4); cod_fish.py + cif2skl.py
  acquisition front-end (9.5); workspace layout +
  general run-reuse cache (mechanism in kaleidoscope,
  key fields per client) (9.6); clients +
  producer-as-client reframing (9.7).  Done 2026-05-20
  (commit 2dc5b30).  The workspace scheme (9.8), left
  open then, is now resolved by D13 / DESIGN 6.2.4 and
  ARCH 9.6/9.8 updated to match (2026-05-21).
- [x] A9. Write ARCHITECTURE §11 (resource & cost guidance
  dataspace, VISION Goal 6).  Done 2026-05-29.  Eight
  subsections paralleling §10: 11.1 relationship to §10 (sibling
  not extension; portability is why they stay separate); 11.2
  layout partitioned by hardware fingerprint + the
  single-execution-observation atomic unit; 11.3 the four
  observation blocks (size signature incl. 1-Schrodinger /
  4-Dirac components + secular dimension, execution config,
  two-layer build config, measured resources); 11.4 capture
  from four sources (dispatch config, CMake build_info, sacct,
  imago self-report) + retained censored failed runs; 11.5
  feature space + physics-informed regressor + provisioning
  consumer; 11.6 module impact (resource_db.py,
  resource_harvest.py, provisioning consumer, future
  resource_migrate.py, CMake build_info hook); 11.7 relationship
  to other prongs; 11.8 open questions (fingerprint + build-knob
  granularity, peak-memory capture, node contention, censored
  data, portable normalization).

---

## DESIGN

- [x] D1. Resolve memory strategy for PDOS projection array
  P(alpha, n, k): store at IBZ k-points only, apply
  atom permutation on-the-fly during tetrahedron
  corner assembly (DESIGN 1.4)
- [x] D2. Verify that electronPopulation_LAT correctly
  replaces electronPopulation for bond order in all
  cases (DESIGN 1.5, 2.5) -- resolved: weight is
  correct but accumulation needs atom permutation
- [x] D3. Design IBZ warning/detection for Gaussian PDOS
  path -- resolved: makeinput mesh mode now writes
  style code 1 (axial counts), so Imago builds the
  full mesh internally for both Gaussian and LAT.
  Style code 0 (explicit list) retained for special
  cases with a prominent warning (DESIGN 2.6)
- [x] D4. Write PSEUDOCODE for atom permutation fix:
  atomPerm table (section 4), fullKPToIBZOpMap save
  (section 5), corrected Q* (section 6), corrected
  bond order (section 7) (DESIGN 2.4, 2.5)
- [x] D7. Design initial SCF potential database: TOML
  schema and validation rules (5.2), in-memory
  representation (5.4), deterministic hand-formatted
  emitter (5.5), lookup with isolated/legacy fallback
  (5.6), layered-reproducibility regeneration pipeline
  (5.7), and 20%-iteration-reduction validation harness
  (5.8) (DESIGN 5; VISION Goal 3; ARCHITECTURE 8)
- [x] D8. Finalize curation manifest schema.  Done 2026-05-11.
  Schema v1 (TOML), validation rules, cache layout, and COD-fetch
  contract are fully specified in DESIGN 5.7.  ARCHITECTURE 8.8's
  open question on file format is closed in favor of TOML, with
  the full spec referenced into DESIGN 5.7.  Key decisions: TOML
  format; two-level `[[reference_solid]]` / `[[reference_solid.entry]]`
  shape with all fields required in v1; manifest lives in
  `share/atomicBDB/`; per-solid source is `cod_id + cod_revision`
  XOR `structure_path`, with structures fetched from the
  Crystallography Open Database at regeneration time using the
  pinned revision; six strict-refusal validation rules;
  direct-comparison cache under `share/atomicBDB/cache/scf/<reference_id>/`
  with byte-compared structure file copies (no hashing) for
  debuggability; strict refusal on COD-fetch failures.
- [x] D9. Phase 2 design: schema v2 with per-entry
  `default` tag and `[[potential.fingerprint]]` sub-blocks
  (DESIGN 5.2); gold sketch update (DESIGN 5.3); in-memory
  `FingerprintRecord` and `default_entry()` (DESIGN 5.4);
  full composition selection algorithm with spatial-scope
  references between `-target`/`-block name=NAME` and
  `-reduce`/`-bispec scope=NAME` (DESIGN 5.6); producer-
  side fingerprint harvest with loen runs cached per
  `(method, sub_spec)` (DESIGN 5.7 v2); two of three
  Phase-2 open questions resolved (DESIGN 5.9, with
  interpolation parked for Phase 3 and element-aware
  bispec carried as a Phase-2 follow-up D10);
  nested-makeinput bootstrap for Fortran-side matchers
  (DESIGN 5.10); matcher protocol architectural section
  inside `makeinput.py` (ARCHITECTURE 8.9).  Done
  2026-05-19, refined the same day to shrink
  ARCHITECTURE 8.2 to a pointer at DESIGN 5.2/5.3,
  rewrite ARCHITECTURE 8.4 around the matcher-driven
  reality, clarify ARCHITECTURE 8.7's
  Fortran-changes scope (Phase-2 base vs follow-up),
  swap "species centroid" for a per-matcher
  `representative` method (DESIGN 5.6.5 / ARCHITECTURE
  8.9), correct DESIGN 5.10.2 step 2 (no `-pot` flag in
  the nested call; default tag flows through), and
  enumerate the full bispec LOEN parameter contract
  (DESIGN 5.10.5).  Element-aware bispectrum design
  is now its own task D10.
- [ ] D10. Design the element-aware bispectrum
  (Phase-2 follow-up).  Specify the `bispecByElement`
  parameter in `O_Input`; the per-neighbor-element
  accumulation in `computeBispectrumComponent`; the
  extended `fort.21` output format (per-element
  vector slices labeled in atomic-number order); the
  `(neighbor_element, vector)` payload shape in
  `[[potential.fingerprint]]` records; the matcher-
  distance semantics that zip by neighbor-element
  symbol with missing-element fallback to L2 of the
  present side; the `by_element = true` sub_spec key
  in `BispecMatcher.to_loen_input` per DESIGN 5.10.5;
  and how the producer harvests the element-aware
  variant.  Adds at least one new sub-section under
  DESIGN 5 (proposed 5.11) and updates DESIGN 5.9 to
  fold this back in.  Must land before C62
  implementation begins.

### Kaleidoscope prong (VISION Goal 4, ARCHITECTURE 9)

- [x] D11. Design the imago.py callable API (ARCH 9.2):
  the result object, the two entry granularities
  (prepared run dir; structure+options driving
  makeinput), the CLI-as-thin-wrapper refactor, and how
  the existing checkpoint and lock-file behavior is
  preserved.  Foundation for D12/D13.

  Done 2026-05-21.  Landed as new DESIGN section 6
  (the kaleidoscope DESIGN home; DESIGN 6 <-> ARCH 9,
  same offset as DESIGN 5 <-> ARCH 8).  Section 6.1
  covers: 6.1.1 the C48.3 client requirements stated up
  front (converge verdict, converged-scfV path, run
  conditions, iteration count); 6.1.2 the `ImagoResult`
  dataclass + `RunStatus` enum (CONVERGED /
  NOT_CONVERGED / FAILED / SKIPPED); 6.1.3 the two entry
  modes (`run_prepared`, `run_structure`) and the CLI
  split into parse-argv / pick-mode / translate-result;
  6.1.4 the private run core and cwd-restore discipline
  for the long-lived kaleidoscope worker; 6.1.5
  lock-file (already per-run-dir) and within-run-dir
  checkpoint preservation, with the boundary against
  kaleidoscope's coarser run-reuse cache; 6.1.6 open
  details deferred to PSEUDOCODE/C63 (output-key
  enumeration, iteration/energy parsing, StructureControl
  input type, API-mode call provenance).  Key judgments:
  run-level failures are returned as statuses while only
  contract failures raise (Principle 10, never sys.exit
  in the API path); cwd restored on exit.  PSEUDOCODE for
  D11 still pending (after D13 lands, then /refine).
- [ ] D12. Design ase_imago.py ImagoCalculator (ARCH
  9.3): implemented_properties (energy / forces / stress
  / charges / ... plus Imago specialties as custom
  result keys), unit conversions, and the ASE-free
  StructureControl factory in structure_control.py plus
  the adapter-layer Atoms-reading glue.
- [x] D13. Design the kaleidoscope flight dispatcher (ARCH
  9.4, 9.6): the Parsl dispatch model, the pluggable
  wingbeat seam, complete-and-report status tracking, the
  workspace layout + id / <calc> tag / status.toml
  scheme (resolves the ARCH 9.8 open item), and the
  run-reuse cache mechanism with client-supplied key
  fields.

  Done 2026-05-21.  Landed as DESIGN section 6.2,
  building on 6.1 (the default wingbeat calls the 6.1 API
  and persists the 6.1.2 ImagoResult).  Subsections:
  6.2.1 the domain-agnostic data model (CalcUnit /
  Flight, with flight.toml as a generated record);
  6.2.2 the pluggable wingbeat seam (Wingbeat.run -> generic
  WingbeatOutcome with an opaque wingbeat-supplied `detail`;
  ImagoWingbeat maps ImagoResult and persists result.toml);
  6.2.3 Parsl dispatch (one python_app per unit; both
  sweep and tight-loop shapes; per-future exception catch
  for Principle 10); 6.2.4 the workspace layout that
  RESOLVES ARCH 9.8 (id charset/uniqueness rule, optional
  <calc> tag + derivation, the five-value status.toml
  schema with convergence carried in `detail` not
  `status`); 6.2.5 the run-reuse cache (mechanism/policy
  split; key = verbatim scalar fields + byte-compared key
  files, generalizing the producer's is_cached_v2;
  resume == re-run); 6.2.6 client-side harvest via the
  run dir + the FlightReport (subsumes C48.3's
  producer-as-client shape); 6.2.7 open details for
  PSEUDOCODE/C68.  Key judgment: kaleidoscope stays
  domain-agnostic (Principle 9) -- it tracks generic
  lifecycle status and records but never interprets the
  wingbeat's `detail`; all harvest is client-side.

  Refined 2026-05-21: ARCH 9.6/9.8 updated to mark the
  workspace scheme resolved by DESIGN 6.2.4 (the open
  item under A7).  PSEUDOCODE follows as P6/P7.
- [ ] D14. Design structure acquisition (ARCH 9.5):
  cod_fish.py COD-fetch contract (urllib, pinned
  cod_revision, strict on failure) and cif2skl.py (ASE
  CIF read -> StructureControl factory -> skl write).
  PSEUDOCODE for D11-D14 follows once each design lands.
- [x] D16. Design the historical-guidance dataspace
  (VISION Goal 5, ARCHITECTURE §10).  Done 2026-05-28
  (categorical signature shape, Jaccard lookup); rewritten
  2026-05-29 (Path B) when the signature changed from
  (elements + stoichiometry + structural_class) to
  (system_type + 13-d composition_vector + 6-axis
  lattice_family) and the lookup became a two-stage k-NN
  regression with variance-based confidence.  Final shape
  (ten subsections): 7.1 motivation + why a dataspace and
  predictor rather than a categorical lookup; 7.2 TOML
  schema v1 (signature / measured / context / verification /
  provenance blocks with 12 validation rules); 7.3 worked
  gold sketch (TiO2-rutile entry); 7.4 in-memory dataclasses
  (Signature, Measured, Context, Verification, Provenance,
  GuidanceEntry, Dataspace, PredictionResult) + the public
  surface (load, save_entry, compute_signature, predict) +
  the elemental_groups.toml element-classification table layout
  per Principle 11; 7.5 deterministic hand-formatted TOML
  emitter + the `<system_type>-<short_sha>` slug derivation;
  7.6 PREDICTOR ALGORITHM (two-stage k-NN with
  inverse-distance weighting for crystalline; canonical
  entry for non-crystalline; sub-models per
  (basis, functional); variance-based confidence;
  is_under_trained flag); 7.7 predict-then-verify flight
  construction + verification-grid widening driven by
  predictor confidence + trust-mode (verify=False) for
  nearly-identical-family flights; 7.8 harvest pipeline
  (reads gap/spin/dos from result.toml; auto-promote rule
  for the seed flight + interactive curator review);
  7.9 bootstrap (canonical entries seeded by hand for
  non-crystalline; wide-grid default for under-trained
  crystalline; non-convergence recovery); 7.10 open
  questions (metalloid assignment, k-NN tuning knobs,
  polytype confusion, AFM total_magnetization, cell-size
  for defects, multi-metric verification, functional/basis
  as sub-models vs features, staleness) + the two
  closed-by-decision records (no link to DESIGN 5 from
  2026-05-28; chemistry not as primary signature axis from
  2026-05-29).
- [x] D17. Patch DESIGN 6.2 for the kaleidoscope <-> guidance
  DB seam.  Done 2026-05-29 in the Path B rewrite session:
  (a) §6.2.4 extended with the **tree-per-varied-axis**
  convention for sweep flights (one directory level per
  varied axis in stable order, with single-tag rules for
  axis names / values / decimal-points, bidirectional path
  parsing, and `flight.toml` recording the axis order +
  fixed axes; chosen over a flat axis-value-string after
  the user flagged filename-length growth); (b) the §6.2.1
  worked example now describes the producer as a
  predict-then-verify client expanding each reference
  solid into a verification sub-grid, including a
  "trust mode for nearly-identical families"
  (`verify=False`) note that builds a length-1 sub-grid
  and skips auto-staging the result; (c) §6.2.8 added as
  a new subsection describing the flight-builder helper
  inside src/scripts/kaleidoscope/ (predict_settings,
  build_verification_grid, build_calc_tag, the
  PredictionRecord shape, cross-references to DESIGN 7.6 /
  7.7 / 7.9 and to 6.2.4's tag convention); (d) the §6.2
  intro paragraph now names Principles 8 / 10 / 12 and
  documents the flight-builder split (option-axis sweeps
  inside kaleidoscope; structure-axis sweeps in
  structure_control / acquisition).
- [x] D15. Design the makeinput callable build API (the
  makeinput-side twin of D11/C63): turn makeinput.py from
  an argv-and-cwd-bound script into one that also exposes a
  callable `build_run_dir(structure, options, wingbeat_dir)`,
  with the CLI a thin wrapper.  This is the design rung
  that was missing under C68 item (a): DESIGN 6.1.3 deferred
  `run_structure`'s structure-and-options mode to "drive
  makeinput," but no in-process makeinput entry point
  existed to drive.

  Done 2026-05-21.  Landed as DESIGN section 6.3 (sibling
  to 6.1 imago API and 6.2 kaleidoscope): 6.3.1 why a
  callable build API + the contract/build error boundary
  (new `MakeinputError`, analog of `ImagoError`); 6.3.2 the
  single `build_run_dir` entry point (skl-path structure for
  now) + CLI-as-thin-wrapper; 6.3.3 the `ScriptSettings`
  split (from_command_line / from_options sharing
  reconcile), mirroring C63; 6.3.4 cwd discipline -- stage
  the skl into wingbeat_dir, chdir in, restore cwd in finally
  (reuse the 6.1.4 pattern rather than rewrite makeinput's
  cwd-relative paths) + factor main()'s body into a callable
  build_inputs(settings, sc); 6.3.5 record_clp and the
  _load_rc sys.exit become CLI-only / worker-safe; 6.3.6
  run_structure becomes build_run_dir -> run_prepared,
  closing 6.1.3; 6.3.7 open details for PSEUDOCODE §14
  (options-dict normalization of multi-valued flags,
  StructureControl structure type per C64, nested-makeinput
  bootstrap left as a subprocess).  Key judgment: cwd
  discipline over path-rewriting -- safe because
  kaleidoscope parallelism is across separate workers, each
  with its own cwd.  PSEUDOCODE follows as P8.
- [x] D18. Write DESIGN §8 (resource & cost dataspace schema,
  VISION 6 / ARCH 11).  Done 2026-05-29.  Parallels DESIGN 7:
  TOML schema with 12 validation rules + four observation
  blocks + provenance (8.2); gold sketch (8.3); in-memory
  dataclasses + the three checked-in registries
  (EXECUTION_KNOB_REGISTRY / BUILD_KNOB_REGISTRY /
  RESOURCE_METRIC_REGISTRY) (8.4); hardware-fingerprint recipe
  with microarch normalization (8.5); physics-informed
  power-law regressor in secular_dimension with k-NN fallback
  (8.6); capture + harvest with censored-data handling (8.7);
  cold-start provisioning + manual seed (8.8); open questions
  (8.9).  Pins the two-layer build record: coarse bucketed
  knobs as features + the full verbatim compile_string as
  provenance.

---

## PSEUDOCODE

- [x] P1. Transcribe exact Bloechl middle-range DOS formula
  for e2 <= E < e3 (PSEUDOCODE 2, Bloechl eqs. 14-16)
- [x] P2. Transcribe cornerIntgWt_LAT formulas for partial
  properties (PSEUDOCODE 3a, derived from first principles)
- [x] P3. Write pseudocode for angle clustering and force
  constant computation: greedy cluster-by-triplet, K from
  geometric mean of arm bond stiffnesses, integration
  into create_lammps_files and normalize_types
  (PSEUDOCODE 10, DESIGN 4.8.3-4.8.8)
- [x] P4. Write pseudocode for initial SCF potential
  database: TOML reader with six-rule validation
  (11.1), deterministic hand-formatted emitter (11.2),
  makeinput.py lookup with fallback chain (11.3),
  regeneration pipeline (11.4), and validation harness
  (11.5) (PSEUDOCODE 11, DESIGN 5)
- [x] P5. Write Phase-2 pseudocode for the matcher
  protocol and the new selection path: matcher base
  class plus `ReduceMatcher` / `BispecMatcher`
  concrete subclasses (ARCHITECTURE 8.9); registry
  lookup; preflight coverage check (DESIGN 5.6.3 step
  4); species-pass composition with `scope=NAME`
  spatial reference resolution (DESIGN 5.6.4);
  manifest-entry pick per species via fingerprint
  match, similarity floor, default-tag fallback
  (DESIGN 5.6.5); type-pass inheritance from parent
  species (DESIGN 5.6.6); nested-makeinput bootstrap
  orchestration with recursion guard (DESIGN 5.10).
  Update existing PSEUDOCODE 11.3 to reflect the
  Phase-2 selection flow rather than the Phase-1
  literal-label flow.  Done 2026-05-19: PSEUDOCODE
  11.3 expanded into seven sub-sections (11.3.a
  matcher protocol; 11.3.b preflight + coverage
  check; 11.3.c species pass with scope resolution;
  11.3.d entry pick with three-step precedence;
  11.3.e type pass with XANES split; 11.3.f bootstrap
  with recursion guard; 11.3.g driver) and 11.4 was
  refreshed to schema v2 in the same pass (manifest
  rules 1-9, is_cached_v2 with field-by-field +
  byte-compare, fingerprint harvest split between
  Python-side and Fortran-side matchers, default tag
  carried into PotentialEntry).

### Kaleidoscope prong (VISION Goal 4, ARCHITECTURE 9)

- [x] P6. Write PSEUDOCODE for the imago.py callable API
  (DESIGN 6.1, D11): the `ImagoResult` dataclass and
  `RunStatus` enum (6.1.2); the `run_prepared` and
  `run_structure` entry points and the CLI wrapper's
  parse-argv / pick-mode / translate-result-to-exit-code
  split (6.1.3); the private run core with cwd-restore
  discipline and the returned-status-vs-raised-error
  boundary (6.1.4); the per-run-dir lock lifecycle and
  the within-run-dir checkpoint short-circuit yielding
  SKIPPED (6.1.5).  Resolve the 6.1.6 open details where
  pseudocode forces the choice: the per-job-type
  `outputs{}` key enumeration factored out of the
  `_manage_*_output` helpers, and the iteration/energy
  parse sources.  Foundation for C63.

  Done 2026-05-21.  Landed as PSEUDOCODE section 12
  (12.1 result/status; 12.2 the `project_home_outputs`
  single-source-of-truth output-name table; 12.3 entry
  points + CLI wrapper with the `ScriptSettings`
  from_command_line/from_options split; 12.4 the
  reentrant `_run_core` with lock lifecycle, cwd restore
  in `finally`, and the raise-vs-return error boundary;
  12.5 harvesting).  Two facts grounded against the code:
  `safe_append` skip is 1-based so the iteration file has
  ONE header line written only on first creation, and
  reruns append data rows.  P6 surfaced a gap DESIGN
  6.1.2 had assumed away -- the driver only checks fort.2
  (ran-without-abort), with no converged-vs-ceiling
  signal -- and resolved it WITHOUT a Fortran change: read
  the iteration file's last data row once for the
  convergence metric (col 4, vs the line after
  CONVERGENCE_TEST in imago.dat), the total energy (col
  5), and the per-run iteration count (col 1, a counter
  that resets each SCF run so appends never inflate it).
  Note for the next refine: DESIGN 6.1.6 correctly
  pre-deferred these to the pseudocode pass (no drift),
  but a refine may want to record in DESIGN 6.1.2 that
  NOT_CONVERGED is now backed by the col-4/CONVERGENCE_TEST
  mechanism.
- [x] P7. Write PSEUDOCODE for the kaleidoscope flight
  wingbeat (DESIGN 6.2, D13): the `CalcUnit` / `Flight`
  data model and `flight.toml` serialization (6.2.1);
  the `Wingbeat` protocol and `WingbeatOutcome`, with
  `ImagoWingbeat` mapping `ImagoResult` and persisting
  `result.toml` (6.2.2); the Parsl per-unit dispatch with
  per-future exception capture for complete-and-report
  (6.2.3); the workspace id/`<calc>`/`status.toml`
  scheme (6.2.4); the cache hit-test (scalar-field
  compare + key-file byte-compare) and resume-as-re-run
  (6.2.5); the `FlightReport` and client-side harvest
  handoff (6.2.6).  Foundation for C68.

  Done 2026-05-21.  Landed as PSEUDOCODE section 13,
  helpers-first then driver: 13.1 data model
  (CalcUnit/KeyFields/Flight) + serialize_flight;
  13.2 Wingbeat protocol + WingbeatOutcome + ImagoWingbeat (calls
  the §12 API, persists result.toml, maps status ->
  ok/detail); 13.3 unit_run_dir + validate_flight (slug
  rule, <calc> derivation, collision abort) + status.toml
  read/write; 13.4 is_cache_hit / cache_key_matches
  (verbatim scalar compare + key-file byte-compare, no
  hashing) / write_cache_key; 13.5 dispatch +
  dispatch_unit (hit -> completed_future; miss -> queued
  + python_app) + execute_wingbeat_task + collect_future (the
  ParslTaskLost -> "lost" vs worker-exception -> "failed"
  split, per-future capture); 13.6 ReportEntry /
  FlightReport / report_entry_from_status + the
  client-side harvest_converged_potentials example fixing
  the C48.3 producer contract.  All five PSEUDOCODE items
  for the kaleidoscope prong (with §12) now done; the
  D11/D13 design + pseudocode chain is complete and
  consistent.  Next: code -- C63 (imago.py API, P6/§12)
  then C68 (kaleidoscope, P7/§13), then C69 + C48.3.
- [x] P8. Write PSEUDOCODE for the makeinput callable build
  API (DESIGN 6.3, D15): the `ScriptSettings`
  from_command_line / from_options split sharing reconcile
  (6.3.3); `build_run_dir(structure, options, wingbeat_dir)` with
  skl staging and the chdir / restore-in-finally cwd
  discipline (6.3.4); `build_inputs(settings, sc)` factored
  out of main() so CLI and API share one build sequence;
  the CLI wrapper's parse-argv / build / translate-error
  split and the `MakeinputError` raise-vs-sys.exit boundary
  (6.3.2, 6.3.5); and `imago.run_structure` as build_run_dir
  -> run_prepared (6.3.6).  Resolve the 6.3.7 open details
  pseudocode forces: the options-dict normalization of the
  multi-valued flags (reduce / target / block / xanes) into
  the args namespace reconcile expects.  Foundation for
  C68(a).

  Done 2026-05-21.  Landed as PSEUDOCODE section 14: 14.1
  the ScriptSettings split + build_args_namespace, with the
  multi-valued-flag normalization pinned (client supplies
  the same list-of-token-lists shape argparse yields, placed
  under args.<dest> verbatim, so a dict-described and a
  flag-described run are byte-identical after reconcile);
  14.2 build_inputs (main()'s body, one shared definition)
  and build_run_dir (skl staging + chdir/restore-in-finally,
  lock-free since the per-run-dir lock is taken later by
  _run_core); 14.3 the cli_main wrapper with record_clp
  CLI-only and MakeinputError -> exit code; 14.4 the
  completed run_structure.  Also reconciled §12.3's
  run_structure (P6 wrote it calling _run_core directly) to
  call run_prepared, matching DESIGN 6.3.6 and §14.4 -- no
  functional change (the dir exists post-build), removes
  chain drift.  Next: code -- C68(a).

### Historical guidance dataspace (VISION 5, ARCH 10, DESIGN 7)

- [x] P9. Write PSEUDOCODE for the historical-guidance
  dataspace library, the predictor, the flight-builder
  helper, and the harvest pipeline (DESIGN 7, D16).  Four
  blocks:
  (a) **Library + I/O.**  load() walking
  entries/<system_type>/ with all 12 validation rules
  yielding clear file/block/field error messages; the
  elemental_groups.toml loader and the element-to-group lookup;
  compute_signature(structure, system_type) producing the
  13-d composition_vector via atom-fraction + the
  lattice_family one-hot; save_entry() with the
  deterministic hand-formatted emitter (16-sig-digit
  floats, fixed block sequence, multi-line composition
  vector layout, comma-trailing array layout) and the
  `<system_type>-<short_sha>` slug derivation; the
  in-memory Dataspace partitioned by system_type.
  (b) **Predictor (DESIGN 7.6).**  The system_type switch
  (canonical entry for non-crystalline; two-stage k-NN for
  crystalline); sub-model selection by (basis, functional)
  with the (basis -> functional-family -> overall pool)
  fallback chain; stage-1 distance d1 over composition +
  lattice_family with inverse-distance weights;
  predicted_gap + predicted_magnetization + confidence_1
  via weighted variance; stage-2 distance d2 over predicted
  electronic character with inverse-distance weights;
  predicted_kpoint_density + confidence_2; combined
  confidence; PredictionResult assembly including the
  is_under_trained flag.
  (c) **Flight-builder helper in
  src/scripts/kaleidoscope/.**  predict_settings(structure,
  options, dataspace, system_type, basis, functional,
  verify, id, extra_axes) per DESIGN 6.2.8; the
  build_verification_grid(center, confidence) widening
  function; the wide-grid fallback for is_under_trained;
  the trust-mode (verify=False) length-1 grid;
  build_calc_tag(calc_axes) emitting the
  tree-per-varied-axis paths per DESIGN 6.2.4;
  PredictionRecord attachment to the flight.
  (d) **Harvest pipeline (guidance_harvest.py).**  Group
  CalcUnits by id; sort each verification sub-grid by
  k-density; parse result.toml for each converged calc
  (gap_ev, gap_kind, spin_polarization,
  total_magnetization, total_energy);
  pick the converged grid point with the two-sided
  delta-below-threshold rule (DESIGN 7.8 step 3c);
  SKIP-and-tag-prediction_mismatch on non-convergence at
  the top; recover predictor_confidence and
  predictor_neighbor_ids from [flight.prediction];
  build a GuidanceEntry; write it to
  staging/<system_type>/ via save_entry().  Plus a
  smaller block describing guidance_promote.py's three
  modes (interactive review default, --auto-promote with
  the middle-60%-of-grid + top-three-energy-variance
  rule, --dry-run).
  Foundation for C70-C73.

  Done 2026-05-28.  Landed as PSEUDOCODE section 15, in
  seven subsections: 15.1 constants + the 7.4 dataclasses
  restated in field order (Verification carries the new
  grid_energies); 15.2 load_elemental_groups + compute_signature
  with the CRYSTAL_SYSTEM_TO_FAMILY map (trigonal lumped
  into hex per 7.10); 15.3 load() + load_entry/
  load_verification enforcing all 12 rules with file/
  block/field messages; 15.4 the hand-formatted emitter
  (save_entry, short_sha slug, format_entry, the float-
  array layout); 15.5 the predictor (predict switch,
  predict_non_crystalline, select_submodel three-tier
  fallback + functional_family, the shared knn_weights,
  stage1/stage2 with d1/d2 + variance confidences);
  15.6 the flight-builder helper (PredictionRecord,
  default_wide_kpoint_density_grid, logspace +
  build_verification_grid, encode_axis_value/
  build_calc_tag, predict_settings, standard_key_fields);
  15.7 harvest_flight + pick_converged and promote +
  auto_promote_ok.  Two design decisions the pseudocode
  forced, both pre-cleaned in DESIGN/ARCH before P9 was
  written (commits 7976c05 + 194b041): the predict() seam
  (free function predict(dataspace, query, basis,
  functional), not a stale db.predict method) and the
  grid_energies array (so auto-promote runs from a staging
  file alone).  One smaller pinning surfaced WHILE writing
  P9 and is recorded inside 15.6: the PredictionRecord
  reaches flight.toml via a generic opaque Flight.metadata
  dict that serialize_flight emits verbatim as
  [flight.<key>], keeping the dispatch core domain-agnostic
  (Principle 9).  The doc-level definition is now in place
  (the metadata field + serialize loop added to DESIGN 6.2.1
  + PSEUDOCODE 13.1, 2026-05-29); only the model.py
  implementation of that field remains, as part of C71's
  model catch-up.

### Resource & cost dataspace (VISION 6, ARCH 11, DESIGN 8)

- [ ] P10. Write PSEUDOCODE (new §16) for the resource & cost
  dataspace, paralleling §15.  Four blocks: (a) **library +
  I/O** -- load() walking entries/<hardware_fingerprint>/ with
  all 12 validation rules and file/block/field error messages;
  the hardware_registry.toml loader + the fingerprint recipe
  (8.5); the three registry validators (execution / build /
  resource) with the unknown-key-fails-loudly rule; save_entry()
  with the deterministic hand-formatted emitter and the
  two-layer build record (coarse knobs + verbatim
  compile_string); the Observation dataclass + the
  ResourceDataspace partitioned by fingerprint.  (b)
  **Predictor (DESIGN 8.6)** -- the per-(fingerprint,
  build-bucket) physics-informed power-law fit in
  secular_dimension (log-log least squares recovering the
  exponent), the parallel + spin corrections, the k-NN fallback
  for thin groups, and the censored (oom / timeout bound)
  handling.  (c) **Provisioning consumer** -- predict
  mem/disk/walltime for a proposed config + safety margin ->
  SLURM request; the cold-start fallback for empty or
  under-populated fingerprints (8.8).  (d) **Harvest** --
  resource_harvest.py joining the four capture sources (8.7),
  one Observation per run dir, censored failed-run staging,
  staging-then-promote.

---

## CODE

### Phase A -- LAT TDOS (eigenvalues only)

- [x] C1a. Save fullKPToIBZKPMap from initializeKPointMesh
  IBZ folding (preserve kPointTracker as module data)
- [x] C1b. Move generateTetrahedra call from readKPoints
  to initializeKPoints (after mesh is built)
- [x] C1c. Complete generateTetrahedra in kpoints.f90
  (PSEUDOCODE 1)
- [x] C2. Compute tetraVol in initializeKPoints after
  lattice initialization (DESIGN 1.2)
- [x] C3. Implement computeTDOS_LAT in dos.F90 using
  fullKPToIBZKPMap for eigenvalue unfolding (PSEUDOCODE 2)
- [x] C4. Validate LAT TDOS against Gaussian broadening at
  high k-point density
- [ ] C4a. Add Bloechl correction terms (eqs. 22-24) to
  computeTDOS_LAT for improved accuracy at lower
  k-point densities (DESIGN 1.3)

### Phase B -- electronPopulation_LAT (integrated properties)

- [x] C5. Implement computeElectronPopulation_LAT
  (PSEUDOCODE 3) -- in populate.F90
- [x] C6. Modify computeBond to use
  electronPopulation_LAT when kPointIntgCode == 1
  (DESIGN 1.5) -- dispatch via statePopulation local
- [x] C7. Modify effective charge in computeBond to use
  electronPopulation_LAT (DESIGN 1.5) -- handled by
  C6 statePopulation dispatch with 2/spin factor

### Phase C -- LAT PDOS (energy-resolved)

- [x] C8. Implement bloechlCornerDOSWt subroutine:
  per-corner DOS density weights, used by both TDOS
  and PDOS (DESIGN 1.3, PSEUDOCODE 2a)
- [x] C8a. Refactor computeTDOS_LAT to call
  bloechlCornerDOSWt + sum, replacing inline
  dosContrib formulas (DESIGN 1.3, PSEUDOCODE 2)
- [x] C8b. Fix deltaDOS * hartree unit in integrated-
  area diagnostic in computeTDOS_LAT and computeDOS
  (DESIGN 1.3)
- [x] C9. Fix integratePDOS_LAT to call
  bloechlCornerDOSWt instead of bloechlCornerWeights
  (DESIGN 1.4, PSEUDOCODE 8.3)
- [x] C9a. Validate: rerun KNbO3 LAT PDOS and verify
  Spin States Calculated ≈ Spin States Expected and
  LAT PDOS matches Gaussian PDOS shape/magnitude --
  validated: Spin States 76.06 vs 76.86 expected
  (1% gap from band-edge effects), per-atom electron
  counts match Gaussian, PDOS shape correct. Minor
  O-atom symmetry spread (<0.1% integrated) noted.

### Phase D -- Density-based k-point input (DESIGN 3)

- [x] C11. Add `-kpd`, `-scfkpd`, `-pscfkpd` CLI options to
  makeinput.py argument parser (DESIGN 3.1)
- [x] C12. Add `_write_density_kp_file()` to makeinput.py that
  writes style-code-2 k-point files directly (DESIGN 3.4)
- [x] C13. Modify `_make_kp()` in makeinput.py to dispatch
  between mesh mode (makeKPoints) and density mode
  (direct write); all-or-nothing (DESIGN 3.2, 3.3)
- [x] C14. Handle summary output for density mode (print
  density value instead of kp count) (DESIGN 3.5)
- [x] C15. Validate end-to-end: makeinput.py -kpd produces
  a file that imago readKPoints correctly parses as
  style code 2
- [x] C16. Embed point group operations in style-code-2
  kpoint file (_extract_point_ops + updated writer)
- [x] C17. Read point group ops in readKPoints style-2
  branch (NUM_POINT_OPS + POINT_OPS labels)
- [x] C18. Port computeRecipPointOps into kpoints.f90
  (convert abc point ops to reciprocal-space)
- [x] C19. Port IBZ mesh folding into initializeKPointMesh
  (foldMesh algorithm from makeKPoints)
- [x] C20. Wire computeRecipPointOps call into
  initializeKPoints style-2 path
- [x] C21. Validate IBZ reduction: compare density-mode
  kpoint count against makeKPoints for same mesh

### Phase E -- Mesh-mode conversion and style-code-0 warning

- [x] C22. Convert makeinput.py mesh mode (`-kp`, `-scfkp`,
  `-pscfkp`) to write style-code-1 k-point files
  (axial counts + shift + point ops) instead of
  calling makeKPoints (ARCHITECTURE 2, DESIGN 2.6)
- [x] C23. Add style-code-1 branch to readKPoints in
  kpoints.f90: read axial counts, shift, and point
  ops; wire computeRecipPointOps into
  initializeKPoints style-1 path (ARCHITECTURE 5)
- [x] C24. Add prominent warning in initializeKPoints
  when kPointStyleCode == 0 that decomposition
  properties may be incorrect (DESIGN 2.6)
- [x] C25. Validate: makeinput.py -kp produces a style-
  code-1 file that imago parses and reduces
  correctly -- validated by user test run

### Phase F -- Atom permutation fix (DESIGN 2.4)

- [x] C26. Save fullKPToIBZOpMap in initializeKPointMesh
  folding loop: store operation index m when a match
  is found, store 1 (identity) for IBZ representatives
  (PSEUDOCODE 5) -- also renamed fullToIBZMap to
  fullKPToIBZKPMap across kpoints.f90, dos.F90, and
  populate.F90 for consistency with design docs
- [x] C27. Implement buildAtomPerm in atomicSites.f90:
  build atomPerm(numPointOps, numAtomSites) from point
  group operations and fractional atom positions
  (PSEUDOCODE 4) -- data and builder both in
  O_AtomicSites; imports abcPointOps from O_KPoints
- [x] C28. Wire buildAtomPerm call into imago.F90 (not
  initializeKPoints, to avoid circular dependency
  between O_AtomicSites and O_KPoints) for style
  codes 1 and 2, in both SCF and PSCF paths
- [x] C29. Corrected Q* accumulation in computeBond:
  buffer per-atom projection per IBZ kpoint, then
  distribute across star via atomPerm (PSEUDOCODE 6).
  Also applies same star distribution to
  atomOrbitalCharge (per-atom, per-QN_l). Unified
  code path for all style codes: style code 0 now
  sets up trivial identity IBZ maps in kpoints.f90
  (numPointOps=1, identity abcPointOps, identity
  fullKPToIBZKPMap/OpMap) so buildAtomPerm and star
  distribution work without special-casing
- [x] C30. Corrected bond order accumulation in
  computeBond: buffer per-pair overlap per IBZ kpoint
  in ibzBondProj, then distribute across star via
  atomPerm with both atom indices permuted
- [x] C31. Validate: run KNbO3 with IBZ and verify that
  symmetry-equivalent atoms produce identical Q* and
  bond order -- confirmed correct

#### Phase F follow-up -- Basis-invariant on-disk operations (DESIGN 2.7)

Surfaced when diamond/prim (Fd-3m reduced to its
primitive rhombohedral cell) tripped
`buildAtomPerm: no atom match found`.  Root cause:
spaceDB operations were stored in conventional-cell-abc
form, but the loaded lattice after primitive reduction
was the primitive cell, so the matrix-vector product
inside buildAtomPerm mixed two bases.  Same latent issue
applied to every non-cubic system in either full or
prim mode.  Resolution: move the basis change to the
kp-file boundary -- makeinput.py emits Cartesian xyz
operations; imago conjugates into whatever cell ended
up loaded.  See DESIGN 2.7 + PSEUDOCODE 4b.

- [x] C31a. Add `full_cell_real_lattice` to
  StructureControl: declared in `__init__` next to
  `full_cell_mag` and friends, documented in the class
  docstring's "Lattice parameters" section, snapshotted
  at the top of `apply_space_group()` before any other
  state changes, and propagated through the existing
  Angstrom-to-Bohr conversion in `_convert_a_to_au` so
  it stays in sync with `real_lattice` units
  (structure_control.py, makeinput.py)
- [x] C31b. Add `_to_cartesian_ops` helper in
  makeinput.py to convert spaceDB operations (rotations
  and translations) from conventional-cell-abc to
  Cartesian xyz via the standard similarity transform.
  Wired into both density-mode and mesh-mode call sites
  in `_make_kp` between `_extract_point_ops` and the
  kp-file writer.  Includes plain-Python 3x3 helpers
  (`_matmul_3x3`, `_matvec_3x3`, `_transpose_3x3`,
  `_inv_3x3`) -- no numpy dependency
  (PSEUDOCODE 4b.1, DESIGN 2.7)
- [x] C31c. Add `computeRealPointOps` to kpoints.f90 as
  real-space sibling of `computeRecipPointOps`.  Emits
  `abcRealPointOps` and `abcRealFracTrans` in the
  basis of the lattice currently in O_Lattice.  Called
  unconditionally from every style-code branch of
  `initializeKPoints` so consumers can stay branch-free
  (PSEUDOCODE 4b.2, DESIGN 2.7)
- [x] C31d. Rename `abcPointOps` -> `xyzPointOps` and
  `abcFracTrans` -> `xyzFracTrans` across kpoints.f90
  to reflect the new on-disk Cartesian convention.
  Tightened docstrings on both compute siblings and on
  the readKPoints style-1 / style-2 branches.
- [x] C31e. Switch `buildAtomPerm` to consume
  `abcRealPointOps` / `abcRealFracTrans` instead of
  the raw `xyzPointOps` / `xyzFracTrans`.  Operations
  and atom positions now agree on basis (whichever cell
  ended up loaded) and the matrix-vector product is
  meaningful (atomicSites.f90)
- [x] C31f. Validate: diamond/prim (227_a / prim) runs
  through buildAtomPerm and completes SCF.  Confirmed
  by user 2026-05-13.

#### Phase F follow-up -- SYBD path bypasses atomPerm (DESIGN 2.6)

Surfaced when diamond/prim was rerun with `-pscfsybd`
and crashed at atomicSites.f90:357 with
`Array bound mismatch for dimension 1 of array
'abcatompos' (3/1)`.  Root cause: the SYBD branch of
`initializeKPoints` calls `makePathKPoints` and skips
all of the point-ops setup (numPointOps assignment,
xyzPointOps/xyzFracTrans allocation,
computeRealPointOps), leaving `abcRealPointOps` and
`abcRealFracTrans` unallocated.  But `setupSCF` and
`intgPSCF` then called `buildAtomPerm` unconditionally
-- a Phase F wiring gap from C28.  Resolution: skip
`buildAtomPerm` / `buildInvAtomPerm` under SYBD.  Band
structure is per-k eigenvalues along a 1-D path; no
shell-summed quantities to unfold; the future partial-
decomposition path is direct (per-atom projection at
the very k-point being plotted) and will not need
atomPerm either.  See DESIGN 2.6.

- [x] C31g. Guard `buildAtomPerm` and `buildInvAtomPerm`
  in setupSCF with `if (doSYBD_SCF /= 1)`, and the
  matching pair in intgPSCF with `if (doSYBD_PSCF /= 1)`.
  Added `doSYBD_SCF` to the `O_CommandLine` use clause
  in setupSCF (already imported in intgPSCF).  Comment
  blocks at both call sites explain why SYBD does not
  need the table (imago.F90).
- [x] C31h. Document the SYBD-path bypass in DESIGN 2.6
  alongside the existing style-code-0 warning, including
  the consumer audit (computeBond / LAT-PDOS channel
  permutation are both already gated by `doBond_*` /
  `doDOS_*`) and the failure mode for the unsupported
  `-sybd + -bond` or `-sybd + -dos` combinations.
- [ ] C31i. Validate: rerun diamond/prim with
  `-pscfsybd` and confirm it completes through
  `bandPSCF` and `printSYBD` (job at
  jobs/c/diamond/prim).
- [ ] C31j. (Follow-up, optional) Add an explicit early
  refusal in `Imago` when `doSYBD_*` is combined with
  `doBond_*` or `doDOS_*`, so the user gets a clear
  message instead of an unallocated-`atomPerm` crash
  downstream.  Defer until a real user actually hits
  that combination; current behavior already fails
  loudly rather than silently producing wrong answers.

#### Phase F follow-up #2 -- Conv-abc on-disk redesign (DESIGN 2.7 revised)

Iterated on the C31 fix.  The Cartesian-xyz on-disk
intermediate worked but pushed a similarity transform onto
the producer side and put numerical values on disk that no
longer traced back to `share/spaceDB/<sg>`.  Reworked the
boundary so the on-disk operations are byte-identical to
spaceDB conv-abc entries and the conjugation happens
entirely on the consumer side via
`C = invRealVectors^T * convLattice`.  Identity shortcut in
`full` mode (M_loaded == M_conv => C = I => copy through)
keeps the common case branch-free.

- [x] C52a. Revise DESIGN 2.7 to make conv-abc the on-disk
  form; document M_conv / CONV_LATTICE / CELL_MODE, the
  consumer-side C and identity shortcut, and a
  "Diagnostic history" paragraph that records the
  Cartesian-xyz iteration this supersedes
- [x] C52b. Update ARCHITECTURE.md kpoints section: rename
  the on-disk arrays to `convAbcPointOps` /
  `convAbcFracTrans`; document the new CONV_LATTICE /
  CELL_MODE storage and the C-based transform with the
  full-mode identity shortcut
- [x] C52c. Rewrite PSEUDOCODE 4b end-to-end: 4b.1 writer
  near-passthrough (spaceDB ops emitted verbatim plus new
  metadata blocks), new 4b.2 reader additions, 4b.3
  consumer-side conjugation with C and identity shortcut
- [x] C52d. makeinput.py: drop `_to_cartesian_ops` and its
  plain-Python 3x3 helpers (`_matmul_3x3`, `_matvec_3x3`,
  `_transpose_3x3`, `_inv_3x3`); pass `conv_lattice` and
  `cell_mode` through `_write_mesh_kp_file` /
  `_write_density_kp_file`; emit the new `CONV_LATTICE` /
  `CELL_MODE` blocks in `_write_point_ops_block`.
  `POINT_OPS` values are now byte-identical to
  `share/spaceDB/<sg>`
- [x] C52e. kpoints.f90: rename `xyzPointOps` ->
  `convAbcPointOps` and `xyzFracTrans` ->
  `convAbcFracTrans` (reverting C31d's direction); add
  module-level `convLattice(3,3)` and `cellMode`; parse
  the new blocks in the style-1 and style-2 branches of
  `readKPoints`; default `cellMode='full'` and
  `convLattice=realVectors` in the style-0 trivial-
  identity setup; rewrite `computeRealPointOps` and
  `computeRecipPointOps` around
  `C = invRealVectors^T * convLattice` with the full-
  mode identity-shortcut copy path; add a private
  `invert3x3` helper for `C^{-1}`
- [x] C52f. Validate: clean build, full/prim agreement on
  decomposition properties (the physics test for the
  shortcut + conjugation paths producing identical
  results).  Confirmed by user 2026-05-18.

### Phase G -- UFF bond parameter database (DESIGN 4)

- [x] C32. Create bond_parameters.dat with UFF per-element
  parameters for Z=1-54 (DESIGN 4.4)
- [x] C33. Rewrite BondData class in condense.py to read
  bond_parameters.dat and compute K_ij, r_ij via
  get_bond_params() (DESIGN 4.6 items 1-3)
- [x] C34. Add bond_parameter_scale to condenserc.py,
  ScriptSettings, Condense.__init__, and
  parse_input_file (DESIGN 4.5, 4.6 item 5)
- [x] C35. Replace linear bond scans in create_lammps_files
  and normalize_types with get_bond_params() calls,
  applying bond_parameter_scale (DESIGN 4.6 items 3-5)
- [x] C36. Update CMakeLists.txt to install
  bond_parameters.dat instead of bonds.dat
  (DESIGN 4.7)
- [ ] C37. Validate: run a condense.py job end-to-end and
  verify that LAMMPS Bond Coeffs contain realistic
  UFF force constants and bond lengths

### Phase H -- Geometry-derived angle parameters (DESIGN 4.8)

#### Step 1: Replace angles.dat with clustering + computed K

- [x] D5. Design geometry-derived angle parameter system:
  cluster observed angles by triplet, compute K from UFF
  bond stiffnesses, add angle_parameter_scale and
  angle_cluster_tolerance keywords (DESIGN 4.8.1-4.8.6)
- [x] C38. Add angle_stiffness_coeff (default 0.15),
  angle_parameter_scale (default 1.0), and
  angle_cluster_tolerance (default 5.0) to
  Condense.__init__ and parse_input_file (DESIGN
  4.8.4-4.8.6, 4.8.8 item 5).  Follows the
  bond_parameter_scale precedent: force-field params
  live in Condense.__init__ with condense.in override,
  not in condenserc.py or ScriptSettings.  Landed in a
  prior session; verified present at lines 767-780
  (defaults) and 917-947 (parse_input_file).
- [x] C39. Implement angle clustering in
  create_lammps_files: collect (Z1, Zv, Z2, theta_obs)
  tuples, cluster within tolerance, compute K_angle from
  get_bond_params(), replace angles.dat lookup
  (DESIGN 4.8.3, 4.8.4, 4.8.8 items 1-3).  Landed in
  commit d51da45 along with a condense.py-wide cleanup
  of cryptic 1-2 letter locals and removal of 6 dead
  variables.  Subtasks C39.1-C39.5 all resolved by the
  same commit (see checkboxes below).

  Plan (captured 2026-04-16 — resolved in commit
  d51da45):

  - [x] C39.1. Add two shared helper functions.  Original
    plan proposed `_cluster_angles` and `_compute_angle_k`
    as methods on Condense; the actually-shipped form (commit
    74805ed) promoted them to module-level functions in
    src/scripts/angle_utils.py to avoid forcing
    make_reactions.py to import the larger condense module:
    - cluster_angles(observations, tolerance): PSEUDOCODE 10a.
      Group observations by canonical (z1, zv, z2) triplet
      (with z1 <= z2), sort each group by observed theta,
      greedy-merge while BOTH |theta - running_mean| <=
      tolerance AND the resulting cluster span stays within
      2 * tolerance (spread cap, consistent with 10e's
      cross-source step).  Return (angle_types,
      angle_type_map), where each angle_type carries z1,
      zv, z2, theta_0, obs_count, and a representative
      base_tag copied from the first observation in the
      cluster.
    - get_angle_k(z1, zv, z2, stiffness, scale,
      get_bond_params): PSEUDOCODE 10b.  Return
      stiffness * sqrt(K_arm1 * K_arm2) * scale, pulling
      K_arm1/K_arm2 from the injected get_bond_params
      callable (bound to Condense.bond_data.get_bond_params
      at the condense.py call site).

  - [x] C39.2. Refactor the angle section of the atom loop
    in create_lammps_files (current lines ~1521-1637) into
    a collect-only pass.  Build each observation dict
    {'z_trip', 'theta_obs', 'base_tag', 'vertex_atom',
    'end_atom_1', 'end_atom_2'} and append to a new
    angle_observations list declared alongside bond_count
    and angle_count near line 1390.  Remove the
    angles.dat lookup, the tag-uniqueness dedup, and the
    in-loop type_id assignment.

  - [x] C39.3. After the atom loop completes, add a
    clustering block that:
    - Calls cluster_angles(angle_observations,
      self.angle_cluster_tolerance) from angle_utils.
    - Sets angle_count = len(angle_observations) and
      num_local_angle_types = len(angle_types).
    - Builds local_angle_tags[t] as the string
      "{base_tag} {theta_0:.4f} {t}" so the final two
      tokens ({theta_0} {t}) match the tag-tail slot
      that normalize_types expects (10c Phase 3, 10d
      Phase 3, 10f Phase B all share this format).
    - Builds local_angle_coeffs[t] = [None, K_angle,
      theta_0] via get_angle_k(...,
      self.bond_data.get_bond_params).
    - Walks angle_observations in collection order and
      populates angle_bonded_atoms and
      ordered_angle_type with local type ids (which
      normalize_types remaps to global ids in 10f
      Phase A) so the per-atom ordering the downstream
      LAMMPS writer depends on is preserved exactly.
    - Exports per-local-type obs_count (slot 5 of
      each angle_types entry from 10a) so 10e can
      weight cross-source clustering by observation
      population.

  - [x] C39.4. Delete the now-unused local alias
    ``ad = self.angle_data`` inside create_lammps_files
    (around line 1338).  The AngleData class itself, and
    the identical alias in normalize_types, stay until
    C41 — normalize_types is updated in C40.

  - [x] C39.5. Sequencing constraint (not a side effect
    to live with): C39 and C39a must land together, and
    the precursor reaction-template DB must be rebuilt
    before any end-to-end test.  Until both producers
    (create_lammps_files in condense.py and the angle
    tag builder in make_reactions.py) emit the new
    "{theta_0} {t}" tag tail, normalize_types will see
    the same physical angle as two distinct types (at4
    string mismatch across lammps.dat vs templates) and
    bond/react type IDs will disagree, even with C40
    already applied.  The fix is to migrate both
    producers, not to paper over it with a hooke-id
    shim in create_lammps_files.
- [x] C39a. Port make_reactions.py's existing angle
  handling off angles.dat.  Remove the duplicated
  AngleData class (lines ~113-193) and the
  self.angle_data instantiation (line ~568).  In the
  _read_angle_data method (around line 2518), replace
  the hooke_angle_coeffs scan with the new collect /
  cluster(tolerance=0) / emit structure from
  PSEUDOCODE 10d so the emitted tag tail matches
  "{theta_0:.4f} {t}" in the same format C39 writes
  in create_lammps_files.  Cross-source unification
  of equivalent angle types is done centrally in
  normalize_types (C40), not by string comparison on
  the producer-emitted tag tail.

  Decision refined (2026-04-18): fixed tolerance = 0
  at the template producer; no K computation at the
  template producer.  The original 2026-04-17
  "hybrid with matching tolerance on both producers"
  plan was revisited once we traced the cross-source
  data flow in detail and realized:

  1. Local clustering at make_reactions.py is an
     optimization, not a correctness step.  Any
     tolerance T_m > condense.py's T_c causes silent
     wrong physics (distinct angles fused at the
     producer cannot be un-fused downstream).  Fixing
     T_m = 0 eliminates the hazard and makes
     reaction templates reusable across any
     condense.py run regardless of its T_c.  See
     DESIGN 4.8.10 for the full hazard analysis and
     rejection of the parameter-manifest alternative.

  2. Reaction template files carry no K values --
     only connectivity, per-atom angle entries, and
     the "{theta_0_local:.4f} {t}" tag tail.
     normalize_types recomputes K authoritatively
     from the triplet in 10f Phase C regardless of
     what any producer computed, so a producer-side
     K in make_reactions.py would be neither written
     nor consumed.  Skipping it keeps
     make_reactions.py independent of BondData and
     avoids plumbing angle_stiffness_coeff and
     angle_parameter_scale into a script that never
     writes their effect to disk.  DESIGN 4.8.8
     item 3 was amended to record this split:
     condense.py computes K, make_reactions.py does
     not.

  Net consequence: the C39a scope is smaller than
  the 2026-04-17 plan called for.  No parameter
  plumbing, no BondData wiring, no shared
  angle_cluster_tolerance handshake.  Just collect,
  cluster at tolerance=0, emit tags.

  Plan (captured 2026-04-17, refined 2026-04-18 --
  pending before writing code):

  - [x] C39a.1. Shared cluster_angles helper.
    Resolved in commit 74805ed: src/scripts/
    angle_utils.py provides AngleType NamedTuple +
    cluster_angles (PSEUDOCODE 10a) + get_angle_k
    (10b), with 10 unit tests in
    src/tests/test_angle_utils.py.  Both producers
    call cluster_angles from there.  make_reactions.py
    will call cluster_angles but NOT get_angle_k --
    get_angle_k is condense.py-only under the
    refined decision.

  - [x] C39a.2. In make_reactions.py, refactor the
    angle loop inside _read_angle_data (around line
    2518) into a collect-only pass per PSEUDOCODE 10d
    Phase 1.  Build each observation as the 8-tuple
    (z1, zv, z2, theta_obs, base_tag, a1, v, a2)
    with z1 <= z2 canonicalization -- same shape as
    condense.py create_lammps_files uses.  Remove
    the inline hooke_angle_coeffs scan, the in-loop
    tag-tail construction, and the existing
    angle_tag uniqueness dedup (cluster_angles does
    the dedup instead).

  - [x] C39a.3. After the template's angle loop
    completes, call cluster_angles(observations, 0.0)
    per PSEUDOCODE 10d Phase 2.  Tolerance = 0 means
    identity-only merge: observations with
    bit-identical theta values (possible because
    _read_angle_data rounds raw angles to 0.5-degree
    resolution before this point) collapse into one
    local type with obs_count > 1; non-identical
    observations each become their own local type.
    For each returned local cluster, build the tag
    tail "{theta_0:.4f} {t}" matching C39.3's
    format.  Do NOT compute or store K -- templates
    carry no angle coefficients (see PSEUDOCODE 10d
    Phase 3 and DESIGN 4.8.8 item 3).  Walk
    observations in collection order to populate
    the per-atom angle_bonded[] and angle_tag_id[]
    structures so downstream template writes see
    the ordering they expect.

  - [x] C39a.4. Delete the self.angle_data
    instantiation at line ~568 and the duplicated
    AngleData class at lines ~113-193.  The
    AngleData class itself lives on in condense.py
    until C41, which removes both copies and
    retires angles.dat.

  - [x] C39a.5. Caller update.  _read_angle_data's
    current return tuple includes unique_angle_tags
    and num_unique_angle_tags, which are no longer
    part of the new export (see PSEUDOCODE 10d
    export list).  Update the caller at
    make_reactions.py:~2171 to unpack only the
    fields normalize_types actually consumes:
    num_bond_angles, angle_bonded, angle_tag_id,
    local_angle_tags, per-local-type obs_count,
    num_angles_total, bond_angles_ext.  Any
    downstream code paths in make_reactions.py that
    read the old unique_angle_tags slot must be
    traced and either removed or updated to read
    local_angle_tags instead.

  - [x] C39a.6. Sequencing note (eased).  Under the
    refined decision, C39a does NOT have to land
    in lockstep with C40 -- the emitted tag-tail
    format already matches what C39 writes and
    what C40 will parse.  The pipeline is still
    inconsistent until C40 is live in
    normalize_types (cross-source clustering still
    uses the old hooke-id-suffix lookup that
    normalize_types currently does), but the C39a
    landing moment is flexible.  Strict ordering:
    C39a then C40 then C41 (teardown) then C42
    (validate) -- C39a in any ordering before C41
    is fine.
- [x] C40. Update normalize_types to own cross-
  source angle unification (DESIGN 4.8.8 item 4).
  Expands beyond the original "replace
  hooke_angle_coeffs scan with get_bond_params()"
  into the authoritative cross-file clustering step
  the hybrid approach requires.  Landed 2026-04-19
  alongside C39a: angle_utils.py gained
  cross_source_cluster + LocalRecord + FinalAngleType
  (with 9 new unit tests), and normalize_types now
  routes every angle through the collect-and-unify
  pipeline rather than the legacy hooke-id lookup.

  Plan (captured 2026-04-17 -- resolved 2026-04-19):

  - [x] C40.1. Collect all angle types emitted by
    every source: the lammps.dat produced by
    create_lammps_files and every reaction template
    produced by make_reactions.py.  Each incoming
    local_record carries (z1, zv, z2, theta_0_local,
    obs_count, base_tag, source, local_type_id) per
    PSEUDOCODE 10e's input layout.  obs_count is the
    observation-population weight used by the
    weighted running-mean merge; base_tag is the
    producer's representative tag prefix
    (element/species/molecule) that 10e carries
    through to the final type's canonical tag.

  - [x] C40.2. Group by canonical triplet (z1, zv,
    z2 with z1 <= z2).  Within each group, sort by
    theta_0_local and greedy-merge any pair within
    angle_cluster_tolerance of the running cluster
    mean.  The running mean is the final canonical
    theta_0 for the merged cluster.  Cap the total
    cluster spread (e.g., max theta span <=
    2 * tolerance) to avoid greedy chaining across
    a wide distribution.

  - [x] C40.3. For each final cluster, compute
    K_angle via get_angle_k() from angle_utils.py
    (PSEUDOCODE 10b; same helper create_lammps_files
    calls for its local step).  Apply
    angle_stiffness_coeff and angle_parameter_scale.

  - [x] C40.4. Rewrite tag tails and remap type IDs
    across every source file.  Each (source,
    source_type_id) pair maps to one final cluster
    id; the lammps.dat Angles section and every
    template angle reference is rewritten to the
    unified id.  Rewrite the tag tail to the final
    canonical theta_0 so any downstream tool that
    inspects the tag sees a consistent value.

  - [x] C40.5. Diagnostic (recommended): write a
    cluster-map file or log section listing, for
    every final cluster: final id, canonical
    theta_0, (z1, zv, z2), and the contributing
    (source, local_theta_0) pairs.  This is the
    main payback for the hybrid approach -- someone
    debugging a bond/react mismatch can open the
    file and see exactly which observations got
    merged where.
- [x] C41. Remove AngleData class and angles.dat from
  build (DESIGN 4.8.8 items 1, 6).  Landed 2026-04-19:
  AngleData class deleted from condense.py (the make_
  reactions.py copy went in C39a), self.angle_data
  init removed, init_env docstring updated, and
  angles.dat removed from the DATABASES list in
  src/data/CMakeLists.txt (per the bonds.dat
  precedent, the Perl AngleData.pm and the on-disk
  angles.dat file remain for historical reference but
  are not installed or read by any active code path).
  Runbook note: any pre-existing on-disk reaction
  templates carrying the old "{rest_angle}
  {hooke_id}" tail must be regenerated by the updated
  make_reactions.py before C42 validation.
- [ ] C42. Validate: run a condense.py job end-to-end
  and verify that LAMMPS Angle Coeffs contain
  physically reasonable force constants and rest angles
  derived from the input geometry

#### Step 2: Look-ahead angles in makeReactions (deferred)

- [ ] D6. Design angle creation in make_reactions.py
  post-reaction templates: port commented-out Perl
  addBondAngle, compute theta_0 from post-reaction
  coordinates, register new angle types
  (DESIGN 4.8.7, see "Empirical confirmation"
  paragraph for the 10-missing-angle decomposition).
  Concrete evidence first observed 2026-04-19 in
  jobs/molecules/b12/3_mol and re-confirmed
  2026-04-25 in jobs/molecules/b12/60_mol: the
  postRxn template has 230 angles (preRxn 240 minus
  the 10 that touched the deleted H's) but contains
  zero angles involving both new-bond endpoints
  (atoms 1 and 18 in the b12h12-b12h12 case).  Ten
  B-B-B angles are expected -- five at each vertex
  (e.g. 2-1-18, 3-1-18, ..., and the symmetric five
  around atom 18).  Without those angles the new
  inter-cage B-B bond is held only by its harmonic
  stretch -- the two cages can rotate freely about
  the bond axis with no restoring torque, and the
  bond axis itself can swing relative to either
  cage's local symmetry.  In an N-mer chain the
  deficit grows as 10N unconstrained angular DoFs
  at the joints.  The existing docstring on
  _build_phase_angles (src/scripts/make_reactions.py
  ~line 2675) already warns "the bonded molecules
  may be too floppy" because of this gap.  Not the
  cause of any reaction-trigger failure observed so
  far; deferred cleanly.

- [ ] C43. In make_reactions.py, the new B-B (or
  more generally new s1-s2) bond written into the
  postRxn template reuses the existing bond type
  for the element pair.  Confirmed by inspection
  for b12h12_1_b-1_b12h12_1_b-1: postRxn bond 61 is
  type 1 (the same B-B harmonic as every other B-B
  bond in the molecule).  That is correct for UFF-
  derived parameters where a new B-B is physically
  the same spring as any other B-B, but the choice
  should be revisited if/when a reaction wants to
  distinguish a "new" bond from a native one (e.g.
  for a different K value or equilibrium length at
  the reaction site).  Revisit alongside D6.

### Phase I -- Initial SCF potential database (DESIGN 5)

> **Documentation requirement (phase-wide).**  Every
> script created or modified in this phase must carry
> a module-level docstring that captures its purpose
> and role in the pipeline (library / producer /
> consumer per DESIGN 5.4) so future students reading
> the source can build the mental model from the file
> itself.  For `build_initial_potentials.py` the
> docstring must additionally include the manifest
> rationale and build analogy from DESIGN 5.7.  This
> requirement reflects the project's documentation
> policy in CLAUDE.md and is not optional.

- [x] C44. Implement src/scripts/initial_potential_db.py
  reader: PotentialEntry and ElementDatabase
  dataclasses; load() enforcing all six DESIGN 5.2
  validation rules (top-level and per-entry presence,
  schema_version == 1, element_symbol vs parent dir,
  coefficients/alphas length, label uniqueness,
  required "isolated" baseline); require_provenance()
  with Imago-source extra fields; lookup() and
  baseline() helpers (PSEUDOCODE 11.1).  Landed
  2026-05-18 alongside C45 and C46 as one commit.
- [x] C45. Implement initial_potential_db.save() --
  deterministic hand-formatted TOML emitter with fixed
  key ordering, %.16e floats, per-block = alignment,
  and multi-line array layout (PSEUDOCODE 11.2,
  DESIGN 5.5).  Landed 2026-05-18 alongside C44/C46.
- [x] C46. Add unit tests for initial_potential_db.py:
  each validation rule fires with the expected error
  message and field name; emitter is bit-deterministic
  for a fixed in-memory database; round-trip
  load(save(db)) preserves all numerical and
  provenance fields.  Landed 2026-05-18 in
  src/tests/test_initial_potential_db.py with 37
  tests covering all six validation rules, emitter
  alignment/format/escapes, %.16e bit-exact
  round-trip, and save idempotency.
- [x] C47. Add -pot CLI argument to makeinput.py and
  integrate the augmented per-element potential database
  per PSEUDOCODE 11.3.0 (reduced flow): load augmented
  per-element database when present; pick the entry via
  -pot LABEL override (FATAL on a missing label) else the
  default-tagged entry; fall back to legacy pot1/coeff1
  when no augmented file exists; emit the potential in the
  current on-the-wire format.

  Done 2026-05-20.  NOTE the original task text above said
  "fall back to 'isolated' on missing label (with
  warning)" -- that was the stale Phase-1 wording.  The
  reconciled spec (PSEUDOCODE 11.3.0, written this session)
  and the programmer's decisions superseded it:
  * Activation: the augmented db is used whenever
    share/atomicPDB/<elem>/s_gaussian_pot.toml is present
    (no files exist yet, so behavior is unchanged until
    curation in C49).
  * -pot LABEL missing from a present db is FATAL (not a
    silent fall-back).  -pot requested for an element with
    NO db warns and uses legacy.
  * Precedence per (element, species): -subpot target wins
    (legacy pot<N>); else augmented db; else legacy pot1.
  Implementation: rather than branch the imago.dat writer
  / Imago's scfV.dat reader, the augmented path
  materializes legacy-format pot/coeff files in
  .inputTemp/ from the chosen PotentialEntry
  (_write_legacy_pot_files_from_entry) and the reduced
  entry pick lives in _select_augmented_pot_entry; both in
  makeinput.py, both behind _obtain_pot_info.  Imago's
  scfV reader consumes only coeff column 1, so the
  generated term lines carry 5 columns with cols 2-5
  harmless (alpha + zeros).  6 new helper tests in
  src/tests/test_makeinput_pot.py; smoke-tested the full
  _obtain_pot_info wiring (default pick, -pot override,
  legacy fallback, -subpot precedence).  Library load
  passes known_methods=None (rule 9 registry arrives with
  C54).
- [ ] C48. Implement src/scripts/build_initial_potentials.py
  per PSEUDOCODE 11.4.  Split into sub-steps:
  - [x] C48.1. Manifest reader: load_manifest_v2 +
    dataclasses (CurationManifest / ReferenceSolid /
    ReferenceEntry / ManifestFingerprint); DESIGN 5.7
    rules 1-8, rule 9 gated on known_methods; missing
    manifest fatal.  Done 2026-05-20 (commit abcccae).
  - [x] C48.2. "isolated" baseline refresh from current
    pot1/coeff1: _parse_pot_file / _parse_coeff_file,
    build_isolated_entry, is_isolated_default_for,
    element_path, refresh_isolated_entries (step 1),
    save_databases (step 3).  Working producer->consumer
    loop with no SCF.  Done 2026-05-20 (commit abcccae).
  - [ ] C48.3. RE-SCOPED by ARCH §9: the producer no
    longer runs SCFs itself.  It becomes a kaleidoscope
    client -- hands kaleidoscope the curated structures
    plus makeinput options, lets it run and track the
    batch, then harvests converged potentials from the
    run dirs (converged scfV matches input scfV; alphas
    from input min/max/number; coeffs and alphas taken
    together).  Carries the C69 doc revisions.
    DEPENDENCY UPDATE 2026-05-28: the original "BLOCKED on
    a usable kaleidoscope slice (C63 + C68)" is satisfied
    -- C63 done, C68(a)+(b)+(b-cont) done and validated
    against real SlurmProvider; only C68(c) lost-vs-failed
    remains and is off the critical path.  NEW dependencies
    in priority order: C69 (DESIGN 5.7 / PSEUDOCODE 11.4 /
    ARCH 8.5 revisions for producer-delegates-SCF); C70
    (guidance_db.py library); C71 (kaleidoscope
    flight-builder helper consuming a guidance entry);
    and -- to actually deliver the acceleration -- C75
    (seed run populating share/historicalGuidanceDB/entries/
    so the producer's predict() calls return useful priors
    rather than the wide-grid fallback).  Without C75 the
    producer still works correctly; it just runs every
    reference solid through the wide-grid sweep, doing as
    much work as if predict-then-verify were not wired in.
- [ ] C49. Curate the first reference solid (e.g.,
  Au fcc) and add its manifest entry; populate the
  first non-trivial "default_solid" potential at
  share/atomicPDB/au/s_gaussian_pot.toml.  Establishes
  the curation pattern for subsequent elements
- [ ] C50. Implement src/scripts/bench_initial_potential.py
  per PSEUDOCODE 11.5: load benchmark manifest with
  held-out sanity check; run each test system under
  both -pot isolated and -pot default_solid; emit
  share/curation/bench_report.md with per-system
  iteration counts, mean_pct, held-out mean, and
  PASS/FAIL verdict; exit 0/1 gated on mean_pct >= 20.0
- [ ] C51. Validate end-to-end: regenerate the
  database via build_initial_potentials.py, run
  bench_initial_potential.py, and confirm mean_pct
  >= 20% reduction per VISION Principle 7

#### Phase 2 -- fingerprint-driven manifest selection

Builds on Phase 1 (C47-C51).  Phase 1 must be at
least functionally complete (literal-label selection
working end-to-end) before this chain lands, because
the matcher protocol and the bootstrap subprocess
both invoke makeinput's Phase-1 emission path
internally.  The chain may be developed in parallel
with Phase 1 but cannot ship until Phase 1 has
shipped.

- [x] C53. initial_potential_db.py: bump schema
  validation from v1 to v2.  Add `default: bool` to
  `PotentialEntry`; add `FingerprintRecord` dataclass
  with `method` / `sub_spec` / `payload`; add
  `fingerprints` list to `PotentialEntry`; add
  `default_entry(db)` and `find_fingerprint(entry,
  method, sub_spec)` public functions; extend the
  reader to enforce the new rules 7/8/9 (DESIGN 5.2);
  extend the emitter to write the new fields
  deterministically (DESIGN 5.5 layout extended to the
  new blocks); update C44/C45/C46 tests for the v2
  schema.  Existing v1 files in `share/atomicPDB/` are
  out-of-tree at the moment (C47 not yet shipped); the
  v1-to-v2 migration is handled by C60 regenerating
  them through the producer.

  Done 2026-05-20.  Library now schema v2: reader
  rejects `schema_version != 2`, requires per-entry
  `default`, enforces rule 7 (exactly one default),
  parses `[[potential.fingerprint]]` records under
  rules 8 (per-entry `(method, sub_spec)` uniqueness
  via the new `canonicalize_sub_spec` freeze) and 9
  (`method` in the optional `known_methods` registry
  arg, skipped when None).  Emitter slots `default`
  after `label` and emits fingerprint sub-blocks
  (`method`, `sub_spec` inline table, payload float-
  vectors as multi-line arrays) via
  `_emit_fingerprint_block` + `_format_inline_table`.
  Added `default_entry` and `find_fingerprint`.  Tests
  grew 37 -> 62 (rules 7/8/9, default_entry,
  find_fingerprint, canonicalize_sub_spec, fingerprint
  emit + round-trip); all green.  Confirmed no on-disk
  augmented db files exist yet, so nothing to migrate.
  Doc reconciliation done in the same session:
  PSEUDOCODE 11.3.0 "Reduced Flow" added (the
  no-environment-matcher path C47 will target), and a
  stale memory claim that "v1 files still load" was
  corrected against DESIGN 5.2's no-v1-compatibility
  rule.
- [ ] C54. makeinput.py: introduce the Matcher
  protocol (ARCHITECTURE 8.9) and the `MATCHERS`
  registry.  Refactor the existing `group_reduce`
  algorithm into `ReduceMatcher` (with
  `needs_loen_run = false`); preserve current
  behavior of `-reduce` as the regression target.
  This lands before the new schemes so the existing
  test surface validates the refactor.
- [ ] C55. makeinput.py: add `BispecMatcher` class
  with `needs_loen_run = true`.  Implement
  `to_loen_input(sub_spec)` returning the full LOEN
  parameter dict per DESIGN 5.10.5 (`loenCode`,
  `twoj1`, `twoj2`, `max_neigh`, `cutoff`,
  `angleSqueeze`; required keys `twoj1`/`twoj2`, the
  others with documented defaults matching today's
  hardcoded values);
  `parse_loen_output(path, sub_spec)` returning the
  per-site real-vector rows of `fort.21`; a
  `distance` metric over those vectors; and
  `representative(members)` returning the
  element-wise mean (ARCHITECTURE 8.9).  Reject the
  optional `by_element` sub_spec key with a clear
  "not yet implemented" error (placeholder for C62).
- [ ] C56. makeinput.py: add `name=NAME` keyword to
  the `-target` and `-block` argparse handlers.
  Validate uniqueness across the run; store the name
  in the existing spatial-flag records so subsequent
  scope references can resolve it.
- [ ] C57. makeinput.py: add `scope=NAME` and
  `scope=~NAME` keyword to the `-reduce` handler;
  introduce the new `-bispec` argparse handler with
  the same scope keyword plus `twoj1=` / `twoj2=`
  sub_spec fields.  Enforce the
  environment-scheme mutual exclusion (DESIGN 5.6.2):
  at most one `-reduce` *or* one `-bispec` per run.
- [ ] C58. makeinput.py: implement the nested-makeinput
  bootstrap of DESIGN 5.10.  Spawn a subprocess into
  `.inputTemp/loen_bootstrap/` with a stripped-down
  argument list (no `-pot` flag -- default tag flows
  through via the per-element preflight; Γ-only
  k-points; no grouping flags), invoke
  `imago.py -loen -scf no` against the produced
  `imago.dat`, and parse `fort.21` for the active
  matcher's per-atom fingerprint vectors.  Thread the
  matcher's `to_loen_input(sub_spec)` parameter dict
  (per DESIGN 5.10.5) into the existing
  `LOEN_INPUT_DATA` block emission in makeinput.py
  (currently hardcoded at lines 4659-4665), so the
  bootstrap's loen call uses the active matcher's
  parameters rather than the hardcoded defaults.  Add
  the `--no-loen-bootstrap` belt-and-suspenders flag
  to the nested call.
- [ ] C59. makeinput.py: implement the Phase-2
  species pass (DESIGN 5.6.4-5.6.7).  Walk
  `settings.methods` in order; dispatch each
  environment-based flag through its matcher with
  the resolved spatial scope; bucket atoms by
  fingerprint distance with the matcher's default
  similarity floor; pick one manifest entry per
  species via fingerprint match + similarity floor
  + default-tag fallback; propagate species choices
  through the type pass (XANES integration intact);
  emit per-type potentials into the Imago input
  file in today's on-the-wire format.
- [ ] C60. build_initial_potentials.py: bump the
  manifest reader to v2 (new `default` per entry,
  new `[[reference_solid.entry.fingerprint]]`
  declarations; new validation rules 7/8/9 per
  DESIGN 5.7).  Add the fingerprint-harvest step
  to the per-entry loop: invoke Python-side
  matchers in-process and Fortran-side matchers
  via cached `imago.py -loen -scf no` runs.
  Attach `FingerprintRecord`s and the `default`
  flag to each produced `PotentialEntry`.
- [ ] C61. Add end-to-end Phase-2 tests: a small
  reference structure runs through the bootstrap +
  loen + bucketing + entry-pick + emit chain
  producing a deterministic Imago input file.
  Include negative tests for the preflight coverage
  check (missing fingerprint family aborts cleanly),
  the mutual exclusion (`-reduce -bispec` together
  errors at parse time), and the similarity floor
  (sub-threshold match falls back to default with a
  warning).

#### Phase 2 follow-up -- element-aware bispectrum (parked)

- [ ] C62. Implement the element-aware bispectrum per
  D10.  Fortran side: add `bispecByElement` to
  `O_Input::readLoEnControl`, default false; extend
  `O_LocalEnv::computeBispectrumComponent` to
  accumulate per-neighbor-element when the flag is
  set; emit the extended `fort.21` format per D10's
  specification.  Python side: update
  `BispecMatcher.parse_loen_output` to return a list
  of `(neighbor_element, vector)` pairs when
  `sub_spec["by_element"]` is true; replace the
  "not yet implemented" rejection in
  `BispecMatcher.to_loen_input` (added under C55)
  with the real handling; update the matcher's
  `distance` per D10's semantics.  Producer side:
  update `build_initial_potentials.py` to harvest
  the element-aware variant when declared.  Library
  side: extend the database schema and reader if
  D10's payload shape requires it.  Parked pending
  both (a) Phase-2 base chain (C53-C61) landing
  first and (b) D10 design landing.

#### CLI override cleanup (DESIGN 5.6.1)

- [ ] C85. Make `-pot` scoped and repeatable in
  `makeinput.py` per DESIGN 5.6.1 / 5.6.5.  Change the
  argparse definition from a single `type=str` to a
  repeatable `action="append"` accepting `LABEL` plus an
  optional `scope=SPEC` / `scope=~SPEC` keyword, where
  `SPEC` is an element (`si`) or species (`si1`) tag.
  Thread the scoped overrides into
  `_select_augmented_pot_entry` / `_obtain_pot_info` so the
  per-element pick becomes a per-(element, species) pick:
  an unscoped `-pot` is the global default, a scoped `-pot`
  applies only to its element/species, and more-specific
  scopes win over broader ones (species > element >
  global), with two equally-specific scopes on the same
  species a parse-time error.  Keep the existing
  hard-error-on-missing-label and legacy-fallback-with-
  warning semantics.  Update the OPTIONS EXPLANATIONS and
  DEFAULTS help text and the selection tests.  Closes the
  gap that today `-pot` can only target the augmented
  database globally while only the legacy `-subpot` could
  target a single species.
- [ ] C86. Remove the `-subpot` option from `makeinput.py`
  (DESIGN 5.6.1 retirement note).  An audit of all 103
  installed element directories found only `pot1`/`coeff1`
  -- no `pot2`-or-higher ever shipped, so `-subpot` had
  nothing to substitute and was never exercised.  Its
  capability is subsumed by the augmented database plus the
  scoped `-pot` from C85.  Delete the argparse definition,
  the `pot_sub_*` parsing, and the path-1 substitution
  branch in `_obtain_pot_info`, collapsing its three-path
  precedence (subpot > augmented > legacy) to two
  (augmented > legacy).  Leave `-subbasis` in place: it
  deprecates together with `-subpot` only once the basis-
  set database gains an augmented, labeled form (a future
  task), preserving the long-standing subbasis/subpot
  symmetry.  Do C86 after C85 so the replacement capability
  lands before the legacy option is withdrawn.

#### Producer-derived entry labels (DESIGN 5.2.1)

- [x] C87 (DONE 2026-06-13). Assemble the augmented-database entry
  `label` at harvest instead of authoring it in the manifest, per
  the DESIGN 5.2.1 scheme
  `<reference_id>-<element><species>-t<type>-a<site>`.  The
  `type` (and `species`) numbers are not known until the
  grouping pass runs, so the label cannot be a manifest
  field; the harvest stage mints it after all executions
  finish.  Three coordinated changes:
  (a) co-opt the existing `datSkl.map` output of makeinput
      (`print_imago` / its sorting helper, makeinput.py:4627)
      rather than adding a new file.  It currently writes two
      columns (DAT# = sorted imago.dat atom number, SKELETON#
      = original imago.skl atom number); add three columns for
      each site's element, `atom_species_id`, and
      `atom_type_id`.  The sorting helper already returns
      `sorted_elem_id` / `sorted_spec_id` / `sorted_type_id`
      at the write point (line 4634), so the data is in hand;
      this just records the verdict so the harvester need not
      re-parse the run's input.  Update the `datSkl.map`
      docstring (makeinput.py:3373) and any reader of the
      file's column layout.
  (b) `build_initial_potentials.py` harvest path
      (`extract_potential` / its caller near line 1158) reads
      `datSkl.map`, finds the SKELETON# row equal to the
      harvested `atom_site`, takes that row's
      `(species, type)`, and builds the label as
      `f"{reference_id}-{element}{species}-t{type}-a{site}"`,
      lowercased, then stores the entry under it (the
      replace-by-label upsert keys on this assembled label).
  (c) `load_manifest_v2` relaxes manifest rule 3 so the entry
      `label` field is OPTIONAL -- when present it overrides
      the derived default (curator escape hatch); when absent
      the producer derives it per (b).  Add a label-safe
      charset check on `reference_id` (rule 5): lowercase
      letters, digits, `-`, `_`; no spaces, since the whole
      assembled label is typed into `-pot`.
  Update the manifest schema docs and the producer tests
  (datSkl.map round-trip + derived-label assembly + optional-
  override path).
  IMPLEMENTED: (a) makeinput `_sort_atoms` writes the five-column
  `datSkl.map` (DAT# SKELETON# ELEMENT SPECIES TYPE); imago.py
  `project_home_outputs` registers it as `outputs["datSkl_map"]`
  (FileNames.dat_skl_map token).  (b) producer
  `read_site_identity_map` + `assemble_entry_label`; the harvest
  loop derives the label when the entry pins none (new injectable
  `identity_fn`).  (c) `load_manifest_v2`: rule 3 no longer
  requires `label`; rule 5 adds the label-safe `reference_id`
  charset check; rule 6 guards derived-label collisions on
  `(reference_id, element, atom_site)`.  Docs: DESIGN 5.7 rules
  3/5/6 + per-entry field list; structure_control reader docstring.
  Tests: 637 pass.  The C74 smoke-test manifest still carries the
  explicit `si_diamond-si1-t1-a1`; it can now be dropped to exercise
  the derive path (label is optional).

### Phase J -- kaleidoscope flight infrastructure (VISION 4, ARCH 9)

The shared submit / track / harvest infrastructure.
Its first client is the C48 potential-DB producer
(C48.3 is blocked on a usable slice of this phase),
but it is general -- convergence sweeps, the C50 bench
harness, and future AIMD / high-throughput screening
are all clients.  Each task wants its DESIGN (D11-D14)
and PSEUDOCODE landed before code.

- [x] C63. Refactor imago.py to expose the callable API
  (ARCH 9.2, D11): the CLI becomes a thin wrapper;
  support both entry modes (prepared dir; structure +
  options); preserve checkpoint and lock behavior.
  Foundation for C65 and C68.

  Done 2026-05-21 (core; run_structure wiring folded into
  C68 per programmer decision).  Implemented PSEUDOCODE §12
  in imago.py -- RunStatus /
  ImagoError / JobIdentity / ImagoResult (12.1); the
  ScriptSettings split into from_command_line /
  from_options sharing reconcile(), with parse_command_line
  taking optional argv (12.3); project_home_outputs (12.2);
  the harvest helpers _read_scf_threshold /
  _last_data_row / _harvest_result reading the iteration
  file's last row for verdict/energy/count (12.5); the
  reentrant _run_core with the per-run-dir lock, cwd
  restore in finally, and the SystemExit->FAILED /
  ImagoError contract boundary (12.4); run_prepared; and
  main() rewritten as the thin CLI wrapper.  15 new tests
  in src/tests/test_imago_api.py; full unit suite green
  (224 passed).
  DEFERRED -- run_structure (structure-and-options mode) is
  an explicit stub that raises: it needs a makeinput
  "build a run dir" API that does not exist yet, which
  lands with C64 (ASE-free StructureControl factory) and
  C68 (kaleidoscope's makeinput->API dispatch, ARCH 9.4).
  Proposal: fold the run_structure wiring into C68 and
  treat the rest of C63 as done.
  Two intentional changes to record: (1) the CLI now exits
  non-zero on SCF NON-convergence (the old driver checked
  only fort.2 = ran-without-abort and exited 0); (2)
  reused_checkpoint / SKIPPED are not yet surfaced (the
  Fortran's within-run-dir checkpointing is preserved
  unchanged, but no Python-visible completion marker exists
  to report it), so they stay False/absent for now.
- [ ] C64. Add the ASE-free StructureControl factory to
  structure_control.py (ARCH 9.3, D12): build a
  StructureControl from (lattice, fractional coords,
  element symbols).  No ASE import; shared by C65 and
  C67.
- [ ] C65. Implement ase_imago.py ImagoCalculator (ARCH
  9.3, D12): calls the C63 API; the Atoms-reading glue
  uses the C64 factory.  Stays a flat module.
- [ ] C66. Implement cod_fish.py (ARCH 9.5, D14): COD
  structure fetch via urllib at a pinned cod_revision;
  strict on failure (no silent revision substitution).
- [ ] C67. Implement cif2skl.py (ARCH 9.5, D14): read
  CIF with ASE -> C64 factory -> skl write.
- [ ] C68. Implement kaleidoscope/ (ARCH 9.4, 9.6, D13):
  Parsl dispatch, the pluggable wingbeat seam, status
  tracking (complete-and-report), the flight
  workspace, and the run-reuse cache mechanism.  Also
  carries the C63 deferral: the run_structure ->
  makeinput "build a run dir" wiring (ARCH 9.4).

  Increment 1, 2026-05-21 (package core + both executors):
  Created src/scripts/kaleidoscope/ implementing PSEUDOCODE
  §13 -- model.py (KeyFields/KeyFile/CalcUnit/Flight/
  WingbeatOutcome/ReportEntry/FlightReport + KaleidoscopeError);
  workspace.py (slug rule, unit_run_dir, <calc> derivation,
  validate_flight, status.toml read/merge-write,
  serialize_flight); cache.py (write_cache_key,
  cache_key_matches with verbatim scalar compare + key-file
  byte-compare, is_cache_hit); wingbeats.py (Wingbeat base,
  ImagoWingbeat mapping ImagoResult + persisting result.toml,
  WINGBEATS registry); dispatch.py (dispatch +
  dispatch_unit + module-level _execute_wingbeat_task + collect,
  the per-future capture, LocalExecutor and a real
  ParslExecutor -- the programmer installed parsl
  2026.05.18, so the Parsl path is implemented and tested
  against a ThreadPoolExecutor Config).  Executor chosen by
  flight.parsl_config (present -> Parsl; absent -> local).
  Tests in src/tests/test_kaleidoscope.py with a fake
  wingbeat (no Imago binary): validate/slug/collision, cache
  hit/miss/byte-compare, dispatch under BOTH executors,
  complete-and-report (one failure does not abort), status
  lifecycle, and ImagoWingbeat mapping.
  Increment 2, 2026-05-21 (item (a): makeinput callable
  build API + run_structure wiring, per D15/P8).  Refactored
  makeinput.py to the §14 shape: added MakeinputError;
  ScriptSettings.__init__ now loads rc defaults only, with
  from_command_line / from_options classmethods +
  _args_from_options (reuses _build_parser().parse_args([])
  for defaults, raises on unknown keys); split
  parse_command_line into _build_parser + parse(argv);
  record_clp takes argv and is CLI-only; factored main()'s
  body into build_inputs(settings, sc); added
  build_run_dir(structure, options, wingbeat_dir, settings=None)
  with skl staging + chdir/restore-in-finally; rewrote main()
  as the thin CLI wrapper (sys.exit(main())).  Wired
  imago.run_structure to build_run_dir -> run_prepared
  (DESIGN 6.3.6).  Tests: test_makeinput_build_api.py (new;
  from_options/_args_from_options, build_run_dir staging +
  cwd discipline on success AND failure, record_clp,
  contract-fault raises -- build_inputs monkeypatched so no
  binaries needed); flipped test_imago_api's stub test into a
  build-then-run delegation test; updated test_makeinput_pot's
  missing-label test (SystemExit -> MakeinputError).  All 55
  affected tests pass.
  WORKER-SAFETY NOTE (refine landed 2026-05-21): 6.3.5
  originally named only the _load_rc sys.exit, but the build
  path held two more (an unsupported reduce op; a -pot override
  naming an absent db entry), each able to kill a worker on one
  bad unit because SystemExit bypasses the dispatcher's
  `except Exception`.  All three were converted to
  MakeinputError.  The refine generalized DESIGN 6.3.5 to "no
  sys.exit may remain on the build path" and tightened 6.3.1's
  fault taxonomy to make the SystemExit-vs-Exception mechanism
  explicit.  The residual helper-module audit came back clean:
  structure_control / initial_potential_db / element_data
  (every in-process module the build reaches) have no sys.exit;
  subprocess execs are exempt.
  C68(b) progress, 2026-05-21 (local multi-process level
  validated; recorded here, no committed test per programmer
  call).  Three cross-process checks on parsl 2026.05.18:
  (1) every arg ParslExecutor ships to a worker (CalcUnit with
  nested KeyFields/KeyFile, wingbeat_dir, default_wingbeat,
  _execute_wingbeat_task, WingbeatOutcome) pickles and round-trips with
  equality; (2) a fresh interpreter auto-registers ['imago'] on
  `import kaleidoscope` (DESIGN 6.2.2 holds cross-process);
  (3) a 4-unit flight run through the REAL dispatch /
  ParslExecutor path on HighThroughputExecutor + LocalProvider
  (max_workers_per_node=2) came back all done+converged across
  TWO distinct worker PIDs, both != main -- proving genuine
  process separation, status written by the worker and
  collected by main, and dfk cleanup.  Finding: a custom wingbeat
  must have its module imported in the worker (PYTHONPATH +
  sitecustomize/worker_init); the default 'imago' auto-registers
  so the real producer path needs nothing extra.  Note the HTEX
  worker-count knob in this parsl is `max_workers_per_node`
  (not `max_workers`); `worker_init` lives on the provider.
  C68(b-cont) DONE -- real SlurmProvider VALIDATED 2026-05-22.
  Hands-on 4-unit flight (silicon/diamond/graphite/silica) from
  jobs/kaleido_slurm/ via HTEX + SlurmProvider(rulisp-lab,
  nodes_per_block=1, cores_per_node=2, exclusive=False,
  max_workers_per_node=2). Proven: submission; worker_init env
  propagation (imago ran, HDF5 from the cpg conda env lib,
  makeinput built inputs); CROSS-NODE TCP (interchange on driver
  c159 <-> worker pool on c083); 2-worker parallelism;
  complete-and-report (2 failures isolated, 2 successes kept);
  status.toml lifecycle worker->main; CACHE HITS on re-run
  (silicon/diamond skipped, only misses re-dispatched). Verified
  worker_init: source conda.sh; conda activate cpg; export
  LD_LIBRARY_PATH=.../mamba/envs/cpg/lib; source the imago venv
  activate; THEN source imagorc (the venv activate clobbers
  PYTHONPATH, so imagorc must come after); OMP_NUM_THREADS=1.
  TWO defects surfaced (kaleidoscope isolated both as failed):
  (1) FIXED makeSGDB.py Perl->Python symlink escaping -- it
  shell-escaped os.symlink targets, leaving 314 broken spaceDB
  links (special-char space groups); now raw name -> 0 broken.
  Committed as bdb5b14. (2) FIXED buildAtomPerm Fortran
  STOP "no atom match found" for hexagonal/trigonal (gamma=120)
  graphite(186)+silica(152); cubic converged all along. Root
  cause was point-op conjugation done in the wrong basis;
  redesigned to T^-1 R_c T with a reciprocal inverse-transpose,
  plus a repeating-decimal input gate.  All 4 demo structures
  now converge.  Committed as ff33f5c.
  REMAINING for C68: (c) lost-vs-failed fine-graining under a
  real worker loss (validate dispatch._is_lost vs this parsl's
  ManagerLost/WorkerLost/BadState names). Plus a diversity-of-
  options study (multi-node, workers/node, exclusive, launcher)
  per the user goal -- off critical path, later cluster session.
- [x] C69. Revise DESIGN 5.7 / PSEUDOCODE 11.4 /
  ARCHITECTURE 8.5 so the producer delegates SCF running
  to kaleidoscope (drops the bespoke run_imago_scf, COD
  fetch, and per-solid cache).  Pairs with C48.3.
  In-scope addition (2026-05-29 refine): the curation
  manifest's per-solid `system_type` field landed in §5.7
  at the same time, so C69's producer rewrite needs to
  thread it through (parse it from the manifest, supply
  it to the guidance-dataspace predictor as the
  `system_type` argument when calling predict(), record
  it on the produced potential-DB entry for forensics).
  DONE (2026-06-11): all three docs rewritten to the
  three-phase build/dispatch/harvest shape.  Settled
  decisions: COD materialization = **Option A** (thin
  `materialize_structure`: local read, or one-shot fetch
  of the pinned `cod_revision` to a plain local file,
  decoupled from any cache; the producer is NOT
  network-free and the strict COD-fetch error contract is
  preserved); `kpoint_spec.density` = optional curator
  override that pins/centres the grid.  Also: validation
  rule 2 now enforces `system_type` presence + the
  four-value domain; `system_type` recorded in the
  produced entry's provenance (not a new 5.2 field); the
  per-solid SCF cache, `is_cached_v2`, and `is_cached_loen`
  are deleted (kaleidoscope run-reuse cache subsumes them);
  Fortran-side fingerprints fold into the dispatched flight
  as structure-only `-loen` units and are read from the
  run dir's `fort.21` at harvest; the converged grid point
  dual-harvests into `historicalGuidanceDB/staging/` via
  C72 `harvest_flight`.  ARCH 9.8 open item marked
  RESOLVED.  Code is C74.

### Phase K -- historical guidance dataspace (VISION 5, ARCH 10, DESIGN 7)

The accumulation prong, post-Path-B rewrite.  Each task
wants P9 landed and the relevant DESIGN 7 subsection
consulted before code lands.  The library + predictor
(C70) is foundational; the flight-builder helper (C71)
is what makes predict-then-verify a real workflow; the
imago result.toml extension (C76) is a small Fortran-side
prerequisite for harvest; the harvest hook (C72) closes
the loop back into the dataspace; the curator helper (C73)
gates staging into canonical entries; the seed flight
(C75) trains the predictor over the chemistry surface;
the C48.3 wiring (C74) is the first major consumer.

- [x] C70. Implement `src/scripts/guidance_db.py` per
  DESIGN 7.2 / 7.4 / 7.5 / 7.6 and PSEUDOCODE P9(a)+(b).
  Two halves rolled into one library because they share
  the in-memory Dataspace:
  **I/O half.** GuidanceDataspaceError; dataclasses
  (Signature, Measured, Context, Verification, Provenance,
  GuidanceEntry, Dataspace, PredictionResult);
  `elemental_groups.toml` loader with element-to-group lookup;
  compute_signature(structure, system_type, group_table)
  building the 13-d composition_vector and the
  lattice_family one-hot; load(root) walking
  entries/<system_type>/ with all 12 validation rules
  yielding clear file/block/field error messages;
  save_entry(entry, root) using the deterministic hand-
  formatted emitter (16-sig-digit floats, fixed block
  sequence, multi-line composition vector layout,
  comma-trailing per-line array layout); the
  `<system_type>-<short_sha>` slug derivation.
  **Predictor half.** predict(dataspace, query, basis,
  functional) returning a PredictionResult; the
  system_type switch (canonical entry for non-crystalline;
  two-stage k-NN for crystalline); sub-model selection
  by (basis, functional) with the fallback chain;
  stage-1 distance d1 over composition + lattice_family
  with inverse-distance weights yielding predicted_gap,
  predicted_magnetization, confidence_1; stage-2 distance d2
  over predicted electronic character yielding
  predicted_kpoint_density and confidence_2; combined
  confidence; the is_under_trained flag whose semantics
  drive the flight-builder's wide-grid fallback.
  Tests: every validation rule fires with the expected
  message; emitter is bit-deterministic; round-trip
  load(save_entry(...)) preserves all fields; the k-NN
  distance and weight formulas produce expected scores
  on a handful of curated pairs; sub-model fallback chain
  exercised under sparse-data conditions; the canonical-
  non-crystalline path returns the expected canonical
  entry.
  Done 2026-05-29 in four tested increments: foundation
  (c7223b2: constants/dataclasses + elemental_groups.toml +
  load_elemental_groups + compute_signature/bravais),
  reader (b82d872: load + 12 rules), emitter (4593d0a:
  save_entry/format_entry, byte-deterministic), predictor
  (cddf68f: predict + two-stage k-NN).  84 tests in
  test_guidance_db.py, all green.  Three §15.2 details
  resolved against the real StructureControl: Bravais via the
  IT space_group_num ranges (hard error on the 0 = "no space
  group"); element symbols matched case-insensitively (SC
  stores them lower-cased).
- [x] C71. Implement the flight-builder helper inside
  `src/scripts/kaleidoscope/` per DESIGN 6.2.8 / 7.7 and
  PSEUDOCODE P9(c).  **Prerequisite (model catch-up):** the
  C68 code shipped `CalcUnit.calc` as `Optional[str]` and
  `Flight` with no `sweep`, but findings 2/3 (2026-05-28
  refine) moved the design to `calc: tuple[str, ...]` (one
  directory component per varied sweep axis) plus a
  `Flight.sweep: SweepRecord | None` field (DESIGN 6.2.1,
  PSEUDOCODE 13.1).  Land that model change first --
  `model.py` (calc -> tuple, add SweepRecord + Flight.sweep,
  plus a generic opaque `Flight.metadata: dict` field per
  PSEUDOCODE 15.6 -- the seam by which the flight-builder
  attaches its PredictionRecord without the dispatch core
  interpreting it, Principle 9),
  `workspace.py` (unit_run_dir splats the tuple onto the path,
  validate_flight slugs each component and keys collisions
  on the tuple, serialize_flight emits calc as a TOML array,
  a `[flight.sweep]` block, and each `metadata[key]` as a
  verbatim `[flight.<key>]` table), and the test suite -- so
  build_calc_tag has a tuple-shaped target to populate and the
  PredictionRecord round-trips through flight.toml.  This
  also implies a one-line DESIGN 6.2.1 / PSEUDOCODE 13.1
  addition (the Flight.metadata field + its serialize line),
  flagged in PSEUDOCODE 15.6.
  predict_settings(structure, options,
  dataspace, system_type, basis, functional, verify, id,
  extra_axes) returning (Flight, PredictionRecord).
  build_verification_grid(center, confidence) lays out
  the logspace grid whose width and point count scale
  inversely with predictor confidence (the 7.7 starting
  heuristics).  is_under_trained falls through to
  default_wide_kpoint_density_grid() (the 8-point bracket
  list in DESIGN 7.9).  trust mode (verify=False)
  collapses to a length-1 grid at the predicted center.
  build_calc_tag(calc_axes) emits the tree-per-varied-
  axis paths per DESIGN 6.2.4.  Tests: prediction record
  shape; high-confidence returns a 3-point tight grid and
  low-confidence returns a 6-point wider grid; log-spacing
  centered on the predicted value; empty dataspace
  returns the wide-grid default with
  policy = "wide_grid_no_prior"; trust mode returns
  length 1; tag derivation matches the 6.2.4 examples.
  **DONE.** Part A model catch-up (72a67b6): calc->tuple,
  SweepRecord, Flight.sweep/metadata, toml_line arrays,
  serialize_flight [flight.sweep]/[flight.<key>] tables.
  Part B builder (8c37c42): organized as a
  `kaleidoscope/builders/` subpackage (one module per
  builder, anticipating XANES/basis-size siblings);
  `builders/predict_verify.py` with predict_settings +
  build_verification_grid + the wide/trust paths +
  build_calc_tag.  16 tests, physics layer monkeypatched
  (no $IMAGO_DATA); full src/tests 576 passed.  C70
  predictor confirmed present (the §15.5 "in progress" note
  was stale).  Unblocks C72.
  **Deferred to C74** (PSEUDOCODE-15.6 elisions resolved with
  documented defaults, none affecting C71's tested logic):
  (a) `imago_commit` is producer-injected via options, not
  yet sourced -- a build-identity concern shared with C78;
  (b) when `structure` is a pre-loaded StructureControl the
  caller must pass an explicit `id` and the key-file source
  falls back to `imago_skl`; (c) `predict_settings` takes an
  optional `root` kwarg (15.6 left the flight root to the
  caller).  Also (from C76, still open for C72): the
  predictor's spin character should key on
  total_magnetization, not spin_polarization.
- [x] C72. Implement `src/scripts/guidance_harvest.py`
  per DESIGN 7.8 (harvest half) and PSEUDOCODE P9(d).
  Walks each structure's verification sub-grid, parses
  each converged calc's result.toml for the measured
  electronic-structure quantities (gap_ev, gap_kind,
  total_magnetization) plus total_energy for the
  convergence test, picks the converged grid point per
  the two-sided delta-below-threshold rule (DESIGN 7.8
  step 3c), SKIPs and tags `prediction_mismatch = true`
  on non-convergence at the top, recovers
  predictor_confidence and predictor_neighbor_ids from
  the flight's [flight.prediction] block, builds a
  rich GuidanceEntry, and writes it to
  staging/<system_type>/ via save_entry().
  **DONE.**  Settled the three-source "Model 1" sourcing
  with the programmer (2026-05-30): flight.toml = the
  plan (units, sweep, prediction; each grid point's kpd
  decoded from its `kpt-density-<int>` calc tag, since
  options are not persisted); result.toml = per-run facts
  (total_energy, gap_ev, gap_kind, total_magnetization,
  and a NEW `scf_threshold` field added to ImagoResult +
  the wingbeat writer); the structure .skl = structural
  facts (cell_atom_count, cell_volume_per_formula_unit in
  Bohr^3, Z=1).  Conventions: `metric_threshold =
  scf_threshold`; `imago_commit` falls back to "unknown".
  `spin_polarization` is recorded as 0.0 (not measured) --
  the predictor's spin character was switched to the
  intensive magnetization `|M|/N_atoms` (guidance_db
  stage1/stage2; `predicted_spin_pol` renamed
  `predicted_magnetization` across guidance_db /
  predict_verify / DESIGN 7.6-7.7 / ARCH 10 / PSEUDOCODE
  15; AFM-blindness noted in DESIGN 7.10).  Added
  `read_flight_toml` + `flight_id_of` to
  kaleidoscope.workspace.  Tests: test_guidance_harvest.py
  (16, synthetic workspaces, no $IMAGO_DATA) covering the
  converged / two-sided-delta / mismatch / record-recovery
  / trust-mode / non-crystalline / no-sweep paths, plus
  read_flight_toml round-trip + 3 intensive-magnetization
  predictor tests.  Full src/tests 594 passed.  Registered
  guidance_harvest.py in scripts CMakeLists.
- [x] C73. Implement `src/scripts/guidance_promote.py`
  per DESIGN 7.8 (curator half).  Four modes:
  interactive review (default), `--auto-promote` (with
  the middle-60%-of-grid + top-three-energy-variance
  rule from DESIGN 7.8), `--all`, `--dry-run`.
  Interactive printed summary covers signature,
  measured, verification, provenance.  Tests exercise
  each mode against a synthetic staging directory.
  **DONE.**  Implemented faithfully to PSEUDOCODE 15.7
  (no design deviation): `auto_promote_ok` (the three
  conditions: mid-60%-of-grid convergence,
  top-three-energy population variance <
  metric_threshold*10, gap_ev/gap_kind consistency, all
  read from the staged file alone); `move_to_entries`
  (a rename, refusing a destination collision);
  `format_summary`; `_ask_choice` (p/s/d, empty ->
  SKIP); the `promote(db_root, mode, *, ask, output)`
  driver returning (entry_id, action) records with
  `ask`/`output` injected for testability.  Reuses
  `guidance_db.load_entry` for schema validation.
  Registered in scripts CMakeLists.  Tests:
  test_guidance_promote.py (14) -- the rule's accept +
  four reject paths, variance, move/collision, and all
  four driver modes against a save_entry-built staging
  dir.  Full src/tests 608 passed.
- [ ] C74. Wire the kaleidoscope flight-builder helper
  (C71) into the C48.3 producer.  Replaces the current
  "user picks settings up front" pattern with predict-
  then-verify: for each reference solid in the curation
  manifest, build_kpoint_convergence(...) returns a
  verification sub-grid Flight, kaleidoscope dispatches it,
  the harvest hook reads back both the converged potential
  (for the C48.3 deliverable) and the rich measured
  quantities (for guidance contribution).  Bundled with
  C48.3 once C70+C71+C72+C76 are in place.
  **BLOCKING OPEN DESIGN QUESTION -- RESOLVED 2026-06-12,
  Option B (per-record sub-model), single-source variant.**
  The (basis, functional, kpoint_integration) sub-model is
  carried ONLY on each per-structure PredictionRecord and is
  NOT duplicated into sweep.fixed_axes (which is now {} for
  builder flights) -- one home, no drift, no reader
  confusion.  The harvest reads the sub-model per-id from the
  record; a structure with NO record is SKIPPED (the record
  is the sole source of both system_type and the sub-model,
  so the old fixed_axes read AND the "default crystalline"
  fallback are both retired).  Propagated this session
  (DESIGN 6.2.8 steps 5/6 / 6.2.9 / 7.7 steps 5/6 / 7.8
  step 3 guard + 3f; PSEUDOCODE 15.6 dataclass + builder
  SweepRecord={} + record, 15.7 harvest skip + per-id reads).
  Matching CODE rides with C74 proper: add
  basis/functional/kpoint_integration to
  predict_verify.PredictionRecord; set them in the builder;
  make the builder's SweepRecord fixed_axes empty; switch
  guidance_harvest's context read onto the per-id record;
  add the no-record skip; drop the crystalline default.
  **Update the test fixtures regardless of size** (record-
  less synthetic flights must now attach a record or be
  recognized as skipped).
  Also clarified (same session) that a sweep-less / one-off
  calc is not blocked: it is a length-1 sweep, harvested for
  the potential deliverable (5.7) and skipped for guidance
  staging; guidance for a known density is seeded manually
  (DESIGN 7.8 intro / 7.9).
  **PROGRESS (2026-06-12): incr 1 (builder, 94ba551) + incr 2
  (harvest, 838878d) DONE -- the fully unit-testable layers,
  612 passed.**  Incr 1 = the #7 renames (predict_verify.py ->
  kpoint_convergence.py, predict_settings ->
  build_kpoint_convergence, predicted_value ->
  predicted_kpoint_density, + import sites), the per-record
  sub-model + empty fixed_axes, the per-id
  metadata["predictions"] stash, and the curator-override
  (center) path.  Incr 2 = the harvest per-id read + no-record
  skip + crystalline-default retired + #8 (gap required via
  _require_field) + len==1 single-point skip + the rts/rt
  rename.  **Incr 3 DONE (3a manifest 55376bb, 3b CalcUnit.kind
  349efc7, 3c producer pipeline 7016ff9):** load_manifest_v2
  rule 2 now requires system_type (4-value domain) +
  basis/functional/kpoint_integration on ReferenceSolid;
  CalcUnit.kind + the harvest kind=="convergence" filter; and
  the full three-phase producer (materialize_structure [Option
  A], make_producer_options, pick_converged_unit,
  extract_potential, the run-log/provenance helpers,
  curation_executor, the CLI) wired predict->dispatch->harvest
  with the guidance contribution.  The ipdb emitter gained
  dict/array provenance support (first Imago-source entry) +
  the system_type forensic extra.  imago_commit injection (C71
  deferral) landed via _git_sha.  **C74 is CODE-COMPLETE,
  full suite 630 passed** -- the toolchain seam (dispatch /
  extract_potential per-site scfV / COD fetch / force cache
  bypass) is injected + mocked in tests and needs a live Imago
  run on the cluster to validate end-to-end.  STILL DEFERRED:
  build_loen_units / harvest_fingerprints stub to [] pending the
  C54 matcher registry + C60 fingerprint harvest.
  **LIVE SMOKE RUN 2026-06-13 surfaced the makeinput/imago
  option-contract gap** (the validation C74 was waiting on).  The
  first live producer run aborted every unit with `unknown
  makeinput option: 'basis'`: the producer emits one options
  dictionary forwarded to both a strict makeinput and a lenient
  imago, but 6 of its 7 keys do not match makeinput's real argparse
  dests, and `basis` is an imago run-time selection (makeinput
  writes all bases into imago.dat), not a makeinput setting at all.
  **Resolved in DESIGN 6.2.10** ("option-contract seam") + a 6.3.6
  pointer; propagated to DESIGN 6.2.2 / 5.7 and PSEUDOCODE 11.4 /
  13.2 / 14.4.  Decisions: the **wingbeat** owns the split; the
  **producer** emits dest-keyed, coded options; SCF convergence
  threads like xccode via a new makeinput `-converg` dest; routing
  is by each tool's recognised-key set, so `basis` can migrate to
  makeinput later without reworking the seam.
  **Refinement (during the code work): the two-dictionary split.**
  Wiring (b) surfaced a conflict between 6.2.10 (producer emits
  dest-keyed options) and 6.2.8 (the builder reads the human
  basis/functional/kpoint_integration FROM options for the
  predictor + record): one shared dict cannot be both coded (for
  the tools) and human (for the builder).  Resolved with the
  programmer as a **two-dictionary** design -- `options` carries
  tool-facing coded keys only, and a separate `submodel` dict
  ({basis, functional, kpoint_integration}, the established name)
  carries the human physics names into `build_kpoint_convergence`.
  Only the basis appears in both channels (as `submodel["basis"]`
  and `scf_basis`), which is intrinsic and benign; `functional` /
  `kpoint_integration` carry different values in each and never
  collide.  Propagated to DESIGN 6.2.8 (signature + predict +
  record reads) / 6.2.9 (input-channel note) / 6.2.10 (decision-2
  "physics names feed the builder through their own channel") / 5.7
  (step 3b) and PSEUDOCODE 11.4 + 15.6.  **Follow-on CODE
  (rides with the live-validation step of C74):** (a) add makeinput
  `-converg` (dest `converg`, overrides rc `converg_main`); (b)
  rewrite `make_producer_options` to emit the dest-keyed vocabulary
  (functional->xccode, kpoint_integration->scfkpint, basis->
  scf_basis, scf_threshold->converg, shift->kpshift) as TOOL-FACING
  keys only, AND give `build_kpoint_convergence` a separate
  `submodel` dict arg ({basis, functional, kpoint_integration}) for
  the predictor + PredictionRecord, so the physics names never enter
  `options` (the two-dictionary split, DESIGN 6.2.8/6.2.9/6.2.10);
  the producer builds `submodel` from the ref and passes it; (c)
  export
  `imago.OPTION_KEYS` + a `CACHE_ONLY_KEYS` set from
  `kaleidoscope.wingbeats`; (d) move the
  partition into `ImagoWingbeat.run`, retire the shared-options
  `run_structure` call; (e) `_KEY_SCALAR_NAMES` ->
  `("converg", "imago_commit")`; (f) harden `pick_converged_unit`
  against a missing `result.toml` (treat a failed / result-less
  unit as non-converged, reading `status.toml`).  Reminder: the
  pinned `kpoint_spec.density` builds a 3-point tight grid, not a
  single point; re-run needs `--force` or a cleared
  `$IMAGO_DATA/curation/workspace/`.
- [ ] C75. Seed `share/historicalGuidanceDB/entries/`
  via a deliberate stratified seed flight.  ~150-250
  calculations covering the chemistry surface
  representatively rather than at random: for each
  pair of element groups (alkali x halide, TM x
  chalcogen, group_iv x chalcogen, etc.) and each common
  stoichiometry pattern (binary AB, A2B, ABO3
  perovskite), pick 3-5 representative COD entries.
  Avoids over-sampling Si compounds and under-sampling
  heavy elements + actinides + lanthanides.  Plus the
  three day-1 canonical entries seeded manually for
  amorphous, nanostructure, and molecular system_types
  per DESIGN 7.9.  Uses kaleidoscope wide-grid sweeps
  (no prior available); harvest via C72; promotes via
  C73 `--auto-promote` so the curator reviews only
  ~20% outliers.  This is the bootstrap that makes the
  predictor non-trivial for the crystalline subtree
  from day-2 on.  Needs cluster time and is a real
  sub-project of its own: stratified-sampling design,
  COD query scripting, allocation budget, post-seed
  calibration of the k-NN tuning knobs per DESIGN 7.10.
- [x] C76. Surface the electronic-structure quantities the
  guidance harvest (C72) needs by making the **iteration
  file the single primary read surface** (programmer
  guidance, 2026-05-29).  It already carries `total_energy`
  and -- spin-polarized runs only -- the magnetization
  column.  Add to the iteration data: `gap_ev` + `gap_kind`
  -- the raw none/direct/indirect electronic-character
  signal that feeds the two-stage predictor.  (dos_at_fermi
  was considered and dropped 2026-05-29: gap/gap_kind alone
  carry the character signal.)  Putting these in the
  iteration data means any
  plain SCF run yields them, so the harvest never has to
  decide whether to run `-scfdos`.  `imago.py` extends its
  iteration-file parser (the primary read), populating
  `result.toml` from it.  Each field optional in the
  guidance schema (a non-spin calc omits magnetization).
  Small Fortran-side change to the iteration-file writer +
  the imago.py parser (DESIGN 6.1).  Prerequisite for C72.
  NOTE: gap / spin depend on the k-point integration
  method, now recorded as Context
  `kpoint_integration` and part of the predictor sub-model
  key (DESIGN 7.2/7.6); the harvest fills it from the flight
  options.
  **Python half DONE (16156a4):** imago.py reads the fixed
  8-column row (col 6 Mag._Mom. -> total_magnetization, col 7
  gap Hartree -> gap_ev in eV, col 8 code -> gap_kind via
  {0:none,1:direct,2:indirect}), length-gated for pre-gap
  files; ImagoResult + the wingbeat result.toml serializer
  carry the fields.
  **Fortran half DONE (f44e6dd):** populate.F90 computes the raw
  gap (CBM-VBM, a.u.) + kind from the sorted spectrum (kpoint
  recovered via the existing 1+(index-1)/(numStates*spin)
  idiom); potentialUpdate.F90 writes the fixed 8-column row in
  both spin branches with explicit 1x field separators (no
  column collisions); imago.F90 emits a single 8-column header.
  Metal detection uses a dedicated `metalGapThresh = 1.0e-3`
  a.u. (~0.027 eV ~ kT, NOT smallThresh) so finite-mesh metal
  artifacts collapse to a zero-gap metal -- see DESIGN 6.1.
  gap_kind code order confirmed against GAP_KIND_BY_CODE
  (0 metal/1 direct/2 indirect).  Verified end to end: diamond
  -> 0.184 a.u. ~ 5.0 eV, kind 2 (indirect); aluminum -> 0.0,
  kind 0 (metal).  C76 COMPLETE; unblocks C72.

### Phase L -- resource & cost dataspace (VISION 6, ARCH 11, DESIGN 8)

Sibling of Phase K.  Near-term consumer is provisioning (C81);
config-optimization / build-comparison / scaling studies ride
on the same data later with no schema change.  Built on P10.

- [ ] C77. Implement src/scripts/resource_db.py: the library +
  predictor.  load/validate (12 rules, file/block/field error
  messages), the hardware_registry.toml loader + fingerprint
  recipe (8.5), the three registry validators, the
  deterministic hand-formatted emitter with the two-layer build
  record, round-trip load/save, and the physics-informed
  power-law predictor with k-NN fallback (8.6).  Modeled on
  guidance_db.py (C70) and initial_potential_db.py.  Tests cover
  the rules, emitter determinism, the regressor on curated
  points, registry rejection of unknown keys, and censored-bound
  handling.  Register in src/scripts/CMakeLists.txt.
- [ ] C78. CMake build-system hook emitting build_info.toml at
  configure/install time -- compiler + full flag string +
  detected HDF5 / ScaLAPACK / BLAS / MPI versions and variants
  -- so the build block (both layers) is captured without
  hand-entry (ARCH 11.4, 11.6; DESIGN 8.7).
- [ ] C79. Capture hooks in the wingbeat / imago.py: record the
  dispatch-time size signature + execution config into the run
  directory, and scrape SLURM sacct (MaxRSS / disk / Elapsed)
  after the run; optional imago self-report of per-phase timings
  into result.toml (ARCH 11.4, DESIGN 8.7).
- [ ] C80. Implement src/scripts/resource_harvest.py: walk a
  finished flight, join the four capture sources into one
  Observation per run dir, stage censored failed runs as bounds,
  write to staging/<fingerprint>/ (DESIGN 8.7).  Sibling of
  guidance_harvest.py (C72).
- [ ] C81. Provisioning consumer in the flight layer (the
  kaleidoscope flight-builder helper or a thin sibling): query
  the predictor with a proposed config + size, apply a safety
  margin, annotate the Parsl provider's SLURM resource request;
  cold-start fallback for empty fingerprints (DESIGN 8.6, 8.8).
  The near-term consumer (VISION Goal 6).
- [ ] C82. Seed the resource dataspace on the local cluster: a
  small manual + flight-harvested set of observations spanning
  the size range on each available hardware fingerprint, so the
  predictor becomes non-trivial.  Needs cluster time.  Sibling
  of the C75 convergence seed.
- [ ] C83. (future) src/scripts/resource_migrate.py: schema
  migration tool mirroring guidance_migrate.py.  Not day-1
  scope (ARCH 11.6).
- [ ] C84. Have Imago stamp its own build commit into its
  output so the guidance/resource harvests record a real
  build identity instead of the `"unknown"` fallback.
  Background: the C74 producer echoes `imago_commit` through
  the run options into the wingbeat-written `result.toml`
  (the near-term fix the harvests read), but that records
  what the *producer* believed it ran, which can drift from
  the binary actually executed.  The robust upgrade is for
  the running binary to report its own build commit (e.g.
  from the C78 `build_info.toml`, or a compiled-in version
  string), which the wingbeat then copies into `result.toml`
  in place of the echoed value.  Pairs with C78 (build
  identity) and C79 (wingbeat/imago.py capture hooks).

---

## TOOLING (lint helpers)

- [ ] T1. Improve `.claude/commands/scripts/rewrap_prose.py`
  heuristics so it stops over-aggressively breaking
  paragraphs.  Three known false positives observed
  while landing C39.1/C39a.1 on 2026-04-18:
    1. Lines ending in `:` always terminate a paragraph
       (lines 657-658, 683-684), which orphans the prose
       that flows into a colon-introduced equation or
       list ("...via the formula:" -> previous line
       ends mid-sentence with no following words).
    2. The `_LABEL_LINE_RE` heuristic
       (`^\s*\w[\w_\. \t]{2,25}:(\s|$)`) catches
       sentence-final phrases like "Default upstream:
       0.15", treating them as field labels and
       refusing to merge them with the surrounding
       prose.
    3. Hyphenated compounds split across lines
       (e.g. "re-clustering", "geometric-mean") are
       sometimes rejoined with an extra space ("re-
       clustering" -> "re- clustering") instead of
       being un-split.
  Each issue has a workaround at the source level
  (rephrase to avoid the trigger pattern), but the
  heuristics should ideally handle these cases so the
  source can read naturally.  See angle_utils.py and
  test_angle_utils.py for the prose patterns that
  initially tripped the script.
- [ ] T2. Consolidate the two parallel fixture trees onto
  `src/tests/fixtures/`, deleting the vestigial top-level
  `tests/fixtures/`.  The whole source tree lives under
  `src/` and the active conftest already reads the
  `src/tests/` copy, so the top-level copy is a drift
  footgun (the same fixture diverging in two places).
  Before deleting, grep for any test or pytest config that
  still points at the top-level `tests/` path and redirect
  it.  Do this as its own small cleanup commit, separate
  from the C74 producer work, so the dedupe is easy to
  review and revert independently.

---

## ARCHIVE

<!-- Resolved items go here. -->
