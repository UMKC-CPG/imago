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
  (5.6), layered-reproducibility build pipeline
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
  sequential loen for Fortran-side matchers (DESIGN
  5.10; later revised from a nested bootstrap to the
  makegroups sequence); matcher protocol architectural
  section
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
- [x] D14. Design structure acquisition (ARCH 9.5):
  cod_fish.py and cif2skl.py.  cod_fish.py is a four-verb
  COD front-end -- get (strict pinned fetch), search
  (composition/author query over COD result.php, exact
  composition via strictmin=strictmax=len(elements), a
  numbered candidate table saved to a session file), pin
  (resolve chosen rows' revisions by index -> manifest
  fragment), and rank (advisory triage by ambient
  conditions / named phase / spacegroup consensus / volume
  outliers).  cif2skl.py PRESERVES the space group: ASE
  parses the CIF (authored asymmetric unit + IT#/setting +
  full expansion), the spaceDB <IT#>_<letter> variant is
  chosen by VERIFICATION (write asymmetric unit + token,
  apply_space_group, match ASE's expansion -- no operation
  parsing, origin falls out), and the skl is written as
  asymmetric unit + token; no match is a hard error with a
  --space override.  No spglib (CIF symmetry taken as
  authored, checked not re-derived).  PSEUDOCODE for
  D11-D14 follows once each design lands.
- [ ] D19. Design the user-facing flight-construction surface
  (ARCH 9 / DESIGN 6).  Today a flight exists only as a
  hand-written Python `Flight` of `CalcUnit`s built on the C71
  `build_kpoint_convergence` helper; there is no way for a user to
  assemble a campaign without writing Python.  Decision 2026-06-17:
  design the front-end BEFORE coding it (design-first).  Settle the
  surface: a CLI, a config/spec file (e.g. a campaign TOML listing
  structures x option axes), or a thin convenience module; how the
  structure list is supplied (skl paths now; CIF/COD once D14 lands);
  how option-axis sweeps are expressed and mapped onto the existing
  SweepRecord; and where this sits relative to kaleidoscope's
  domain-agnostic core (Principle 9 -- the builder is a client of
  the dispatcher, never inside it).  Pairs later with a C-level
  implementation task once the surface is fixed.  Relates to D14
  (structure intake) and the C92 single-structure guidance path.
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
  StructureControl structure type per C64; the loen flow
  is now the makegroups sequence, not a makeinput
  subprocess).  Key judgment: cwd
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
  build pipeline (11.4), and validation harness
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
  species (DESIGN 5.6.6); the Fortran-side loen flow
  (DESIGN 5.10; since revised to the makegroups
  sequence -- no nested bootstrap or recursion guard).
  Update existing PSEUDOCODE 11.3 to reflect the
  Phase-2 selection flow rather than the Phase-1
  literal-label flow.  Done 2026-05-19: PSEUDOCODE
  11.3 expanded into seven sub-sections (11.3.a
  matcher protocol; 11.3.b preflight + coverage
  check; 11.3.c species pass with scope resolution;
  11.3.d entry pick with three-step precedence;
  11.3.e type pass with XANES split; 11.3.f the loen
  flow -- since revised to makegroups; 11.3.g driver)
  and 11.4 was
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
  the curation pattern for subsequent elements.  A curation
  run of the shipped producer, not new machinery.  CODE;
  DESIGN 5.7
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
the matcher protocol and the makegroups sequential
loen flow both invoke makeinput's Phase-1 emission
path.  The chain may be developed in parallel
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
- [x] C54. makeinput.py: introduce the Matcher
  protocol (ARCHITECTURE 8.9) and the `MATCHERS`
  registry.  Refactor the existing `group_reduce`
  algorithm into `ReduceMatcher` (with
  `needs_loen_run = false`); preserve current
  behavior of `-reduce` as the regression target.
  This lands before the new schemes so the existing
  test surface validates the refactor.  DONE
  (f772484): Matcher base + ReduceShellCode/
  ReduceShellLevel + ReduceStructureView +
  ReduceMatcher + BispecMatcher placeholder + MATCHERS;
  group_reduce reduced to the species-assignment
  workflow.  New test_makeinput_reduce.py adds two
  hand-traced regression fixtures plus matcher-surface
  unit tests; full suite 648 passed.
- [x] C55. makeinput.py: add `BispecMatcher` class
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
  DONE: `to_loen_input` rejects `by_element` (the
  guard C62 will replace), requires `twoj1`/`twoj2`,
  defaults the other three; `distance` is L2 (rejects
  unequal length); `representative` is the element-wise
  mean (rejects empty); also added `build_payload`/
  `extract_query_vector` (the `values` field) to
  complete the protocol surface.  14 unit tests in
  `test_makeinput_bispec.py`.  `compute_query` is NOT
  implemented for bispectrum -- the makegroups flow
  (C58) reads `fort.21` via `parse_loen_output` instead.
  KNOWN BUG: `parse_loen_output` was written against an
  incomplete `fort.21` spec (no header, no `site#`
  column); C89 revises it to the real enriched format.
- [ ] C56. makeinput.py: add `name=NAME` keyword to
  the `-target` and `-block` argparse handlers.
  Validate uniqueness across the run; store the name
  in the existing spatial-flag records so subsequent
  scope references can resolve it.
- [~] C57. makeinput.py: add `scope=NAME` /
  `scope=~NAME` to the `-reduce` handler.  Add a
  LOEN-block flag (e.g. `-loeninput`) that fills the
  `LOEN_INPUT_DATA` block from
  a `to_loen_input(sub_spec)` parameter dict instead of
  the hardcoded `4 4` -- this is a plain input-writer
  flag, NOT a grouping flag, used by makegroups' first
  makeinput call (DESIGN 5.10.2).  The `-bispec`
  *grouping* handler does NOT live in makeinput; it
  belongs to makegroups (C58).
  LOEN HALF DONE: `-loeninput CODE TWOJ1 TWOJ2 MAXNEIGH
  CUTOFF ANGLESQUEEZE` (six values in LOEN-block order)
  reconciles into six `loen_*` settings whose defaults
  reproduce the old hardcoded `1 / 4 4 / 50 9.0 0.85`
  and match the descriptor contract; the writer now
  emits them.  Flows through `_args_from_options` too,
  so the producer can drive it in-process.  REMAINING:
  the `scope=NAME` reduce handler (belongs with the C59
  reduce species pass, not the bispectrum path).
- [x] C58. makegroups.py (new; DESIGN 5.10, PSEUDOCODE
  11.3.f): the sequential bispectrum grouping helper,
  dual-mode (importable `group_by_bispectrum` + a
  `__main__` CLI; shebang + exec bit).  Steps: run
  makeinput on the skeleton with the LOEN params (C57)
  and no grouping; run `imago.py -loen -scf no`; read
  the enriched `fort.21` via
  `BispecMatcher.parse_loen_output` (C89); bucket atoms
  by fingerprint distance within the floor; rewrite the
  skeleton with explicit per-element species tags
  (`Si1,Si2,...,O1,O2,...` restarting at 1 per element,
  DESIGN 5.10.4).  Round-trip test on the per-element
  numbering.  No makeinput self-invocation, no recursion
  guard.  Replaces the retired nested-bootstrap design.
  DONE: pure core (P1 guard `_require_p1`; per-element
  bucketing `_assign_species_from_rows` via the shared
  `bucket_by_fingerprint`; in-place tag rewrite
  `_rewrite_skeleton_species`) + subprocess
  orchestration (`_run_makeinput`/`_run_loen`, fort.21
  discovery, `<skeleton>.orig` backup, scratch workdir).
  Row->atom mapping keys on (element, first-pass
  species) off the self-describing fort.21 -- no
  datSkl.map.  17 tests (round-trip via StructureControl,
  crystal/supercell rejection, CLI).  LIVE-VALIDATED
  2026-06-16 against the real engine (4-carbon P1 cell:
  c1 alone, symmetric c2/c3/c4 merged -> 2 species).  The
  live run fixed two seam bugs: makeinput takes no
  skeleton positional (reads ./imago.skl in cwd), and the
  loen descriptor is renamed by imago.py to
  `<edge>_loen<basis>.plot` (e.g. `gs_loen-fb.plot`), NOT
  `fort.21` -- `_find_loen_descriptor` globs `*loen*.plot`
  and raises if absent (absence = the loen run failed).
  Also confirms the C89 enriched fort.21 format live.
- [x] C59. makeinput.py: implement the Phase-2 species
  pass (DESIGN 5.6.4-5.6.7) for the schemes makeinput
  still owns -- position-based (`-target`, `-block`) and
  reduce.  Bucket reduce atoms by fingerprint distance
  (in-Python) with the matcher's similarity floor; pick
  one manifest entry per reduce species via fingerprint
  match + floor + default-tag fallback; propagate
  species choices through the type pass (XANES intact);
  emit per-type potentials in today's on-the-wire format.
  Bispectrum grouping and bispectrum potential-picking
  are NOT here -- makegroups (C58) / the orchestrator do
  them and hand makeinput explicit types (and `-pot`).
  DONE 2026-06-16 in two commits.  Part 1 (08d87ca): the
  matcher protocol gained `match_distance` (reduce
  element-only per DESIGN 5.2; bispectrum L2);
  `bucket_by_fingerprint` fixed to compare
  distance(representative, candidate) so the asymmetric
  reduce distance keeps the seed as reference; precedence
  2 (fingerprint match) added to `_select_augmented_pot_
  entry`; `_obtain_pot_info` picks per (element, species);
  `group_reduce` stores per-atom fingerprints + sub_spec
  and now buckets via the shared helper (retiring its own
  loop).  Part 2 (this commit): `scope=NAME`/`scope=~NAME`
  on `-reduce` + `name=NAME` on `-target`/`-block` with a
  `settings.named_regions` registry, so reduce regroups
  only inside/outside a position region (new species
  offset past preserved out-of-scope ones; relax_group
  compresses).  This also finishes C57's deferred
  `scope=NAME` reduce handler.  REMAINING (-> C61): a full
  end-to-end "small reference -> deterministic imago.dat"
  test; the unit-level pieces are covered.
- [~] C60. build_initial_potentials.py: bump the
  manifest reader to v2 (new `default` per entry,
  new `[[reference_solid.entry.fingerprint]]`
  declarations; new validation rules 7/8/9 per
  DESIGN 5.7).  Add the fingerprint-harvest step
  to the per-entry loop: invoke Python-side
  matchers in-process and Fortran-side matchers
  via cached `imago.py -loen -scf no` runs.
  Attach `FingerprintRecord`s and the `default`
  flag to each produced `PotentialEntry`.
  LIGHT HALF DONE (b88ba01 docs; a9be84f matcher
  build_payload; 9ba7f31 outputs["structure"];
  d761cec producer harvest): the manifest reader
  (already v2 with rules 7/8/9; rule 9 now wired to
  MATCHERS) plus the Python-side (reduce) harvest --
  it reads the run's expanded imago.fract-mi, maps
  atom_site->dat via datSkl.map, computes the shell
  code with ReduceMatcher, and stores an element-only
  `shell_code` (DESIGN 5.2, species dropped as
  non-transferable).  FORTRAN-SIDE (bispectrum) HARVEST
  DONE 2026-06-16: `build_loen_units(ref, struct, options)`
  now builds one structure-only `-loen -scf no` CalcUnit
  per distinct (method, sub_spec) -- kind="fingerprint"
  (convergence harvest skips it), job=loen/scf_basis=no,
  LOEN block via makegroups.loen_input_values, calc tag
  `loen-<method>-<sub_spec slug>` (slug-safe, %.6g floats).
  `harvest_loen_fingerprint` reconstructs that run dir from
  the same tag, finds the `*loen*.plot` descriptor
  (makegroups.find_loen_descriptor, live-validated),
  maps atom_site->dat via datSkl.map, guards the row's
  self-describing element, and wraps the vector via
  build_payload.  harvest_fingerprints now dispatches
  per-declaration (reduce in-process / bispectrum from the
  loen run) sharing one skeleton_to_dat.  This is the
  WITNESS path (per-entry atom_site = one fingerprint per
  declared type, no grouping) -- correct for crystalline
  references.  REMAINING (follow-on): a NON-crystalline
  reference would need the producer to call
  makegroups.group_by_bispectrum first (rewrite skeleton
  -> re-materialize) before harvest; not wired, since the
  current manifest references are crystalline.  C55 + C58
  + C89 all DONE; loen seam live-validated via makegroups.
  NOTE (2026-06-26): the DESIGN 5.7 harvest model is now
  *auto-discover one representative per distinct environment*,
  with manifest entries as optional customizations and the
  fingerprint recipe lifted to a `[characterization]` block
  (C101).  The per-entry witness path above is the interim
  crystalline-only implementation; it is superseded by C101
  plus the C88 bispectrum-keyed insert-or-skip.
- [x] C93. Revise the design chain to DECOUPLE the
  fingerprint-based initial-potential pick from the grouping
  scheme.  Decision 2026-06-17: precedence 2 (the environment
  fingerprint match, DESIGN 5.6.5) was wired to fire only when a
  `-reduce` scheme populated per-atom fingerprints, so a crystalline
  run using the normal crystallographic/space-group type grouping
  never matched and always fell through to the default isolated-atom
  entry.  That conflates two independent concerns -- how species are
  GROUPED vs which DB entry each species PICKS (the phase-i
  potential/type unification note).  Fix: the pick runs for every
  species, default-on, with a `-nofingerprint` opt-out that disables
  it (forcing the default-tagged entry; `-pot` still applies).
  TWO REGIMES settled across this session:
  (a) USER chose `-reduce`/`-bispec`: that family + the user's
  sub_spec drive both grouping and the match (grouping descriptors
  reused for the representative); the database NEVER overrules the
  user -- if it carries no fingerprint at that sub_spec the species
  misses to the default SILENTLY (best-effort, NOT a hard error;
  this SUPERSEDES the earlier "require-to-match / abort on mismatch"
  decision).
  (b) FILE-DICTATED species (crystalline / pre-assigned, no env
  flag): the database decides via the `preferred` record per family;
  the bispectrum-then-reduce priority applies ONLY here (prefer a
  preferred bispectrum record, else reduce).  Exactly one family,
  one sub_spec, one query -- a simple best-effort lookup, never a
  search across multiple sub_specs, never more than one loen run.
  The bispectrum (loen) query is fast, so the file-dictated case
  pays it up front; no loen run when the element has no preferred
  bispectrum record.
  PREFERRED-RECORD MECHANISM (the linchpin): the schema keeps
  allowing many sub_specs per family, but exactly one per family is
  flagged `preferred = true`, and the preferred sub_spec for a
  family is UNIFORM database-wide (set once in the manifest, stamped
  on every element file).  A non-preferred record at another
  sub_spec is always storable; only a divergent preferred flag is
  rejected.  This is what makes the file-dictated lookup
  unambiguous and one-loen-run.
  DESIGN DONE this session: 5.6.1 (`-nofingerprint`), 5.6.3 step 4
  (coverage note generalized, never fatal, no step 5), 5.6.4
  (grouping/pick decoupling + the two regimes + group->
  representative->match order + descriptor reuse), 5.6.5 step 2
  (two-regime single-sub_spec best-effort pick), 5.2 (FingerprintRecord
  `preferred` field + dataclass + per-element rule 10), 5.7 (manifest
  `preferred` field + canonical-sub_spec convention note + rules
  10/11).  PSEUDOCODE DONE: 11.3 intro + 11.3.0 (reduced flow now
  the -nofingerprint / no-preferred path), 11.3.a (consumer-side
  loen query via loen_descriptors), 11.3.b (noteCoverage
  generalized, fingerprinting_enabled), 11.3.c (env_species_ids
  pick-gate dropped), 11.3.d (two-regime preferred-driven pick +
  find_preferred), 11.3.g driver (new signatures, step renumber).
  ARCH DONE: 8.9 ("the pick is a per-species step" -- regime split,
  -nofingerprint, preferred/convention, consumer loen seam).
  DESIGN-CHAIN COMPLETE for C93; only the C94 code remains.  C61's
  negative tests follow the new default-on, never-abort behavior.
- [x] C94. Implement the C93 decoupling in code.
  makeinput.py: lift the reduce-only gate in
  `_species_query_fingerprint` (today it returns None unless
  `group_reduce` populated `atom_reduce_fingerprint`).  Implement
  the two-regime precedence-2 pick: when the user chose an env
  scheme, match that family at the user's sub_spec (reuse the
  grouping descriptors), miss -> default SILENTLY, never error on a
  DB that lacks the sub_spec; when species are file-dictated, read
  the element db's `preferred` record per family and use bispectrum
  if a preferred bispectrum record exists else reduce -- one family,
  one sub_spec, one query, accept misses.  Add the `-nofingerprint`
  CLI flag (disables precedence 2; `-pot` still wins).  Keep
  precedence 1 (`-pot` override) and the similarity-floor warn-and-
  fall-through unchanged.  initial_potential_db.py: parse, validate
  (per-element rule 10: exactly one preferred per PRESENT method --
  a family with records but no preferred is a load error), expose,
  and `save()` the `preferred` field on FingerprintRecord.
  build_initial_potentials.py (ties into C60): the harvest stamps
  `preferred` per manifest rules 10/11 (one preferred per
  (element, method); preferred sub_spec uniform database-wide).
  Depends on C93 (design); most meaningful once C91 seeds real
  fingerprints, but the wiring is independent of the seed.
- [x] C95. `--materialize-only` structure pre-flight in
  build_initial_potentials.py.  Expose Phase 1 (fetch + cif2skl
  convert) as a standalone CLI mode so a curator can get every
  reference structure materializing cleanly before filling in
  the run and harvest fields.  Adds `structure_cache_dir`, the
  relaxed reader `load_structure_sources` (only the rule 1/4/5
  checks materializing needs), `materialize_only` (per-solid
  report, continues past failures), and the `--materialize-only`
  / `--materialize-dir` flags.  CODE; PSEUDOCODE 11.4.
- [x] C96. Wire the `kpoint_integration` Gaussian smearing end
  to end.  The `gaussian-<sigma>` token's width was parsed but
  discarded (only the integration code reached makeinput); now
  the producer's `_thermsmear_for` extracts the sigma (eV) and
  `make_producer_options` forwards it as the `thermsmear`
  option, makeinput gains a `-thermsmear` dest that overrides
  the rc `therm_smear_main` and writes `THERMAL_SMEARING_SIGMA`.
  A bare `gaussian` keeps the rc default (0.0).  CODE; DESIGN
  5.7; PSEUDOCODE 11.4.  Also captured in VISION Goal 3: use the
  database to benchmark non-Wigner functionals (and a possible
  future libxc/xclib integration).
- [x] C97. Manifest authoring pipeline ("Approach C").  Extract
  the manifest schema (the dataclasses, `load_manifest_v2`, and
  `load_structure_sources`) out of the producer into a leaf
  library `curation_manifest.py` (which imports only
  initial_potential_db / guidance_db, nothing in the tool
  layer); add the manifest writer (`format_manifest` /
  `write_manifest`, human-readable TOML that round-trips through
  `load_manifest_v2`); and add `expand_manifest.py`, the
  sketch-to-manifest authoring tool (interactive `-i` and
  mechanical modes, with an injected prompt for testability).
  cod_fish stays a pure discovery tool.  CODE; DESIGN 5.7; ARCH
  8.5/8.7; PSEUDOCODE 11.4/11.6.
- [x] C98. `cod_fish pin` emits a ready-to-use sketch. It now
  prints a complete sketch -- a `schema_version = 2` header plus
  one `[[reference_solid]]` stub per pinned structure -- and
  auto-derives each `reference_id` from the CIF metadata it already
  downloads: `<formula>_<H-M symbol>_<IT number>_<year>` (e.g.
  `si_fd-3m_227_2010`). The year dates the entry and separates
  phases sharing a space group (4H vs hcp Si, both P6_3/mmc); a
  residual clash gets a trailing counter (rule 5). Removes the
  hand-prepend of schema_version and the per-entry rename from the
  seeding workflow, so `cod_fish.py pin <ids> > sketch.toml` is a
  ready sketch. Each stub also carries two CIF-read discovery hints
  -- `elements` (composition) and `source_description` (chemical
  name + space group + year) -- which `expand_manifest` auto-fills
  into each entry's element and description, so the curator invents
  neither (they are not schema fields; the producer ignores them and
  the finished manifest omits them). CODE; DESIGN 5.7; ARCH 9.5.
- [x] C99. cif2skl tolerates near-full site occupancy. Experimental
  CIFs report a fully-occupied site as a refined value just below 1
  (e.g. BC8 silicon, COD 4350826, lists `Si 0.9999`); the old check
  refused anything below `1 - 1e-6`, spuriously rejecting it. Now a
  new `_FULL_OCCUPANCY_TOLERANCE` (1e-2) treats occupancy within 1%
  of full as full, while genuine partial occupancy (disorder /
  vacancies, well below the band) is still refused. CODE; ARCH 9.5.
- [x] C91. **Side-quest: populate the augmented potential
  database with real fingerprint records.**  Today every
  `share/atomicPDB/<elem>/s_gaussian_pot.toml` is
  Phase-1-shaped -- a single `default = true` entry and zero
  `[[potential.fingerprint]]` records (an audit of all 103
  element files on 2026-06-16 found no fingerprints anywhere
  in the installed database).  The Phase-2 consumer path
  (the `-reduce` / `-bispec` per-species entry pick,
  precedence 2 in DESIGN 5.6.5) therefore has nothing to
  match against: grouping still runs, but every species
  falls through to the default entry, so the fingerprint
  pick cannot be exercised against live data at all.  This
  task runs the producer chain end-to-end to seed real
  fingerprints: curate a small reference set (building on
  C49's first solid), run `build_initial_potentials.py` as a
  kaleidoscope client (C48.3 / C74) so each reference solid
  converges its SCF and the C60 harvest attaches both the
  reduce (Python-side, in-process) and bispectrum (loen)
  `FingerprintRecord`s, and land the updated per-element
  TOML files back under `share/atomicPDB/`.  Deliverable: at
  least the benchmark elements carry one or more
  fingerprint-bearing entries so the Phase-2 selection (C61)
  and the later dedup harvest (C88) have genuine data to
  work with.  Dependencies: C49 (curation pattern), C74 (the
  one remaining live producer-as-client validation seam),
  and cluster time for the SCF + loen runs.  Note: C61's
  synthetic-DB fixture for its positive case still stands
  (deterministic, no cluster needed), but seeding real data
  lets C61's positive path also be confirmed against the
  wild database and surfaces any producer<->consumer schema
  drift that hand-built fixtures would miss.

  **DONE 2026-07-20**, by the clean production run (orchestrator
  job 15032332, build `61429b9`, smearing off).  The eight-solid
  Si manifest ran end to end and
  `share/atomicPDB/si/s_gaussian_pot.toml` now carries seven
  reference-derived entries beside `isolated`, each bearing a
  `reduce` and a `bispectrum` record -- fourteen fingerprints in
  all -- so the Phase-2 consumer path finally has live data to
  match against.  Seven and not eight because
  `si_fd-3m_227_1962_3` fell within the similarity floor of an
  already-stored diamond entry and was dropped by
  `insert_or_skip`: the DESIGN 5.2.3 dedup behaving exactly as
  specified.  That the other five diamond allotropes were each
  judged novel is a floor-calibration observation for C88/C103,
  not a defect -- small lattice-constant differences currently
  read as distinct environments.
- [ ] C61. Add end-to-end Phase-2 tests.  For the
  in-makeinput path: a small reference runs through
  reduce bucketing + entry-pick + emit producing a
  deterministic Imago input file, with negative tests
  for the missing-fingerprint-family case (the run does
  NOT abort -- it emits the coverage info note and every
  species falls through to the default entry, DESIGN
  5.6.3 step 4 / PSEUDOCODE 11.3.b `noteCoverage`) and
  the similarity floor (sub-threshold match falls back to
  default with a warning).  For the makegroups path
  (C58): the sequential loen -> bucket -> skeleton-rewrite
  chain produces the expected per-element species tags.
- [~] C89. fort.21 enrichment + parser fix (DESIGN
  5.10.3).  CODE DONE, pending a live recompile + loen
  run to validate the Fortran format.  Fortran (loen.f90,
  the `open(unit=21,...)` block): each row now leads with
  `site#`, `element`, `species`, `type_in_species`,
  `type_flat` (looked up from
  `potTypes(potSites(i)%potTypeAssn)`), and the header
  names them.  Python (`BispecMatcher.parse_loen_output`):
  skips the header, reads the 5 identity columns + the
  `twoj2+1` components, returns a `LoenSite` record per
  site.  Also fixed the second bug: the component count is
  `twoj2+1` (coupling channels j in `|j1-j2| <= j <=
  j1+j2`), NOT `2*twoj2+1` -- corrected in code, DESIGN
  (5.2/5.3/5.10), PSEUDOCODE, ARCHITECTURE 8.9, and the
  tests.  `test_makeinput_bispec.py` updated to the real
  format.  REMAINING: the live recompile/run check.
  Prerequisite for C58 and the C60 loen harvest.
- [x] C89.1. Extract the matcher protocol from
  makeinput.py into a neutral `src/scripts/matchers.py`
  library (Matcher / Reduce* / Bispec* / LoenSite /
  MATCHERS).  makeinput imports ReduceMatcher +
  ReduceStructureView inside `group_reduce`;
  build_initial_potentials and the bispec/reduce tests
  import from `matchers`.  Done so `makegroups.py` (C58)
  can import `BispecMatcher` without a makeinput<->makegroups
  import cycle.  ARCHITECTURE 8.9 "Location" + layout
  updated.  669 tests pass.
- [ ] C90. Full grouping extraction (the makeinput /
  makegroups split proper).  Move the geometric species
  pass -- `group_reduce`, the `-target` / `-block`
  handlers, and `scope=` resolution -- out of makeinput
  into makegroups, so makeinput becomes a pure
  input-writer (reading explicit skeleton types + the
  xanes electronic pass).  Deferred deliberately: reduce
  grouping composes with the position-based flags via the
  `scope=` feature (a `-reduce scope=region` confined to a
  `-target`/`-block` region, in CLI order, DESIGN 5.6.4),
  so this is a real refactor with the reduce regression
  tests (C54 fixtures) to carry, not a mechanical move.
- [ ] C88. Dedup storage model + native/witness
  fingerprints (DESIGN 5.2.2-5.2.4).  The database
  stores *distinct environments, not atoms*, so a
  tens-of-thousands-atom model adds entries only for
  its genuinely distinct environments.  Work:
  (a) add the `type_assignment` provenance field
  (already in 5.2's table) and derive each
  fingerprint's native/witness role from it
  (`M == type_assignment`); the producer must record
  the run's grouping scheme.  Reduce is always a
  witness, never the assigning method (5.2.2), so it
  stays a droppable column.  (b) Symmetric dual
  harvest: compute *both* registered methods (native +
  witness) for every harvested environment -- reduce
  is free given the structure, a bispectrum loen pass
  is cheap next to the SCF.  (c) Make the build
  INCREMENTAL: load each existing per-element file (or
  seed an empty one), refresh only the `"isolated"`
  baseline, and append this run's harvest -- never
  reset -- so several manifests accrete into one
  database (DESIGN 5.7).  (d) Make the harvest an
  *insert-or-skip* keyed on the **bispectrum**
  descriptor at the preferred sub_spec (the
  transferable one every entry carries -- NOT a
  conservative all-methods union and NOT symmetry,
  which cannot compare across structures; 5.2.3): a
  duplicate is SKIPPED (the first representative's
  potential stands), a novel environment appends.  The
  dedup tolerance is the producer-side mirror of C61's
  similarity floor.  The leaner model stores NOTHING
  extra per entry -- no counts, no spread, no
  contributor list; `coefficients` is the
  representative's potential verbatim.  Optional: carry
  the origin label onto records so the per-atom cross-
  method pairing survives normalization.  Lands after
  the Phase-2 base chain (C53-C61); the
  bundled->normalized (potential pool + per-method
  index) migration and the learned-predictor training
  pipeline are separate, later items this sets up.  The
  statistical merge (mean/spread/counts) is the
  deferred C103 upgrade.
- [ ] C101. Manifest characterization + customization model +
  auto-harvest (DESIGN 5.2.1/5.2.2/5.7).  Recast the
  manifest and harvest from "one entry per declared atom
  site" to "auto-discover one representative per distinct
  environment; entries are optional customizations."  Work:
  (a) curation_manifest.py: add the top-level
  `[characterization]` block (the database-wide preferred
  fingerprint recipe, one sub_spec per method; makes rule
  11 structural), add a persisted per-solid
  `source_description` field, and relax `ReferenceEntry`
  so `atom_site`, `element`, `default`, and `description`
  are all optional customizations (label already optional via
  C87).  Update validation rules 2/3/6/7/8/10/11 to the
  revised DESIGN 5.7.  (b) build_initial_potentials.py:
  harvest one order-independent representative per
  distinct environment (the assigning method's partition;
  5.6.5), compute the `[characterization]` fingerprints
  for every environment (preferred), apply per-entry
  overrides (rare, non-preferred), and auto-compose each
  environment's description from `source_description` +
  species/site when no customization supplies one.
  (c) expand_manifest.py: emit the `[characterization]`
  block and per-solid settings; entry customizations optional.
  Pairs with C88 (the insert-or-skip the auto-harvest
  feeds) and reframes the C60 per-entry witness path.
  CODE; DESIGN 5.2.1/5.2.2/5.7; PSEUDOCODE 11.4; ARCH 8.5.
- [ ] C103. Statistical merge + exact rebuild (deferred
  upgrade of C88's skip-on-match; DESIGN 5.2.3
  "Deferred").  Replace insert-or-skip with a merge that
  folds EVERY atom mapping to an environment into the
  stored entry, adding: (a) per-coefficient `coefficients`
  as a running **mean** plus a per-coefficient
  `coefficient_std` spread (the empirical test of the
  fingerprint's fidelity, and a per-entry confidence /
  predictor variance, 5.2.4); (b) an atom `multiplicity`
  (sample weight) and a `model_count` (distinct-solid
  corroboration); (c) one **contribution record per
  reference solid** (`reference_id`, `atom_count`,
  `coeff_sum`, `coeff_sumsq`) as the order-free source of
  truth the summaries derive from -- which also unlocks
  **removal** of a solid (drop its record, re-derive) and
  an **exact wholesale rebuild** to a byte-identical
  file, the two properties skip-on-match gives up.
  Thread the new fields through the field list +
  validation (5.2), the `PotentialEntry` dataclass (5.4),
  the reader/emitter (5.5 / PSEUDOCODE 11.1-11.2), and
  the harvest (5.7); the merge asserts equal alpha SETS
  before averaging.  CODE; DESIGN 5.2.3; PSEUDOCODE 11.4.
- [x] C104. Manifest `[defaults]` hoist + cod_fish complete
  manifest (DESIGN 5.7; ARCH 8.5).  Design DONE (this commit).
  Code work: (a) curation_manifest.py: add a top-level
  `[defaults]` block holding the five shared run settings
  (`basis`, `functional`, `kpoint_integration`, `kpoint_spec`,
  `scf_threshold`); make those per-solid fields optional; the
  reader parses `[defaults]` and RESOLVES each solid's omitted
  setting from it (so the producer reads fully-populated solids
  unchanged), validating that every setting is resolvable
  (rule 2); the writer emits a `[defaults]` block and per-solid
  settings only when they differ from the default.  Move the
  default-value constants, `default_run_settings()`,
  `default_characterization()`, and the customization-comment
  template into curation_manifest.py as the shared source.
  (b) cod_fish.py: emit a complete, runnable manifest by default
  (the recipe + `[defaults]` + stubs, values from the shared
  library), with `--sketch-only` for the bare stubs.
  (c) expand_manifest.py: read the sketch from stdin when no
  file is named; emit the `[defaults]` block instead of
  repeating settings per solid; source defaults from the shared
  library.  `system_type` stays per-solid (structure metadata,
  not defaulted).  CODE; DESIGN 5.7; ARCH 8.5.
  DONE: all three sub-items shipped -- `curation_manifest`
  parses, resolves, and writes `[defaults]`; `cod_fish` emits a
  complete runnable manifest with `--sketch-only` for bare stubs;
  `expand_manifest` reads a sketch from stdin and hoists the
  shared settings via `shared_defaults()` + `sparse_solid()`.
  The box stayed unchecked past the work because the matching
  PSEUDOCODE 11.6 sync was still outstanding; that sync is C114,
  which is now done as well.
- [x] C105. On-disk potential-file schema migration / version
  guard (initial_potential_db).  Surfaced 2026-07-01 during the
  C91 Si seed: the incremental producer loads each existing
  `share/atomicPDB/<elem>/s_gaussian_pot.toml` through the strict
  reader, so a file written before a now-required field was added
  is rejected with a bare `missing required field` ValueError and
  no recovery but hand-deletion.  The concrete trigger was a
  pre-B3 `si/s_gaussian_pot.toml` whose Imago provenance predated
  the required `type_assignment` field (added by C88(a)); deleting
  the stale file let the producer re-seed and re-harvest cleanly.
  This is acceptable now while the databases are disposable seed
  data, but once a curated set accretes across many manifests that
  we do NOT want to rebuild, an older file must be handled
  gracefully.  Work: have `initial_potential_db.load()` compare
  the file's `schema_version` against the current one and EITHER
  migrate a known-older file forward (fill newly required fields
  with a documented default, e.g. `type_assignment` from the
  provenance already present) OR raise a clear, actionable error
  ("this file predates schema vN; regenerate it with
  build_initial_potentials.py, or migrate with <tool>") instead of
  the low-level missing-field message.  Pairs with the guidance /
  resource dataspaces, which have their own `*_migrate.py` tools
  (ARCH 10 / 11) -- the initial-potential DB should grow the same
  story.  CODE; DESIGN 5.2/5.5.

  **DONE 2026-07-21.**  Reading the code first reframed the task:
  the version guard already existed (rule 1) and did NOT fire,
  because `type_assignment` was added to the required Imago
  provenance set without moving `schema_version` off 2.  The file
  was a legal-looking v2 that failed a v2 field check, so no guard
  could have caught it.  The fix therefore starts with a policy,
  not a migrator: **the required-field set IS the version**, and
  changing it is a bump.  Written top-down as new DESIGN 5.2.5
  (versioning rule, the reader's four outcomes, the honesty test
  for a derivation, the error contract, where the machinery
  lives), rippled into DESIGN 5.2's key table and 5.4's rule 1,
  ARCHITECTURE 8.7 (which now names the future
  `potential_migrate.py` beside `guidance_migrate` /
  `resource_migrate`), then PSEUDOCODE 11.1, then the code.

  Policy chosen: **migrate when honest, else refuse.**  Each bump
  declares, per newly required field, either a derivation from
  what an older file already carries or `NOT_DERIVABLE`; a
  plausible default is not a derivation, because an invented value
  leaves the file well-formed and wrong where no later check can
  catch it.  `type_assignment` is the worked example of
  not-derivable -- nothing else in an older file records which
  scheme drew the type partition, so filling it would silently
  mislabel every fingerprint's native/witness role.

  Code (`initial_potential_db.py`): `CURRENT_SCHEMA_VERSION`, the
  `NOT_DERIVABLE` sentinel, a `SchemaMigration` dataclass and the
  `SCHEMA_MIGRATIONS` table (empty -- version 2 is current), and
  the `apply_schema_migrations` gate.  `load` checks
  `schema_version`'s presence on its own and runs the gate BEFORE
  the required-field sweep, which is the ordering the whole fix
  turns on: an out-of-date file is diagnosed as out of date rather
  than as missing a field.  Every refusal names the file, both
  versions, and a recovery; a newer-than-us file is told to update
  Imago and explicitly warned NOT to regenerate.  Seven tests
  cover the four outcomes plus the two ordering guarantees; 997
  non-integration tests pass, and all 103 installed element
  databases load through the gate unchanged.

  Deferred, as ARCHITECTURE 8.7 now records: the bulk
  `potential_migrate.py`.  It is a convenience, not a correctness
  requirement -- in-memory migration plus the producer's next save
  already carries forward every file a run touches -- so it earns
  its keep only for files no producer run will revisit.
- [x] C106. Single-source the producer's fingerprint-declaration
  set so the build and harvest sides cannot drift (DESIGN
  5.7/5.10).  Back-burner.  Scenario that surfaced it (2026-07-02,
  C91 Si seed live run): the producer dispatched every SCF
  convergence run for a solid, then died at harvest with
  `FileNotFoundError: no loen descriptor ... the loen run did not
  complete`.  Root cause was a split source of truth --
  `harvest_fingerprints` reads the declaration set as
  `characterization + entry.overrides` (the database-wide recipe
  plus per-entry overrides), but `build_loen_units` only iterated
  `entry.fingerprints`.  After the C93/C94 decoupling moved the
  bispectrum recipe into the database-wide `[characterization]`
  block (with no per-entry overrides in a Si default manifest),
  the build side built no loen unit at all, so the harvest read a
  descriptor that was never dispatched.  The immediate bug is
  fixed by passing `characterization` into `build_loen_units` and
  unioning both sources (calc-tag deduped), but nothing ENFORCES
  that the two sides stay in agreement: a future third path (a new
  fingerprint family, an override-precedence rule) could
  reintroduce the same drift.  Work: compute the
  `(method, sub_spec)` declaration set ONCE (e.g. a
  `producer_fingerprint_declarations(ref, characterization)`
  helper) and have BOTH the loen-unit build and the harvest
  consume that single list, so the build set is the harvest set by
  construction.  CODE; DESIGN 5.7/5.10.

  **DONE 2026-07-21**, with C107 (they share a design rung and a
  commit).  New DESIGN 5.10.6 states the rule and why the drift is
  structural rather than incidental; PSEUDOCODE 11.4 gains
  `fingerprintDeclarations` and `producerFingerprintDeclarations`;
  the code follows.  The shapes differ between the two sides, which
  is what made the single source non-obvious: the harvest composes
  the set for ONE environment, while the build cannot know the
  environments yet (they are discovered from the converged run).
  So the rule is defined per-environment and the build applies it
  to every case that could arise -- the override-less environment,
  plus each manifest entry -- and unions.  The build set is then a
  superset of any harvest set by construction.  Two consequences
  documented: a site-less customization may build one spare
  geometry-only run (the right direction to err), and the calc-tag
  dedup still collapses the repeated recipe to one unit per
  `(method, sub_spec)`.

  A prerequisite surfaced on the way in and was done first, as its
  own commit (aa6f341): PSEUDOCODE 11.4's `harvestFingerprints`
  still carried the interim C55/C58 guard refusing every loen-side
  declaration, and `harvestLoenFingerprint` was spec'd as an unbuilt
  path with naive positional indexing, though the code had
  implemented the finished state since C60 and DESIGN 5.10 described
  it correctly.  Writing the single-source rule into a spec that
  still refused loen declarations would have been incoherent.  That
  repair was a legitimate upward edit (code verified faithful to
  DESIGN first, per the chain rules), kept separate so a later
  reader can tell catching-up from new specification.

- [x] C107. Fail fast before dispatch when a fingerprint the
  harvest will read has no dispatched unit (DESIGN 5.10/6.2).
  Back-burner.  Same scenario as C106: the missing loen unit was
  not detected until the HARVEST phase, after the producer had
  already spent minutes of cluster SCF time on all eight solids'
  convergence sweeps -- the failure was silent-until-harvest.  A
  cheap pre-dispatch invariant would catch it before any run
  launches: after assembling `all_units`, assert that every
  Fortran-side declaration the harvest will read (per solid, over
  `characterization + overrides`) has a matching loen unit in the
  flight, and raise a clear error naming the solid and sub_spec if
  not.  This is a defensive backstop that makes the C106 class of
  bug loud and free instead of expensive and late; it is worth
  keeping even after C106 single-sources the declaration set,
  since it also guards against a unit that was built but dropped
  during dispatch assembly.  CODE; DESIGN 5.10/6.2.

  **DONE 2026-07-21**, with C106.  `assert_loen_coverage` runs in
  the producer main once the units are assembled and before any is
  sent; it is a set comparison over calc tags, so it costs nothing
  next to the minutes of cluster SCF time it protects.  The error
  names the solid AND the sub_spec, since the curator needs to know
  which declaration went unrun.  Kept even though C106 should make
  it unreachable, exactly as this entry argued -- it guards what
  the rule does not, and a cheap invariant that can only fire when
  something else is already broken is the kind worth asserting.
  Four tests, including one that a convergence-tagged unit sharing
  the calc tag does NOT satisfy the check.

- [ ] C108. Intermediate-scratch cleanup for the producer (and a
  reusable cleanup subsystem).  Motivation: the Si seed run
  (2026-07-02) left 3.7 GB of per-calc scratch under each run
  directory's `intermediate ->` symlink (roughly 20 MB per calc
  dir), against only 18 MB of kept home-side artifacts
  (status/result/scfV/descriptor).  Scratch of this kind fills up
  fast and is tedious to locate and remove by hand -- and stale
  units from earlier manifests linger in the shared workspace (the
  seed run's workspace still held `si_diamond`, `si_fd-3m_227_2010`,
  and a half-finished `si_p63mmc_194_2018` from prior experiments).
  Three layers, built in this order:
    (c) A STANDALONE cleanup script -- the eventual home of the
        logic: a generic, selective find/remove over intermediate
        scratch, with options that let the user target what to
        prune (by workspace, by unit, by calc kind, by age, dry-run
        preview, keep-the-harvestable-artifacts).  This is the real
        deliverable and can be written later.
    (a) Once (c) exists, the build/harvest script gains a
        `--clean-after` option that simply invokes the standalone
        script with the known-good options a just-finished harvest
        can supply (which units harvested, what to keep) -- so the
        post-harvest cleanup and the standalone tool are one code
        path, not two.
    (b) NEAR-TERM, before (c) lands: no cleanup option implies NO
        cleaning (today's behaviour, unchanged).  A `--tidy-run`
        option turns on prune-as-you-go, discarding a unit's
        superseded intermediate scratch as the flight advances.  The
        pruning ACTION must be written generically -- a
        builder-supplied policy/hook -- so builders other than the
        k-point convergence producer can define what "safe to prune
        now" means for their own units.
  Level note: this introduces a small cleanup subsystem, so the
  standalone-script boundary and the generic-pruning-hook belong in
  ARCHITECTURE/DESIGN before (c) is coded.  CODE + DESIGN.

- [x] C109. Decide the default cell -- full (conventional) vs
  primitive -- for the materialized `imago.skl`.  Surfaced by the
  Si seed run (2026-07-02): a primitive cell has fewer atoms, so a
  smaller secular equation and a faster SCF, and the harvested
  quantity is a per-species / per-environment POTENTIAL, which is
  cell-invariant -- so the choice changes run COST, not the
  science.  Two wrinkles to work through before picking a default:
  (1) COD CIFs arrive as the conventional or as-published cell, so
  going primitive means a symmetry reduction (spglib / ASE
  `find_primitive`) inside `cif2skl`, and the reference_id <->
  structure mapping and `datSkl.map` type assignment must stay
  consistent through that reduction; (2) the `kpt-density-N` sweep
  metric must be defined so a primitive and a conventional cell
  receive EQUIVALENT k-meshes -- automatic if the density is a
  reciprocal-space length target, but not if it is a per-axis
  count.  Decide the default and whether it becomes a manifest
  knob.  DESIGN 5.7 (materialize_structure) / ARCH 9.5 (cif2skl).

  **IN PROGRESS 2026-07-21.  Both wrinkles above turned out to be
  already resolved -- neither the way this entry assumed.**

  (1) The reduction needs no spglib and no `cif2skl` work: the
  skeleton carries a `full` / `prim` token on a line of its own,
  `structure_control` reads it (`do_full_cell`), and the reduction
  is done with the space-group operations already in hand.  It runs
  one way only (`full -> prim`), which is the direction needed,
  since `cif2skl` writes the conventional cell.  The whole climb
  path was already prim-ready too: `axis_classes_for_cell` has an
  explicit `"prim"` branch conjugating the conventional-abc point
  ops into the loaded basis (DESIGN 2.7), and the producer already
  computes `cell_mode` from the loaded skeleton and recomputes the
  reciprocal lattice from the final cell.  Only ONE line forced the
  conventional cell: a hardcoded `"full\n"` in `cif2skl`.

  (2) C118 settled the k-mesh question by construction.  The
  density ladder that prompted the worry is gone, and the climb
  picks meshes from a reciprocal-space SPACING (`counts_at_spacing`:
  `|b_i| / h`) derived from a volume density.  A primitive cell has
  an n-fold larger reciprocal cell and so receives n-fold more
  k-points at the same density -- the physically equivalent
  sampling, with nothing to decide.

  **The knob is built** (this commit): `cell` joins `[defaults]` as
  a sixth run setting, `"full"` (default) or `"prim"`, per-solid
  overridable.  It is a COST setting, not a physics one -- the
  harvested potential and every fingerprint are cell-invariant --
  so DESIGN 5.7 records that it selects no predictor sub-model,
  that a recorded `converged_mesh` is NOT comparable across cells,
  and that a curator's own `structure_path` skeleton keeps its own
  token.  Alone among the run settings it is exempt from rule 2's
  resolvability requirement, with the exemption's expiry stated:
  the moment `cell` is recorded on an entry it becomes emitted
  knowledge (VISION Principle 11) and rejoins the rule.

  **A cache trap was found and fixed with it.**  The materialized
  skeleton was cached as `<reference_id>.skl` and reused whenever
  present, so changing the cell and re-running would silently hand
  back the earlier cell's file -- no error, a well-formed skeleton,
  the wrong answer reported as success.  A full-vs-prim comparison
  would then have compared full against full and read as a
  confirmation.  Cached skeletons are now `<reference_id>-<cell>.skl`,
  and DESIGN 5.7 states the general rule: the cached name carries
  every manifest setting that changes what the conversion writes.
  The relaxed `--materialize-only` reader resolves `cell` for real
  rather than leaving it at a placeholder, so a pre-flight and the
  run that follows it still share one cache.

  **First evidence, no cluster time needed.**  New `cif2skl` tests
  drive the reduction through `structure_control` end to end: fcc
  gold goes 4 -> 1 atom and diamond silicon 8 -> 2, with the stored
  asymmetric unit identical either way.  So the reduction is
  correct at least for F-centred cubic.  A regression test also
  pins the ordering the conversion depends on: the space-group
  candidates must be built and verified as `full` (the CIF's atom
  list is a conventional list), and only the winner is rebuilt with
  the caller's cell.

  **NEXT: the live comparison.**  Same manifest twice, `cell =
  "full"` and `cell = "prim"`, into two `--pdb-root`s (which
  separates the structure cache, the workspace, and the database in
  one flag).  The eight seed solids already span three centrings --
  F (six diamond allotropes), I (si_ia-3), C (si_cmce) -- so they
  exercise the reduction broadly with nothing new to curate.
  Compare fingerprints and coefficients for cell-invariance, and
  converged mesh and wall time for the payoff, which scales as
  `1/n^2` (both dominant cost terms carry atom count cubed against
  one factor of k-point count): about 16x for the diamond cells and
  4x for the other two.  THEN decide the default.

  **DECIDED 2026-07-21: the default is `prim`.**  Both live runs
  completed 8/8 (jobs 15165434 full, 15167972 prim) from a wiped
  workspace, so nothing was served from the run-reuse cache, and
  they were run SEQUENTIALLY rather than side by side -- the
  earlier attempt had both on one node, which made its timings
  worthless.

  Cost, per converged calculation:

      diamond Si   (8 -> 2 atoms)   5.1 s -> 2.7 s    1.9x
      Si III BC8   (16 -> 8 atoms)  43.6 s -> 22.8 s  1.9x
      si_cmce      (metal)          13.8 s -> 12.1 s  1.1x

  So about **twice**, not the 16x the `1/n^2` scaling argument
  predicted -- at these sizes the cubic term simply does not
  dominate, and the 8x more k-points a primitive F-centred cell
  needs ([6,6,6] -> [12,12,12]) nearly cancels the cheaper
  diagonalization.  The whole-campaign figure is smaller again,
  1.33x (492 s -> 369 s), because the primitive climb walks
  eleven rungs to the conventional seven.  The campaign number is
  the pessimistic one and NOT the one that governs: a climb is a
  one-time seeding cost per material, while a production run reads
  its density from the guidance dataspace and pays only the single
  converged calculation.  At thousands of simulations, 1.9x is the
  figure that compounds.

  Correctness cost: none measurable.  Converged energies agree to
  **0.002 meV/atom** on every insulator (0.8 meV/atom on si_cmce,
  the metal, which settles at a deliberately rough mesh anyway).
  Reduce and bispectrum fingerprints match across cells in both
  directions for all three centrings, once C126 made the reduce
  descriptor transferable.

  The one genuine risk -- the prim-only conjugation in
  `axis_classes_for_cell` (DESIGN 2.7), the same class of math as
  the `buildAtomPerm` hex/trig bug -- is now directly validated
  rather than merely unexercised.  si_cmce's reduction GENUINELY
  changes its axis classes (`1 2 3` conventional -> `1 1 3`
  primitive, because a C-centred lattice's primitive vectors
  `(a +/- b)/2` have equal length), and the Python port and imago's
  own runtime `computeAxisClasses` agree on the partition for every
  cell and centring tested.  Two independent implementations
  agreeing on a non-trivial answer is much stronger evidence than
  an absence of errors.

  Changed: `DEFAULT_CELL` -> `"prim"` with the measurement recorded
  beside it, `default_run_settings()` emits it, DESIGN 5.7 rewritten
  around the measured figures (replacing the `1/n^2` estimate),
  PSEUDOCODE constant updated.  A manifest naming no cell now gets
  the primitive reduction; `cell = "full"` per solid or in
  `[defaults]` restores the old behaviour for any structure that
  needs it.

  **Follow-up:** recording `cell` in each entry's provenance, so a
  run can be reconstructed from its entry and the cost dataspace
  can tell an 8-atom run from a 2-atom one.  Tracked as C127.  Note
  it does NOT affect dedup -- an earlier draft of this entry claimed
  it did, and that was wrong; see C127 for the measurement that
  corrected it.

- [ ] C127. Record `cell` in each entry's provenance, and make it a
  required, resolvable run setting.  Raised when C109 moved the
  default to `prim` (2026-07-21), so a database can now hold entries
  harvested under either cell.  LOW priority -- see the correction
  below, which removed the urgency the first draft claimed.  DESIGN
  + CODE; DESIGN 5.2 (provenance block) / 5.2.5 (the version gate) /
  5.7 (the run setting) / 8.2 (the cost size signature).

  **A correction, recorded because the reasoning is the useful
  part.**  This entry was first written claiming that mixing cells
  would leave DUPLICATE environments the dedup could not collapse,
  because two harvests of one environment in different cells differ
  in the last digits of their stored distances.  That was wrong, and
  measuring it took one script.  The dedup keys on the preferred
  BISPECTRUM descriptor (5.2.3), not on the reduce shell code: the
  engine computes it from a periodic neighbour list, emits seven
  significant figures, and both cells produce the SAME seven --
  bitwise identical across all nine channels, L2 distance 0.0
  against a similarity floor of 0.10.  Inserting every prim entry
  into a full-harvested database under deliberately different labels
  (so the label-replace path could not mask the question) skipped
  all seven as duplicates: 8 entries before, 8 after.  The last-digit
  differences that prompted the worry are in the REDUCE distances,
  computed in Python from coordinates at double precision -- a
  different fingerprint family, and one the dedup never consults.

  **`cell` must never enter the dedup or the match.**  This is the
  trap the correction exposes, and it is worth stating as a rule
  rather than leaving implicit.  Folding `cell` into either key
  would manufacture distinctions the physics does not have: one
  environment stored twice because a curator drew a different cell,
  inflating the database with redundancy and teaching the learned
  predictor (5.2.4) that cell choice is a property of an atom's
  surroundings.  C126 and the C109 validation exist precisely to
  establish that it is not.  `cell` belongs beside `commit` and
  `generated_at` -- provenance describing where a number came from,
  never a key describing what it means.

  **What recording it actually buys.**  Two things, both real but
  neither urgent:
    (a) RECONSTRUCTION.  The Imago provenance fields exist so an
        entry's originating SCF can be identified and re-run.  The
        cell is now part of what defines that run -- it fixes the
        atom count and, through the reciprocal cell, the mesh the
        climb converges on -- so an entry without it cannot be
        reproduced from what it records.
    (b) THE COST DATASPACE (section 8), whose size signature is
        built on `atom_count` and `secular_dimension`.  A primitive
        cell halves both.  A resource model fitted across mixed
        cells with no field distinguishing them would be fitting
        the cell choice as unexplained scatter.  This is the
        stronger of the two, and it only bites once the cost
        dataspace is actually being trained (C77-C82).

  **The work.**  Add `cell` to the required Imago provenance set
  beside `type_assignment`, which by the C105 rule (the
  required-field set IS the version) makes it a v2 -> v3 schema
  bump.  Register the migration: `cell` IS derivable for every file
  written before C109, because the default was `full` throughout, so
  the derivation fills `"full"` honestly.  Then drop `cell` from
  `EXEMPT_RUN_SETTING_KEYS` so a manifest must resolve it -- the
  exemption in DESIGN 5.7 lasts exactly as long as `cell` is
  recorded nowhere, since the rule it is exempt from exists so that
  nothing EMITTED rides on an implicit default (VISION Principle
  11).  Recording it ends the exemption by definition.

  **Also the first real exercise of the C105 machinery**, which so
  far has an empty `SCHEMA_MIGRATIONS` table and a version gate that
  only ever refuses.  Doing it here would prove the migrate-when-
  honest path with a derivation that is genuinely derivable.  The
  one time-sensitive aspect: that derivation is honest only while
  the pre-C109 history is entirely `full`, so it gets murkier the
  longer mixed-cell harvesting runs -- though since every database
  can simply be regenerated at this stage, that is a mild argument
  rather than a deadline.

#### Seed-run refinement -- producer code tasks (design settled)

The four code tasks that follow implement the DESIGN decisions
taken in the 2026-07 seed-run refinement pass (companions to the
C108 cleanup and C109 cell-choice design items above).  Each
DESIGN rung is already written; these are the code that follows
it.  Reminder: `src/scripts` edits sync to `bin/` (the producer
imports its neighbours from `$IMAGO_BIN`).

- [x] C110. Make the default wingbeat re-apply the unit's
  imago-side settings on every launch.  Motivation: the Si seed
  run showed a re-dispatched `-loen -scf no` unit run a full
  ground-state SCF ("SCF after loen").  Root cause: `ImagoWingbeat`
  (`kaleidoscope/wingbeats.py`) reaches `imago.run_prepared(dir)`
  with NO settings, so `job` / `edge` / `scf_basis` -- imago
  *runtime* options that do not live in the staged `imago.dat`
  (DESIGN 6.2.10) -- are lost.  Fix, per the now-written PSEUDOCODE
  13.2 (Model A): the wingbeat always rebuilds the imago-side
  settings from the unit's options (`{k: v for ... if k in
  imago.OPTION_KEYS}` -> `ScriptSettings.from_options`) and passes
  them to `run_prepared(dir, settings=...)`.  The plan here
  assumed Model A would leave the wingbeat a SINGLE path (commit
  the staged inputs, then run), since the driver prepares every
  unit.  C111 found otherwise and did not build it: `ImagoWingbeat`
  is generic, and a client that never prepares arrives with no
  staging area, so the shipped `_stage_inputs` keeps three cases
  (commit a staged copy / nothing to do / build with makeinput).
  Settings are re-applied on every path regardless, which is what
  this item is actually about.  Merges with C111 at the wingbeat
  (C111 moves the build out to the driver's prepare pass; this
  makes the surviving run always carry its settings).  DESIGN
  6.2.2 / 6.1; PSEUDOCODE 13.2.  CODE.
  DONE (settings fix): shipped the always-pass-settings half --
  settings are built once and passed to `run_prepared` on both the
  prepared and build paths, keeping the current two-path structure.
  The single-path Model-A consolidation (drop the `_is_prepared`
  branch, add `commit_prepared_inputs`) lands with C111.  Test:
  `test_imago_runner_prepared_reapplies_imago_settings`.

- [x] C111. Key the run-reuse cache on `structure.dat` and run
  the prepare step in the driver (Model A).  Motivation: the seed
  run's "cache" never hit -- every re-run re-executed the SCF for
  one warm-restart iteration.  Root cause: the k-point convergence
  builder names its key file `"structure"`
  (`kaleidoscope/builders/kpoint_convergence.py`,
  `standard_key_fields` `KeyFile(name="structure", ...)`) but the
  wingbeat stages it on disk as `imago.skl`, so `cache_key_matches`
  never finds the file and reports a MISS.  Fix, per the now-written
  PSEUDOCODE (F2, Model A -- 15.6 / 13.4 / 13.2 and
  buildInitialPotentials Phase 1b; ARCH 8.5): (a) the key file
  becomes `structure.dat`, makeinput's *resolved* output, which
  bakes in the type/species assignment, basis, functional, and
  potential so an unlisted option cannot silently reuse a stale
  result; (b) `CalcUnit` gains a `prepared_dir` field; (c) the
  producer runs a driver-side PREPARE PASS before dispatch
  (partition the makeinput-side options, `build_run_dir` into a
  per-unit `prepare/<id>/*calc` staging area SEPARATE from the run
  dir -- the "must not clobber" rule -- plus a fast `imago -loen`
  when a solid's species/type assignment needs one), setting
  `prepared_dir`; (d) define `key_file_source` = the staged copy
  (`prepared_dir/<name>`) and `key_file_staged` = the run dir's
  copy, so the domain-agnostic cache core is UNCHANGED (it only
  byte-compares the two); (e) on a miss the wingbeat commits the
  staged inputs into the run dir (`commit_prepared_inputs`) and
  runs `run_prepared` -- no rebuild (merges with C110).  DESIGN
  6.2.5 / 5.7; ARCH 8.5; PSEUDOCODE 15.6 / 13.2 / 13.4.  CODE.
  DONE: standard_key_fields key file -> `structure.dat`;
  `CalcUnit.prepared_dir`; a driver-side `prepare_units` pass
  (Phase 1b) wired via a new `prepare_fn` injectable seam; the
  wingbeat's `_stage_inputs` (commit prepared / skip if staged /
  build) + `_commit_prepared_inputs` -- which also lands C110's
  deferred single-path consolidation.  First reconciled the
  PSEUDOCODE cache model to the code's `KeyFile(name, source)` and
  the 13.2 wingbeat to `stage_inputs`.  Tests: producer end-to-end
  stubs `prepare_fn`; new `test_imago_runner_commits_prepared_inputs`.

- [x] C112. Source the k-point grid-flatness threshold from the
  manifest `[harvest]` block, per atom, in eV.  Motivation: the
  seed run reported ia-3 and cmce non-converged though every SCF
  converged cleanly -- "not all converged."  Root cause:
  `guidance_harvest.py` reuses `scf_threshold` (1e-6 hartree, ~
  1e-8 relative) as the grid-flatness `metric_threshold`, far
  below real k-point sampling noise.  Fix, per DESIGN 7.8: divide
  each consecutive-pair total-energy delta by `cell_atom_count`
  (deltas in eV) and compare against `metric_threshold`, now
  resolved from the solid's `kpoint_convergence_threshold` (its
  own value, else the manifest `[harvest]` block, else the
  built-in 5e-4 eV/atom default; DESIGN 5.7).  Retire the
  `metric_threshold = scf_threshold` v1 convention (its docstring
  note in `guidance_harvest.py`) and thread the resolved harvest
  threshold through `pick_converged_unit`
  (`build_initial_potentials.py`), which today passes
  `scf_threshold`.  Also: teach the manifest loader to parse and
  per-solid-resolve `[harvest].kpoint_convergence_threshold`.
  Storage (Option B, decided 2026-07-08): keep the raw total-cell
  energies (hartree) in each guidance entry exactly as the run
  produced them, and do the per-atom eV conversion at every site
  that compares them against the threshold -- both `pick_converged`
  and `guidance_promote.py`'s `auto_promote_ok` -- reading
  `cell_atom_count` from the entry to normalize.  Single-source
  that conversion so the two comparison sites cannot drift.  While
  there, fix a pre-existing dimensional mismatch in
  `auto_promote_ok`: its flatness bar tests a statistical variance
  (an energy squared) against the plain threshold (an energy), so
  the two sides do not share dimensions regardless of units --
  make it a like-for-like comparison (a per-atom eV spread against
  the per-atom eV threshold).  Resolve target (PSEUDOCODE 11.4,
  now written): extend the shipped 2-arg `resolve_run_settings`
  into the 3-arg `resolve_settings(solid, defaults, harvest)`; add
  a `harvest` field on `CurationManifest` and a
  `kpoint_convergence_threshold` field on `ReferenceSolid`; and
  have `apply_manifest_defaults` pass `manifest.harvest` so the
  harvest arm resolves the solid's own value, else `[harvest]`,
  else built-in `DEFAULT_KPOINT_CONVERGENCE_THRESHOLD`.
  DESIGN 7.8 / 5.7.  CODE (+ PSEUDOCODE).
  DONE: `curation_manifest` (DEFAULT_KPOINT_CONVERGENCE_THRESHOLD,
  HARVEST_SETTING_KEYS, `[harvest]` parse, `ReferenceSolid`.
  `kpoint_convergence_threshold`, `CurationManifest.harvest`,
  `resolve_run_settings` -> 3-arg `resolve_settings`);
  `build_initial_potentials` (`apply_manifest_defaults` 3-arg,
  prediction-record stamp, `pick_converged_unit` per-atom via a
  `load_structure` cell-atom-count read); `guidance_harvest`
  (`per_atom_ev` helper on `HARTREE`, `pick_converged` 3-arg,
  `_load_structure` -> public `load_structure`, `build_entry` takes
  the loaded structure + `kpoint_threshold`, scf_threshold guarded
  as a required context fact); `guidance_promote.auto_promote_ok`
  (per-atom eV SPREAD, not variance).  333 tests pass across the
  affected suites.

- [x] C113. Let the producer/orchestrator run as its own batch
  job, with a separate orchestrator resource block.  Motivation:
  the driver now does real per-unit prepare work (a makeinput
  build, plus a fast `imago -loen` when assignment needs it, once
  per unit including cache hits, C111), which at scale would tie
  up a login node's terminal for a whole flight.  Fix, per DESIGN
  6.2.11 ("Driver location"): (a) support wrapping the run in its
  own `sbatch` job; (b) add a separate orchestrator resource
  block -- sized to the dispatch shape (tiny when it only fans
  out to worker jobs, compute-sized under `--dispatch local`) --
  as a new `clusterrc.py` / per-run setting distinct from
  `memory_per_worker`; (c) materialize-then-submit: run
  `--materialize-only` on the login node first (the only
  network-touching step; compute nodes may lack internet), THEN
  submit the batch job whose prepare + dispatch touch no network;
  (d) keep `--dispatch` a per-run flag (`local` for seed scale
  now, `slurm-per-job` / `slurm-pooled` later).  Deferred
  sub-item, turned on only when the serial prepare cost bites:
  move prepare-and-hit-test onto dispatched worker units (a hit
  then costs a cheap worker slot instead of being decided
  driver-local).  DESIGN 6.2.11.  CODE.
  DONE: clusterrc gains a grouped `orchestrator` block (cores /
  memory / walltime) -- the driver job class, distinct from the
  per-worker sizing (ARCH 9.4 note added; two job classes, not one
  per builder); `cluster_config.build_orchestrator_sbatch` renders
  the sbatch script (one node, the orchestrator shape, account +
  partition, worker_init, then the command);
  `build_initial_potentials --submit` materializes on the login node
  then submits a batch job that re-runs the producer minus --submit
  (structures already cached).  Tests cover the generator, the
  clusterrc block, and the submit path (sbatch mocked).  The
  deferred prepare-on-workers sub-item stays deferred.  Actual
  sbatch submission + running-as-batch-job is validated by a live
  cluster run (rides with C75).
  AMENDED 2026-07-09: C113 originally shipped code with NO
  PSEUDOCODE above it -- its TODO entry was detailed enough to
  implement from, so the level was skipped (the trap now named in
  CLAUDE.md, "Chain Discipline").  A `/refine` caught it.  The
  back-fill wrote PSEUDOCODE 13.7 (the `orchestrator` block, the
  derived memory request, `build_orchestrator_sbatch`) and 11.4
  (`main_submit_mode`), and verifying the code against DESIGN
  6.2.11 first turned up two defects the original pass had missed:
  `build_orchestrator_sbatch` prefixed a second `#SBATCH` onto each
  already-complete `extra_scheduler_options` line, and
  `cluster_probe._starter_schema` never learned the `orchestrator`
  key, leaving `test_cluster_probe` RED on `main` since 10c1741
  (the C113 test run covered only the suites it assumed affected).
  Both fixed, each with a test; full non-integration suite green
  (844 passed).  The `memory_per_node` / `memory_per_worker` split
  (shipped in 6e8db47, DESIGN + code only) was back-filled into
  13.7 in the same pass.

- [x] C114. Sync the PSEUDOCODE writer to the shipped `[defaults]`
  manifest.  Surfaced 2026-07-08 during the seed-run refinement
  `/refine`: the `[defaults]` block shipped in code (C104 --
  `curation_manifest.py` resolution and the `expand_manifest` /
  `cod_fish` writers), and PSEUDOCODE 11.4 (the reader) is now
  updated to model both `[defaults]` and `[harvest]` resolution,
  but PSEUDOCODE 11.6 `format_manifest` (the *writer*) still emits
  every run setting per solid with no `[defaults]` block -- so the
  writer pseudocode lags its own shipped code.  Work: update
  `format_manifest` to emit the `[defaults]` block plus the compact
  per-solid overrides the shipped writer produces, so the read and
  write sides of the schema library agree in pseudocode as they
  already do in code.  Open sub-question, decide when C112 lands:
  whether the authoring tools also emit an explicit `[harvest]`
  block (making the tolerance visible and editable) or omit it and
  lean on the built-in default -- the harvest setting is exempt
  from the write-it-down rule (it has a default and the resolved
  value is recorded on each guidance entry), so either is valid.
  PSEUDOCODE 11.6; DESIGN 5.7.
  DONE: PSEUDOCODE 11.6 rewritten against DESIGN 5.7 -- the four
  top-level blocks each gain a writing counterpart, solids carry
  only the run settings they override, `preferred` is never
  serialized (it is structural, recovered from the block a
  declaration lands in), and `expand_manifest` is re-specified as
  the `shared_defaults()` + `sparse_solid()` hoist the code
  actually uses.  The open sub-question is SETTLED: the authoring
  tools populate no harvest dict and emit no `[harvest]` block
  (leaning on the built-in default), but a curator-authored block
  must survive -- and did not.  Writing the pseudocode exposed a
  live bug: `format_manifest` never emitted `[harvest]` at all, so
  an authored `kpoint_convergence_threshold` was silently dropped
  on rewrite and the next load fell back to the built-in, doubling
  the convergence tolerance with nothing on screen.  Fixed with
  three tests (one fails as `assert 0.0005 == 0.00025`).

- [x] C116. The k-density ladder is not a refinement sequence, and
  `pick_converged` assumes it is.  A requested k-point density does
  not map monotonically onto the mesh imago actually integrates
  over: raising the density can leave the mesh unchanged, and can
  even coarsen it.  The two-sided flatness test in DESIGN 7.8 step
  3c reads an unchanged mesh as a converged energy, because two
  runs of the same calculation differ by exactly zero.  This is a
  DESIGN defect in the acceptance rule, not a bug in any script.

  Measured 2026-07-09/10 from the first successful end-to-end seed
  run (orchestrator job 14882026).  `gs_scf-fb.out` prints the
  resolved irreducible k-point list under the heading `Kpoints in
  x,y,z cartesian form in recip space cell are:`, so the mesh each
  rung actually used can simply be counted.  Across the eight
  ladder rungs (densities 25, 50, 100, 150, 200, 250, 300, 400)
  the irreducible k-point counts are:

      si_fd-3m_227_*      3    6   12   24   20   30   60   45
      si_ia-3_206_2016    2    4    7   18   18   24   32   48
      si_cmce_64_1999     1    4    2    4    8   12   12    9

  Read the diamond row across: density 150 gives 24 points and
  density 200 gives 20.  Asking for a finer sampling returned a
  coarser one.  The same reversal appears from 300 to 400 (60 down
  to 45), and in the Cmce row from 300 to 400 (12 down to 9).  The
  duplicated rungs -- `si_ia-3` at 150 and 200, `si_cmce` at 250
  and 300 -- are the degenerate case of this, not a separate
  phenomenon: in each pair the printed k-point lists are identical
  point for point, so the two runs are one calculation performed
  twice.  Every energy delta up such a ladder therefore mixes real
  convergence with mesh-shape noise, which is why the deltas
  alternate in sign instead of decaying.

  This has already corrupted a shipped potential.  `pick_converged`
  accepts the smallest interior rung whose energy is within the
  threshold of BOTH neighbours.  `si_ia-3_206_2016` was accepted at
  density 200 on a left-hand delta of exactly 0.000 meV/atom -- the
  duplicate mesh compared against itself.  Its right-hand deltas
  then GROW (0.131, 0.110, 0.732 meV/atom) and its gap is still
  falling hard, 0.507 -> 0.287 eV between densities 200 and 400.
  The solid is not converged; the rule declared it converged on an
  arithmetic tie, and its potential is in
  `share/atomicPDB/si/s_gaussian_pot.toml` today.  DESIGN 7.8
  reasons carefully about false plateaus from a single numerical
  dip, and the two-sided test defeats exactly that -- but it
  assumes each grid point is a distinct calculation, and a repeated
  rung manufactures a perfect zero on one side for free.

  Ruled out by the measurement, and recorded so they are not
  re-litigated: `MIN_KP_LINE_DENSITY` is neither ignored, clamped,
  nor saturated (the collapse happens at DIFFERENT density pairs
  for different lattices -- 150/200 for Ia-3, 250/300 for Cmce --
  whereas a clamp would bite at the same density for every solid);
  and no run directory was reused or stale (same reason, plus the
  fd-3m rows show eight distinct meshes).  What remains, and does
  NOT need settling before the DESIGN work, is whether the collapse
  happens in the density-to-mesh map or in the IBZ reduction that
  follows it.  Point-op counts differ widely across these solids
  (Fd-3m 48, Ia-3 24, Cmce 8), so both remain live.

  The DESIGN question, which is where this item lives.  A ladder
  indexed by requested density asks a question the code cannot
  answer, because density is not what the calculation consumes.
  The candidate fix is to index the ladder by the resolved mesh
  instead: generate the rungs, resolve each to its k-point set,
  drop or merge rungs that resolve alike, and only then apply the
  flatness test -- so that consecutive rungs are guaranteed to be
  genuinely different calculations and a zero delta means a
  converged energy rather than a repeated one.  Whether the ladder
  should additionally require monotonically increasing k-point
  counts, and what to do with a rung that coarsens, are open.  Do
  not write code before DESIGN 7.8 answers this.

  Adjacent, and now cheap: DESIGN 8.2's size signature already
  wants `kpoint_count`, "number of k-points actually computed."
  That number is the same one this item needs, and it is already
  printed -- harvesting it into `result.toml` would let the
  acceptance rule see the mesh directly instead of inferring it
  from the requested density.  DESIGN 7.8 (defect) and 8.2
  (consumer); then PSEUDOCODE 11.4 and code.  Cf. C117, whose
  near-metal oscillation rides on top of this same noise.

  **CLOSED 2026-07-20 -- superseded, not patched.**  The candidate
  fix above (index the ladder by resolved mesh) was overtaken by
  the answer C118 gave: retire the requested-density ladder
  altogether and climb through symmetry-compatible MESHES, which
  removes the defect at its root because every rung is a mesh by
  construction and no two consecutive rungs can resolve alike.
  The adjacent ask landed too -- `result.toml` now carries the
  mesh and the gap the calculation actually used, so the
  acceptance rule reads them directly rather than inferring
  anything from a requested density.  Two live seed runs confirm
  it end to end (jobs 15026798 and 15032332).

- [x] C117. Decide the k-point integration scheme for near-metallic
  reference solids.  Observed 2026-07-09 in the same seed run:
  `si_cmce_64_1999` (16-atom Cmce cell) is the one solid of eight
  that did not converge, and its gap oscillates between metal and
  small-gap semiconductor along the k-density ladder -- 0.328,
  0.000, 0.294, 0.093, 0.000, 0.093, 0.093, 0.086 eV at densities
  25 through 400.  The manifest's `[defaults]` set
  `kpoint_integration = "gaussian"`, and a bare `gaussian` token
  carries no width, so makeinput applies NO smearing (DESIGN 5.7).
  With no smearing and a near-zero gap, occupations jump each time
  a band crosses the Fermi level, which is a plausible source of
  the 100+ meV/atom energy swings at low density.  Not the only
  source, though: C116 shows this solid's mesh does not refine
  monotonically with the requested density (1, 4, 2, 4, 8, 12, 12,
  9 irreducible points up the ladder), so part of the swing is the
  mesh changing shape rather than the occupations jumping.  The two
  effects compound, and C116 must be understood first, or a smearing
  width will be chosen to damp noise that smearing did not cause.

  This is a curation question, not a bug: the schema already
  supports a per-solid override, so `linear-tetrahedral`
  (parameter-free) or `gaussian-0.1` can be named on this solid
  alone.  Open: whether the seed manifest should instead default
  to `linear-tetrahedral`, which DESIGN 5.7 calls the producer's
  default, and why the Si defaults were settled on `gaussian`.
  Settle before re-running the seed, since the choice selects the
  guidance predictor's sub-model and a database must not mix
  sub-models silently.  DESIGN 5.7 / 7.6; CURATION.

  **DECIDED 2026-07-20: the default stays bare `gaussian` -- no
  smearing -- and near-metals are handled by classification, not
  by broadening.**  C116's half of the compounding was removed
  first, exactly as this entry demanded, so what remained could be
  attributed cleanly.  Smearing was then tested directly on
  si_cmce rather than argued about: three live runs at widths
  `0.0`, `0.026` (room-temperature kT) and `0.1` eV.  The wobble
  barely moved -- quadrupling the width from `0.026` to `0.1`
  reproduced nearly the same per-step energy changes, including
  the same `+7.4x` reversal at mesh `[4,6,7]` -- so a wider
  Gaussian is not the lever, and the occupation-jump hypothesis
  above is not the dominant term.  Smearing also costs accuracy
  where it is not needed: against the clean run, wide-gap diamond
  Si shifted by `-0.02` meV/atom (noise) but small-gap si_ia-3
  carried a real `+1.40` meV/atom bias at `0.1` eV.  Leaving the
  default off therefore keeps every insulator exact, and C125's
  gap test disposes of the metals the smearing was meant to tame.
  The per-solid `gaussian-<width>` override remains available for
  a curator who wants it; the seed manifest uses none.

- [x] C118. Implement the adaptive mesh climb (DESIGN 3.12 /
  PSEUDOCODE 4e).  The convergence search moves from a fixed
  density grid to a climb through symmetry-compatible meshes,
  seeded by the guidance prediction and stopped when the energy
  is flat.  This is the fix for the 0/8 seed regression: a fixed
  density ladder collapses onto too few distinct meshes for
  high-symmetry cells (cubic Si tops out at [5,5,5] still moving
  ~1.7 meV/atom), while the climb keeps going until flat.  Scope:
    - Producer (`build_initial_potentials.py`): the round-based
      climb loop (`converge_by_climb`, 4e.5) -- serial within a
      material, parallel across -- replacing the one-shot verify
      grid.  Rung mechanics (`climbOneRung` / `descendOneRung`,
      4e.1) reuse `selectAxialCounts` / `spacingSpread` (4c.2);
      the stop test with confidence-scaled persistence (4e.2) and
      the ceiling; the two dispatch modes gated by confidence
      (4e.4).
    - Iterative dispatch: the producer drives kaleidoscope round
      by round, reading energies between rounds (ARCH 9.7).  The
      dispatch core stays dumb (Principle 12) -- each round is a
      flat CalcUnit list; no dispatch-side change.
    - Guidance schema: `Verification` gains `converged_mesh` (the
      resolved axial counts); the reader (15.3), emitter (15.4),
      and `build_entry` (15.7) carry it; the recorded density is
      the converged mesh's full-mesh volume density (4e.6).
    - Config, not constants (Principle 11): the ceiling, the
      confidence thresholds, the grid width, the start offset,
      and `flat_needed` live in the manifest characterization
      block (like `metric_threshold`) or the site rc for the cost
      budget -- never hardcoded (DESIGN 3.12.6).
  The numeric values of those knobs are to be fixed by experiment
  on the seed set (3.12.6), so a first pass can carry documented
  provisional defaults in the config.  Presupposes the Stage 1-4
  mesh rework (done).  After it lands, re-run the seed to close
  C116 and unblock C117.  CODE; DESIGN 3.12; PSEUDOCODE 4e;
  ARCH 9.7 / 10.

  Increments (implementation ordering; each governed by the named
  PSEUDOCODE section, which is where the specification lives):
    - [x] Inc 0. Chain gate: PSEUDOCODE 4c.7 (producer-side
      axis-class sourcing) and the 4d.5 amendment (imago emits
      `RESOLVED_KP_CLASSES` as the port's validation hook).
    - [x] Inc 1. `src/scripts/mesh_climb.py` primitives (axis
      classes 4c.1 / 4c.7, count selection 4c.2, rung mechanics
      4e.1) plus their unit tests (`test_mesh_climb.py`, 27 tests).
    - [x] Inc 2. The pure decision helpers.  Stop test with
      confidence-scaled persistence (`pick_converged_climb`, 4e.2)
      beside `per_atom_ev`/`pick_converged` in `guidance_harvest`
      (single-sourced); the per-axis ceiling (`at_ceiling`, 4e.2),
      first-round seeding (`initial_meshes`, 4e.4), and the
      confidence-to-mode policy (`resolve_climb_policy`, 4e.4) in
      `mesh_climb`, with provisional threshold defaults.  Tests:
      +9 in `test_mesh_climb.py`, +4 in `test_guidance_harvest.py`.
      `climbAction` moved to Inc 3: ARCH 9.7 places the energy-
      reading / next-mesh decision in the producer, with the loop.
    - [x] Inc 3a. The producer's climb control loop, in
      `build_initial_potentials.py` (ARCH 9.7): `climb_action`
      (4e.3) and the `converge_by_climb` round-based loop (4e.5),
      with the `dispatch_round` runner INJECTED, plus the `Rung` /
      `ClimbConfig` / `ClimbAction` types and `_sort_by_mesh` /
      `_merge_distinct` helpers.  The ceiling-tag call (7.8 step
      3d) is an injected `on_non_converged` hook.  Tested with a
      synthetic `dispatch_round` (+9 in
      `test_build_initial_potentials.py`).
    - [x] Inc 3b-design. Mesh-dispatch design + pseudocode (the
      chain gate for 3b-code).  DESIGN 7.7 gains the mesh->run
      mechanics (explicit `scfkp` mesh, `kpt-mesh` calc tag,
      `total_energy`/`kpoint_mesh` read-back, the round adapter,
      the builder split, and the fail-fast rule); PSEUDOCODE 4e.7
      adds `encodeMeshValue`/`decodeMeshValue`, `build_mesh_unit`,
      `predict_kpoint_density`, and `make_dispatch_round`; 4e.5 is
      amended for fail-fast (round-0 empty + a requested rung that
      does not return) and reconciled to `seed_densities` /
      `on_non_converged`; 4e.4 reconciled to seed from a density.
    - [x] Inc 3b-code. Mesh-dispatch adapter per 4e.7 / DESIGN 7.7.
      Builder split in `kpoint_convergence.py`: `encode_mesh_value`
      / `decode_mesh_value`, `build_mesh_unit` (`scfkp` option +
      `kpt-mesh` tag), `predict_kpoint_density` (predict-only).
      `make_dispatch_round` in the producer (injected prepare /
      dispatch / completed / read; omits non-completed units;
      asserts the mesh is honoured exactly).  `converge_by_climb`
      gained the fail-fast guard (round-0 empty + missing-rung ->
      NON_CONVERGED).  +7 builder tests, +5 producer tests; 221
      affected-suite tests green.
    - [x] Inc 4. `converged_mesh` through the guidance schema.
      Propagated the field from DESIGN 7.2 into PSEUDOCODE (the
      `Verification` dataclass, the reader 15.3, the emitter 15.4,
      `build_entry` 15.7) and simplified `record_converged` (4e.6)
      to `(rung, rungs, config)`.  Code: `Verification.converged_mesh`
      + reader (optional, three-count check) + emitter (inline int
      array) in `guidance_db.py`; `build_entry` reads it from the
      chosen rung's `result.toml` in `guidance_harvest.py`;
      `record_converged` in the producer.  +7 tests.
    - [x] Inc 5. The climb-policy knobs sourced from the manifest
      (Principle 11).  DESIGN 5.7 gains the optional
      `[harvest.kpoint_climb]` sub-table (7 knobs, database-wide,
      each with a provisional default); PSEUDOCODE carries the
      `KPOINT_CLIMB_KEYS` validation and the confidence-policy
      prose.  Code: `curation_manifest` parses / validates / emits
      the sub-table (`HARVEST_SETTING_KEYS` + `KPOINT_CLIMB_KEYS`);
      `mesh_climb.climb_policy_from_manifest` merges it over the
      provisional defaults into `(PolicyThresholds, max_count)`.
      A test pins `KPOINT_CLIMB_KEYS` to mesh_climb so they cannot
      drift.  Wiring into the producer main is Inc 6.
    - [ ] Inc 6. Wire the producer main to the climb.  The 6a-6c code
      is DONE; only the live seed re-run remains.  **6a** (DONE):
      DESIGN 5.7 + PSEUDOCODE 11.4 rewritten for the climb; the
      chosen-facts `build_entry` (Q1/Q2 resolved); an A'' diversion
      added `symmetry.py` -- the shared `share/spaceDB` point-ops
      reader both the kp writer and the producer read through (ARCH
      2/3/7 started the Section-7 symmetry split in the small).
      **6b** (DONE): rewired `build_initial_potentials` main to the
      three-phase climb (build / converge / harvest) with `build_-
      climb_config`; FULL-SWEEP retirement of `build_kpoint_-
      convergence` + its helper tail + `pick_converged_unit` across
      code / tests / PSEUDOCODE 15.6 / DESIGN; orchestration tests
      rewired; the make_run_log_entry mesh-vs-density fix.  **6c**
      (DONE): imago emits `RESOLVED_KP_CLASSES` (4d.5) from its own
      `computeAxisClasses` call in *both* mesh style codes -- every
      climb run is an explicit mesh (style 1, `-scfkp`), and imago
      must compute the classes independently of the producer or the
      cross-check against `axis_classes_for_cell` is circular.  The
      separate axis-class self-test (4c.7) is descoped: the live seed
      re-run exercises the port on real cells.  **6d** (DONE): two
      blockers the first live run surfaced.  The script install list
      omitted `mesh_climb.py` and `symmetry.py`, so the installed
      `structure_control` raised `ModuleNotFoundError`.  And the
      producer reloaded one Parsl `Config` per climb round, which a
      single-use executor forbids: PSEUDOCODE 13.5 gained
      `make_executor`, and 11.4 / 4e now build one executor per run,
      pin it to every dispatch, and close it once in a `finally`
      (DESIGN 6.2.11's pooled shape; 6.2.3 <-> 3.12.5 now cross-
      referenced, since the climb is a producer-side control loop by
      Principle 12 rather than a wingbeat inner loop).  NEXT: the
      live seed re-run against the rebuilt binary closes C116 and
      unblocks C117.

  **DONE.**  The live seed re-runs that this entry was waiting on
  have all been made, and the climb has since been tuned (C123),
  floored (C124), and taught to recognise metals (C125) against
  what they showed.  The final clean production run (job 15032332)
  converged or settled all eight reference solids and harvested
  every one, which is the outcome the climb was built to reach
  after the fixed density grid returned 0/8.  C116 is closed above
  and C117 is decided above, as this entry predicted they would
  be.

- [x] C119. A block asks for its slice, not for the node
  (DESIGN 6.2.11 / PSEUDOCODE 13.7 / code).  Surfaced by the
  C118 seed re-run, which held a 128-core node for 41 minutes
  and used one core of it (0.8%), and left seven of its eight
  blocks queued behind whole-node requests they never needed.
  Two causes, both ours: the provider never named `exclusive`,
  so Parsl's whole-node default stood; and `scheduler_options`
  requested memory but never cores, which only "worked"
  because exclusivity handed us the node anyway.  The two are
  one request -- declining the node without naming the cores
  takes SLURM's one-core default and would starve a packed
  pool -- so DESIGN 6.2.11 now states both, deriving the cores
  as `cores_per_worker x workers_on_the_node` exactly as the
  memory is derived, and pointing a site that truly needs
  whole nodes at `extra_scheduler_options` rather than adding
  a knob.  Verified against the real `clusterrc`: per-job asks
  `--mem=10G --cpus-per-task=1`, pooled asks `--mem=400G
  --cpus-per-task=40`, both non-exclusive.  The partition
  permits sharing (`OverSubscribe=NO`, `ExclusiveUser=NO`
  bound cores, not nodes), so the seven blocks that queued had
  nothing to wait for.  Note this is throughput and courtesy,
  not correctness: the seed run's answers were right, and its
  41 minutes were 40 minutes of real compute -- the waste was
  the 127 idle cores no one else could use.

- [x] C120. Retire the climb's round barrier -- CODE the
  wait-for-any climb (option a).  The design is SETTLED and
  written top-down; code against these sections, NOT against
  this entry: DESIGN 3.12.5 (concurrent, no chain waits on
  another), 6.2.3 (two-phase public surface -- send-off /
  collect, `done()` on futures), 7.7 (the climb dispatcher);
  PSEUDOCODE 13.5 (`send_off` / `collect` / `collect_next`,
  `dispatch` the wrapper), 4e.5 (`converge_by_climb` rewritten),
  4e.7 (`make_climb_dispatcher`), 11.4 (producer main);
  ARCHITECTURE 9.7.  Work to do (a tracking checklist, not a
  spec):
  - `dispatch.py`: expose `send_off` and `collect`, add
    `collect_next` (poll `done()`, take whichever lands first),
    make `dispatch` the send-off-then-collect-all wrapper; add
    `done()` to `_LocalFuture` (always True) and `_ParslFuture`
    (delegate to the AppFuture); a cache hit becomes an
    already-done future so hits and misses sit uniformly.
  - `build_initial_potentials.py`: rewrite `converge_by_climb`
    to the concurrent shape (per-material `in_air` + `opening`
    sets; judge a chain the moment its own rung lands); replace
    `make_dispatch_round` with `make_climb_dispatcher` (send /
    next_rung over ONE flight whose unit list accretes, so
    `flight.toml` records every rung and the per-send race is
    gone).
  - tests: the 4 orchestration tests + the dispatch-layer tests.
  Decisions locked: one accreting flight; a dispatcher object
  owning the in-flight set; hits as done futures; `on_outcome`
  now fires in LANDING order.

  Reclamation -- freeing a retired chain's workers for the
  chains still climbing -- is DESIGNED but NOT built here; it
  waits on a parallel imago (DESIGN 6.2.11 forward note;
  ARCHITECTURE 9.8).  Option (a) only keeps the door open: the
  producer tracks chains not slots (the `in_air` map does this),
  a unit MAY carry a resource field later (not added now), and
  the rung width must never enter the cache key.

  Measured payoff, so nobody re-derives it: the barrier cost
  ~2 min of ~34 in the C118 seed re-run -- small only because
  `si_cmce` was 82% of compute (1949 s of 2385 s, 28 serial
  rungs against 8-9 for the others).  It will bite when
  materials are heterogeneous rather than seven-of-eight
  identical.  **The step-size rule is the larger prize** (C118
  inc 6, C116).

  **DONE, committed 5cb6e0d.**  Every item on the checklist above
  shipped as written, and the concurrent shape has since carried
  five live seed runs without a barrier stall.  Reclamation is
  still designed-but-not-built, as this entry scoped it; the door
  stayed open, since the producer tracks chains rather than slots
  and the rung width never entered the cache key.

- [x] C121. A loen-descriptor build must skip the fingerprint
  match (makeinput recursion fix).  DESIGN 5.6.5 + 5.10.2 first,
  then PSEUDOCODE (the potential-resolution function), then the
  code in `makeinput.py` (`_obtain_pot_info`).  Chain-compliant:
  DESIGN + PSEUDOCODE before the guard is coded.

  **The failure mode.**  makeinput builds the loen input -- the
  provisional `imago.dat` a loen run reads to compute a
  bispectrum descriptor -- by running `makeinput -loeninput`.
  That build runs the 5.6.5 fingerprint match like any other, and
  in the file-dictated regime the match asks the database for its
  preferred fingerprint; when that is bispectrum, computing the
  query runs a loen descriptor computation
  (`loen_site_descriptors` in makegroups), whose first step is
  another `makeinput -loeninput`.  So the build invokes itself,
  each level nesting a `loen_pick_work/` inside the last until the
  path overflows (`OSError Errno 36, File name too long`).  The
  per-process descriptor cache never helps because every level is
  a fresh subprocess.

  **Why it is a second-run failure.**  An empty database has no
  preferred bispectrum record, so the match returns early and no
  loen run is triggered -- the first seed run into a fresh
  database succeeds and harvests potentials WITH `preferred =
  true` bispectrum fingerprints.  Every later run reads that
  populated database and takes the recursive path.  (Confirmed:
  the first success stamped `share/atomicPDB/si/s_gaussian_pot`
  with preferred bispectrum records; the next run recursed at the
  loen pre-flight.)

  **The fix.**  A loen-descriptor build -- a makeinput run given
  `-loeninput` -- always skips the fingerprint match and takes the
  default-tagged entry (5.6.5 step 3).  This is correct, not a
  compromise: the bispectrum is geometric, so the potential the
  loen run sees is irrelevant to the descriptor it produces
  (5.10.2).  The guard belongs in `_obtain_pot_info`, keyed on
  `-loeninput`, alongside the existing `-nofingerprint` and `-pot`
  skips.  `-loeninput` has exactly two callers and both are
  loen-descriptor builds -- the producer's loen pre-flight unit
  (`build_initial_potentials.py`) and the makegroups descriptor
  sub-run (`makegroups._run_makeinput`) -- so one guard on that
  flag closes both entry points at the single place the recursion
  is born.

  **DONE, committed f3ec42a**, chain-compliant (DESIGN then
  PSEUDOCODE then the guard).  The second-run failure it describes
  is gone: the Si database now holds preferred bispectrum records
  from the first seed run, and every run since has read that
  populated database and built its loen inputs without recursing.

- [x] C123. Code the bracket-then-refine mesh climb (the
  step-size rule, option 3).  DONE against DESIGN 3.12.2/3/5/6 and
  PSEUDOCODE 4e.1-4e.6, /refine-clean.  Shipped in three
  increments:
  - `mesh_climb.py`: the stride/fill/ceiling/trace primitives
    (`climb_n_rungs`, `next_fill_mesh`, `ceiling_mesh`,
    `consecutive_block`, `rung_at`); the three modes
    (PARALLEL_GRID / BRACKET_REFINE default / UNIT_STEP option);
    `max_stride` + `climb_shape` policy knobs.  The energy tests
    (`stride_is_flat`) live in `guidance_harvest` beside
    `pick_converged_climb`, single-sourced on `per_atom_ev`.
  - `build_initial_potentials.py`: the bracket-refine state
    machine (`bracket_refine_next` + `_stride_up` + `_enter_refine`
    + `new_bracket_refine_state` / `new_search_state`),
    `climb_next` dispatch, `climb_action` returns CONVERGED(rung),
    `converge_by_climb` threads the per-material search state, and
    `record_converged` records the consecutive block.  Lives in the
    producer (not `mesh_climb`) because it reads energies
    (ARCHITECTURE 9.7).
  - `curation_manifest.py`: `max_stride` + `climb_shape` join
    `KPOINT_CLIMB_KEYS`; `climb_shape` value validated against the
    known shapes.

  DESIGN fix found while coding: the refine must fill `flat_needed
  + 1` rungs above the flat stride's bottom, not one -- the bottom
  rung need not itself be settled (its lower neighbour may still be
  moving), so the first confirmable rung is often one higher.
  Without it a `flat_needed = 2` climb never converged (endless
  three-rung brackets).  Corrected across DESIGN 3.12.3 +
  PSEUDOCODE 4e.3 + code, with a regression test.

  Remaining (not code): provisional knob defaults (`max_stride`,
  the confidence thresholds) are still to be fixed by the seed
  experiment (3.12.6); the live seed re-run validates the climb
  end-to-end.  The payoff is the cold/seeding case -- si_cmce's
  ~28 rungs should fall to ~8-10; a warm seed barely brackets.

- [x] C124. Code the crystalline opening floor and retire the
  coarse-mesh gate on the near-metal bail, per DESIGN 3.12.4 /
  PSEUDOCODE 4c.2 + 4e.3 + 4e.4.  Found from the seed re-run:
  si_ia-3 (BC8 insulator) false-bailed NON_CONVERGED because its
  cold climb opened at `[1,1,1]` and the `[1,1,1] -> [2,2,2]`
  Gamma-sampling rise tripped the near-metal bail; the point-count
  gate meant to suppress it did not.  The fix floors a crystalline
  climb's opening at a per-axis cap (`crystalline_floor_axis_count`,
  default 4: densest axis gets the cap, others scale down), so the
  climb never visits the unreliable coarse regime and the bail
  needs no gate.  Surface (all bounded, inline):
  - `mesh_climb.py`: factor `counts_at_spacing` out of
    `select_axial_counts`; add `crystalline_floor_mesh`; swap the
    `metallic_min_points` PolicyThresholds field / default /
    validation for `crystalline_floor_axis_count`; `initial_meshes`
    gains an `opening_floor` arg and applies the max in the climb
    branch.
  - `build_initial_potentials.py`: `ClimbConfig` field
    `metallic_min_points` -> `opening_floor`; `build_climb_config`
    computes it (crystalline -> `crystalline_floor_mesh`, else
    None); drop the coarse-mesh guard in `bracket_refine_next`; the
    seeding call passes `config.opening_floor`.
  - `curation_manifest.py`: `KPOINT_CLIMB_KEYS` swaps the knob.
  - tests: replace the coarse-gate test with a floor test (keep
    `_gamma_noise_energy`); `_cubic_config` + monkeypatched
    `ClimbConfig` builders take `opening_floor`; assert
    `build_climb_config` floors a cubic cell at `[4,4,4]`, an
    anisotropic one lower per axis, and a non-crystalline one to
    None.  Committed e896eac.  The live seed re-run (2026-07-20,
    job 15019795) confirmed the floor's primary purpose: every
    insulator opened at its floor and converged -- si_ia-3 back
    to ~`[6,6,6]` with the false bail gone, and all six `fd-3m`
    diamond-Si allotropes at ~`[6,6,6]`.  It also overturned the
    prediction above that si_cmce "still bails cheaply": removing
    the coarse gate exposed a near-metal *runaway* instead,
    captured and analysed as C125 below.

- [x] C125. Close the near-metal dead zone the C124 seed re-run
  exposed -- RESOLVED, by classifying metals directly on the gap.
  CODE-COMPLETE down the whole chain (DESIGN 3.12.3/3.12.4/3.12.6/5.7
  -> PSEUDOCODE 4e.2/4e.3 + knobs -> code); 991 non-integration tests
  pass.  Code: `guidance_harvest.is_gapless` (replaces `stride_rose`);
  producer `_ACTION_METAL` settles + RECORDS the rung (only `ceiling`
  is non-converged now), gap check atop `bracket_refine_next`, `Rung`
  gains `gap` (from `gap_ev`), `build_climb_config` passes the knob
  directly; `metallic_rise_multiple` -> `metal_gap_threshold` (eV,
  default 0.05, checked `> 0`) across mesh_climb / curation_manifest.
  Committed 61429b9.  DEFERRED (minimal-for-now): a metal is recorded
  as a plain rung, no distinguishing flag -- the harvest (7.8) /
  dataspace (7.2) metal-flag + predictor-learning question is left
  for later, and should be settled before the guidance predictor is
  ever trained on a database that contains metals.

  **VALIDATED LIVE 2026-07-20**, twice.  Job 15026798 confirmed the
  fix directly: si_cmce read a gap of `0.0` at its very first rung,
  was declared a metal, and settled at the floor mesh `[2,4,4]` after
  a single rung instead of the 21 it had run away through -- and it
  was RECORDED, taking the producer from 6/8 to 8/8 harvested.  The
  seven insulators were untouched, converging at `[6,6,6]` as before,
  which is the outcome the gap threshold was chosen to give.  Job
  15032332 then repeated it as a clean production run with smearing
  off, harvesting 8/8 into `share/atomicPDB/` (see C91).
  C124 worked as intended and, in doing so, removed the
  *accidental* early stop that used to terminate the si_cmce
  near-metal -- the coarse `[1,1,1] -> [2,2,2]` Gamma-artifact
  rise, which tripped the near-metal bail for a reason unrelated
  to metallicity.  With that gone, si_cmce (42-atom Cmce, floor
  `[2,4,4]`) climbed 21 rungs to `[6,10,11]` without ever
  converging or bailing, and was cancelled only because it was
  still short of the per-axis ceiling (`DEFAULT_MAX_COUNT = 20`).

  **Root cause -- a dead zone between the two stop tests.**  Along
  the refining (finer-mesh) direction si_cmce's energy oscillates
  with amplitude ~1-3x the convergence threshold (`5e-4`
  eV/atom): too rough for the flatness test to latch (it needs a
  stride below 1x) yet too gentle for the near-metal bail to fire
  (it needs a finer-raises-energy stride above 50x = `0.025`
  eV/atom).  The large `+135x` / `+25x` strides in its trace are
  bracket-fills to *coarser* meshes, which the one-sided bail
  rightly ignores.  So neither test bites and the climb runs to
  the count ceiling, dragging through exactly the expensive
  high-mesh rungs the bail is meant to avoid.

  **This falsifies a specific DESIGN 3.12.3 claim** (~lines
  1898-1902: the bail "saves the expensive high meshes a
  near-metal would otherwise be dragged through -- si_cmce's
  ladder ran to thousands of k-points per rung before the count
  ceiling bit").  The bail never fires for si_cmce, so the
  section's si_cmce example must be revised as part of the fix.

  **Smearing was tried and ruled out (2026-07-20).**  Thermal
  smearing is wired end to end (engine
  `THERMAL_SMEARING_SIGMA` / makeinput `-thermsmear` + rc
  `therm_smear_main` / producer `gaussian-<width>` token), so
  turning it on is a value choice, not construction.  Three live
  si_cmce runs settle it: `0.0` (off) -> runaway, 21 rungs;
  `0.026` (room-temperature kT) -> still runaway, 16+ rungs, the
  oscillation only slightly damped; `0.1` (~4x room T) -> still
  not converging, its per-step energy changes NEARLY IDENTICAL to
  the `0.026` run (same `+7.4x` reversal at `4-6-7`).  Quadrupling
  the width barely moved the wobble, so a wider Gaussian smear is
  not the lever.  Smearing stays available via the token, but the
  default stays OFF, which keeps every insulator exact.

  **Resolution -- classify on the gap, then settle.**  The gap is
  already computed and stored (`gap_ev` in every `result.toml`)
  and separates the cases cleanly: si_cmce reads ~`0.000` at
  nearly every mesh, while si_ia-3 (`0.24`-`0.92`) and diamond Si
  (`1.1`-`1.9`) stay well above.  It is a DIRECT metal signal,
  better than the indirect rising-energy proxy the bail uses.  And
  the deliverable is a ROUGH starting potential, not a converged
  energy: a metal's energy never settles cleanly with mesh, but
  its potential at a modest mesh is a fine starting guess -- far
  better than the isolated-atom potential we use now.  So chasing
  convergence on a metal is wasted effort.  The rule:
  - As the climb proceeds, check the gap at each mesh.  The FIRST
    mesh whose gap is essentially zero (below a small threshold)
    declares the system a metal.  This is a LIVE trigger, not a
    pre-scan: a near-metal may show a small non-zero gap at coarse
    meshes and collapse to zero at a finer one.  The collapse to
    zero always eventually returns (or is already there) for a
    true metal, and that returning zero is the trigger -- so a
    material that starts at, say, `0.15` and later reads `0`
    settles at that later mesh.
  - On the trigger, STOP the convergence chase and settle at the
    current mesh (or one step further for a slightly better
    sampling), recording the result as a metal / rough potential.
    This guarantees cheap termination -- no runaway, no smearing,
    no stall heuristic.  si_cmce's gap is ~0 already at its floor
    `[2,4,4]`, so it would settle within a rung or two.
  - Insulators are unchanged: they never read gap ~0 at floor-and-
    above meshes, so they never enter the metal path and converge
    by the existing flatness rule.

  **What this simplifies away.**  The gap classifier can RETIRE
  the indirect rising-energy near-metal bail, and it makes the
  stall / two-sided-oscillation ideas we weighed earlier
  unnecessary: branch once on the gap and treat metals and
  insulators each simply.  (C124 already removed the coarse-mesh
  gate.)

  **Two small decisions to settle in DESIGN.**  (i) the gap
  threshold for "essentially zero" -- ~`0.05` eV cleanly separates
  a real metal from a small-gap insulator (si_ia-3's converged
  ~`0.25` eV sits safely above); (ii) settle at the triggering
  mesh, or `+1` step.

  **Chain / next step.**  Rewrite DESIGN 3.12.3 around gap-based
  metal detection -- replacing the falsified bail narrative -- then
  PSEUDOCODE 4e.3, then code.  Evidence preserved under
  `share/curation/workspace/wingbeats/si_cmce_64_1999/`.

- [x] C126. The reduce shell code counts atoms IN THE CELL, so it is
  not the transferable descriptor DESIGN says it is.  Found
  2026-07-21 by the C109 full-vs-prim comparison, which the defect
  crashed outright.  RESOLVED 2026-07-21 by option A (make it
  transferable).  DESIGN + CODE; DESIGN 5.2 (the `shell_code`
  record) / 5.6.5 (the consumer's match) / 5.10 (the family split).

  **What was measured.**  The same material, the same reduce recipe,
  two cells of it -- reading the neighbour geometry
  `ReduceMatcher.compute_query` actually sees (`min_dist` after
  `set_limit_dist(5.0)`), for diamond silicon `si_fd-3m_227_2001`:

      conventional (8 atoms):  2.33 x4,  3.80 x3
      primitive    (2 atoms):  2.33 x1

  `structure.min_dist` is an `N x N` MINIMUM-IMAGE matrix over the
  atoms of the cell.  `set_limit_dist` sizes the periodic search so
  each PAIR gets its shortest periodic distance, but the matrix
  still has one row per atom in the cell: no periodic neighbour list
  is ever built.  So the shell walk enumerates CELL CONTENTS, and
  the descriptor it produces is a property of the cell as much as of
  the environment.

  **Why that is a defect and not a documented limitation.**  DESIGN
  5.2 introduces the `shell_code` as the cross-structure descriptor
  and justifies its element-only neighbour multiset precisely so
  that it "would transfer to the query structures this stored
  fingerprint is later matched against (5.6.5)".  Transferability is
  the stated purpose.  A descriptor that reports four neighbours or
  one for the same atom, depending on how the curator drew the cell,
  does not have it.

  **It is already wrong in the cell we ship.**  Diamond silicon has
  4 nearest neighbours and 12 second neighbours.  The conventional
  cell reports 4 and 3: the twelve second neighbours are periodic
  images of only three distinct cell atoms, so the second shell is
  truncated by a factor of four.  The first shell is right only
  because an 8-atom cell happens to contain all four of them.  Every
  reduce record harvested so far is therefore a cell artifact that
  agrees with the physics by luck of cell size, not a coordination
  environment.

  **The bispectrum is NOT affected, and this was verified rather
  than assumed.**  The loen descriptors of the two cells are
  bit-identical across all nine channels
  (`0.4249608E+01 ... 0.3259187E+01`), because the Fortran engine
  builds a real neighbour list under periodic boundary conditions
  out to the sub_spec cutoff.  So the dedup key (5.2.3 keys on the
  preferred bispectrum) and the transferable descriptor every
  harvested entry carries are sound; the damage is confined to the
  reduce family.

  **Interim state (done, commit with this entry).**  The exhaustion
  case now REFUSES with a message naming the level, the cell's atom
  count, and the two ways out, instead of crashing with
  `'>=' not supported between instances of 'int' and 'NoneType'`
  (the exhausted search left the atom index at 0, and
  `min_dist[atom][0]` is the 1-indexed padding slot).  Note the
  crash was never prim-specific: ANY cell with `num_atoms <= level`
  hit it, including a 2-atom conventional cell.  Refusing rather
  than emitting an empty shell is deliberate -- an empty shell is a
  value the walk did not find, and inventing one would put a
  descriptor in the database that no structure produced.

  **The work.**  Decide what the reduce descriptor is FOR, then make
  it that.  If it is meant to be transferable (DESIGN's claim), the
  shell walk must run over a periodic neighbour list built to the
  sub_spec cutoff -- the same geometry the bispectrum already uses
  -- rather than over cell atoms, and every stored reduce record is
  then invalidated and must be re-harvested.  If it is meant to stay
  a cheap within-one-structure grouping key (which is what
  `group_reduce` uses it for, and where cell-dependence is harmless
  because the comparison never leaves one structure), then DESIGN
  5.2 must stop calling it transferable and 5.6.5 must stop matching
  it across structures.  Those are different descriptors and the
  chain currently claims both.

  **Chain note.**  The shell walk was specified at NO level -- neither
  DESIGN nor PSEUDOCODE described the closest-atom / thick-band
  algorithm; it was ported from the historical `group_reduce` and
  PSEUDOCODE 11.3 only delegated to `run_reduce_in_python`.  That
  absence is why a cell-dependence this basic went unreviewed.

  **RESOLUTION -- option A, make it transferable.**  Chosen after
  measuring that the change does NOT disturb grouping, which was the
  only reason to hesitate.  The `compare_walks.py` study ran both the
  old cell-atom walk and a periodic neighbour-list walk over a
  1296-atom amorphous silica model and four small crystals: on the
  glass the two agree atom-for-atom on the shells themselves, and on
  every structure -- including an 8-atom diamond cell where the shells
  differ completely (second-level counts of 3 against 12) -- the
  species PARTITION is identical.  Grouping compares atoms WITHIN one
  structure, and the old cap truncated every symmetry-equivalent atom
  the same way, so the relative comparison survives even where the
  absolute counts were wrong.  So option A corrects the descriptor
  without changing what grouping does with it, and A's walk is
  strictly more correct for grouping too -- no reason left to keep two
  descriptors (option C) or merely re-document the defect (option B).

  Written down the chain: new DESIGN 5.11 (what the descriptor IS --
  a periodic neighbour list, why that makes it transferable, the
  walk, exhaustion) with 5.11.1 recording why it replaced the
  cell-atom walk and that grouping was measured unchanged; DESIGN 5.2
  updated where it introduces `shell_code`; PSEUDOCODE 11.3's
  `compute_query` now spells out the `shellCode` walk instead of
  delegating to a named-but-unwritten helper -- closing the chain gap
  above.  Code: `ReduceMatcher.compute_query` walks the extended-cell
  images via a new `_neighbor_list` helper (a neighbour is a periodic
  IMAGE, counted once per image, the central atom's own images
  included; only distance zero is excluded); `ReduceStructureView`
  swaps `min_dist` for the `direct_xyz` / `num_atoms_ext` /
  `ext_direct_xyz_list` / `ext_to_central_item_map` the walk reads,
  all already populated by `create_min_dist_matrix`; both hand-built
  views in makeinput (`group_reduce`, `_file_reduce_query`) updated.
  No new geometry is computed -- the extended arrays already existed,
  which is what made A a day's work rather than a rewrite.

  The exhaustion guard from the earlier crash-only commit survives in
  new form: a cutoff too short to seed every level is refused (a
  cutoff problem now, not a cell-size one, since a neighbour list is
  not bounded by the atom count).  Verified on the real full-vs-prim
  diamond cells: both now report identical shells (4 nearest + 12
  second), where before they gave 4+3 and 1.  A real-structure
  regression test reads the shipped diamond fixture through the
  production matcher and asserts the physical 4-and-12 multiplicities
  -- the test that would have caught this originally.  1021
  non-integration tests pass.

  **Consequence for C109.**  The blocker is gone: adopting `cell =
  "prim"` no longer ships reduce records that disagree with the
  conventional ones, because they no longer disagree.  C109's default
  decision can now be made on cost alone.

  **Left for later (small):** the ~1.3x slower per-atom walk (it scans
  the extended list rather than a row of a matrix) is invisible at
  seed scale but could be cached per structure if a large campaign
  ever makes it matter.

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
- [x] C64. Add the ASE-free StructureControl factory to
  structure_control.py (ARCH 9.3, D12): build a
  StructureControl from (lattice, fractional coords,
  element symbols).  No ASE import; shared by C65 and
  C67.
- [ ] C65. Implement ase_imago.py ImagoCalculator (ARCH
  9.3, D12): calls the C63 API; the Atoms-reading glue
  uses the C64 factory.  Stays a flat module.
- [x] C66. Implement cod_fish.py (ARCH 9.5, D14): the
  four-verb COD front-end -- get (strict pinned urllib
  fetch, revision-verified), search (result.php query by
  composition/author, exact-composition element-count
  bounds, numbered table saved to a cwd session file),
  pin (index -> resolved revision -> manifest fragment),
  rank (advisory triage).  build_initial_potentials
  imports its get for the canonical fetch (DRY with
  _fetch_cod_structure).
- [x] C67. Implement cif2skl.py (ARCH 9.5, D14): ASE CIF
  read -> authored asymmetric unit + IT#/setting; resolve
  the spaceDB variant by apply_space_group verification
  against ASE's expansion; write asymmetric unit + token;
  hard error + --space override on no match.  No spglib.
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
  options study (multi-node, workers/node, launcher) per the
  user goal -- off critical path, later cluster session.
  Whether to claim whole nodes is outside this study: C119
  settles it in DESIGN 6.2.11.
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
  **Follow-on code (a)-(f) DONE + tested (2026-06-13, full suite 640
  passed).**  Includes the two-dictionary `submodel` refinement
  (DESIGN 6.2.8/6.2.10).  makeinput grew `-converg`;
  `make_producer_options` now emits dest-keyed coded tool settings
  only (functional->xccode, kpoint_integration->scfkpint via a
  smearing-aware mapper, basis->scf_basis, scf_threshold->converg,
  shift->kpshift) while the builder takes the human sub-model in its
  own `submodel` dict; `imago.OPTION_KEYS` +
  `wingbeats.CACHE_ONLY_KEYS` exported; `ImagoWingbeat.run`
  partitions options across makeinput/imago and drops the cache-only
  build identity (retiring the shared-options `run_structure` call);
  `_KEY_SCALAR_NAMES` keys on `converg`; `pick_converged_unit` reads
  `status.toml` first and skips failed/result-less units.  **What
  remains for C74: the live cluster smoke run** to validate the seam
  end to end (the `jobs/c74_si_test/` manifest, re-run with `--force`
  per the reminder above).  C54 and the C60 light half (the reduce
  fingerprint harvest) are done; only the C55/C58 bispectrum half of
  the harvest remains stubbed.
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
- [ ] C92. Make historical guidance reachable from an INDIVIDUAL
  Imago run, not only through a flight.  Today the predictor
  (`guidance_db.predict`) is consumed exclusively by the
  flight-builder helper (`build_kpoint_convergence`, C71), so a
  user who wants to run one structure must hand-write a flight to
  borrow the cluster's accumulated guidance.  Provide a direct
  path: EITHER a `-guidance` switch on `makeinput.py` (load the
  dataspace, build the Signature query from the structure +
  settings, call `predict`, and fill the k-point density from the
  prediction unless the user pinned one explicitly) OR a thin
  standalone tool (`predict_settings.py`) that takes a structure +
  the calculation settings and prints the recommended settings for
  the user to apply by hand.  Decision pending on switch-vs-tool
  (the makeinput switch is the lower-friction default; the standalone
  tool keeps makeinput free of a dataspace dependency).  Must reuse
  the same predict-then-fallback contract as the flight path
  (under-trained -> wide-grid default, DESIGN 7.9) so individual and
  flight runs never disagree.  Needs a seeded dataspace (C75) to do
  anything useful, but the wiring is independent of the seed.

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
  handling.  Register in src/scripts/CMakeLists.txt.  CODE;
  DESIGN 8.5/8.6; PSEUDOCODE 16; ARCH 11.
  PSEUDOCODE 16 now exists, and splits this item in two.  The
  library half -- constants and registries (16.1), fingerprint
  and registry loader (16.2), the twelve-rule reader (16.3), the
  deterministic emitter and its round-trip (16.4) -- is fully
  specified and can be built now.  The predictor half (16.5) is
  BLOCKED on C115: `fit_group` cannot be written faithfully
  until DESIGN settles how censored observations enter the fit,
  and `feature_row`'s correction terms are a seed-calibrated
  form.  Build the library first; do not guess the predictor.
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
- [x] C100. Wire the producer (and other clients) onto SLURM
  dispatch.  Today `curation_executor` returns a `LocalExecutor`,
  so flights run locally and one-at-a-time even on a cluster login
  node (ARCH 9.7).  Build the dispatch-config story: a per-site
  resource-control file (queues, account, per-node cores and
  accelerators, worker bring-up), a few per-run choices (topology,
  partition, nodes, walltime), and a generator that assembles a
  Parsl `Config` for either topology -- one shared pooled
  allocation, or one scheduler job per unit -- with local kept as
  the default.  Change every client (not just the producer): set
  `flight.parsl_config` from the generator and let `dispatch`
  auto-select Local vs Parsl -- no client builds its own
  executor.  The re-run/cache-bypass switch moves to a
  `dispatch` argument (DESIGN 6.2.5) -- DONE in `dispatch.py`.
  The command-line default is `slurm-per-job` (a missing
  settings file is a config error, not a quiet local
  fall-back); `local` is the explicit opt-out the tests and
  laptops request, and the library entry point defaults to
  `local`.  Ship `cluster_probe.py` (a separate tool; the
  `clusterrc.py` settings file stays pure data) that reads the
  discoverable tiers off the scheduler/node (`sinfo`,
  `scontrol`, `lscpu`) and writes a starter `clusterrc.py`
  with the required core left blank.  The design points
  (site-config home, per-run UX, deferral of per-unit
  right-sizing, generator home, per-job default, discovery
  tool) are settled in DESIGN 6.2.11; ARCH 9.8 RESOLVED.  C81
  layers predictive sizing on top of this.
  CODE; VISION Goals 4/6/7, ARCH 9.4/9.7/9.8, DESIGN 6.2.11.
- [ ] C81. Provisioning consumer in the flight layer (the
  kaleidoscope flight-builder helper or a thin sibling): query
  the predictor with a proposed config + size, apply a safety
  margin, annotate the Parsl provider's SLURM resource request;
  cold-start fallback for empty fingerprints (DESIGN 8.6, 8.8).
  The near-term consumer (VISION Goal 6).  Builds on C100.
- [ ] C82. Seed the resource dataspace on the local cluster: a
  small manual + flight-harvested set of observations spanning
  the size range on each available hardware fingerprint, so the
  predictor becomes non-trivial.  Needs cluster time.  Sibling
  of the C75 convergence seed.  A seeding run of shipped tools,
  not new machinery.  CODE; DESIGN 8.
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
  CODE (Fortran + wingbeat); DESIGN 7.2; ARCH 11.

- [ ] C115. Settle the DESIGN 8.9 open questions, plus two gaps
  the PSEUDOCODE 16 pass surfaced.  Writing section 8's
  pseudocode (2026-07-09) showed that the library half of the
  dataspace is fully specifiable today and the predictor half is
  not.  Blocks C77's predictor and C81's provisioning consumer;
  does NOT block C77's library half.  Work:
  (a) **How censored observations enter the fit.**  DESIGN 8.7
  says an OOM memory figure is a lower bound and a timeout
  walltime an upper bound, never a point; 8.9 leaves open how
  such bounds reach the least squares.  PSEUDOCODE 16.5's
  `fit_group` currently EXCLUDES them and says so -- a
  placeholder, since discarding an OOM wastes real evidence
  while feeding it in as a point biases every exponent.  Decide:
  a censored (Tobit-style) regression, or a bound-weighting
  scheme.
  (b) **The correction-term form.**  16.5's `feature_row` gives
  the v1 terms DESIGN 8.6 names by quantity, single-sourced for
  re-derivation from the seed (C82).
  (c) **Aggregate vs per-rank memory** (8.9), which decides what
  the parallel correction even means.
  (d) **Build effects on numerics** (8.9): whether the build
  block is ever referenced from the section-7 convergence side,
  against the no-cross-reference boundary.
  (e) **NEW: no objective auto-promote rule for a cost
  observation.**  DESIGN 8.7 borrows 7.8's curation discipline,
  but 7.8's `auto_promote_ok` rests on a flatness test with no
  cost analogue.  Either design a criterion or drop
  `auto-promote` from the promoter's modes and have the curator
  review every observation.
  (f) **NEW: `SAFETY_MARGIN` is a safety parameter with no
  value.**  DESIGN 8.6 rightly defers it to calibration, but the
  provisioning consumer must not ship with it unset -- so the
  seed (C82) is a hard prerequisite for C81, not merely a source
  of accuracy.  Say so in DESIGN 8.8.
  DESIGN 8.6/8.7/8.8/8.9; PSEUDOCODE 16.9.

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
