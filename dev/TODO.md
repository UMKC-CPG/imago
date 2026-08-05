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
- [x] V4. Add Principle 15 (attribution is a non-negotiable,
  not a formality) to VISION.  Done 2026-07-27 alongside the
  licensing chain landing (A10): the architecture had grown a
  convention governing every shipped source file while the top
  of the chain recorded no principle motivating it.  States the
  two obligations that follow -- cite the work Imago derives
  from, and build so that Imago's own credit survives
  redistribution -- and fixes the copyright-versus-authorship
  distinction the file header convention rests on
  (VISION Principles.15).

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
- [x] A10. Licensing and attribution machinery (ARCHITECTURE
  1.1).  Done 2026-07-27.  Imago shipped with no license at all,
  which made it technically all-rights-reserved and left nothing
  for a downstream reader to attribute to.  Settled on ECL-2.0,
  continuing OLCAO's choice: verbatim LICENSE, NOTICE carrying
  the OLCAO lineage and the Birkbeck College space group
  provenance, CITATION.cff with ORCID, CONTRIBUTING.md stating
  the copyright-versus-authorship policy, README.md, and a
  two-line SPDX + copyright header on all 259 shipped source
  files (256 under src/ plus rappture/ and middleware/).  The
  tree is REUSE 3.3 compliant, verified by `reuse lint` at
  313/313 files, with REUSE.toml covering data and documentation
  that cannot hold a comment.  Governed at ARCHITECTURE rather
  than PSEUDOCODE because file headers carry no algorithmic
  content; see the "Why this level" note in 1.1.
- [ ] A11. Decide the licensing status of src/data/spaceDB.tgz
  (ARCHITECTURE 1.1).  The REUSE.toml blanket rule declares
  src/data/** as ECL-2.0, which sweeps in spaceDB.tgz -- whose
  contents derive from the Birkbeck College crystallographic
  tables credited in NOTICE.  The compilation and encoding are
  ours; the underlying tables are not, so the blanket
  declaration asserts a license over material that did not
  originate here.  Crystallographic tables are largely
  uncopyrightable facts, so this may well be fine, but it is a
  claim rather than a formality and should not rest on a glob
  pattern.  Resolve by either confirming the blanket rule or
  splitting spaceDB.tgz into its own annotation.  A NOTE in
  REUSE.toml flags the question at the point of decision.
- [ ] A12. Mint a DOI and complete CITATION.cff (ARCHITECTURE
  1.1).  Connect the repository to Zenodo, tag a release, and
  backfill the three fields left commented in CITATION.cff:
  version, date-released, and doi.  The DOI is what a reference
  list can actually point at -- a repository URL is not a
  citable object and does not accrue credit through the
  indexing services -- so this is the step that converts the
  citation metadata from a statement of intent into something
  a reader can cite.  Feeds D20, which needs the DOI for the
  banner text.

- [x] A13. Establish the runtime output control facility
  (ARCHITECTURE 12).  Written 2026-07-27 alongside D20, which
  is its first client.  IMAGO_VERBOSENESS takes a
  comma-separated list of category NAMES; the bitmask is
  private.  A numeric form was considered and deliberately
  rejected: accepting one publishes the bit assignment as soon
  as anyone writes it into a job script, and this project
  records every invocation into a `command` file, so those
  settings persist -- renumbering afterward silently changes
  what old scripts do.  Names can be reordered freely and an
  unrecognized one can be reported; a numeric form can be added
  later without breaking anything, but cannot be withdrawn.
  Three behaviours fixed: unset means `normal` not silent,
  unknown names warn and continue, `none` is explicit.  Only
  the `banner` category is defined; the categories that matter
  for the debugging and parallelization campaigns are
  deliberately left unenumerated so that work names them.
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
- [ ] D21. The convergence threshold sits BELOW the ladder's own
  noise floor, so passing it is luck rather than convergence.
  This one sentence accounts for the whole cluster of symptoms
  chased on 2026-07-31: si_cmce's Gaussian ladder "converging" on
  scatter, three separate one-rung-short near-misses, and four of
  five elemental metals running to the ceiling.  A bar beneath the
  noise can only be cleared by a coin landing the same way twice.
  **The measurement.**  Top-of-ladder rung-to-rung scatter, all in
  eV per atom: Al 0.00082 under LAT and 0.00469 under Gaussian,
  Fe ~0.0015, the transition metals ~0.001.  The threshold is
  5e-4.
  **What relaxing it does.**  Converged mesh by threshold, mesh
  order, failed rungs excluded throughout:

        solid              5e-4        1e-3        2e-3        5e-3
        Al                 --          [16,16,16]  [16,16,16]  [16,16,16]
        Co                 [16,16,16]  [16,16,16]  [16,16,16]  [12,12,12]
        Cu                 --          [12,12,12]  [10,10,10]  [8,8,8]
        Fe                 --          --          [13,13,13]  [12,12,12]
        Ni                 --          --          [16,16,16]  [12,12,12]
        si_cmce            [14,14,13]  [11,11,10]  [8,8,7]     [8,8,7]
        si_fd-3m (six)     [12,12,12]  [10,10,10]  [10,10,10]  [8,8,8]
        si_ia-3            [11,11,11]  [9,9,9]     [7,7,7]     [6,6,6]

  **The si_ia-3 row is CORRECTED, 2026-08-04.**  It previously read
  `--  --  [11,11,11]  [9,9,9]`, i.e. failing at the two tightest
  bars.  That is wrong in a specific and instructive way: the whole
  row was shifted one column right, so every entry was attributed
  to the next looser threshold.  D22's numbers for the same solid
  disagreed with it and turned out to be the correct ones.
  Re-measured from scratch on a clean 17-rung ladder (meshes 4 to
  20, one binary, one integration scheme, every rung converged,
  computed in a single run) under `jobs/si_ia3_remeasure/`.  The
  new energies match the old ones digit for digit, so the shift was
  in the SCORING, not the physics.
  Consequence: si_ia-3 converges at every bar, including 5e-4.  The
  claim that tightening the threshold would cost us this solid --
  which nearly blocked C146 -- was never true.

  At 2e-3, ALL THIRTEEN converge.  The insulators degrade
  gracefully -- [12,12,12] to [10,10,10] -- while si_cmce, a
  metal, pays most ([14,14,13] to [8,8,7]), which is precisely the
  case DESIGN 3.12.3 already concedes should be rough.
  The Al row was first recorded as failing at every threshold
  below 5e-3.  That row had been read off the GAUSSIAN ladder --
  the al_parity comparison ran its two arms into the same run
  directory and Gaussian went second, so one row of an otherwise
  all-LAT table was a different scheme.  For the record, the two
  differ sharply and the difference is the point:

        Al LAT       --   [16,16,16]  [16,16,16]  [16,16,16]
        Al GAUSSIAN  --   --          --          [18,18,18]

  Under LAT Al converges at 1e-3; under unsmeared Gaussian it
  needs 5e-3.  "Al resists relaxation" was an artifact of the
  mislabelled row, not a property of aluminium.
  **Why this suits the deliverable.**  3.12.3 justifies the metal
  short-circuit by observing that "the initial-potential database
  wants a *rough* good starting point for a later self-consistent
  calculation, not a converged energy".  The threshold never got
  the same scrutiny.  5e-4 eV/atom is 0.5 meV/atom -- a
  publication-grade bar applied to a starting guess that the next
  calculation re-converges anyway, and to guidance a curator reads
  as advice rather than as an answer.
  **Two properties that make this cheap and honest.**  It needs no
  code: `kpoint_convergence_threshold` is already a per-solid
  manifest key (`HARVEST_SETTING_KEYS`, curation_manifest.py), so
  a curator wanting a tighter bar on one material can still ask
  for it.  And it is self-describing: a staged entry carries
  `metric_threshold` in `entry.verification`, so a looser entry is
  labelled as such and the dataspace does not degrade silently.
  **Proposed:** raise the default to sit above the measured noise
  floor -- 2e-3 on this evidence -- and say in DESIGN that the
  threshold MUST sit above it, with these numbers as the reason.
  A later reader should not tighten it back without re-measuring
  the floor.
  **Confirmed live, 2026-07-31.**  Both campaigns re-run with
  `kpoint_convergence_threshold = 2e-3` in the top-level
  `[harvest]` block.  Every rung hit the cache -- the threshold is
  a harvest setting and is not part of the cache key -- so this
  cost minutes and tested the rule as the climb GROWS a ladder
  rather than as a re-score of a finished one.  Al [16,16,16], Cu
  [10,10,10] and all six si_fd-3m [10,10,10] matched the table
  exactly, and staged entries carry
  `metric_threshold = 2.0e-03` as intended.
  Two rows did not match, both for understood reasons.  si_cmce
  settled at [5,5,5] because the Si manifest carries no
  `climb_shape`, so the automatic climb's metal short-circuit
  fired -- a different question, not a different answer.  si_ia-3
  gave [7,7,7] against a predicted [11,11,11] because it computed
  a NEW rung at [3,3,3]: `stride_threshold` is
  `threshold * stride_flatness_multiple`, so relaxing the
  threshold also loosens the bracket and the climb takes a
  different PATH, not merely a different verdict on fixed rungs.
  Offline re-scoring is therefore exact for `unit_step` solids
  (all five metals matched) and only indicative for bracket-refine
  ones.  Worth remembering before trusting any future re-score.
  **What this does NOT resolve.**  It makes the ladder's
  non-monotonicity stop mattering; it does not remove it.  That
  property is real, measured, and recorded below, and this entry
  should not be closed as though the ladder were fixed.

  ---

  **The recorded property: the mesh ladder is not a monotone
  refinement.**  DESIGN 3.12.1 separates the roles cleanly --
  density is the currency, the mesh is the step -- and 3.12.2
  justifies the step as "each is the *finest* mesh the symmetry
  permits".  Both hold.  Neither establishes what the two-sided
  flatness test assumes: that a later rung is a BETTER sampling
  than the one before it.  For fcc at the standard shift it is
  not.  Irreducible point counts by parity:

        odd  n:  19  44  85 146 231 344 489 670
        even n:   8  16  29  47  72 104 145 195 256

  At comparable n the odd meshes carry 2-3x the sampling, so
  [20,20,20] (256 points) is a step BACKWARDS from [19,19,19]
  (670).  The ladder alternates between a good sampling and a poor
  one and the energy zigzags accordingly.
  How much of the metals' non-convergence this explains is
  PARTIAL, and the distinction matters: Cu converges once the
  rungs are ordered by sampling, so for Cu the ladder was the
  whole problem.  Al and Fe converge under NEITHER ordering, so
  something further is at work for them -- Fe especially, which is
  bcc and has no parity split to begin with (counts 6, 10, 14 ...
  250, already monotone) yet still scatters at ~1.5e-3.  Do not
  close this expecting the ladder fix to convert every metal.
  **Not a LAT defect -- measured, not assumed.**  fcc Al, same
  manifest, one token changed, run directories/prepare/scratch
  cleared before each arm, `atomicPDB` verified to hold only
  "isolated" both times so neither arm warm-started from the
  other's potential, both arms opening on [4,4,4] and agreeing on
  the irreducible count at every rung.  Top-four rung-to-rung
  scatter (eV/atom):

        family      LAT       gaussian
        full      0.00082     0.00469
        odd-only  0.00120     0.00127
        even-only 0.00943     0.01803

  BOTH schemes zigzag and Gaussian zigzags worse, so the
  alternation belongs to the k-mesh sequence.  `generateTetrahedra`
  and the `fullKPToIBZKPMap` unfolding are exonerated; no code bug
  is in play.  A by-product worth keeping: on a parity-consistent
  ladder LAT reaches 1e-4 flatness at [17,17,17] while unsmeared
  Gaussian is still scattering at 1.3e-3 -- LAT converges Al and
  Gaussian does not.  Ladders kept at
  `jobs/al_parity/ladders/{linear-tetrahedral,gaussian}/`.
  **Secondary options, kept because they were measured, not
  because they are the fix.**  Raising the threshold above the
  noise floor buys twelve of thirteen; ordering buys two of five.
  Reach for these only if the threshold change proves
  insufficient, or in combination where a solid still resists.
  (a) Walk a sampling-consistent ladder.  Two candidates: fix a
  parity and step cubic n by 2, or order rungs by IRREDUCIBLE
  COUNT.  Re-scoring the existing metal ladders both ways favours
  the second, and not marginally:

        solid  mesh order      irreducible-count order
        Al     NOT converged   NOT converged
        Cu     NOT converged   converges at n=16
        Ni     converges n=18  converges at n=13
        Co     converges n=16  converges at n=13
        Fe     NOT converged   NOT converged

  Three of five against two, and Ni and Co settle at n=13 rather
  than 18 and 16 -- cheaper meshes, not merely more passes.
  Ordering also generalises: Fe is bcc and its counts are already
  monotone (6, 10, 14, ... 250) with no parity split at all, and a
  low-symmetry cell increments one axis class at a time so
  "parity" is not even defined for it.  Fixing a parity would also
  commit every rung to the EXPENSIVE family, since for fcc the odd
  meshes are the well-sampled ones.  Read the parity split as the
  symptom that exposed the problem in fcc, not as the thing to
  fix.
  (b) Extend by one rung when the ceiling is reached with the
  PENULTIMATE rung two-sided flat.  Endpoints are never eligible
  (3.12.3), so a ladder can stop exactly one rung short of being
  able to confirm what it has already found.  Observed three
  times and re-verified against `pick_converged_climb`: si_cmce
  uncorrected LAT at `max_count = 18` (penultimate [17,17,15]), Al
  re-sorted by irreducible count, and Al's odd-only family here.
  The trigger is the CONJUNCTION -- ceiling reached, verdict not
  converged, penultimate rung two-sided flat.  Penultimate
  flatness alone is not the signal; converged ladders show it too,
  which is unsurprising since they settle near the top.  Cheap to
  detect, and it converts near-misses into results.
  (d) `gap_ev` is measured on a mesh chosen to converge the
  ENERGY, which mis-classifies metals and degrades insulator
  entries.  Split out as **D22** -- it concerns what a guidance
  entry measures, not how the climb stops.  It belongs here only
  in that raising the threshold makes it worse, by moving the rung
  the gap is read from.
  (c) Exclude rungs whose SCF did not converge from the flatness
  test.  NOT secondary -- this one is a plain defect and survives
  every threshold.  Today such rungs are read as ordinary energies
  and can be SELECTED: ni_fm-3m_225_2006 reported
  `converged_mesh = [18,18,18]`, a rung whose own `result.toml`
  says `status = "not_converged"`.  Cu and Fe carry failed rungs
  inside their ladders too.  Do this regardless of what happens to
  the threshold.
  DESIGN 3.12.1 / 3.12.2 / 3.12.3.

- [ ] D22. A gap read at the energy-converged mesh is not a
  converged gap.  `gap_ev` is recorded from whichever rung the
  climb happens to stop on -- a mesh chosen to converge the
  ENERGY -- and it is then used as a predictor key (DESIGN 7.6)
  and as the metal test's sole input (3.12.3).  Neither use was
  ever argued for.  Split out of D21(d) on 2026-07-31.
  **Failure 1: metals recorded as gapped insulators.**  The 2e-3
  confirmation staged these:

        al_fm-3m_225       converged_mesh=[16,16,16]  gap_ev=0.124
        cu_fm-3m_225_2011  converged_mesh=[10,10,10]  gap_ev=0.185

  both with `gap_kind = "indirect"`, against
  `metal_gap_threshold = 0.05`.  Aluminium is the textbook
  free-electron metal.  The 0.124 eV is the finite-mesh artifact
  of C138(e) -- a level SPACING in the globally sorted spectrum,
  not a gap at E_F -- and at the rung this climb stopped on it
  cleared the metal cutoff by 2.5x.  Ni, Co and Fe were gated
  correctly on the same run, so this is not a blanket failure but
  a silent, material-dependent one, which is worse.  Both entries
  were deleted from staging before promotion could fix them in
  place; the underlying defect is untouched.
  **Failure 2: the recorded gap moves with the convergence bar.**
  si_ia-3, same material and same scheme, differing only in the
  bar (numbers confirmed 2026-08-04 on a clean re-measured ladder):

        thr=5e-4   mesh=[11,11,11]   gap_ev=0.193
        thr=1e-3   mesh=[9,9,9]      gap_ev=0.257
        thr=2e-3   mesh=[7,7,7]      gap_ev=0.370
        thr=5e-3   mesh=[6,6,6]      gap_ev=0.429

  **This entry's reading of that was WRONG, and the truth is
  worse.**  It said "si_ia-3 really is an insulator" and called the
  move a degraded gap.  The premise is false.  Scored over 17 rungs
  from [4,4,4] to [20,20,20], si_ia-3's gap goes to ZERO as 1/n^2:

        n     gap      gap*n^2
        10    0.1682   16.82
        14    0.0900   17.64
        18    0.0564   18.27
        20    0.0464   18.56

  `gap*n^2` is flat -- it rises by a factor of 1.10 while n
  DOUBLES.  A real gap approaches a constant, so that product would
  have grown by a factor of 4 over the same span.  This is not a
  gap converging; it is the level spacing at the Fermi surface
  shrinking with the mesh, the artifact `populate.F90:253-265`
  describes in advance.  Si III in the BC8 structure is reported as
  semimetallic, which fits.
  So the 0.193-to-0.370 move is not a real gap degrading.  It is
  the SAME artifact read at two different meshes, and every value
  in the table above is fictitious.
  At [20,20,20] the gap reads 0.0464 eV -- BELOW the 0.05 eV metal
  cut.  Push the ladder a little further and this solid classifies
  as a metal.
  **Failure 2 as originally stated is REFUTED, and the real defect
  is a different one.**  Measured 2026-08-04 on diamond silicon,
  the unambiguous insulator of the seed set: 21 rungs from [4,4,4]
  to [24,24,24], one build, one scheme, every rung converged
  (`jobs/si_fd3m_gapcheck/`).  Its gap IS real -- `gap*n^2` grew
  30.45x from n=4 to n=24 against 36x predicted for a constant gap,
  as against si_ia-3's 1.10x -- and it settles near 0.805 eV.
  The threshold costs that gap almost nothing:

        thr=5e-4   mesh=[12,12,12]   gap_ev=0.8046
        thr=1e-3   mesh=[10,10,10]   gap_ev=0.8061
        thr=2e-3   mesh=[10,10,10]   gap_ev=0.8061
        thr=5e-3   mesh=[8,8,8]      gap_ev=0.8285

  A TENFOLD change in the bar moves the recorded gap by 3% -- 24
  meV.  The 0.18 eV move this entry was built on came entirely from
  si_ia-3, where the gap was fictitious at every mesh.  On a
  material with a real gap, the convergence bar is not the problem.
  **What IS the problem: mesh parity.**  Consecutive rungs of the
  same ladder disagree far more than any two thresholds do:

        [11,11,11]  0.9572
        [12,12,12]  0.8046
        [13,13,13]  0.9145

  19% between neighbours.  The ladder carries a strong parity
  sawtooth -- even meshes settle to ~0.805 by n=12 while odd meshes
  are still 4.6% high at n=23 (0.8419) -- so the two families
  approach the same limit at very different rates.  The recorded
  gap therefore depends on the PARITY of whichever rung the climb
  stops at.  This run landed on 12, 10, 10 and 8, all even, and
  agreed with itself by luck; a bar settling on [11,11,11] would
  have recorded 0.957 against 0.806.
  So the headline claim survives intact -- a gap read at the
  energy-converged mesh is not a converged gap -- but the mechanism
  is mesh-family sensitivity of order 19% between neighbours, not
  threshold sensitivity of order 2x.  That is a smaller error and a
  much more erratic one: it does not scale with anything a curator
  controls, and it cannot be reasoned about from the threshold.
  The energy on that same ladder is settled to the eighth decimal,
  -7.77652760 repeated across the top nine rungs.
  Caveat on provenance: the si_fd-3m ladder was built by commit
  9e5a936b and the si_ia-3 one by a5fc5bc9, a POPTC commit having
  landed between them.  Each ladder is internally consistent, which
  is what the scaling test needs, and the change was to the optical
  path rather than the ground-state SCF -- but the two solids were
  not measured by one binary, and that is worth saying rather than
  assuming away.
  **And this is a live miss by the metal machinery.**  At the
  meshes the climb actually stops on -- 7, 9, 11 -- the gap reads
  0.37, 0.26, 0.19 and looks solidly insulating.  Nothing C142,
  C143 or C144 added catches it, because all of those read gaps at
  meshes the ENERGY chose.  The classification's reliability
  therefore depends on the energy threshold, which was never meant
  to carry that weight: a tighter bar climbs higher and is more
  likely to see the collapse.  This is the same defect as Failure
  1, reached from the other end.
  Note the contrast with the energy on that same ladder, which is
  genuinely settled: scatter across meshes 14-20 is 9.8e-5 eV/atom,
  five times below even the old 5e-4 bar.  The energy converges and
  the gap does not, on one ladder, which is this item in one line.
  **Why this is not just a D21 side-effect.**  The case for a
  looser threshold (D21) was argued about the ENERGY, and holds
  there: the deliverable is a rough starting potential and a rough
  guide, and the SCF that consumes the potential re-converges it.
  Nothing re-converges `gap_ev`.  It is stored, predicted from,
  and compared across materials.  So the trade D21 makes is energy
  roughness we accept in exchange for gap roughness nobody agreed
  to -- and the defect predates D21, which only widened it.
  **Candidates** -- three possible fixes, in increasing cost:
  (1) **Record which mesh the gap came from**, in the entry, so a
      consumer can discount it.  Bookkeeping only: it does not
      improve the number, it makes the number inspectable.
  (2) **Require the gap to persist as the mesh densifies** -- the
      same two-sided, flat-over-consecutive-rungs test the ENERGY
      already gets, applied to the gap.  A gap that has not gone
      flat is not recorded as though it had.
  (3) **Converge the gap on its own terms**, climbing until the gap
      settles regardless of where the energy settled.  Most
      correct, most expensive: it can demand meshes the energy
      never needed.
  **DECIDED 2026-08-04: (1) with the flatness measurement recorded.
  (3) is out** -- not as a controller for initial-potential
  generation, nor for guidance, for now.
  **(1) as literally worded was already done.**  The entry has
  carried `verification.converged_mesh` since C118, and `gap_ev` is
  read off that same rung, so the mesh the gap came from was always
  recorded.  What was missing is any statement of how far that gap
  could be trusted.  So the work became: measure (2)'s test, record
  the number, act on nothing.  Landed as C147.
  **A number, not a boolean.**  Storing a verdict would freeze a
  tolerance nobody has chosen, and entries written under different
  tolerances would then disagree silently -- exactly what
  `metric_threshold` exists to prevent for the energy.  The stored
  `gap_spread` is the raw relative movement, so a consumer picks
  its own bar and every entry stays comparable.  It also means the
  later decision about unsettled gaps can be taken from data
  gathered during the rebuild rather than guessed at now.
  **The 2026-08-04 measurements point at (2), and (1) is worth
  doing anyway.**  The defect turned out to be a parity sawtooth
  rather than a slow drift, and a sawtooth is precisely what a
  two-sided flatness test rejects: at [12,12,12] the gap sits 19%
  below both its neighbours, so it cannot read flat, while a
  genuinely settled gap passes.  The test would also have caught
  si_ia-3, whose gap falls monotonically and never goes flat at
  all -- so one rule covers both failures.  (3) is not obviously
  needed: on diamond silicon the even-mesh family settles by n=12,
  which is where the energy converges anyway, so the cost may be
  ladder-shape rather than extra rungs.  (1) is cheap, independent,
  and would have made both of these findings visible from the
  stored entries instead of requiring a re-measure.
  The metal test needs the same treatment either way -- reading
  one rung's gap is what let a free-electron metal past it.
  Relates to C138(e) (what `gap_ev` measures, and why a coarse
  mesh reports a spacing as a gap) and to D21 (which rung the
  climb stops on).
  DESIGN 3.12.3 / 7.6; guidance schema in DESIGN 7.2.

- [x] D23. A declared cache key file that a unit's job never stages
  makes that unit permanently uncacheable.  Found 2026-08-05 while
  checking whether the producer works at its own defaults
  (`jobs/defaults_live_check/`).  Two consecutive campaigns under a
  byte-identical manifest each reported the fingerprint calculation
  as `run (no usable result)` while all ten ladder rungs reused, so
  the miss is deterministic rather than a stale-workspace accident.
  **The mechanism.**  `standard_key_fields` declares the same two
  key files for every unit
  (`kaleidoscope/builders/kpoint_convergence.py:276`):

        _KEY_FILE_NAMES = ("structure.dat", "kp-scf.dat")

  `cache_key_matches` requires each of them to exist at the run
  directory's ROOT and to byte-equal the freshly prepared copy, and
  returns a miss when the staged side is absent (`cache.py:90-93`).
  The fingerprint unit is a `loen` job built with `scf_basis = "no"`
  (`build_initial_potentials.py:1067-1069`): it runs no SCF, so the
  flattening step puts `kp-pscf.dat` at its run-directory root and
  never `kp-scf.dat`.  The staged side is therefore absent forever
  and the unit can never be reused.  Read off disk, the fingerprint
  unit beside a mesh rung of the same solid:

        unit             structure.dat       kp-scf.dat
        loen-bispectrum  present, identical  ABSENT from run dir
        kpt-mesh-10      present, identical  present, identical

  **Introduced 2026-07-30 by C135** (`f7bfc96`), which added
  `kp-scf.dat` as a second key file so that a change of integration
  scheme could not silently return another scheme's answer.  Its own
  comment states the assumption that broke: "A key FILE costs
  nothing, because every run directory already stages this one."
  True of the SCF rungs it was reasoned about; false of the
  fingerprint job, which was not in view.
  **What it costs.**  Every fingerprint calculation in every
  campaign since 2026-07-30 has been recomputed from scratch.  On
  the one-atom silicon cell that is 0.28 s, which is why it went
  unseen, but the cost scales with cell size and with the
  fingerprint cutoff.  Second, for that unit the discrimination C135
  bought does not operate at all, because the file it compares is
  not there to compare.  It fails SAFE -- a permanent miss rather
  than a wrong hit -- but that is not the designed behaviour and no
  document says it is what should happen.
  **Candidates** -- three, in increasing damage:
  (1) **Compare the `inputs/` copies rather than the flattened root
      copies.**  Both sides keep an `inputs/` carrying
      `structure.dat`, `kp-scf.dat` and `kp-pscf.dat` for every
      unit -- verified on the fingerprint unit itself -- so the
      comparison becomes symmetric and every declared key file
      exists for every unit.  This also dissolves the prepare step's
      asymmetry, which `build_initial_potentials.py:1896-1918`
      spends a paragraph explaining.
  **DECIDED: (1).  Its preconditions were CHECKED 2026-08-05, all
  three clear.**
  *Nothing reclaims a run directory's `inputs/`.*  The whole path
  holds two deletion sites.  `tidy_scratch.py:703` removes
  `scratch_target(run_dir, scratch_root)`, which must lie under
  `$IMAGO_TEMP` -- the tree the `intermediate` symlink points at,
  never the run directory; both `--clean-after` and `--tidy-run`
  route here and both skip when `$IMAGO_TEMP` is unset.
  `wingbeats.py:265` removes the flattened ROOT copies on every
  commit, by design, so `imago.py` stays their only writer.  That
  second site corroborates the defect rather than threatening the
  fix: the root copies are cleared on each commit and refilled only
  for names the unit's job reads.
  *Every wingbeat stages an `inputs/`.*  `ImagoWingbeat` is the only
  one registered (`wingbeats.py:361`); `Wingbeat` is a
  documentation-only base class.  The suite's fakes build no
  `inputs/`, but the generic cache tests pass names with no
  directory component and are indifferent -- which confirms the
  shape: `KeyFile.name` should carry a path RELATIVE to the run
  directory, so `cache.py` only joins it and never learns about
  `inputs/` or imports `makeinput`.
  *No existing cache is invalidated.*  Over every run directory in
  the workspace, comparing each key file's root copy against its
  `inputs/` copy: 411 identical, 0 differing, 0 missing an
  `inputs/` copy, and 13 missing a root copy -- exactly the
  fingerprint units, ONE PER SOLID, all thirteen solids.  The
  defect is universal, not a silicon quirk.  Forward-checking the
  proposed rule for the one solid with prepare directories on disk:
  22 of 22 units would hit, the fingerprint among them.  The
  load-bearing row is `kpt-mesh-10-10-10`, whose prepare copy was
  built under HEAD `76eff98` and whose run copy was staged a day
  earlier under `9e5a936b` -- different campaigns, different builds,
  byte-identical, which is the cross-build reuse the fix must
  preserve.  Weaker for the fingerprint unit itself, whose two
  copies are same-campaign: that shows the structural blocker is
  removed, not that its inputs are reproducible across campaigns.
  Only the live re-run closes that.
  **Test impact, bounded:** `test_kpoint_convergence.py:184-186`
  asserts the literal name pair and
  `test_build_initial_potentials.py:4232-4301` stages copies at the
  root, so both change; the generic cache tests
  (`test_kaleidoscope.py:396-526`) and the C139 root-copy tests
  (`1253-1361`) are untouched.
  (2) **Declare key files per job type**, so the fingerprint unit
      keys on `kp-pscf.dat`.  Correct, but it restores the
      hand-maintained list of "options that matter" that
      `standard_key_fields` was written to avoid.
  (3) **Treat a declared-but-unstaged key file as not applicable**
      rather than as a miss.  Cheapest and worst: it silently
      weakens the key in exactly the way C135 existed to prevent.
  Whichever is chosen, DESIGN 6.2.5 must first say what a declared
  key file means for a unit whose job does not produce it.  The
  contract is silent on that today, which is how the two halves of
  the mechanism came to disagree without either one looking wrong on
  its own.
  Relates to C135 (which added the second key file) and to C134
  (which made the cache key physics-only, so a surviving workspace
  is now expected to be reused across builds and this miss is felt
  on every campaign rather than only after a rebuild).
  **DONE 2026-08-05, candidate (1), down the chain in order:**
  DESIGN 6.2.5 (`095517e`), then PSEUDOCODE 13.1 / 13.4 / 11.4
  Phase 1b (`1d277b8`), then the code (`b2babdc`), then the tests
  (`512baef`).  `KeyFile.name` became `KeyFile.path` and holds a
  path relative to the unit's directory; the producer declares
  `inputs/structure.dat` and `inputs/kp-scf.dat`; `prepare_units`
  collapsed to the same join the hit-test makes.
  **The agreement test was added rather than the property being
  handed to C139.**  Moving identity off the run-directory root
  would have transferred a guarantee: today the cache compares the
  very file the engine reads, so the key necessarily describes what
  runs.  So a root copy that EXISTS must byte-equal its `inputs/`
  copy, which separates the two cases one test could not -- absent
  means this job does not read that file, present-and-disagreeing
  means the engine would run inputs the key does not describe.
  C139 still prevents the second at the source; this is the
  backstop, and 6.2.5 says neither may be dropped for the other.
  **A second, pre-existing defect fell out of the tests.**  The
  agreement test passed, then failed, then passed.  `filecmp.cmp`
  memoizes on a (mode, size, mtime) signature, which is blind to
  precisely what this cache must catch: a same-size rewrite in
  place.  Two writes microseconds apart routinely share one mtime
  tick -- 169 of 200 measured here -- and in EVERY one of those the
  comparison reported differing files as equal.  That is a false
  HIT, a stored result returned for inputs that changed, and it is
  the failure with no escape valve, since nobody knows to reach for
  `--force`.  It predates this work: `filecmp` was already how
  bytes were compared, and the agreement test is simply the first
  test that compares one pair twice with a change in between.
  PSEUDOCODE 13.4 already specified an unconditional read, so the
  constraint was written there and `cache_key_matches` now clears
  the memo first.  A test forces the condition with `os.utime`
  rather than leaving it to luck.
  **Confirmed live 2026-08-05**, campaign `15725038` after
  `cmake --install`: 11 of 11 units reused, 0 run, the fingerprint
  line reading `reuse` for the first time since C135 landed.  The
  ten rungs still reuse from `9e5a936b` while the fingerprint
  reuses from `76eff98`, so the physics-only key still works across
  builds and no existing hit was lost.  The staged entry diffs
  IDENTICAL against the pre-fix run apart from `entry_id` and
  `generated_at` -- same mesh, same `gap_ev`, same `gap_spread` --
  so the cache change moved no science.
  Each campaign stages a fresh entry file for the same solid at the
  same settings (three now, byte-identical in content), which is
  the REGENERATE model working as designed: promotion resolves
  duplicates into one entry plus superseded ones.  Noted only
  because promotion has been seen to handle two copies and not yet
  three.
  DESIGN 6.2.5; PSEUDOCODE 13.1 / 13.4 / 11.4.

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
- [x] D20. Design the runtime citation banner (VISION
  Principles.15, ARCHITECTURE 1.1 and 12).  Written 2026-07-27
  as DESIGN 10.  Imago prints a "if you use these results,
  please cite" block, the way LAMMPS, VASP, and Quantum
  ESPRESSO do -- the highest-leverage attribution mechanism
  available, because unlike a license header or CITATION.cff it
  reaches the person at the moment they have results in hand
  and are writing the paper.  Every open question is now
  settled.  SPLIT IN TWO (10.2): an identity block at startup
  (art, version, Imago's own citation) and a methods block at
  the end, because at startup the program does not yet know
  which methods a run will exercise.  Art is LITERAL TEXT, not
  generated (10.3) -- it was hand-kerned from the project logo
  and no script reproduces that; 51 columns to match the
  character(len=51) opLabels.  Art AND citation text live in
  src/data/banner.txt (10.3, 10.4), installed to share and
  found through IMAGO_DATA exactly as elementData.f90 and
  potential.f90 already find elements.dat -- no new mechanism
  and no new failure mode, since a broken IMAGO_DATA kills
  elementData before the banner is ever reached.  An earlier
  draft argued for compiling the citation in; that rested on a
  false premise about the cost of a data file and was corrected
  once the existing mechanism was checked.  Destination (10.6) is
  the log and nothing else: not the tabular DOS/bond/optc
  outputs, which are parsed positionally, and not the three
  HDF5 files either -- an attribute was considered and rejected,
  because citation guidance belongs in the human-readable
  output and cluttering data files to restate it buys nothing a
  reader sees.  Methods block (10.5) is a
  registry pairing each DESIGN References entry with a
  predicate reading state the engine already holds
  (kPointIntgCode and the like).  Suppression via
  ARCHITECTURE 12 rather than a flag (10.7): the engine's
  arguments are positional, so a flag would change the argv
  contract in Fortran and Python together.
---

## PSEUDOCODE

- [x] P11. Write pseudocode for the runtime output control
  facility and the citation banner (ARCHITECTURE 12, DESIGN
  10).  Written 2026-07-28 as PSEUDOCODE 17.  THREE modules,
  not the two ARCHITECTURE 12.4 planned: one module holding
  both the identity block and the methods registry closes the
  Fortran cycle O_CommandLine -> O_Banner -> O_KPoints ->
  O_CommandLine, because kpoints.f90:1093 already uses
  O_CommandLine.  The split falls on the seam DESIGN 10.2 drew
  -- a block that prints before the work begins can depend on
  nothing the work produces -- so O_Banner keeps the identity
  block and depends only on O_Verboseness, while a new
  O_MethodCitations holds the registry and may read engine
  state freely.  ARCHITECTURE 12.4 corrected to match; that is
  the legitimate upward edit, since the constraint is factual
  and was only visible at this level.  Also settled here: the
  category table pairs name with bit implicitly (bit = row
  index - 1), so no second column can drift; `normal` and
  `none` are reserved set-valued aliases, not categories;
  tokens combine by union and nothing subtracts; matching is
  case-insensitive; every failure path warns to unit 20 and
  continues, unlike elementData's fatal IMAGO_DATA check.  The
  read and write formats must be '(a)' and never list-directed,
  or every line of art shifts one column right of the len=51
  timestamp rules; and the read buffer is 132, not 51, because
  the citation lines run to 64 today and would truncate.  Two
  registry entries, not three: UFF appears nowhere in
  src/imago/, so the Rappe reference belongs to the Python
  force-field path and an entry for it would be dead on
  arrival.  Feeds C133.
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
- [x] C4a. Add Bloechl's curvature correction to the corner
  integration weights (DESIGN 1.3.1; PSEUDOCODE 3a).
  DONE 2026-07-31.
  `bloechlCornerCorrection` in `mathSubs.f90` returns the four
  `dw_i`; `computeElectronPopulation_LAT` adds them to the
  weights at the accumulation step.  One insertion point covers
  both consumers, because `populateLAT` and the
  bond/effective-charge path both fill the same array through
  that routine.
  **Scope, corrected twice.**  The entry originally named
  `computeTDOS_LAT`, the one consumer the correction does NOT
  reach: the DOS uses `cornerDOSWt_LAT`, the energy derivative
  dw/dE, and eq. 22 corrects `w`.  It was then widened here to
  "all three consumers", which was also wrong for the same
  reason.  It reaches the two that evaluate cumulative
  integration weights.  A DOS correction is a separate question
  nothing has answered.
  **Measured, not asserted.**  Three si_cmce ladders, same
  manifest, one variable at a time (README beside the ladders):
  Gaussian converges at [15,15,13] / 448 k-points; uncorrected
  LAT does not converge and hits the ceiling at 720; corrected
  LAT converges at [14,14,13] / 392 -- a LOWER density than
  Gaussian, which is the accuracy-at-coarse-meshes claim the
  correction exists to make.
  The verdict is the weaker evidence, since all three rest on a
  flatness test over a plateau whose scatter is near the
  threshold.  The strong evidence is cross-scheme agreement:
  from ~48 k-points up, corrected LAT tracks Gaussian to within
  +/-0.0006 eV/atom, while uncorrected LAT sat +0.0084 away at
  504 and was still moving.  Two approximations to the same BZ
  integral must converge to each other; with the correction they
  do, and that test is independent of the flatness rule.
  The correction's size falls from -0.040 eV/atom at 12 k-points
  to -0.0094 at 504, as a finite-mesh correction should, and it
  resolves most of the C138(e) gapped-mesh puzzle.
  **Verified invariants.**  A throwaway harness over 7
  tetrahedron shapes (including two-fold, three-fold and total
  degeneracy) x 201 energies: the four terms sum to zero to
  5.6e-17, so the electron count and the Fermi search are
  untouched; the correction is exactly zero at all 609
  out-of-range samples, which is what leaves gapped systems
  alone; the sign follows `(epsBar - eps_i)` everywhere.  The
  harness is not in the build -- see C138(d), there is still no
  Fortran test suite to put it in.
  **The citation, settled.**  This entry and DESIGN's reference
  list both said "eqs. 22-24", which reads as three equations
  owed and cost a round of hedging in the code and the design.
  Checked against the paper 2026-07-31: the correction is eq. 22
  alone.  Equations 23 and 24 compare the true Fermi surface
  against the interpolated polyhedral one -- an assessment the
  paper makes, not a formula Imago computes.  The bibliography
  now says so, so the range cannot mislead the next reader.
  **Left open.**  The correction on a genuine insulator, where
  it should vanish identically rather than merely become small.
  si_fd-3m at a COARSE mesh is the cheap check, and it doubles
  as the diagnostic C138(e) still wants.

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

- [x] C108. Intermediate-scratch cleanup for the producer (and a
  reusable cleanup subsystem).  Motivation: the Si seed run
  (2026-07-02) left 3.7 GB of per-calc scratch under each run
  directory's `intermediate ->` symlink (roughly 20 MB per calc
  dir), against only 18 MB of kept home-side artifacts
  (status/result/scfV/descriptor).  Scratch of this kind fills up
  fast and is tedious to locate and remove by hand -- and stale
  units from earlier manifests linger in the shared workspace (the
  seed run's workspace still held `si_diamond`, `si_fd-3m_227_2010`,
  and a half-finished `si_p63mmc_194_2018` from prior experiments).
  Three layers.  The build order below was originally written two
  ways -- (b) is labelled NEAR-TERM "before (c) lands" yet listed
  after it -- and was settled 2026-07-21 in favour of (c) first:
  measurement showed scratch is 99.7% of the bytes and that
  pruning it is provably safe (DESIGN 6.2.12), so the standalone
  tool recovers essentially the whole saving at once, while (b)
  only matters for a campaign large enough to exhaust scratch
  MID-flight, which nothing has yet reached.  Order: (c), (a),
  then (b) when that pressure arrives.
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

  **DONE.  (c) and (a) built 2026-07-21; (b) built 2026-07-23.**

  Measured first, which reshaped the task.  Scratch is **99.7%** of
  the bytes a run leaves behind -- 3.17 GB of 3.2 GB on a seed-scale
  run, essentially all the HDF5 carrying the wavefunctions, about 25
  MB per calculation -- against 222 KB of kept files per calc.  So
  reclaiming scratch is nearly all of the available saving, and the
  motivation's "roughly 20 MB per calc dir" was about right.

  **Why it is safe, established rather than assumed.**  Every path
  in a `result.toml` `outputs` table resolves INSIDE the run
  directory; none points through `intermediate`.  And `is_cache_hit`
  decides from `status.toml` plus `cache_key.toml` alone.  So a
  reclaimed run still harvests and still counts as a cache hit --
  verified on the live workspace after reclaiming it: the potential
  still extracted (16 coefficients), every output still present, the
  cache key and status intact.

  Chain: ARCHITECTURE 9.6 gains the two-tier kept/scratch model and
  names `tidy_scratch.py`; DESIGN 6.2.12 defines the subsystem --
  the mechanism/policy split (mirroring the cache's split in 6.2.5,
  and for the same reason: only the client knows when a finished run
  is finished WITH), the refusals, dry-run-by-default, and all
  three layers; PSEUDOCODE 13.8 gives the walk.  (The tool began as
  `tidy_workspace.py` with four refusals; the hand-run block below
  records its growth to two roots and seven refusals.)

  Code: `src/scripts/tidy_scratch.py` (layer (c)) -- selective by
  id, calc glob, and age, previewing by default and removing nothing
  without `--apply`.  Producer `--clean-after` (layer (a)) calls
  that same planner rather than reimplementing the walk, so the two
  cannot drift about what is safe to delete.  Reclamation runs after
  the databases and run log are written, and swallows failures: a
  stuck directory must not turn a successful build into a failed
  one.  28 tests; 1049 non-integration tests pass.  Proven live:
  1.6 GB reclaimed from the C109 full workspace, 6 KB left, every
  kept file intact.

  **(b) PRUNE-AS-YOU-GO BUILT 2026-07-23.**  There is no fixed
  headroom to name as its trigger: `$IMAGO_TEMP` has no per-user
  quota, so the ceiling is whatever a SHARED filesystem happens to
  have free at the time -- unguaranteed, and reduced by everyone
  else's jobs.  (An earlier note here claimed 37 TB free and a
  ~40,000-calculation threshold; that was the whole shared
  filesystem, not an allowance, and the figure derived from it is
  withdrawn.)  A campaign cannot know in advance how much room it
  will get, which argues for pruning in flight rather than
  against it.

  **The seam, decided 2026-07-23.**  Layer (b) adds NOTHING to
  kaleidoscope.  A flight already fires `on_outcome` once per unit,
  in landing order, carrying that unit's run directory -- exactly
  the moment and the fact a prune needs -- so (b) is a producer
  wiring: `make_prune_callback` builds the callback, and both the
  loen pre-flight and the climb flight carry it.  The alternative,
  a reclaim-policy field on the flight that the dispatcher acts on,
  was rejected: reclamation reads imago's own names, and inside
  kaleidoscope the one place engine knowledge belongs is the
  wingbeat.  The dispatch core beneath it names no imago file at
  all, and this must not be the change that teaches it one.

  That also corrected a claim the earlier documents carried:
  DESIGN 6.2.12 and ARCHITECTURE 9.6 both said the reclamation
  MECHANISM was kaleidoscope's.  It never was -- it lives in
  `tidy_scratch.py`, a sibling tool -- and both now say so and say
  why.

  **A prune that fails is contained but never hidden.**  The
  callback runs inside `collect`, which does NOT guard the hook, so
  an escaping exception would abandon the whole campaign over
  housekeeping.  Nothing may escape.  But a REFUSAL (unfinished,
  nested, too recent) is the mechanism working and is reported
  quietly, while a FAILURE (a removal attempted and declined, or
  the mechanism raising) means a permission/mount/lock assumption
  has broken -- so it is printed when it happens AND carried to an
  end-of-run summary, because a line an hour deep in a log is lost.
  Deliberately not fatal: the databases and run log are intact.

  Code: `plan_one_dir` / `reclaim_one_dir` in `tidy_scratch.py`
  (the per-directory decision the whole-tree sweep now also calls,
  so in-flight and after-the-fact cannot diverge); the nesting
  refusal takes its comparison set as an argument, since it is the
  one refusal a single directory cannot judge alone.  Producer
  `--tidy-run`, `make_prune_callback`, `report_prune_problems`.
  17 new tests; 1093 non-integration tests pass.  NOT yet exercised
  on a live campaign.

  **HAND-RUN GAP CLOSED 2026-07-22 (the actual daily cost).**  The
  tool only ever saw kaleidoscope workspaces, but an ordinary
  `imago.py` run plants the same `intermediate` link -- and those
  are the common case, not the exception.  On this machine they
  held ~1.0 GB the tool could not reach: no `wingbeats/` (so the
  CLI refused the root) and no `status.toml` (so the default
  policy would have refused each run anyway).

  Renamed `tidy_workspace.py` -> **`tidy_scratch.py`**, since a
  workspace is now one of two roots it recognizes.  The root's
  kind is DETECTED, and one call handles exactly one kind, so a
  single report never gathers two safety contracts under one set
  of totals.  Two refusals were added for the job tree (a run must
  hold no `imagoLock` AND end its `runtime` log with the
  completion marker; a workspace nested in a job tree is named but
  never descended into).

  Two findings worth keeping.  First, `runtime` is opened in
  APPEND mode, so a directory run four times holds four completion
  markers -- only the TAIL is truthful, and `c/diamond/full2` is
  the live proof: four markers, a log ending mid-run, and a stale
  lock.  Searching the file would have deleted a killed run's
  scratch.  Second, the marker is written from a `finally`, so it
  means the driver reached cleanup, NOT that the run succeeded; a
  job tree has no `result.toml` with which to preserve failures
  the way the workspace contract does.

  **A safety hole the live run exposed, now refusal 5.**  Scratch
  mirrors the run directory's path, so a run nested inside another
  (`knbo3/cubic/debug` inside `knbo3/cubic`) has its scratch
  nested too.  Removing the outer tree would have taken the inner
  one with it -- double-counting the saving, and deleting the
  working files of a run the other refusals had just declined to
  touch, which would make "never touch an unfinished run" a
  formality.  An outer tree holding another run's scratch is now
  deferred and named; a second pass takes it once the inner ones
  are gone.  Containment is tested against every run the walk
  found, not only the selected ones, since a filtered-out inner
  run is exactly the collateral case.  On the live tree this
  corrected the reported saving from 895.9 MB to 787.3 MB.

  Chain: DESIGN 6.2.12 (two roots, refusals 5-7), PSEUDOCODE 13.8
  (`detect_root_kind`, `find_job_run_dirs`, `hand_run_policy`,
  resolve-then-judge planner), ARCHITECTURE 9.6.  1075 tests pass.

  Layer (b)'s chain: DESIGN 6.2.12 (mechanism placement corrected,
  how (b) hooks, contained-not-hidden failures), PSEUDOCODE 13.8
  (`plan_one_dir`, `reclaim_one_dir`) and 11.4 (`make_prune_callback`
  plus the producer's `clean_after` / `tidy_run` wiring, which had
  never been specified for layer (a) either), ARCHITECTURE 9.6.

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
- [x] C128. Retire condense.py's hardcoded SLURM template in
  favour of the md job class.  The template in
  `create_lammps_files` asks for 125 tasks on a partition
  whose nodes hold 48, hardcodes the account, and predates
  the site-config machinery entirely.  Five edits: (a) add
  the `md` block to `clusterrc.py` -- ranks/walltime/memory
  plus `init`, the latter shipped blank and flagged REQUIRED;
  (b) add the same block to `cluster_probe.py`'s own copy of
  the schema and `md.init` to its blanked keys -- the drift
  test compares whole dictionaries, so a block added to one
  file and not the other fails it; (c) add `build_md_sbatch`
  to `kaleidoscope/cluster_config.py`, beside the
  orchestrator generator; (d) call it from
  `create_lammps_files`; (e) tests.
  Three behaviours to get right, each settled in DESIGN
  6.2.11.  An unfilled required core makes `load_site_config`
  raise and condense.py inherits that refusal rather than
  writing a file built from guesses.  An empty `md.init` is
  refused by the generator and NOT by `_require_core`, so a
  site that flies calculations and never condenses is not
  refused a flight over a setting no flight reads.  And the
  generator writes `OMP_NUM_THREADS=1` after the bring-up,
  because ranks sized to fill a node must each hold one core.
  Note (c) also changes `build_orchestrator_sbatch` to open
  with `#!/bin/bash -l`: the login-shell rule belongs to the
  act of generating a submission file, not to either class.
  Validated by C37/C42, which already call for an end-to-end
  condense.py run.
  CODE; ARCH 9.4, DESIGN 6.2.11, PSEUDOCODE 10c/13.7.
  DONE.  All five edits landed, and the refine pass that
  followed them changed two things it found.  The site file is
  now read at *settings time*, in `ScriptSettings.__init__`
  beside the `condenserc.py` read, and handed down to
  `condense_write_submission(site)`: refusing an unconfigured
  site at the moment the file is written would have discarded a
  whole run's bonds, angles and LAMMPS input first, and the
  loader searches the current directory, which
  `create_lammps_files` has by then left for `lammps/`.  And
  the queue-override typo guard now descends one level, into a
  block, matching the merge -- a pre-existing gap that let
  `{"md": {"rank": 999}}` merge to a stray key while the job
  ran at the site's width.  Still to be validated live by
  C37/C42.
- [x] C129. Move the producer's structure cache out of the
  basis database.  It sat at `share/atomicBDB/cache/structures/`
  because an earlier draft put it beside a per-solid SCF cache
  at `.../cache/scf/`; DESIGN 5.7 later dropped that cache in
  favour of kaleidoscope's run-reuse cache, and the survivor
  left downloaded CIFs inside the atomic *basis* database,
  where they are neither per-element nor basis data.  The
  layout level had never named the location at all -- ARCH 8.1
  described only the atomicPDB tree -- which is how it drifted
  there unnoticed.  New home `share/curation/structures/`,
  beside the flight workspace and run log, so
  `structure_cache_dir` and `curation_workspace_root` now
  derive from one root and a campaign's whole footprint clears
  in one gesture.  The split is along reconstructibility: a
  deleted structure is re-fetched by the next run, a harvested
  potential entry cost cluster time and cannot be.
  ARCHITECTURE 8.1 (new layout block + the rationale),
  DESIGN 5.7 (names the location; the dropped-SCF-cache
  paragraph no longer cites the dead path), PSEUDOCODE 11.4
  (`structure_cache_dir` now specified, not just referenced).
  CODE.  DONE.  Old cache cleared on disk first -- 74 files,
  including two stray `makeinput -cif` run directories left
  inside it and skeletons under the pre-`<cell>` naming rule.
  Nothing migrated: the cache is rebuilt by the next run.
- [x] C130. Dedup re-runs at promotion.  The guidance dataspace
  had no dedup at any stage, and nothing else could supply one:
  the harvest writes one entry per converged solid with no view
  of the promoted corpus, and a re-run mints a fresh `entry_id`
  by construction (the slug hashes flight id, structure, and
  timestamp), so no collision check ever fired.  Found while
  clearing staging before the seed re-run -- 67 staged files
  describing 8 solids, one of them (`si_cmce`) recorded at two
  wildly different densities either side of the C125 metal fix.
  With `neighbor_count = 5` the copies would have filled the
  whole neighbour set, collapsed the variance to zero, and
  reported a single measurement as near-certainty, which drives
  the climb into its narrowest search and past the crystalline
  opening floor.
  Claim key is `(system_type, basis, functional,
  kpoint_integration, basename(source_structure))` -- settings
  in, because they are the sub-model partition; basename not
  path, because the cache moved in C129; `imago_commit` out,
  because it is what the comparison examines.  Agreement is the
  `converged_mesh`, an exact integer test with no tolerance
  knob.  Three outcomes: redundant retires the staged copy to a
  new `superseded/` area (promotion only ever ADDS to
  `entries/`); a differing or uncomparable mesh is a CONFLICT
  that is reported and left in staging, never resolved
  automatically; anything else takes the ordinary per-mode path.
  Every mode applies it, `--all` included -- it waives the
  quality rule, not the correctness guard.  Batch duplicates
  resolve among themselves first, which retires the standing
  claim that staging is not a uniqueness namespace.
  ARCHITECTURE 10.1 (`superseded/`) + 10.5, DESIGN 7.8,
  PSEUDOCODE 15.7.  CODE.  DONE.  +9 tests, 1123 pass.
  Not yet exercised on a live corpus: the seed re-run is the
  first campaign whose staging this will judge.
  **Superseded in design by C134 (2026-07-29), before ever running
  live.**  The mesh-agreement test, the three outcomes, and the
  batch pre-pass described above are all gone: the check is now a
  plain existence test with two outcomes, the mesh is printed
  rather than compared, and a live in-memory index makes the
  within-batch case fall out of the ordinary rule.  Read C134 for
  what is actually being built; this entry is kept as the record of
  what was tried and why it was cut.
- [x] C131. Report outcomes and problems, not progress.  The
  producer narrated every reference solid it fetched and every
  scratch tree it pruned, so a clean eight-solid pre-flight
  filled the screen with lines saying that nothing was wrong,
  and a 79-unit campaign added one prune line per unit.  That is
  where a real failure goes to hide, which is VISION 10 read at
  the level of one run: a failure surfaced among two hundred
  lines of success has not been surfaced.
  The rule: a run says what it achieved and what went wrong, and
  does not narrate what it is doing.  Silent on a clean
  pre-flight; every fetch failure named with its reason whether
  or not anyone asked; the pre-flight tally printed only when
  something failed (it then says how much of the set survived)
  or on request.  Prune *failures* stay loud, prune *refusals*
  are the mechanism working and became narration.  `--verbose`
  restores the old behaviour verbatim, which is what makes it
  safe to have taken away.
  Verbosity is one module-level setting established by `main`
  before any work, NOT an argument: it describes how the process
  talks to its user, not how a flight converges, and threading
  it would have put a reporting concern inside four functions
  about physics and dispatch.
  DESIGN 5.7 (Reporting + the `--verbose` flag), PSEUDOCODE
  (`set_verbosity` / `narrate` / `print_materialize_report`, and
  the prune callback's narration branches).  CODE.  DONE.
  +6 tests, 1128 pass.
- [x] C132. The run-reuse cache could not survive a second run.
  Found live: the seed was auto-promoted and re-run, and the
  campaign died eight seconds in, before dispatching a unit,
  with `FileNotFoundError` out of `filecmp` inside the
  loen pre-flight's hit-test.  Two defects, and NEITHER was code
  drifting from a spec -- both were holes in the spec itself,
  which is why the chain had not caught them.
  (a) `prepare_units` pointed the `structure.dat` KeyFile source
  at `<staging>/structure.dat`, but makeinput writes
  `<staging>/inputs/structure.dat`.  The two sides of the
  byte-compare are NOT symmetric: a run directory carries the
  file both at its root and under `inputs/`, because imago reads
  it flattened at the root when the unit runs, while a prepare
  directory is never run and so is never flattened.  PSEUDOCODE
  referenced `prepare_units` twice and never specified it -- the
  function had no governing section at all -- so the spec was
  WRITTEN, then the code fixed to match.
  (b) `cache_key_matches` guarded the staged side but not the
  source, and the code implemented that faithfully, so the SPEC
  was wrong.  `filecmp.cmp` stats its arguments and raises rather
  than returning False, so one absent file aborted a whole
  campaign -- Principle 10 inverted, a per-unit doubt failing the
  flight.  The rule is now stated: a key file unreadable on
  EITHER side is a miss, never an error.
  Latent since the cache was written, because the compare is
  reached only once a prior `cache_key.toml` exists and its
  scalars match.  Every seed run until now began from a wiped
  workspace, so no run had ever executed it.
  DESIGN 6.2.5, PSEUDOCODE 13.4 (`cache_key_matches`) + 11.4
  (new `prepare_units`).  CODE.  DONE.  +4 tests, 1132 pass;
  each verified to FAIL against the restored defects.
- [x] C133. Code the output control facility and the citation
  banner (PSEUDOCODE 17).  DONE 2026-07-28.  Three new Fortran
  modules -- verboseness.f90, banner.f90, methodCitations.f90 --
  in both engine source lists; initVerboseness then
  printIdentityBlock in parseCommandLine; printMethodsBlock at
  the end of Imago; dead `banner` variable removed from
  timeStamps.f90.  Verified on a live silicon SCF: the art
  lands byte for byte identical to banner.txt, so it is flush
  with the len=51 timestamp rules; all seven IMAGO_VERBOSENESS
  cases behave (unset, none, banner, BANNER, banner+typo,
  none,banner, and ",,banner,"), with the typo warning and the
  run completing; and both registry entries print when their
  predicates hold.
  THE FIRST RUN DESTROYED ITS OWN LOG, and the spec was what
  was wrong.  cleanUpSCF, cleanUpPSCF, and loen each closed
  unit 20 before the end of Imago, so printMethodsBlock wrote
  to a closed unit -- which does not fail.  Fortran reconnects
  it to fort.20 and truncates, leaving an 8-line file holding
  only the citation that had just erased 900 lines of results.
  Nothing in the chain could have caught it: every level agreed
  with the level above it and the error was in all of them.
  Only running it found it.
  The same three routines also raised the fort.2 success
  signal, which DESIGN 6.1.2 says certifies the binary ran
  without an abortive error and which imago.py:1598 treats as
  its sole success gate.  Since every 200-series job runs both
  stages in one invocation, cleanUpSCF certified success before
  the post-SCF stage had begun -- and gfortran exits 0 on both
  STOP and STOP 'message', so the return code did not catch it
  either.  A post-SCF abort was reported as success.  Both are
  one defect: a whole-run event raised from a routine that
  knows about a single stage.  The tail of the run is now
  ordered once at the outermost level -- citations, close (20),
  fort.2 -- and printMethodsBlock inquires whether the log is
  open, complaining to stdout rather than writing if it is not,
  so a future early close costs the citations instead of
  everything.  DESIGN 10.6 and PSEUDOCODE 17.7/17.8 corrected.
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
  Background: the build identity is *recorded, never compared*
  (DESIGN 6.2.5).  The producer hangs it on each unit's `record`
  mapping, the driver stamps it into `status.toml` for the reuse
  plan, and the wingbeat echoes it into `result.toml` where a
  guidance entry's provenance reads it (C134).  It reached
  `result.toml` through neither path before C134, which is why
  every entry on disk today records `"unknown"`.  Even once it
  does, it records what the *producer* believed it ran, which can
  drift from the binary actually executed.  The robust upgrade is
  for the running binary to report its own build commit (e.g.
  from the C78 `build_info.toml`, or a compiled-in version
  string), which the wingbeat then prefers over the recorded
  value -- PSEUDOCODE 13.2 already writes the echo as a
  fallback, so this lands as a substitution in one field of one
  file rather than new plumbing.  Pairs with C78 (build
  identity) and C79 (wingbeat/imago.py capture hooks).
  CODE (Fortran + wingbeat); DESIGN 6.2.2, 7.2; ARCH 11.
  Its priority rises with C134: the cache no longer compares the
  build at all, so a *recorded* identity is the only thing telling
  a curator what produced a reused result.

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

- [x] C134. Code the cache-key and guidance-uniqueness
  simplification.  Both mechanisms had grown to ask a question
  neither could answer: the cache asked whether a stored result was
  still *good*, and promotion asked whether two entries *agreed*.
  Both hold starting guesses that downstream machinery
  re-converges, so the honest question in each place is the same
  and much smaller -- does a result for this already exist? -- with
  two outcomes and no quality judgment.  Anything expensive or
  irreversible becomes a verb a person types.  Design landed
  2026-07-29 across VISION 16, ARCH 9.6 / 10.1 / 10.5, DESIGN 6.2.2
  / 6.2.4 / 6.2.5 / 6.2.10 / 7.8, PSEUDOCODE 4e.2 / 4e.7 / 11.4 /
  13.1 / 13.2 / 13.3 / 13.5 / 15.6 / 15.7.  Three code steps, in
  this order:
  (a) **The key, and the record that replaces the guard.**
  `_KEY_SCALAR_NAMES` becomes `("converg",)`; the producer stops
  putting `imago_commit` in `options` and hangs `{imago_commit:
  <sha>}` on each unit's new `record` instead; `dispatch_unit`
  stamps it into `status.toml` at launch and `write_status`
  preserves the `[record]` table across every later rewrite; the
  wingbeat echoes it into `result.toml`; `CACHE_ONLY_KEYS` and its
  partition branch are deleted.  **The two halves of that last
  point must land together** -- with the bucket gone a key reaching
  neither tool raises, so a `make_producer_options` still emitting
  `imago_commit` aborts every unit with `unknown makeinput option`,
  the very C74 failure the seam was written to fix.
  (b) **The reuse plan and preview.**  `reuse_plan` /
  `print_reuse_plan` / `dispatch(preview=True)`; counts always, the
  per-unit lines under the driver's own module-level verbosity
  switch or in a preview.  DESIGN 5.7's rule holds here: the climb
  calls `send_off` once per round, so unconditional per-unit lines
  would refill the screen C131 cleared.
  (c) **Promotion, the metal skip, and the prediction record.**
  `promote` loses the batch pre-pass and `_by_generated_at`
  entirely and keeps ONE live index shaped `key -> (path, entry)`;
  the occupied branch reports both entries and retires the
  newcomer, with REPLACE offered in interactive mode only.
  `build_entry` returns None for a gapless run via a new scalar
  `is_gapless_value` that the rung-shaped `is_gapless` also calls,
  and both harvest paths skip on None.  `PredictionRecord` gains
  the two resolved knobs the standalone harvest cannot look up,
  `kpoint_convergence_threshold` (today stamped onto the dict after
  the fact) and `metal_gap_threshold`.
  Also: the two staged metal entries `crystalline-fbd1b7` and
  `crystalline-214bb2` are si_cmce runs the new rule would never
  have produced, and want deleting by hand.
  CODE.
  **Status 2026-07-30: complete.  (a), (b), and (c) are written and
  the suite is green at 1252 passed / 23 skipped**, +23 tests over the
  1229 the code alone reached.  Seventeen existing tests had been
  written against the retired behaviour and were moved onto the new
  spec; three changed meaning rather than shape and were renamed to
  say so -- `test_a_disagreeing_mesh_takes_the_same_branch` and
  `test_a_mesh_that_cannot_be_compared_is_not_special`, plus
  `test_a_duplicate_within_one_batch_takes_the_ordinary_branch`, which
  had claimed a `generated_at` tie-break that no longer exists.
  The new behaviour itself is now pinned: the `[record]` table
  surviving every lifecycle rewrite, and a HIT leaving it describing
  the run that produced the result; a differing recorded build still
  HITTING, which is the whole point of the change; the reuse plan's
  counts-always / lines-when-asked, and a preview that spends nothing
  and reports nothing; the wingbeat's echo and its preference for an
  engine-reported build (the C84 seam); `is_gapless_value` on both
  sides of the cut and on a missing gap; `build_entry` returning None
  for a metal, both harvest paths skipping, and the producer still
  writing the potential; and in promotion the within-batch duplicate
  through the ordinary branch, REPLACE with a SECOND replace against
  one claim still working, the unattended modes never prompting, and
  `dry-run` modelling the index.
  Two gaps the tests surfaced, both fixed.  `set_verbose` had NO
  caller, so `--verbose` never reached the driver and the per-unit
  reuse lines could not print at all -- PSEUDOCODE 13.5 already
  specified the call, so `main` now makes it and the four driver
  reporting names are exported from the package.  And
  `wingbeats._stage_inputs` still documented dropping "the cache-only
  build identity", which `CACHE_ONLY_KEYS`' deletion had retired.
  **Validated live 2026-07-30.**  A seed run from scratch, a re-run at
  the same commit (87/87 hits), then an empty commit `82e178b` and a
  third run that still hit all 87 -- its only misses seven genuinely
  new meshes the now-guided climb chose.  `cache_key.toml` is written
  on a MISS only, so its mtimes are the evidence.  `imago_commit`
  reached entry provenance as a real sha for the first time, and
  `[record]` survived the full lifecycle.  Promotion then met its
  motivating case unaided: two runs had staged 14 files, one claim
  each twice, and it promoted 7 and superseded 7.  The output also
  showed the `generated_at` tie-break really is gone -- sorted
  filename picked the LATER entry for one solid and the EARLIER for
  another.

- [x] C135. The run-reuse cache cannot tell two k-point integration
  schemes apart.  The producer's key is `_KEY_SCALAR_NAMES =
  ("converg",)` plus one key file, `structure.dat`.  The integration
  scheme travels as makeinput's `scfkpint` -> `kp_intg_code` and is
  written into `kp-scf.dat` as `KPOINT_INTG_CODE`, a file the key
  never reads, and it does NOT appear in `structure.dat`.  So one
  solid at one mesh under two different schemes resolves to the same
  run directory with a matching key and a `done` status -- a HIT that
  returns the other scheme's energy.
  **This is a different class of fault from the one C134 weighed.**
  That argument rested on a false hit being recoverable: an older
  engine computing the same physics, with `--force` as the escape
  valve.  Here a false hit returns *different physics* under the name
  of the physics that was asked for, and nothing in the output says
  so.  DESIGN 6.2.5 currently justifies keying on `structure.dat`
  because it "bakes in every input that changes the result"; that
  claim is false for the integration scheme, so it must be CORRECTED
  rather than merely supplemented.
  Found live 2026-07-30 setting up the si_cmce tetrahedral trial: the
  trial reused ten cached `gaussian` rungs and computed only two.  It
  has stayed harmless until now only because `kpoint_integration` has
  had one value in every run to date; the metals work is what turns it
  into a real axis.
  **Latent today, and the reason says when it goes live.**  LAT is
  implemented for the post-SCF properties only -- DOS/PDOS and bond
  consume the tetrahedron weights, while `valeCharge` accumulates the
  SCF valence charge from the gaussian/histogram `electronPopulation`
  with no branch on `kPointIntgCode`.  So a ground-state SCF returns
  bit-identical total energies under either code today, and a wrong
  cache hit currently costs nothing.  It starts costing correctness
  the moment LAT reaches the SCF occupation path, which is exactly
  what the metals work needs -- so this must land WITH that change,
  not after it.
  **Fix it with a second key FILE, not a key scalar.**  Adding
  `scfkpint` to `_KEY_SCALAR_NAMES` invalidates every existing
  `cache_key.toml` at once, since `cache_key_matches` compares the
  saved `scalars` table verbatim and no stored file carries the new
  name.  Adding `kp-scf.dat` to `key_fields` costs nothing: that
  function does NOT compare the saved `files` list -- it byte-compares
  each DECLARED file against its staged copy, and every existing run
  directory already stages one -- so a same-scheme re-run still hits
  and only a scheme change misses.  It also keeps 6.2.5's stated
  preference for byte-compared files over hashing: two `kp-scf.dat`
  files diff to the one `KPOINT_INTG_CODE` line that differs.
  Wiring: the driver's prepare step re-points the new KeyFile's
  `source` at `<prepare_dir>/inputs/kp-scf.dat`, the same move it
  already makes for `structure.dat`.
  A sibling gap, deliberately NOT folded in here because it wants its
  own decision: the initial-potential database has no notion of the
  sub-model at all.  `insert_or_skip` matches on LABEL
  (`<reference_id>-<element><species>-t<type>-a<site>`, no scheme in
  it), so a second-scheme run REPLACES the stored potential rather
  than sitting beside it, and `make_imago_provenance` records no
  basis, functional, or kpoint_integration -- so nothing on disk says
  which scheme produced what is stored.  The guidance dataspace treats
  the sub-model as identity; the potential database does not know the
  concept exists.
  DESIGN 6.2.5 (the corrected claim) and PSEUDOCODE 15.6
  (`KEY_FILE_NAMES`); PSEUDOCODE 11.4 needed NO change -- it already
  re-pointed every key file, and only the code had narrowed that to
  `structure.dat`, so the code was what disagreed with the spec.
  CODE.  DONE 2026-07-30.  +6 tests, 1257 pass.  Not yet exercised
  live: the si_cmce workspace is being cleared, so the LAT trial is
  the first run whose scheme change this will distinguish.

- [x] C147. Record how settled a guidance entry's gap was
  (`verification.gap_spread`).  DONE 2026-08-04.  DESIGN 7.2 led,
  then PSEUDOCODE 15.2 / 15.4 / 15.7, then the code.  D22 option
  (1), with the flatness measurement of option (2) recorded but NOT
  acted on.
  **What was actually missing.**  Option (1) said "record which
  mesh the gap came from"; the entry has carried
  `verification.converged_mesh` all along, and `gap_ev` is read off
  that same rung, so that half was already done.  The gap in the
  record was any statement of how far the gap could be TRUSTED.
  **The measurement.**  `measure_gap_spread` returns the largest
  relative change between the converged rung's gap and the rungs
  two ladder positions either side, as a fraction of that gap.
  Two design choices, both forced by measurement rather than taste:
  *Two positions, not one.*  A ladder carries a strong parity
  sawtooth in the gap.  Diamond silicon's adjacent rungs disagree
  by 19% -- [11,11,11] reads 0.9572 eV against [12,12,12]'s
  0.8046 -- even where the gap has settled to ~1% within one parity
  family.  Comparing immediate neighbours calls every ladder
  unsettled and discriminates nothing.
  *Relative, not absolute.*  Near the top of its ladder si_ia-3's
  gap moves 0.010-0.014 eV per two rungs, SMALLER in absolute terms
  than settled silicon's mid-ladder movement -- yet si_ia-3's gap
  is collapsing to zero.  As fractions: ~20% against ~1%.
  **A number, not a verdict.**  A stored boolean would freeze a
  tolerance nobody has chosen, and entries written under different
  tolerances would disagree silently -- the failure
  `metric_threshold` exists to prevent for the energy.  The raw
  fraction lets a consumer pick its own bar and keeps entries
  comparable.  Nothing acts on it today, by decision.
  **None means NOT MEASURED, never "settled".**  Absent when the
  ladder reaches two positions on neither side, when a gap was
  never read, or when the gap is zero (a metal, whose relative
  change is undefined and which stages no entry anyway).  The field
  is omitted from the file rather than written as 0.0, because a
  stored zero would claim a perfectly settled gap -- the one thing
  an unmeasured gap must never assert.
  Optional in the schema and therefore last and defaulted in
  `Verification`, so a hand-written entry need not carry one.
  Both harvest paths feed it: the producer from
  `record_converged`'s new `grid_gaps`, the standalone sweep from
  its own collapsed ladder.  The chosen rung's position is found by
  matching the converged density in `grid_values` -- exact, since
  both callers derive the two from the same numbers -- and a miss
  measures nothing rather than the wrong rung.
  +7 tests, suite green at 1278.
  Does NOT close D22.  It makes the defect visible in the stored
  entries instead of requiring a re-measure to find; whether an
  unsettled gap should suppress an entry, or the metal test should
  read this, are both still open.
  CODE; DESIGN 7.2, PSEUDOCODE 15.2 / 15.4 / 15.7.

- [x] C146. Lower `DEFAULT_KPOINT_CONVERGENCE_THRESHOLD` from 2e-3
  to 1e-3.  DONE 2026-08-04.  DESIGN 3.12.3 led, then PSEUDOCODE
  11.4, then the code -- the order C140 and C141 took.
  **Not a reversal of C140; the population changed.**  The floor
  principle is untouched: a threshold must sit above the scatter of
  the ladders it judges.  What changed is WHICH ladders those are.
  Every scatter figure that argued for 2e-3 -- Al 0.0047 under
  Gaussian, Fe ~0.0015, the transition metals ~0.001 -- was measured
  on a METAL, and that is not incidental: a metal's energy
  oscillates as the mesh crosses the Fermi surface, which is both
  why its ladder is noisy and why it cannot converge in k-points at
  all.  Those ladders no longer reach the flatness test.  C142 and
  C143 stop a metal on the gap, on every search shape, before any
  convergence work.  So 2e-3 was a floor set by ladders the test
  never sees, and "all thirteen converge" counted five solids that
  now stop for an entirely different reason.
  **What it costs the insulators: nothing measurable.**  For the six
  ordinary si_fd-3m seeds, 1e-3 and 2e-3 pick the IDENTICAL mesh,
  [10,10,10].
  **What it buys: gap quality.**  `gap_ev` is read off whichever
  rung the climb stopped on, is a predictor key, and is the one
  recorded quantity nothing downstream re-converges -- unlike the
  potential, which the consuming SCF re-converges by construction.
  A looser bar stops earlier and records a coarser-mesh gap.  This
  does not fix D22, but it stops widening it.
  Swept the stale 5e-4 text with it: five DESIGN sites plus the
  "0.5 meV/atom is looser than the textbook 1 meV/atom" paragraph,
  which now reads the other way round -- 1 meV/atom IS the textbook
  bar, reached from below by measurement rather than adopted for
  being familiar.  One `mesh_climb` comment cited 5e-4 as "the
  default" while explaining how `stride_flatness_multiple` was
  tuned; it now names 5e-4 as the condition that experiment ran
  under, and flags the multiple for re-checking if the bar moves far.
  Suite green at 1271, no test changed -- the one test that asserts
  the default reads the CONSTANT, which is what C140 fixed it to do.
  **The one case that needed watching is CLEAR, 2026-08-04.**
  `si_ia-3` was the worry: the seed closest to metallic behaviour,
  and the two tables contradicted each other about it.  Re-measured
  from scratch on a clean 17-rung ladder
  (`jobs/si_ia3_remeasure/`), it converges at EVERY bar, 1e-3
  included, at [9,9,9] against [7,7,7] at 2e-3.  So the tightening
  costs it two rungs and loses nothing.  D22's numbers were the
  correct ones; D21's row was shifted a column and is now fixed.
  Nothing else in the seed set is at risk from the change: the six
  si_fd-3m insulators pick the same mesh at 1e-3 and 2e-3, and the
  metals no longer reach this test at all.
  The re-measure turned up something larger, carried in D22:
  si_ia-3 is not an insulator.  Its gap goes to zero as 1/n^2, so
  every gap value recorded for it is a finite-mesh artifact.  That
  does not affect this item -- the ENERGY on that ladder is
  genuinely settled, to 9.8e-5 eV/atom at the top.
  CODE; DESIGN 3.12.3, PSEUDOCODE 11.4.

- [x] C145. Reject a rung whose SCF did not converge (D21(c)).
  DONE 2026-08-04.  DESIGN 5.7 + 7.8 step 3b led, then PSEUDOCODE
  4e.7 and 15.7, then the code.
  **Two different questions, answered in two different places.**
  The flight entry's status says whether the JOB completed;
  `converged` in the run's own result.toml says whether the SCF
  reached its fixed point.  A run that hits its iteration ceiling
  does both -- exits cleanly AND writes a total energy -- so it
  passed the only check either path made, and from that point its
  energy was indistinguishable from a real one.  DESIGN 6.2 already
  said such a run's outputs "must not be harvested as a reference
  potential"; neither path enforced it.
  **Why an unconverged energy is worse than useless on a ladder.**
  It is wherever the iteration happened to stop, which inverts what
  the flatness test asks -- the test wants an energy that has
  stopped moving with the MESH, and this one stopped moving for a
  reason that has nothing to do with the mesh.  It can read flat by
  coincidence and it can break a plateau that was real.
  **The two paths need opposite treatments, and get them.**  The
  climb treats it as a rung that did not run, which stops the
  material: the climb picks its next mesh FROM the ladder, so a
  ladder that does not grow would re-request the same mesh forever.
  The standalone density harvest DROPS the point and carries on --
  its grid is a fixed set of points, not a sequence that chooses
  its next member, so removing one cannot stall anything.
  **Missing data never discards.**  Only an EXPLICIT `converged =
  false` rejects.  A result.toml with no such field cannot be
  judged and is kept -- the same side taken on a missing `gap_ev`,
  and for the same reason.  An older workspace behaves exactly as
  before rather than stopping every material in the campaign.
  Neither drop is silent: the climb prints the material and mesh,
  the harvest adds a summary line naming the dropped densities.
  Also fixed while there: the harvest's single-point guard counted
  the REQUESTED points, so after filtering it could hand an empty
  grid downstream.  It now counts survivors and catches zero.
  +4 tests, suite green at 1271.  The sharpest is the harvest one:
  four converged points hold a real plateau and one unconverged
  point sits in the middle of it, so nothing converges at all while
  it is on the ladder.
  CODE; DESIGN 5.7 / 7.8, PSEUDOCODE 4e.7 / 15.7.

- [x] C144. Carry the climb's verdict forward instead of discarding
  it.  DONE 2026-08-04.  DESIGN 5.7 + 7.8 led, then PSEUDOCODE 4e.5 /
  11.4 / 15.7, then the code.
  **The defect.** `converge_by_climb` computed WHY a material stopped
  and threw it away one line later: CONVERGED and METAL both wrote
  `outcomes[m] = action.rung` and nothing else.  Two different stops,
  one indistinguishable record -- a settled rung carries no trace of
  which produced it.  So every later stage had to work the
  classification out again from whatever evidence it held, and the
  harvest holds ONE rung.
  **Three things now carry it.**  (a) `converge_by_climb` returns a
  third dict, `verdicts[m]`, written by `retire()` alongside the
  outcome so an outcome can never be recorded without its reason.
  (b) The run log records `verdict` verbatim and DERIVES `converged`
  from it -- so `converged` finally means k-point converged and reads
  false for a metal.  Both fields are kept because they answer
  different questions: `converged` gates the 5.8 harness, `verdict`
  separates the two false cases (a metal row names a mesh and yielded
  a potential; a not_converged row yielded nothing).  (c)
  `build_entry` takes `ladder_is_metal`.
  **Either reading says metal (decision, 2026-08-04).**
  `ladder_is_metal` is the CALLER's multi-rung reading and the chosen
  rung's own gap is the second; either suffices, neither can overrule
  the other.  Default false = "no evidence offered", never "known not
  to be a metal", so nothing previously caught stops being caught.
  **The standalone harvest re-derives nothing either.**  It has no
  climb to read a verdict from, so it makes the multi-rung reading
  itself over its own grid's gaps -- same rule, on the evidence that
  path holds.  Those gaps were always parsed there; only the chosen
  one had ever been consulted.  This closes the same hole on the
  second path.
  +4 tests, suite green at 1267.  The sharpest is
  `test_a_sweep_is_metallic_if_any_point_is`: chosen point 0.124 eV
  against a 0.05 cut, a coarser point at zero -- fcc Al's ladder, and
  the shape that recorded Al and Cu as gapped insulators.
  **Fixed a defect C142 introduced.**
  `climb_policy_from_manifest` rejected `metal_gap_threshold <= 0`, so
  the negative-disables escape hatch C142 documented could not be set
  from a manifest at all -- only by building a ClimbConfig directly,
  which is what its test did.  The range check is gone: every real
  value is meaningful, and a `> 0` check rejected exactly the setting
  DESIGN 3.12.6 tells a curator to use.  Zero is meaningful too (a
  true metal's gap collapses to exactly zero).
  Deferred, by decision: marking the harvested POTENTIAL itself as a
  metal's.  Nothing in the potential database distinguishes a
  deliberately rough metal potential from a converged one, but fixing
  it means a schema change, and not now.
  CODE; DESIGN 5.7 / 7.8, PSEUDOCODE 4e.4 / 4e.5 / 11.4 / 15.7.

- [x] C143. Settle a metal on the rung that READ gapless, not on the
  densest rung computed.  DONE 2026-08-04.
  **This was a drift, not a design change.**  DESIGN 3.12.3 has said
  since it was written that the trigger is "the first rung that
  actually reads zero, wherever on the climb that falls" and that the
  material "settles there".  PSEUDOCODE 4e.3 said `METAL(rungs[-1])`,
  glossed as "the densest rung reached so far", and the code followed
  the pseudocode.  So the two lower levels had drifted from DESIGN and
  agreed with each other, which is the shape of drift that survives
  review.  Corrected downward, per the chain: DESIGN was the level in
  the right, so PSEUDOCODE and the code moved to it.
  **The two readings coincide on a plain upward walk**, which is how
  the distinction was lost.  They part wherever more than one rung is
  in hand when the test fires: a confident opening GRID resolves
  several rungs at once, and the bracket-refine REFINE phase fills
  from the bottom, landing rungs BELOW ones already computed.  There
  `rungs[-1]` can be a rung that read a GAP while a lower rung read
  zero.
  **What that cost.**  The harvest re-reads the ONE settled rung to
  decide whether to stage a guidance entry.  Handed a gapped rung it
  sees an insulator where the climb saw a metal, and stages a k-point
  convergence claim for a material that has none -- with the
  finite-mesh artifact recorded as a real gap, which is a predictor
  key (DESIGN 7.6).  This is the Al/Cu failure of 3.12.3 reached by a
  second route.  Settling on the gapless rung makes the climb and the
  harvest agree BY CONSTRUCTION rather than by coincidence, which is
  why it was worth doing before the larger verdict-carrying work.
  Also the cheaper rung, and roughness is the intent.
  DESIGN 3.12.3 gained a short paragraph naming which rung is taken
  when several are in hand, so the same drift cannot recur silently.
  +1 test (`test_metal_settles_on_the_gapless_rung_not_the_densest`),
  covering both search shapes since the rule sits in `climb_next`
  above the dispatch.  Suite green at 1263.
  Does NOT close D22, and does not remove the need to carry the
  climb's verdict forward: the run log's `"converged": True` at
  `make_run_log_entry` is a hardcoded literal that no amount of
  climb/harvest agreement can fix.  Settled follow-on decisions,
  2026-08-04: a material counts as a metal if EITHER the climb or the
  harvest says so; the standalone harvest READS the verdict off the
  prediction record rather than re-deriving it; and marking the
  harvested potential itself as a metal's is deferred -- no potential
  database schema change for now.
  CODE; DESIGN 3.12.3, PSEUDOCODE 4e.3.

- [x] C142. Apply the metal test to every climb shape, not just the
  automatic one.  DONE 2026-08-04.  DESIGN 3.12.3 led, then
  PSEUDOCODE 4e.3, then the code.
  **The old scoping confused two different things.**  DESIGN 3.12.3
  put the gap test in the bracket-refine climb alone, reasoning that
  the fine unit-step walk is the conservative shape a curator pins on
  purpose and should compute every rung.  But a STOPPING RULE is
  properly a matter of search shape -- the shapes exist to disagree
  about which rungs are worth computing -- while a CLASSIFICATION is
  not.  Recognising a metal is the second kind.
  **And the omission cost a result, not just time.**  Without the
  test the unit-step climb walks a metal to the ceiling hunting a
  flatness a metal does not have, stops non-converged, and harvests
  NO potential at all -- strictly worse than the rough floor-level
  potential the metal path exists to produce, and paid for in full.
  It also left the harvest's single-rung gap reading as the only
  metal judgment anywhere in the system, which is the reading that
  recorded Al and Cu as insulators (3.12.3, and D22).
  **One test, one place.**  It moved OUT of `bracket_refine_next` and
  UP into `climb_next`, above the shape dispatch, so all three shapes
  -- grid, bracket-refine, unit-step -- share one copy of the rule
  and none can grow a variant.  The recursive resume at the foot of
  `bracket_refine_next` may skip it safely: the rungs it re-judges
  are the ones `climb_next` already read.
  **The diagnostic ladder survives, by an existing dial.**  A band
  gap cannot be negative, so `metal_gap_threshold = -1.0` is a test
  no rung can trigger -- how a curator asks for every rung of a KNOWN
  metal to be computed, as the seed-run gap ladders did.  No special
  case was needed; it falls out of the comparison.  DESIGN 3.12.6
  now states it, and a test pins it.
  +2 tests (`test_unit_step_settles_a_metal_on_its_gap`,
  `test_negative_gap_threshold_disables_the_metal_test`); suite green
  at 1262.  The pair is the demonstration: same energy model, same
  zero gap, and the only difference is the threshold's sign.
  Does NOT close D22.  This fixes which climbs ASK the question; D22
  is that the climb and the harvest ask it of different evidence (any
  rung versus the one settled rung), which item 3 addresses.
  CODE; DESIGN 3.12.3 / 3.12.6, PSEUDOCODE 4e.3.

- [x] C141. Make linear tetrahedral integration the authored
  default for the initial-potential producer.  DONE 2026-08-04.
  DESIGN 5.7 led, then PSEUDOCODE 11.6, then the code constant --
  the same order C140 took.
  **The reason is that the choice has to be made before the answer
  is known.**  The producer cannot tell a metal from an insulator
  until the ladder has climbed several rungs, so the integration
  scheme is fixed while the question is still open.  That rules out
  "pick the one that suits the commoner case" and leaves "pick the
  one that is safe under both".  Unsmeared Gaussian integration
  moves whole states across the Fermi level as the mesh refines and
  rattles the energy by amounts that do not shrink with the mesh
  spacing -- the noise floor C140's threshold had to be raised
  above.  Tetrahedral integration varies the occupied volume
  continuously instead.
  **It costs the insulators nothing**, which is what makes it a
  free choice rather than a trade: in a gapped system every
  tetrahedron is wholly occupied or wholly empty, the Bloechl
  weights reduce to a quarter per corner, and the scheme returns
  the Gaussian answer exactly.  Measured on `si_fd-3m_227_2001` at
  mesh 6-6-6, where the total energy is unchanged between schemes
  (C138's insulator check).
  **Scoped to the ground state.**  DESIGN 5.7 now says so
  explicitly: a XANES/ELNES run names `"gaussian"` and gets it,
  because `populateLAT` refuses the core-hole path deliberately
  (C138(b)).  The default is the producer's and the ground-state
  SCF's, not a global one.
  Three touches only: DESIGN 5.7's `kpoint_integration` bullet
  gains the decision and its reasons; PSEUDOCODE 11.6 gains the
  four `[defaults]` run-setting constants and `default_run_settings`,
  which it had referenced but never specified; and
  `curation_manifest.py` changes the one token.  These constants are
  what a NEWLY AUTHORED manifest says, not resolve-time fallbacks --
  rule 2 makes every run setting resolvable from the solid or
  `[defaults]`, so a manifest naming `gaussian` is still honoured as
  written, and the existing seed manifests are unaffected.
  One test changed, `test_default_helpers_match_authoring_values`.
  Suite green at 1260.
  CODE; DESIGN 5.7, PSEUDOCODE 11.6.

- [x] C140. Raise `DEFAULT_KPOINT_CONVERGENCE_THRESHOLD` from
  5.0e-4 to 2e-3 (`curation_manifest.py`), per DESIGN 3.12.3.
  DONE 2026-07-31.  PSEUDOCODE led the code, as the chain requires:
  it pinned 5e-4 in three places (the constant, the
  resolvability-exemption comment, and the ClimbConfig field note)
  and all three were amended first.  One test failed on the
  change, `test_build_initial_potentials_resolves_defaults`,
  asserting the literal 5.0e-4; it now asserts against the
  CONSTANT, since what it cares about is that the resolved value
  reaches the predictor and the default itself is a floor DESIGN
  may move again.  Suite green at 1260.
  The design now states the bar as a FLOOR set by the ladder's own
  rung-to-rung scatter rather than as a preference: measured
  scatter is 0.0008 to 0.0047 eV/atom depending on solid and
  scheme, so 5.0e-4 sat beneath the noise and could be cleared
  only by coincidence.  Evidence and the confirmation run are in
  D21; the design change is committed, this is the code catching
  up.
  Small change, three things to get right.  The docstrings at
  lines 173 and 381 quote the value's role and should say WHY the
  number is what it is, so the next reader does not tune it back
  down as a taste.  Manifests that pin
  `kpoint_convergence_threshold` explicitly are unaffected, which
  includes the two seed manifests already pinned at 2e-3 -- so
  this changes behaviour only for manifests that say nothing.
  And nothing needs regenerating: entries already staged at 5e-4
  carry their own `metric_threshold` (DESIGN 7.2), so the
  dataspace stays self-describing across the change rather than
  silently mixing bars.
  Does NOT close D21.  Raising the threshold makes the ladder's
  non-monotonicity stop mattering; D21(a) (order rungs by
  irreducible count), D21(b) (extend by one rung at the ceiling)
  and D21(c) (exclude failed-SCF rungs) all survive it, and (c)
  is a plain defect at any threshold.
  CODE; DESIGN 3.12.3.

- [x] C139. Make a commit clear the run-directory root copies it
  supersedes.  **VERIFIED LIVE 2026-08-04 -- the acceptance below was
  run and it PASSES.**  A cache MISS over a surviving run directory
  had been re-running the unit against the PREVIOUS calculation's
  inputs.
  `commit_prepared_inputs` refreshes the staged `inputs/`, but the
  flattened copies at the run-directory root are what `imago.py`
  actually reads, and it only ever populates those when they are
  absent -- so they keep the old contents and the engine runs the old
  physics while the key file, the run's own `summary`, and the flight
  report all describe the new one.  Nothing prints.
  **Why this is not a stale hit.**  The hit-test byte-compares the
  ROOT copy, the same file the engine reads, so the directory misses
  CORRECTLY and forever: full price paid on every re-run, old answer
  returned every time.  `--force` cannot help -- it only turns hits
  into misses, and this is already a miss.  Deleting the run directory
  is the only recovery, which is exactly the local-file reuse the
  cache exists to provide.
  **How it surfaced.**  2026-07-31, si_cmce_64_1999 re-run with
  `kpoint_integration` changed from `linear-tetrahedral` to
  `gaussian` for C138(a)'s baseline.  All 27 rungs missed on
  `kp-scf.dat` exactly as C135 intends, re-ran, and every root
  `kp-scf.dat` still read `KPOINT_INTG_CODE 1`.  All 27 energies
  reproduced the tetrahedral ladder to the eighth decimal, the SCF
  output still printed its tetrahedron population line, and the
  `summary` beside it said `SCF KP Integration = Gaussian`.  Eight
  minutes of compute returned a verbatim copy of the ladder it was
  meant to be compared against.
  **Note what C135 did and did not buy.**  The cache key was right --
  it distinguished the schemes and refused the hit.  It bought a
  correct DECISION and no correct result, because nothing between
  that decision and the engine acted on it.  A key that misses into
  the wrong inputs is not a smaller version of a key that hits
  wrongly; it is the same wrong answer with the compute paid for.
  **Scope.**  Only the staged input names are cleared.  A prior run's
  outputs stay -- including the converged potential, which is a
  starting point every later SCF re-converges (DESIGN 6.2.5), and
  which is absent from `inputs/` so it survives without a carve-out.
  This is not "wipe the run directory": 6.1's within-directory
  checkpointing depends on it surviving.  Nor is it makeinput's job;
  makeinput never owns the root, and on the producer path it builds a
  prepare directory and never sees the run directory at all.
  A `--reset` for makeinput that guarantees a clean build directory
  is a fair convenience for hand-run rebuilds and for the wingbeat's
  build-in-place branch, but it does not reach these files and must
  not be logged as this fix.
  CODE (`src/scripts/kaleidoscope/wingbeats.py`);
  DESIGN 6.2.5 ("What a commit owes a surviving run directory");
  PSEUDOCODE 13.2 (`clear_superseded_root_copies`).
  **Acceptance:** run one solid at one mesh, change only
  `kpoint_integration`, re-run over the surviving workspace, and
  confirm the root `kp-scf.dat` carries the new `KPOINT_INTG_CODE`
  and the energy moves.
  **RESULT, 2026-08-04.  All three observables PASS.**  The check is
  kept at `jobs/c139_check/` -- manifest, README, and a
  `before_state.txt` recording the directory as it stood at launch, so
  the comparison rests on a reading rather than on memory.  It reused
  the si_cmce_64_1999 run at [4,4,3] left by the corrected-LAT ladder
  and asked for `gaussian` over the top of it, which IS the
  re-run-over-a-surviving-directory this item asks for.
  Both candidate energies were already measured on this binary at
  this mesh, so the numbers were PREDICTED, not read after the fact:
        observable            before      predicted   actual
        root KPOINT_INTG_CODE    1            0          0
        total_energy         -30.8009762  -30.798312  -30.79831196
        summary SCF KP Intg     LAT       Gaussian    Gaussian
  The energy is the Gaussian value to every predicted digit -- the
  other scheme's answer exactly, not a generously-read near miss.
  Two corroborations point the same way: `scf_iterations` 12 -> 5,
  and `gap_ev` 0.05881958 -> 0.05880924, the small shift two
  integration schemes should give at one mesh.  `inputs/kp-scf.dat`
  and the root copy now BOTH read 0, so the commit refreshed the
  staged inputs and cleared the superseded root copy.
  The cache behaved as predicted too -- the driver logged `0 to
  reuse, 1 to run`.  The miss was correct, which is why `--force`
  could never have rescued this.
  **The three unit tests had been passing against a mocked
  filesystem the whole time.**  Every re-run since the fix landed was
  handed an EMPTIED run directory, which is precisely the workaround
  that hides the defect, so none of them had touched the real path.
  That is the general lesson worth keeping: a fix whose only
  exercise is a test that mocks the thing it fixes is unverified.

- [ ] C138. Finish and validate LAT in the SCF occupation path.
  The path itself is CODED and building (DESIGN 1.6, PSEUDOCODE 3a,
  commit `9da3000`): `populateStates` dispatches to `populateLAT`
  when the integration code selects it, `populateLAT` finds the
  Fermi level from the tetrahedron integral by safeguarded Newton,
  `latElectronCount` returns the count and its derivative from one
  corner sort, and `valeCharge` branches at the unpack applying the
  same `2/spin` conversion `computeBond` uses.  The Bloechl corner
  routines moved from `O_DOS` to `O_MathSubs` unchanged, because
  `O_DOS` already uses `O_Populate` and the reverse import would
  have been a module cycle.
  **Verified so far (insulator only).**  `si_fd-3m_227_2001` at
  mesh 6-6-6 under `scfkpint = 1`: the LAT electron count comes out
  at exactly `NUM_ELECTRONS` (8.0 -- the calibration DESIGN 1.6d
  asks for, and the one a wrong `spinFactor` would have failed),
  the run converges in 5 iterations, and the total energy is
  UNCHANGED from the Gaussian run at the same mesh.  That last is
  the correct result, not a null one: in a gapped system every
  tetrahedron is wholly occupied or wholly empty, so LAT must be
  neutral there.  It is evidence the substitution is sound where it
  should change nothing.
  **What remains:**
  (a) **The metal comparison -- DONE 2026-07-31, and it answers
  the question.**  Three ladders for `si_cmce_64_1999` are kept
  under `jobs/si_fingerprint/seed/ladders/` -- `gaussian`,
  `linear-tetrahedral`, `linear-tetrahedral-corrected` -- with
  the analysis in the README beside them.  Same manifest, one
  variable at a time.

        scheme                 verdict      mesh          nk
        gaussian (unsmeared)   converged    [15,15,13]   448
        LAT uncorrected        NOT conv.    [18,18,15]*  720
        LAT corrected (C4a)    converged    [14,14,13]   392
        * the per-axis ceiling, not convergence

  **The schemes return different verdicts on the same solid**,
  which is what this item was waiting for: LAT changes an answer
  it ought to be able to change, at the level of the verdict
  rather than the sixth decimal.  Gaussian's convergence is
  declared on scatter -- from [10,10,8] up it moves +0.0005
  eV/atom in total across thirteen strides while individual
  strides reach +/-0.0028, having hit its noise floor near 120
  k-points.  The [12,12,11] false positive reproduces EXACTLY
  (-0.00042 down, -0.00012 up, flat both sides at 252 k-points,
  which `flat_needed = 1` would have taken); neither LAT ladder
  is catchable by it.  The electron count reads 16.000000 at
  every mesh and iteration in both LAT runs.
  **Read the agreement, not the verdict.**  All three verdicts
  rest on a flatness test over a plateau whose scatter is near
  the threshold, so none is strong evidence by itself.  What is
  strong: from ~48 k-points up, corrected LAT tracks Gaussian to
  within +/-0.0006 eV/atom, while uncorrected LAT sat +0.0084
  away at 504 and was still moving.  Two approximations to one BZ
  integral must converge to each other; with the correction they
  do.  That test does not use the flatness rule at all.
  Caveat on the baseline: the Gaussian ladder is unsmeared, and
  smearing is a Gaussian-path setting only (DESIGN 1.6e --
  `populate` takes the LAT branch before it tests
  `thermalSigma`), so a smeared-versus-LAT run is a DIFFERENT
  comparison and not a refinement of this one.
  A first attempt at the Gaussian baseline returned the
  tetrahedral ladder verbatim while reporting Gaussian; that
  defect is C139.
  (b) **XANES/ELNES under LAT is refused, deliberately.**
  `populateLAT` stops with a message naming Gaussian integration as
  the supported path.  The Gaussian core-hole correction addresses
  the flat sorted occupation array through
  `indexEnergyEigenValues`, and its band arithmetic does not
  transfer to the `(band, kpoint, spin)` array unexamined --
  `numOrbitalStates` is scaled by `spin` at initialization, so the
  mapping is not the obvious one.  A wrong correction misplaces
  exactly one electron, which reads as a convergence problem and
  never gets questioned.  PSEUDOCODE 3a keeps the shape as the
  spec and marks the band range as the unsettled part.
  (c) DESIGN 1.6's two remaining open questions: whether SCF and
  post-SCF integration schemes may legitimately differ, and tying
  the Fermi root-find tolerance to the SCF convergence criterion
  rather than leaving it fixed at `smallThresh`.
  (d) No automated coverage exists for any of this -- the test
  suite is Python and this is Fortran.  The insulator electron-count
  check is currently a thing a person reads out of `gs_scf-fb.out`.
  (e) **The gapped-mesh discrepancy -- RESOLVED 2026-07-31, the
  premise was false.**  As posed: at `si_cmce_64_1999` mesh
  [4,4,3] both schemes report the SAME gap (0.0588 eV) and yet
  differed by 0.031 eV/atom, already in the FIRST SCF iteration
  (-30.581 vs -30.448 hartree) from the same starting potential,
  with Fermi levels 0.099 eV apart -- more than the reported gap
  is wide.  The argument for why that should be impossible: with
  a true gap at E_F every tetrahedron is wholly occupied or
  wholly empty, the Blochl weights reduce to 1/4 per corner
  (Case 0b), summing over the tetrahedra sharing each k-point
  gives back the uniform mesh weight, and LAT must return the
  Gaussian answer -- as the si_fd-3m insulator at 6-6-6 showed.
  **C4a settled it, by a route this entry got backwards.**  This
  item originally ruled the correction out as an explanation,
  reasoning that it vanishes wherever the neutrality argument
  applies.  True, and that is exactly what makes it decisive: the
  correction moved [4,4,3] by 0.040 eV/atom, and since it is
  proportional to a tetrahedron's DOS at the Fermi level, a
  nonzero shift PROVES tetrahedra straddle E_F at that mesh.  So
  there is no true gap at the Fermi level there, the neutrality
  argument never applied, and the premise of the puzzle was
  false.  A term that vanishes under a hypothesis is a test OF
  that hypothesis, not merely something the hypothesis excludes.
  With the correction the same-mesh disagreement falls from
  +0.031 to -0.009 eV/atom, ordinary for two approximations at 12
  k-points, and the `fullKPToIBZKPMap` misweighting branch is
  correspondingly weakened.
  **`gap_ev` is not implicated, and this entry briefly said it
  was.**  It reports exactly what it documents:
  `populate.F90:482` takes the spacing between the highest
  occupied and lowest unoccupied entries of the GLOBALLY sorted
  eigenvalue list, and collapses it to zero below
  `metalGapThresh = 1.0e-3` a.u.  The comment at
  `populate.F90:253-265` describes this case in advance -- a true
  metal on a discrete mesh shows a small artificial "gap" of
  order the level spacing at the Fermi energy, ~1e-4 to 1e-2
  a.u. depending on mesh density, and one exceeding the cutoff is
  "a k-point convergence problem to be cured with a denser mesh".
  0.0588 eV is 0.00216 a.u.: twice the cutoff, inside the
  predicted band, far below the >= 0.5 eV named for a real
  semiconductor gap.  The prescribed cure works -- the reported
  gap is 0.0000 from [5,5,5] up in all three ladders.
  So the coarse-mesh "gap" is a documented finite-mesh artifact,
  not a gap at E_F, which is precisely why the neutrality
  argument did not apply.  Two independent lines agree on that:
  the correction being nonzero, and `metalGapThresh`'s own
  account of what a sub-0.03 eV gap means at 12 k-points.
  Also retracted: "both schemes cannot sit inside the same gap"
  confused a SPACING with a position -- two runs converging to
  different potentials can report the same spacing at different
  absolute energies, which is all the 0.099 eV Fermi-level
  difference says.  And C125's metal classification is not at
  risk, since it reads the settled mesh, by which point the
  artifact is long gone.
  CODE (Fortran); DESIGN 1.6; PSEUDOCODE 3a.

- [ ] C136. Let a run suppress the k-mesh reduction while keeping
  its atomic symmetry.  **The invariant both this and C137 serve:**
  the point group used to reduce the k-point mesh must be a symmetry
  of the Hamiltonian ACTUALLY being solved -- the geometry, the type
  assignment, AND the electronic configuration.  Today the reduction
  group is a pure function of the declared space group:
  `makeinput._extract_point_ops` reads it via
  `symmetry.read_conv_abc_point_ops(space_db, space_group_name)` and
  writes it into `kp-scf.dat`, and `kpoints.f90` reads it back
  (line ~290) without ever consulting atom types, potential types, or
  the electronic state.  So a Hamiltonian that has LESS symmetry than
  its geometry is reduced by operations it does not possess.
  The motivating case is XANES/ELNES.  A core hole on one atom leaves
  the geometry's full space group intact while the true electronic
  symmetry is only the subgroup fixing the excited site, so
  eigenvalues at `k` and `R.k` genuinely differ.  The reduction is
  then invalid at the EIGENVALUE level -- even a TDOS is wrong -- and
  no permutation bookkeeping repairs it.  A full mesh is the fix.
  **What the option does.**  The `.skl` still declares a crystal with
  its space group, and the structure is produced exactly as now,
  INCLUDING symmetry-based assignment of types to equivalent atoms.
  Only the k-mesh side is changed: the writer emits the identity
  alone into the `NUM_POINT_OPS` / `POINT_OPS` block, so the run
  integrates a P1-effective full mesh.  This separates two uses of
  the space group that are currently one act, and lets a calculation
  say "this structure HAS symmetry that this calculation must not
  use" -- a statement declaring `P1` cannot make, because that also
  changes the structure and the types.
  Note the cost is real (the full mesh is |G| times the IBZ) and
  should be reported, not silent.  Note also that the option changes
  `kp-scf.dat` and NOTHING else, so without C135 a symmetry-off run
  is a cache HIT against a symmetry-on one -- the same fault, and a
  second reason C135 must land first.
  DESIGN 2 (the invariant) and 3.2 (mesh reduction); PSEUDOCODE 4b.1
  (writer additions); then CODE (makeinput only -- imago simply reads
  fewer operations).

- [ ] C137. Refuse to run when the type assignment splits a symmetry
  orbit the k-mesh reduction relies on.  The second way to violate
  C136's invariant, and the silent one.  Giving two symmetry-
  equivalent atoms different types gives them different potentials,
  which IS a symmetry breaking of the Hamiltonian -- so a mesh
  reduced by the operation exchanging them is invalid for exactly the
  reason the core hole is.  Types COARSER than symmetry stay
  perfectly safe: a type spanning a union of whole orbits leaves the
  orbit sums invariant, which is what makes the SCF's `potRho`
  accumulation (indexed by potential TYPE, not by atom) exact under
  the IBZ reduction in the first place.  Only splitting an orbit
  hurts.
  **The check.**  `makeinput` holds both halves already -- the
  operations and the type assignment -- so it belongs there, at build
  time, before anything dispatches (the fail-fast shape of DESIGN
  5.10.6).  For each operation `R` in the reduction set and each atom
  `A`, if `type(A) != type(R.A)` the orbit is split: refuse, naming
  the two atoms, the operation, and the two types.  The atom
  permutation under each operation is already built (DESIGN 2.4 /
  2.7).
  **With an override**, because the size of the error is worth
  measuring rather than assuming: a flag that downgrades the refusal
  to a loud warning and runs anyway, so the same system can be run
  split-and-reduced against split-and-full (C136) and the difference
  read off.  The refusal message should name C136's option as the
  correct fix, since suppressing the reduction is what actually makes
  a split-orbit run sound.
  **Not a complete implementation of the invariant.**  This catches
  the type-assignment route only; the electronic route (a core hole
  with types intact) is invisible to it and is C136's to handle by
  hand.  Say so where the guard is documented, so a later reader does
  not mistake a passing check for a guarantee.
  DESIGN 2; PSEUDOCODE 4b.1 plus the section-4 permutation table;
  then CODE.

---

## OPTICAL PROPERTIES (imported OLCAO code)

Items in this section sit outside the document chain.  The
optical properties code (`optc.F90`, `optcPrint.F90`,
`imagoKKc.f90`, `processPOPTC.py`) came across from OLCAO and
has never had a DESIGN or PSEUDOCODE section written above it.
Bringing it into the chain is itself a future task; until then
these entries record what is known to be open so the knowledge
is not carried only in conversation.

- [x] O1. Settle the Kramers-Kronig prefactor in
  `imagoKKc.f90`.  RESOLVED 2026-08-04 by inspection, with no
  measurement needed: **the prefactor is correct.**  No code
  changed; the two factors it multiplies together are now
  documented where it is assigned and where it is used.

  `kramersKronig` sets `multFactor = 2.0/3.0/pi` and applies
  it to each Cartesian component of eps1 separately, while
  `averageFunctions` divides the three components by 3 again
  to form `totalEps1`.  Read as a single number that is a
  factor of three applied twice.  It is not a single number.

  **The `1/3` is Simpson's third, not a directional average.**
  The quadrature multiplies by `fineEnergyDiff`, which
  `getFineEnergyDiff` computes as the mean fine-grid spacing
  -- the bare `h`, not `h/3`.  Composite Simpson needs `h/3`,
  so the missing third is carried in `multFactor` alongside
  the `2/pi` that is the actual Kramers-Kronig prefactor.
  Each component therefore receives `2/pi` exactly once, and
  the only average over x, y and z is the one in
  `averageFunctions`.  The historical agreement of computed
  eps1(0) with measured refractive indices is not a
  coincidence in need of explanation; it is what a correct
  calculation produces.

  Also settled while reading, since it is the other unlabelled
  constant in the same expression: the trailing `pOptcFactor`
  is the additive 1 of `eps1 = 1 + (2/pi)*Int[...]`.  For a
  total run it is literally 1.  For a partial run `optc.F90`
  writes each pair a share, `partialsIndex(j)*partialsIndex(k)
  / valeDimIndex**2`, which sums over all pairs to exactly 1
  -- so a complete set of partials carries the constant once
  between them rather than once each.  `partialsIndex` is
  `real`, so this is real division and not the integer
  division it resembles.

  What this does NOT settle is whether the quadrature that
  the `1/3` belongs to is itself right; see O6.

- [ ] O2. Apply IBZ atom-permutation unfolding to the
  atom-resolved POPTC detail codes.  DESIGN 2.5 already
  settled the governing rule for PDOS, and it carries over to
  POPTC once each detail code is matched to its PDOS analogue.

  **Specified in PSEUDOCODE 7a**, which carries the detail
  code table, the correction itself, the proof that the
  component-summed pair matrix transforms by pure index
  permutation, the guards, and the cost.  Code against that
  section, not against the notes below; what follows is the
  reasoning that produced it and the state of the task.

  **State.** The isotropic half is CODED, in
  `computePOPTCPairs` immediately after the transition double
  loop, and VERIFIED live on 2026-08-05.  The style code 0
  warning in `kpoints.f90` now names partial optical
  properties among the decompositions it covers.  The entry
  stays open for the per-axis half, which O3 governs.

  **The verification, and why it is believable.**  Cubic
  KNbO3 (Pm-3m, 5 atoms), whose three oxygens are one type on
  the 3c site -- symmetry equivalent, so their isotropic
  spectra must match, but each on a different Cartesian axis,
  so the cubic operations permute them non-trivially.  K and
  Nb are alone in their types and act as controls.  Two runs
  of a 4x4x4 shifted mesh differing in one thing only: the
  post-SCF k-points folded to 4 points against the same mesh
  left whole at 64.  The unfolded reference was produced by
  trimming the operation list in `kp-pscf.dat` to the
  identity, which leaves the fold nothing to fold with and
  keeps both runs on the same code path.  Both took the same
  4-point SCF; the potential is expanded per type, so the
  three oxygens share one by construction and no potential
  asymmetry can leak into the comparison.  Job directory
  `jobs/knbo3/o2_poptc_unfold` (gitignored), 83 and 86
  seconds.

  Isotropic eps2, worst relative disagreement over all 25
  atom pairs: 1.2e-8, which is the print precision of the
  `.raw` files -- agreement to the last written digit.

  Three findings make that a pass rather than a coincidence:

  1. **The comparison is demonstrably sensitive.**  The
     per-axis columns of the same files, read by the same
     parser, disagree by 1.08, 1.29 and 1.45 relative.  Those
     columns are the ones the correction deliberately does
     not repair, so a comparison that could not see a
     difference would have shown agreement there too.
  2. **The averaging did real work.**  In the reduced run the
     three oxygens agree to 5.8e-17 -- exactly, because they
     were averaged over their orbit.  In the full run they
     agree only to 5.8e-10, the roundoff of summing 64
     independent k-points.  Had the raw IBZ decomposition
     already been permutation-symmetric, the correction would
     have been a no-op and the reduced oxygens would agree at
     the raw data's own level, not at machine epsilon.
  3. **The totals did not move**, as predicted, so the sum
     rule confirmed nothing -- which is the point.

  Not done: a direct A/B against a binary with the block
  disabled.  The three findings above make that redundant
  rather than merely inconvenient, but it is the one check
  that would settle it by construction rather than by
  inference.

  Incidental, found while setting this up and worth its own
  fix: the user-facing POPTC legend disagrees with the code
  on detail code 4.  `makeinput.py` writes
  `0N;1t;2a;3enl;4enlm` into every `imago.dat` and
  `makeinputrc.py` carries `0=NONE,1=elem,2=a tot,3=e nl,4=e
  nlm`, both of which read as though code 4 were element or
  type resolved.  It is atom resolved: `pOptcIndex(vdi) =
  vdi` makes every basis function its own partial, and
  `printSpectrumPOPTC` walks `numAtomSites`.  The labels also
  hide the alternation that matters here -- type, atom,
  type+nl, atom+nlm -- which is why it is codes 2 and 4 that
  need the permutation rather than 3 and 4.  This is O5's
  finding surviving in the strings a user actually reads;
  O5 corrected only the comments inside the Fortran.

  Scope is isotropic now, per-axis blocked behind O3.  The
  atom permutation makes the per-atom-pair isotropic column
  correct; the per-axis columns need the Cartesian rotation,
  which is O3's.  O2 stays open with that half deferred.

  The decisive property is that a sum taken over all atoms of
  a type is invariant under the operations used for IBZ
  reduction, while anything resolving individual atoms is not.
  Note carefully that this does *not* follow from type meaning
  "symmetry equivalent", because in Imago it often does not:
  in an amorphous cell a type is a bin of locally similar
  environments with no symmetry content at all, and in a point
  defect supercell types are deliberately assigned from the
  *pre-defect* symmetry to keep the cost down even though the
  defect has destroyed that symmetry.  DESIGN 2.3 states the
  rule as though species implied equivalence; that premise is
  false in general and should be corrected there.

  What actually guarantees the invariance is a runtime check.
  `buildAtomPerm` in `atomicSites.f90` requires every point
  operation to carry each atom onto an atom carrying the same
  `atomTypeAssn`, and stops with a fatal error if no such
  partner exists.  So whatever the types happen to mean, the
  operations Imago actually reduces by are guaranteed to
  permute atoms within a type -- which is exactly the closure
  property a type-level sum needs.  The two awkward cases both
  land safely for this reason rather than by assumption: an
  amorphous cell carries no operations to reduce by (P1, so
  the permutation is the identity), and a defect supercell
  typed from the original symmetry is *coarser* than the true
  orbits of the reduced group, which is the safe direction --
  a union of orbits is still closed.  Typing *finer* than the
  orbits is the unsafe direction, and that is the case
  `buildAtomPerm` aborts on.

  The gap is style code 0.  With an explicit k-point list
  Imago cannot build the symmetry maps, so neither the
  unfolding nor the `buildAtomPerm` consistency check runs --
  only the warning noted in DESIGN.  A hand-supplied reduced
  k-point list for a defect supercell is therefore the one
  configuration where a wrong type-level decomposition can
  pass silently.  Worth restating in the O2 write-up.

  Mapping the four detail codes onto the rule:

  ```
  code  grouping        PDOS analogue   needs
  ------------------------------------------------------
  1     type            mode 0          nothing extra
  3     type + QN_nl    mode 0          nothing extra
  2     atom            mode 1          atomPerm/invAtomPerm
  4     atom + QN_nlm   mode 3          D^l(R); deferred
  ```

  So codes 1 and 3 are already correct on a symmetry-reduced
  mesh and need no work.  Code 2 needs the same treatment bond
  order and Q* received in Phase F.  Code 4 is blocked behind
  the same deferred D^l(R) representation matrices that block
  PDOS mode 3, and is in any case unaffordable at present (see
  the scaling note below).

  Code 3 is grouped by *type*, not by atom, and every comment
  that said otherwise has been corrected (see O5).  Three
  independent places agree on the reading: `pOptcIndex` is
  built from `cumulNumPartials(currentType)`, the partial count
  is `cumulNumPartials(numAtomTypes + 1)`, and
  `printSpectrumPOPTC` walks `numAtomTypes` and writes a
  `TYPE_1`/`TYPE_2` header.  Only the labels were wrong.  The
  wrong reading would have put code 3 in the unsafe column and
  bought it an unfolding correction it does not need.

  Two things make POPTC harder than PDOS, neither of which
  needs new infrastructure:

  1. The decomposed quantity carries *two* group indices (an
     initial-state group and a final-state group), so an
     operation permutes both at once -- the pair matrix is
     conjugated by the permutation rather than re-indexed
     along a single axis:
       M(a,b) at k_full = M(invAtomPerm(R,a), invAtomPerm(R,b))
                          at k_IBZ
     The existing `invAtomPerm` table supplies this; it is
     simply applied twice.
  2. **The partials-sum-to-total identity cannot be used to
     check this.**  That identity (`optc.F90`, the
     `Re*ReSum + Im*ImSum` construction) holds under *any*
     permutation of the two indices, because permuting the
     summands does not change their sum.  It is exactly blind
     to a missing or wrong unfolding, so a passing sum rule
     proves nothing here.  Check instead against a structure
     with symmetry-equivalent-but-inequivalently-oriented
     atoms, comparing a full-mesh run against an IBZ-reduced
     run of the same system; those must agree per atom.

  Do this after O1 and after the `optc.F90` correctness fixes
  have landed, so the baseline being compared against is
  trustworthy.

- [ ] O3. Decide whether the x, y, and z resolved columns of
  the optical output are meaningful on a symmetry-reduced
  k-point mesh.  This is independent of POPTC and affects the
  *total* spectra that Imago already produces.
  The momentum operator is a vector: under a point group
  operation its Cartesian components mix, P_i(Rk) = sum_j
  R_ij P_j(k).  `getOptcCond` accumulates
  `transitionProb(:,j,i,h) * kPointFactor(i)` over the IBZ
  points only, so every member of a k-point's star is credited
  with the *representative's* orientation, scaled by the star
  multiplicity carried in `kPointWeight`, with no rotation of
  the components.
  The isotropic average is safe: summing x, y and z is a trace
  and is invariant under the rotation, so the "total" column
  that `printSpectrum` and `printSpectrumPOPTC` write (the
  sum over the three components divided by three) is correct.
  The three individual columns are not, unless the crystal is
  cubic or the mesh was never reduced.
  Decide between rotating the components during accumulation
  and documenting the per-axis columns as valid only for
  unreduced meshes.  Until then treat anisotropic per-axis
  eps2, conductivity, n, k, R and alpha as unverified.

- [ ] O4. Repair or remove the serial XYZ path.  `imago.py`
  exposes `serialxyz` as a user setting and passes it through
  to `serialXYZ` in `commandLine.f90`, but the branch it
  selects in `computeTransitions` is entirely commented out.
  A run with the flag set therefore computes no transitions
  at all, prints zero-valued spectra, and then fails in
  cleanup trying to deallocate a `transitionProbPOPTC` that
  was never allocated.  Either restore the branch (its intent
  was to do the three Cartesian components one at a time to
  reduce peak memory, which is directly relevant to the POPTC
  scaling problem) or remove the option so it cannot be
  selected.  Silently producing zeros is the worst of the
  three outcomes.

- [x] O5. Correctness fixes along the partial optical
  properties path.  DONE 2026-08-04.  These are the fixes O2
  waits on: until they land, a full-mesh run and a reduced
  run cannot be compared, because neither one is right.
  Recompiled clean.

  **The one that was silently wrong, not loudly wrong.**
  `finalStateIndex` indexes `conjWaveMomSum`, which is filled
  for *every* final state in the `firstFin` to `lastFin`
  range.  It was a counter incremented only on *accepted*
  pairs.  With thermal smearing a state can be both initial
  and final, and such a pair is skipped -- but its slot in
  `conjWaveMomSum` still exists, so from the first skip
  onward the counter read a different final state's momentum
  sum, for the rest of that initial state, with no symptom.
  It is now derived, `j - firstFin + 1`.  Both copies of the
  transition loop had it: `computePairs` and
  `computePOPTCPairs`.

  **Job IDs 106 and 107 were crossed.**  The `doOPTC` codes
  and the job ID number the last two optical properties in
  opposite order, so 106 (non-linear) was setting code 3
  (sigma(E)) and 107 the reverse.  Both SCF and PSCF.  The
  mapping is now correct and the crossover is commented at
  both sites, since it looks like a typo either way round.

  **Non-linear optical properties now stop rather than
  pretend.**  They have an input block (`NLOP_INPUT_DATA`)
  and a job ID that reaches `getEnergyStatistics`, but no
  routine anywhere computes them -- there is no counterpart
  to `computePairs` for the second-order response.  Falling
  through emitted whichever spectrum the surrounding code
  happened to produce, labelled as a non-linear result.

  Also: `maxPairs` is initialized before being built up with
  `max()`, so the result no longer depends on how a compiler
  pre-fills module storage; the detail-code-1 basis-function
  counting loop no longer forms the out-of-range subscript
  `pOptcIndex(0)` on its first iteration; `transitionProb` is
  released for partial runs too (it is allocated and used
  unconditionally, so the old branch leaked it whenever a
  decomposition was requested); and `valeValeXMom` is
  released alongside `conjWaveMomSum`, which it never was.

  **Two spellings that meant the post-processing simply did
  not run.**  `processPOPTC.py` invoked `makePDOS` and
  `imagoKkc`; the installed names are `makePDOS.py` and
  `imagoKKc`, so on a case-sensitive filesystem neither was
  ever found.  `imago.py` called `processPOPTC` (no suffix)
  and discarded its output, which is how the whole thing
  stayed hidden -- the total spectra still appeared and only
  the partial files quietly went missing.  It now calls
  `processPOPTC.py -s <channel>` and writes the result into
  the runtime log.  Every comment naming the executable has
  been normalized to `imagoKKc`.

  **One more Perl-to-Python port defect.**  `copy_data`
  stripped two leading columns from each `imagoKKc` output
  line.  The Perl original indexed from 2 because Perl's
  `split(/\s+/)` leaves an empty leading token on these
  Fortran-formatted lines and Python's does not, so the port
  was discarding the first real spectral column of every
  partial.  Now strips one.

  **Comment-only, but the point of O2:** code 3 was labelled
  "each QN_nl resolved atom" in four places in `optc.F90` and
  as "Decompose by atom and QN_nl" in `optcPrint.F90`, while
  the code groups by *type*.  Corrected, with the reason the
  distinction matters recorded at the index construction:
  a type-level sum is invariant under the operations used to
  reduce the k-point mesh and an atom-level one is not.

  **Three things noted and deliberately not changed**, each
  because touching it buys nothing and costs a re-verified
  build.  `cumulNumPartials` is never read for detail code 1,
  and the running sum built there adds `i` rather than 1 --
  inert, now commented as such rather than repaired, since
  what it *ought* to hold is undecided.  For code 3 the same
  array is allocated to `numAtomSites + 1` while only
  `numAtomTypes + 1` entries are used; harmless, since the
  site count is always the larger.  And the program unit
  inside `imagoKKc.f90` is still spelled `program imagoKkc`,
  which affects nothing (CMake's target name sets the
  executable name) but should be renamed the next time that
  file is opened for real work.

- [x] O6. Repair the integration in `kramersKronig`.  DONE
  2026-08-04.  Opened out of the O1 reading and closed the same
  day: deciding that the `1/3` in `multFactor` belongs to
  Simpson's rule means the sum it multiplies has to actually be
  a Simpson sum, and reading it closely to check that turned up
  five separate faults, two of them more serious than the one
  that prompted the look.  Line numbers below are as the file
  stood at commit `3bb73f2`, before the repair.

  **How this was measured, since none of it is visible in a
  single run.**  The routine integrates a straight-line
  interpolation of eps2 between the computed energy points, and
  that integral can be written in closed form, so the yardstick
  used here is an exact formula rather than a finer version of
  the same quadrature.  (An earlier attempt used a finer grid
  and was discarded: its own error near the pole was larger
  than the effects being measured, which made a refinement
  study read as though the code stopped converging.)  A Python
  copy of the routine, defect for defect, reproduces the
  compiled program to six figures, so each fault could be
  switched off on its own and its cost read directly.  The test
  spectrum is silicon-like: absorption beginning at a 1.1 eV
  gap, two peaks, 201 energy points to 30 eV.

  **Fault 1, the worst one: the last point of the range was
  evaluated at the wrong energy.**  Line 517 read
  `original_e=energy(length)` where every other point in the
  sum uses `energy(i)`, the energy the eps1 value is being
  computed for.  The consequence is not a small inaccuracy.
  That node sits exactly one fine step below the top of the
  range, so the difference in its denominator was always that
  one step, making the node evaluate to about
  `-eps2(top)/(2h)`.  The `h` in that expression then cancelled
  against the `h` the whole sum is multiplied by, leaving a
  constant `-eps2(top)/(3*pi)` subtracted from **every** eps1
  value.  Being independent of the grid, it was the one error
  that survived any amount of refinement -- refining to a grain
  of 320 left it untouched at 0.067.  Predicted -0.06366,
  measured -0.06386.  It vanishes only when the spectrum is
  carried far enough out that eps2 has decayed to nothing,
  which is why it has never been caught: a well-converged
  spectrum hides it completely, and a truncated one shifts by
  about 0.4 percent of eps1(0).

  **Fault 2: the first point of the range was read but never
  written.**  `makeFineGrainEnergy` sets `fine_energy(1)`
  exactly equal to `energy(1)`, and line 482 compared it
  against `energy(1)` as well, so the guard on line 483 asked
  whether zero exceeds 1e-5 and answered no on every pass of
  every iteration.  The assignment on line 484 therefore never
  executed -- yet line 535 read `integrand(:,1)` into the sum.
  That is a read of memory that was allocated and never
  written.  Confirmed by filling the array with -12345 before
  use, which moved eps1(0) from 16.954 to -22.276, matching the
  predicted shift.  It has been harmless in practice only
  because a fresh allocation usually hands back zeroed memory,
  and zero happens to be the right contribution for a material
  whose energy range starts below its gap.  For a metal the
  right value is not zero, though even there the term is worth
  only about 5e-5.  The danger was never the size of the error;
  it was that the answer depended on memory nobody had set.

  **Fault 3: the Simpson weights were the wrong way round.**
  Counting the first node as number one, Simpson's rule weights
  the even-numbered nodes by four and the odd-numbered interior
  ones by two.  The test on line 503 sent even-numbered nodes
  to `evenSum`, which line 535 multiplied by two, and odd ones
  to `oddSum`, which it multiplied by four.  Exactly doubling
  the error of the integration: measured largest error 0.193
  against 0.097 for the corrected weighting, on a spectrum with
  no truncation, and the same factor of two at every grid
  refinement tested.

  **Fault 4, recorded but NOT repaired: the number of nodes is
  even.**  Simpson's rule needs an odd number, so that the
  points pair up into complete spans.  The sum runs over
  `numValues - grain` nodes, which is `grain*(length-1)` and
  therefore even for the default grain of 10.  Repairing this
  means changing how the fine grid is built, not adjusting a
  weight, and it is why the corrected weighting still converges
  only in proportion to the step size rather than far faster as
  a correct Simpson rule would.  Left alone deliberately: the
  gain is accuracy the calculation does not currently need, and
  the change is to a data structure rather than an expression.

  **Fault 5: two accumulators computed and never read.**
  `totalSum` and `totalSumi` were built up on lines 473, 486,
  489, 510, 511, 521 and 524 and never used.  Deleted.

  **What the repair is worth, compiled program against the
  exact formula, on a spectrum truncated while eps2 is still
  0.60:** eps1(0) moves from 16.95412 to 17.01786 against an
  exact 17.01812, so its error falls from 0.064 to 0.00026.
  Across the spectrum the largest error falls from 0.258 to
  0.098 and the average error from -0.067 to -0.002.  Note what
  the remaining 0.098 is: it is the ordinary error of a
  step-size-limited rule near the pole, it shrinks with the
  grid, and it is no longer a constant bias.

  **This changes every optical spectrum Imago produces.**  The
  shift is largest for spectra cut off while absorption is
  still appreciable, and negligible for spectra carried out to
  where it has died away.  Anyone comparing new output against
  a previously published figure should expect a small upward
  move in eps1 and in everything derived from it.

  Unrelated but in the same expression, and cheap to fix
  whenever this is opened: `valeDimIndex**2` in the
  `POPTC_KKC_FACTOR` write in `optc.F90` is evaluated in
  default integer, so a valence basis above about 46000
  functions overflows it.  Remote, but silent if reached.

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
