"""Tests for the flight -> guidance-entry harvest
(guidance_harvest.py; DESIGN 7.8 / PSEUDOCODE 15.7).

The harvest reads three things off disk -- the flight.toml plan,
each run's result.toml, and the structure .skl -- and stages one
GuidanceEntry per converged structure sweep.  These tests build
synthetic flight workspaces with tmp_path (no Imago binary, no
real run) and monkeypatch the two physics-layer touch points so
the harvest's own logic is what is under test:

* ``compute_signature`` -> a fixed Signature (so no elements.dat
  / $IMAGO_DATA is needed), recording the system_type it was
  asked for;
* ``load_structure`` -> a tiny stand-in carrying ``num_atoms``
  and ``real_cell_volume`` (so no real StructureControl read);
* ``save_entry`` -> a capture that records the GuidanceEntry and
  returns a stub path (so the byte-deterministic emitter, already
  tested in test_guidance_db.py, is not re-exercised and no real
  dataspace tree is needed).

What remains under test is exactly the harvest contract: reading
the swept k-density out of each calc tag, the two-sided
convergence rule, the trust-mode and prediction-mismatch skips,
and the field-by-field assembly of the staged entry.
"""

import os
import types

import pytest

import guidance_harvest as gh
from guidance_db import Signature
from kaleidoscope import CalcUnit, Flight, SweepRecord
from kaleidoscope.workspace import serialize_flight, toml_line


_DATASPACE = types.SimpleNamespace(group_table={})


# --------------------------------------------------------------
#  Workspace builder + physics-layer patch
# --------------------------------------------------------------

def _make_workspace(tmp_path, kpds, energies, *,
                    gaps=None, kinds=None, mags=None, meshes=None,
                    scf_threshold=1.0, write_scf_threshold=True,
                    kpoint_convergence_threshold=5.0e-4,
                    write_gap=True, add_loen=False,
                    policy="verify_around_prediction",
                    system_type="crystalline", confidence=0.9,
                    neighbor_ids=("mp-1", "mp-2"),
                    unit_id="si", structure="si.skl",
                    write_prediction=True):
    """Lay out a one-structure flight workspace under tmp_path: a
    flight.toml (via the real serialize_flight) plus one
    result.toml per grid point.  Returns the workspace root.  The
    knobs let each test shape the grid, the energies, and the
    per-structure [flight.predictions.<id>] block."""
    root = str(tmp_path / "flight")
    units = [CalcUnit(id=unit_id, structure=structure,
                      calc=(f"kpt-density-{k}",)) for k in kpds]
    if add_loen:
        # A structure-only "fingerprint" loen unit sharing the
        #   structure id (DESIGN 6.2.9): the convergence harvest must
        #   filter it out.  It gets no result.toml below, so if the
        #   harvest did NOT skip it the run would fail loudly.
        units.append(CalcUnit(id=unit_id, structure=structure,
                              calc=("loen",), kind="fingerprint"))
    # fixed_axes is empty: the sub-model rides on the per-structure
    #   prediction record (DESIGN 6.2.9), which the harvest reads.
    sweep = SweepRecord(varied_axes=("kpt-density",), fixed_axes={})
    metadata = {}
    if write_prediction:
        # Per-id predictions mapping (DESIGN 6.2.9); this record is
        #   the SOLE home of system_type AND the (basis, functional,
        #   kpoint_integration) sub-model the harvest reads.
        metadata["predictions"] = {unit_id: {
            "policy": policy, "system_type": system_type,
            "confidence": confidence,
            "neighbor_entry_ids": list(neighbor_ids),
            "predicted_kpoint_density": 100.0,
            "is_under_trained": False,
            "basis": "fb", "functional": "gga-pbe",
            "kpoint_integration": "gaussian-0.1",
            # The producer stamps the resolved per-atom k-point
            #   flatness tolerance on the record; the harvest reads
            #   it as the convergence metric_threshold (DESIGN 7.8).
            "kpoint_convergence_threshold":
                kpoint_convergence_threshold}}
    serialize_flight(Flight(root=root, units=units, sweep=sweep,
                            metadata=metadata))

    for index, kpd in enumerate(kpds):
        run_dir = os.path.join(root, "wingbeats", unit_id,
                               f"kpt-density-{kpd}")
        os.makedirs(run_dir, exist_ok=True)
        gap = gaps[index] if gaps is not None else 1.5
        kind = kinds[index] if kinds is not None else "indirect"
        mag = mags[index] if mags is not None else 0.0
        path = os.path.join(run_dir, "result.toml")
        with open(path, "w") as result_file:
            result_file.write(
                toml_line("total_energy", energies[index]))
            if write_gap:
                result_file.write(toml_line("gap_ev", gap))
                result_file.write(toml_line("gap_kind", kind))
            result_file.write(
                toml_line("total_magnetization", mag))
            if meshes is not None:
                # The resolved axial counts imago records for this
                #   grid point (DESIGN 6.1.2), which build_entry
                #   stores as verification.converged_mesh.
                a_count, b_count, c_count = meshes[index]
                result_file.write(
                    f"kpoint_mesh = [{a_count}, {b_count}, "
                    f"{c_count}]\n")
            if write_scf_threshold:
                result_file.write(
                    toml_line("scf_threshold", scf_threshold))
            result_file.write(toml_line("imago_commit", "abc123"))
    return root


@pytest.fixture
def patched(monkeypatch):
    """Patch the harvest's physics-layer touch points and capture
    every staged entry.  Returns a dict with the captured entries
    and the system_types compute_signature was asked for."""
    captured = {"entries": [], "system_types": []}

    def fake_save(entry, db_root):
        captured["entries"].append(entry)
        return os.path.join(
            db_root, "staging", entry.signature.system_type,
            "stub.toml")

    def fake_signature(structure, system_type, group_table,
                       label="<structure>"):
        captured["system_types"].append(system_type)
        crystalline = (system_type == "crystalline")
        return Signature(
            system_type=system_type,
            composition_vector=(1.0,) + (0.0,) * 12,
            lattice_family="cubic" if crystalline else "",
            lattice_onehot=((1.0,) + (0.0,) * 5 if crystalline
                            else (0.0,) * 6))

    def fake_load_structure(path):
        return types.SimpleNamespace(
            num_atoms=8, real_cell_volume=100.0)

    monkeypatch.setattr(gh, "save_entry", fake_save)
    monkeypatch.setattr(gh, "compute_signature", fake_signature)
    monkeypatch.setattr(gh, "load_structure", fake_load_structure)
    return captured


# --------------------------------------------------------------
#  Pure helpers (no workspace needed)
# --------------------------------------------------------------

def test_decode_axis_value_inverts_the_encoding():
    """decode_axis_value undoes encode_axis_value: plain ints,
    'p' as the decimal point, leading 'm' as a minus."""
    assert gh.decode_axis_value("100") == 100
    assert gh.decode_axis_value("1p5") == 1.5
    assert gh.decode_axis_value("m2") == -2
    assert isinstance(gh.decode_axis_value("100"), int)


def test_swept_value_of_reads_the_axis_component():
    """swept_value_of pulls the value from the calc component
    whose prefix matches the swept axis, even though the axis name
    itself contains a hyphen."""
    unit = CalcUnit(id="si", structure="si.skl",
                    calc=("kpt-density-150",))
    assert gh.swept_value_of(unit, "kpt-density") == 150


def test_swept_value_of_missing_axis_raises():
    """A unit whose calc tag has no component for the swept axis is
    malformed -- the harvest refuses rather than guess."""
    unit = CalcUnit(id="si", structure="si.skl",
                    calc=("basis-size-3",))
    with pytest.raises(ValueError):
        gh.swept_value_of(unit, "kpt-density")


def test_pick_converged_two_sided_and_no_endpoints():
    """pick_converged needs BOTH neighbour PER-ATOM-eV deltas below
    threshold and never returns an endpoint.  Energies are raw
    total-cell hartree; with cell_atom_count=1 each delta becomes
    |dE| * HARTREE eV/atom, so a 1.0 eV/atom threshold tolerates the
    small tail deltas but not the big step down to 3.0 (DESIGN 7.8)."""
    # index 1 is flat on the high side but its down-delta to 3.0 is
    #   huge; index 2 is flat on both sides -> chosen.
    energies = [3.0, 2.0, 1.99, 1.985]
    assert gh.pick_converged(energies, 1, 1.0) == 2
    # Still moving everywhere -> None.
    assert gh.pick_converged([3.0, 2.0, 1.0], 1, 1.0) is None


def test_pick_converged_is_the_flat_needed_one_climb():
    """pick_converged is exactly pick_converged_climb with
    flat_needed=1, so the two agree on the same ladder (they share
    one two-sided rule and cannot drift, DESIGN 3.12.3)."""
    energies = [3.0, 2.0, 1.99, 1.985]
    assert gh.pick_converged(energies, 1, 1.0) == \
        gh.pick_converged_climb(energies, 1, 1.0, 1)


def test_climb_demands_flatness_persist_over_two_rungs():
    """With flat_needed=2 the returned index and the next interior
    rung must BOTH be two-sided flat.  Here indices 2 and 3 are both
    flat, so index 2 is the first of a persistent run."""
    # Per-atom deltas (atom_count=1): the 3.0->2.0 step is huge, the
    #   rest are small, so interior rungs 2 and 3 are flat.
    energies = [3.0, 2.0, 1.99, 1.985, 1.983]
    assert gh.pick_converged_climb(energies, 1, 1.0, 2) == 2


def test_climb_rejects_a_single_lucky_flat_rung():
    """A lone flat interior rung satisfies flat_needed=1 but not
    flat_needed=2, so a persistence-demanding climb keeps going
    where a single-grid pick would have stopped."""
    # Index 2 is flat both sides; index 3 jumps up to 5.0, so the
    #   flat run has length one.
    energies = [3.0, 2.0, 1.99, 1.985, 5.0]
    assert gh.pick_converged_climb(energies, 1, 1.0, 1) == 2
    assert gh.pick_converged_climb(energies, 1, 1.0, 2) is None


def test_climb_needs_enough_rungs_above_to_confirm():
    """A flat rung cannot be confirmed at flat_needed=2 until a
    further interior rung exists above it; with only one interior
    rung the climb returns None even though that rung is flat."""
    energies = [3.0, 2.99, 2.985]        # single interior rung (i=1)
    assert gh.pick_converged_climb(energies, 1, 1.0, 1) == 1
    assert gh.pick_converged_climb(energies, 1, 1.0, 2) is None


# --------------------------------------------------------------
#  collapse_by_mesh -- the duplicate-rung guard (DESIGN 7.8 3c)
# --------------------------------------------------------------

def test_collapse_by_mesh_merges_contiguous_duplicates():
    """Two rungs resolving to the same mesh collapse to one --
    the lowest-density member -- and kept maps the survivor back
    to its original index."""
    densities = [100, 150, 200, 250]
    energies = [-10.0, -10.5, -10.5, -10.7]
    meshes = [[2, 2, 2], [3, 3, 3], [3, 3, 3], [4, 4, 4]]
    dens, ergs, kept = gh.collapse_by_mesh(
        densities, energies, meshes)
    assert dens == [100, 150, 250]      # 200 (dup of 150) dropped
    assert ergs == [-10.0, -10.5, -10.7]
    assert kept == [0, 1, 3]            # survivors' original idx


def test_collapse_by_mesh_none_is_inert():
    """If any mesh is None (older result.toml / pre-emit binary)
    the guard cannot act: the grid returns unchanged with
    identity indices, so behavior matches the pre-guard code."""
    densities = [100, 150, 200]
    energies = [-10.0, -10.5, -10.5]
    meshes = [[2, 2, 2], None, [3, 3, 3]]
    dens, ergs, kept = gh.collapse_by_mesh(
        densities, energies, meshes)
    assert dens == densities
    assert ergs == energies
    assert kept == [0, 1, 2]


def test_collapse_by_mesh_energy_mismatch_raises():
    """Equal mesh MUST give equal energy; a disagreement means the
    runs were not identical and is surfaced, not averaged."""
    with pytest.raises(ValueError, match="same mesh"):
        gh.collapse_by_mesh(
            [100, 150], [-10.0, -10.4],
            [[3, 3, 3], [3, 3, 3]])


def test_guard_corrects_duplicate_mesh_false_plateau():
    """The si_ia-3 seed failure, with its real ladder (16-atom
    cell).  Densities 150 and 200 resolved to the SAME mesh
    [2,3,3], so their energies are bit-identical and the raw
    ladder carries a zero delta at 200.  Without the guard the
    two-sided test accepts 200 on that manufactured zero; with
    the guard, 200 is merged into 150 and the test lands on 250 --
    a genuinely two-sided-flat point on distinct meshes (DESIGN
    7.8 step 3c)."""
    densities = [25, 50, 100, 150, 200, 250, 300, 400]
    energies = [-62.0875583, -62.08399468, -62.08286333,
                -62.07818894, -62.07818894, -62.0781121,
                -62.07804731, -62.07761701]
    meshes = [[1, 2, 1], [1, 2, 2], [3, 3, 1], [2, 3, 3],
              [2, 3, 3], [4, 2, 3], [4, 4, 2], [3, 4, 4]]
    threshold = 5e-4                            # 0.5 meV/atom

    # Raw: the duplicate zero at 200 is read as convergence.
    raw = gh.pick_converged(energies, 16, threshold)
    assert densities[raw] == 200                # the false accept

    # Guarded: 200 is merged into 150; the survivor at 250 is
    #   flat on both sides over genuinely distinct meshes.
    dens, ergs, kept = gh.collapse_by_mesh(
        densities, energies, meshes)
    assert 200 not in dens                      # duplicate removed
    guarded = gh.pick_converged(ergs, 16, threshold)
    assert dens[guarded] == 250                 # corrected point


# --------------------------------------------------------------
#  Converged harvest -- the entry's fields
# --------------------------------------------------------------

def test_converged_sweep_stages_one_entry(patched, tmp_path):
    """A flat interior grid point yields exactly one staged entry,
    at the converged k-density."""
    root = _make_workspace(
        tmp_path, [50, 100, 200], [0.5, 0.5, 0.5])
    summaries = gh.harvest_flight(root, str(tmp_path / "db"),
                                  _DATASPACE)
    assert len(patched["entries"]) == 1
    assert len(summaries) == 1 and "staged" in summaries[0]
    entry = patched["entries"][0]
    assert entry.measured.kpoint_density == 100
    assert entry.verification.converged_at == 100


def test_converged_entry_measured_and_context(patched, tmp_path):
    """The staged entry's measured + context fields come from the
    chosen grid point's result.toml and this structure's prediction
    record (the sub-model's sole home, DESIGN 6.2.9); cell info
    comes from the loaded structure (Bohr^3)."""
    root = _make_workspace(
        tmp_path, [50, 100, 200], [0.5, 0.5, 0.5],
        gaps=[5.0, 5.0, 5.0], kinds=["indirect"] * 3,
        mags=[0.0, 0.0, 0.0], scf_threshold=1.0)
    gh.harvest_flight(root, str(tmp_path / "db"), _DATASPACE)
    entry = patched["entries"][0]
    assert entry.measured.gap_ev == 5.0
    assert entry.measured.gap_kind == "indirect"
    assert entry.measured.total_magnetization == 0.0
    # spin_polarization is not measured -> honest 0.0 placeholder.
    assert entry.measured.spin_polarization == 0.0
    assert entry.context.basis == "fb"
    assert entry.context.functional == "gga-pbe"
    assert entry.context.kpoint_integration == "gaussian-0.1"
    assert entry.context.scf_threshold == 1.0
    assert entry.context.cell_atom_count == 8
    # 100 Angstrom^3 -> Bohr^3 via the module's own factor.
    assert entry.context.cell_volume_per_formula_unit == \
        pytest.approx(100.0 * gh._ANGSTROM3_TO_BOHR3)


def test_converged_entry_records_the_mesh(patched, tmp_path):
    """build_entry records the chosen rung's resolved axial counts
    as verification.converged_mesh, read from its result.toml
    (DESIGN 3.12.4 / 7.2)."""
    root = _make_workspace(
        tmp_path, [25, 50, 100], [0.5, 0.5, 0.5],
        meshes=[[3, 3, 3], [4, 4, 4], [5, 5, 5]])
    gh.harvest_flight(root, str(tmp_path / "db"), _DATASPACE)
    entry = patched["entries"][0]
    # Flat energies converge at the interior [4,4,4] rung.
    assert entry.verification.converged_at == 50
    assert entry.verification.converged_mesh == (4, 4, 4)


def test_converged_entry_mesh_absent_is_none(patched, tmp_path):
    """A run whose result.toml carries no kpoint_mesh (an older
    binary) records converged_mesh as None, not a failure."""
    root = _make_workspace(
        tmp_path, [25, 50, 100], [0.5, 0.5, 0.5])   # no meshes
    gh.harvest_flight(root, str(tmp_path / "db"), _DATASPACE)
    assert patched["entries"][0].verification.converged_mesh is None


def test_converged_entry_verification_and_provenance(patched,
                                                     tmp_path):
    """The verification block carries the full grid + the recovered
    predictor confidence/neighbours; provenance carries the
    flight_id (workspace basename), structure, and commit."""
    root = _make_workspace(
        tmp_path, [50, 100, 200], [0.5, 0.5, 0.5],
        confidence=0.83, neighbor_ids=("mp-1", "mp-2"))
    gh.harvest_flight(root, str(tmp_path / "db"), _DATASPACE)
    entry = patched["entries"][0]
    v = entry.verification
    assert v.grid_values == (50, 100, 200)
    assert v.grid_energies == (0.5, 0.5, 0.5)
    assert v.metric == "total_energy"
    # metric_threshold is the resolved per-atom kpoint tolerance
    #   from the prediction record, NOT the run's scf_threshold.
    assert v.metric_threshold == 5.0e-4
    assert v.predictor_confidence == 0.83
    assert v.predictor_neighbor_ids == ("mp-1", "mp-2")
    assert entry.source == "flight"
    assert entry.provenance.flight_id == "flight"   # root basename
    assert entry.provenance.source_structure == "si.skl"
    assert entry.provenance.imago_commit == "abc123"
    assert entry.provenance.curator == "guidance_harvest.py"


def test_imago_commit_falls_back_to_unknown(patched, tmp_path):
    """When result.toml carries no commit, provenance records the
    non-empty 'unknown' stand-in so the entry still satisfies the
    schema's rule-11 check (v1 convention)."""
    root = _make_workspace(tmp_path, [50, 100, 200],
                           [0.5, 0.5, 0.5])
    # Strip imago_commit from each result.toml.
    for kpd in (50, 100, 200):
        path = os.path.join(root, "wingbeats", "si",
                            f"kpt-density-{kpd}", "result.toml")
        kept = [line for line in open(path)
                if not line.startswith("imago_commit")]
        with open(path, "w") as rt:
            rt.writelines(kept)
    gh.harvest_flight(root, str(tmp_path / "db"), _DATASPACE)
    assert patched["entries"][0].provenance.imago_commit == \
        gh._UNKNOWN_COMMIT


# --------------------------------------------------------------
#  The skip paths
# --------------------------------------------------------------

def test_single_point_grid_stages_nothing(patched, tmp_path):
    """A single-point grid (trust mode OR a single-point curator
    override) harvests no entry -- one calc is weaker evidence than
    a grid, and the skip keys on len(grid)==1, not the policy
    string (DESIGN 6.2.9 / 7.7)."""
    root = _make_workspace(tmp_path, [137], [0.5],
                           policy="trust_no_verify")
    summaries = gh.harvest_flight(root, str(tmp_path / "db"),
                                  _DATASPACE)
    assert patched["entries"] == []
    assert "single point" in summaries[0]


def test_non_converged_sweep_skips_and_tags(patched, tmp_path):
    """A sweep whose energy is still moving at the top of the grid
    stages no entry, is reported as skipped, and tags the flight
    with prediction_mismatch (DESIGN 7.8 step 3d)."""
    root = _make_workspace(tmp_path, [50, 100, 200],
                           [3.0, 2.0, 1.0], scf_threshold=0.1)
    summaries = gh.harvest_flight(root, str(tmp_path / "db"),
                                  _DATASPACE)
    assert patched["entries"] == []
    assert "skipped" in summaries[0]
    marker = os.path.join(root, "prediction_mismatch.toml")
    assert os.path.exists(marker)
    assert "si = true" in open(marker).read()


def test_missing_scf_threshold_raises(patched, tmp_path):
    """result.toml with no scf_threshold cannot supply the
    convergence metric_threshold, so the harvest aborts clearly."""
    root = _make_workspace(tmp_path, [50, 100, 200],
                           [0.5, 0.5, 0.5],
                           write_scf_threshold=False)
    with pytest.raises(ValueError):
        gh.harvest_flight(root, str(tmp_path / "db"), _DATASPACE)


def test_fingerprint_loen_unit_is_ignored(patched, tmp_path):
    """A structure-only 'fingerprint' loen unit sharing the
    structure id is NOT a grid point: the convergence harvest filters
    it out by kind (DESIGN 6.2.9 / 7.8 step 2) and still stages
    exactly one entry from the convergence sweep.  (The loen unit has
    no result.toml, so a failure to skip it would raise.)"""
    root = _make_workspace(tmp_path, [50, 100, 200],
                           [0.5, 0.5, 0.5], add_loen=True)
    summaries = gh.harvest_flight(root, str(tmp_path / "db"),
                                  _DATASPACE)
    assert len(patched["entries"]) == 1
    assert patched["entries"][0].measured.kpoint_density == 100
    # one summary line for the single structure (the loen unit did
    #   not create a second group).
    assert len(summaries) == 1 and "staged" in summaries[0]


def test_no_prediction_record_skips(patched, tmp_path):
    """A structure with no [flight.predictions.<id>] record is not
    guidance-harvestable -- the record is the sole source of
    system_type and the sub-model -- so it is skipped, not staged
    (DESIGN 6.2.9 / 7.8 step 3)."""
    root = _make_workspace(tmp_path, [50, 100, 200],
                           [0.5, 0.5, 0.5],
                           write_prediction=False)
    summaries = gh.harvest_flight(root, str(tmp_path / "db"),
                                  _DATASPACE)
    assert patched["entries"] == []
    assert "no prediction record" in summaries[0]


def test_missing_gap_raises(patched, tmp_path):
    """A converged run whose result.toml omits the (required)
    electronic-character gap fails loudly rather than fabricating a
    value (DESIGN 7.8 / group-A #8).  Unlike total_magnetization, an
    absent gap means the run never surfaced what the entry exists to
    record."""
    root = _make_workspace(tmp_path, [50, 100, 200],
                           [0.5, 0.5, 0.5], write_gap=False)
    with pytest.raises(ValueError):
        gh.harvest_flight(root, str(tmp_path / "db"), _DATASPACE)


# --------------------------------------------------------------
#  Non-crystalline path
# --------------------------------------------------------------

def test_non_crystalline_signature_uses_prediction_type(patched,
                                                        tmp_path):
    """system_type rides on this structure's
    [flight.predictions.<id>] record, so a molecular flight signs
    (and stages) under 'molecular'."""
    root = _make_workspace(tmp_path, [25, 50, 100],
                           [0.5, 0.5, 0.5],
                           system_type="molecular")
    gh.harvest_flight(root, str(tmp_path / "db"), _DATASPACE)
    assert patched["system_types"] == ["molecular"]
    assert patched["entries"][0].signature.system_type == \
        "molecular"


# --------------------------------------------------------------
#  No-sweep guard
# --------------------------------------------------------------

def test_flight_without_sweep_raises(patched, tmp_path):
    """A hand-built flight with no [flight.sweep] cannot be
    harvested into the k-density dataspace."""
    root = str(tmp_path / "flight")
    serialize_flight(Flight(
        root=root,
        units=[CalcUnit(id="si", structure="si.skl",
                        calc=("kpt-density-50",))]))
    with pytest.raises(ValueError):
        gh.harvest_flight(root, str(tmp_path / "db"), _DATASPACE)
