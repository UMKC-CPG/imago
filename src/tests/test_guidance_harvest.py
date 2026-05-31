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
* ``_load_structure`` -> a tiny stand-in carrying ``num_atoms``
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
                    gaps=None, kinds=None, mags=None,
                    scf_threshold=1.0, write_scf_threshold=True,
                    policy="verify_around_prediction",
                    system_type="crystalline", confidence=0.9,
                    neighbor_ids=("mp-1", "mp-2"),
                    unit_id="si", structure="si.skl",
                    write_prediction=True):
    """Lay out a one-structure flight workspace under tmp_path: a
    flight.toml (via the real serialize_flight) plus one
    result.toml per grid point.  Returns the workspace root.  The
    knobs let each test shape the grid, the energies, and the
    [flight.prediction] block."""
    root = str(tmp_path / "flight")
    units = [CalcUnit(id=unit_id, structure=structure,
                      calc=(f"kpt-density-{k}",)) for k in kpds]
    sweep = SweepRecord(
        varied_axes=("kpt-density",),
        fixed_axes={"basis": "fb", "functional": "gga-pbe",
                    "kpoint_integration": "gaussian-0.1"})
    metadata = {}
    if write_prediction:
        metadata["prediction"] = {
            "policy": policy, "system_type": system_type,
            "confidence": confidence,
            "neighbor_entry_ids": list(neighbor_ids),
            "predicted_value": 100.0, "is_under_trained": False}
    serialize_flight(Flight(root=root, units=units, sweep=sweep,
                            metadata=metadata))

    for index, kpd in enumerate(kpds):
        run_dir = os.path.join(root, "wingbeats", unit_id,
                               f"kpt-density-{kpd}")
        os.makedirs(run_dir, exist_ok=True)
        gap = gaps[index] if gaps is not None else 1.5
        kind = kinds[index] if kinds is not None else "indirect"
        mag = mags[index] if mags is not None else 0.0
        with open(os.path.join(run_dir, "result.toml"), "w") as rt:
            rt.write(toml_line("total_energy", energies[index]))
            rt.write(toml_line("gap_ev", gap))
            rt.write(toml_line("gap_kind", kind))
            rt.write(toml_line("total_magnetization", mag))
            if write_scf_threshold:
                rt.write(toml_line("scf_threshold", scf_threshold))
            rt.write(toml_line("imago_commit", "abc123"))
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
    monkeypatch.setattr(gh, "_load_structure", fake_load_structure)
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
    """pick_converged needs BOTH neighbour deltas below threshold
    and never returns an endpoint."""
    # Flat only on the high side of index 1 -> not converged there;
    #   index 2 is flat on both sides -> chosen.
    energies = [3.0, 2.0, 1.99, 1.985]
    assert gh.pick_converged(energies, 0.1) == 2
    # Still moving everywhere -> None.
    assert gh.pick_converged([3.0, 2.0, 1.0], 0.1) is None


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
    chosen grid point's result.toml and the sweep's fixed axes;
    cell info comes from the loaded structure (Bohr^3)."""
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
    assert v.metric_threshold == 1.0          # == scf_threshold
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

def test_trust_mode_stages_nothing(patched, tmp_path):
    """A trust-mode flight (single point, policy trust_no_verify)
    harvests no entry -- one calc is weaker evidence than a grid
    (DESIGN 7.7)."""
    root = _make_workspace(tmp_path, [137], [0.5],
                           policy="trust_no_verify")
    summaries = gh.harvest_flight(root, str(tmp_path / "db"),
                                  _DATASPACE)
    assert patched["entries"] == []
    assert "trusted" in summaries[0]


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


# --------------------------------------------------------------
#  Non-crystalline path
# --------------------------------------------------------------

def test_non_crystalline_signature_uses_prediction_type(patched,
                                                        tmp_path):
    """system_type rides on the [flight.prediction] block, so a
    molecular flight signs (and stages) under 'molecular'."""
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
