"""Tests for the predict-then-verify flight builder
(kaleidoscope.builders.kpoint_convergence; DESIGN 6.2.8 / 7.7;
PSEUDOCODE 15.6).

The builder's job is purely to turn a prediction into a
verification-grid Flight: choose the grid of k-densities and the
policy, tag each grid point, and assemble the PredictionRecord.
The prediction itself (compute_signature + predict) belongs to
guidance_db and is tested there, so here it is monkeypatched to a
caller-chosen PredictionResult.  That isolates the grid/policy/tag
logic and lets every test run with no $IMAGO_DATA and no real
dataspace: a non-path structure object means StructureControl is
never instantiated.
"""

import math
import tomllib
import types

import pytest

from kaleidoscope import CalcUnit, Flight, KaleidoscopeError
from kaleidoscope.builders import kpoint_convergence as kc
from kaleidoscope.workspace import serialize_flight
from guidance_db import PredictionResult, Signature


# --------------------------------------------------------------
#  Fixtures / helpers
# --------------------------------------------------------------

# A non-path structure: _load_structure returns it unchanged, so
#   no StructureControl (hence no elements.dat) is ever touched.
_STRUCTURE = object()

# Tool-facing run settings (dest-keyed, coded): copied verbatim
#   into every unit and never inspected by the builder.
_OPTIONS = {"scf_basis": "fb", "xccode": 200,
            "scfkpint": 0, "converg": 1.0e-6}

# The human sub-model the predictor + record read, in its own dict
#   (DESIGN 6.2.8): kept separate from the tool-facing _OPTIONS.
_SUBMODEL = {"basis": "fb", "functional": "gga-pbe",
             "kpoint_integration": "gaussian-0.1"}

_DATASPACE = types.SimpleNamespace(group_table={})


def _signature():
    """A stand-in crystalline Signature (the patched
    compute_signature returns this; its exact values do not matter
    to the builder)."""
    return Signature(
        system_type="crystalline",
        composition_vector=(1.0,) + (0.0,) * 12,
        lattice_family="cubic",
        lattice_onehot=(1.0,) + (0.0,) * 5)


def _result(kpd=100.0, confidence=0.9, under_trained=False,
            neighbor_ids=("mp-1", "mp-2"), gap=1.2,
            magnetization=0.0):
    """Build a PredictionResult for the patched predictor."""
    return PredictionResult(
        predicted_kpoint_density=kpd, confidence=confidence,
        is_under_trained=under_trained,
        neighbor_entry_ids=neighbor_ids,
        predicted_gap=gap, predicted_magnetization=magnetization)


def _no_predict(*args, **kwargs):
    """A stand-in predictor that fails if called -- used by the
    curator-override tests to prove the predictor is bypassed."""
    raise AssertionError(
        "predict() must not run in curator-override mode")


@pytest.fixture
def patched(monkeypatch):
    """Return a function that patches the physics layer to a fixed
    signature and a caller-supplied PredictionResult."""
    def _apply(result):
        monkeypatch.setattr(kc, "compute_signature",
                            lambda *args, **kw: _signature())
        monkeypatch.setattr(kc, "predict",
                            lambda *args, **kw: result)
    return _apply


@pytest.fixture
def patched_no_predict(monkeypatch):
    """Patch only compute_signature; install a predictor that
    raises if consulted (for the override-bypass tests)."""
    monkeypatch.setattr(kc, "compute_signature",
                        lambda *args, **kw: _signature())
    monkeypatch.setattr(kc, "predict", _no_predict)


def _build(**kwargs):
    """Run build_kpoint_convergence with the standard stand-in
    inputs and a caller-chosen id (a non-path structure cannot
    derive one).  The predictor itself is already installed by the
    ``patched`` fixture, so only the builder's own arguments are
    passed here."""
    return kc.build_kpoint_convergence(
        _STRUCTURE, _OPTIONS, _DATASPACE, "crystalline", _SUBMODEL,
        id="si", **kwargs)


# --------------------------------------------------------------
#  Pure grid / tag helpers (no patching needed)
# --------------------------------------------------------------

def test_logspace_endpoints_and_geometric_center():
    """logspace includes both endpoints and is geometrically
    spaced; a single-point request is the geometric midpoint."""
    grid = kc.logspace(10.0, 40.0, 3)
    assert grid[0] == pytest.approx(10.0)
    assert grid[-1] == pytest.approx(40.0)
    assert grid[1] == pytest.approx(20.0)        # sqrt(10*40)
    assert kc.logspace(10.0, 40.0, 1) == \
        pytest.approx([math.sqrt(400.0)])


def test_build_verification_grid_tight_when_confident():
    """confidence=1 -> a tight 3-point span [c/1.2, c*1.2]."""
    grid = kc.build_verification_grid(100.0, 1.0)
    assert len(grid) == 3
    assert grid[0] == pytest.approx(100.0 / 1.2)
    assert grid[-1] == pytest.approx(100.0 * 1.2)


def test_build_verification_grid_wide_when_unsure():
    """confidence=0.3 -> ~6 points over a wider [c/2.25, c*2.25]
    span (width = 1.2 + 1.5*0.7, n = round(3 + 4*0.7))."""
    grid = kc.build_verification_grid(100.0, 0.3)
    assert len(grid) == 6
    assert grid[0] == pytest.approx(100.0 / 2.25)
    assert grid[-1] == pytest.approx(100.0 * 2.25)


def test_encode_axis_value_examples():
    """Integer-valued floats render as plain integers; a decimal
    uses 'p' for '.' and a leading 'm' for a negative
    (DESIGN 6.2.4 rule 3)."""
    assert kc.encode_axis_value(50.0) == "50"
    assert kc.encode_axis_value(1.5) == "1p5"
    assert kc.encode_axis_value(-2.0) == "m2"
    assert kc.encode_axis_value(0.1) == "0p1"


def test_build_calc_tag_examples_and_order():
    """build_calc_tag returns one '<axis>-<value>' component per
    axis, in mapping order (DESIGN 6.2.4)."""
    assert kc.build_calc_tag({"kpt-density": 50}) == \
        ("kpt-density-50",)
    assert kc.build_calc_tag(
        {"kpt-density": 50, "basis-size": 3}) == \
        ("kpt-density-50", "basis-size-3")


def test_build_calc_tag_rejects_non_slug_axis():
    """An axis name that is not a slug aborts -- it would be an
    unsafe directory level."""
    with pytest.raises(KaleidoscopeError):
        kc.build_calc_tag({"Bad Axis": 5})


# --------------------------------------------------------------
#  build_kpoint_convergence: grid + policy selection
# --------------------------------------------------------------

def test_high_confidence_gives_tight_three_point_grid(patched):
    """A confident prediction verifies with a tight 3-point grid
    centered on the predicted value."""
    patched(_result(kpd=100.0, confidence=1.0))
    flight, record = _build(verify=True)
    kpds = [unit.options["kpd"] for unit in flight.units]
    assert len(kpds) == 3
    assert kpds == sorted(kpds)
    assert 100 in kpds                       # center is preserved
    assert record.policy == "verify_around_prediction"


def test_low_confidence_gives_wider_six_point_grid(patched):
    """A low-confidence prediction widens to a 6-point grid."""
    patched(_result(kpd=100.0, confidence=0.3))
    flight, record = _build(verify=True)
    assert len(flight.units) == 6
    assert record.policy == "verify_around_prediction"


def test_under_trained_falls_back_to_wide_grid(patched):
    """When the predictor is under-trained the builder ignores the
    (meaningless) predicted value and lays out the fixed wide-grid
    default with the no-prior policy (DESIGN 7.9)."""
    patched(_result(kpd=0.0, confidence=0.0, under_trained=True))
    flight, record = _build(verify=True)
    kpds = [unit.options["kpd"] for unit in flight.units]
    assert kpds == [25, 50, 100, 150, 200, 250, 300, 400]
    assert record.policy == "wide_grid_no_prior"


def test_trust_mode_gives_single_point(patched):
    """Trust mode (verify=False) collapses to one unit at the
    predicted value, regardless of confidence."""
    patched(_result(kpd=137.0, confidence=0.4))
    flight, record = _build(verify=False)
    assert len(flight.units) == 1
    assert flight.units[0].options["kpd"] == 137
    assert record.policy == "trust_no_verify"


# --------------------------------------------------------------
#  build_kpoint_convergence: curator override (center given)
# --------------------------------------------------------------

def test_curator_override_pins_grid_and_bypasses_predictor(
        patched_no_predict):
    """A pinned ``center`` (the 5.7 kpoint_spec override) lays out a
    tight verify grid around the curator value WITHOUT consulting
    the predictor; policy is curator_override (DESIGN 6.2.9)."""
    flight, record = kc.build_kpoint_convergence(
        _STRUCTURE, _OPTIONS, _DATASPACE, "crystalline", _SUBMODEL,
        id="si", center=120.0)
    kpds = [unit.options["kpd"] for unit in flight.units]
    assert record.policy == "curator_override"
    assert record.predicted_kpoint_density == 120.0
    assert record.confidence == 1.0
    assert record.is_under_trained is False
    assert record.neighbor_entry_ids == ()
    assert record.predicted_gap is None
    assert len(kpds) == 3                     # tight verify grid
    assert 120 in kpds                        # centred on the pin
    # The sub-model is still recorded on the override record.
    assert record.basis == "fb"
    assert record.kpoint_integration == "gaussian-0.1"


def test_curator_override_single_point_when_not_verifying(
        patched_no_predict):
    """center with verify=False yields a single pinned point and
    still bypasses the predictor."""
    flight, record = kc.build_kpoint_convergence(
        _STRUCTURE, _OPTIONS, _DATASPACE, "crystalline", _SUBMODEL,
        id="si", center=120.0, verify=False)
    assert [u.options["kpd"] for u in flight.units] == [120]
    assert record.policy == "curator_override"


# --------------------------------------------------------------
#  build_kpoint_convergence: units and the prediction record
# --------------------------------------------------------------

def test_units_carry_kpd_tag_options_and_key_fields(patched):
    """Each grid unit pins its k-density in both the makeinput
    'kpd' option and the calc tag, carries the structure and the
    imago wingbeat, and gets cache key fields built from options."""
    patched(_result(kpd=100.0, confidence=1.0))
    flight, _ = _build(verify=True)
    for unit in flight.units:
        kpd = unit.options["kpd"]
        assert unit.calc == (f"kpt-density-{kpd}",)
        assert unit.structure is _STRUCTURE
        assert unit.wingbeat == "imago"
        assert unit.id == "si"
        # converg (the SCF limit) flows into the cache scalars;
        #   the structure is the single byte-compared key file.
        assert unit.key_fields.scalars["converg"] == 1.0e-6
        assert [f.name for f in unit.key_fields.files] == \
            ["structure"]


def test_sweep_record_varied_axis_and_empty_fixed_axes(patched):
    """The flight's SweepRecord names the swept axis; fixed_axes is
    EMPTY because the (basis, functional, kpoint_integration)
    sub-model rides on the per-structure record, never duplicated
    here (DESIGN 6.2.9)."""
    patched(_result(kpd=100.0, confidence=1.0))
    flight, _ = _build(verify=True)
    assert flight.sweep.varied_axes == ("kpt-density",)
    assert flight.sweep.fixed_axes == {}


def test_prediction_record_shape_and_stash(patched):
    """The returned record carries the full provenance -- including
    the sub-model (the sole home now that fixed_axes is empty) --
    and the same data is stashed (as a plain dict) under
    flight.metadata["predictions"][id] so the harvest recovers it
    per structure (PSEUDOCODE 15.6 / DESIGN 6.2.9)."""
    patched(_result(kpd=100.0, confidence=0.83,
                    neighbor_ids=("mp-1", "mp-2"), gap=1.2,
                    magnetization=0.0))
    flight, record = _build(verify=True)
    assert record.predicted_kpoint_density == 100.0
    assert record.confidence == 0.83
    assert record.is_under_trained is False
    assert record.neighbor_entry_ids == ("mp-1", "mp-2")
    assert record.predicted_gap == 1.2
    assert record.system_type == "crystalline"
    # The sub-model rides on the record (DESIGN 6.2.9).
    assert record.basis == "fb"
    assert record.functional == "gga-pbe"
    assert record.kpoint_integration == "gaussian-0.1"
    # Keyed by structure id so one flight can carry many.
    stashed = flight.metadata["predictions"]["si"]
    assert stashed["policy"] == "verify_around_prediction"
    # asdict keeps tuples in memory; they become TOML arrays only
    #   once serialize_flight writes them (see the round-trip test).
    assert stashed["neighbor_entry_ids"] == ("mp-1", "mp-2")
    assert stashed["kpoint_integration"] == "gaussian-0.1"
    # feature_vector flattened to a nested dict (the Signature).
    assert stashed["feature_vector"]["lattice_family"] == "cubic"


# --------------------------------------------------------------
#  Error paths and the Part-A serialization round-trip
# --------------------------------------------------------------

def test_missing_required_submodel_key_raises(patched):
    """A submodel dict missing a sub-model-selecting key aborts."""
    patched(_result())
    bad = {"functional": "gga-pbe",
           "kpoint_integration": "gaussian-0.1"}   # no basis
    with pytest.raises(KaleidoscopeError):
        kc.build_kpoint_convergence(
            _STRUCTURE, _OPTIONS, _DATASPACE, "crystalline", bad,
            id="si")


def test_non_path_structure_without_id_raises(patched):
    """A non-path structure with no explicit id cannot yield a
    stable slug, so the builder refuses rather than guess."""
    patched(_result())
    with pytest.raises(KaleidoscopeError):
        kc.build_kpoint_convergence(
            _STRUCTURE, _OPTIONS, _DATASPACE, "crystalline",
            _SUBMODEL)


def test_flight_serializes_prediction_and_sweep(patched, tmp_path):
    """End to end with Part A: a built flight serializes its
    [flight.sweep] and per-id [flight.predictions.<id>] blocks, and
    they read back through tomllib unchanged."""
    patched(_result(kpd=100.0, confidence=1.0))
    flight, _ = _build(verify=True, root=str(tmp_path))
    serialize_flight(flight)
    with open(tmp_path / "flight.toml", "rb") as flight_file:
        data = tomllib.load(flight_file)
    assert data["flight"]["sweep"]["varied_axes"] == ["kpt-density"]
    pred = data["flight"]["predictions"]["si"]
    assert pred["policy"] == "verify_around_prediction"
    assert pred["feature_vector"]["lattice_family"] == "cubic"
