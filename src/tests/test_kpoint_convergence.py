## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

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


# --------------------------------------------------------------
#  Pure calc-tag helpers (no patching needed)
# --------------------------------------------------------------

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
#  The mesh-dispatch split (DESIGN 7.7; PSEUDOCODE 4e.7)
# --------------------------------------------------------------

def test_encode_mesh_value_joins_counts_with_hyphens():
    """A mesh renders as its three counts hyphen-joined -- slug-safe
    because axial counts are positive integers."""
    assert kc.encode_mesh_value([4, 4, 4]) == "4-4-4"
    assert kc.encode_mesh_value([10, 10, 4]) == "10-10-4"
    assert kc.encode_mesh_value([3, 3, 2]) == "3-3-2"


def test_decode_mesh_value_inverts_encode():
    """decode_mesh_value round-trips encode_mesh_value back to the
    integer count triple."""
    for mesh in ([4, 4, 4], [10, 10, 4], [3, 3, 2]):
        assert kc.decode_mesh_value(kc.encode_mesh_value(mesh)) \
            == mesh


def test_build_mesh_unit_sets_scfkp_and_kpt_mesh_tag():
    """An explicit-mesh unit carries the mesh as the makeinput
    `scfkp` option and tags itself `kpt-mesh-<a>-<b>-<c>`, keeping
    the caller's options dict untouched."""
    unit = kc.build_mesh_unit(_STRUCTURE, _OPTIONS, [4, 4, 4], "si")
    assert unit.id == "si"
    assert unit.wingbeat == "imago"
    assert unit.calc == ("kpt-mesh-4-4-4",)
    assert unit.options["scfkp"] == [4, 4, 4]
    # The shared options dict must not have gained a mesh.
    assert "scfkp" not in _OPTIONS


def test_build_mesh_unit_tags_an_anisotropic_mesh():
    """A non-cubic mesh tags each axis count in order."""
    unit = kc.build_mesh_unit(_STRUCTURE, _OPTIONS, [5, 5, 2], "gr")
    assert unit.calc == ("kpt-mesh-5-5-2",)
    assert unit.options["scfkp"] == [5, 5, 2]


def test_key_fields_distinguish_the_integration_scheme():
    """The cache identity byte-compares TWO makeinput outputs, and
    the second is not optional polish (DESIGN 6.2.5).

    ``structure.dat`` bakes in the type/species assignment, basis,
    functional, and potential -- but NOT the k-point integration
    scheme, which reaches ``kp-scf.dat`` as KPOINT_INTG_CODE.  Keyed
    on ``structure.dat`` alone, one solid at one mesh under two
    schemes shares a run directory and HITS, returning the other
    scheme's answer under the name of the one asked for.  That is
    wrong physics reported silently, not merely a stale result."""
    fields = kc.standard_key_fields(_STRUCTURE, _OPTIONS)
    assert [key_file.name for key_file in fields.files] == [
        "structure.dat", "kp-scf.dat"]


def test_key_scalars_stay_the_convergence_limit_alone():
    """The scheme is a key FILE, deliberately not a key scalar.  The
    scalars are compared as a whole table, so a name no stored
    ``cache_key.toml`` carries would invalidate every cached unit in
    every surviving workspace at once -- a mass false miss, which is
    the failure the cache design works hardest to avoid.  A key file
    costs nothing, because every run directory already stages it."""
    fields = kc.standard_key_fields(
        _STRUCTURE, dict(_OPTIONS, converg=1.0e-6, scfkpint=1))
    assert fields.scalars == {"converg": 1.0e-6}


def test_key_file_sources_start_provisional():
    """Every key file's source starts at the skeleton and is
    re-pointed by the driver's prepare step once makeinput has
    actually written the files (DESIGN 6.2.5, Model A).  Pinned for
    both files, since a second file added without the matching
    re-point would name a path that never exists."""
    fields = kc.standard_key_fields("/cache/si-prim.skl", _OPTIONS)
    assert {key_file.source for key_file in fields.files} == {
        "/cache/si-prim.skl"}


def test_predict_kpoint_density_returns_prediction(patched):
    """The predict-only builder returns the predicted density,
    confidence, under-trained flag, and a `predict_then_climb`
    record -- and lays no grid."""
    patched(_result(kpd=150.0, confidence=0.8, under_trained=False))
    density, confidence, under_trained, record = \
        kc.predict_kpoint_density(
            _STRUCTURE, _DATASPACE, "crystalline", _SUBMODEL)
    assert density == 150.0
    assert confidence == 0.8
    assert under_trained is False
    assert record.policy == "predict_then_climb"
    assert record.predicted_kpoint_density == 150.0
    assert record.basis == "fb"


def test_predict_kpoint_density_curator_override(patched_no_predict):
    """A pinned `center` bypasses the predictor: the density is the
    pinned value at full confidence with a `curator_override`
    record."""
    density, confidence, under_trained, record = \
        kc.predict_kpoint_density(
            _STRUCTURE, _DATASPACE, "crystalline", _SUBMODEL,
            center=250.0)
    assert density == 250.0
    assert confidence == 1.0
    assert under_trained is False
    assert record.policy == "curator_override"
    assert record.predicted_kpoint_density == 250.0


def test_predict_kpoint_density_requires_full_submodel(patched):
    """A submodel missing one of the three sub-model-selecting names
    is rejected, exactly as the density builder rejects it."""
    patched(_result())
    with pytest.raises(KaleidoscopeError):
        kc.predict_kpoint_density(
            _STRUCTURE, _DATASPACE, "crystalline",
            {"basis": "fb", "functional": "gga-pbe"})
