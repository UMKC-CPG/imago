"""Tests for the curation-manifest schema library's writer.

These exercise :func:`curation_manifest.format_manifest` and its
file-producing counterpart :func:`write_manifest`, with the central
guarantee being a round-trip: a manifest written by the emitter reads
back through :func:`load_manifest_v2` to the same reference solids.
The readers themselves are covered alongside the producer in
test_build_initial_potentials.py.
"""

import pytest

from curation_manifest import (
    CurationManifest,
    ReferenceSolid,
    ReferenceEntry,
    ManifestFingerprint,
    format_manifest,
    write_manifest,
    load_manifest_v2,
)

pytestmark = pytest.mark.unit


def _si_solid(**overrides) -> ReferenceSolid:
    """A valid single-element Si reference solid (cod-sourced) with
    one default entry carrying a preferred reduce and a preferred
    bispectrum fingerprint -- enough to satisfy rules 7, 10, 11."""

    base = dict(
        reference_id="si_diamond",
        system_type="crystalline",
        basis="fb",
        functional="wigner",
        kpoint_integration="gaussian",
        kpoint_spec={"shift": [0.0, 0.0, 0.0]},
        scf_threshold=1.0e-6,
        cod_id=2104737,
        cod_revision="201401",
        structure_path=None,
        entries=[ReferenceEntry(
            element="Si", atom_site=1, default=True,
            description="Si diamond, site 1", label=None,
            fingerprints=[
                ManifestFingerprint(
                    method="reduce",
                    sub_spec={"level": 2, "thick": 0.5,
                              "cutoff": 5.0, "tolerance": 0.05},
                    preferred=True),
                ManifestFingerprint(
                    method="bispectrum",
                    sub_spec={"twoj1": 8, "twoj2": 8,
                              "cutoff": 9.0},
                    preferred=True)])])
    base.update(overrides)
    return ReferenceSolid(**base)


def _manifest(*solids) -> CurationManifest:
    return CurationManifest(
        schema_version=2, manifest_path="(in memory)",
        reference_solids=list(solids))


def _write_and_load(manifest, tmp_path):
    """Write a manifest to a temp file and read it back."""
    path = tmp_path / "manifest.toml"
    write_manifest(manifest, str(path))
    return load_manifest_v2(str(path))


# ==============================================================
#  Round-trip: written manifest reads back to the same solids
# ==============================================================

def test_cod_solid_round_trips(tmp_path):
    """A cod-sourced solid with entries and fingerprints survives a
    format -> write -> load_manifest_v2 round-trip unchanged."""
    original = _si_solid()
    loaded = _write_and_load(_manifest(original), tmp_path)
    assert loaded.schema_version == 2
    assert loaded.reference_solids == [original]


def test_structure_path_solid_round_trips(tmp_path):
    """A structure_path solid round-trips too.  The path must
    resolve to an existing file (rule 4), so the skeleton is
    created beside the manifest before loading."""
    (tmp_path / "a_si.skl").write_text("# placeholder skeleton\n")
    original = _si_solid(
        reference_id="a_si", system_type="amorphous",
        cod_id=None, cod_revision=None,
        structure_path="a_si.skl")
    loaded = _write_and_load(_manifest(original), tmp_path)
    assert loaded.reference_solids == [original]


def test_multiple_solids_round_trip(tmp_path):
    """Two reference solids in one manifest both survive.  Only the
    first carries the element's single default entry and its
    preferred fingerprints (rules 7, 10); the second contributes a
    non-default, non-fingerprinted Si site."""
    first = _si_solid()
    second = _si_solid(
        reference_id="si_bc8", cod_id=4350826,
        cod_revision="297572",
        entries=[ReferenceEntry(
            element="Si", atom_site=1, default=False,
            description="Si BC8, site 1", label=None,
            fingerprints=[])])
    loaded = _write_and_load(_manifest(first, second), tmp_path)
    assert loaded.reference_solids == [first, second]


# ==============================================================
#  Layout choices the emitter makes
# ==============================================================

def test_floats_use_readable_repr_not_padded_scientific():
    """Manifest floats read naturally -- shortest round-trippable
    repr, not the database writer's 17-digit padded form."""
    text = format_manifest(_manifest(_si_solid()))
    assert "scf_threshold = 1e-06" in text
    assert "thick = 0.5" in text
    assert "tolerance = 0.05" in text
    # The padded binary64 form must NOT appear.
    assert "e-06" not in text.replace("1e-06", "")


def test_inline_table_keeps_authored_key_order():
    """A sub_spec keeps the order it was authored in, rather than
    being alphabetised, so it reads in its natural sequence."""
    text = format_manifest(_manifest(_si_solid()))
    assert ("sub_spec = { level = 2, thick = 0.5, "
            "cutoff = 5.0, tolerance = 0.05 }") in text


def test_label_omitted_when_absent_emitted_when_present(tmp_path):
    """An absent label is left out (derived at harvest); an explicit
    label is emitted and round-trips."""
    without = format_manifest(_manifest(_si_solid()))
    assert "label =" not in without

    entry = ReferenceEntry(
        element="Si", atom_site=1, default=True,
        description="Si diamond", label="si_diamond-Si1",
        fingerprints=[])
    labelled = _si_solid(entries=[entry])
    text = format_manifest(_manifest(labelled))
    assert 'label = "si_diamond-Si1"' in text
    loaded = _write_and_load(_manifest(labelled), tmp_path)
    assert loaded.reference_solids[0].entries[0].label == \
        "si_diamond-Si1"


def test_preferred_emitted_only_when_true():
    """preferred = true is emitted for a preferred record; a
    non-preferred record omits the flag (the reader defaults it to
    false)."""
    entry = ReferenceEntry(
        element="Si", atom_site=1, default=True,
        description="Si diamond", label=None,
        fingerprints=[
            ManifestFingerprint(
                method="reduce",
                sub_spec={"level": 2, "thick": 0.5,
                          "cutoff": 5.0, "tolerance": 0.05},
                preferred=True),
            ManifestFingerprint(
                method="reduce",
                sub_spec={"level": 3, "thick": 0.5,
                          "cutoff": 5.0, "tolerance": 0.05},
                preferred=False)])
    text = format_manifest(_manifest(_si_solid(entries=[entry])))
    assert text.count("preferred = true") == 1
    assert "preferred = false" not in text


def test_empty_kpoint_spec_renders_braces(tmp_path):
    """An empty kpoint_spec renders as {} and round-trips (the
    density is then left to predict-then-verify)."""
    solid = _si_solid(kpoint_spec={})
    text = format_manifest(_manifest(solid))
    assert "kpoint_spec = {}" in text
    loaded = _write_and_load(_manifest(solid), tmp_path)
    assert loaded.reference_solids[0].kpoint_spec == {}


def test_description_with_special_characters_round_trips(tmp_path):
    """A description carrying quotes and a tab is escaped and reads
    back byte-for-byte."""
    entry = ReferenceEntry(
        element="Si", atom_site=1, default=True,
        description='Si "diamond"\tphase', label=None,
        fingerprints=[])
    solid = _si_solid(entries=[entry])
    loaded = _write_and_load(_manifest(solid), tmp_path)
    assert loaded.reference_solids[0].entries[0].description == \
        'Si "diamond"\tphase'
