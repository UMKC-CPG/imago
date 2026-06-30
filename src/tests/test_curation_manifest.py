"""Tests for the curation-manifest schema library's writer.

These exercise :func:`curation_manifest.format_manifest` and its
file-producing counterpart :func:`write_manifest`, with the central
guarantee being a round-trip: a manifest written by the emitter reads
back through :func:`load_manifest_v2` to the same manifest -- the
database-wide ``[characterization]`` recipe and the reference solids
both survive unchanged.  The readers themselves are covered alongside
the producer in test_build_initial_potentials.py.
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
    default_run_settings,
    default_characterization,
    resolve_run_settings,
)

pytestmark = pytest.mark.unit


def _characterization() -> list[ManifestFingerprint]:
    """The database-wide preferred recipe: one reduce and one
    bispectrum declaration, each the family's single preferred
    record (DESIGN 5.7).  The reader stamps ``preferred=True`` on
    every ``[characterization]`` record, so that is the shape a
    round-trip recovers."""

    return [
        ManifestFingerprint(
            method="reduce",
            sub_spec={"level": 2, "thick": 0.5,
                      "cutoff": 5.0, "tolerance": 0.05},
            preferred=True),
        ManifestFingerprint(
            method="bispectrum",
            sub_spec={"twoj1": 8, "twoj2": 8, "cutoff": 9.0},
            preferred=True)]


def _si_solid(**overrides) -> ReferenceSolid:
    """A valid single-element Si reference solid (cod-sourced) with
    one default customization.  The customization carries a single
    NON-preferred per-entry fingerprint override (a coarser
    bispectrum), exercising the per-entry path; the preferred recipe
    lives in the manifest's [characterization] block, not here."""

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
                    method="bispectrum",
                    sub_spec={"twoj1": 6, "twoj2": 4},
                    preferred=False)])])
    base.update(overrides)
    return ReferenceSolid(**base)


def _manifest(*solids, characterization=None,
              defaults=None) -> CurationManifest:
    """Wrap reference solids in a manifest, defaulting to the
    standard preferred recipe so every round-trip exercises the
    [characterization] block unless a test overrides it.
    ``defaults`` populates the optional [defaults] run settings."""

    return CurationManifest(
        schema_version=2, manifest_path="(in memory)",
        characterization=(characterization
                          if characterization is not None
                          else _characterization()),
        defaults=defaults if defaults is not None else {},
        reference_solids=list(solids))


def _write_and_load(manifest, tmp_path):
    """Write a manifest to a temp file and read it back."""
    path = tmp_path / "manifest.toml"
    write_manifest(manifest, str(path))
    return load_manifest_v2(str(path))


# ==============================================================
#  Round-trip: written manifest reads back to the same manifest
# ==============================================================

def test_cod_solid_round_trips(tmp_path):
    """A cod-sourced solid with a customization and a per-entry
    fingerprint override survives a format -> write ->
    load_manifest_v2 round-trip, and the database-wide preferred
    recipe round-trips alongside it."""
    original = _si_solid()
    loaded = _write_and_load(_manifest(original), tmp_path)
    assert loaded.schema_version == 2
    assert loaded.characterization == _characterization()
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
    first carries the element's single default customization (rule
    7); the second contributes a non-default, non-fingerprinted Si
    site."""
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


def test_source_description_round_trips(tmp_path):
    """A per-solid source_description is emitted when present and
    reads back verbatim; an absent one stays None and is omitted."""
    with_hint = _si_solid(
        source_description="Si in the diamond structure (Fd-3m).")
    text = format_manifest(_manifest(with_hint))
    assert ('source_description = "Si in the diamond '
            'structure (Fd-3m)."') in text
    loaded = _write_and_load(_manifest(with_hint), tmp_path)
    assert loaded.reference_solids[0].source_description == \
        "Si in the diamond structure (Fd-3m)."

    without = format_manifest(_manifest(_si_solid()))
    assert "source_description =" not in without


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
    being alphabetised, so it reads in its natural sequence.  The
    preferred reduce recipe lives in the [characterization] block."""
    text = format_manifest(_manifest(_si_solid()))
    assert ("sub_spec = { level = 2, thick = 0.5, "
            "cutoff = 5.0, tolerance = 0.05 }") in text


def test_characterization_block_emitted_without_preferred_key(
        tmp_path):
    """The preferred recipe is emitted as a [characterization]
    block, one [[characterization.fingerprint]] per method.  The
    preferred flag is structural, so it is NEVER serialized -- not
    on the characterization records and not on per-entry overrides
    -- yet the reader recovers it from the block each record lands
    in."""
    text = format_manifest(_manifest(_si_solid()))
    assert "[characterization]" in text
    assert text.count("[[characterization.fingerprint]]") == 2
    assert "preferred" not in text

    loaded = _write_and_load(_manifest(_si_solid()), tmp_path)
    # Characterization records read back preferred; the per-entry
    #   override reads back non-preferred -- both derived from the
    #   block, not from any written flag.
    assert all(fp.preferred for fp in loaded.characterization)
    entry = loaded.reference_solids[0].entries[0]
    assert all(not fp.preferred for fp in entry.fingerprints)


def test_no_characterization_block_when_recipe_empty():
    """A manifest with no preferred recipe emits no
    [characterization] block, and the reader tolerates its
    absence (an empty recipe)."""
    text = format_manifest(
        _manifest(_si_solid(), characterization=[]))
    assert "[characterization]" not in text


# ==============================================================
#  Optional customization fields: emitted only when set
# ==============================================================

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


def test_optional_customization_fields_omitted_when_absent(
        tmp_path):
    """Every customization field is optional (DESIGN 5.2.2).  An
    entry that sets none emits a bare [[reference_solid.entry]]
    header with no field lines, and reads back all-None / False."""
    bare = _si_solid(entries=[ReferenceEntry()])
    text = format_manifest(_manifest(bare))
    assert "[[reference_solid.entry]]" in text
    assert "element =" not in text
    assert "atom_site =" not in text
    assert "description =" not in text
    # default=False is the absent reading, so no default line.
    assert "default =" not in text

    loaded = _write_and_load(_manifest(bare), tmp_path)
    entry = loaded.reference_solids[0].entries[0]
    assert entry.element is None
    assert entry.atom_site is None
    assert entry.default is False
    assert entry.description is None
    assert entry.label is None


def test_default_emitted_only_when_true(tmp_path):
    """A default = true customization emits its flag and round-trips;
    a default = false one omits it (the absent flag reads as
    false)."""
    text = format_manifest(_manifest(_si_solid()))
    assert "default = true" in text
    assert "default = false" not in text


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


# ==============================================================
#  [defaults] block + run-setting resolution (DESIGN 5.7)
# ==============================================================

_DEFAULTS = {
    "basis": "fb", "functional": "wigner",
    "kpoint_integration": "gaussian",
    "kpoint_spec": {"shift": [0.0, 0.0, 0.0]},
    "scf_threshold": 1.0e-6}


def _sparse_si(**overrides) -> ReferenceSolid:
    """A Si solid that inherits every run setting (all None),
    unless an override pins one."""
    sparse = dict(basis=None, functional=None,
                  kpoint_integration=None, kpoint_spec=None,
                  scf_threshold=None)
    sparse.update(overrides)
    return _si_solid(**sparse)


def test_defaults_block_round_trips(tmp_path):
    # A [defaults] block plus a solid that inherits every run
    #   setting round-trips: the loaded solid stays sparse (None
    #   where inherited) and the [defaults] survive.
    original = _sparse_si()
    loaded = _write_and_load(
        _manifest(original, defaults=dict(_DEFAULTS)), tmp_path)
    assert loaded.defaults == _DEFAULTS
    assert loaded.reference_solids == [original]
    assert loaded.reference_solids[0].basis is None


def test_writer_omits_inherited_settings():
    # The emitted TOML carries a [defaults] table and the solid
    #   names none of the inherited run settings.
    text = format_manifest(
        _manifest(_sparse_si(), defaults=dict(_DEFAULTS)))
    assert "[defaults]" in text
    solid_block = text.split("[[reference_solid]]", 1)[1]
    for key in ("basis", "functional", "kpoint_integration",
                "kpoint_spec", "scf_threshold"):
        assert key not in solid_block


def test_per_solid_override_is_emitted(tmp_path):
    # A solid that overrides one setting keeps just that one; the
    #   rest still inherit (stay None).
    loaded = _write_and_load(
        _manifest(_sparse_si(basis="eb"),
                  defaults=dict(_DEFAULTS)), tmp_path)
    solid = loaded.reference_solids[0]
    assert solid.basis == "eb"
    assert solid.functional is None


def test_resolve_run_settings_fills_from_defaults():
    # An inherited setting resolves from [defaults]; an override
    #   survives resolution.
    resolved = resolve_run_settings(_sparse_si(basis="eb"), _DEFAULTS)
    assert resolved.basis == "eb"
    assert resolved.functional == "wigner"
    assert resolved.kpoint_spec == {"shift": [0.0, 0.0, 0.0]}
    assert resolved.scf_threshold == 1.0e-6


def test_unresolvable_setting_raises(tmp_path):
    # A solid that omits a setting with no [defaults] for it is a
    #   hard error (rule 2): the value is not resolvable.
    partial = dict(_DEFAULTS)
    del partial["scf_threshold"]
    with pytest.raises(ValueError, match="not resolvable"):
        _write_and_load(
            _manifest(_sparse_si(), defaults=partial), tmp_path)


def test_default_helpers_match_authoring_values():
    # The shared defaults the authoring tools emit (DESIGN 5.7).
    settings = default_run_settings()
    assert settings["basis"] == "fb"
    assert settings["functional"] == "wigner"
    assert settings["kpoint_integration"] == "gaussian"
    assert settings["kpoint_spec"] == {}
    assert settings["scf_threshold"] == 1.0e-6
    recipe = default_characterization()
    assert [fp.method for fp in recipe] == ["reduce", "bispectrum"]
    assert all(fp.preferred for fp in recipe)
