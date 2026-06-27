"""Tests for expand_manifest, the sketch-to-manifest authoring tool.

The pure builders (:func:`stamp_shared_defaults`, :func:`build_
mechanical`, :func:`build_interactive`) are exercised directly; the
interactive flow is driven by an injected scripted ``ask`` rather than
real prompts, so it is fully deterministic.  The central guarantee is
that what either mode produces reads back through
``load_manifest_v2``.
"""

import pytest

import expand_manifest as em
from expand_manifest import (
    stamp_shared_defaults,
    build_mechanical,
    build_interactive,
)
from curation_manifest import load_structure_sources, load_manifest_v2

pytestmark = pytest.mark.unit


# A sketch as the curator would assemble it from cod_fish stubs: only
#   identity and structure source per solid, no run or harvest fields.
_SKETCH = """\
schema_version = 2

[[reference_solid]]
reference_id = "si_diamond"
system_type = "crystalline"
cod_id = 2104737
cod_revision = "201401"
"""


def _write_sketch(tmp_path, text=_SKETCH):
    path = tmp_path / "sketch.toml"
    path.write_text(text)
    return str(path)


def _scripted(answers):
    """An ``ask`` callable that returns successive scripted answers,
    ignoring the prompt.  An answer of ``None`` means "press Enter",
    so the call returns the ``default`` -- matching the real prompt's
    accept-the-default behaviour, which lets a test exercise the
    hint-supplied element and description defaults."""
    answer_iter = iter(answers)

    def ask(prompt, default):
        answer = next(answer_iter)
        return default if answer is None else answer

    return ask


_SHARED = dict(
    basis="fb", functional="wigner",
    kpoint_integration="gaussian", scf_threshold=1.0e-6,
    system_type_default="crystalline")


# ==============================================================
#  stamp_shared_defaults
# ==============================================================

def test_stamp_fills_run_settings_and_keeps_source(tmp_path):
    source = load_structure_sources(_write_sketch(tmp_path))[0]
    solid = stamp_shared_defaults(source, **_SHARED)
    assert solid.reference_id == "si_diamond"
    assert solid.cod_id == 2104737
    assert solid.cod_revision == "201401"
    assert solid.basis == "fb"
    assert solid.functional == "wigner"
    assert solid.kpoint_integration == "gaussian"
    assert solid.scf_threshold == pytest.approx(1.0e-6)
    assert solid.system_type == "crystalline"   # from the sketch
    # k-density is left to predict-then-verify; entries unfilled.
    assert solid.kpoint_spec == {}
    assert solid.entries == []


def test_stamp_uses_system_type_default_when_sketch_blank(tmp_path):
    # A sketch with no system_type falls back to the supplied default.
    sketch = ('schema_version = 2\n\n[[reference_solid]]\n'
              'reference_id = "x"\nstructure_path = "x.skl"\n')
    source = load_structure_sources(_write_sketch(tmp_path, sketch))[0]
    solid = stamp_shared_defaults(
        source, **{**_SHARED, "system_type_default": "amorphous"})
    assert solid.system_type == "amorphous"


# ==============================================================
#  Mechanical mode
# ==============================================================

def test_mechanical_stamps_defaults_and_loads(tmp_path):
    sources = load_structure_sources(_write_sketch(tmp_path))
    manifest = build_mechanical(sources, **_SHARED)
    out = tmp_path / "manifest.toml"
    from curation_manifest import format_manifest
    out.write_text(format_manifest(manifest) + "\n" + em.ENTRY_TEMPLATE)
    # The stamped manifest carries the required [characterization]
    #   recipe, so it loads with no entries (the harvest fills them).
    loaded = load_manifest_v2(str(out))
    solid = loaded.reference_solids[0]
    assert solid.basis == "fb"
    assert solid.functional == "wigner"
    assert solid.entries == []
    methods = {fp.method for fp in loaded.characterization}
    assert methods == {"reduce", "bispectrum"}
    assert all(fp.preferred for fp in loaded.characterization)


def test_mechanical_output_carries_fill_in_template():
    text = em.ENTRY_TEMPLATE
    assert "Per-structure customizations" in text
    assert "[[reference_solid.entry]]" in text
    # The template shows a RARE per-entry override (non-preferred);
    #   the preferred recipe lives in the [characterization] block.
    assert "preferred = true" not in text
    assert "twoj1 = 6, twoj2 = 4" in text


# ==============================================================
#  Interactive mode (scripted ask)
# ==============================================================

def test_interactive_builds_complete_manifest(tmp_path):
    answers = [
        "fb", "wigner", "gaussian", "1e-6",   # shared defaults
        "crystalline",                        # system_type
        "y",                                  # add a customization?
        "Si", "1", "y", "Si diamond site 1", "",   # entry fields
        "n",                                  # add another?
    ]
    sources = load_structure_sources(_write_sketch(tmp_path))
    manifest = build_interactive(
        sources, _scripted(answers), **_SHARED)

    solid = manifest.reference_solids[0]
    assert solid.basis == "fb"
    assert solid.scf_threshold == pytest.approx(1.0e-6)
    entry = solid.entries[0]
    assert entry.element == "Si"
    assert entry.atom_site == 1
    assert entry.default is True
    assert entry.label is None            # blank reply -> derived
    # No per-entry fingerprints are authored: the preferred recipe
    #   is the database-wide [characterization] block.
    assert entry.fingerprints == []
    methods = {fp.method for fp in manifest.characterization}
    assert methods == {"reduce", "bispectrum"}
    assert all(fp.preferred for fp in manifest.characterization)

    # And it round-trips through the strict loader.
    from curation_manifest import write_manifest
    out = tmp_path / "manifest.toml"
    write_manifest(manifest, str(out))
    reloaded = load_manifest_v2(str(out))
    assert reloaded.reference_solids == manifest.reference_solids
    assert reloaded.characterization == manifest.characterization


def test_interactive_sets_database_wide_recipe_once(tmp_path):
    # Two Si customizations: neither carries a per-entry fingerprint,
    #   and the manifest carries the one database-wide preferred
    #   recipe, so the strict loader (which enforces rules 2 and 10)
    #   accepts it.
    answers = [
        "fb", "wigner", "gaussian", "1e-6",
        "crystalline",
        "y", "Si", "1", "y", "site 1", "",   # entry 1 (default)
        "y", "Si", "2", "n", "site 2", "",   # entry 2
        "n",
    ]
    sources = load_structure_sources(_write_sketch(tmp_path))
    manifest = build_interactive(
        sources, _scripted(answers), **_SHARED)
    first, second = manifest.reference_solids[0].entries
    assert first.fingerprints == []
    assert second.fingerprints == []
    assert len(manifest.characterization) == 2
    # rules 2 and 10 hold, so the strict loader accepts the file.
    from curation_manifest import write_manifest
    out = tmp_path / "manifest.toml"
    write_manifest(manifest, str(out))
    load_manifest_v2(str(out))   # would raise on a rule violation


# ==============================================================
#  Discovery hints (cod_fish-written) auto-fill element + desc
# ==============================================================

_HINTED_SKETCH = """\
schema_version = 2

[[reference_solid]]
reference_id = "si_imma_74_1993"
system_type = "crystalline"
cod_id = 9011656
cod_revision = "291877"
elements = ["Si"]
source_description = "Silicon, I m m a (74), 1993"
"""


def test_source_description_persists_through_expand(tmp_path):
    # The CIF-derived source_description is a persisted reference_solid
    #   field (DESIGN 5.7): stamp_shared_defaults carries it from the
    #   sketch onto the finished solid, and it round-trips through the
    #   strict loader.  (`elements`, by contrast, is a transient hint
    #   the finished manifest omits.)
    sources = load_structure_sources(
        _write_sketch(tmp_path, _HINTED_SKETCH))
    manifest = build_mechanical(sources, **_SHARED)
    solid = manifest.reference_solids[0]
    assert solid.source_description == "Silicon, I m m a (74), 1993"

    from curation_manifest import write_manifest
    out = tmp_path / "manifest.toml"
    write_manifest(manifest, str(out))
    reloaded = load_manifest_v2(str(out))
    assert reloaded.reference_solids[0].source_description == \
        "Silicon, I m m a (74), 1993"


def test_load_sketch_hints_reads_elements_and_description(tmp_path):
    path = _write_sketch(tmp_path, _HINTED_SKETCH)
    hints = em.load_sketch_hints(path)
    assert hints["si_imma_74_1993"]["elements"] == ["Si"]
    assert hints["si_imma_74_1993"]["description"] == \
        "Silicon, I m m a (74), 1993"


def test_load_sketch_hints_empty_without_fields(tmp_path):
    # A hand-written sketch lacking the hint fields yields empties,
    #   so the interactive flow just prompts with blank defaults.
    path = _write_sketch(tmp_path)            # the plain _SKETCH
    hints = em.load_sketch_hints(path)
    assert hints["si_diamond"] == {"elements": [], "description": ""}


def test_interactive_accepts_hinted_element_and_description(tmp_path):
    # None means "press Enter": element and description are left to
    #   their hint-supplied defaults.
    answers = [
        None, None, None, None,   # shared defaults
        None,                     # system_type
        "y",                      # add a customization?
        None, None, "y", None, None,   # element, atom_site, default,
                                       #   description, label (all Enter)
        "n",                      # add another?
    ]
    sources = load_structure_sources(_write_sketch(tmp_path))
    hints = {"si_diamond": {
        "elements": ["Si"],
        "description": "Silicon, F d -3 m (227), 2010"}}
    manifest = build_interactive(
        sources, _scripted(answers), hints=hints, **_SHARED)
    entry = manifest.reference_solids[0].entries[0]
    assert entry.element == "Si"
    assert entry.description == "Silicon, F d -3 m (227), 2010"


# ==============================================================
#  CLI
# ==============================================================

def test_main_mechanical_writes_loadable_manifest(tmp_path):
    sketch = _write_sketch(tmp_path)
    out = tmp_path / "out.toml"
    rc = em.main([sketch, "-o", str(out)])
    assert rc == 0
    loaded = load_manifest_v2(str(out))
    assert loaded.reference_solids[0].functional == "wigner"


def test_main_interactive_requires_output(tmp_path):
    sketch = _write_sketch(tmp_path)
    with pytest.raises(SystemExit):
        em.main([sketch, "-i"])
