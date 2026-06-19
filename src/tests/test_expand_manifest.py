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
    ignoring the prompt and default -- the test supplies exact text."""
    answer_iter = iter(answers)

    def ask(prompt, default):
        return next(answer_iter)

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
    # The stamped-but-entry-less manifest is structurally loadable.
    loaded = load_manifest_v2(str(out))
    solid = loaded.reference_solids[0]
    assert solid.basis == "fb"
    assert solid.functional == "wigner"
    assert solid.entries == []


def test_mechanical_output_carries_fill_in_template():
    text = em.ENTRY_TEMPLATE
    assert "Per-structure curation" in text
    assert "[[reference_solid.entry]]" in text
    # The canonical sub-specs are shown so the curator can copy them.
    assert "level = 2" in text and "twoj1 = 8" in text


# ==============================================================
#  Interactive mode (scripted ask)
# ==============================================================

def test_interactive_builds_complete_manifest(tmp_path):
    answers = [
        "fb", "wigner", "gaussian", "1e-6",   # shared defaults
        "crystalline",                        # system_type
        "y",                                  # add an entry?
        "Si", "1", "y", "Si diamond site 1", "",   # entry fields
        "y",                                  # attach fingerprints?
        "n",                                  # add another entry?
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
    methods = [fp.method for fp in entry.fingerprints]
    assert methods == ["reduce", "bispectrum"]
    assert all(fp.preferred for fp in entry.fingerprints)

    # And it round-trips through the strict loader.
    from curation_manifest import write_manifest
    out = tmp_path / "manifest.toml"
    write_manifest(manifest, str(out))
    assert load_manifest_v2(str(out)).reference_solids == \
        manifest.reference_solids


def test_interactive_marks_preferred_only_once_per_family(tmp_path):
    # Two Si entries, both fingerprinted: only the FIRST is preferred
    #   for each (element, method), so rule 10 holds.
    answers = [
        "fb", "wigner", "gaussian", "1e-6",
        "crystalline",
        "y", "Si", "1", "y", "site 1", "", "y",   # entry 1 (default)
        "y", "Si", "2", "n", "site 2", "", "y",   # entry 2
        "n",
    ]
    sources = load_structure_sources(_write_sketch(tmp_path))
    manifest = build_interactive(
        sources, _scripted(answers), **_SHARED)
    first, second = manifest.reference_solids[0].entries
    assert all(fp.preferred for fp in first.fingerprints)
    assert not any(fp.preferred for fp in second.fingerprints)
    # rule 10 (exactly one preferred per family) is satisfied:
    from curation_manifest import write_manifest
    out = tmp_path / "manifest.toml"
    write_manifest(manifest, str(out))
    load_manifest_v2(str(out))   # would raise on a rule-10 violation


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
