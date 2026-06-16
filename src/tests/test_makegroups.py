"""test_makegroups.py -- Unit and integration tests for makegroups.py,
the bispectrum grouping orchestrator (DESIGN 5.10; PSEUDOCODE 11.3.f;
TODO C58).

makegroups runs a sequence of real programs (makeinput.py, then
``imago.py -loen -scf no``) that need the Fortran engine, so the full
``group_by_bispectrum`` flow has no pure-unit harness.  These tests pin
the pieces that carry the grouping *logic* -- the P1 guard, the
per-element species bucketing, and the skeleton rewrite -- directly,
and exercise the orchestration with the two engine runs monkeypatched
out (a fake loen step drops a hand-written fort.21 in the work dir).

conftest.py's ``SCRIPTS_DIR`` insertion lets us import makegroups
directly without installing the package.  Tests that read a skeleton
through StructureControl need the space-group database and so take the
``imago_data_dir`` fixture, which skips them when $IMAGO_DATA is unset.
"""

import os
import types

import pytest

import makegroups
from makegroups import (
    _parse_atom_tag,
    _assign_species_from_rows,
    _rewrite_skeleton_species,
    _require_p1,
    _loen_input_values,
    group_by_bispectrum,
)
from matchers import MATCHERS, BispecMatcher, LoenSite


# A P1 skeleton with every atom its own species -- the case-1 input the
# grouping contract assumes (DESIGN 5.10.2): two Si and two O, distinctly
# tagged, in a 1x1x1 P1 cell.  Written into tmp_path by the helper below.
_P1_SKELETON = """\
title
Grouping test cell
end
cell
   5 5 5 90 90 90
space 1_a
supercell 1 1 1
full
frac 4
Si1  0.00000000  0.00000000  0.00000000
Si2  0.10000000  0.10000000  0.10000000
O1   0.20000000  0.20000000  0.20000000
O2   0.40000000  0.40000000  0.40000000
"""


def _write_skeleton(tmp_path, text=_P1_SKELETON, name="imago.skl"):
    """Write a skeleton file into tmp_path and return its path."""
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def _atom_lines(skeleton_path):
    """Return just the four atom-tag tokens of the test skeleton, in
    file order, for asserting on a rewrite."""
    lines = open(skeleton_path).read().splitlines()
    start = next(i for i, ln in enumerate(lines)
                 if ln.lower().startswith("frac")) + 1
    return [lines[i].split()[0] for i in range(start, start + 4)]


# ==============================================================
#  _parse_atom_tag
# ==============================================================

def test_parse_atom_tag_splits_element_and_species():
    """An element symbol plus an integer id splits into the two,
    preserving the element's letter case for a faithful rewrite."""
    assert _parse_atom_tag("Si12") == ("Si", 12)
    assert _parse_atom_tag("o2") == ("o", 2)


def test_parse_atom_tag_defaults_missing_number_to_one():
    """A tag with no trailing digit means species 1, matching the
    engine's own convention of appending a '1'."""
    assert _parse_atom_tag("C") == ("C", 1)


def test_parse_atom_tag_rejects_garbage():
    """A token that is not element-plus-digits is a fault, not a
    silent miss."""
    with pytest.raises(ValueError):
        _parse_atom_tag("12x!")


# ==============================================================
#  _assign_species_from_rows -- per-element bucketing
# ==============================================================

def test_assign_species_buckets_each_element_from_one():
    """Atoms of each element are clustered on their own and numbered
    from 1.  Two near-identical Si collapse to species 1; two O atoms
    far apart split into species 1 and 2 -- the O numbering restarts at
    1, independent of Si (DESIGN 5.10.4)."""
    matcher = BispecMatcher()
    rows = [
        LoenSite(1, "Si", 1, 1, 1, [0.0, 0.0]),
        LoenSite(2, "Si", 2, 1, 2, [0.05, 0.0]),
        LoenSite(3, "O", 1, 1, 3, [10.0, 0.0]),
        LoenSite(4, "O", 2, 1, 4, [99.0, 0.0]),
    ]
    assignment = _assign_species_from_rows(rows, matcher, 0.5)
    assert assignment == {
        ("si", 1): 1, ("si", 2): 1,
        ("o", 1): 1, ("o", 2): 2,
    }


def test_assign_species_never_merges_across_elements():
    """The bispectrum distance is element-blind, so even identical
    fingerprints on different elements must not share a species: each
    element gets its own species 1."""
    matcher = BispecMatcher()
    rows = [
        LoenSite(1, "Si", 1, 1, 1, [1.0, 2.0]),
        LoenSite(2, "O", 1, 1, 2, [1.0, 2.0]),
    ]
    assignment = _assign_species_from_rows(rows, matcher, 0.5)
    assert assignment == {("si", 1): 1, ("o", 1): 1}


# ==============================================================
#  _rewrite_skeleton_species -- the in-place tag rewrite
# ==============================================================

def test_rewrite_restarts_species_per_element(tmp_path):
    """Rewriting collapses Si2 -> Si1 and renumbers O2 -> O1 only when
    grouped together; here Si1,Si2 merge and O1,O2 stay distinct, so
    the tags become Si1, Si1, O1, O2 with each element starting at 1."""
    path = _write_skeleton(tmp_path)
    key = {("si", 1): 1, ("si", 2): 1, ("o", 1): 1, ("o", 2): 2}
    species_of = _rewrite_skeleton_species(path, key)
    assert species_of == [1, 1, 1, 2]
    assert _atom_lines(path) == ["Si1", "Si1", "O1", "O2"]


def test_rewrite_preserves_coordinates_and_other_lines(tmp_path):
    """Only the tag token changes: coordinates, the title, cell, space
    group, and supercell lines all survive verbatim."""
    path = _write_skeleton(tmp_path)
    key = {("si", 1): 1, ("si", 2): 1, ("o", 1): 1, ("o", 2): 1}
    _rewrite_skeleton_species(path, key)
    text = open(path).read()
    assert "0.10000000  0.10000000  0.10000000" in text
    assert "space 1_a" in text
    assert "supercell 1 1 1" in text
    assert "Grouping test cell" in text


def test_rewrite_missing_fingerprint_key_raises(tmp_path):
    """An atom whose (element, species) has no fingerprint row is a
    fault: the first ungrouped pass must give every atom its own
    species (DESIGN 5.10.2)."""
    path = _write_skeleton(tmp_path)
    # O2's key is absent from the map.
    key = {("si", 1): 1, ("si", 2): 1, ("o", 1): 1}
    with pytest.raises(ValueError):
        _rewrite_skeleton_species(path, key)


def test_rewrite_round_trip_through_structurecontrol(
        tmp_path, imago_data_dir):
    """The DESIGN 5.10.4 round-trip: group, write the skeleton, reread
    it through StructureControl, and recover the same per-element
    grouping.  After merging Si1,Si2 and keeping O1,O2 distinct, the
    reread tags are si1, si1, o1, o2 (atom_tag is 1-indexed, slot 0 is a
    placeholder)."""
    from structure_control import StructureControl
    path = _write_skeleton(tmp_path)
    key = {("si", 1): 1, ("si", 2): 1, ("o", 1): 1, ("o", 2): 2}
    _rewrite_skeleton_species(path, key)

    reread = StructureControl()
    reread.read_imago_skl(path, use_file_species=True)
    assert reread.atom_tag == [None, "si1", "si1", "o1", "o2"]


# ==============================================================
#  _require_p1 -- the crystal guard
# ==============================================================

def test_require_p1_accepts_a_p1_skeleton(sc_c2_molecule):
    """A genuine P1 cell (space group 1, 1x1x1 supercell) passes the
    guard without raising."""
    _require_p1(sc_c2_molecule, "c2.skl")


def test_require_p1_rejects_a_crystal(sc_si_diamond):
    """A symmetry-bearing skeleton is refused: grouping in P1 would drop
    the space group the k-point folding needs (DESIGN 5.10.1)."""
    with pytest.raises(ValueError, match="not P1"):
        _require_p1(sc_si_diamond, "si.skl")


def test_require_p1_rejects_a_nonunit_supercell():
    """Even at space group 1, a non-unit supercell is not the single
    cell grouping rewrites, so it is refused."""
    structure = types.SimpleNamespace(
        space_group_num=1, space_group="1_a",
        supercell=[None, 2, 1, 1])
    with pytest.raises(ValueError, match="not P1"):
        _require_p1(structure, "big.skl")


# ==============================================================
#  _loen_input_values -- the sub_spec -> CLI mapping
# ==============================================================

def test_loen_input_values_in_block_order():
    """The sub_spec maps through to_loen_input into the six -loeninput
    strings in LOEN-block order (DESIGN 5.10.5): method code, 2j1, 2j2,
    max_neigh, cutoff, angleSqueeze.  Unspecified parameters take the
    descriptor-contract defaults."""
    matcher = BispecMatcher()
    values = _loen_input_values(matcher, {"twoj1": 6, "twoj2": 4})
    assert values == ["1", "6", "4", "50", "9.0", "0.85"]


# ==============================================================
#  group_by_bispectrum -- the orchestration, engine runs faked
# ==============================================================

# A loen descriptor matching the four-atom P1 skeleton above, in the
# enriched self-describing format parse_loen_output reads (DESIGN
# 5.10.3): a header line, then one row per site in dat (element-sorted)
# order with five identity columns, five components (twoj2 + 1 = 5), and
# a trailing sum.  The two Si vectors are near-identical (they will
# merge); the two O vectors are far apart (they will split).
_FAKE_LOEN_DESCRIPTOR = """\
site# element species type_in_species type_flat  c0 c1 c2 c3 c4  sum
1 Si 1 1 1   0.00 0 0 0 0   0.0
2 Si 2 1 2   0.05 0 0 0 0   0.0
3 O 1 1 3   10.00 0 0 0 0   0.0
4 O 2 1 4   99.00 0 0 0 0   0.0
"""


def _fake_loen_writing_descriptor(work_dir):
    """Stand-in for _run_loen: drop the descriptor file the real engine
    leaves behind.  imago.py renames the loen fort.21 to the
    ``<edge>_loen<basis>.plot`` 1D-profile output, so the fake uses that
    naming (``gs_loen-fb.plot``) -- the exact name _find_loen_descriptor
    globs for, confirmed against a live run."""
    with open(os.path.join(work_dir, "gs_loen-fb.plot"), "w") as handle:
        handle.write(_FAKE_LOEN_DESCRIPTOR)


def test_group_by_bispectrum_round_trip(
        tmp_path, imago_data_dir, monkeypatch):
    """End to end with the two engine runs faked: the P1 guard passes,
    the faked loen drops a fort.21, the atoms bucket per element, and the
    skeleton is rewritten to Si1, Si1, O1, O2.  The original is preserved
    as <skeleton>.orig and the scratch dir is cleaned up."""
    path = _write_skeleton(tmp_path)
    original_text = open(path).read()

    monkeypatch.setattr(makegroups, "_run_makeinput",
                        lambda work_dir, loen_values: None)
    monkeypatch.setattr(makegroups, "_run_loen",
                        _fake_loen_writing_descriptor)

    species_of = group_by_bispectrum(
        path, {"twoj1": 4, "twoj2": 4}, similarity_floor=0.5)

    assert species_of == [1, 1, 1, 2]
    assert _atom_lines(path) == ["Si1", "Si1", "O1", "O2"]
    # The original skeleton is backed up untouched.
    assert open(path + ".orig").read() == original_text
    # The scratch directory is removed by default.
    assert not os.path.isdir(os.path.join(tmp_path, "makegroups_work"))


def test_group_by_bispectrum_refuses_crystal(
        tmp_path, imago_data_dir, monkeypatch):
    """The guard fires before any engine run: a crystalline skeleton is
    refused and no makeinput/loen subprocess is launched."""
    # Copy the diamond fixture (a real space group) into tmp_path.
    from conftest import STRUCTURES_DIR
    crystal = open(
        os.path.join(STRUCTURES_DIR, "si_diamond.skl")).read()
    path = _write_skeleton(tmp_path, text=crystal, name="si.skl")

    def _should_not_run(*args, **kwargs):
        raise AssertionError("engine run attempted on a crystal")

    monkeypatch.setattr(makegroups, "_run_makeinput", _should_not_run)
    monkeypatch.setattr(makegroups, "_run_loen", _should_not_run)

    with pytest.raises(ValueError, match="not P1"):
        group_by_bispectrum(path, {"twoj1": 4, "twoj2": 4})


# ==============================================================
#  CLI
# ==============================================================

def test_cli_defaults_match_descriptor_contract():
    """The CLI's bispectrum defaults are the descriptor-contract values,
    so a bare invocation groups with the same parameters the producer
    uses."""
    parser = makegroups._build_parser()
    args = parser.parse_args(["some.skl"])
    assert args.twoj1 == 4
    assert args.twoj2 == 4
    assert args.maxneigh == 50
    assert args.cutoff == 9.0
    assert args.anglesqueeze == 0.85
    assert args.floor is None


def test_cli_main_drives_group_by_bispectrum(monkeypatch, capsys):
    """main parses argv into a sub_spec and calls group_by_bispectrum,
    then reports the species count.  The grouping itself is stubbed."""
    captured = {}

    def _fake_group(skeleton_path, sub_spec, **kwargs):
        captured["skeleton"] = skeleton_path
        captured["sub_spec"] = sub_spec
        captured["kwargs"] = kwargs
        return [1, 1, 2]

    monkeypatch.setattr(makegroups, "group_by_bispectrum", _fake_group)
    makegroups.main(["my.skl", "-twoj1", "6", "-floor", "0.2"])

    assert captured["skeleton"] == "my.skl"
    assert captured["sub_spec"]["twoj1"] == 6
    assert captured["kwargs"]["similarity_floor"] == 0.2
    assert "2 species" in capsys.readouterr().out
