"""test_makeinput_build_api.py -- Unit tests for makeinput.py's
callable build API (C68(a); DESIGN 6.3; PSEUDOCODE 14).

C68(a) turns makeinput.py from an argv-and-cwd-bound script into
one that also exposes a callable ``build_run_dir(structure,
options, run_dir)``, so that ``imago.run_structure`` can drive it
in-process from a kaleidoscope worker (DESIGN 6.3.6).  A full
build runs the makeinput workflow end to end -- it reads a real
structure, the contract/basis databases, and invokes the
makeKPoints/contract executables -- so it has no pure-unit
harness.  These tests instead pin the *wrapper* logic the API
adds, with the heavy ``build_inputs`` workflow monkeypatched out:

* ``ScriptSettings.from_options`` / ``_args_from_options`` -- the
  options-mapping path that mirrors the C63 imago split: every
  argparse default is reproduced and the caller's options are
  overlaid, with an unknown key raising ``MakeinputError``
  (PSEUDOCODE 14.1).
* ``build_run_dir`` -- structure staging into ``run_dir`` and the
  cwd discipline (chdir in, restore on every exit including
  failure), plus the contract-fault raises (PSEUDOCODE 14.2 /
  DESIGN 6.3.4).
* ``record_clp`` -- the CLI-only ``command``-file append, now
  taking an explicit argv (DESIGN 6.3.5).
* ``MakeinputError`` -- raised, not ``sys.exit``, so a build
  fault propagates out of a worker instead of killing it
  (DESIGN 6.3.1).

conftest.py's ``SCRIPTS_DIR`` insertion lets us import
``makeinput`` directly without installing the package.
"""

import os
import types

import pytest

import makeinput
from makeinput import MakeinputError, ScriptSettings


# Pure / temp-file unit tests: the heavy makeinput workflow is
# monkeypatched out, so no Fortran binaries, no database reads,
# and no full build are needed.
pytestmark = pytest.mark.unit


def _bare_settings():
    """A ScriptSettings whose __init__ (and its rc load) is
    bypassed, for testing the methods that need an instance but
    not the resource-control defaults (``_args_from_options``,
    ``record_clp``)."""
    return ScriptSettings.__new__(ScriptSettings)


# ==============================================================
#  _args_from_options (PSEUDOCODE 14.1)
# ==============================================================

def test_args_from_options_overlays_known_dest():
    """A known argparse dest is overlaid onto the default
    namespace, and untouched dests keep their argparse defaults
    (so an option dict reproduces a bare invocation plus the
    overrides)."""
    settings = _bare_settings()
    args = settings._args_from_options({"basisdb": "/some/bdb"})
    assert args.basisdb == "/some/bdb"
    # potdb was not supplied, so it keeps its argparse default.
    assert args.potdb is None


def test_args_from_options_empty_is_all_defaults():
    """An empty options mapping yields exactly the argparse
    defaults -- the same namespace a bare ``makeinput`` sees."""
    settings = _bare_settings()
    args = settings._args_from_options({})
    assert args.basisdb is None
    assert args.potdb is None


def test_args_from_options_unknown_key_raises():
    """An option key that is not an argparse dest is a contract
    fault and raises MakeinputError, catching typos rather than
    silently ignoring them."""
    settings = _bare_settings()
    with pytest.raises(MakeinputError):
        settings._args_from_options({"not_a_real_option": 1})


def test_args_from_options_carries_loeninput_list():
    """The -loeninput dest (C57) is a registered argparse option, so
    the producer can drive the LOEN block in-process by passing a
    six-value list through the options mapping -- the same set that
    BispecMatcher.to_loen_input yields (DESIGN 5.10.5), in LOEN-block
    order."""
    settings = _bare_settings()
    args = settings._args_from_options(
        {"loeninput": ["1", "6", "4", "40", "7.5", "0.9"]})
    assert args.loeninput == ["1", "6", "4", "40", "7.5", "0.9"]
    # Omitting it leaves the dest at its argparse default of None, so a
    #   bare invocation falls back to the built-in contract values.
    assert _bare_settings()._args_from_options({}).loeninput is None


# ==============================================================
#  record_clp -- CLI-only command-file append (DESIGN 6.3.5)
# ==============================================================

def test_record_clp_writes_given_argv(tmp_path, monkeypatch):
    """record_clp appends the *passed* argv (not sys.argv) to the
    ``command`` file in the cwd, so the CLI controls what is
    recorded and the API path can decline to call it."""
    monkeypatch.chdir(tmp_path)
    _bare_settings().record_clp(["makeinput", "-scfkp", "4",
                                 "4", "4"])
    contents = (tmp_path / "command").read_text()
    assert "makeinput -scfkp 4 4 4" in contents


# ==============================================================
#  build_run_dir -- staging, cwd discipline, contract faults
# ==============================================================

def test_build_run_dir_missing_structure_raises(tmp_path):
    """A structure path that does not exist is a contract fault.
    An explicit (dummy) settings object is passed so the
    from_options default path is not exercised; the missing-file
    check happens before settings are otherwise used."""
    with pytest.raises(MakeinputError):
        makeinput.build_run_dir(
            str(tmp_path / "no_such.skl"), {},
            str(tmp_path / "run"), settings=object())


def test_build_run_dir_stages_skl_and_runs_in_run_dir(
        tmp_path, monkeypatch):
    """build_run_dir copies the structure into run_dir as
    imago.skl, runs build_inputs *with the cwd set to run_dir*,
    and restores the original cwd afterward."""
    source = tmp_path / "mystruct.skl"
    source.write_text("TITLE\nfake skeleton\n")
    run_dir = tmp_path / "run"
    seen = {}

    def fake_build_inputs(settings, sc):
        # Prove the workflow runs inside run_dir and that the
        #   skeleton was staged under the fixed name it reads.
        seen["cwd"] = os.getcwd()
        seen["staged"] = os.path.exists("imago.skl")

    monkeypatch.setattr(makeinput, "build_inputs",
                        fake_build_inputs)
    original_cwd = os.getcwd()
    returned = makeinput.build_run_dir(
        str(source), {}, str(run_dir), settings=object())

    assert returned == str(run_dir)
    assert seen["cwd"] == str(run_dir)        # ran in run_dir
    assert seen["staged"] is True             # skl was staged
    assert os.getcwd() == original_cwd        # cwd restored
    assert (run_dir / "imago.skl").read_text() == \
        "TITLE\nfake skeleton\n"


def test_build_run_dir_restores_cwd_on_failure(
        tmp_path, monkeypatch):
    """The cwd is restored even when the build raises -- the
    reentrancy guarantee that keeps one failed build from
    stranding a long-lived worker (DESIGN 6.3.4)."""
    source = tmp_path / "s.skl"
    source.write_text("x\n")
    run_dir = tmp_path / "run"

    def boom(settings, sc):
        raise MakeinputError("build blew up")

    monkeypatch.setattr(makeinput, "build_inputs", boom)
    original_cwd = os.getcwd()
    with pytest.raises(MakeinputError):
        makeinput.build_run_dir(str(source), {}, str(run_dir),
                                settings=object())
    assert os.getcwd() == original_cwd        # restored on fail


def test_build_run_dir_no_copy_when_source_is_staged(
        tmp_path, monkeypatch):
    """When the structure already IS run_dir/imago.skl, staging
    is a no-op (no self-copy) and the build still runs."""
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    staged = run_dir / "imago.skl"
    staged.write_text("already here\n")
    monkeypatch.setattr(makeinput, "build_inputs",
                        lambda settings, sc: None)
    makeinput.build_run_dir(str(staged), {}, str(run_dir),
                            settings=object())
    assert staged.read_text() == "already here\n"


# ==============================================================
#  from_options end-to-end (needs $IMAGO_RC for the rc defaults)
# ==============================================================

def test_from_options_builds_reconciled_settings(imago_data_dir):
    """from_options loads rc defaults and reconciles an options
    mapping into a usable settings object.  Skipped when the
    environment is not configured (the rc load needs $IMAGO_RC,
    asserted via the imago_data_dir fixture's session guard)."""
    if not os.getenv("IMAGO_RC"):
        pytest.skip("$IMAGO_RC not set")
    settings = ScriptSettings.from_options({})
    # The rc defaults were loaded by the constructor: a
    #   representative rc-backed attribute exists.
    assert hasattr(settings, "makekpoints_exec")


def test_loeninput_overrides_all_six_loen_settings(imago_data_dir):
    """-loeninput (C57) reconciles all six LOEN-block values, replacing
    the built-in bispectrum-contract defaults.  The first four parse as
    integers and the last two as reals, matching the LOEN-block order
    the engine reads (DESIGN 5.10.2 / 5.10.5)."""
    if not os.getenv("IMAGO_RC"):
        pytest.skip("$IMAGO_RC not set")
    # With no override the constructor leaves the contract defaults in
    #   place: code 1, the (4, 4) angular-momentum pair, 50 neighbors,
    #   a 9.0-Bohr cutoff, and an angleSqueeze of 0.85.
    plain = ScriptSettings.from_options({})
    assert (plain.loen_code, plain.loen_twoj1, plain.loen_twoj2,
            plain.loen_max_neigh, plain.loen_cutoff,
            plain.loen_angle_squeeze) == (1, 4, 4, 50, 9.0, 0.85)
    # An explicit -loeninput overlays every one of the six.
    overridden = ScriptSettings.from_options(
        {"loeninput": ["1", "6", "4", "40", "7.5", "0.9"]})
    assert overridden.loen_code == 1
    assert overridden.loen_twoj1 == 6
    assert overridden.loen_twoj2 == 4
    assert overridden.loen_max_neigh == 40
    assert overridden.loen_cutoff == 7.5
    assert overridden.loen_angle_squeeze == 0.9


# ==============================================================
#  _sort_atoms: the datSkl.map writer (DESIGN 5.2.1 / C87)
# ==============================================================

def test_sort_atoms_writes_datskl_map_with_identity_columns(
        tmp_path, monkeypatch):
    """``_sort_atoms`` writes ``inputs/datSkl.map`` with the five
    columns the C87 label-derivation reads: DAT#, SKELETON#, and
    each site's element symbol, species number, and type number.

    Two atoms (Si then O), already in element order, drive the
    writer directly with a synthetic settings/sc so the species and
    type columns are exercised without a full makeinput build."""

    settings = types.SimpleNamespace(
        num_atoms=2,
        atom_element_id=[0, 1, 2],          # 1-indexed
        atom_species_id=[0, 1, 1],
        atom_type_id=[[0, 1, 1]],           # [file_set][atom]
        atom_element_name=[None, "Si", "O"],
        xanes_atoms=[])
    sc = types.SimpleNamespace(
        fract_abc=[[None, 0, 0, 0], [None, 0.0, 0.0, 0.0],
                   [None, 0.5, 0.5, 0.5]],
        direct_abc=[[None, 0, 0, 0], [None, 0.0, 0.0, 0.0],
                    [None, 0.5, 0.5, 0.5]],
        direct_xyz=[[None, 0, 0, 0], [None, 0.0, 0.0, 0.0],
                    [None, 0.5, 0.5, 0.5]])

    monkeypatch.chdir(tmp_path)
    os.makedirs("inputs")
    makeinput._sort_atoms(settings, sc, file_set=0,
                          cumulative_num_types=None)

    rows = [line.split() for line in
            (tmp_path / "inputs" / "datSkl.map").read_text()
            .splitlines() if line.strip()]
    # Header carries the new column names; data rows carry element,
    #   species, and type alongside the dat<->skl numbering.
    assert rows[0] == ["DAT#", "SKELETON#", "ELEMENT",
                       "SPECIES", "TYPE"]
    assert rows[1] == ["1", "1", "si", "1", "1"]
    assert rows[2] == ["2", "2", "o", "1", "1"]
