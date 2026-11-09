## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

"""test_makeinput_subfile_flavor.py -- Unit tests for carrying the
active build flavor into the generated submission file (DESIGN 14.6,
PSEUDOCODE 35.6).

A batch job begins by sourcing ``imagorc``, which points ``IMAGO_BIN``
at the production build.  A user working under a build flavor
(``imago_env mpi``, say) therefore needs the submission file to select
that flavor again, or the job runs a different build from the one in
their shell.  makeinput does this by writing ``imago_env <flavor> ||
exit 1`` directly after the ``source`` line.

Two things are tested, matching the two pieces of the pseudocode:

* ``active_flavor_name`` -- the pure function that reads the flavor
  out of ``IMAGO_DIR`` / ``IMAGO_BIN`` (production, a flavor, an
  unrecognized directory, and names that must never reach a shell
  script).
* ``_make_sub_file`` -- the writer, checked for where the line lands
  in each of the three job bodies and for the unchanged production
  output.

The check numbers in the test names are those of PSEUDOCODE 35.6.
Checks 11 and 12 need a scheduler and a built engine and are run by
hand, not here.  No Fortran binary and no database is touched.
"""

from types import SimpleNamespace

import pytest

import makeinput
from makeinput import active_flavor_name


pytestmark = pytest.mark.unit


# The install root every test pretends Imago lives under.  Nothing is
#   read from it; the function under test only compares path text.
INSTALL_ROOT = "/opt/imago"


def slurm_settings(ranks=1, xanes=0):
    """Return the handful of settings fields the writer reads.

    ``_make_sub_file`` consults only the queue type, the XANES switch
    and the scheduler resources, so a plain namespace stands in for
    the full settings object and no command line has to be parsed.
    """
    return SimpleNamespace(
        queue_type=3, xanes=xanes, ranks=ranks, nodes=1,
        cores_per_node=64, partition="some-partition", account="",
        time="01:00:00", memory="10G")


def write_sub_file(tmp_path, monkeypatch, imago_bin, file_set=0,
                   ranks=1, xanes=0):
    """Run the writer under a chosen IMAGO_BIN; return the file lines.

    ``imago_bin`` of None leaves the variable unset.  ``file_set`` is
    the writer's ``h`` argument: 0 is the base file set and a positive
    value is an XANES target.  The XANES control list is 1-indexed and
    holds one "n l" orbital entry, enough to produce one command.
    """
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("IMAGO_RC", f"{INSTALL_ROOT}/.imago")
    monkeypatch.setenv("IMAGO_DIR", INSTALL_ROOT)
    if imago_bin is None:
        monkeypatch.delenv("IMAGO_BIN", raising=False)
    else:
        monkeypatch.setenv("IMAGO_BIN", imago_bin)
    makeinput._make_sub_file(
        slurm_settings(ranks=ranks, xanes=xanes), file_set, "target",
        str(tmp_path), "deck", 1, [None, "1 0"])
    return (tmp_path / makeinput.SLURM_SUB).read_text().splitlines()


# ==============================================================
# active_flavor_name -- the resolver
# ==============================================================

@pytest.mark.parametrize("imago_dir, imago_bin", [
    (INSTALL_ROOT, f"{INSTALL_ROOT}/bin"),
    (INSTALL_ROOT + "/", f"{INSTALL_ROOT}//bin/"),
    (INSTALL_ROOT, ""),
    ("", f"{INSTALL_ROOT}/envs/mpi/bin"),
])
def test_check7_production_or_unset_names_no_flavor(imago_dir,
                                                    imago_bin):
    """Production, or nothing to compare, is the empty flavor."""
    assert active_flavor_name(imago_dir, imago_bin) == ("", False)


@pytest.mark.parametrize("imago_dir, imago_bin, expected_name", [
    (INSTALL_ROOT, f"{INSTALL_ROOT}/envs/mpi/bin", "mpi"),
    (INSTALL_ROOT + "/", f"{INSTALL_ROOT}//envs/mpi/bin/", "mpi"),
    (INSTALL_ROOT, f"{INSTALL_ROOT}/envs/asan-O0_v2/bin",
     "asan-O0_v2"),
])
def test_check8_flavor_directory_is_named(imago_dir, imago_bin,
                                          expected_name):
    """A flavor bin is recognized, slashes notwithstanding."""
    assert active_flavor_name(imago_dir, imago_bin) \
        == (expected_name, False)


@pytest.mark.parametrize("imago_bin", [
    "/some/overlay/bin",                     # outside the install
    f"{INSTALL_ROOT}/jobs/pa7/bin_mpi",      # ad-hoc overlay inside
    f"{INSTALL_ROOT}/envs/mpi",              # a flavor, but not its bin
    f"{INSTALL_ROOT}/envs/mpi/lib",          # wrong leaf directory
    f"{INSTALL_ROOT}/envs/a/b/bin",          # one level too deep
])
def test_check9_other_directories_are_unrecognized(imago_bin):
    """Anything that is not production or a flavor bin is flagged."""
    assert active_flavor_name(INSTALL_ROOT, imago_bin) == ("", True)


@pytest.mark.parametrize("hostile_name", [
    "a;rm -rf x",
    "two words",
    "$HOME",
    "$(touch pwned)",
    "`id`",
    "mpi\n",
    "mpi\nrm -rf x",
    "mpi||true",
])
def test_check10_names_outside_the_alphabet_are_refused(hostile_name):
    """A name that could act as shell text never becomes a flavor."""
    imago_bin = f"{INSTALL_ROOT}/envs/{hostile_name}/bin"
    assert active_flavor_name(INSTALL_ROOT, imago_bin) == ("", True)


# ==============================================================
# _make_sub_file -- the writer
# ==============================================================

def test_check7_production_file_has_no_flavor_line(
        tmp_path, monkeypatch, capsys):
    """Under production the file is silent about flavors."""
    lines = write_sub_file(tmp_path, monkeypatch,
                           f"{INSTALL_ROOT}/bin")
    assert not any("imago_env" in line for line in lines)
    source_at = lines.index(f"source {INSTALL_ROOT}/.imago/imagorc")
    assert lines[source_at + 1] == "export OMP_NUM_THREADS=1"
    assert capsys.readouterr().err == ""


def test_check7_unset_imago_bin_has_no_flavor_line(
        tmp_path, monkeypatch, capsys):
    """With IMAGO_BIN unset there is nothing to carry and no noise."""
    lines = write_sub_file(tmp_path, monkeypatch, None)
    assert not any("imago_env" in line for line in lines)
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize("file_set, xanes", [
    (0, 0),      # the plain body
    (0, 1),      # the XANES base body
    (1, 1),      # an XANES target body
])
def test_check8_flavor_line_follows_source_in_every_body(
        tmp_path, monkeypatch, file_set, xanes):
    """Each body selects the flavor directly after sourcing imagorc."""
    lines = write_sub_file(
        tmp_path, monkeypatch, f"{INSTALL_ROOT}/envs/mpi/bin",
        file_set=file_set, ranks=4, xanes=xanes)
    source_at = lines.index(f"source {INSTALL_ROOT}/.imago/imagorc")
    assert lines[source_at + 1] == "imago_env mpi || exit 1"
    assert sum("imago_env" in line for line in lines) == 1
    # The launcher export must still come after the source line,
    #   because imagorc resets IMAGO_EXE to empty.
    launcher_at = lines.index(
        'export IMAGO_EXE="mpirun -np $SLURM_NTASKS"')
    assert launcher_at > source_at


def test_check9_overlay_directory_warns_and_writes_no_line(
        tmp_path, monkeypatch, capsys):
    """An overlay bin cannot be carried over; say so, do not abort."""
    lines = write_sub_file(tmp_path, monkeypatch, "/some/overlay/bin")
    assert not any("imago_env" in line for line in lines)
    warning = capsys.readouterr().err
    assert "WARNING" in warning
    assert "/some/overlay/bin" in warning


def test_check10_hostile_name_never_reaches_the_file(
        tmp_path, monkeypatch, capsys):
    """Shell text in the flavor position is not written anywhere."""
    lines = write_sub_file(
        tmp_path, monkeypatch,
        f"{INSTALL_ROOT}/envs/a;touch pwned/bin")
    assert not any("imago_env" in line for line in lines)
    assert not any("pwned" in line for line in lines)
    assert "WARNING" in capsys.readouterr().err
