"""test_imago_api.py -- Unit tests for imago.py's callable API
(C63; DESIGN 6.1; PSEUDOCODE 12).

C63 refactors imago.py from a command-line-only driver into a
callable Python API, with the CLI reduced to a thin wrapper.
The orchestration that actually stages inputs, locks the run
directory, runs the Fortran binary, and collects outputs needs
a real structure, the input databases, and the Imago binaries,
so it has no pure-unit harness.  These tests instead pin the
*pure* pieces of the API that carry the new logic and need
neither a binary nor $IMAGO_RC:

* ``project_home_outputs`` -- the single-source-of-truth table
  of project-home output filenames (PSEUDOCODE 12.2).
* ``_harvest_result`` -- builds the ImagoResult by reading the
  iteration file's last data row for the convergence verdict
  (column 4 vs the CONVERGENCE_TEST criterion), the total
  energy (column 5), and the per-run iteration count (column 1)
  (PSEUDOCODE 12.5).
* ``_read_convergence_threshold`` / ``_last_data_row`` -- the
  small parse helpers behind the harvest.
* ``run_structure`` -- confirmed to be an explicit
  not-yet-wired stub (C64 / C68).
* ``ImagoResult.success`` / ``_build_result`` -- the result
  object's derived state.

conftest.py's ``SCRIPTS_DIR`` insertion lets us import
``imago`` directly without installing the package.
"""

import os
from types import SimpleNamespace

import pytest

import imago
from imago import RunStatus, ImagoError, ImagoResult


# Pure-computation / temp-file unit tests -- no Fortran
# binaries, no full imago run, no $IMAGO_RC.
pytestmark = pytest.mark.unit


def _settings(edge="gs", job_name="scf", job_id=0,
              basis_scf="fb", basis_pscf="no"):
    """A stand-in for ScriptSettings carrying only the fields
    the pure helpers read.  Avoids needing $IMAGO_RC, which a
    real ScriptSettings construction requires."""
    return SimpleNamespace(
        edge=edge, job_name=job_name, job_id=job_id,
        basis_scf=basis_scf, basis_pscf=basis_pscf,
    )


# --------------------------------------------------------------
#  project_home_outputs (PSEUDOCODE 12.2)
# --------------------------------------------------------------

def test_outputs_scf_only():
    """A plain SCF emits the three SCF-always keys plus out,
    using the SCF basis tag, and no property key."""
    out = imago.project_home_outputs(_settings())
    assert out["scfV"] == "gs_scfV-fb.dat"
    assert out["energy"] == "gs_enrg-fb.dat"
    assert out["iteration"] == "gs_iter-fb.dat"
    assert out["out"] == "gs_scf-fb.out"
    assert "property" not in out
    assert "loen" not in out


def test_outputs_postscf_bond_uses_pscf_basis():
    """A 200-series job tags its files with the post-SCF basis
    and adds the property key for its property tag."""
    s = _settings(job_name="pscf+bond", job_id=202,
                  basis_scf="fb", basis_pscf="mb")
    out = imago.project_home_outputs(s)
    assert out["scfV"] == "gs_scfV-mb.dat"
    assert out["property"] == "gs_bond-mb.t"


def test_outputs_scf_no_omits_iteration():
    """A post-SCF property pass with no SCF (basis_scf == 'no')
    must not advertise the SCF-always keys, so the harvester
    will not read a stale iteration file."""
    s = _settings(job_name="pscf+dos", job_id=201,
                  basis_scf="no", basis_pscf="fb")
    out = imago.project_home_outputs(s)
    assert "iteration" not in out
    assert "scfV" not in out
    assert out["out"] == "gs_pscf+dos-fb.out"


def test_outputs_loen_key():
    """The local-environment job (311) adds a loen key and uses
    the fixed -fb basis tag of the 300-series."""
    s = _settings(job_name="loen", job_id=311,
                  basis_scf="no", basis_pscf="fb")
    out = imago.project_home_outputs(s)
    assert out["loen"] == "gs_loen-fb.plot"


# --------------------------------------------------------------
#  _read_convergence_threshold / _last_data_row (12.5)
# --------------------------------------------------------------

def test_read_convergence_threshold(tmp_path):
    """The criterion is the value on the line right after the
    CONVERGENCE_TEST label."""
    p = tmp_path / "imago.dat"
    p.write_text("OTHER\nCONVERGENCE_TEST\n  1.0e-4\nMORE\n")
    assert imago._read_convergence_threshold(str(p)) == 1.0e-4


def test_read_convergence_threshold_missing(tmp_path):
    """A missing label raises ImagoError naming the file."""
    p = tmp_path / "imago.dat"
    p.write_text("NO LABEL HERE\n")
    with pytest.raises(ImagoError):
        imago._read_convergence_threshold(str(p))


def test_last_data_row_skips_trailing_blanks(tmp_path):
    """The last *non-empty* line wins, split into fields."""
    p = tmp_path / "iter.dat"
    p.write_text("a b\n1 2 3\n\n  \n")
    assert imago._last_data_row(str(p)) == ["1", "2", "3"]


def test_last_data_row_empty(tmp_path):
    """An all-blank file yields None."""
    p = tmp_path / "iter.dat"
    p.write_text("\n  \n")
    assert imago._last_data_row(str(p)) is None


# --------------------------------------------------------------
#  _harvest_result (PSEUDOCODE 12.5)
# --------------------------------------------------------------

def _make_run_dir(tmp_path, last_delta, threshold=1.0e-4):
    """Build a fake converged-or-not run directory with an
    imago.dat carrying the criterion and a gs_iter-fb.dat whose
    final row has the given convergence delta in column 4."""
    (tmp_path / "imago.dat").write_text(
        f"CONVERGENCE_TEST\n  {threshold}\n")
    (tmp_path / "gs_iter-fb.dat").write_text(
        "# iter c2 c3 delta energy\n"
        "1 0.0 0.0 0.01 -10.0\n"
        f"2 0.0 0.0 {last_delta} -10.5\n")
    return tmp_path


def test_harvest_converged(tmp_path):
    """Final delta below the criterion => CONVERGED, with the
    iteration count (col 1) and total energy (col 5) read off
    the same last row."""
    d = _make_run_dir(tmp_path, last_delta="0.00005")
    r = imago._harvest_result(str(d), "/tmp/x", _settings(),
                              1.23, reused=False)
    assert r.status is RunStatus.CONVERGED
    assert r.success is True
    assert r.converged is True
    assert r.scf_iterations == 2
    assert abs(r.total_energy - (-10.5)) < 1e-12
    assert r.outputs["iteration"].endswith("gs_iter-fb.dat")
    assert r.runtime_seconds == 1.23


def test_harvest_not_converged(tmp_path):
    """Final delta above the criterion => NOT_CONVERGED, but
    the run still completed so the count/energy are read."""
    d = _make_run_dir(tmp_path, last_delta="0.005")
    r = imago._harvest_result(str(d), "/tmp/x", _settings(),
                              1.0, reused=False)
    assert r.status is RunStatus.NOT_CONVERGED
    assert r.success is False
    assert r.converged is False
    assert r.scf_iterations == 2


def test_harvest_no_scf_is_converged(tmp_path):
    """With no iteration output (a -scf no property pass) there
    is nothing to converge, so a clean run reports CONVERGED
    with no iteration count."""
    s = _settings(job_name="pscf+dos", job_id=201,
                  basis_scf="no", basis_pscf="fb")
    r = imago._harvest_result(str(tmp_path), "/tmp/x", s,
                              0.5, reused=False)
    assert r.status is RunStatus.CONVERGED
    assert r.scf_iterations is None
    assert r.total_energy is None


# --------------------------------------------------------------
#  Result object and entry-point contracts
# --------------------------------------------------------------

def test_build_result_failed():
    """_build_result yields a FAILED result that is not a
    success and carries the message."""
    r = imago._build_result(
        RunStatus.FAILED, "/run", "/tmp", _settings(), 2.0,
        message="boom")
    assert isinstance(r, ImagoResult)
    assert r.status is RunStatus.FAILED
    assert r.success is False
    assert r.message == "boom"
    assert r.scf_iterations is None


def test_success_only_for_converged():
    """success is True only for CONVERGED, not SKIPPED."""
    base = dict(run_dir="/r", temp_dir="/t",
                job=imago.JobIdentity("gs", "scf", "fb", "no"),
                runtime_seconds=0.0)
    assert ImagoResult(status=RunStatus.CONVERGED,
                       **base).success is True
    assert ImagoResult(status=RunStatus.SKIPPED,
                       **base).success is False
    assert ImagoResult(status=RunStatus.NOT_CONVERGED,
                       **base).success is False


def test_run_structure_is_stub():
    """structure-and-options mode is an explicit not-yet-wired
    stub (C64 / C68): it raises rather than silently failing."""
    with pytest.raises(ImagoError):
        imago.run_structure("x.skl", {}, "/tmp")


def test_run_prepared_missing_dir_raises():
    """A missing run directory is a contract fault.  An
    explicit settings object is passed so the check is
    exercised without needing $IMAGO_RC (the default
    settings=None path builds a real ScriptSettings)."""
    with pytest.raises(ImagoError):
        imago.run_prepared("/no/such/run/dir/at/all",
                           settings=_settings())
