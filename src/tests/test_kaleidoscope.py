"""test_kaleidoscope.py -- Unit tests for the kaleidoscope
flight dispatcher (C68; DESIGN 6.2; PSEUDOCODE 13).

kaleidoscope drives a *set* of Imago calculations: it dispatches
the per-structure work, tracks each one's outcome, caches
completed runs so a flight resumes by re-running, and surfaces
a report the client harvests.  Per VISION Principle 9 it is
domain-agnostic -- it never interprets what a run computed.  That
is exactly what lets these tests exercise the whole machinery
*without* an Imago binary or $IMAGO_RC: a fake wingbeat stands in
for the real one and reports whatever generic outcome a test
wants.  The pieces pinned here, helpers-first then the driver:

* ``validate_flight`` -- the slug rule, the derived ``<calc>``
  tag when one id hosts several units, and the run-directory
  collision guard (PSEUDOCODE 13.3).
* ``unit_run_dir`` -- the ``<root>/wingbeats/<id>[/<calc>]`` layout.
* the cache hit-test -- verbatim scalar compare plus key-file
  byte-comparison, and the done-status precondition
  (PSEUDOCODE 13.4).
* ``dispatch`` under *both* executors (LocalExecutor and a
  Parsl ThreadPoolExecutor config) -- dispatch, the status
  lifecycle, complete-and-report (one failure never aborts the
  batch), resume-skips-done, and the ``on_outcome`` hook
  (PSEUDOCODE 13.5).
* ``ImagoWingbeat`` -- the ImagoResult -> WingbeatOutcome mapping and
  the persisted ``result.toml`` handoff, with imago's run entry
  points monkeypatched so no binary is needed (PSEUDOCODE 13.2).
* the report views ``by_status`` / ``with_detail`` / ``failures``
  the client selects on (PSEUDOCODE 13.6).

conftest.py's ``SCRIPTS_DIR`` insertion lets us import
``kaleidoscope`` (and ``imago``) directly without installing the
package.
"""

import os
import tomllib

import pytest

from kaleidoscope import (
    Flight, CalcUnit, KeyFields, KeyFile, WingbeatOutcome,
    FlightReport, ReportEntry, KaleidoscopeError, SweepRecord,
    dispatch, register_wingbeat, resolve_wingbeat,
    validate_flight, unit_run_dir,
    is_cache_hit, cache_key_matches, write_cache_key,
    read_status, write_status, ImagoWingbeat,
)
from kaleidoscope.wingbeats import Wingbeat
from kaleidoscope.workspace import derive_calc_tag, serialize_flight


# Pure-computation / temp-file tests -- they read no fixture
# files beyond what tmp_path holds, need no Imago binary, and
# never touch $IMAGO_RC, so the whole module is a unit suite.
pytestmark = pytest.mark.unit


# ==============================================================
#  Fake wingbeats (stand-ins for the real ImagoWingbeat)
# ==============================================================

class CountingRunner(Wingbeat):
    """A wingbeat that always completes and counts how many times
    it actually executed.  Used to prove the cache skips a unit
    on resume: a cache hit must NOT call the wingbeat again, so the
    count stays put across a second ``dispatch``."""

    def __init__(self):
        self.calls = 0

    def run(self, unit, wingbeat_dir):
        self.calls += 1
        return WingbeatOutcome(ok=True, detail="converged",
                          runtime_seconds=0.1, message="")


class ModalRunner(Wingbeat):
    """A wingbeat whose behavior a unit selects through its
    options, so one flight can mix outcomes.  ``fake_mode`` is
    one of: ``ok`` (completes, detail "converged"); ``not_ok``
    (completes but WingbeatOutcome.ok is False -> status "failed");
    ``raise`` (raises on the worker, exercising the dispatcher's
    per-future capture -> status "failed")."""

    def run(self, unit, wingbeat_dir):
        mode = unit.options.get("fake_mode", "ok")
        if mode == "raise":
            raise RuntimeError("kaboom from the fake wingbeat")
        if mode == "not_ok":
            return WingbeatOutcome(ok=False, detail="not_converged",
                              runtime_seconds=0.0, message="hit ceiling")
        return WingbeatOutcome(ok=True, detail="converged",
                          runtime_seconds=0.2, message="")


# ==============================================================
#  Executor selection -- run each driver test on both paths
# ==============================================================

def _parsl_thread_config():
    """A minimal Parsl ``Config`` backed by a thread pool -- the
    laptop deployment.  Threads run in-process, so the wingbeat
    instance registered in this process is the very one a worker
    resolves (no pickling), and no cluster is needed."""
    from parsl.config import Config
    from parsl.executors.threads import ThreadPoolExecutor
    return Config(executors=[ThreadPoolExecutor(
        max_threads=2, label="kaleidoscope_test_threads")])


@pytest.fixture(params=["local", "parsl"])
def parsl_config(request):
    """Parametrize driver tests over both executors.  ``local``
    yields None (LocalExecutor, synchronous in-process); ``parsl``
    yields a thread-pool Config (ParslExecutor), skipped when
    Parsl is not installed.  dispatch chooses the executor
    from the presence of this config."""
    if request.param == "parsl":
        pytest.importorskip("parsl")
        return _parsl_thread_config()
    return None


# ==============================================================
#  13.3 -- workspace: slugs, run dirs, validate_flight
# ==============================================================

def test_unit_run_dir_with_and_without_calc(tmp_path):
    """The run directory is ``<root>/wingbeats/<id>``; the ``<calc>``
    level appears only when the unit carries a calc tag."""
    flight = Flight(root=str(tmp_path), units=[])
    plain = CalcUnit(id="s1", structure="a.skl")
    assert unit_run_dir(flight, plain) == os.path.join(
        str(tmp_path), "wingbeats", "s1")
    tagged = CalcUnit(id="s1", structure="a.skl", calc=("v2",))
    assert unit_run_dir(flight, tagged) == os.path.join(
        str(tmp_path), "wingbeats", "s1", "v2")
    # A multi-axis calc tuple nests one directory level per
    #   component, in tuple order (DESIGN 6.2.1/6.2.4).
    multi = CalcUnit(id="s1", structure="a.skl",
                     calc=("kpt-density-50", "smear-gauss"))
    assert unit_run_dir(flight, multi) == os.path.join(
        str(tmp_path), "wingbeats", "s1",
        "kpt-density-50", "smear-gauss")


def test_validate_rejects_non_slug_id():
    """An id that is not a filesystem-safe slug aborts the
    flight rather than being silently rewritten -- a rewrite
    would break the cache hit-test (different directory)."""
    flight = Flight(root="/tmp",
                        units=[CalcUnit(id="Bad ID!", structure="a")])
    with pytest.raises(KaleidoscopeError):
        validate_flight(flight)


def test_validate_derives_calc_for_shared_id():
    """When two units share one id but neither names a calc, a
    distinguishing ``<job>-<scf_basis>`` tag is derived in place
    so the two no longer collide on one run directory."""
    first = CalcUnit(id="s1", structure="a",
                     options={"job": "scf", "scf_basis": "fb"})
    second = CalcUnit(id="s1", structure="a",
                      options={"job": "pscf", "scf_basis": "mb"})
    validate_flight(Flight(root="/tmp",
                               units=[first, second]))
    assert first.calc == ("scf-fb",)
    assert second.calc == ("pscf-mb",)


def test_derive_calc_tag_defaults():
    """With no makeinput options the derived tag falls back to
    the documented ``scf``/``fb`` defaults, as a one-element
    tuple matching the CalcUnit.calc shape."""
    assert derive_calc_tag(CalcUnit(id="s1", structure="a")) == \
        ("scf-fb",)


def test_validate_duplicate_run_dir_raises():
    """Two units that resolve to the same id+calc would clobber
    one run directory, so validation aborts and names them."""
    first = CalcUnit(id="s1", structure="a", calc=("v1",))
    second = CalcUnit(id="s1", structure="a", calc=("v1",))
    with pytest.raises(KaleidoscopeError):
        validate_flight(Flight(root="/tmp",
                                   units=[first, second]))


def test_serialize_flight_calc_array_and_no_calc(tmp_path):
    """Each unit's calc tuple serializes as a TOML array that
    parses back to a list; a unit with no calc emits ``calc =
    []`` (DESIGN 6.2.1)."""
    swept = CalcUnit(id="s1", structure="a.skl",
                     calc=("kpt-density-50",))
    plain = CalcUnit(id="s2", structure="b.skl")
    serialize_flight(Flight(root=str(tmp_path),
                            units=[swept, plain]))
    with open(os.path.join(str(tmp_path), "flight.toml"),
              "rb") as flight_file:
        data = tomllib.load(flight_file)
    assert data["unit"][0]["calc"] == ["kpt-density-50"]
    assert data["unit"][1]["calc"] == []


def test_serialize_flight_sweep_and_metadata_round_trip(tmp_path):
    """A predict-then-verify flight emits a [flight.sweep] block
    (with its nested fixed_axes sub-table) and each metadata key
    as a verbatim [flight.<key>] table that round-trips through
    tomllib unchanged (DESIGN 6.2.8; the harvest recovers these
    without parsing run-dir paths)."""
    sweep = SweepRecord(
        varied_axes=("kpt-density",),
        fixed_axes={"basis": "fb", "functional": "ldau"})
    # A flat prediction table standing in for the 6.2.8 helper's
    #   PredictionRecord: scalars plus a string array.
    prediction = {
        "policy": "verify_around_prediction",
        "predicted_value": 50.0,
        "confidence": 0.83,
        "is_under_trained": False,
        "neighbor_entry_ids": ["mp-1", "mp-2"],
    }
    unit = CalcUnit(id="s1", structure="a.skl",
                    calc=("kpt-density-50",))
    flight = Flight(root=str(tmp_path), units=[unit], sweep=sweep,
                    metadata={"prediction": prediction})
    serialize_flight(flight)
    with open(os.path.join(str(tmp_path), "flight.toml"),
              "rb") as flight_file:
        data = tomllib.load(flight_file)

    assert data["flight"]["sweep"]["varied_axes"] == ["kpt-density"]
    assert data["flight"]["sweep"]["fixed_axes"] == {
        "basis": "fb", "functional": "ldau"}
    assert data["flight"]["prediction"] == prediction


def test_resolve_unknown_runner_raises():
    """Asking for a wingbeat that was never registered is a
    flight-construction fault, not a silent default."""
    with pytest.raises(KaleidoscopeError):
        resolve_wingbeat("no-such-wingbeat-name")


# ==============================================================
#  13.3 -- status.toml merge/read lifecycle
# ==============================================================

def test_write_status_merges_across_lifecycle(tmp_path):
    """A later terminal write must preserve fields accumulated
    earlier (``submitted_at`` from queue time survives the
    ``done`` write).  None-valued fields are skipped."""
    write_status(str(tmp_path), status="queued",
                 submitted_at="t0", calc=None)
    write_status(str(tmp_path), status="done",
                 detail="converged", finished_at="t1")
    status = read_status(str(tmp_path))
    assert status["submitted_at"] == "t0"
    assert status["status"] == "done"
    assert status["detail"] == "converged"
    assert status["finished_at"] == "t1"
    assert "calc" not in status        # None was skipped


def test_read_status_absent_is_none(tmp_path):
    """A run directory with no status.toml reads back as None."""
    assert read_status(str(tmp_path)) is None


# ==============================================================
#  13.4 -- cache hit-test (scalars verbatim, files byte-compare)
# ==============================================================

def _staged_unit(tmp_path, source_text, staged_text,
                 scalars):
    """Build a (wingbeat_dir, unit) pair for a cache test: write the
    current source file, stage a copy under the run directory
    (as a prior run would have left it), snapshot the key, and
    return both.  ``staged_text`` may differ from ``source_text``
    to exercise a byte mismatch."""
    wingbeat_dir = tmp_path / "run"
    wingbeat_dir.mkdir()
    source = tmp_path / "structure.skl"
    source.write_text(source_text)
    (wingbeat_dir / "structure.skl").write_text(staged_text)
    unit = CalcUnit(
        id="s1", structure=str(source),
        key_fields=KeyFields(
            scalars=scalars,
            files=[KeyFile(name="structure.skl",
                           source=str(source))]))
    write_cache_key(str(wingbeat_dir), unit)
    return str(wingbeat_dir), unit


def test_cache_matches_when_scalars_and_files_agree(tmp_path):
    """A unit whose scalars equal the snapshot and whose key file
    byte-equals its staged copy is a key match."""
    wingbeat_dir, unit = _staged_unit(
        tmp_path, "LATTICE 1 2 3\n", "LATTICE 1 2 3\n",
        {"kpoints": "4x4x4", "threshold": 0.0001})
    assert cache_key_matches(unit, wingbeat_dir) is True


def test_cache_misses_on_changed_scalar(tmp_path):
    """A single differing scalar field is a miss -- the key is
    compared verbatim, field by field."""
    wingbeat_dir, _ = _staged_unit(
        tmp_path, "LATTICE 1 2 3\n", "LATTICE 1 2 3\n",
        {"kpoints": "4x4x4"})
    changed = CalcUnit(
        id="s1", structure="x",
        key_fields=KeyFields(scalars={"kpoints": "6x6x6"}))
    assert cache_key_matches(changed, wingbeat_dir) is False


def test_cache_misses_on_byte_differing_key_file(tmp_path):
    """When the current source no longer byte-equals the staged
    copy the cache misses, even though the names match."""
    wingbeat_dir, unit = _staged_unit(
        tmp_path, "LATTICE 9 9 9\n", "LATTICE 1 2 3\n",
        {"kpoints": "4x4x4"})
    assert cache_key_matches(unit, wingbeat_dir) is False


def test_is_cache_hit_requires_done_status(tmp_path):
    """A matching key is necessary but not sufficient: the run
    must also have reached the ``done`` status.  A still-running
    directory is a miss, so the unit is relaunched."""
    wingbeat_dir, unit = _staged_unit(
        tmp_path, "X\n", "X\n", {"v": 1})
    # No status.toml yet -> miss.
    assert is_cache_hit(unit, wingbeat_dir) is False
    write_status(wingbeat_dir, status="running")
    assert is_cache_hit(unit, wingbeat_dir) is False
    write_status(wingbeat_dir, status="done")
    assert is_cache_hit(unit, wingbeat_dir) is True


# ==============================================================
#  13.5 -- the dispatch driver, on BOTH executors
# ==============================================================

def test_flight_runs_and_reports_done(tmp_path, parsl_config):
    """A clean unit runs to ``done`` with the wingbeat's ``detail``
    recorded, flight.toml is written, and the report carries
    the entry in unit order."""
    register_wingbeat("fake_ok", CountingRunner())
    unit = CalcUnit(id="u1", structure="s.skl", wingbeat="fake_ok",
                    key_fields=KeyFields(scalars={"v": 1}))
    flight = Flight(root=str(tmp_path), units=[unit],
                        parsl_config=parsl_config)
    report = dispatch(flight)

    assert len(report.entries) == 1
    entry = report.entries[0]
    assert entry.id == "u1"
    assert entry.status == "done"
    assert entry.detail == "converged"
    assert os.path.exists(os.path.join(str(tmp_path),
                                       "flight.toml"))


def test_status_lifecycle_fields_present(tmp_path, parsl_config):
    """After a successful run the run directory's status.toml
    carries the full lifecycle: a queued-time ``submitted_at``,
    a worker-time ``started_at``, and the terminal ``done`` plus
    ``detail`` / ``finished_at`` / ``runtime_seconds``."""
    register_wingbeat("fake_ok", CountingRunner())
    unit = CalcUnit(id="u1", structure="s.skl", wingbeat="fake_ok",
                    key_fields=KeyFields(scalars={"v": 1}))
    flight = Flight(root=str(tmp_path), units=[unit],
                        parsl_config=parsl_config)
    dispatch(flight)

    status = read_status(unit_run_dir(flight, unit))
    assert status["status"] == "done"
    assert status["detail"] == "converged"
    for field in ("submitted_at", "started_at", "finished_at",
                  "runtime_seconds"):
        assert field in status


def test_one_failure_does_not_abort_batch(tmp_path, parsl_config):
    """Complete-and-report (Principle 10): a unit that raises on
    the worker becomes ``failed`` while its siblings still reach
    ``done``.  A unit that completes-but-not-ok is also
    ``failed``, and both land in ``failures()``."""
    register_wingbeat("fake_modal", ModalRunner())
    units = [
        CalcUnit(id="ok1", structure="s", wingbeat="fake_modal",
                 options={"fake_mode": "ok"},
                 key_fields=KeyFields(scalars={"v": 1})),
        CalcUnit(id="boom", structure="s", wingbeat="fake_modal",
                 options={"fake_mode": "raise"},
                 key_fields=KeyFields(scalars={"v": 1})),
        CalcUnit(id="notok", structure="s", wingbeat="fake_modal",
                 options={"fake_mode": "not_ok"},
                 key_fields=KeyFields(scalars={"v": 1})),
    ]
    flight = Flight(root=str(tmp_path), units=units,
                        parsl_config=parsl_config)
    report = dispatch(flight)

    by_id = {e.id: e for e in report.entries}
    assert by_id["ok1"].status == "done"
    assert by_id["boom"].status == "failed"
    assert by_id["notok"].status == "failed"
    assert {e.id for e in report.failures()} == {"boom", "notok"}
    # The raised exception's message is captured for the report.
    assert "kaboom" in (by_id["boom"].message or "")


def test_resume_skips_done_units(tmp_path):
    """Re-running a flight is its resume: a unit already
    ``done`` with a still-matching key is a cache hit and the
    wingbeat is NOT called again (LocalExecutor path)."""
    wingbeat = CountingRunner()
    register_wingbeat("fake_count", wingbeat)
    unit = CalcUnit(id="u1", structure="s.skl",
                    wingbeat="fake_count",
                    key_fields=KeyFields(scalars={"v": 1}))
    flight = Flight(root=str(tmp_path), units=[unit])

    first = dispatch(flight)
    assert first.entries[0].status == "done"
    assert wingbeat.calls == 1

    second = dispatch(flight)        # resume == re-run
    assert second.entries[0].status == "done"
    assert wingbeat.calls == 1               # hit: not re-run


def test_on_outcome_callback_fires_per_unit(tmp_path):
    """The optional streaming hook is invoked once per unit with
    its terminal ReportEntry."""
    register_wingbeat("fake_ok", CountingRunner())
    seen = []
    units = [CalcUnit(id=f"u{i}", structure="s",
                      wingbeat="fake_ok",
                      key_fields=KeyFields(scalars={"v": i}))
             for i in range(3)]
    flight = Flight(root=str(tmp_path), units=units,
                        on_outcome=seen.append)
    dispatch(flight)
    assert [e.id for e in seen] == ["u0", "u1", "u2"]


# ==============================================================
#  13.2 -- ImagoWingbeat: ImagoResult -> WingbeatOutcome + result.toml
# ==============================================================

def _imago_result(status, **overrides):
    """Fabricate an ImagoResult for the mapping tests, with the
    fields ImagoWingbeat reads (status, runtime, message) and the
    ones _persist_result echoes (outputs, job)."""
    import imago
    fields = dict(
        run_dir="/r", temp_dir="/t",
        job=imago.JobIdentity("gs", "scf", "fb", "no"),
        runtime_seconds=2.5,
        outputs={"scfV": "/r/gs_scfV-fb.dat"},
        message="")
    fields.update(overrides)
    return imago.ImagoResult(status=status, **fields)


def test_imago_runner_maps_converged(tmp_path, monkeypatch):
    """A prepared run directory (holds imago.dat) is run as-is;
    a CONVERGED ImagoResult maps to ok=True, detail="converged",
    and the native result is persisted to result.toml for the
    client's harvest."""
    import imago
    (tmp_path / "imago.dat").write_text("CONVERGENCE_TEST\n 1e-4\n")
    captured = {}

    def fake_run_prepared(wingbeat_dir, **kwargs):
        captured["wingbeat_dir"] = wingbeat_dir
        return _imago_result(imago.RunStatus.CONVERGED)

    monkeypatch.setattr(imago, "run_prepared", fake_run_prepared)
    unit = CalcUnit(id="x", structure="s.skl")
    outcome = ImagoWingbeat().run(unit, str(tmp_path))

    assert captured["wingbeat_dir"] == str(tmp_path)   # prepared mode
    assert outcome.ok is True
    assert outcome.detail == "converged"
    assert outcome.runtime_seconds == 2.5
    assert (tmp_path / "result.toml").exists()


def test_imago_runner_maps_not_converged(tmp_path, monkeypatch):
    """A directory with no staged imago.dat goes through the
    structure-and-options build path; NOT_CONVERGED still
    *completed*, so it is ok=True with detail="not_converged"."""
    import imago

    def fake_run_structure(structure, options, wingbeat_dir):
        return _imago_result(imago.RunStatus.NOT_CONVERGED)

    monkeypatch.setattr(imago, "run_structure", fake_run_structure)
    unit = CalcUnit(id="x", structure="s.skl")
    outcome = ImagoWingbeat().run(unit, str(tmp_path))
    assert outcome.ok is True
    assert outcome.detail == "not_converged"


def test_imago_runner_maps_failed_to_not_ok(tmp_path, monkeypatch):
    """A hard FAILED is the only status that maps to ok=False --
    the unit did not complete."""
    import imago

    def fake_run_structure(structure, options, wingbeat_dir):
        return _imago_result(imago.RunStatus.FAILED,
                             message="fortran abort")

    monkeypatch.setattr(imago, "run_structure", fake_run_structure)
    outcome = ImagoWingbeat().run(CalcUnit(id="x", structure="s"),
                                str(tmp_path))
    assert outcome.ok is False
    assert outcome.detail == "failed"
    assert outcome.message == "fortran abort"


def test_imago_runner_prepared_detection_under_inputs(tmp_path,
                                                      monkeypatch):
    """A staged ``inputs/imago.dat`` also marks a directory as
    prepared, so the run-as-is path is taken."""
    import imago
    (tmp_path / "inputs").mkdir()
    (tmp_path / "inputs" / "imago.dat").write_text("X\n")
    used = {}

    def fake_run_prepared(wingbeat_dir, **kwargs):
        used["prepared"] = True
        return _imago_result(imago.RunStatus.CONVERGED)

    def fake_run_structure(structure, options, wingbeat_dir):
        used["structure"] = True
        return _imago_result(imago.RunStatus.CONVERGED)

    monkeypatch.setattr(imago, "run_prepared", fake_run_prepared)
    monkeypatch.setattr(imago, "run_structure", fake_run_structure)
    ImagoWingbeat().run(CalcUnit(id="x", structure="s"),
                      str(tmp_path))
    assert used == {"prepared": True}


# ==============================================================
#  13.6 -- report views the client selects on
# ==============================================================

def _entry(id, status, detail):
    """A minimal ReportEntry for the view tests."""
    return ReportEntry(id=id, calc=(), status=status,
                       detail=detail, wingbeat_dir=f"/wingbeats/{id}",
                       runtime_seconds=0.0, message="")


def test_report_views_select_correctly():
    """``by_status`` filters on the generic lifecycle status,
    ``with_detail`` on the wingbeat-supplied detail (how a client
    selects converged units), and ``failures`` collects the
    failed/lost entries."""
    report = FlightReport(entries=[
        _entry("a", "done", "converged"),
        _entry("b", "done", "not_converged"),
        _entry("c", "failed", None),
        _entry("d", "lost", None),
    ])
    assert {e.id for e in report.by_status("done")} == {"a", "b"}
    assert [e.id for e in report.with_detail("converged")] == ["a"]
    assert {e.id for e in report.failures()} == {"c", "d"}
