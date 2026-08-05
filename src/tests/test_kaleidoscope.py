## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

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
    dispatch, send_off, collect_next,
    make_executor, register_wingbeat, resolve_wingbeat,
    validate_flight, unit_run_dir,
    is_cache_hit, cache_key_matches, write_cache_key,
    read_status, write_status, ImagoWingbeat,
)
from kaleidoscope.wingbeats import Wingbeat
from kaleidoscope.workspace import (
    derive_calc_tag, serialize_flight, read_flight_toml, flight_id_of)


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
    # A flat stand-in payload (scalars plus a string array) that
    #   the dispatch core stores and round-trips without ever
    #   reading its contents.  This exercises that verbatim
    #   write-then-read cycle, so the exact keys are just example
    #   data, not the literal builder output.
    prediction = {
        "policy": "verify_around_prediction",
        "predicted_kpoint_density": 50.0,
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


def test_read_flight_toml_round_trips_units_sweep_metadata(tmp_path):
    """read_flight_toml is the disk-side inverse of
    serialize_flight: it restores each unit's identity (id,
    structure, calc tuple, wingbeat), the SweepRecord, and the
    opaque metadata tables -- the fields the harvest reads back
    (DESIGN 7.8).  options/key_fields are deliberately NOT
    persisted, so they come back empty."""
    sweep = SweepRecord(
        varied_axes=("kpt-density",),
        fixed_axes={"basis": "fb", "functional": "gga-pbe",
                    "kpoint_integration": "gaussian-0.1"})
    prediction = {
        "policy": "verify_around_prediction",
        "confidence": 0.9,
        "is_under_trained": False,
        "neighbor_entry_ids": ["mp-1", "mp-2"],
        "system_type": "crystalline",
    }
    units = [
        CalcUnit(id="si", structure="si.skl",
                 calc=("kpt-density-50",), wingbeat="imago"),
        CalcUnit(id="si", structure="si.skl",
                 calc=("kpt-density-100",), wingbeat="imago"),
    ]
    serialize_flight(Flight(root=str(tmp_path), units=units,
                            sweep=sweep,
                            metadata={"prediction": prediction}))

    flight = read_flight_toml(
        os.path.join(str(tmp_path), "flight.toml"))
    assert flight.root == str(tmp_path)
    assert [u.calc for u in flight.units] == [
        ("kpt-density-50",), ("kpt-density-100",)]
    assert all(u.id == "si" and u.structure == "si.skl"
               and u.wingbeat == "imago" for u in flight.units)
    # calc came back as a tuple (not a list) so unit_run_dir can
    #   splat it; options were not persisted.
    assert isinstance(flight.units[0].calc, tuple)
    assert flight.units[0].options == {}
    assert flight.sweep.varied_axes == ("kpt-density",)
    assert flight.sweep.fixed_axes["functional"] == "gga-pbe"
    assert flight.metadata["prediction"] == prediction


def test_calc_unit_kind_round_trips_and_defaults(tmp_path):
    """CalcUnit.kind (DESIGN 6.2.9) serializes and restores: a unit
    that does not set it defaults to 'convergence', and an explicit
    'fingerprint' loen unit round-trips as such -- this is what lets
    the convergence harvest filter out the loen runs that share a
    structure id."""
    units = [
        CalcUnit(id="si", structure="si.skl",
                 calc=("kpt-density-50",), wingbeat="imago"),
        CalcUnit(id="si", structure="si.skl", calc=("loen",),
                 wingbeat="imago", kind="fingerprint"),
    ]
    serialize_flight(Flight(root=str(tmp_path), units=units))
    flight = read_flight_toml(
        os.path.join(str(tmp_path), "flight.toml"))
    assert [u.kind for u in flight.units] == [
        "convergence", "fingerprint"]


def test_flight_id_of_is_the_workspace_basename():
    """The provenance flight_id is the workspace root's basename,
    trailing slash and all (DESIGN 7.8)."""
    assert flight_id_of("/work/flights/diamond-seed") == \
        "diamond-seed"
    assert flight_id_of("/work/flights/diamond-seed/") == \
        "diamond-seed"


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


def test_record_table_survives_every_lifecycle_rewrite(tmp_path):
    """The ``[record]`` table is stamped once, at launch, and describes
    the run rather than the lifecycle (DESIGN 6.2.4).  status.toml is
    rewritten several times after that -- running, then terminal -- and
    every rewrite must carry the table forward, or the one fact written
    to be read months later would be erased at the first transition
    after the write that put it there.

    It round-trips as a TOML sub-table, which is why ``write_status``
    holds dict-valued fields back until every scalar line is out: a
    sub-table header must follow the bare keys of the table that
    contains it, so emitting them in dict order would produce a file
    that no longer parses."""
    write_status(str(tmp_path), status="queued", submitted_at="t0",
                 record={"imago_commit": "abc123"})
    write_status(str(tmp_path), status="running", started_at="t1")
    write_status(str(tmp_path), status="done", detail="converged",
                 finished_at="t2", runtime_seconds=1.5)
    status = read_status(str(tmp_path))
    assert status["record"] == {"imago_commit": "abc123"}
    assert status["status"] == "done"
    assert status["submitted_at"] == "t0"


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
            files=[KeyFile(path="structure.skl",
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


def test_cache_misses_when_the_staged_copy_is_gone(tmp_path):
    """A run directory left half-written -- the key snapshot
    present but the staged file missing -- is a miss, not an
    error."""
    wingbeat_dir, unit = _staged_unit(
        tmp_path, "LATTICE 1 2 3\n", "LATTICE 1 2 3\n",
        {"kpoints": "4x4x4"})
    os.remove(os.path.join(wingbeat_dir, "structure.skl"))
    assert cache_key_matches(unit, wingbeat_dir) is False


def test_cache_misses_when_the_source_is_gone(tmp_path):
    """The other side of the same rule (DESIGN 6.2.5).  A source
    can vanish between runs -- a prepare directory reclaimed as
    scratch, a structure cache that moved -- and that must mean
    "re-run this unit", never an exception.

    This is the case that cost a live campaign: filecmp raises
    rather than returning False on a missing file, so an
    unguarded source turned one absent file into a crash that
    abandoned a whole flight before its first unit dispatched.
    """
    wingbeat_dir, unit = _staged_unit(
        tmp_path, "LATTICE 1 2 3\n", "LATTICE 1 2 3\n",
        {"kpoints": "4x4x4"})
    os.remove(unit.key_fields.files[0].source)
    assert cache_key_matches(unit, wingbeat_dir) is False


def _two_file_unit(tmp_path, structure_text, kpoint_text,
                   staged_kpoint_text, root_copies=True):
    """A unit keyed on TWO files, in the shape the producer declares
    (DESIGN 6.2.5): the resolved structure and the resolved k-point
    file, both named under ``inputs/`` because that is the one
    surface makeinput writes for every unit whatever its job reads.

    The directories mirror a real flight.  makeinput builds into
    ``inputs/`` on both sides, so the source is
    ``<prepare>/inputs/<name>`` and the staged copy is
    ``<run>/inputs/<name>`` -- the same relative path.  A run
    directory additionally carries the flattened root copies the
    engine reads, which is what ``root_copies`` writes.

    ``staged_kpoint_text`` may differ from ``kpoint_text`` to model a
    prior run made under a different integration scheme.  Setting
    ``root_copies=False`` models a job that reads neither file at its
    root -- the fingerprint unit of TODO D23."""
    wingbeat_dir = tmp_path / "run"
    (wingbeat_dir / "inputs").mkdir(parents=True)
    sources = tmp_path / "prepare"
    (sources / "inputs").mkdir(parents=True)
    files = []
    for name, current, staged in (
            ("structure.dat", structure_text, structure_text),
            ("kp-scf.dat", kpoint_text, staged_kpoint_text)):
        (sources / "inputs" / name).write_text(current)
        (wingbeat_dir / "inputs" / name).write_text(staged)
        if root_copies:
            (wingbeat_dir / name).write_text(staged)
        files.append(KeyFile(
            path=f"inputs/{name}",
            source=str(sources / "inputs" / name)))
    unit = CalcUnit(id="s1", structure=str(sources),
                    key_fields=KeyFields(scalars={}, files=files))
    write_cache_key(str(wingbeat_dir), unit)
    write_status(str(wingbeat_dir), status="done")
    return str(wingbeat_dir), unit


def test_a_changed_integration_scheme_misses(tmp_path):
    """The fault C135 closes.  The k-point integration scheme lives
    in ``kp-scf.dat`` as KPOINT_INTG_CODE and appears nowhere in
    ``structure.dat``, so a key naming only the structure would call
    these the same calculation and return the stored answer.

    That is a different failure from a stale hit: it reports the
    physics of one integration scheme under the name of another, and
    prints nothing.  Declaring the k-point file as a second key file
    is what makes the two distinguishable."""
    wingbeat_dir, unit = _two_file_unit(
        tmp_path, "cell 1 2 3\n",
        "KPOINT_INTG_CODE\n1\n",          # tetrahedron now
        "KPOINT_INTG_CODE\n0\n")          # histogram before
    assert is_cache_hit(unit, wingbeat_dir) is False


def test_the_same_scheme_still_hits(tmp_path):
    """The property that decided key FILE over key scalar.  Adding a
    name to the scalar table would mismatch every stored
    ``cache_key.toml`` at once, since the scalars are compared whole
    -- a mass false miss, the failure with no escape valve.  A key
    file costs nothing provided its path names a file every unit
    has, which is what ``inputs/`` guarantees, so an unchanged
    scheme reuses its result exactly as before."""
    wingbeat_dir, unit = _two_file_unit(
        tmp_path, "cell 1 2 3\n",
        "KPOINT_INTG_CODE\n0\n", "KPOINT_INTG_CODE\n0\n")
    assert is_cache_hit(unit, wingbeat_dir) is True


def test_a_unit_whose_job_stages_no_root_copy_still_hits(tmp_path):
    """The fault D23 closes.  A run directory carries flattened root
    copies only of the names its own job reads, so a fingerprint
    unit -- which runs no SCF -- has no ``kp-scf.dat`` at its root
    however many times it has been computed.

    Keyed on the root, such a unit can never match: it is not stale,
    it is uncacheable, and it recomputes on every campaign at full
    price while the log says only 'no usable result'.  Keyed under
    ``inputs/``, which makeinput writes for every unit, the same
    directory is recognised."""
    wingbeat_dir, unit = _two_file_unit(
        tmp_path, "cell 1 2 3\n",
        "KPOINT_INTG_CODE\n0\n", "KPOINT_INTG_CODE\n0\n",
        root_copies=False)
    assert is_cache_hit(unit, wingbeat_dir) is True


def test_declaring_a_bare_name_is_what_made_it_uncacheable(tmp_path):
    """The same directory, judged under the two declarations, so the
    defect and its fix are visible side by side rather than asserted.

    Nothing about the stored run changes here.  Only the declared
    path changes -- a bare ``kp-scf.dat``, which resolves to the
    run-directory root, against ``inputs/kp-scf.dat``, which resolves
    to the surface every unit has.  The first cannot match and the
    second does, which is the whole of D23 in one comparison."""
    wingbeat_dir, unit = _two_file_unit(
        tmp_path, "cell 1 2 3\n",
        "KPOINT_INTG_CODE\n0\n", "KPOINT_INTG_CODE\n0\n",
        root_copies=False)
    assert is_cache_hit(unit, wingbeat_dir) is True

    for key_file in unit.key_fields.files:
        key_file.path = os.path.basename(key_file.path)
    assert is_cache_hit(unit, wingbeat_dir) is False


def test_a_root_copy_that_disagrees_misses(tmp_path):
    """The agreement test, and the backstop under C139.

    The root copy is the file the engine actually reads.  A commit
    that refreshed ``inputs/`` while leaving a previous
    calculation's root copy in place would run the OLD physics while
    the key file, the run's summary and the flight report all
    described the new -- silently, and on a hit, so the wrong answer
    would be returned for free rather than merely recomputed.

    C139's rule stops that arising by clearing the root copies a
    commit supersedes.  This test is why the cache does not simply
    rely on it: the property that the key describes what the engine
    reads is checked where a reader checking the cache can see it,
    and it holds for a directory built before that rule existed."""
    wingbeat_dir, unit = _two_file_unit(
        tmp_path, "cell 1 2 3\n",
        "KPOINT_INTG_CODE\n0\n", "KPOINT_INTG_CODE\n0\n")
    assert is_cache_hit(unit, wingbeat_dir) is True

    # The staged inputs stay correct; only the flattened copy the
    #   engine would read goes stale.
    with open(os.path.join(wingbeat_dir, "kp-scf.dat"), "w") as f:
        f.write("KPOINT_INTG_CODE\n1\n")
    assert is_cache_hit(unit, wingbeat_dir) is False


def test_a_same_size_rewrite_is_caught_despite_the_memo(tmp_path):
    """The comparison must read both files, every call.

    ``filecmp`` memoizes on a (mode, size, mtime) signature, which
    is blind to exactly the change this cache exists to catch: a
    rewrite in place that keeps the length.  mtime resolution is
    coarse enough that two writes microseconds apart routinely
    share one tick, so the memo answers "equal" for files that
    differ -- a false HIT, returning a stored result for inputs
    that have changed.  That is the failure with no escape valve:
    ``--force`` turns hits into misses, but nobody knows to reach
    for it.

    The test above meets this by luck, since it depends on whether
    two writes happened to land in the same tick.  This one removes
    the luck: the mtime is restored explicitly, so the memoized
    signature is guaranteed identical and only a real re-read of
    the bytes can answer correctly."""
    wingbeat_dir, unit = _two_file_unit(
        tmp_path, "cell 1 2 3\n",
        "KPOINT_INTG_CODE\n0\n", "KPOINT_INTG_CODE\n0\n")
    assert is_cache_hit(unit, wingbeat_dir) is True   # warms it

    root_copy = os.path.join(wingbeat_dir, "kp-scf.dat")
    was = os.stat(root_copy)
    with open(root_copy, "w") as handle:
        handle.write("KPOINT_INTG_CODE\n1\n")         # same length
    os.utime(root_copy, ns=(was.st_atime_ns, was.st_mtime_ns))

    assert os.path.getsize(root_copy) == was.st_size
    assert os.stat(root_copy).st_mtime_ns == was.st_mtime_ns
    assert is_cache_hit(unit, wingbeat_dir) is False


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


def test_force_bypasses_cache_and_reruns(tmp_path):
    """``force=True`` bypasses the run-reuse cache (DESIGN 6.2.5):
    a unit already ``done`` with a still-matching key is NOT
    treated as a hit, so the wingbeat runs again.  The switch
    lives on the driver because the cache it governs does."""
    wingbeat = CountingRunner()
    register_wingbeat("fake_force", wingbeat)
    unit = CalcUnit(id="u1", structure="s.skl",
                    wingbeat="fake_force",
                    key_fields=KeyFields(scalars={"v": 1}))
    flight = Flight(root=str(tmp_path), units=[unit])

    first = dispatch(flight)
    assert first.entries[0].status == "done"
    assert wingbeat.calls == 1

    # Without force this second run would be a cache hit (the
    #   resume test above); force makes it re-run regardless.
    second = dispatch(flight, force=True)
    assert second.entries[0].status == "done"
    assert wingbeat.calls == 2               # forced: re-run


# ==============================================================
#  13.5 -- the recorded build, and the key that no longer holds
#  it (DESIGN 6.2.5).  The cache asks whether this is the same
#  CALCULATION, not whether its result is still good, so the
#  engine build rides on the unit's ``record``: written down at
#  launch, printed in the reuse plan, never compared.
# ==============================================================

def _recorded_build_flight(root, wingbeat_name, build):
    """A one-unit flight whose recorded build is ``build``.  The cache
    identity is identical in every call, so two flights differing only
    in what they record are the same calculation as far as the key is
    concerned -- which is precisely what the tests below check."""
    unit = CalcUnit(id="u1", structure="s.skl",
                    wingbeat=wingbeat_name,
                    key_fields=KeyFields(scalars={"v": 1}),
                    record={"imago_commit": build})
    return Flight(root=str(root), units=[unit])


def test_a_rebuilt_engine_is_still_a_cache_hit(tmp_path):
    """The whole point of the change (DESIGN 6.2.5; VISION 16).  A
    second run whose recorded build differs -- an ordinary development
    commit, which happens constantly and changes the physics almost
    never -- must still HIT.  A rebuilt engine does not make a stored
    potential wrong: that potential is a starting point every later SCF
    re-converges.

    The two ways of being wrong are not symmetric, which is why the old
    guard was dropped rather than tightened.  A false hit has an escape
    valve in ``force``; a false miss has none, and the hours are simply
    spent again."""
    wingbeat = CountingRunner()
    register_wingbeat("fake_build", wingbeat)
    dispatch(_recorded_build_flight(tmp_path, "fake_build", "old111"))
    assert wingbeat.calls == 1

    report = dispatch(
        _recorded_build_flight(tmp_path, "fake_build", "new222"))
    assert report.entries[0].status == "done"
    assert wingbeat.calls == 1                # hit: not re-run


def test_a_hit_keeps_the_record_of_the_run_that_produced_it(tmp_path):
    """``record`` is written on the MISS only, so a later hit leaves it
    describing the run that produced the stored result rather than the
    flight that reused it (DESIGN 6.2.4).  That is what makes it worth
    printing: the build named beside a reused result is the build that
    result actually came out of."""
    register_wingbeat("fake_keep", CountingRunner())
    dispatch(_recorded_build_flight(tmp_path, "fake_keep", "old111"))
    dispatch(_recorded_build_flight(tmp_path, "fake_keep", "new222"))
    status = read_status(str(tmp_path / "wingbeats" / "u1"))
    assert status["record"] == {"imago_commit": "old111"}


def test_a_unit_that_records_nothing_writes_no_record_table(tmp_path):
    """The mapping is free-form and optional: a client that hangs
    nothing on its units gets no bare ``[record]`` header at all
    (DESIGN 6.2.4)."""
    register_wingbeat("fake_bare", CountingRunner())
    unit = CalcUnit(id="u1", structure="s.skl", wingbeat="fake_bare",
                    key_fields=KeyFields(scalars={"v": 1}))
    dispatch(Flight(root=str(tmp_path), units=[unit]))
    status = read_status(str(tmp_path / "wingbeats" / "u1"))
    assert "record" not in status


# ==============================================================
#  13.5 -- the reuse plan and the preview (DESIGN 6.2.5).  This
#  is what stands in for the automatic staleness guard: the
#  build behind a reused result is REPORTED to a person who can
#  act on it, instead of being silently compared.
# ==============================================================

def test_reuse_plan_names_the_hits_the_misses_and_the_reason(tmp_path):
    """The plan is decided from local files and touches nothing.  A
    unit with a completed run directory reuses, and carries the facts a
    judgment would want -- when the result finished and the build
    recorded behind it.  Under ``force`` the same unit runs, and the
    plan says *why* rather than leaving a reader to guess."""
    from kaleidoscope import reuse_plan
    register_wingbeat("fake_plan", CountingRunner())
    flight = _recorded_build_flight(tmp_path, "fake_plan", "old111")
    dispatch(flight)

    (_, action, detail), = reuse_plan(flight, flight.units)
    assert action == "reuse"
    assert detail["record"] == {"imago_commit": "old111"}
    assert detail["finished_at"]

    (_, forced, why), = reuse_plan(flight, flight.units, force=True)
    assert forced == "run"
    assert why["reason"] == "forced"


def test_the_counts_always_print_and_the_lines_wait_to_be_asked(
        tmp_path, capsys):
    """The counts are the decision being announced, so they always
    print.  The per-unit lines are that decision's evidence and are
    held back until asked for (DESIGN 5.7 / 6.2.5): the climb calls
    ``send_off`` once per round, so an unconditional line per unit
    would refill the screen the reporting rule cleared, on the very
    path it cleared."""
    from kaleidoscope import print_reuse_plan, reuse_plan
    register_wingbeat("fake_quiet", CountingRunner())
    flight = _recorded_build_flight(tmp_path, "fake_quiet", "old111")
    dispatch(flight)
    capsys.readouterr()                  # discard the first run's own

    plan = reuse_plan(flight, flight.units)
    print_reuse_plan(plan)
    quiet = capsys.readouterr().out
    assert "1 to reuse, 0 to run" in quiet
    assert "u1" not in quiet

    print_reuse_plan(plan, per_unit=True)
    loud = capsys.readouterr().out
    assert "1 to reuse, 0 to run" in loud
    assert "u1" in loud and "old111" in loud


def test_send_off_takes_its_per_unit_lines_from_the_switch(
        tmp_path, capsys):
    """Verbosity is a module-level switch a client sets once from its
    entry point, NOT an argument threaded through ``send_off``: it
    describes how the process talks to its user, not how a flight
    dispatches, so threading it would put a reporting concern into the
    signature of every function between a client's main and the printer
    (PSEUDOCODE 13.5)."""
    from kaleidoscope import set_verbose
    register_wingbeat("fake_switch", CountingRunner())
    flight = _recorded_build_flight(tmp_path, "fake_switch", "old111")
    executor = make_executor(None)
    try:
        set_verbose(True)
        send_off(flight, flight.units, executor)
        assert "u1" in capsys.readouterr().out

        set_verbose(False)
        send_off(flight, flight.units, executor)
        assert "u1" not in capsys.readouterr().out
    finally:
        set_verbose(False)               # never leak the switch
        executor.close()


def test_preview_prints_the_lines_spends_nothing_and_reports_nothing(
        tmp_path, capsys):
    """A preview answers "what will this cost?" BEFORE a flight starts
    rather than while it goes past (DESIGN 6.2.5).  No executor is
    built and no unit runs; the per-unit lines print whether or not
    verbosity is on, since reading them one by one is the whole purpose
    of a preview; and the report comes back EMPTY rather than partial,
    because nothing ran to report on."""
    from kaleidoscope import set_verbose
    wingbeat = CountingRunner()
    register_wingbeat("fake_preview", wingbeat)
    flight = _recorded_build_flight(tmp_path, "fake_preview", "old111")
    set_verbose(False)

    report = dispatch(flight, preview=True)

    printed = capsys.readouterr().out
    assert "u1" in printed                    # the lines, not just
    assert "0 to reuse, 1 to run" in printed  #   the counts
    assert wingbeat.calls == 0                # nothing dispatched
    assert report.entries == []
    # Nothing was touched either: a preview builds no run directory.
    assert not os.path.exists(str(tmp_path / "wingbeats"))


def test_shared_executor_serves_repeated_dispatches(tmp_path,
                                                    parsl_config):
    """The producer's climb dispatches MANY flights under one
    executor (a pre-flight batch, then one flight per round): it
    builds the executor once with make_executor and pins it to
    every dispatch, so a Parsl config -- single-use once closed --
    is never reloaded (PSEUDOCODE 13.5).  A caller-pinned executor
    must survive one dispatch and serve the next; the caller, not
    dispatch, closes it.  This guards the multi-round regression
    that reloading the same config triggered."""
    register_wingbeat("fake_ok", CountingRunner())
    executor = make_executor(parsl_config)       # one warm pool
    try:
        # Two rounds, distinct unit ids -> distinct run dirs, so the
        #   second is real work, not a cache hit: it can only reach
        #   `done` if the pinned executor was NOT torn down after the
        #   first dispatch.
        for round_id in ("round1", "round2"):
            unit = CalcUnit(
                id=round_id, structure="s.skl", wingbeat="fake_ok",
                key_fields=KeyFields(scalars={"v": 1}))
            flight = Flight(root=str(tmp_path), units=[unit],
                            parsl_config=parsl_config)
            report = dispatch(flight, executor=executor)
            assert report.entries[0].status == "done"
    finally:
        executor.close()


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
#  13.5 -- the two public phases: send_off / collect_next
#  (DESIGN 6.2.3).  These are what a control-loop client (the
#  k-point climb) drives directly instead of the one-shot
#  dispatch wrapper.
# ==============================================================

def test_send_off_then_collect_next_drains_all(tmp_path,
                                               parsl_config):
    """The two phases compose: send_off launches the units and returns
    one (unit, future) pair each without waiting, and collect_next
    takes whichever has landed until the outstanding list empties.
    Every unit reaches the terminal 'done' status, on both
    executors."""
    register_wingbeat("fake_ok", CountingRunner())
    units = [CalcUnit(id=f"u{i}", structure="s", wingbeat="fake_ok",
                      key_fields=KeyFields(scalars={"v": i}))
             for i in range(3)]
    flight = Flight(root=str(tmp_path), units=units,
                    parsl_config=parsl_config)
    executor = make_executor(parsl_config)
    try:
        outstanding = send_off(flight, flight.units, executor,
                               force=False)
        assert len(outstanding) == 3

        collected = []
        while outstanding:
            unit, entry, outstanding = collect_next(
                flight, outstanding)
            collected.append((unit.id, entry.status))
    finally:
        executor.close()

    assert sorted(collected) == [("u0", "done"), ("u1", "done"),
                                 ("u2", "done")]


def test_send_off_futures_are_done_on_the_local_executor(tmp_path):
    """The local executor runs each unit synchronously in send_off, so
    every future it returns is already done() -- which is what lets a
    local climb never reach collect_next's poll sleep."""
    register_wingbeat("fake_ok", CountingRunner())
    units = [CalcUnit(id=f"u{i}", structure="s", wingbeat="fake_ok",
                      key_fields=KeyFields(scalars={"v": i}))
             for i in range(2)]
    flight = Flight(root=str(tmp_path), units=units)
    executor = make_executor(None)              # LocalExecutor
    try:
        outstanding = send_off(flight, flight.units, executor,
                               force=False)
        assert all(future.done() for _unit, future in outstanding)
    finally:
        executor.close()


def test_send_off_returns_a_cache_hit_as_a_done_future(tmp_path):
    """A unit already 'done' with a matching key is a cache hit:
    send_off submits no task and hands back an already-done future, so
    collect_next reports it from the existing status.toml without
    re-running the wingbeat (hits and misses sit uniformly in the
    outstanding set)."""
    wingbeat = CountingRunner()
    register_wingbeat("fake_count", wingbeat)
    unit = CalcUnit(id="u1", structure="s.skl", wingbeat="fake_count",
                    key_fields=KeyFields(scalars={"v": 1}))
    flight = Flight(root=str(tmp_path), units=[unit])
    dispatch(flight)                            # first run fills cache
    assert wingbeat.calls == 1

    executor = make_executor(None)
    try:
        outstanding = send_off(flight, flight.units, executor,
                               force=False)
        assert outstanding[0][1].done() is True
        _unit, entry, _rest = collect_next(flight, outstanding)
    finally:
        executor.close()
    assert entry.status == "done"
    assert wingbeat.calls == 1                  # hit: not re-run


def test_collect_next_streams_the_outcome_hook(tmp_path):
    """collect_next fires the flight's on_outcome hook as each unit is
    collected (in landing order), so a control-loop consumer sees each
    outcome the moment it lands rather than after the whole batch."""
    register_wingbeat("fake_ok", CountingRunner())
    seen = []
    units = [CalcUnit(id=f"u{i}", structure="s", wingbeat="fake_ok",
                      key_fields=KeyFields(scalars={"v": i}))
             for i in range(3)]
    flight = Flight(root=str(tmp_path), units=units,
                    on_outcome=seen.append)
    executor = make_executor(None)
    try:
        outstanding = send_off(flight, flight.units, executor,
                               force=False)
        while outstanding:
            _unit, _entry, outstanding = collect_next(
                flight, outstanding)
    finally:
        executor.close()
    assert sorted(e.id for e in seen) == ["u0", "u1", "u2"]


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
    structure-and-options build path (the wingbeat builds the deck
    with makeinput, then runs it); NOT_CONVERGED still *completed*,
    so it is ok=True with detail="not_converged"."""
    import imago
    import makeinput

    monkeypatch.setattr(makeinput, "build_run_dir",
                        lambda *a, **k: None)
    monkeypatch.setattr(
        imago, "run_prepared",
        lambda wingbeat_dir, **kwargs: _imago_result(
            imago.RunStatus.NOT_CONVERGED))
    unit = CalcUnit(id="x", structure="s.skl")
    outcome = ImagoWingbeat().run(unit, str(tmp_path))
    assert outcome.ok is True
    assert outcome.detail == "not_converged"


def test_imago_runner_maps_failed_to_not_ok(tmp_path, monkeypatch):
    """A hard FAILED is the only status that maps to ok=False --
    the unit did not complete."""
    import imago
    import makeinput

    monkeypatch.setattr(makeinput, "build_run_dir",
                        lambda *a, **k: None)
    monkeypatch.setattr(
        imago, "run_prepared",
        lambda wingbeat_dir, **kwargs: _imago_result(
            imago.RunStatus.FAILED, message="fortran abort"))
    outcome = ImagoWingbeat().run(CalcUnit(id="x", structure="s"),
                                str(tmp_path))
    assert outcome.ok is False
    assert outcome.detail == "failed"
    assert outcome.message == "fortran abort"


def test_persist_result_carries_resolved_mesh(tmp_path):
    """The resolved mesh reaches result.toml (DESIGN 6.1.2):
    kpoint_mesh as a TOML array, kpoint_count as a scalar, so the
    k-density guard (PSEUDOCODE 15.7) can read them back."""
    import tomllib
    import imago
    result = _imago_result(
        imago.RunStatus.CONVERGED,
        total_energy=-31.1, kpoint_mesh=[4, 2, 3], kpoint_count=24)
    ImagoWingbeat._persist_result(str(tmp_path), result)
    with open(tmp_path / "result.toml", "rb") as handle:
        loaded = tomllib.load(handle)
    assert loaded["kpoint_mesh"] == [4, 2, 3]
    assert loaded["kpoint_count"] == 24


def test_persist_result_omits_absent_mesh(tmp_path):
    """An explicit-list run or older binary emits no mesh; the
    None fields are omitted from result.toml (so the guard reads
    them back as absent and stays inert, PSEUDOCODE 15.7)."""
    import tomllib
    import imago
    result = _imago_result(imago.RunStatus.CONVERGED)
    ImagoWingbeat._persist_result(str(tmp_path), result)
    with open(tmp_path / "result.toml", "rb") as handle:
        loaded = tomllib.load(handle)
    assert "kpoint_mesh" not in loaded
    assert "kpoint_count" not in loaded


def test_the_wingbeat_echoes_the_recorded_build(tmp_path, monkeypatch):
    """One RECORDED fact rides into result.toml beside the measured
    ones: the build identity out of the unit's ``record`` mapping
    (DESIGN 6.2.2 / 6.2.4).  A guidance entry's provenance reads it
    there, which is what keeps that harvest on the three per-run
    sources it already has instead of opening the dispatch core's
    status.toml.  The fact lands in two files that cannot disagree,
    because both are copied from the same mapping at launch."""
    import imago
    (tmp_path / "imago.dat").write_text("X\n")
    monkeypatch.setattr(
        imago, "run_prepared",
        lambda wingbeat_dir, **kwargs: _imago_result(
            imago.RunStatus.CONVERGED))
    unit = CalcUnit(id="x", structure="s.skl",
                    record={"imago_commit": "abc123"})
    ImagoWingbeat().run(unit, str(tmp_path))
    with open(tmp_path / "result.toml", "rb") as handle:
        loaded = tomllib.load(handle)
    assert loaded["imago_commit"] == "abc123"


def test_persist_result_prefers_the_build_the_engine_reported(tmp_path):
    """The seam TODO C84 lands in.  The recorded value is what the
    *producer believed* it launched, which can drift from the binary
    that actually ran; the engine's own word is worth strictly more.
    imago does not report its build yet, so the echo is written as a
    FALLBACK -- once it does, preferring it is the whole change, one
    substitution in one field of one file rather than new plumbing."""
    import imago
    result = _imago_result(imago.RunStatus.CONVERGED)
    result.imago_commit = "from-the-binary"
    ImagoWingbeat._persist_result(str(tmp_path), result,
                                  {"imago_commit": "what-we-thought"})
    with open(tmp_path / "result.toml", "rb") as handle:
        loaded = tomllib.load(handle)
    assert loaded["imago_commit"] == "from-the-binary"


def test_persist_result_omits_an_unrecorded_build(tmp_path):
    """A client that hangs nothing on its units writes no
    ``imago_commit`` line at all; the harvest's ``"unknown"`` floor
    covers that case, and it stays non-empty so the schema's rule-11
    check passes and a curator can spot it on review (DESIGN 7.8)."""
    import imago
    ImagoWingbeat._persist_result(
        str(tmp_path), _imago_result(imago.RunStatus.CONVERGED))
    with open(tmp_path / "result.toml", "rb") as handle:
        loaded = tomllib.load(handle)
    assert "imago_commit" not in loaded


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

    def fake_build_run_dir(structure, options, wingbeat_dir):
        # The structure-build path would call this; in prepared
        #   mode it must be skipped entirely.
        used["structure"] = True

    import makeinput
    monkeypatch.setattr(imago, "run_prepared", fake_run_prepared)
    monkeypatch.setattr(makeinput, "build_run_dir",
                        fake_build_run_dir)
    ImagoWingbeat().run(CalcUnit(id="x", structure="s"),
                      str(tmp_path))
    assert used == {"prepared": True}


def test_imago_runner_prepared_reapplies_imago_settings(tmp_path,
                                                        monkeypatch):
    """C110: a re-run of an already-prepared directory must STILL
    receive the unit's imago-side settings.  The job type and the
    SCF suppression live only in these settings (not in the staged
    imago.dat, DESIGN 6.2.10), so if they are dropped imago no
    longer sees the unit's ``-loen -scf no`` request and falls back
    to its default job, a ground-state SCF ("SCF after loen").  The
    prepared path must build the settings from the imago-side
    options and pass them to run_prepared, as the build path does."""
    import imago
    (tmp_path / "imago.dat").write_text("X\n")   # prepared directory
    captured = {}
    settings_sentinel = object()

    def fake_from_options(options):
        captured["imago_options"] = options
        return settings_sentinel

    def fake_run_prepared(wingbeat_dir, settings=None, **kwargs):
        captured["settings"] = settings
        return _imago_result(imago.RunStatus.CONVERGED)

    monkeypatch.setattr(imago.ScriptSettings, "from_options",
                        fake_from_options)
    monkeypatch.setattr(imago, "run_prepared", fake_run_prepared)

    unit = CalcUnit(id="x", structure="s.skl",
                    options={"job": "loen", "scf_basis": "no"})
    ImagoWingbeat().run(unit, str(tmp_path))

    # The built settings reached run_prepared (not the pre-fix
    #   None), and they were built from the imago-side options only.
    assert captured["settings"] is settings_sentinel
    assert captured["imago_options"] == {"job": "loen",
                                         "scf_basis": "no"}


def test_imago_runner_commits_prepared_inputs(tmp_path, monkeypatch):
    """Model A (C111): when the driver prepared the unit's inputs
    into unit.prepared_dir, the wingbeat COMMITS them into the run
    dir (no rebuild) and runs the prepared directory.  makeinput is
    never called on this path, and the imago-side settings are still
    passed to run_prepared."""
    import imago
    import makeinput

    staging = tmp_path / "staging"
    staging.mkdir()
    (staging / "structure.dat").write_text("STRUCT\n")
    (staging / "imago.dat").write_text("DECK\n")
    wingbeat_dir = tmp_path / "run"
    used = {}

    def fake_run_prepared(wb_dir, settings=None, **kwargs):
        used["settings"] = settings
        return _imago_result(imago.RunStatus.CONVERGED)

    def boom_build(*a, **k):
        used["built"] = True   # must NOT happen on the commit path

    monkeypatch.setattr(imago, "run_prepared", fake_run_prepared)
    monkeypatch.setattr(makeinput, "build_run_dir", boom_build)

    unit = CalcUnit(id="x", structure="s.skl",
                    options={"job": "loen", "scf_basis": "no"},
                    prepared_dir=str(staging))
    ImagoWingbeat().run(unit, str(wingbeat_dir))

    # The staged inputs were committed into the run dir, the deck
    #   was NOT rebuilt, and the settings still reached run_prepared.
    assert (wingbeat_dir / "structure.dat").read_text() == "STRUCT\n"
    assert (wingbeat_dir / "imago.dat").read_text() == "DECK\n"
    assert "built" not in used
    assert used["settings"] is not None


def _staged_dir_with_inputs(root, **files):
    """Build a prepare-style staging directory: makeinput's outputs
    under ``inputs/`` and nothing flattened at the root, which is the
    shape ``prepare_units`` produces (DESIGN 6.2.5)."""
    root.mkdir(parents=True, exist_ok=True)
    inputs = root / "inputs"
    inputs.mkdir()
    for name, text in files.items():
        (inputs / name).write_text(text)
    return root


def test_commit_clears_superseded_root_copies(tmp_path, monkeypatch):
    """A commit onto a SURVIVING run directory must not leave the
    previous calculation's flattened root copies in place (DESIGN
    6.2.5, "What a commit owes a surviving run directory").

    imago.py reads the root copy and only ever creates one when it is
    absent, so a root ``kp-scf.dat`` left behind means the engine runs
    the OLD physics on a cache MISS -- the compute paid in full for
    the previous answer, with nothing printed.  The commit therefore
    REMOVES each root copy it supersedes and lets imago.py refill it
    from the staged file."""
    import imago
    import makeinput

    staging = _staged_dir_with_inputs(
        tmp_path / "staging",
        **{"kp-scf.dat": "KPOINT_INTG_CODE\n0\n",     # gaussian now
           "structure.dat": "STRUCT-NEW\n"})
    wingbeat_dir = tmp_path / "run"
    (wingbeat_dir / "inputs").mkdir(parents=True)
    # The prior calculation's flattened copies, tetrahedral.
    (wingbeat_dir / "kp-scf.dat").write_text(
        "KPOINT_INTG_CODE\n1\n")
    (wingbeat_dir / "structure.dat").write_text("STRUCT-OLD\n")

    monkeypatch.setattr(
        imago, "run_prepared",
        lambda wb_dir, settings=None, **kw: _imago_result(
            imago.RunStatus.CONVERGED))
    monkeypatch.setattr(makeinput, "build_run_dir",
                        lambda *a, **k: None)

    unit = CalcUnit(id="x", structure="s.skl", options={},
                    prepared_dir=str(staging))
    ImagoWingbeat().run(unit, str(wingbeat_dir))

    # The stale root copies are GONE, not merely still stale, so
    #   imago.py's copy-up refills them from the staged files.
    assert not (wingbeat_dir / "kp-scf.dat").exists()
    assert not (wingbeat_dir / "structure.dat").exists()
    # And the authoritative staged copies did land under inputs/.
    assert (wingbeat_dir / "inputs" / "kp-scf.dat").read_text() == (
        "KPOINT_INTG_CODE\n0\n")
    assert (wingbeat_dir / "inputs"
            / "structure.dat").read_text() == "STRUCT-NEW\n"


def test_commit_keeps_a_prior_runs_outputs(tmp_path, monkeypatch):
    """The clearing is scoped to the staged INPUT names, so a prior
    run's outputs survive -- above all the converged potential, which
    is a starting point every later SCF re-converges (DESIGN 6.2.5).

    It survives by construction rather than by a carve-out: an output
    name is simply absent from the staged ``inputs/`` listing.  This
    is what separates the rule from "wipe the run directory", which
    would also destroy DESIGN 6.1's within-directory checkpointing."""
    import imago
    import makeinput

    staging = _staged_dir_with_inputs(
        tmp_path / "staging",
        **{"kp-scf.dat": "KPOINT_INTG_CODE\n0\n",
           "scfV.dat": "DATABASE-POTENTIAL\n"})
    wingbeat_dir = tmp_path / "run"
    wingbeat_dir.mkdir()
    (wingbeat_dir / "kp-scf.dat").write_text(
        "KPOINT_INTG_CODE\n1\n")            # a staged name: cleared
    (wingbeat_dir / "gs_scfV-fb.dat").write_text(
        "CONVERGED-POTENTIAL\n")            # an output name: kept
    (wingbeat_dir / "gs_scf-fb.out").write_text("LOG\n")
    (wingbeat_dir / "fort.15").write_text("UNIT\n")

    monkeypatch.setattr(
        imago, "run_prepared",
        lambda wb_dir, settings=None, **kw: _imago_result(
            imago.RunStatus.CONVERGED))
    monkeypatch.setattr(makeinput, "build_run_dir",
                        lambda *a, **k: None)

    unit = CalcUnit(id="x", structure="s.skl", options={},
                    prepared_dir=str(staging))
    ImagoWingbeat().run(unit, str(wingbeat_dir))

    assert not (wingbeat_dir / "kp-scf.dat").exists()
    assert (wingbeat_dir / "gs_scfV-fb.dat").read_text() == (
        "CONVERGED-POTENTIAL\n")
    assert (wingbeat_dir / "gs_scf-fb.out").exists()
    assert (wingbeat_dir / "fort.15").exists()


def test_commit_onto_a_clean_run_dir_removes_nothing(
        tmp_path, monkeypatch):
    """The first run over a clean workspace has no root copies to
    clear, so the pass removes nothing and raises nothing.  A staged
    name with no root copy is the ordinary case, not an error."""
    import imago
    import makeinput

    staging = _staged_dir_with_inputs(
        tmp_path / "staging",
        **{"kp-scf.dat": "KPOINT_INTG_CODE\n0\n"})
    wingbeat_dir = tmp_path / "run"

    monkeypatch.setattr(
        imago, "run_prepared",
        lambda wb_dir, settings=None, **kw: _imago_result(
            imago.RunStatus.CONVERGED))
    monkeypatch.setattr(makeinput, "build_run_dir",
                        lambda *a, **k: None)

    unit = CalcUnit(id="x", structure="s.skl", options={},
                    prepared_dir=str(staging))
    ImagoWingbeat().run(unit, str(wingbeat_dir))

    assert (wingbeat_dir / "inputs" / "kp-scf.dat").read_text() == (
        "KPOINT_INTG_CODE\n0\n")


def test_partition_options_routes_by_recognised_key_set():
    """The wingbeat splits a unit's options TWO ways and no more
    (DESIGN 6.2.10): imago run-time selections go to imago, and
    everything else goes to the strict makeinput build.  There is no
    third "dropped before forwarding" bucket -- every key in
    ``options`` is a real tool input."""
    from kaleidoscope.wingbeats import _partition_options
    options = {
        "scf_basis": "fb",         # imago run-time selection
        "job": "scf",              # imago
        "xccode": 100,             # makeinput
        "scfkpint": 1,             # makeinput
        "converg": 1.0e-6,         # makeinput
    }
    makeinput_options, imago_options = _partition_options(options)
    assert imago_options == {"scf_basis": "fb", "job": "scf"}
    assert makeinput_options == {
        "xccode": 100, "scfkpint": 1, "converg": 1.0e-6}
    # Every key landed somewhere: nothing is silently swallowed, which
    #   is what keeps makeinput's strict check a pure typo backstop.
    assert (set(makeinput_options) | set(imago_options)
            == set(options))


def test_bookkeeping_is_not_an_option():
    """A fact ABOUT a run does not ride in ``options`` (DESIGN
    6.2.10).  With no third bucket to swallow it, an unrecognised key
    now falls through to makeinput -- which is the point: it is
    forwarded so makeinput's strict check can name it as the typo it
    almost certainly is, rather than being dropped in silence."""
    from kaleidoscope.wingbeats import _partition_options
    makeinput_options, imago_options = _partition_options(
        {"imago_commit": "abc123"})
    assert imago_options == {}
    assert makeinput_options == {"imago_commit": "abc123"}


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
