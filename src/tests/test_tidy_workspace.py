"""test_tidy_workspace.py -- unit tests for the scratch reclamation
tool (``src/scripts/tidy_workspace.py``; DESIGN 6.2.12, PSEUDOCODE
13.8).

The tool deletes, so most of what is worth testing is what it
REFUSES to delete.  Three refusals define it -- only scratch is
removed and never the run directory, an unfinished unit is never
touched, and a link pointing outside the scratch area is never
followed -- and each gets a dedicated test here, alongside the
selection filters and the dry-run-by-default rule.

The fixtures build a miniature workspace on tmp_path: a
``wingbeats/<id>/<calc>/`` run directory holding the kept tier
(``status.toml``, ``result.toml``) plus an ``intermediate`` symlink
into a separate "scratch" tree.  That mirrors the real layout
closely enough to exercise every path without needing a flight.
"""

import os

import pytest

import tidy_workspace


pytestmark = pytest.mark.unit


# ==============================================================
#  Fixture builders
# ==============================================================

def _make_run(root, scratch, unit_id, calc, *, status="done",
              result=True, scratch_bytes=1024, link_to=None):
    """Create one run directory and its scratch, and return the
    run directory path.

    ``status`` and ``result`` control the two conditions the default
    policy tests; ``link_to`` overrides where ``intermediate``
    points, so a test can aim it outside the scratch area.
    """

    run_dir = os.path.join(root, "wingbeats", unit_id, calc)
    os.makedirs(run_dir, exist_ok=True)

    if status is not None:
        with open(os.path.join(run_dir, "status.toml"), "w") as f:
            f.write(f'status = "{status}"\n')
    if result:
        with open(os.path.join(run_dir, "result.toml"), "w") as f:
            f.write('total_energy = -1.0\n')

    target = link_to
    if target is None:
        target = os.path.join(scratch, unit_id, calc)
        os.makedirs(target, exist_ok=True)
        # One "HDF5"-sized file, standing in for the working files
        #   that make up almost all of a real run's scratch.
        with open(os.path.join(target, "gs_scf-fb.hdf5"), "wb") as f:
            f.write(b"\0" * scratch_bytes)
    os.symlink(target, os.path.join(run_dir, "intermediate"))
    return run_dir


@pytest.fixture
def workspace(tmp_path):
    """A workspace root and its scratch root, as two sibling trees
    (the real ones live on different filesystems)."""

    root = tmp_path / "workspace"
    scratch = tmp_path / "scratch"
    (root / "wingbeats").mkdir(parents=True)
    scratch.mkdir()
    return str(root), str(scratch)


# ==============================================================
#  The default policy: when is a unit done with its scratch?
# ==============================================================

class TestDefaultPolicy:
    """DESIGN 6.2.12: a unit is spent once it finished AND left a
    result.  ``done`` alone is deliberately not enough."""

    def test_done_with_a_result_is_spent(self, workspace):
        root, scratch = workspace
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        spent, reason = tidy_workspace.default_reclaim_policy(run)
        assert spent is True
        assert reason == ""

    @pytest.mark.parametrize("status", ["running", "queued",
                                        "failed", "lost"])
    def test_an_unfinished_unit_is_never_spent(self, workspace,
                                               status):
        # The engine may still be writing, so its working files are
        #   not ours to remove.
        root, scratch = workspace
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2",
                        status=status)
        spent, reason = tidy_workspace.default_reclaim_policy(run)
        assert spent is False
        assert status in reason

    def test_done_without_a_result_is_preserved(self, workspace):
        # The most important half of the rule: a run that finished
        #   but produced nothing is usually exactly the state a
        #   curator wants to investigate, so its working files stay.
        root, scratch = workspace
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2",
                        result=False)
        spent, reason = tidy_workspace.default_reclaim_policy(run)
        assert spent is False
        assert "result.toml" in reason

    def test_a_missing_status_is_preserved(self, workspace):
        root, scratch = workspace
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2",
                        status=None)
        spent, reason = tidy_workspace.default_reclaim_policy(run)
        assert spent is False
        assert "status.toml" in reason


# ==============================================================
#  The refusal that matters most: never follow a link out
# ==============================================================

class TestScratchTargetRefusals:
    """DESIGN 6.2.12 refusal 3.  A cleanup tool a symlink can
    redirect is a hazard, not a convenience."""

    def test_a_link_outside_the_scratch_root_is_refused(
            self, workspace, tmp_path):
        root, scratch = workspace
        elsewhere = tmp_path / "somewhere-else"
        elsewhere.mkdir()
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2",
                        link_to=str(elsewhere))
        target, reason = tidy_workspace.scratch_target(run, scratch)
        assert target is None
        assert "not under scratch root" in reason

    def test_a_sibling_prefix_is_not_mistaken_for_the_root(
            self, workspace, tmp_path):
        # "/scratch-other" must not pass a check for "/scratch".
        #   A plain startswith test would accept it; commonpath
        #   does not.
        root, scratch = workspace
        sibling = str(scratch) + "-other"
        os.makedirs(sibling, exist_ok=True)
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2",
                        link_to=sibling)
        target, reason = tidy_workspace.scratch_target(run, scratch)
        assert target is None
        assert "not under scratch root" in reason

    def test_a_missing_link_is_refused(self, workspace):
        root, scratch = workspace
        run = os.path.join(root, "wingbeats", "si", "kpt-mesh-2-2-2")
        os.makedirs(run)
        target, reason = tidy_workspace.scratch_target(run, scratch)
        assert target is None
        assert "no intermediate link" in reason

    def test_an_already_reclaimed_unit_is_recognised(
            self, workspace):
        # A dangling link is what a previous reclamation LEAVES
        #   behind on purpose, so meeting one again is normal and
        #   must not be an error.
        root, scratch = workspace
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        target, _ = tidy_workspace.scratch_target(run, scratch)
        import shutil
        shutil.rmtree(target)
        again, reason = tidy_workspace.scratch_target(run, scratch)
        assert again is None
        assert "already reclaimed" in reason


# ==============================================================
#  The plan, the filters, and applying it
# ==============================================================

class TestPlanAndApply:

    def test_a_plan_finds_every_run_and_sizes_it(self, workspace):
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2",
                  scratch_bytes=4096)
        _make_run(root, scratch, "si", "kpt-mesh-4-4-4",
                  scratch_bytes=8192)
        plan = tidy_workspace.plan_reclamation(root, scratch)
        assert len(plan) == 2
        assert all(item["ok"] for item in plan)
        assert sum(item["bytes"] for item in plan) == 4096 + 8192

    def test_planning_removes_nothing(self, workspace):
        # Dry run is the default and the plan is only a report:
        #   an operation whose purpose is deletion must not delete
        #   until asked (DESIGN 6.2.12).
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        plan = tidy_workspace.plan_reclamation(root, scratch)
        assert os.path.isdir(plan[0]["target"])

    def test_apply_removes_scratch_and_keeps_the_run_dir(
            self, workspace):
        # Refusals 1 and 2 together: the scratch tree goes, the run
        #   directory and every kept file stay, and the dangling
        #   intermediate link is left in place on purpose so the run
        #   still records where its scratch was.
        root, scratch = workspace
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        plan = tidy_workspace.plan_reclamation(root, scratch)
        target = plan[0]["target"]

        removed, freed, failures = tidy_workspace.apply_reclamation(
            plan)

        assert (removed, failures) == (1, [])
        assert freed > 0
        assert not os.path.exists(target)
        assert os.path.isdir(run)
        assert os.path.isfile(os.path.join(run, "result.toml"))
        assert os.path.isfile(os.path.join(run, "status.toml"))
        assert os.path.islink(os.path.join(run, "intermediate"))

    def test_an_unfinished_unit_survives_apply(self, workspace):
        root, scratch = workspace
        done = _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        busy = _make_run(root, scratch, "si", "kpt-mesh-4-4-4",
                         status="running")
        plan = tidy_workspace.plan_reclamation(root, scratch)
        tidy_workspace.apply_reclamation(plan)

        assert not os.path.exists(
            os.path.realpath(os.path.join(done, "intermediate")))
        assert os.path.isdir(
            os.path.realpath(os.path.join(busy, "intermediate")))

    def test_the_id_filter_selects_one_unit(self, workspace):
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        _make_run(root, scratch, "ge", "kpt-mesh-2-2-2")
        plan = tidy_workspace.plan_reclamation(root, scratch,
                                               ids={"si"})
        assert [item["unit"].split(os.sep)[0] for item in plan] \
            == ["si"]

    def test_the_calc_filter_selects_by_kind(self, workspace):
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        _make_run(root, scratch, "si", "loen-bispectrum-x")
        plan = tidy_workspace.plan_reclamation(
            root, scratch, calc_pattern="loen-*")
        assert len(plan) == 1
        assert plan[0]["unit"].endswith("loen-bispectrum-x")

    def test_the_age_filter_spares_recent_scratch(self, workspace):
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        plan = tidy_workspace.plan_reclamation(root, scratch,
                                               older_than=7.0)
        assert plan[0]["ok"] is False
        assert "days old" in plan[0]["reason"]

    def test_a_client_policy_overrides_the_default(self, workspace):
        # The mechanism/policy split (DESIGN 6.2.12): a client that
        #   needs its working files for longer says so, and the walk
        #   honours it without knowing why.
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        plan = tidy_workspace.plan_reclamation(
            root, scratch,
            policy=lambda run_dir: (False, "client says keep"))
        assert plan[0]["ok"] is False
        assert plan[0]["reason"] == "client says keep"

    def test_the_walk_does_not_descend_into_scratch(
            self, workspace):
        # `intermediate` is a symlink into another tree; walking
        #   through it would find whatever lives there and, worse,
        #   could plan a removal outside the workspace.
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        target = os.path.join(scratch, "si", "kpt-mesh-2-2-2")
        # Plant something that would look like a run directory if
        #   the walk ever followed the link into scratch.
        with open(os.path.join(target, "status.toml"), "w") as f:
            f.write('status = "done"\n')
        found = list(tidy_workspace.find_run_dirs(root))
        assert len(found) == 1
        assert found[0].startswith(os.path.join(root, "wingbeats"))


# ==============================================================
#  CLI surface
# ==============================================================

class TestCli:

    def test_dry_run_is_the_default(self, workspace, capsys):
        root, scratch = workspace
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        rc = tidy_workspace.main([root, "--scratch-root", scratch])
        assert rc == 0
        assert os.path.isdir(
            os.path.realpath(os.path.join(run, "intermediate")))
        out = capsys.readouterr().out
        assert "would free" in out
        assert "--apply" in out

    def test_apply_reports_what_it_freed(self, workspace, capsys):
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        rc = tidy_workspace.main([root, "--scratch-root", scratch,
                                  "--apply"])
        assert rc == 0
        out = capsys.readouterr().out
        assert "reclaimed 1" in out

    def test_a_non_workspace_is_refused(self, tmp_path, capsys):
        rc = tidy_workspace.main([str(tmp_path), "--scratch-root",
                                  str(tmp_path)])
        assert rc == 2
        assert "not a flight workspace" in capsys.readouterr().err

    def test_no_scratch_root_is_refused(self, workspace, capsys,
                                        monkeypatch):
        # Without a scratch root the containment check cannot run,
        #   and that check is the tool's main safety property -- so
        #   refuse rather than proceed without it.
        root, _ = workspace
        monkeypatch.delenv("IMAGO_TEMP", raising=False)
        rc = tidy_workspace.main([root, "--scratch-root", ""])
        assert rc == 2
        assert "no scratch root" in capsys.readouterr().err
