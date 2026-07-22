"""test_tidy_scratch.py -- unit tests for the scratch reclamation
tool (``src/scripts/tidy_scratch.py``; DESIGN 6.2.12, PSEUDOCODE
13.8).

The tool deletes, so most of what is worth testing is what it
REFUSES to delete.  Five refusals hold for every root -- only
scratch is removed and never the run directory, an unfinished
unit is never touched, a link pointing outside the scratch area
is never followed, a symlink is never descended while walking,
and a tree holding another run's scratch is deferred rather than
taken -- and a job tree adds two more: a run that has not
declared itself finished is never reclaimed, and a workspace is
never descended into from a job tree.  Each gets a dedicated test
here, alongside the selection filters and the dry-run default.

Two fixtures build miniature roots on tmp_path, one per kind:

* ``workspace`` -- a ``wingbeats/<id>/<calc>/`` run directory
  holding the kept tier (``status.toml``, ``result.toml``);
* ``job_tree`` -- ordinary ``imago.py`` run directories at
  whatever depth, holding a ``runtime`` log and no status at all.

Both plant an ``intermediate`` symlink into a separate "scratch"
tree, mirroring the real two-filesystem layout closely enough to
exercise every path without needing a flight.
"""

import os
from datetime import datetime

import pytest

import tidy_scratch


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


def _make_hand_run(root, scratch, relative, *, complete=True,
                   locked=False, runtime=True, extra_runs=0,
                   scratch_bytes=1024):
    """Create one ordinary ``imago.py`` run directory inside a job
    tree, and return its path.

    A hand run has no ``status.toml`` and no ``result.toml``; it
    proves it finished by leaving no ``imagoLock`` in its scratch
    and ending its ``runtime`` log with the completion marker.
    ``extra_runs`` prepends earlier completed runs, so a test can
    exercise the append behaviour that makes reading the log's
    TAIL correct and searching it wrong.
    """

    run_dir = os.path.join(root, relative)
    os.makedirs(run_dir, exist_ok=True)
    target = os.path.join(scratch, relative)
    os.makedirs(target, exist_ok=True)
    with open(os.path.join(target, "gs_scf-fb.hdf5"), "wb") as f:
        f.write(b"\0" * scratch_bytes)
    if locked:
        with open(os.path.join(target, "imagoLock"), "w") as f:
            f.write("held\n")

    if runtime:
        with open(os.path.join(run_dir, "runtime"), "w") as f:
            for _ in range(extra_runs):
                f.write("Start: earlier run\n")
                f.write("Program Sequence Complete.\n")
            f.write("Start: this run\n")
            if complete:
                f.write("Program Sequence Complete.\n")
            else:
                f.write("sys\t0m0.052s\n")

    os.symlink(target, os.path.join(run_dir, "intermediate"))
    return run_dir


@pytest.fixture
def job_tree(tmp_path):
    """A job-tree root and its scratch root: ordinary run
    directories with no flight above them."""

    root = tmp_path / "jobs"
    scratch = tmp_path / "scratch"
    root.mkdir()
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
        spent, reason = tidy_scratch.default_reclaim_policy(run)
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
        spent, reason = tidy_scratch.default_reclaim_policy(run)
        assert spent is False
        assert status in reason

    def test_done_without_a_result_is_preserved(self, workspace):
        # The most important half of the rule: a run that finished
        #   but produced nothing is usually exactly the state a
        #   curator wants to investigate, so its working files stay.
        root, scratch = workspace
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2",
                        result=False)
        spent, reason = tidy_scratch.default_reclaim_policy(run)
        assert spent is False
        assert "result.toml" in reason

    def test_a_missing_status_is_preserved(self, workspace):
        root, scratch = workspace
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2",
                        status=None)
        spent, reason = tidy_scratch.default_reclaim_policy(run)
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
        target, reason = tidy_scratch.scratch_target(run, scratch)
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
        target, reason = tidy_scratch.scratch_target(run, scratch)
        assert target is None
        assert "not under scratch root" in reason

    def test_a_missing_link_is_refused(self, workspace):
        root, scratch = workspace
        run = os.path.join(root, "wingbeats", "si", "kpt-mesh-2-2-2")
        os.makedirs(run)
        target, reason = tidy_scratch.scratch_target(run, scratch)
        assert target is None
        assert "no intermediate link" in reason

    def test_an_already_reclaimed_unit_is_recognised(
            self, workspace):
        # A dangling link is what a previous reclamation LEAVES
        #   behind on purpose, so meeting one again is normal and
        #   must not be an error.
        root, scratch = workspace
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        target, _ = tidy_scratch.scratch_target(run, scratch)
        import shutil
        shutil.rmtree(target)
        again, reason = tidy_scratch.scratch_target(run, scratch)
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
        plan = tidy_scratch.plan_reclamation(root, scratch)
        assert len(plan) == 2
        assert all(item["ok"] for item in plan)
        assert sum(item["bytes"] for item in plan) == 4096 + 8192

    def test_planning_removes_nothing(self, workspace):
        # Dry run is the default and the plan is only a report:
        #   an operation whose purpose is deletion must not delete
        #   until asked (DESIGN 6.2.12).
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        plan = tidy_scratch.plan_reclamation(root, scratch)
        assert os.path.isdir(plan[0]["target"])

    def test_apply_removes_scratch_and_keeps_the_run_dir(
            self, workspace):
        # Refusals 1 and 2 together: the scratch tree goes, the run
        #   directory and every kept file stay, and the dangling
        #   intermediate link is left in place on purpose so the run
        #   still records where its scratch was.
        root, scratch = workspace
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        plan = tidy_scratch.plan_reclamation(root, scratch)
        target = plan[0]["target"]

        removed, freed, failures = tidy_scratch.apply_reclamation(
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
        plan = tidy_scratch.plan_reclamation(root, scratch)
        tidy_scratch.apply_reclamation(plan)

        assert not os.path.exists(
            os.path.realpath(os.path.join(done, "intermediate")))
        assert os.path.isdir(
            os.path.realpath(os.path.join(busy, "intermediate")))

    def test_the_id_filter_selects_one_unit(self, workspace):
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        _make_run(root, scratch, "ge", "kpt-mesh-2-2-2")
        plan = tidy_scratch.plan_reclamation(root, scratch,
                                               ids={"si"})
        assert [item["unit"].split(os.sep)[0] for item in plan] \
            == ["si"]

    def test_the_calc_filter_selects_by_kind(self, workspace):
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        _make_run(root, scratch, "si", "loen-bispectrum-x")
        plan = tidy_scratch.plan_reclamation(
            root, scratch, calc_pattern="loen-*")
        assert len(plan) == 1
        assert plan[0]["unit"].endswith("loen-bispectrum-x")

    def test_the_age_filter_spares_recent_scratch(self, workspace):
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        plan = tidy_scratch.plan_reclamation(root, scratch,
                                               older_than=7.0)
        assert plan[0]["ok"] is False
        assert "days old" in plan[0]["reason"]

    def test_a_client_policy_overrides_the_default(self, workspace):
        # The mechanism/policy split (DESIGN 6.2.12): a client that
        #   needs its working files for longer says so, and the walk
        #   honours it without knowing why.
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        plan = tidy_scratch.plan_reclamation(
            root, scratch,
            policy=lambda run_dir, target: (False,
                                            "client says keep"))
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
        found = list(tidy_scratch.find_run_dirs(root))
        assert len(found) == 1
        assert found[0].startswith(os.path.join(root, "wingbeats"))


# ==============================================================
#  CLI surface
# ==============================================================

class TestCli:

    def test_dry_run_is_the_default(self, workspace, capsys):
        root, scratch = workspace
        run = _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        rc = tidy_scratch.main([root, "--scratch-root", scratch])
        assert rc == 0
        assert os.path.isdir(
            os.path.realpath(os.path.join(run, "intermediate")))
        out = capsys.readouterr().out
        assert "would free" in out
        assert "--apply" in out

    def test_apply_reports_what_it_freed(self, workspace, capsys):
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        rc = tidy_scratch.main([root, "--scratch-root", scratch,
                                  "--apply"])
        assert rc == 0
        out = capsys.readouterr().out
        assert "reclaimed 1" in out

    def test_a_root_of_neither_kind_is_refused(self, tmp_path,
                                               capsys):
        # Neither a wingbeats/ workspace nor a job tree holding
        #   any run directory, so there is nothing to reclaim and
        #   no contract under which to try.
        rc = tidy_scratch.main([str(tmp_path), "--scratch-root",
                                  str(tmp_path)])
        assert rc == 2
        assert "holds neither" in capsys.readouterr().err

    def test_no_scratch_root_is_refused(self, workspace, capsys,
                                        monkeypatch):
        # Without a scratch root the containment check cannot run,
        #   and that check is the tool's main safety property -- so
        #   refuse rather than proceed without it.
        root, _ = workspace
        monkeypatch.delenv("IMAGO_TEMP", raising=False)
        rc = tidy_scratch.main([root, "--scratch-root", ""])
        assert rc == 2
        assert "no scratch root" in capsys.readouterr().err


# ==============================================================
#  Root detection: one call handles one kind
# ==============================================================

class TestRootDetection:
    """DESIGN 6.2.12: the tool decides what kind of thing it was
    pointed at by looking, and that decision fixes the contract
    for the whole call."""

    def test_a_wingbeats_root_is_a_workspace(self, workspace):
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        assert tidy_scratch.detect_root_kind(root) == \
            tidy_scratch.WORKSPACE_ROOT

    def test_a_tree_of_run_dirs_is_a_job_tree(self, job_tree):
        root, scratch = job_tree
        _make_hand_run(root, scratch, "c/diamond/full")
        assert tidy_scratch.detect_root_kind(root) == \
            tidy_scratch.JOB_TREE_ROOT

    def test_a_root_with_neither_is_no_target(self, tmp_path):
        assert tidy_scratch.detect_root_kind(str(tmp_path)) is None

    def test_an_empty_workspace_is_still_a_workspace(self,
                                                     workspace):
        # wingbeats/ decides it, not the presence of units: an
        #   emptied workspace must not be re-read as a job tree
        #   and judged by the other contract.
        root, _ = workspace
        assert tidy_scratch.detect_root_kind(root) == \
            tidy_scratch.WORKSPACE_ROOT


# ==============================================================
#  The job-tree policy: refusal 5
# ==============================================================

class TestHandRunPolicy:
    """DESIGN 6.2.12 refusal 5: a hand run is reclaimable only
    when its scratch holds no ``imagoLock`` AND its ``runtime``
    log ends with the completion marker.  Both, not either."""

    def _policy(self, root, scratch, relative, **kwargs):
        run = _make_hand_run(root, scratch, relative, **kwargs)
        return tidy_scratch.hand_run_policy(
            run, os.path.join(scratch, relative))

    def test_a_finished_run_is_spent(self, job_tree):
        root, scratch = job_tree
        spent, reason = self._policy(root, scratch, "al")
        assert spent is True
        assert reason == ""

    def test_a_held_lock_refuses(self, job_tree):
        # The lock is taken before any work begins, so its
        #   presence means the run owns the directory now or died
        #   without releasing it.  Either way it is not ours.
        root, scratch = job_tree
        spent, reason = self._policy(root, scratch, "al",
                                     locked=True)
        assert spent is False
        assert "imagoLock" in reason

    def test_a_lock_refuses_even_when_the_marker_is_present(
            self, job_tree):
        # BOTH conditions are required.  A stale lock beside a
        #   completed log is exactly the interrupted-run state
        #   that must be preserved for diagnosis.
        root, scratch = job_tree
        spent, _ = self._policy(root, scratch, "al", locked=True,
                                complete=True)
        assert spent is False

    def test_a_log_ending_mid_run_refuses(self, job_tree):
        root, scratch = job_tree
        spent, reason = self._policy(root, scratch, "al",
                                     complete=False)
        assert spent is False
        assert "mid-run" in reason

    def test_the_tail_decides_not_the_presence_of_a_marker(
            self, job_tree):
        # The runtime log is opened in APPEND mode, so a directory
        #   run four times holds four markers.  A fifth run that
        #   was interrupted leaves all four in place -- searching
        #   the file would call this complete, and delete the
        #   scratch of a killed run.  Only the tail is truthful.
        root, scratch = job_tree
        run = _make_hand_run(root, scratch, "c/diamond/full2",
                             complete=False, extra_runs=4)
        with open(os.path.join(run, "runtime")) as handle:
            assert handle.read().count(
                "Program Sequence Complete.") == 4
        spent, reason = tidy_scratch.hand_run_policy(
            run, os.path.join(scratch, "c/diamond/full2"))
        assert spent is False
        assert "mid-run" in reason

    def test_no_runtime_log_refuses(self, job_tree):
        # Absence of evidence is refusal, not permission.
        root, scratch = job_tree
        spent, reason = self._policy(root, scratch, "al",
                                     runtime=False)
        assert spent is False
        assert "no runtime log" in reason


# ==============================================================
#  The job-tree walk: refusal 6
# ==============================================================

class TestJobTreeWalk:
    """DESIGN 6.2.12 refusal 6: a job-tree walk never descends
    into a workspace, and names each one it declined so the bytes
    left behind stay visible."""

    def test_run_dirs_are_found_at_any_depth(self, job_tree):
        # Job trees are organized however their author liked, so
        #   the walk keys on the `intermediate` link rather than
        #   on a fixed depth.
        root, scratch = job_tree
        _make_hand_run(root, scratch, "al")
        _make_hand_run(root, scratch, "c/diamond/full")
        found = sorted(path for kind, path
                       in tidy_scratch.find_job_run_dirs(root)
                       if kind == "run")
        assert [os.path.relpath(p, root) for p in found] == \
            ["al", "c/diamond/full"]

    def test_a_nested_workspace_is_declined_not_entered(
            self, job_tree, tmp_path):
        root, scratch = job_tree
        _make_hand_run(root, scratch, "al")
        nested = os.path.join(root, "campaign", "workspace")
        os.makedirs(nested)
        _make_run(nested, str(tmp_path / "wsscratch"), "si",
                  "kpt-mesh-2-2-2")

        found = list(tidy_scratch.find_job_run_dirs(root))
        workspaces = [p for kind, p in found if kind == "workspace"]
        runs = [p for kind, p in found if kind == "run"]

        assert [os.path.relpath(p, root) for p in workspaces] == \
            ["campaign/workspace"]
        # The workspace's own unit is inside it, so declining to
        #   descend means it never appears as a job-tree run.
        assert [os.path.relpath(p, root) for p in runs] == ["al"]

    def test_a_declined_workspace_is_reported(self, job_tree,
                                              tmp_path, capsys):
        root, scratch = job_tree
        _make_hand_run(root, scratch, "al")
        nested = os.path.join(root, "campaign", "workspace")
        os.makedirs(nested)
        _make_run(nested, str(tmp_path / "wsscratch"), "si",
                  "kpt-mesh-2-2-2")

        plan = tidy_scratch.plan_reclamation(root, scratch)
        tidy_scratch.print_report(plan, applied=False)
        out = capsys.readouterr().out
        assert "campaign/workspace" in out
        assert "not descended" in out


# ==============================================================
#  Ageing scratch by its newest file
# ==============================================================

class TestAgeing:
    """DESIGN 6.2.12: age is the newest mtime anywhere in the
    tree, never the scratch directory's own."""

    def test_a_recent_file_keeps_an_old_directory_young(
            self, job_tree):
        # A directory's mtime moves only when entries are added or
        #   removed, so a job that has been writing into an
        #   already-created HDF5 for a week still presents a
        #   week-old directory.  Ageing by the directory would
        #   call this stale and delete a live run's scratch.
        root, scratch = job_tree
        _make_hand_run(root, scratch, "al")
        target = os.path.join(scratch, "al")
        old = datetime.now().timestamp() - 30 * 86400.0
        os.utime(target, (old, old))          # directory: ancient
        # The file inside it, though, was written moments ago.

        plan = tidy_scratch.plan_reclamation(root, scratch,
                                             older_than=7)
        assert plan[0]["ok"] is False
        assert "days old" in plan[0]["reason"]

    def test_a_genuinely_stale_tree_passes_the_age_filter(
            self, job_tree):
        root, scratch = job_tree
        _make_hand_run(root, scratch, "al")
        target = os.path.join(scratch, "al")
        old = datetime.now().timestamp() - 30 * 86400.0
        for name in os.listdir(target):
            os.utime(os.path.join(target, name), (old, old))
        os.utime(target, (old, old))

        plan = tidy_scratch.plan_reclamation(root, scratch,
                                             older_than=7)
        assert plan[0]["ok"] is True


# ==============================================================
#  The job tree end to end
# ==============================================================

class TestJobTreeCli:
    """The whole contract through the command line."""

    def test_a_job_tree_reclaims_only_finished_runs(self,
                                                    job_tree):
        root, scratch = job_tree
        _make_hand_run(root, scratch, "al")
        _make_hand_run(root, scratch, "knbo3/cubic",
                       complete=False)
        _make_hand_run(root, scratch, "c/diamond/full2",
                       locked=True)

        plan = tidy_scratch.plan_reclamation(root, scratch)
        reclaimable = sorted(item["unit"] for item in plan
                             if item["ok"])
        assert reclaimable == ["al"]

        tidy_scratch.apply_reclamation(plan)
        assert not os.path.isdir(os.path.join(scratch, "al"))
        # The two refused runs keep every working file.
        assert os.path.isfile(os.path.join(
            scratch, "knbo3/cubic", "gs_scf-fb.hdf5"))
        assert os.path.isfile(os.path.join(
            scratch, "c/diamond/full2", "gs_scf-fb.hdf5"))

    def test_the_run_directory_itself_is_never_touched(
            self, job_tree):
        root, scratch = job_tree
        run = _make_hand_run(root, scratch, "al")
        plan = tidy_scratch.plan_reclamation(root, scratch)
        tidy_scratch.apply_reclamation(plan)
        assert os.path.isdir(run)
        assert os.path.isfile(os.path.join(run, "runtime"))
        # The link is left dangling on purpose, so the run still
        #   records where its scratch was.
        assert os.path.islink(os.path.join(run, "intermediate"))

    def test_the_match_filter_selects_by_path(self, job_tree):
        root, scratch = job_tree
        _make_hand_run(root, scratch, "al")
        _make_hand_run(root, scratch, "c/diamond/full")
        _make_hand_run(root, scratch, "c/graphite")

        plan = tidy_scratch.plan_reclamation(root, scratch,
                                             match="c/*")
        assert sorted(item["unit"] for item in plan) == \
            ["c/diamond/full", "c/graphite"]

    def test_a_workspace_filter_on_a_job_tree_is_refused(
            self, job_tree, capsys):
        # A filter that silently did nothing would widen a
        #   deleting run to everything it was meant to narrow.
        root, scratch = job_tree
        _make_hand_run(root, scratch, "al")
        rc = tidy_scratch.main([root, "--scratch-root", scratch,
                                "--id", "si"])
        assert rc == 2
        assert "--id" in capsys.readouterr().err

    def test_a_job_tree_filter_on_a_workspace_is_refused(
            self, workspace, capsys):
        root, scratch = workspace
        _make_run(root, scratch, "si", "kpt-mesh-2-2-2")
        rc = tidy_scratch.main([root, "--scratch-root", scratch,
                                "--match", "si/*"])
        assert rc == 2
        assert "--match" in capsys.readouterr().err

    def test_the_dry_run_default_removes_nothing(self, job_tree):
        root, scratch = job_tree
        _make_hand_run(root, scratch, "al")
        rc = tidy_scratch.main([root, "--scratch-root", scratch])
        assert rc == 0
        assert os.path.isfile(os.path.join(scratch, "al",
                                           "gs_scf-fb.hdf5"))


# ==============================================================
#  Nested scratch: refusal 5
# ==============================================================

class TestNestedScratch:
    """DESIGN 6.2.12 refusal 5.  Scratch mirrors the run
    directory's path, so a run nested inside another has its
    scratch nested too, and removing the outer tree would take
    the inner one with it."""

    def test_an_outer_tree_holding_another_run_is_deferred(
            self, job_tree):
        root, scratch = job_tree
        _make_hand_run(root, scratch, "knbo3/cubic")
        _make_hand_run(root, scratch, "knbo3/cubic/debug")
        # The layout the refusal is about: one scratch inside the
        #   other, because the run directories nest.
        assert os.path.isdir(os.path.join(scratch, "knbo3/cubic",
                                          "debug"))

        plan = tidy_scratch.plan_reclamation(root, scratch)
        by_unit = {item["unit"]: item for item in plan}
        assert by_unit["knbo3/cubic"]["ok"] is False
        assert "nested" in by_unit["knbo3/cubic"]["reason"]
        # The inner run is perfectly reclaimable on its own.
        assert by_unit["knbo3/cubic/debug"]["ok"] is True

    def test_a_second_pass_reclaims_the_outer_tree(self,
                                                   job_tree):
        # Nothing is lost by deferring: once the inner scratch is
        #   gone the outer tree is an ordinary candidate again.
        root, scratch = job_tree
        _make_hand_run(root, scratch, "knbo3/cubic")
        _make_hand_run(root, scratch, "knbo3/cubic/debug")

        first = tidy_scratch.plan_reclamation(root, scratch)
        tidy_scratch.apply_reclamation(first)
        second = tidy_scratch.plan_reclamation(root, scratch)
        by_unit = {item["unit"]: item for item in second}
        assert by_unit["knbo3/cubic"]["ok"] is True

        tidy_scratch.apply_reclamation(second)
        assert not os.path.isdir(os.path.join(scratch,
                                              "knbo3/cubic"))

    def test_a_running_inner_run_is_not_taken_as_collateral(
            self, job_tree):
        # The reason the refusal is not merely about accounting.
        #   The inner run holds its lock, so refusal 6 declines
        #   it -- and without refusal 5 the outer removal would
        #   have deleted its working files anyway.
        root, scratch = job_tree
        _make_hand_run(root, scratch, "knbo3/cubic")
        _make_hand_run(root, scratch, "knbo3/cubic/debug",
                       locked=True)

        plan = tidy_scratch.plan_reclamation(root, scratch)
        tidy_scratch.apply_reclamation(plan)
        assert os.path.isfile(os.path.join(
            scratch, "knbo3/cubic/debug", "gs_scf-fb.hdf5"))

    def test_containment_is_checked_against_filtered_out_runs(
            self, job_tree):
        # The filter narrows what is REMOVED, never what is
        #   protected: an inner run excluded by --match must
        #   still stop the outer removal that would delete it.
        root, scratch = job_tree
        _make_hand_run(root, scratch, "knbo3/cubic")
        _make_hand_run(root, scratch, "knbo3/cubic/debug")

        plan = tidy_scratch.plan_reclamation(
            root, scratch, match="knbo3/cubic")
        assert [item["unit"] for item in plan] == ["knbo3/cubic"]
        assert plan[0]["ok"] is False
        assert "nested" in plan[0]["reason"]

    def test_a_sibling_tree_is_not_mistaken_for_a_nested_one(
            self, job_tree):
        # "cubic2" starts with "cubic" but is not inside it; a
        #   prefix test rather than a path test would defer both.
        root, scratch = job_tree
        _make_hand_run(root, scratch, "knbo3/cubic")
        _make_hand_run(root, scratch, "knbo3/cubic2")

        plan = tidy_scratch.plan_reclamation(root, scratch)
        assert all(item["ok"] for item in plan)

    def test_the_refusal_is_universal_and_fires_in_a_workspace(
            self, workspace):
        # Refusal 5 holds for every root, not only job trees.  A
        #   workspace CAN nest run directories: DESIGN 6.2.4 makes
        #   the <calc> level optional, so a unit may sit directly
        #   under its id AND carry a sub-run below that.  Scratch
        #   mirrors the path either way, so the outer must defer.
        root, scratch = workspace
        outer = os.path.join(root, "wingbeats", "si")
        inner = os.path.join(outer, "kpt-mesh-2-2-2")
        os.makedirs(inner)
        for run, rel in ((outer, "si"),
                         (inner, "si/kpt-mesh-2-2-2")):
            with open(os.path.join(run, "status.toml"), "w") as f:
                f.write('status = "done"\n')
            with open(os.path.join(run, "result.toml"), "w") as f:
                f.write("total_energy = -1.0\n")
            target = os.path.join(scratch, rel)
            os.makedirs(target)
            with open(os.path.join(target, "gs_scf-fb.hdf5"),
                      "wb") as f:
                f.write(b"\0" * 1024)
            os.symlink(target, os.path.join(run, "intermediate"))

        plan = tidy_scratch.plan_reclamation(root, scratch)
        by_unit = {item["unit"]: item for item in plan}
        assert by_unit["si"]["ok"] is False
        assert "nested" in by_unit["si"]["reason"]
        assert by_unit["si/kpt-mesh-2-2-2"]["ok"] is True
