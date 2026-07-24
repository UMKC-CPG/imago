"""Tests for the curator promotion helper (guidance_promote.py;
DESIGN 7.8 curator half / PSEUDOCODE 15.7).

These build a synthetic dataspace under tmp_path -- staging entries
written with the real ``guidance_db.save_entry`` so they round-trip
through ``load_entry`` exactly as the live curator would see them --
and exercise each of the four modes plus the objective acceptance
rule.  No $IMAGO_DATA and no real Imago run are needed: promotion
reads only the staged TOML files.

The acceptance-rule tests that need a schema-*invalid* entry (the
gap/gap_kind inconsistency, a missing verification block) call
``auto_promote_ok`` on an in-memory entry directly, because such an
entry could never be staged through ``save_entry`` +
``load_entry`` in the first place -- the rule is a defense that
stands on its own.
"""

import os

import pytest

import guidance_promote as gp
from guidance_db import (
    CANONICAL_GROUP_ORDER,
    CANONICAL_LATTICE_ORDER,
    Context,
    GuidanceEntry,
    Measured,
    Provenance,
    Signature,
    Verification,
    save_entry,
)


# --------------------------------------------------------------
#  Entry builders
# --------------------------------------------------------------

def _make_entry(*, converged_at=100.0,
                grid_values=(25.0, 50.0, 100.0, 200.0, 400.0),
                grid_energies=(3.0, 2.0, 1.0, 1.0, 1.0),
                gap_ev=2.0, gap_kind="direct", structure="si.skl",
                generated_at="2026-01-01T00:00:00Z",
                metric_threshold=1.0, system_type="crystalline",
                lattice="cubic", verification="default",
                converged_mesh=None, commit="abc123",
                kpoint_integration="gaussian-0.1"):
    """Build an in-memory GuidanceEntry for the rule tests.  The
    defaults describe a clean, auto-promotable crystalline sweep:
    converged at the 0.2 position with a perfectly flat top of
    grid and a self-consistent gap.  ``verification=None`` drops
    the block entirely."""

    composition = tuple(1.0 if group == "group_iv" else 0.0
                        for group in CANONICAL_GROUP_ORDER)
    if system_type == "crystalline":
        family = lattice
        onehot = tuple(1.0 if name == lattice else 0.0
                       for name in CANONICAL_LATTICE_ORDER)
    else:
        family = ""
        onehot = tuple(0.0 for _ in CANONICAL_LATTICE_ORDER)

    if verification == "default":
        verification = Verification(
            grid_values=tuple(grid_values),
            grid_energies=(tuple(grid_energies)
                           if grid_energies is not None else None),
            converged_at=converged_at,
            converged_mesh=(tuple(converged_mesh)
                            if converged_mesh is not None else None),
            metric="total_energy",
            metric_threshold=metric_threshold,
            predictor_confidence=0.9,
            predictor_neighbor_ids=("mp-1",))

    return GuidanceEntry(
        entry_id="", generated_at=generated_at, source="flight",
        signature=Signature(system_type, composition, family, onehot),
        measured=Measured(gap_ev, gap_kind, 0.0, 0.0, converged_at),
        context=Context("fb", "gga-pbe", kpoint_integration, 1.0e-6,
                        8, 100.0),
        verification=verification,
        provenance=Provenance("flight", structure, commit,
                              "guidance_harvest.py"))


def _stage(db_root, **kwargs):
    """Build an entry and write it into staging/ via save_entry,
    returning the staged path.  Distinct ``structure`` values give
    distinct provenance-hash slugs, so several entries can be
    staged without colliding."""

    return save_entry(_make_entry(**kwargs), db_root)


def _count(db_root, area, system_type="crystalline"):
    """Number of *.toml files under db_root/<area>/<system_type>/."""

    import glob
    return len(glob.glob(os.path.join(
        db_root, area, system_type, "*.toml")))


# --------------------------------------------------------------
#  The objective acceptance rule
# --------------------------------------------------------------

def test_auto_promote_ok_accepts_clean_sweep():
    """A mid-grid convergence, a flat top of grid, and a
    consistent gap clear all three conditions."""
    assert gp.auto_promote_ok(_make_entry()) is True


def test_auto_promote_ok_rejects_endpoint_convergence():
    """Converging near the low end of the grid (position < 0.2)
    is suspect -- the grid may have been too narrow."""
    entry = _make_entry(converged_at=50.0)     # (50-25)/375 ~ 0.07
    assert gp.auto_promote_ok(entry) is False


def test_auto_promote_ok_rejects_unflat_top_of_grid():
    """A top-of-grid whose per-atom eV energy SPREAD exceeds
    metric_threshold*10 is not convincingly converged."""
    entry = _make_entry(grid_energies=(3.0, 2.0, 1.0, 3.0, 1.0),
                        metric_threshold=0.01)
    assert gp.auto_promote_ok(entry) is False


def test_auto_promote_ok_rejects_gap_inconsistency():
    """gap_kind must be 'none' iff gap_ev == 0.0 (a defense even
    though the loader also enforces it -- tested in-memory because
    such an entry could never be staged)."""
    bad = _make_entry(gap_ev=2.0, gap_kind="none")
    assert gp.auto_promote_ok(bad) is False


def test_auto_promote_ok_rejects_missing_verification():
    """An entry with no verification block (a manual entry) is
    never auto-promoted -- it has no flatness evidence."""
    assert gp.auto_promote_ok(
        _make_entry(verification=None)) is False


def test_auto_promote_ok_rejects_missing_grid_energies():
    """Verification without grid_energies cannot satisfy the
    flatness test, so it falls to interactive review."""
    entry = _make_entry(grid_energies=None)
    assert gp.auto_promote_ok(entry) is False


# --------------------------------------------------------------
#  move_to_entries
# --------------------------------------------------------------

def test_move_to_entries_moves_and_refuses_collision(tmp_path):
    """A promote is a rename from staging/ to entries/; a
    pre-existing destination is refused, not overwritten."""
    db_root = str(tmp_path / "db")
    staged = _stage(db_root, structure="a.skl")
    moved = gp.move_to_entries(staged, db_root, "crystalline")
    assert os.path.exists(moved)
    assert not os.path.exists(staged)
    assert os.path.dirname(moved).endswith(
        os.path.join("entries", "crystalline"))
    # A second move onto the same destination name must refuse.
    again = _stage(db_root, structure="a.skl")   # same slug
    with pytest.raises(ValueError):
        gp.move_to_entries(again, db_root, "crystalline")


# --------------------------------------------------------------
#  Driver modes
# --------------------------------------------------------------

def test_promote_all_moves_everything(tmp_path):
    """--all promotes every staged entry regardless of the rule
    (including one that fails the auto-promote test)."""
    db_root = str(tmp_path / "db")
    _stage(db_root, structure="good.skl")              # passes rule
    _stage(db_root, structure="edge.skl", converged_at=50.0)  # fails
    results = gp.promote(db_root, "all", output=lambda *a: None)
    assert sorted(a for _, a in results) == ["promoted", "promoted"]
    assert _count(db_root, "staging") == 0
    assert _count(db_root, "entries") == 2


def test_promote_auto_promote_moves_only_passing(tmp_path):
    """--auto-promote moves the rule-passing entry and leaves the
    endpoint-converged one in staging for review."""
    db_root = str(tmp_path / "db")
    _stage(db_root, structure="good.skl")
    _stage(db_root, structure="edge.skl", converged_at=50.0)
    results = dict(gp.promote(db_root, "auto-promote",
                              output=lambda *a: None))
    assert sorted(results.values()) == ["promoted", "skipped"]
    assert _count(db_root, "staging") == 1            # the failing one
    assert _count(db_root, "entries") == 1            # the passing one


def test_promote_dry_run_moves_nothing(tmp_path):
    """--dry-run reports would-promote / would-skip and touches no
    files."""
    db_root = str(tmp_path / "db")
    _stage(db_root, structure="good.skl")
    _stage(db_root, structure="edge.skl", converged_at=50.0)
    lines = []
    results = dict(gp.promote(db_root, "dry-run",
                              output=lines.append))
    assert sorted(results.values()) == ["would-promote", "would-skip"]
    assert _count(db_root, "staging") == 2            # untouched
    assert _count(db_root, "entries") == 0
    assert any("WOULD PROMOTE" in line for line in lines)


def test_promote_interactive_promote_skip_delete(tmp_path):
    """Interactive review honors the curator's PROMOTE / SKIP /
    DELETE choices: promoted files move to entries/, skipped files
    stay in staging/, deleted files are removed entirely."""
    db_root = str(tmp_path / "db")
    _stage(db_root, structure="a.skl")
    _stage(db_root, structure="b.skl")
    _stage(db_root, structure="c.skl")
    answers = iter(["promote", "skip", "delete"])
    results = gp.promote(
        db_root, "interactive",
        ask=lambda prompt: next(answers),
        output=lambda *a: None)
    actions = sorted(action for _, action in results)
    assert actions == ["deleted", "promoted", "skipped"]
    assert _count(db_root, "entries") == 1            # promoted
    assert _count(db_root, "staging") == 1            # skipped (1 deleted)


def test_promote_rejects_unknown_mode(tmp_path):
    """An unrecognized mode aborts rather than silently doing
    nothing."""
    with pytest.raises(ValueError):
        gp.promote(str(tmp_path / "db"), "frobnicate")


def test_ask_choice_normalizes_and_defaults_to_skip():
    """_ask_choice accepts the initials and full words, re-prompts
    on garbage, and treats an empty answer as the safe SKIP."""
    assert gp._ask_choice(lambda p: "p") == "PROMOTE"
    assert gp._ask_choice(lambda p: "DELETE") == "DELETE"
    assert gp._ask_choice(lambda p: "") == "SKIP"
    bad_then_good = iter(["huh?", "s"])
    assert gp._ask_choice(lambda p: next(bad_then_good)) == "SKIP"


# --------------------------------------------------------------
#  Re-run dedup (DESIGN 7.8; PSEUDOCODE 15.7)
# --------------------------------------------------------------
# A re-run of an already-promoted solid is invisible to every
# other guard: the harvest cannot see the promoted corpus, and the
# entry_id slug hashes the timestamp, so a second run of one solid
# never collides.  Promotion is the only stage that can catch it,
# and these tests pin the three outcomes -- redundant, conflict,
# and genuinely new -- plus the batch and dry-run behaviour.

def _promote_one(db_root, **kwargs):
    """Stage an entry and promote it into entries/, returning the
    entry so a test can compare a later re-run against it."""

    path = _stage(db_root, **kwargs)
    gp.move_to_entries(path, db_root, "crystalline")
    return path


def test_dedup_key_ignores_where_the_structure_cache_sits():
    """The key uses the structure's basename, not its path: the
    cache has moved once already (ARCHITECTURE 8.1) and a path
    change must not make a re-run look like a new solid."""
    old = _make_entry(structure="/share/atomicBDB/cache/si.skl")
    new = _make_entry(structure="/share/curation/structures/si.skl")
    assert gp.dedup_key(old) == gp.dedup_key(new)


def test_dedup_key_separates_integration_sub_models():
    """A gaussian and a gaussian-0.1 run of one solid are
    different physics -- the predictor keeps them in separate
    sub-models, so the dedup must not merge them."""
    plain = _make_entry(kpoint_integration="gaussian")
    smeared = _make_entry(kpoint_integration="gaussian-0.1")
    assert gp.dedup_key(plain) != gp.dedup_key(smeared)


def test_rerun_agreeing_on_mesh_is_retired_not_promoted(tmp_path):
    """Same claim, same converged mesh: the promoted entry stands
    untouched and the staged copy moves to superseded/."""
    db_root = str(tmp_path / "db")
    _promote_one(db_root, converged_mesh=(6, 6, 6),
                 generated_at="2026-01-01T00:00:00Z")
    _stage(db_root, converged_mesh=(6, 6, 6),
           generated_at="2026-02-02T00:00:00Z")

    results = gp.promote(db_root, "auto-promote", output=lambda m: None)

    assert [action for _, action in results] == ["superseded"]
    assert _count(db_root, "entries") == 1
    assert _count(db_root, "staging") == 0
    assert _count(db_root, "superseded") == 1


def test_rerun_disagreeing_on_mesh_is_a_conflict(tmp_path):
    """Same claim, different converged mesh -- the shape a code
    change produces.  Nothing moves, and the staged file stays in
    staging for the curator to adjudicate."""
    db_root = str(tmp_path / "db")
    _promote_one(db_root, converged_mesh=(6, 6, 6), commit="old111",
                 generated_at="2026-01-01T00:00:00Z")
    _stage(db_root, converged_mesh=(2, 4, 4), commit="new222",
           generated_at="2026-02-02T00:00:00Z")

    said = []
    results = gp.promote(db_root, "auto-promote", output=said.append)

    assert [action for _, action in results] == ["conflicted"]
    assert _count(db_root, "entries") == 1            # untouched
    assert _count(db_root, "staging") == 1            # left alone
    assert _count(db_root, "superseded") == 0
    # The report has to name both sides well enough to compare
    # them by hand -- both meshes and both commits.
    report = "\n".join(said)
    assert "CONFLICT" in report
    assert "old111" in report and "new222" in report
    assert "[6, 6, 6]" in report and "[2, 4, 4]" in report


def test_a_mesh_that_cannot_be_compared_conflicts(tmp_path):
    """converged_mesh is optional (a manual entry has none).  An
    absent mesh cannot be shown to agree, so the pair goes to the
    curator rather than being waved through."""
    db_root = str(tmp_path / "db")
    _promote_one(db_root, converged_mesh=None,
                 generated_at="2026-01-01T00:00:00Z")
    _stage(db_root, converged_mesh=(6, 6, 6),
           generated_at="2026-02-02T00:00:00Z")

    results = gp.promote(db_root, "auto-promote", output=lambda m: None)

    assert [action for _, action in results] == ["conflicted"]
    assert _count(db_root, "staging") == 1
    assert _count(db_root, "entries") == 1


def test_all_mode_does_not_bypass_the_dedup(tmp_path):
    """--all waives the quality rule, not the correctness guard:
    it must still refuse to store one claim twice."""
    db_root = str(tmp_path / "db")
    _promote_one(db_root, converged_mesh=(6, 6, 6),
                 generated_at="2026-01-01T00:00:00Z")
    _stage(db_root, converged_mesh=(6, 6, 6),
           generated_at="2026-02-02T00:00:00Z")

    results = gp.promote(db_root, "all", output=lambda m: None)

    assert [action for _, action in results] == ["superseded"]
    assert _count(db_root, "entries") == 1
    assert _count(db_root, "superseded") == 1


def test_duplicates_within_one_batch_resolve_to_the_newest(tmp_path):
    """Two staged files can share a claim before either is
    promoted -- the ordinary shape of a solid harvested twice.
    The later run wins and the earlier is retired."""
    db_root = str(tmp_path / "db")
    _stage(db_root, converged_mesh=(6, 6, 6),
           generated_at="2026-01-01T00:00:00Z")
    _stage(db_root, converged_mesh=(6, 6, 6),
           generated_at="2026-03-03T00:00:00Z")

    results = gp.promote(db_root, "all", output=lambda m: None)

    actions = sorted(action for _, action in results)
    assert actions == ["promoted", "superseded"]
    assert _count(db_root, "entries") == 1
    assert _count(db_root, "superseded") == 1
    assert _count(db_root, "staging") == 0


def test_dry_run_reports_the_dedup_but_moves_nothing(tmp_path):
    """A curator must be able to see every outcome -- including
    retirements -- before a single file is touched."""
    db_root = str(tmp_path / "db")
    _promote_one(db_root, converged_mesh=(6, 6, 6),
                 generated_at="2026-01-01T00:00:00Z")
    _stage(db_root, converged_mesh=(6, 6, 6),
           generated_at="2026-02-02T00:00:00Z")

    results = gp.promote(db_root, "dry-run", output=lambda m: None)

    assert [action for _, action in results] == ["would-supersede"]
    assert _count(db_root, "staging") == 1            # nothing moved
    assert _count(db_root, "superseded") == 0
    assert _count(db_root, "entries") == 1


def test_distinct_structures_are_not_duplicates(tmp_path):
    """The six diamond-silicon COD entries are near-identical in
    the feature space but are genuinely different structures.  The
    dedup keys on the structure, so all of them promote."""
    db_root = str(tmp_path / "db")
    for name in ("si_a.skl", "si_b.skl", "si_c.skl"):
        _stage(db_root, structure=name, converged_mesh=(6, 6, 6))

    results = gp.promote(db_root, "all", output=lambda m: None)

    assert [action for _, action in results] == ["promoted"] * 3
    assert _count(db_root, "entries") == 3
    assert _count(db_root, "superseded") == 0
