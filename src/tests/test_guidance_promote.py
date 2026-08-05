## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

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


def test_ask_occupied_choice_offers_replace_and_never_promote():
    """On an OCCUPIED claim the verbs change: REPLACE joins SKIP and
    DELETE, and PROMOTE is deliberately absent -- promoting beside the
    existing record is the one outcome the uniqueness rule exists to
    prevent.  An empty answer still defaults to the safe SKIP."""
    assert gp._ask_occupied_choice(lambda p: "r") == "REPLACE"
    assert gp._ask_occupied_choice(lambda p: "DELETE") == "DELETE"
    assert gp._ask_occupied_choice(lambda p: "") == "SKIP"
    # "p" is not a verb here, so it re-prompts rather than guessing.
    bad_then_good = iter(["p", "s"])
    assert gp._ask_occupied_choice(
        lambda p: next(bad_then_good)) == "SKIP"


# --------------------------------------------------------------
#  Re-run dedup (DESIGN 7.8; PSEUDOCODE 15.7)
# --------------------------------------------------------------
# A re-run of an already-promoted solid is invisible to every
# other guard: the harvest cannot see the promoted corpus, and the
# entry_id slug hashes the timestamp, so a second run of one solid
# never collides.  Promotion is the only stage that can catch it.
#
# The check is an EXISTENCE test, not a comparison: either the
# collection holds a record for the claim or it does not, and that
# has two answers, not three.  The converged meshes are printed and
# never tested, because the action is the same whether they agree
# or not -- promotion has no verb for retracting a reviewed entry
# unasked, so all a comparison could decide is whether the newcomer
# was archived or left in staging to be re-reported by every later
# pass.  These tests pin both answers, the curator's REPLACE, and
# the batch and dry-run behaviour that falls out of applying the
# one rule uniformly.

def _promote_one(db_root, **kwargs):
    """Stage an entry and promote it into entries/, returning the
    entry so a test can compare a later re-run against it."""

    path = _stage(db_root, **kwargs)
    gp.move_to_entries(path, db_root, "crystalline")
    return path


def _promoted_commits(db_root):
    """The build recorded on every promoted entry, read back through
    the promoter's own index so a test sees exactly what a later
    promotion pass would find in ``entries/``."""

    return sorted(entry.provenance.imago_commit for _, entry
                  in gp.load_promoted_entries(db_root).values())


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


def test_a_disagreeing_mesh_takes_the_same_branch(tmp_path):
    """Same claim, DIFFERENT converged mesh -- the shape a code change
    produces.  The mesh is not a branch: the outcome is the ordinary
    occupied one, because promotion has no verb for retracting a
    reviewed entry unasked, so all a comparison could ever decide is
    whether the newcomer was archived or left in staging to be
    re-reported by every later pass (DESIGN 7.8)."""
    db_root = str(tmp_path / "db")
    _promote_one(db_root, converged_mesh=(6, 6, 6), commit="old111",
                 generated_at="2026-01-01T00:00:00Z")
    _stage(db_root, converged_mesh=(2, 4, 4), commit="new222",
           generated_at="2026-02-02T00:00:00Z")

    said = []
    results = gp.promote(db_root, "auto-promote", output=said.append)

    assert [action for _, action in results] == ["superseded"]
    assert _count(db_root, "entries") == 1            # untouched
    assert _count(db_root, "staging") == 0            # archived
    assert _count(db_root, "superseded") == 1
    # The report still has to name both sides well enough to compare
    # them by hand -- both meshes and both builds.  That is where a
    # disagreement gets resolved: by a person, not by a branch.
    report = "\n".join(said)
    assert "OCCUPIED" in report
    assert "old111" in report and "new222" in report
    assert "[6, 6, 6]" in report and "[2, 4, 4]" in report


def test_a_mesh_that_cannot_be_compared_is_not_special(tmp_path):
    """converged_mesh is optional (a manual entry has none).  An
    absent mesh used to count as disagreement; now it is simply
    printed as ``<none recorded>`` and the claim is occupied either
    way."""
    db_root = str(tmp_path / "db")
    _promote_one(db_root, converged_mesh=None,
                 generated_at="2026-01-01T00:00:00Z")
    _stage(db_root, converged_mesh=(6, 6, 6),
           generated_at="2026-02-02T00:00:00Z")

    said = []
    results = gp.promote(db_root, "auto-promote", output=said.append)

    assert [action for _, action in results] == ["superseded"]
    assert _count(db_root, "staging") == 0
    assert _count(db_root, "entries") == 1
    assert "<none recorded>" in "\n".join(said)


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


def test_a_duplicate_within_one_batch_takes_the_ordinary_branch(
        tmp_path):
    """Two staged files can share a claim before either is promoted --
    the ordinary shape of a solid harvested twice.  Nothing special
    happens to them, and that is the point: the promoted index is
    updated as each entry is promoted, so the second file finds the
    claim occupied and takes the same branch any re-run takes.

    There is no separate batch-resolution pass and no ``generated_at``
    tie-break.  Both existed only to save one rule from having to apply
    twice in a row, and one consequence is worth stating because it
    reverses an older contract: promotion no longer judges each staged
    file in isolation, so ``staging/`` IS a uniqueness namespace."""
    db_root = str(tmp_path / "db")
    _stage(db_root, converged_mesh=(6, 6, 6),
           generated_at="2026-01-01T00:00:00Z")
    _stage(db_root, converged_mesh=(6, 6, 6),
           generated_at="2026-03-03T00:00:00Z")

    results = gp.promote(db_root, "all", output=lambda m: None)

    # Which of the two is promoted is decided by the sorted file
    #   order alone -- neither date nor mesh is consulted -- so the
    #   claim under test is that exactly one lands and one is retired.
    actions = sorted(action for _, action in results)
    assert actions == ["promoted", "superseded"]
    assert _count(db_root, "entries") == 1
    assert _count(db_root, "superseded") == 1
    assert _count(db_root, "staging") == 0


def test_a_claim_staged_many_times_still_lands_once(tmp_path):
    """The rule applies once per file after the first, not once.

    The test above stages a pair, which is the shallowest case and
    the only one that was covered.  Depth is not hypothetical: the
    live staging area reached SIX files for one solid, because every
    campaign harvests afresh and the REGENERATE model leaves the
    duplicates for promotion to resolve.  Six is used here for that
    reason rather than picked.

    Nothing in the mechanism should care how deep a group runs --
    the promoted index is live, so the second file finds the claim
    occupied and so does the sixth -- but "should not care" is what
    a test is for.  A batch-resolution pass or a ``generated_at``
    tie-break, both of which this design removed, are exactly the
    shapes that work for two and go wrong for more."""
    db_root = str(tmp_path / "db")
    for month in range(1, 7):
        _stage(db_root, converged_mesh=(6, 6, 6),
               generated_at=f"2026-{month:02d}-01T00:00:00Z")
    assert _count(db_root, "staging") == 6

    results = gp.promote(db_root, "all", output=lambda m: None)

    actions = sorted(action for _, action in results)
    assert actions == ["promoted"] + ["superseded"] * 5
    assert _count(db_root, "entries") == 1
    assert _count(db_root, "superseded") == 5
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


def test_replace_retires_the_promoted_entry_for_the_newcomer(tmp_path):
    """REPLACE is the only route by which anything ever leaves
    ``entries/``, and it exists solely behind a per-record prompt a
    person answers by hand (DESIGN 7.8).  The standing objection --
    that a tool able to retract a reviewed entry can do so by accident
    -- is an argument against an AUTOMATIC retraction and does not
    reach an explicit one.  Without the verb, the report names a
    situation the curator could act on only by moving files
    themselves, which is a worse place to leave them than a prompt."""
    db_root = str(tmp_path / "db")
    _promote_one(db_root, converged_mesh=(6, 6, 6), commit="old111",
                 generated_at="2026-01-01T00:00:00Z")
    _stage(db_root, converged_mesh=(4, 4, 4), commit="new222",
           generated_at="2026-02-02T00:00:00Z")

    results = gp.promote(db_root, "interactive",
                         ask=lambda prompt: "r",
                         output=lambda m: None)

    assert [action for _, action in results] == ["replaced"]
    assert _count(db_root, "entries") == 1
    assert _count(db_root, "superseded") == 1     # the retired one
    assert _count(db_root, "staging") == 0
    # The collection now holds the newcomer, not what it displaced.
    assert _promoted_commits(db_root) == ["new222"]


def test_a_second_replace_against_one_claim_still_works(tmp_path):
    """Why the promoted index carries a PATH beside each entry instead
    of a bare entry.  After the first REPLACE the slot holds the
    newcomer, and retiring THAT one requires the file it now occupies
    -- the destination the move returned, not the path the entry it
    displaced used to sit at.  Store bare entries and the second
    replace has nothing to retire."""
    db_root = str(tmp_path / "db")
    _promote_one(db_root, converged_mesh=(6, 6, 6), commit="first",
                 generated_at="2026-01-01T00:00:00Z")
    _stage(db_root, converged_mesh=(4, 4, 4), commit="second",
           generated_at="2026-02-02T00:00:00Z")
    _stage(db_root, converged_mesh=(2, 2, 2), commit="third",
           generated_at="2026-03-03T00:00:00Z")

    results = gp.promote(db_root, "interactive",
                         ask=lambda prompt: "r",
                         output=lambda m: None)

    assert [action for _, action in results] == ["replaced",
                                                 "replaced"]
    assert _count(db_root, "entries") == 1        # still exactly one
    assert _count(db_root, "superseded") == 2     # both displaced
    assert _count(db_root, "staging") == 0


def test_the_unattended_modes_never_retract(tmp_path):
    """``--all`` and ``--auto-promote`` apply the existence test --
    refusing to store one claim twice is a correctness guard, not a
    quality judgment -- but neither is ever offered REPLACE.  ``--all``
    means "I have reviewed these," not "store them however many times
    they appear," and least of all "retract what a curator already
    accepted."  The prompt is not merely declined here but never
    reached: ``ask`` raises if it is called at all."""
    def never_asked(prompt):
        raise AssertionError(
            "an unattended mode must not prompt: " + prompt)

    for mode in ("all", "auto-promote"):
        db_root = str(tmp_path / mode)
        _promote_one(db_root, converged_mesh=(6, 6, 6),
                     commit="reviewed-and-promoted",
                     generated_at="2026-01-01T00:00:00Z")
        _stage(db_root, converged_mesh=(4, 4, 4), commit="newcomer",
               generated_at="2026-02-02T00:00:00Z")

        results = gp.promote(db_root, mode, ask=never_asked,
                             output=lambda m: None)

        assert [action for _, action in results] == ["superseded"]
        # The reviewed entry stands, exactly as it was accepted.
        assert _promoted_commits(db_root) == ["reviewed-and-promoted"]


def test_dry_run_models_the_index_so_a_second_file_is_not_free(
        tmp_path):
    """A dry run moves nothing, so it has to record where each file
    WOULD land.  Without that, a second staged file claiming the slot
    the first just filled is reported free when a real run would find
    it taken -- and a preview that disagrees with the run it previews
    is worse than no preview at all."""
    db_root = str(tmp_path / "db")
    _stage(db_root, converged_mesh=(6, 6, 6),
           generated_at="2026-01-01T00:00:00Z")
    _stage(db_root, converged_mesh=(6, 6, 6),
           generated_at="2026-03-03T00:00:00Z")

    results = gp.promote(db_root, "dry-run", output=lambda m: None)

    assert sorted(action for _, action in results) == [
        "would-promote", "would-supersede"]
    assert _count(db_root, "staging") == 2        # nothing moved
    assert _count(db_root, "entries") == 0
    assert _count(db_root, "superseded") == 0


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
