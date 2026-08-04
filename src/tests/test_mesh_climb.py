## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

"""Tests for the producer-side k-point mesh primitives
(``mesh_climb``; DESIGN 3.12 / 2.7; PSEUDOCODE 4c / 4e.1).

These exercise the pure count-vector arithmetic the historical-
guidance producer uses to reason about meshes before it dispatches
anything: the axis classes that couple reciprocal axes, the
isotropic count selection that meets a density floor, and the
one-step rung moves the climb takes.  Everything here is numpy-
free pure computation, so no ``$IMAGO_DATA`` and no
``StructureControl`` are needed -- the module never touches disk.

The load-bearing behavioural checks reproduce the worked mesh
sequences of DESIGN 3.12.2 by hand:

  - a cubic cell climbs in lockstep, [2,2,2] -> [3,3,3] -> ...,
    because all three axes share one class;
  - a hexagonal cell steps its coupled in-plane pair several
    times per single c-axis bump.

The orthorhombic example of 3.12.2 is deliberately NOT pinned to a
literal sequence here: which class the rung rule picks at each
step depends on the specific reciprocal magnitudes of the real
cell, so a generic set of magnitudes need not reproduce the
document's illustrative numbers.  It is checked instead by the
structural invariant the rule guarantees for any cell -- each rung
moves exactly one class and never lowers a count.
"""

import math
from collections import namedtuple

import pytest

import mesh_climb

# A minimal stand-in for the producer's ``Rung(mesh, energy, gap)``:
#   the bracket-refine mesh helpers read only ``.mesh``, so the tests
#   carry a nominal energy the geometry never touches (and no gap).
_Rung = namedtuple("_Rung", ["mesh", "energy"])


def _rung(mesh):
    """A rung at ``mesh`` with a placeholder energy (unused by the
    pure mesh helpers under test here)."""
    return _Rung(mesh, 0.0)


# ------------------------------------------------------------------
#  Axis-class labels used throughout
#
#  The labels are union-find roots, so their absolute values are
#  arbitrary; only the equality pattern matters.  A cubic cell
#  couples all three axes into one class; a hexagonal cell couples
#  the two in-plane axes and leaves c on its own; an orthorhombic
#  cell couples nothing.
# ------------------------------------------------------------------
CUBIC_CLASSES = [0, 0, 0]
HEX_CLASSES = [0, 0, 1]
ORTHO_CLASSES = [0, 1, 2]

# Graphite-derived reciprocal magnitudes (a = 2.46, c = 6.70 Ang),
#   the values hand-traced against DESIGN 3.12.2.  The in-plane
#   pair share a magnitude; the c* axis is much shorter.
HEX_RECIP_MAG = [2.949, 2.949, 0.938]

# Equal magnitudes stand in for a cubic cell; the actual number is
#   irrelevant because every axis shares it.
CUBIC_RECIP_MAG = [1.0, 1.0, 1.0]

# Three distinct magnitudes stand in for an orthorhombic cell.
ORTHO_RECIP_MAG = [1.0, 1.8, 2.6]


# ==================================================================
#  Small 3x3 linear algebra
# ==================================================================

def test_inverse_of_identity_is_identity():
    """Inverting the identity returns the identity.  (Compared as a
    flat sequence because pytest.approx does not nest.)"""
    identity = mesh_climb._identity3x3()
    inverse = mesh_climb._inverse3x3(identity)
    flat = [inverse[row][column]
            for row in range(3) for column in range(3)]
    assert flat == pytest.approx(
        [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])


def test_inverse_times_original_is_identity():
    """A * A^-1 is the identity for a non-trivial invertible
    matrix, confirming the adjugate/determinant formula."""
    matrix = [[2.0, 1.0, 0.0],
              [1.0, 3.0, 1.0],
              [0.0, 1.0, 2.0]]
    product = mesh_climb._mat_mul(
        matrix, mesh_climb._inverse3x3(matrix))
    for row in range(3):
        for column in range(3):
            expected = 1.0 if row == column else 0.0
            assert product[row][column] == pytest.approx(expected)


def test_inverse_of_singular_matrix_raises():
    """A singular matrix (a repeated row) has no inverse, so the
    routine raises rather than dividing by a zero determinant."""
    singular = [[1.0, 2.0, 3.0],
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0]]
    with pytest.raises(ValueError):
        mesh_climb._inverse3x3(singular)


def test_determinant_matches_known_value():
    """The determinant of a small integer matrix matches the
    value computed by hand along the first row."""
    matrix = [[1, 2, 3], [0, 1, 4], [5, 6, 0]]
    assert mesh_climb._determinant3x3(matrix) == pytest.approx(1.0)


def test_round_to_integer_cleans_roundoff():
    """Entries a hair away from integers are snapped to those
    integers when they are all within tolerance."""
    almost = [[1.0000001, -0.0000002, 0.0],
              [0.0, 1.0, 0.9999999],
              [-1.0000001, 0.0, 1.0]]
    cleaned = mesh_climb._round_to_integer_matrix(almost)
    assert cleaned == [[1, 0, 0], [0, 1, 1], [-1, 0, 1]]


def test_round_to_integer_rejects_non_integer():
    """A genuinely fractional entry means the inputs are not a
    lattice automorphism, so the routine raises."""
    fractional = [[1.0, 0.5, 0.0],
                  [0.0, 1.0, 0.0],
                  [0.0, 0.0, 1.0]]
    with pytest.raises(ValueError):
        mesh_climb._round_to_integer_matrix(fractional)


# ==================================================================
#  Axis classes (PSEUDOCODE 4c.1)
# ==================================================================

def test_identity_couples_nothing():
    """With only the identity operation, no off-diagonal coupling
    exists, so all three axes land in distinct classes."""
    classes = mesh_climb.compute_axis_classes(
        [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]])
    assert len(mesh_climb.distinct_classes(classes)) == 3


def test_fourfold_z_couples_the_in_plane_pair():
    """A 90-degree rotation about z maps x <-> y off the diagonal,
    coupling axes 0 and 1 while c stays free -- the hexagonal /
    tetragonal class pattern."""
    rot_z = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    classes = mesh_climb.compute_axis_classes([rot_z])
    assert classes[0] == classes[1]
    assert classes[2] != classes[0]


def test_two_operations_chain_all_axes_together():
    """One operation coupling {0,1} and another coupling {1,2}
    merge transitively, so all three axes become one class -- the
    cubic pattern."""
    couple_xy = [[0, 1, 0], [1, 0, 0], [0, 0, 1]]
    couple_yz = [[1, 0, 0], [0, 0, 1], [0, 1, 0]]
    classes = mesh_climb.compute_axis_classes(
        [couple_xy, couple_yz])
    assert classes[0] == classes[1] == classes[2]


def test_full_mode_collapses_to_recip_of_ops():
    """In "full" cell mode the change of basis is the identity, so
    the resolved classes must equal those of the reciprocal twins
    ``transpose(inverse(op))`` of the given operations -- exactly
    what the code computes with no conjugation."""
    identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    rot_z = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    ops = [identity, rot_z]
    lattice = mesh_climb._identity3x3()

    resolved = mesh_climb.axis_classes_for_cell(
        ops, lattice, lattice, "full")
    expected = mesh_climb.compute_axis_classes(
        [mesh_climb._transpose(mesh_climb._inverse3x3(op))
         for op in ops])
    assert resolved == expected
    # The fourfold-z twin still couples the in-plane pair.
    assert resolved[0] == resolved[1] != resolved[2]


def test_prim_mode_with_equal_cells_matches_full():
    """When the loaded cell equals the conventional cell the change
    of basis T = Lc^-1 * L is the identity, so "prim" conjugation
    must reproduce the "full" result -- a check that the full
    conjugation path collapses correctly in the degenerate case."""
    rot_z = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    ops = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]], rot_z]
    lattice = [[2.0, 0.0, 0.0],
               [0.0, 2.0, 0.0],
               [0.0, 0.0, 3.0]]

    via_prim = mesh_climb.axis_classes_for_cell(
        ops, lattice, lattice, "prim")
    via_full = mesh_climb.axis_classes_for_cell(
        ops, lattice, lattice, "full")
    assert via_prim == via_full


def test_axis_classes_rejects_unknown_cell_mode():
    """Only "full" and "prim" are meaningful cell modes; anything
    else is a programming error and raises."""
    identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    lattice = mesh_climb._identity3x3()
    with pytest.raises(ValueError):
        mesh_climb.axis_classes_for_cell(
            [identity], lattice, lattice, "sideways")


# ==================================================================
#  distinct_classes / spacing_spread helpers
# ==================================================================

def test_distinct_classes_preserves_first_seen_order():
    """Distinct labels come back in the order they first appear, so
    class iteration is deterministic."""
    assert mesh_climb.distinct_classes([2, 2, 5]) == [2, 5]
    assert mesh_climb.distinct_classes([7, 3, 7]) == [7, 3]


def test_spacing_spread_is_max_minus_min():
    """The spread is the largest minus the smallest of the three
    per-axis spacings ``|b_i| / n_i``."""
    # spacings = [2/2, 4/1, 6/3] = [1, 4, 2]; spread = 4 - 1 = 3.
    spread = mesh_climb.spacing_spread([2, 1, 3], [2.0, 4.0, 6.0])
    assert spread == pytest.approx(3.0)


def test_spacing_spread_zero_for_isotropic_mesh():
    """Equal spacings on every axis give a spread of zero, the
    perfectly isotropic mesh the selection drives toward."""
    spread = mesh_climb.spacing_spread([2, 2, 2], [1.0, 1.0, 1.0])
    assert spread == pytest.approx(0.0)


# ==================================================================
#  Axial count selection (PSEUDOCODE 4c.2)
# ==================================================================

def test_gamma_sentinel_returns_single_point():
    """A non-positive density is the Gamma sentinel (DESIGN 3.6)
    and resolves to a single k-point regardless of the cell."""
    assert mesh_climb.select_axial_counts(
        0.0, HEX_RECIP_MAG, 1.0, HEX_CLASSES) == [1, 1, 1]
    assert mesh_climb.select_axial_counts(
        -5.0, CUBIC_RECIP_MAG, 1.0, CUBIC_CLASSES) == [1, 1, 1]


def test_cubic_selection_is_isotropic():
    """For a cubic cell with equal reciprocal magnitudes the three
    counts must come out equal -- a single shared class cannot be
    split -- and a density of 27 per unit volume lands squarely on
    the 3x3x3 mesh."""
    counts = mesh_climb.select_axial_counts(
        27.0, CUBIC_RECIP_MAG, 1.0, CUBIC_CLASSES)
    assert counts == [3, 3, 3]


def test_selection_meets_the_density_floor():
    """Whatever the cell, the chosen full-mesh point count must
    reach the requested floor ``density * recip_cell_volume``."""
    density = 200.0
    recip_cell_volume = 1.0
    counts = mesh_climb.select_axial_counts(
        density, HEX_RECIP_MAG, recip_cell_volume, HEX_CLASSES)
    floor = density * recip_cell_volume
    assert counts[0] * counts[1] * counts[2] >= floor


def test_selection_keeps_coupled_axes_equal():
    """A hexagonal cell's two in-plane axes share a class, so the
    selection must never give them different counts, no matter the
    density."""
    for density in (10.0, 55.0, 130.0, 400.0):
        counts = mesh_climb.select_axial_counts(
            density, HEX_RECIP_MAG, 1.0, HEX_CLASSES)
        assert counts[0] == counts[1]


# ==================================================================
#  bump -- moving a whole class (PSEUDOCODE 4e.1)
# ==================================================================

def test_bump_moves_only_the_named_class():
    """Bumping a class adds the step to every axis in that class
    and leaves the others untouched."""
    assert mesh_climb.bump([2, 2, 1], HEX_CLASSES, 0, +1) \
        == [3, 3, 1]
    assert mesh_climb.bump([2, 2, 1], HEX_CLASSES, 1, +1) \
        == [2, 2, 2]


def test_bump_does_not_mutate_its_input():
    """``bump`` returns a fresh list and leaves the caller's counts
    unchanged, so callers can trial many bumps safely."""
    counts = [2, 2, 1]
    mesh_climb.bump(counts, HEX_CLASSES, 0, +1)
    assert counts == [2, 2, 1]


# ==================================================================
#  Rung mechanics -- the worked sequences of DESIGN 3.12.2
# ==================================================================

def test_cubic_climbs_in_lockstep():
    """All three axes share one class, so each rung raises every
    axis together: [2,2,2] -> [3,3,3] -> [4,4,4]."""
    step_one = mesh_climb.climb_one_rung(
        [2, 2, 2], CUBIC_CLASSES, CUBIC_RECIP_MAG)
    assert step_one == [3, 3, 3]
    step_two = mesh_climb.climb_one_rung(
        step_one, CUBIC_CLASSES, CUBIC_RECIP_MAG)
    assert step_two == [4, 4, 4]


def test_hexagonal_steps_in_plane_pair_then_c_axis():
    """The graphite-like cell reproduces DESIGN 3.12.2 exactly: the
    coupled in-plane pair climbs three times before the short c*
    axis is bumped, then the pair resumes."""
    expected = [[2, 2, 1], [3, 3, 1], [4, 4, 1], [4, 4, 2],
                [5, 5, 2]]
    counts = expected[0]
    for next_counts in expected[1:]:
        counts = mesh_climb.climb_one_rung(
            counts, HEX_CLASSES, HEX_RECIP_MAG)
        assert counts == next_counts


def test_orthorhombic_moves_one_class_and_never_descends():
    """With three independent axes the rung rule's choice depends
    on the real cell, so rather than pin a literal sequence we
    assert the invariant it guarantees for any cell: each rung
    raises exactly one axis by one and lowers none."""
    counts = [2, 2, 2]
    for _ in range(8):
        next_counts = mesh_climb.climb_one_rung(
            counts, ORTHO_CLASSES, ORTHO_RECIP_MAG)
        differences = [next_counts[axis] - counts[axis]
                       for axis in range(3)]
        assert sorted(differences) == [0, 0, 1]
        counts = next_counts


def test_descend_produces_a_lower_neighbour():
    """Descending a rung yields a mesh that climbs straight back to
    the original: ``climb(descend(m)) == m``.  It need not be the
    same predecessor a forward climb arrived from -- a mesh can have
    more than one lower neighbour (e.g. hexagonal [4,4,2] descends
    from both [4,4,1] and [3,3,2]) -- but it is always a genuine
    one, which is what seeding the climb below the prediction needs
    (4e.4)."""
    cases = [
        (CUBIC_CLASSES, CUBIC_RECIP_MAG, [3, 3, 3]),
        (HEX_CLASSES, HEX_RECIP_MAG, [4, 4, 2]),
        (HEX_CLASSES, HEX_RECIP_MAG, [5, 5, 2]),
        (ORTHO_CLASSES, ORTHO_RECIP_MAG, [3, 2, 4]),
    ]
    for classes, recip_mag, mesh in cases:
        lower = mesh_climb.descend_one_rung(
            mesh, classes, recip_mag)
        assert lower != mesh          # each has a lower neighbour
        assert mesh_climb.climb_one_rung(
            lower, classes, recip_mag) == mesh


def test_descend_at_minimum_mesh_stays_put():
    """The [1,1,1] mesh has no lower neighbour, so descending it
    returns it unchanged rather than going below one point."""
    recovered = mesh_climb.descend_one_rung(
        [1, 1, 1], CUBIC_CLASSES, CUBIC_RECIP_MAG)
    assert recovered == [1, 1, 1]


def test_descend_known_hexagonal_rung():
    """The explicit inverse of the hexagonal sequence's last step:
    [5,5,2] descends to [4,4,2]."""
    recovered = mesh_climb.descend_one_rung(
        [5, 5, 2], HEX_CLASSES, HEX_RECIP_MAG)
    assert recovered == [4, 4, 2]


# ==================================================================
#  The ceiling (PSEUDOCODE 4e.2)
# ==================================================================

def test_at_ceiling_fires_when_any_axis_reaches_the_backstop():
    """The backstop bites when the LARGEST axis count reaches
    max_count, whatever the other axes are doing."""
    assert mesh_climb.at_ceiling([20, 4, 4], 20) is True
    assert mesh_climb.at_ceiling([21, 1, 1], 20) is True


def test_at_ceiling_false_below_the_backstop():
    """A mesh whose every axis is under max_count has not hit the
    ceiling and the climb may continue."""
    assert mesh_climb.at_ceiling([19, 19, 19], 20) is False
    assert mesh_climb.at_ceiling([5, 5, 5], 20) is False


# ==================================================================
#  Confidence-to-policy map (PSEUDOCODE 4e.4)
# ==================================================================

def test_confident_prediction_selects_a_parallel_grid():
    """A confidence at or above the threshold, not under-trained,
    yields a parallel grid with single-rung confirmation and no
    downward start offset."""
    policy = mesh_climb.resolve_climb_policy(
        confidence=0.9, under_trained=False)
    assert policy.mode == mesh_climb.PARALLEL_GRID
    assert policy.flat_needed == 1
    assert policy.grid_width == 1
    assert policy.start_offset == 0


def test_confidence_threshold_is_inclusive():
    """A confidence exactly at ``confidence_high`` counts as
    confident (the comparison is >=), so the boundary lands on the
    parallel grid."""
    thresholds = mesh_climb.DEFAULT_POLICY_THRESHOLDS
    policy = mesh_climb.resolve_climb_policy(
        confidence=thresholds.confidence_high, under_trained=False)
    assert policy.mode == mesh_climb.PARALLEL_GRID


def test_weak_prediction_selects_the_bracket_refine_climb():
    """A low confidence yields the default serial climb -- the
    bracket-refine shape: two flat rungs demanded, starting one rung
    below the seed, and carrying the stride cap the bracket phase
    reads."""
    policy = mesh_climb.resolve_climb_policy(
        confidence=0.4, under_trained=False)
    assert policy.mode == mesh_climb.BRACKET_REFINE
    assert policy.flat_needed == 2
    assert policy.grid_width == 0
    assert policy.start_offset == 1
    assert policy.max_stride == \
        mesh_climb.DEFAULT_POLICY_THRESHOLDS.max_stride


def test_curator_can_pin_the_unit_step_climb():
    """A curator who sets ``climb_shape`` to UNIT_STEP gets the fine
    walk-every-rung climb for a non-confident prediction instead of
    the bracketing default; a confident prediction still grids."""
    unit_step = mesh_climb.DEFAULT_POLICY_THRESHOLDS._replace(
        climb_shape=mesh_climb.UNIT_STEP)
    weak = mesh_climb.resolve_climb_policy(
        confidence=0.4, under_trained=False, thresholds=unit_step)
    assert weak.mode == mesh_climb.UNIT_STEP
    # The shape choice governs only the serial climb; a confident
    #   prediction still lays a parallel grid.
    strong = mesh_climb.resolve_climb_policy(
        confidence=0.95, under_trained=False, thresholds=unit_step)
    assert strong.mode == mesh_climb.PARALLEL_GRID


def test_climb_policy_from_empty_manifest_is_the_defaults():
    """An empty [harvest.kpoint_climb] sub-table yields the built-in
    provisional policy: DEFAULT_POLICY_THRESHOLDS and
    DEFAULT_MAX_COUNT unchanged."""
    thresholds, max_count = mesh_climb.climb_policy_from_manifest({})
    assert thresholds == mesh_climb.DEFAULT_POLICY_THRESHOLDS
    assert max_count == mesh_climb.DEFAULT_MAX_COUNT


def test_climb_policy_from_manifest_overrides_named_knobs():
    """A sub-table overrides only the knobs it names; the rest keep
    their provisional defaults."""
    thresholds, max_count = mesh_climb.climb_policy_from_manifest(
        {"confidence_high": 0.9, "max_count": 30})
    assert thresholds.confidence_high == 0.9          # overridden
    assert max_count == 30                            # overridden
    # An unnamed knob keeps its default.
    assert thresholds.flat_needed_cold == \
        mesh_climb.DEFAULT_POLICY_THRESHOLDS.flat_needed_cold


def test_climb_policy_from_manifest_overrides_stride_and_shape():
    """The bracket phase's ``max_stride`` cap and the curator's
    ``climb_shape`` choice are manifest knobs like the rest, merged
    over the provisional defaults."""
    thresholds, _ = mesh_climb.climb_policy_from_manifest(
        {"max_stride": 4, "climb_shape": mesh_climb.UNIT_STEP})
    assert thresholds.max_stride == 4
    assert thresholds.climb_shape == mesh_climb.UNIT_STEP
    # An unnamed knob still keeps its default.
    assert thresholds.confidence_high == \
        mesh_climb.DEFAULT_POLICY_THRESHOLDS.confidence_high


def test_climb_policy_rejects_an_unknown_climb_shape():
    """``climb_shape`` is the one knob with a restricted value, so a
    typo (a value that is not a known climb shape) fails loudly at
    resolve time rather than silently falling through to a default
    shape once the producer dispatches on it."""
    with pytest.raises(ValueError, match="climb_shape"):
        mesh_climb.climb_policy_from_manifest(
            {"climb_shape": "unit-step"})       # hyphen, not the enum
    # A valid shape is accepted.
    thresholds, _ = mesh_climb.climb_policy_from_manifest(
        {"climb_shape": mesh_climb.UNIT_STEP})
    assert thresholds.climb_shape == mesh_climb.UNIT_STEP


def test_climb_policy_merges_the_stride_flatness_multiple():
    """The bracket phase's looseness knob merges over the default like
    any other, and the default is a multiple greater than one (the
    bracket test is looser than convergence)."""
    assert mesh_climb.DEFAULT_POLICY_THRESHOLDS \
        .stride_flatness_multiple > 1
    thresholds, _ = mesh_climb.climb_policy_from_manifest(
        {"stride_flatness_multiple": 5.0})
    assert thresholds.stride_flatness_multiple == 5.0


def test_climb_policy_rejects_a_stride_multiple_below_one():
    """A multiple below one would make the bracket test STRICTER than
    convergence -- surely a mistake -- so it fails loudly, while a
    multiple of exactly one (bracket == convergence) is allowed."""
    with pytest.raises(ValueError, match="stride_flatness_multiple"):
        mesh_climb.climb_policy_from_manifest(
            {"stride_flatness_multiple": 0.5})
    thresholds, _ = mesh_climb.climb_policy_from_manifest(
        {"stride_flatness_multiple": 1.0})
    assert thresholds.stride_flatness_multiple == 1.0


def test_climb_policy_merges_the_metal_gap_threshold():
    """The metal test's gap threshold merges over the default, which
    is a small positive band gap in eV (above a true metal's near-zero
    gap, below any real insulator's)."""
    assert mesh_climb.DEFAULT_POLICY_THRESHOLDS \
        .metal_gap_threshold > 0
    thresholds, _ = mesh_climb.climb_policy_from_manifest(
        {"metal_gap_threshold": 0.1})
    assert thresholds.metal_gap_threshold == 0.1


def test_climb_policy_accepts_any_metal_gap_threshold():
    """Every real gap threshold is meaningful, so none is rejected.

    A NEGATIVE one is the documented way to disable the metal test for
    a diagnostic ladder -- no band gap can be negative, so the test can
    never fire (DESIGN 3.12.3 / 3.12.6).  A range check here would
    reject exactly the setting the design tells a curator to use.  Zero
    is meaningful too: a true metal's gap collapses to exactly zero, so
    a zero threshold is the strictest test that still fires on one."""
    thresholds, _ = mesh_climb.climb_policy_from_manifest(
        {"metal_gap_threshold": -1.0})
    assert thresholds.metal_gap_threshold == -1.0

    thresholds, _ = mesh_climb.climb_policy_from_manifest(
        {"metal_gap_threshold": 0.0})
    assert thresholds.metal_gap_threshold == 0.0


def test_climb_policy_merges_the_crystalline_floor_axis_count():
    """The crystalline floor cap merges over the default, which opens
    a crystalline climb above the ultra-coarse (Gamma) region -- more
    than a single k-point on its densest axis."""
    assert (mesh_climb.DEFAULT_POLICY_THRESHOLDS
            .crystalline_floor_axis_count > 1)
    thresholds, _ = mesh_climb.climb_policy_from_manifest(
        {"crystalline_floor_axis_count": 8})
    assert thresholds.crystalline_floor_axis_count == 8


def test_climb_policy_rejects_a_floor_axis_count_below_one():
    """The floor cap is a per-axis k-point count, so a value below
    one is a mistake and fails loudly."""
    with pytest.raises(ValueError,
                       match="crystalline_floor_axis_count"):
        mesh_climb.climb_policy_from_manifest(
            {"crystalline_floor_axis_count": 0})


def test_climb_policy_from_manifest_full_override():
    """Every knob named is carried into the resolved policy."""
    thresholds, max_count = mesh_climb.climb_policy_from_manifest(
        {"confidence_high": 0.7, "grid_width": 2,
         "start_offset_moderate": 1, "start_offset_cold": 3,
         "flat_needed_confident": 1, "flat_needed_cold": 2,
         "max_stride": 16, "climb_shape": mesh_climb.BRACKET_REFINE,
         "stride_flatness_multiple": 5.0,
         "metal_gap_threshold": 0.1,
         "crystalline_floor_axis_count": 8,
         "max_count": 24})
    assert thresholds == mesh_climb.PolicyThresholds(
        confidence_high=0.7, grid_width=2,
        start_offset_moderate=1, start_offset_cold=3,
        flat_needed_confident=1, flat_needed_cold=2,
        max_stride=16, climb_shape=mesh_climb.BRACKET_REFINE,
        stride_flatness_multiple=5.0, metal_gap_threshold=0.1,
        crystalline_floor_axis_count=8)
    assert max_count == 24


def test_under_trained_forces_a_lower_cold_climb():
    """An under-trained prediction climbs even when its confidence
    number is high, and starts lower than a merely-weak one (the
    bootstrap regime, DESIGN 7.9)."""
    policy = mesh_climb.resolve_climb_policy(
        confidence=0.95, under_trained=True)
    assert policy.mode == mesh_climb.BRACKET_REFINE
    assert policy.flat_needed == 2
    assert policy.start_offset == 2


# ==================================================================
#  First-round seeding (PSEUDOCODE 4e.4)
# ==================================================================

def test_climb_mode_seeds_one_mesh_below_the_prediction():
    """In a serial climb the first round is a single mesh, sitting
    ``start_offset`` rungs below the predicted seed, so it already
    has room to climb upward.  Density 125 on the cubic cell seeds
    [5,5,5]; one rung down is [4,4,4].  Both climb shapes open the
    same way, so this stands for BRACKET_REFINE and UNIT_STEP
    alike."""
    policy = mesh_climb.ClimbPolicy(
        mode=mesh_climb.BRACKET_REFINE, flat_needed=2, grid_width=0,
        start_offset=1, max_stride=8)
    meshes = mesh_climb.initial_meshes(
        125.0, policy, CUBIC_CLASSES, CUBIC_RECIP_MAG, 1.0)
    assert meshes == [[4, 4, 4]]
    # Climbing back up the start_offset returns to the seed.
    assert mesh_climb.climb_one_rung(
        meshes[0], CUBIC_CLASSES, CUBIC_RECIP_MAG) == [5, 5, 5]


def test_grid_mode_lays_a_symmetric_ladder_around_the_seed():
    """In PARALLEL_GRID mode the first round is the seed plus
    ``grid_width`` rungs each side, ascending, and each adjacent
    pair is climb-connected so the stop test reads across it."""
    policy = mesh_climb.ClimbPolicy(
        mode=mesh_climb.PARALLEL_GRID, flat_needed=1, grid_width=1,
        start_offset=0, max_stride=8)
    meshes = mesh_climb.initial_meshes(
        125.0, policy, CUBIC_CLASSES, CUBIC_RECIP_MAG, 1.0)
    assert meshes == [[4, 4, 4], [5, 5, 5], [6, 6, 6]]
    for lower, higher in zip(meshes, meshes[1:]):
        assert mesh_climb.climb_one_rung(
            lower, CUBIC_CLASSES, CUBIC_RECIP_MAG) == higher


def test_grid_collapses_duplicate_floor_meshes():
    """A grid whose seed sits at the minimal [1,1,1] mesh descends
    onto that same floor, and the duplicate is removed so no
    manufactured zero-delta pair reaches the stop test."""
    policy = mesh_climb.ClimbPolicy(
        mode=mesh_climb.PARALLEL_GRID, flat_needed=1, grid_width=1,
        start_offset=0, max_stride=8)
    meshes = mesh_climb.initial_meshes(
        1.0, policy, CUBIC_CLASSES, CUBIC_RECIP_MAG, 1.0)
    assert meshes == [[1, 1, 1], [2, 2, 2]]
    assert len(meshes) == len({tuple(m) for m in meshes})


# ==================================================================
#  The stride -- climb_n_rungs (PSEUDOCODE 4e.1)
# ==================================================================

def test_climb_n_rungs_is_repeated_single_steps():
    """A stride of ``n`` is exactly ``n`` single climbs: on the cubic
    cell three rungs from [2,2,2] reach [5,5,5], and a stride of one
    matches climb_one_rung."""
    assert mesh_climb.climb_n_rungs(
        [2, 2, 2], 3, CUBIC_CLASSES, CUBIC_RECIP_MAG) == [5, 5, 5]
    assert mesh_climb.climb_n_rungs(
        [2, 2, 2], 1, CUBIC_CLASSES, CUBIC_RECIP_MAG) \
        == mesh_climb.climb_one_rung(
            [2, 2, 2], CUBIC_CLASSES, CUBIC_RECIP_MAG)


def test_climb_n_rungs_zero_is_the_identity():
    """A stride of zero computes no new mesh and returns the input
    unchanged, so a degenerate stride cannot move the climb."""
    assert mesh_climb.climb_n_rungs(
        [3, 3, 1], 0, HEX_CLASSES, HEX_RECIP_MAG) == [3, 3, 1]


def test_climb_n_rungs_crosses_the_hexagonal_ladder():
    """On the anisotropic hexagonal ladder a stride crosses several
    distinct rungs at once: from [2,2,1] a stride of four lands on
    [5,5,2], skipping the three computed-only-at-endpoints meshes
    between (DESIGN 3.12.2)."""
    assert mesh_climb.climb_n_rungs(
        [2, 2, 1], 4, HEX_CLASSES, HEX_RECIP_MAG) == [5, 5, 2]


# ==================================================================
#  The up-to-ceiling walk -- ceiling_mesh (PSEUDOCODE 4e.3)
# ==================================================================

def test_ceiling_mesh_climbs_to_the_first_capped_mesh():
    """From below the cap, ceiling_mesh climbs to the lowest ladder
    position whose largest axis first reaches max_count: on the cubic
    cell with cap 5, both [2,2,2] and [3,3,3] resolve to [5,5,5]."""
    assert mesh_climb.ceiling_mesh(
        [2, 2, 2], CUBIC_CLASSES, CUBIC_RECIP_MAG, 5) == [5, 5, 5]
    assert mesh_climb.ceiling_mesh(
        [3, 3, 3], CUBIC_CLASSES, CUBIC_RECIP_MAG, 5) == [5, 5, 5]


def test_ceiling_mesh_returns_an_already_capped_mesh_unchanged():
    """A mesh already at the cap is its own ceiling mesh -- the walk
    stops before taking any step."""
    assert mesh_climb.ceiling_mesh(
        [5, 5, 5], CUBIC_CLASSES, CUBIC_RECIP_MAG, 5) == [5, 5, 5]


def test_ceiling_mesh_on_the_anisotropic_ladder():
    """On the hexagonal ladder the cap bites on the in-plane pair:
    from [4,4,2] with cap 5 the next rung [5,5,2] first reaches it."""
    assert mesh_climb.ceiling_mesh(
        [4, 4, 2], HEX_CLASSES, HEX_RECIP_MAG, 5) == [5, 5, 2]


# ==================================================================
#  Rung lookup -- rung_at (PSEUDOCODE 4e.3)
# ==================================================================

def test_rung_at_finds_the_computed_rung():
    """rung_at returns the ladder entry whose mesh matches, so a
    bracket endpoint's energy can be read back for the stride test."""
    rungs = [_rung([2, 2, 2]), _rung([4, 4, 4]), _rung([8, 8, 8])]
    assert mesh_climb.rung_at(rungs, [4, 4, 4]) is rungs[1]


def test_rung_at_raises_when_the_mesh_is_absent():
    """Asking for a mesh that has not been computed is a search-state
    bug (a stride tested before its endpoint landed), so rung_at
    raises loudly rather than returning None."""
    rungs = [_rung([2, 2, 2]), _rung([4, 4, 4])]
    with pytest.raises(ValueError):
        mesh_climb.rung_at(rungs, [3, 3, 3])


# ==================================================================
#  Filling the bracket -- next_fill_mesh (PSEUDOCODE 4e.3)
# ==================================================================

def test_next_fill_mesh_returns_the_lowest_uncomputed_rung():
    """With only the bracket endpoints [2,2,2] and [5,5,5] computed,
    filling walks up from the low end and asks for the lowest gap
    first: [3,3,3], then [4,4,4] once [3,3,3] lands, then None."""
    endpoints = [_rung([2, 2, 2]), _rung([5, 5, 5])]
    assert mesh_climb.next_fill_mesh(
        endpoints, [2, 2, 2], [5, 5, 5],
        CUBIC_CLASSES, CUBIC_RECIP_MAG) == [3, 3, 3]

    with_three = endpoints + [_rung([3, 3, 3])]
    assert mesh_climb.next_fill_mesh(
        with_three, [2, 2, 2], [5, 5, 5],
        CUBIC_CLASSES, CUBIC_RECIP_MAG) == [4, 4, 4]


def test_next_fill_mesh_none_when_the_interval_is_full():
    """When every ladder position in [lo, hi] is present the interval
    is filled and next_fill_mesh returns None, ending the fill."""
    filled = [_rung([2, 2, 2]), _rung([3, 3, 3]), _rung([4, 4, 4])]
    assert mesh_climb.next_fill_mesh(
        filled, [2, 2, 2], [4, 4, 4],
        CUBIC_CLASSES, CUBIC_RECIP_MAG) is None


# ==================================================================
#  The recorded flatness trace -- consecutive_block (PSEUDOCODE 4e.6)
# ==================================================================

def test_consecutive_block_drops_sparse_bracket_endpoints():
    """A bracket-refine ladder holds sparse stride endpoints below the
    filled bracket.  With a gap at [2,2,2] the block around the
    converged [4,4,4] keeps only the consecutive run [3,3,3], [4,4,4],
    [5,5,5] and drops the stray low [1,1,1] endpoint (DESIGN
    3.12.3)."""
    ladder = [_rung([1, 1, 1]), _rung([3, 3, 3]),
              _rung([4, 4, 4]), _rung([5, 5, 5])]
    block = mesh_climb.consecutive_block(
        ladder, _rung([4, 4, 4]), CUBIC_CLASSES, CUBIC_RECIP_MAG)
    assert [one.mesh for one in block] \
        == [[3, 3, 3], [4, 4, 4], [5, 5, 5]]


def test_consecutive_block_is_the_whole_ladder_when_consecutive():
    """A unit-step climb (or a grid) has no gaps, so the block around
    any interior rung is the entire ladder -- including the minimal
    [1,1,1], where the downward walk stops."""
    ladder = [_rung([1, 1, 1]), _rung([2, 2, 2]),
              _rung([3, 3, 3]), _rung([4, 4, 4])]
    block = mesh_climb.consecutive_block(
        ladder, _rung([3, 3, 3]), CUBIC_CLASSES, CUBIC_RECIP_MAG)
    assert [one.mesh for one in block] \
        == [[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]]


def test_consecutive_block_of_a_lone_rung_is_itself():
    """A converged rung with neither neighbour present is its own
    block, so a one-rung trace records just that mesh."""
    ladder = [_rung([4, 4, 4])]
    block = mesh_climb.consecutive_block(
        ladder, _rung([4, 4, 4]), CUBIC_CLASSES, CUBIC_RECIP_MAG)
    assert [one.mesh for one in block] == [[4, 4, 4]]
