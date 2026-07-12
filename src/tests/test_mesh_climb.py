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

import pytest

import mesh_climb


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


def test_weak_prediction_selects_a_single_rung_climb():
    """A low confidence yields the CLIMB mode: one rung at a time,
    two flat rungs demanded, starting one rung below the seed."""
    policy = mesh_climb.resolve_climb_policy(
        confidence=0.4, under_trained=False)
    assert policy.mode == mesh_climb.CLIMB
    assert policy.flat_needed == 2
    assert policy.grid_width == 0
    assert policy.start_offset == 1


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


def test_climb_policy_from_manifest_full_override():
    """Every knob named is carried into the resolved policy."""
    thresholds, max_count = mesh_climb.climb_policy_from_manifest(
        {"confidence_high": 0.7, "grid_width": 2,
         "start_offset_moderate": 1, "start_offset_cold": 3,
         "flat_needed_confident": 1, "flat_needed_cold": 2,
         "max_count": 24})
    assert thresholds == mesh_climb.PolicyThresholds(
        confidence_high=0.7, grid_width=2,
        start_offset_moderate=1, start_offset_cold=3,
        flat_needed_confident=1, flat_needed_cold=2)
    assert max_count == 24


def test_under_trained_forces_a_lower_cold_climb():
    """An under-trained prediction climbs even when its confidence
    number is high, and starts lower than a merely-weak one (the
    bootstrap regime, DESIGN 7.9)."""
    policy = mesh_climb.resolve_climb_policy(
        confidence=0.95, under_trained=True)
    assert policy.mode == mesh_climb.CLIMB
    assert policy.flat_needed == 2
    assert policy.start_offset == 2


# ==================================================================
#  First-round seeding (PSEUDOCODE 4e.4)
# ==================================================================

def test_climb_mode_seeds_one_mesh_below_the_prediction():
    """In CLIMB mode the first round is a single mesh, sitting
    ``start_offset`` rungs below the predicted seed, so it already
    has room to climb upward.  Density 125 on the cubic cell seeds
    [5,5,5]; one rung down is [4,4,4]."""
    policy = mesh_climb.ClimbPolicy(
        mode=mesh_climb.CLIMB, flat_needed=2, grid_width=0,
        start_offset=1)
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
        start_offset=0)
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
        start_offset=0)
    meshes = mesh_climb.initial_meshes(
        1.0, policy, CUBIC_CLASSES, CUBIC_RECIP_MAG, 1.0)
    assert meshes == [[1, 1, 1], [2, 2, 2]]
    assert len(meshes) == len({tuple(m) for m in meshes})
