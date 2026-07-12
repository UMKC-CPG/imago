"""mesh_climb.py -- the producer-side k-point mesh primitives and
the adaptive-climb search over them (DESIGN 3.12; PSEUDOCODE 4c /
4e).

Why this module exists
----------------------
The historical-guidance producer (``build_initial_potentials.py``)
converges each reference solid's k-point sampling before it
harvests a potential.  DESIGN 3.12 replaces the old fixed grid of
requested *densities* with an adaptive climb through symmetry-
compatible *meshes*: the density is the transferable currency the
guidance database is keyed on, but the mesh is the natural unit
for *searching*, because for a high-symmetry cell only a handful of
distinct meshes exist and each is a genuinely different
calculation.

imago (the Fortran engine) already knows how to turn a density
into a mesh and reduce it to the irreducible wedge -- but at
runtime, one cell at a time.  The producer reasons about meshes
for many materials *before* it dispatches anything (to seed and
step the climb), so it needs the same mesh arithmetic in Python.
This module is that arithmetic, re-expressed on plain count
vectors:

- **Axis classes** (4c.1 / 4c.7): which reciprocal axes the point
  group couples, so they must share a count.  Sourced from the
  space-group operations and the loaded cell, mirroring imago's
  ``computeRealPointOps -> computeRecipPointOps ->
  computeAxisClasses`` (DESIGN 2.7).
- **Axial count selection** (4c.2): the most isotropic integer
  mesh, per class, that meets a density floor.
- **Rung mechanics** (4e.1): the one-step-up / one-step-down
  moves the climb takes through the distinct meshes.

The climb's control loop (``converge_by_climb``) and the two
dispatch modes live in this module too, but they take the actual
run-a-round action as an injected callback, so nothing here
imports the dispatch layer -- these functions stay pure and
unit-testable (VISION Principle 9/12).

Conventions
-----------
Count vectors, reciprocal magnitudes, and axis-class labels are
plain 0-indexed length-3 Python lists ``[a, b, c]`` -- NOT the
1-indexed sentinel-slot layout the Fortran uses.  A 3x3 matrix is
a list of three rows, each a list of three floats.  Lattices are
passed with their vectors as ROWS (the ``StructureControl``
layout).
"""

from collections import namedtuple


# ==================================================================
#  Small 3x3 linear algebra
#
#  Hand-rolled to keep this module free of a numpy dependency, in
#  the same spirit as ``structure_control.make_inv_or_recip_lattice``
#  (which inverts 3x3 lattices by cofactor expansion).  The guidance
#  producer chain is deliberately import-light so the test suite
#  starts fast on the networked venv.
# ==================================================================

def _transpose(matrix):
    """Return the transpose of a 3x3 matrix (rows <-> columns)."""
    return [[matrix[column][row] for column in range(3)]
            for row in range(3)]


def _mat_mul(left, right):
    """Return the 3x3 matrix product ``left * right``."""
    return [[sum(left[row][k] * right[k][column] for k in range(3))
             for column in range(3)]
            for row in range(3)]


def _determinant3x3(matrix):
    """Return the determinant of a 3x3 matrix by cofactor
    expansion along the first row."""
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2]
                        - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2]
                          - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1]
                          - matrix[1][1] * matrix[2][0]))


def _inverse3x3(matrix):
    """Return the inverse of a 3x3 matrix via the adjugate divided
    by the determinant.  Raises ``ValueError`` on a singular
    matrix -- for the lattices and point operations this module
    handles, that signals a degenerate cell or a corrupt space
    group rather than a recoverable condition."""
    determinant = _determinant3x3(matrix)
    if abs(determinant) < 1.0e-12:
        raise ValueError(
            "cannot invert a (near-)singular 3x3 matrix; the cell "
            "or space-group operation is degenerate")
    # Cofactor matrix, then transpose (the adjugate), then scale.
    cofactor = [
        [(matrix[(row + 1) % 3][(col + 1) % 3]
          * matrix[(row + 2) % 3][(col + 2) % 3]
          - matrix[(row + 1) % 3][(col + 2) % 3]
          * matrix[(row + 2) % 3][(col + 1) % 3])
         for col in range(3)]
        for row in range(3)]
    adjugate = _transpose(cofactor)
    return [[adjugate[row][col] / determinant for col in range(3)]
            for row in range(3)]


def _round_to_integer_matrix(matrix, tolerance=1.0e-6):
    """Round a 3x3 matrix to the nearest integer entries, but only
    when every entry is already within ``tolerance`` of an integer.

    A valid lattice automorphism (a point-group operation expressed
    in a lattice basis) is an integer matrix of determinant +/-1
    (DESIGN 2.7).  The conjugation ``T^-1 * R * T`` therefore lands
    on integers up to floating-point round-off, which this cleans.
    An entry that is NOT near an integer after conjugation means the
    inputs do not describe a genuine lattice automorphism -- a
    corrupt cell or space group -- so we raise rather than silently
    round a bogus operation into the axis-class union-find."""
    cleaned = [[0] * 3 for _ in range(3)]
    for row in range(3):
        for col in range(3):
            nearest = round(matrix[row][col])
            if abs(matrix[row][col] - nearest) > tolerance:
                raise ValueError(
                    "point-group operation is not integer in the "
                    "lattice basis after conjugation "
                    f"(entry {matrix[row][col]!r}); the cell or "
                    "space group is inconsistent")
            cleaned[row][col] = int(nearest)
    return cleaned


# ==================================================================
#  Axis classes (DESIGN 3.8; PSEUDOCODE 4c.1, 4c.7)
# ==================================================================

def _union_find_root(parent, node):
    """Return the representative root of ``node`` with path
    compression, over a 0-indexed ``parent`` array."""
    while parent[node] != node:
        parent[node] = parent[parent[node]]
        node = parent[node]
    return node


def compute_axis_classes(recip_point_ops):
    """Return the axis-class label of each reciprocal axis
    (PSEUDOCODE 4c.1; DESIGN 3.8).

    Two reciprocal axes are *coupled* when some operation connects
    them off the diagonal; the transitive closure of coupling gives
    the classes that must share a k-point count.  The reciprocal-
    space operations are used because the k-point mesh transforms
    under them.  The returned list has three entries, each the
    union-find root of its axis; axes with equal labels are one
    class.  (The labels themselves are arbitrary; only their
    equality pattern is meaningful.)

    ``recip_point_ops`` is a sequence of 3x3 integer matrices.
    """
    parent = [0, 1, 2]
    for operation in recip_point_ops:
        for row in range(3):
            for column in range(3):
                if row != column and operation[row][column] != 0:
                    root_row = _union_find_root(parent, row)
                    root_column = _union_find_root(parent, column)
                    if root_row != root_column:
                        parent[root_column] = root_row
    return [_union_find_root(parent, axis) for axis in range(3)]


def axis_classes_for_cell(conv_abc_point_ops, loaded_lattice,
                          conv_lattice, cell_mode):
    """Return a cell's axis classes, sourced the way the producer
    must (PSEUDOCODE 4c.7; DESIGN 2.7).

    This is the Python mirror of imago's runtime chain
    ``computeRealPointOps -> computeRecipPointOps ->
    computeAxisClasses``.  The producer needs a cell's axis classes
    *before* it dispatches any run (to seed and step the climb), so
    it cannot read them back from a run; it recomputes them from the
    same two ingredients imago uses.

    Parameters
    ----------
    conv_abc_point_ops
        The space-group rotation matrices as stored on disk in
        ``share/spaceDB/<sg>`` -- conventional-cell-abc fractional
        form, each acting as ``r' = R*r`` (the producer already
        reads these to write the kp file).  A sequence of 3x3
        matrices.
    loaded_lattice, conv_lattice
        The loaded and conventional cells, 3x3, with their lattice
        vectors as ROWS (the ``StructureControl`` layout), in a
        common length unit.  Only their ratio enters, so the unit
        cancels.
    cell_mode
        ``"full"`` when the loaded cell IS the conventional cell
        (the common case) -- the change of basis is the identity
        and the conjugation collapses to a copy.  ``"prim"`` when
        the structure loaded a primitive reduction, where the
        conventional-abc operations must be conjugated into the
        loaded basis before they mean anything (DESIGN 2.7).

    Returns
    -------
    list of int
        The three axis-class labels (see ``compute_axis_classes``).
    """
    if cell_mode == "full":
        change_of_basis = _identity3x3()
        inverse_change = _identity3x3()
    elif cell_mode == "prim":
        # Vectors-as-columns are the transpose of the row-major
        #   storage.  T = Lc^-1 * L carries loaded-fractional to
        #   conventional-fractional coordinates (DESIGN 2.7).
        conv_columns = _transpose(conv_lattice)
        loaded_columns = _transpose(loaded_lattice)
        change_of_basis = _mat_mul(
            _inverse3x3(conv_columns), loaded_columns)
        inverse_change = _inverse3x3(change_of_basis)
    else:
        raise ValueError(
            f"cell_mode must be 'full' or 'prim', got "
            f"{cell_mode!r}")

    recip_point_ops = []
    for conventional_op in conv_abc_point_ops:
        # Direct-space op in the loaded basis: an ordinary
        #   similarity, because direct fractional coordinates are
        #   covariant (DESIGN 2.7).
        loaded_op = _mat_mul(
            inverse_change, _mat_mul(conventional_op,
                                     change_of_basis))
        # Reciprocal-space twin: the inverse transpose, because
        #   k-point fractional coordinates are contravariant.  This
        #   is the operation the mesh actually folds under.
        recip_op = _transpose(_inverse3x3(loaded_op))
        recip_point_ops.append(_round_to_integer_matrix(recip_op))

    return compute_axis_classes(recip_point_ops)


def _identity3x3():
    """Return a fresh 3x3 identity matrix."""
    return [[1 if row == column else 0 for column in range(3)]
            for row in range(3)]


# ==================================================================
#  Axial count selection (DESIGN 3.7; PSEUDOCODE 4c.2)
# ==================================================================

def distinct_classes(classes):
    """Return the distinct axis-class labels in first-seen order.

    Order-stable so callers that iterate classes (the density-floor
    loop, the climb's rung choice) behave deterministically."""
    seen = []
    for label in classes:
        if label not in seen:
            seen.append(label)
    return seen


def _class_members(classes, label):
    """Return the axis indices belonging to class ``label``."""
    return [axis for axis in range(3) if classes[axis] == label]


def spacing_spread(counts, recip_mag):
    """Return how uneven the three inter-point spacings are for a
    mesh (PSEUDOCODE 4c.2).

    The spacing along reciprocal axis ``i`` is ``|b_i| / n_i``; the
    spread is ``max - min`` over the three.  A smaller spread is a
    more isotropic mesh, which is what the count selection and the
    climb both drive toward when they choose which class to bump."""
    spacings = [recip_mag[axis] / counts[axis] for axis in range(3)]
    return max(spacings) - min(spacings)


def select_axial_counts(density, recip_mag, recip_cell_volume,
                        classes):
    """Return the most isotropic symmetry-compatible integer mesh
    that meets a density floor (PSEUDOCODE 4c.2; DESIGN 3.7).

    ``density`` is the target volume density D (full-mesh points per
    unit reciprocal volume).  Axes in one class share a count by
    construction, so the result never breaks symmetry compatibility
    (DESIGN 3.8).  A non-positive density is the Gamma sentinel
    (DESIGN 3.6) and resolves to a single point.

    Parameters
    ----------
    density
        The target volume density D.
    recip_mag
        The three reciprocal-axis magnitudes ``|b_i|``.
    recip_cell_volume
        The reciprocal cell volume, in units consistent with
        ``recip_mag`` (so ``product(|b_i|) / recip_cell_volume`` is
        dimensionless-ish and the density floor is comparable).
    classes
        The axis-class labels (from ``axis_classes_for_cell`` /
        ``compute_axis_classes``).
    """
    if density <= 0:
        return [1, 1, 1]

    # Continuous isotropic counts at a common spacing h:
    #   h = (prod|b_i| / (recipCellVolume * D))^(1/3),
    #   x_i = |b_i| / h.
    spacing = (recip_mag[0] * recip_mag[1] * recip_mag[2]
               / (recip_cell_volume * density)) ** (1.0 / 3.0)
    continuous = [recip_mag[axis] / spacing for axis in range(3)]

    # Force one shared real count per class.  Coupled axes already
    #   have equal |b_i| (hence equal x); the mean guards against
    #   round-off before rounding.
    for label in distinct_classes(classes):
        members = _class_members(classes, label)
        shared = sum(continuous[axis] for axis in members) \
            / len(members)
        for axis in members:
            continuous[axis] = shared

    # Nearest positive integer, per class (already equal within a
    #   class, so the class stays uniform).
    counts = [max(1, round(continuous[axis])) for axis in range(3)]

    # Raise WHOLE classes until the full-mesh product meets the
    #   floor.  Never raise a single axis inside a multi-axis class
    #   -- that would break symmetry compatibility.
    floor = density * recip_cell_volume
    while (counts[0] * counts[1] * counts[2]) < floor:
        best_label = min(
            distinct_classes(classes),
            key=lambda label: spacing_spread(
                bump(counts, classes, label, +1), recip_mag))
        counts = bump(counts, classes, best_label, +1)

    return counts


# ==================================================================
#  Rung mechanics (DESIGN 3.12.2; PSEUDOCODE 4e.1)
# ==================================================================

def bump(counts, classes, label, step):
    """Return ``counts`` with ``step`` added to every axis in class
    ``label`` (PSEUDOCODE 4e.1).

    Moving a whole class at once keeps its axes equal, so the mesh
    stays symmetry-compatible (DESIGN 3.8).  Does not mutate the
    input."""
    return [count + step if classes[axis] == label else count
            for axis, count in enumerate(counts)]


def climb_one_rung(counts, classes, recip_mag):
    """Return the next distinct mesh up the climb (PSEUDOCODE 4e.1;
    DESIGN 3.12.2).

    One step of the count-selection floor loop: increment the axis
    class whose bump most evens the three inter-point spacings
    (smallest ``spacing_spread``).  Because the class structure
    comes from the point group, this single rule produces the right
    sequence for every crystal system -- cubic in lockstep, a
    hexagonal cell's in-plane pair several steps per c-axis bump, an
    orthorhombic cell one axis at a time -- with no per-material
    table."""
    best_label = min(
        distinct_classes(classes),
        key=lambda label: spacing_spread(
            bump(counts, classes, label, +1), recip_mag))
    return bump(counts, classes, best_label, +1)


def descend_one_rung(counts, classes, recip_mag):
    """Return a mesh one rung below ``counts`` on the climb
    (PSEUDOCODE 4e.1; DESIGN 3.12.4).

    A lower mesh that ``climb_one_rung`` steps back up to ``counts``
    from.  A mesh can have more than one such lower neighbour -- a
    hexagonal [4,4,2] is reached both from [4,4,1] (bumping the c
    axis) and from [3,3,2] (bumping the in-plane pair) -- so this
    returns the first in class order.  That guarantees
    ``climb_one_rung(descend_one_rung(counts)) == counts``, but it
    does not promise to undo any particular climb path.  Used to
    seed the climb a rung or more below the predicted mesh (so it
    acquires a lower neighbour) and to lay out the confident mode's
    grid.  Returns ``counts`` unchanged when it is already the
    minimal mesh the cell admits."""
    for label in distinct_classes(classes):
        members = _class_members(classes, label)
        # Cannot take any axis of this class below one point.
        if any(counts[axis] == 1 for axis in members):
            continue
        trial = bump(counts, classes, label, -1)
        if climb_one_rung(trial, classes, recip_mag) == counts:
            return trial
    return counts


# ==================================================================
#  Stopping, seeding, and the dispatch-mode policy
#  (DESIGN 3.12.3-3.12.6; PSEUDOCODE 4e.2 / 4e.4)
#
#  These decide the SHAPE of the climb: when a mesh is high enough
#  that the search must stop (the ceiling), what to run in the very
#  first round, and how the predictor's confidence chooses between
#  a single climbing rung and a small parallel grid.  The energy-
#  based stop test itself (``pick_converged_climb``) lives beside
#  ``per_atom_ev`` in ``guidance_harvest``, single-sourced against
#  the harvest's own convergence rule so the two cannot drift; the
#  round loop that reads those energies and calls ``climb_one_rung``
#  lives in the producer (ARCHITECTURE 9.7).  What stays here is the
#  pure mesh-and-confidence arithmetic those two callers lean on.
# ==================================================================

# The two dispatch modes (DESIGN 3.12.5).  A confident prediction
#   lays a small grid around the seed and judges it in one round; a
#   cold or moderate one climbs a single rung at a time.
PARALLEL_GRID = "parallel_grid"
CLIMB = "climb"


## The confidence-derived shape of one material's climb, produced by
##   resolve_climb_policy and read by initial_meshes (4e.4) and the
##   producer's per-material decision (climbAction, 4e.3):
##     mode          PARALLEL_GRID or CLIMB -- which first round to
##                   run (a small grid, or a single climbing rung)
##     flat_needed   consecutive flat interior rungs the stop test
##                   must see before it accepts convergence (4e.2)
##     grid_width    rungs laid on EACH side of the seed in a
##                   parallel grid (ignored in CLIMB mode)
##     start_offset  rungs BELOW the seed a climb begins from, so the
##                   first rung gains a lower neighbour (ignored for
##                   a parallel grid)
ClimbPolicy = namedtuple(
    "ClimbPolicy",
    ["mode", "flat_needed", "grid_width", "start_offset"])


## The tunable knobs the confidence-to-policy map reads.  Their
##   numeric values are to be fixed by the seed experiment (DESIGN
##   3.12.6); DEFAULT_POLICY_THRESHOLDS carries documented
##   provisional values, and the manifest characterisation block
##   will override them (that config wiring is later work):
##     confidence_high        at or above this predictor confidence
##                            the search is confident (parallel grid)
##     grid_width             rungs each side of the seed when
##                            confident
##     start_offset_moderate  rungs below the seed for a low-but-
##                            trained prediction's climb
##     start_offset_cold      rungs below for an under-trained /
##                            bootstrap climb (starts lower, wider)
##     flat_needed_confident  persistence demanded when confident (1)
##     flat_needed_cold       persistence demanded when cold (2)
PolicyThresholds = namedtuple(
    "PolicyThresholds",
    ["confidence_high", "grid_width", "start_offset_moderate",
     "start_offset_cold", "flat_needed_confident",
     "flat_needed_cold"])


# Provisional defaults (DESIGN 3.12.6): placeholders until the seed
#   experiment tunes them, and the fallback the manifest layer will
#   override.  A confident search lays a +/-1-rung grid and accepts a
#   single flat interior rung; a cold climb starts two rungs low and
#   demands two consecutive flat rungs.
DEFAULT_POLICY_THRESHOLDS = PolicyThresholds(
    confidence_high=0.75,
    grid_width=1,
    start_offset_moderate=1,
    start_offset_cold=2,
    flat_needed_confident=1,
    flat_needed_cold=2)


# The fixed per-axis backstop the climb never exceeds (DESIGN
#   3.12.3): a provisional ceiling pending the seed experiment, and
#   ahead of the resource-dataspace cost ceiling that will layer on
#   later.  The climb stops at whichever ceiling bites first.
DEFAULT_MAX_COUNT = 20


def at_ceiling(mesh, max_count):
    """Return whether a mesh has reached the per-axis backstop
    (PSEUDOCODE 4e.2; DESIGN 3.12.3).

    The climb stops -- non-converged -- once any single axis count
    reaches ``max_count``.  This fixed backstop guarantees the
    search terminates even for a cell whose energy never quite goes
    flat; the resource dataspace's cost ceiling layers on later, and
    the climb halts at whichever ceiling bites first."""
    return max(mesh) >= max_count


def resolve_climb_policy(confidence, under_trained,
                         thresholds=DEFAULT_POLICY_THRESHOLDS):
    """Turn a prediction's confidence into the shape of its climb
    (PSEUDOCODE 4e.4; DESIGN 3.12.4-3.12.6).

    A confident prediction (``confidence`` at or above
    ``thresholds.confidence_high`` and not flagged under-trained)
    warrants a short, tight search: lay a small parallel grid around
    the seed and accept a single flat interior rung.  A weak or
    under-trained prediction warrants a wider one: climb a single
    rung at a time, beginning below the seed and demanding the
    flatness persist, so one lucky flat step cannot end it early.
    An under-trained prediction (the bootstrap regime, DESIGN 7.9)
    begins the climb lower still; the producer additionally seeds it
    from the wide-grid floor rather than a predicted density, but
    that choice of seed density is the producer's, not this
    policy's."""
    confident = (not under_trained
                 and confidence >= thresholds.confidence_high)
    if confident:
        return ClimbPolicy(
            mode=PARALLEL_GRID,
            flat_needed=thresholds.flat_needed_confident,
            grid_width=thresholds.grid_width,
            start_offset=0)

    start_offset = (thresholds.start_offset_cold if under_trained
                    else thresholds.start_offset_moderate)
    return ClimbPolicy(
        mode=CLIMB,
        flat_needed=thresholds.flat_needed_cold,
        grid_width=0,
        start_offset=start_offset)


def climb_policy_from_manifest(climb_settings):
    """Merge a manifest ``[harvest.kpoint_climb]`` sub-table over the
    provisional defaults (PSEUDOCODE 4e.4; DESIGN 5.7 / 3.12.6).

    Returns ``(thresholds, max_count)``: a ``PolicyThresholds`` the
    confidence-to-mode policy reads (``resolve_climb_policy``) and
    the per-axis ceiling ``at_ceiling`` reads.  Every knob is
    optional -- an omitted one keeps its provisional default
    (``DEFAULT_POLICY_THRESHOLDS`` / ``DEFAULT_MAX_COUNT``, whose
    values are still to be fixed by the seed experiment, 3.12.6) --
    so an empty ``climb_settings`` yields the built-in policy and a
    partial one overrides only the knobs it names.

    ``climb_settings`` is the plain dict the manifest reader parsed
    from ``[harvest.kpoint_climb]``; its keys were already validated
    against the known knob names at load
    (``curation_manifest.KPOINT_CLIMB_KEYS``), so only the merge
    remains here."""
    threshold_overrides = {
        field: climb_settings[field]
        for field in PolicyThresholds._fields
        if field in climb_settings}
    thresholds = DEFAULT_POLICY_THRESHOLDS._replace(
        **threshold_overrides)
    max_count = climb_settings.get("max_count", DEFAULT_MAX_COUNT)
    return thresholds, max_count


def _distinct_meshes(meshes):
    """Return ``meshes`` with exact duplicates removed and first-
    seen order preserved.  A parallel grid whose seed sits near the
    minimal mesh can descend onto the same floor mesh more than once
    (repeated [1,1,1], say); those collapse to a single rung so the
    grid does not carry a manufactured zero-energy-delta pair into
    the stop test."""
    unique = []
    for mesh in meshes:
        if mesh not in unique:
            unique.append(mesh)
    return unique


def initial_meshes(density, policy, classes, recip_mag,
                   recip_cell_volume):
    """Return the mesh or meshes to run in the climb's first round
    (PSEUDOCODE 4e.4; DESIGN 3.12.4-3.12.5).

    ``density`` is the seed density -- the guidance prediction, or
    the wide-grid floor for an under-trained bootstrap (DESIGN 7.9).
    It is converted to a starting mesh by ``select_axial_counts``
    (4c.2), and ``policy`` (from ``resolve_climb_policy``) sets the
    round's shape:

    - ``PARALLEL_GRID``: the seed plus ``policy.grid_width`` rungs on
      each side, laid down together and judged as one grid.  The grid
      is a genuine climb ladder -- ``climb_one_rung`` maps each lower
      rung to the next (4e.1) -- so the two-sided stop test reads
      straight across it.
    - ``CLIMB``: a single mesh ``policy.start_offset`` rungs below
      the seed, so the first rung already has room to climb upward
      and acquire the upper neighbour the stop test needs.

    Returns a list of axial-count meshes in ascending order, with any
    duplicate floor meshes collapsed (``_distinct_meshes``)."""
    seed = select_axial_counts(density, recip_mag,
                               recip_cell_volume, classes)

    if policy.mode == PARALLEL_GRID:
        meshes = [seed]
        lower = seed
        upper = seed
        for _ in range(policy.grid_width):
            lower = descend_one_rung(lower, classes, recip_mag)
            upper = climb_one_rung(upper, classes, recip_mag)
            meshes = [lower] + meshes + [upper]
        return _distinct_meshes(meshes)

    start = seed
    for _ in range(policy.start_offset):
        start = descend_one_rung(start, classes, recip_mag)
    return [start]
