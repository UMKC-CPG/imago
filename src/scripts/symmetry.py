## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

"""symmetry.py -- the structure toolchain's space-group operations
(ARCHITECTURE 2, 7; PSEUDOCODE 4b.4).

Why this module exists
----------------------
Reading a space group's point operations off disk is needed in two
places that must agree: the k-point-file writer in ``makeinput.py``
(which emits the operations and carries their fractional
translations for ``buildAtomPerm``) and the initial-potential
producer in ``build_initial_potentials.py`` (whose adaptive mesh
climb, DESIGN 3.12, must compute a cell's k-point axis classes in
Python *before* it dispatches any run).  If each read the operation
file with its own parser, the two could silently drift, and the
climb's later self-check -- do the producer's Python axis classes
match the operations imago actually ran under (DESIGN 2.7)? -- would
be meaningless.  So the read lives here, once.

This is the first focused module split out of the large
``structure_control.py`` at the symmetry seam that ARCHITECTURE 7
earmarks.  Its scope for now is only the on-disk operation reader;
the rest of the symmetry domain (supercell construction, atom
permutation tables) migrates here as that split proceeds.
"""

import os


def read_conv_abc_point_ops(space_group_db, space_group_name):
    """Read a space group's point group operations from the space
    group database (PSEUDOCODE 4b.4; DESIGN 2.7 / 3.2).

    The operations are returned in the on-disk *conventional-cell-abc
    fractional* form exactly as ``share/spaceDB/<sg>`` stores them --
    no change of basis is applied here.  Consumers rebase them as
    they need: imago conjugates them into the loaded cell at runtime
    (``computeRecipPointOps`` / ``computeRealPointOps``), and the
    producer's axis-class port mirrors that in Python
    (:func:`mesh_climb.axis_classes_for_cell`).

    The space group database file has this format:
      - Line 1: description (e.g. ``"F Fm3~m"``)
      - Line 2: root space group number, sub-number
      - Line 3: ``numSpaceOps numShifts``
      - Then for each of the ``numSpaceOps`` operations:
          a blank line, the 3x3 rotation matrix (3 lines), and the
          fractional translation vector (1 line)

    The first ``numSpaceOps / numShifts`` operations are the pure
    point group operations (no centering translations); those are
    the ones read.  For symmorphic space groups the translations are
    all zero; for non-symmorphic groups (e.g. Fd-3m) some operations
    carry non-zero translations (screw axes, glide planes), which the
    k-point writer needs for ``buildAtomPerm``'s real-space atom
    mapping.

    Parameters
    ----------
    space_group_db : str
        The space group database directory (``share/spaceDB``).  A
        loaded ``StructureControl`` exposes this as ``space_group_db``
        and ``makeinput``'s settings as ``space_db`` -- both callers
        adapt their own name to this argument.
    space_group_name : str
        The operation file within that directory (the canonical space
        group name, e.g. ``"Fm3~m"``).

    Returns
    -------
    tuple of (list, list)
        ``point_ops`` -- each element is a 3x3 rotation matrix (a
        list of 3 rows, each row a list of 3 floats).
        ``frac_trans`` -- each element is a list of 3 floats giving
        the fractional translation for that operation.
    """
    operation_file_path = os.path.join(
        space_group_db, space_group_name)
    with open(operation_file_path, "r") as operation_file:
        # Skip the description line and the root space-group number
        #   line; neither carries an operation.
        operation_file.readline()
        operation_file.readline()

        # The count line gives the total space operations and the
        #   number of centering shifts.  The pure point group is the
        #   first block of size (total / shifts) -- the operations
        #   before centering translations repeat them.
        count_fields = operation_file.readline().split()
        num_space_ops = int(count_fields[0])
        num_shifts = int(count_fields[1])
        num_point_ops = num_space_ops // num_shifts

        point_ops = []
        frac_trans = []
        for _operation_index in range(num_point_ops):
            # Each operation is preceded by a blank separator line.
            operation_file.readline()
            # Read the 3x3 rotation matrix, one row per line.
            rotation_matrix = []
            for _matrix_row in range(3):
                row_fields = operation_file.readline().split()
                rotation_matrix.append(
                    [float(value) for value in row_fields[:3]])
            # Read the fractional translation.  Non-symmorphic
            #   operations (screw axes, glide planes) carry non-zero
            #   values here; the k-point writer needs them for
            #   buildAtomPerm's real-space atom mapping.
            translation_fields = operation_file.readline().split()
            translation_vector = [
                float(value) for value in translation_fields[:3]]
            point_ops.append(rotation_matrix)
            frac_trans.append(translation_vector)

    return point_ops, frac_trans
