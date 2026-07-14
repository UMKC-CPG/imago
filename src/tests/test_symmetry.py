"""Tests for the shared space-group operation reader
(``symmetry.read_conv_abc_point_ops``; PSEUDOCODE 4b.4; DESIGN
2.7 / 3.2) and the ``StructureControl.point_ops`` accessor over it.

The reader is the single parser of a ``share/spaceDB/<sg>``
operation file, shared by the k-point-file writer (``makeinput``)
and the initial-potential producer (whose adaptive mesh climb
needs a cell's point operations in Python before it dispatches any
run).  These tests pin the parse itself: the header is skipped,
only the pure point group -- the first ``numSpaceOps / numShifts``
operations -- is read, each rotation lands as 3x3 rows, and
non-symmorphic fractional translations are preserved.

The self-contained cases write a tiny operation file in the spaceDB
format and need no ``$IMAGO_DATA``; the integration cases read a
real database file and exercise the accessor through a loaded
``StructureControl``.
"""

import os

import pytest

import symmetry


# The identity rotation and a distinct four-fold rotation about the
#   c axis, so a test can tell operation 2 apart from operation 1.
_IDENTITY = [[1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0],
             [0.0, 0.0, 1.0]]
_ROT4_Z = [[0.0, -1.0, 0.0],
           [1.0, 0.0, 0.0],
           [0.0, 0.0, 1.0]]


def _write_space_group_file(path, num_space_ops, num_shifts,
                            operations):
    """Write one operation file in the spaceDB text format.

    ``operations`` is a list of ``(rotation_rows, translation)``
    pairs, each rotation a 3x3 list of rows.  The description and
    root-number header lines carry placeholder values the reader
    skips; the count line carries ``num_space_ops num_shifts``.
    Numbers are formatted the way the real database writes them, but
    the reader only splits on whitespace, so the exact spacing does
    not matter.
    """
    lines = ["P Test", "1 1", f"{num_space_ops} {num_shifts}"]
    for rotation_rows, translation in operations:
        lines.append("")                     # blank operation separator
        for row in rotation_rows:
            lines.append(" ".join(f"{value:.8f}" for value in row))
        lines.append(
            " ".join(f"{value:.8f}" for value in translation))
    path.write_text("\n".join(lines) + "\n")


def test_reads_single_identity_operation(tmp_path):
    """A one-operation, symmorphic (P1-like) file yields exactly the
    identity rotation and a zero translation."""
    _write_space_group_file(
        tmp_path / "P1", num_space_ops=1, num_shifts=1,
        operations=[(_IDENTITY, [0.0, 0.0, 0.0])])

    point_ops, frac_trans = symmetry.read_conv_abc_point_ops(
        str(tmp_path), "P1")

    assert point_ops == [_IDENTITY]
    assert frac_trans == [[0.0, 0.0, 0.0]]


def test_reads_only_the_pure_point_group_block(tmp_path):
    """With ``numSpaceOps / numShifts`` giving two, only the first
    two operations are read; the trailing centering copies below
    them are not, and a non-symmorphic translation is preserved."""
    # Four total operations, two centering shifts -> two pure point
    #   ops.  Write four and expect the reader to stop after two.
    centering_copy = (_IDENTITY, [0.5, 0.5, 0.0])
    _write_space_group_file(
        tmp_path / "Ctest", num_space_ops=4, num_shifts=2,
        operations=[(_IDENTITY, [0.0, 0.0, 0.0]),
                    (_ROT4_Z, [0.5, 0.0, 0.0]),
                    centering_copy, centering_copy])

    point_ops, frac_trans = symmetry.read_conv_abc_point_ops(
        str(tmp_path), "Ctest")

    assert point_ops == [_IDENTITY, _ROT4_Z]
    # Operation 2's screw-axis translation survives the read; the
    #   two centering copies are never reached.
    assert frac_trans == [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]


@pytest.mark.integration
def test_reads_a_real_database_file(imago_data_dir):
    """The reader parses a real ``share/spaceDB`` file: space group
    1 (``1_a``, P1) has a single identity operation."""
    space_group_db = os.path.join(imago_data_dir, "spaceDB")

    point_ops, frac_trans = symmetry.read_conv_abc_point_ops(
        space_group_db, "1_a")

    assert len(point_ops) == 1
    assert point_ops[0] == _IDENTITY
    assert frac_trans[0] == [0.0, 0.0, 0.0]


@pytest.mark.integration
def test_structure_control_accessor_returns_rotations(sc_c2_molecule):
    """``StructureControl.point_ops`` delegates to the shared reader,
    keyed on the loaded cell's own space group, and returns just the
    rotations.  The C2 molecule loads in P1, so its only operation is
    the identity."""
    rotations = sc_c2_molecule.point_ops()

    assert rotations == [_IDENTITY]
