## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

"""test_imago_gamma_kp.py -- imago.py's gamma k-point detector.

Locks ``check_gamma_kp`` against the CURRENT k-point file grammar,
the one ``readKPoints`` (src/imago/kpoints.f90) actually reads:
every style shares four leading label/value pairs
(``KPOINT_STYLE_CODE``, ``KPOINT_INTG_CODE``,
``NUM_TETRA_DIAGONALS``, ``SYMMETRIZE_LAT_PARTIALS``) before the
style-specific fields.  The detector must locate every value by its
LABEL, never by a fixed line offset: offset parsing silently broke
when the two tetrahedron fields were added for all styles, and the
gamma executable stopped being selected without any visible error.
These tests are the tripwire for the next grammar change -- if a
field moves, the label lookup follows it; if a label disappears,
the detector must fail loudly rather than misread the file.

The companion file ``test_makeinput_gamma.py`` locks the OTHER side
of the contract: that makeinput writes a gamma request canonically
as a 1x1x1 style-1 mesh with zero shift.  This file locks that
imago.py recognizes such a file (and the other styles) correctly.

conftest.py's ``SCRIPTS_DIR`` insertion lets us import ``imago``
directly; each test writes a small k-point file into ``tmp_path``
and calls ``check_gamma_kp`` on it, the same call ``init_exes``
makes when choosing between imagoG and the general executable.
"""

import pytest

import imago


def _write_kp(tmp_path, body):
    """Write ``body`` as kp-scf.dat under tmp_path and return the
    (file name, directory) pair check_gamma_kp expects."""
    kp_file = tmp_path / "kp-scf.dat"
    kp_file.write_text(body)
    return "kp-scf.dat", str(tmp_path)


## The shared four-pair header every style code carries, with the
##   default tetrahedron settings makeinput writes:
_HEADER = """KPOINT_STYLE_CODE
{style}
KPOINT_INTG_CODE
0
NUM_TETRA_DIAGONALS
4
SYMMETRIZE_LAT_PARTIALS
0
"""


def test_style1_gamma_mesh_is_gamma(tmp_path):
    """A 1x1x1 style-1 mesh with zero shift -- the canonical gamma
    request makeinput writes -- must select the gamma executable."""
    body = _HEADER.format(style=1) + (
        "NUM_KP_A_B_C\n1 1 1\nKP_SHIFT_A_B_C\n0 0 0\n"
        "NUM_POINT_OPS\n48\n")
    assert imago.check_gamma_kp(*_write_kp(tmp_path, body)) is True


def test_style1_multipoint_mesh_is_not_gamma(tmp_path):
    """Any multi-point style-1 mesh must stay on the general
    executable."""
    body = _HEADER.format(style=1) + (
        "NUM_KP_A_B_C\n4 4 4\nKP_SHIFT_A_B_C\n0 0 0\n")
    assert imago.check_gamma_kp(*_write_kp(tmp_path, body)) is False


def test_style1_shifted_single_point_is_not_gamma(tmp_path):
    """A 1x1x1 mesh with a nonzero shift is one point NOT at the
    origin (a mean-value sample), so it is not gamma."""
    body = _HEADER.format(style=1) + (
        "NUM_KP_A_B_C\n1 1 1\nKP_SHIFT_A_B_C\n0.5 0.5 0.5\n")
    assert imago.check_gamma_kp(*_write_kp(tmp_path, body)) is False


def test_style0_single_origin_kpoint_is_gamma(tmp_path):
    """A style-0 explicit list holding exactly one kpoint at the
    origin is gamma."""
    body = _HEADER.format(style=0) + (
        "NUM_BLOCH_VECTORS\n1\nNUM_WEIGHT_KA_KB_KC\n"
        "1 1.0 0.0 0.0 0.0\n")
    assert imago.check_gamma_kp(*_write_kp(tmp_path, body)) is True


def test_style0_off_origin_kpoint_is_not_gamma(tmp_path):
    """One explicit kpoint away from the origin is not gamma."""
    body = _HEADER.format(style=0) + (
        "NUM_BLOCH_VECTORS\n1\nNUM_WEIGHT_KA_KB_KC\n"
        "1 1.0 0.5 0.0 0.0\n")
    assert imago.check_gamma_kp(*_write_kp(tmp_path, body)) is False


def test_style0_multiple_kpoints_are_not_gamma(tmp_path):
    """A multi-kpoint explicit list is not gamma, and the count
    test must decide that WITHOUT touching the kpoint lines."""
    body = _HEADER.format(style=0) + (
        "NUM_BLOCH_VECTORS\n2\nNUM_WEIGHT_KA_KB_KC\n"
        "1 0.5 0.0 0.0 0.0\n2 0.5 0.5 0.5 0.5\n")
    assert imago.check_gamma_kp(*_write_kp(tmp_path, body)) is False


def test_style1_nondefault_tetra_settings_do_not_confuse(tmp_path):
    """The tetrahedron fields sit between the style code and the
    mesh fields for EVERY style, and their values vary with the
    -scftetradiag option.  Changing them must not change gamma
    detection -- this is exactly the axis along which fixed-offset
    parsing broke."""
    body = _HEADER.format(style=1).replace(
        "NUM_TETRA_DIAGONALS\n4", "NUM_TETRA_DIAGONALS\n1").replace(
        "SYMMETRIZE_LAT_PARTIALS\n0", "SYMMETRIZE_LAT_PARTIALS\n1"
    ) + "NUM_KP_A_B_C\n1 1 1\nKP_SHIFT_A_B_C\n0 0 0\n"
    assert imago.check_gamma_kp(*_write_kp(tmp_path, body)) is True


def test_style2_density_mode_is_never_gamma(tmp_path):
    """Density mode defers mesh sizing to Imago at runtime, so the
    script routes it to the general executable unconditionally."""
    body = _HEADER.format(style=2) + (
        "MIN_KP_LINE_DENSITY\n200.0\nKP_SHIFT_A_B_C\n0 0 0\n")
    assert imago.check_gamma_kp(*_write_kp(tmp_path, body)) is False


def test_missing_label_fails_loudly(tmp_path):
    """A file missing an expected label does not follow the grammar
    this script understands.  It must exit with an error naming the
    label, never silently return False -- a silent False is the
    failure mode that hid the offset drift."""
    body = _HEADER.format(style=1) + "KP_SHIFT_A_B_C\n0 0 0\n"
    with pytest.raises(SystemExit):
        imago.check_gamma_kp(*_write_kp(tmp_path, body))
