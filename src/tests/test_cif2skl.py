"""Tests for cif2skl.py -- the space-group-preserving CIF to imago.skl
converter (ARCHITECTURE 9.5).

These exercise the real conversion path: ASE parses a committed CIF
fixture, cif2skl resolves the spaceDB setting by running the structure
through ``apply_space_group`` (the applySpaceGroup binary) and matching
ASE's symmetry expansion, and writes the asymmetric-unit skeleton.  They
therefore need ASE (skipped if absent) and the same IMAGO_DATA /
applySpaceGroup binary the StructureControl symmetry tests already
require.

The two silicon fixtures are the same phase (diamond-cubic, Fd-3m) in
the two different origin settings COD publishes -- origin 1 and origin 2.
They are the discrimination test: a correct converter must resolve them
to *different* spaceDB tokens, because only the matching origin's
operations regenerate each one's coordinates.
"""

import os

import pytest

pytest.importorskip("ase")

import cif2skl
from cif2skl import convert, cif_to_skeleton, CifConversionError

CIF_DIR = os.path.join(os.path.dirname(__file__), "fixtures", "cif")


def _cif(name):
    return os.path.join(CIF_DIR, name)


def _count_asymmetric_atoms(skl_path):
    """Return the atom count declared on the skeleton's ``fract`` line --
    the size of the stored asymmetric unit."""

    with open(skl_path) as handle:
        for line in handle:
            tokens = line.split()
            if tokens and tokens[0] == "fract":
                return int(tokens[1])
    raise AssertionError(f"no 'fract' line in {skl_path}")


class TestSpaceGroupPreserved:
    """A converted crystal keeps its space group: the skeleton stores the
    asymmetric unit plus a space-group token, not a flattened P1 cell."""

    def test_au_fcc_asymmetric_unit_expands(self, tmp_path,
                                            _sc_import):
        out = str(tmp_path / "au.skl")
        token = convert(_cif("au_fcc_225.cif"), out)
        # Fm-3m: a single authored atom is stored ...
        assert _count_asymmetric_atoms(out) == 1
        assert token.startswith("225")
        # ... and Imago regenerates the four-atom fcc cell on read.
        structure = _sc_import()
        structure.read_input_file(out)
        assert structure.num_atoms == 4

    def test_asymmetric_unit_smaller_than_full_cell(self, tmp_path,
                                                    _sc_import):
        # The silicon diamond cell holds more atoms than its asymmetric
        # unit, proving the converter did not flatten to P1.
        out = str(tmp_path / "si.skl")
        convert(_cif("si_fd3m_origin1.cif"), out)
        stored = _count_asymmetric_atoms(out)
        structure = _sc_import()
        structure.read_input_file(out)
        assert stored < structure.num_atoms


class TestOriginDiscrimination:
    """The two origin settings of one phase must resolve to different
    spaceDB tokens -- the verification picks the origin whose operations
    reproduce the CIF, not a fixed default."""

    def test_two_origins_resolve_differently(self, tmp_path):
        _, token1 = cif_to_skeleton(_cif("si_fd3m_origin1.cif"))
        _, token2 = cif_to_skeleton(_cif("si_fd3m_origin2.cif"))
        assert token1.startswith("227")
        assert token2.startswith("227")
        assert token1 != token2

    def test_both_origins_verify_round_trip(self, tmp_path, _sc_import):
        # Each origin's skeleton, read back and expanded, reproduces the
        # same number of atoms (the verification that gated the write).
        counts = []
        for name in ("si_fd3m_origin1.cif", "si_fd3m_origin2.cif"):
            out = str(tmp_path / (name + ".skl"))
            convert(_cif(name), out)
            structure = _sc_import()
            structure.read_input_file(out)
            counts.append(structure.num_atoms)
        assert counts[0] == counts[1]      # same phase, same atom count


class TestRefusalsAndOverride:
    """Unsupported inputs fail loudly; a forced setting is still
    verified."""

    def test_partial_occupancy_refused(self, tmp_path):
        cif = tmp_path / "disordered.cif"
        cif.write_text(
            "data_test\n"
            "_cell_length_a 5.0\n_cell_length_b 5.0\n"
            "_cell_length_c 5.0\n"
            "_cell_angle_alpha 90\n_cell_angle_beta 90\n"
            "_cell_angle_gamma 90\n"
            "_space_group_IT_number 1\n"
            "loop_\n"
            "_atom_site_label\n_atom_site_fract_x\n"
            "_atom_site_fract_y\n_atom_site_fract_z\n"
            "_atom_site_occupancy\n"
            "Si1 0.0 0.0 0.0 0.5\n")
        with pytest.raises(CifConversionError, match="occupancy"):
            cif_to_skeleton(str(cif))

    def test_wrong_forced_setting_errors(self, tmp_path):
        # Forcing P1 on a structure that needs Fm-3m cannot reproduce the
        # four-atom expansion, so the verified override fails loudly.
        with pytest.raises(CifConversionError, match="reproduce"):
            cif_to_skeleton(_cif("au_fcc_225.cif"), space_override="1_a")

    def test_forced_setting_accepts_correct_token(self, tmp_path):
        # Forcing the token auto-resolution would have chosen succeeds
        # (and still verifies) -- without hard-coding which letter it is.
        _, auto = cif_to_skeleton(_cif("si_fd3m_origin1.cif"))
        _, forced = cif_to_skeleton(
            _cif("si_fd3m_origin1.cif"), space_override=auto)
        assert forced == auto
