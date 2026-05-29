"""test_guidance_db.py -- Tests for the historical-guidance
dataspace library (DESIGN 7 / PSEUDOCODE 15), foundation layer.

This file covers the increment-1 surface of ``guidance_db.py``:
the element-group table loader, the space-group-number to
lattice-family mapping, and ``compute_signature``.  The reader,
emitter, and predictor are covered as those increments land.

conftest.py puts ``src/scripts`` on the path, so the module
imports directly.  The pure tests use a lightweight fake
structure (no ``IMAGO_DATA`` needed); the integration tests use
the real StructureControl fixtures and auto-skip when
``IMAGO_DATA`` is unset.
"""

import os

import pytest

from guidance_db import (
    CANONICAL_GROUP_ORDER,
    CANONICAL_LATTICE_ORDER,
    bravais_family_of,
    compute_signature,
    load_elemental_groups,
    _crystal_system_of,
)


# Path to the checked-in element-group table (src/data/).
GROUPS_PATH = os.path.join(
    os.path.dirname(__file__), "..", "data",
    "elemental_groups.toml")


class FakeStructure:
    """Minimal stand-in carrying just the attributes
    ``compute_signature`` duck-types: ``num_atoms``,
    ``atom_element_name`` (1-indexed, lower-cased, index 0
    unused -- mirroring StructureControl), and
    ``space_group_num``.
    """

    def __init__(self, symbols, space_group_num=0):
        self.num_atoms = len(symbols)
        self.atom_element_name = (
            [None] + [s.lower() for s in symbols])
        self.space_group_num = space_group_num


def _write_groups(tmp_path, overrides=None, schema_version=1,
                  omit=None, drop_groups_table=False):
    """Write a synthetic elemental_groups.toml.  Starts from all
    13 groups empty, applies ``overrides`` (group -> symbol
    list), optionally omits one group key, and returns the path.
    """

    groups = {group: [] for group in CANONICAL_GROUP_ORDER}
    if overrides:
        groups.update(overrides)
    if omit:
        groups.pop(omit, None)
    lines = ["schema_version = " + str(schema_version), ""]
    if not drop_groups_table:
        lines.append("[groups]")
        for group, symbols in groups.items():
            rendered = ", ".join('"' + s + '"' for s in symbols)
            lines.append(group + " = [" + rendered + "]")
    path = tmp_path / "elemental_groups.toml"
    path.write_text("\n".join(lines) + "\n")
    return str(path)


# ============================================================
#  load_elemental_groups (PSEUDOCODE 15.2)
# ============================================================

class TestLoadElementalGroups:
    def test_loads_real_checked_in_table(self):
        table = load_elemental_groups(GROUPS_PATH)
        # The day-1 table classifies the whole periodic table;
        # spot-check one element from several groups.
        assert table["si"] == "group_iv"
        assert table["o"] == "chalcogen"
        assert table["fe"] == "transition_metal"
        assert table["h"] == "hydrogen"
        assert table["la"] == "lanthanide"
        assert table["u"] == "actinide"
        assert table["cs"] == "alkali"
        assert table["xe"] == "noble_gas"
        assert table["b"] == "group_iii"
        assert table["as"] == "pnictogen"

    def test_metalloid_ships_empty(self):
        # v1 is column-based: nothing maps to metalloid yet
        # (DESIGN 7.4 / 7.10).
        table = load_elemental_groups(GROUPS_PATH)
        assert "metalloid" not in set(table.values())

    def test_lookup_is_case_insensitive(self):
        # The file stores "Fe"; lookups are keyed lower-cased so
        # StructureControl's lower-cased symbols match.
        table = load_elemental_groups(GROUPS_PATH)
        assert "fe" in table
        assert "Fe" not in table

    def test_rejects_double_assignment(self, tmp_path):
        # Si placed in two groups must fail loudly, not let the
        # last assignment silently win.
        path = _write_groups(tmp_path, overrides={
            "group_iv": ["Si"], "group_iii": ["Si"]})
        with pytest.raises(ValueError, match="assigned to"):
            load_elemental_groups(path)

    def test_rejects_wrong_schema_version(self, tmp_path):
        path = _write_groups(tmp_path, schema_version=2)
        with pytest.raises(ValueError, match="schema_version"):
            load_elemental_groups(path)

    def test_rejects_missing_group(self, tmp_path):
        path = _write_groups(tmp_path, omit="noble_gas")
        with pytest.raises(ValueError, match="missing group"):
            load_elemental_groups(path)

    def test_rejects_missing_groups_table(self, tmp_path):
        path = _write_groups(tmp_path, drop_groups_table=True)
        with pytest.raises(ValueError, match=r"\[groups\]"):
            load_elemental_groups(path)


# ============================================================
#  Space-group-number -> lattice family (PSEUDOCODE 15.2)
# ============================================================

class TestCrystalSystemAndFamily:
    @pytest.mark.parametrize("sg, system, family", [
        (1, "triclinic", "tri"),
        (2, "triclinic", "tri"),
        (5, "monoclinic", "mono"),
        (60, "orthorhombic", "ortho"),
        (136, "tetragonal", "tet"),
        (160, "trigonal", "hex"),     # trigonal lumps into hex
        (186, "hexagonal", "hex"),
        (225, "cubic", "cubic"),
        (230, "cubic", "cubic"),
    ])
    def test_ranges_and_family(self, sg, system, family):
        assert _crystal_system_of(sg) == system
        assert bravais_family_of(FakeStructure(["Si"], sg)) == family

    def test_sg_zero_raises(self):
        # 0 = no space group recorded -> cannot be signed.
        with pytest.raises(ValueError, match="1..230"):
            bravais_family_of(FakeStructure(["Si"], 0))

    def test_sg_out_of_range_raises(self):
        with pytest.raises(ValueError, match="1..230"):
            bravais_family_of(FakeStructure(["Si"], 231))


# ============================================================
#  compute_signature (PSEUDOCODE 15.2)
# ============================================================

class TestComputeSignature:
    def setup_method(self):
        self.table = load_elemental_groups(GROUPS_PATH)

    def _comp(self, signature, group):
        # Read one group's atom fraction out of the vector.
        return signature.composition_vector[
            CANONICAL_GROUP_ORDER.index(group)]

    def test_composition_tio2(self):
        # Rutile TiO2: Ti (transition_metal) + 2 O (chalcogen).
        struct = FakeStructure(["Ti", "O", "O"], 136)
        sig = compute_signature(struct, "crystalline", self.table)
        assert self._comp(sig, "transition_metal") == pytest.approx(1/3)
        assert self._comp(sig, "chalcogen") == pytest.approx(2/3)
        assert self._comp(sig, "alkali") == 0.0

    def test_composition_sums_to_one(self):
        struct = FakeStructure(["Ti", "O", "O"], 136)
        sig = compute_signature(struct, "crystalline", self.table)
        assert sum(sig.composition_vector) == pytest.approx(1.0)

    def test_crystalline_lattice_onehot(self):
        struct = FakeStructure(["Ti", "O", "O"], 136)
        sig = compute_signature(struct, "crystalline", self.table)
        assert sig.lattice_family == "tet"
        expected = tuple(
            1.0 if name == "tet" else 0.0
            for name in CANONICAL_LATTICE_ORDER)
        assert sig.lattice_onehot == expected
        assert sum(sig.lattice_onehot) == 1.0

    def test_noncrystalline_has_empty_lattice(self):
        # Even with a space_group_num set, a non-crystalline
        # system gets the empty family + zero one-hot.
        struct = FakeStructure(["Si", "O", "O"], 136)
        sig = compute_signature(struct, "amorphous", self.table)
        assert sig.lattice_family == ""
        assert sig.lattice_onehot == tuple(
            0.0 for _ in CANONICAL_LATTICE_ORDER)
        # Composition is still computed for non-crystalline.
        assert sum(sig.composition_vector) == pytest.approx(1.0)

    def test_unknown_element_raises(self):
        struct = FakeStructure(["Xx"], 225)
        with pytest.raises(ValueError, match="not in"):
            compute_signature(struct, "crystalline", self.table,
                              label="mystery.skl")

    def test_unknown_system_type_raises(self):
        struct = FakeStructure(["Si"], 225)
        with pytest.raises(ValueError, match="system_type"):
            compute_signature(struct, "liquid_crystal", self.table)

    def test_crystalline_without_space_group_raises(self):
        struct = FakeStructure(["Si"], 0)
        with pytest.raises(ValueError, match="1..230"):
            compute_signature(struct, "crystalline", self.table)

    def test_none_atom_name_raises(self):
        struct = FakeStructure(["Si"], 225)
        struct.atom_element_name[1] = None
        with pytest.raises(ValueError, match="no element name"):
            compute_signature(struct, "crystalline", self.table)


# ============================================================
#  Integration against real StructureControl fixtures
#  (auto-skipped when IMAGO_DATA is unset; conftest.py)
# ============================================================

class TestRealStructureIntegration:
    def test_si_diamond_signature(self, sc_si_diamond):
        # Si diamond, space group 227 -> cubic; pure group_iv.
        table = load_elemental_groups(GROUPS_PATH)
        sig = compute_signature(
            sc_si_diamond, "crystalline", table, "Si")
        assert sig.composition_vector[
            CANONICAL_GROUP_ORDER.index("group_iv")] == pytest.approx(1.0)
        assert sig.lattice_family == "cubic"

    def test_beo_hex_signature(self, sc_beo_hex):
        # BeO wurtzite, space group 186 -> hex; Be + O.
        table = load_elemental_groups(GROUPS_PATH)
        sig = compute_signature(
            sc_beo_hex, "crystalline", table, "BeO")
        assert sig.composition_vector[
            CANONICAL_GROUP_ORDER.index("alkali_earth")] == pytest.approx(0.5)
        assert sig.composition_vector[
            CANONICAL_GROUP_ORDER.index("chalcogen")] == pytest.approx(0.5)
        assert sig.lattice_family == "hex"
