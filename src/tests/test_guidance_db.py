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
import shutil

import pytest

from guidance_db import (
    CANONICAL_GROUP_ORDER,
    CANONICAL_LATTICE_ORDER,
    VALID_SYSTEM_TYPES,
    Context,
    Dataspace,
    GuidanceEntry,
    Measured,
    Provenance,
    Signature,
    bravais_family_of,
    compute_signature,
    format_entry,
    functional_family,
    k_neighbors,
    knn_weights,
    load,
    load_elemental_groups,
    load_entry,
    predict,
    save_entry,
    select_submodel,
    short_sha,
    slug_for,
    stage1,
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


# ============================================================
#  Reader: load() and load_entry() (PSEUDOCODE 15.3)
# ============================================================
#
# Default building blocks for a valid flight entry.  Each rule
# test copies one block and tweaks a single field, so a failure
# isolates exactly the rule under test.

DEF_TOP = {
    "schema_version": 1, "entry_id": "crystalline-aaa111",
    "generated_at": "2026-05-29T00:00:00Z", "source": "flight"}
DEF_SIG = {"system_type": "crystalline", "lattice_family": "tet"}
DEF_COMP = {group: 0.0 for group in CANONICAL_GROUP_ORDER}
DEF_COMP["transition_metal"] = 1.0 / 3.0
DEF_COMP["chalcogen"] = 2.0 / 3.0
DEF_MEAS = {
    "gap_ev": 3.0, "gap_kind": "direct", "spin_polarization": 0.0,
    "total_magnetization": 0.0, "kpoint_density": 50.0,
    "dos_at_fermi": 0.0}
DEF_CTX = {
    "basis": "fb", "functional": "gga-pbe", "scf_threshold": 1.0e-6,
    "cell_atom_count": 6, "cell_volume_per_formula_unit": 100.0}
DEF_VER = {
    "grid_values": [25.0, 50.0, 100.0],
    "grid_energies": [-1.0, -1.5, -1.51], "converged_at": 50.0,
    "metric": "total_energy", "metric_threshold": 1.0e-4,
    "predictor_confidence": 0.8, "predictor_neighbor_ids": ["a", "b"]}
DEF_PROV = {
    "flight_id": "seed_2026", "source_structure": "COD-1",
    "imago_commit": "abc1234", "curator": "harvest.py"}


def _val(value):
    """Render one Python value as a TOML scalar / array."""

    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return '"' + value + '"'
    if isinstance(value, list):
        return "[" + ", ".join(_val(item) for item in value) + "]"
    return repr(value)            # int -> "6"; float -> "1e-06"


def _entry_text(top=None, signature=None, composition=None,
                measured=None, context=None, verification=DEF_VER,
                provenance=None):
    """Build entry-file TOML text from the default blocks, with
    any block replaced wholesale.  ``verification=None`` omits
    the verification block entirely.
    """

    top = DEF_TOP if top is None else top
    signature = DEF_SIG if signature is None else signature
    composition = DEF_COMP if composition is None else composition
    measured = DEF_MEAS if measured is None else measured
    context = DEF_CTX if context is None else context
    provenance = DEF_PROV if provenance is None else provenance

    lines = [k + " = " + _val(v) for k, v in top.items()]
    lines += ["", "[entry.signature]"]
    lines += [k + " = " + _val(v) for k, v in signature.items()]
    lines += ["", "[entry.signature.composition_vector]"]
    lines += [g + " = " + _val(composition[g]) for g in composition]
    lines += ["", "[entry.measured]"]
    lines += [k + " = " + _val(v) for k, v in measured.items()]
    lines += ["", "[entry.context]"]
    lines += [k + " = " + _val(v) for k, v in context.items()]
    if verification is not None:
        lines += ["", "[entry.verification]"]
        lines += [k + " = " + _val(v)
                  for k, v in verification.items()]
    lines += ["", "[entry.provenance]"]
    lines += [k + " = " + _val(v) for k, v in provenance.items()]
    return "\n".join(lines) + "\n"


def _write_dataspace(tmp_path, entry_text=None,
                     system_type="crystalline", entry_name="e1.toml",
                     marker="1", entries=None):
    """Lay out a minimal dataspace under ``tmp_path``: the marker
    file, a copy of the real elemental_groups.toml, and one or
    more entry files.  ``entries`` is a list of
    ``(system_type, name, text)`` triples; if omitted a single
    entry is written from ``entry_text`` (default valid).
    """

    (tmp_path / "SCHEMA_VERSION").write_text(marker + "\n")
    shutil.copy(GROUPS_PATH,
                str(tmp_path / "elemental_groups.toml"))
    if entries is None:
        entries = [(system_type, entry_name,
                    entry_text if entry_text is not None
                    else _entry_text())]
    for system_type_dir, name, text in entries:
        directory = tmp_path / "entries" / system_type_dir
        directory.mkdir(parents=True, exist_ok=True)
        (directory / name).write_text(text)
    return str(tmp_path)


def _load_entry_text(tmp_path, text, system_type="crystalline"):
    """Write one entry file and validate it directly via
    load_entry (skipping the marker/groups scaffolding)."""

    path = tmp_path / "entry.toml"
    path.write_text(text)
    return load_entry(str(path), system_type, {})


class TestLoadDataspace:
    def test_loads_valid_flight_entry(self, tmp_path):
        space = load(_write_dataspace(tmp_path))
        assert space.schema_version == 1
        entries = space.entries_by_system_type["crystalline"]
        assert len(entries) == 1
        entry = entries[0]
        assert entry.source == "flight"
        assert entry.signature.lattice_family == "tet"
        assert entry.measured.kpoint_density == 50.0
        assert entry.verification is not None
        assert entry.verification.converged_at == 50.0
        assert entry.provenance.flight_id == "seed_2026"
        # Untouched partitions are empty, not missing.
        assert space.entries_by_system_type["amorphous"] == []

    def test_empty_dataspace_is_valid(self, tmp_path):
        (tmp_path / "SCHEMA_VERSION").write_text("1\n")
        shutil.copy(GROUPS_PATH,
                    str(tmp_path / "elemental_groups.toml"))
        space = load(str(tmp_path))
        assert all(entries == [] for entries
                   in space.entries_by_system_type.values())

    def test_manual_entry_without_verification(self, tmp_path):
        text = _entry_text(
            top={**DEF_TOP, "source": "manual"}, verification=None)
        space = load(_write_dataspace(tmp_path, text))
        entry = space.entries_by_system_type["crystalline"][0]
        assert entry.source == "manual"
        assert entry.verification is None

    def test_bad_marker_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="marker"):
            load(_write_dataspace(tmp_path, marker="2"))

    def test_duplicate_entry_id_across_files(self, tmp_path):
        text = _entry_text()
        space = _write_dataspace(tmp_path, entries=[
            ("crystalline", "a.toml", text),
            ("crystalline", "b.toml", text)])
        with pytest.raises(ValueError, match="duplicate entry_id"):
            load(space)

    def test_system_type_directory_mismatch(self, tmp_path):
        # Entry declares crystalline but lives under amorphous/.
        text = _entry_text()
        space = _write_dataspace(tmp_path, entries=[
            ("amorphous", "a.toml", text)])
        with pytest.raises(ValueError,
                           match="under entries/amorphous"):
            load(space)


class TestEntryValidationRules:
    def test_missing_top_level_source(self, tmp_path):
        top = {k: v for k, v in DEF_TOP.items() if k != "source"}
        with pytest.raises(ValueError,
                           match="missing top-level: source"):
            _load_entry_text(tmp_path, _entry_text(top=top))

    def test_entry_schema_version_mismatch(self, tmp_path):
        with pytest.raises(ValueError, match="schema_version"):
            _load_entry_text(tmp_path, _entry_text(
                top={**DEF_TOP, "schema_version": 2}))

    def test_invalid_source(self, tmp_path):
        with pytest.raises(ValueError, match="flight"):
            _load_entry_text(tmp_path, _entry_text(
                top={**DEF_TOP, "source": "bogus"}))

    def test_invalid_system_type(self, tmp_path):
        with pytest.raises(ValueError, match="invalid system_type"):
            _load_entry_text(tmp_path, _entry_text(
                signature={**DEF_SIG, "system_type": "bogus"}))

    def test_composition_keys_wrong(self, tmp_path):
        comp = dict(DEF_COMP)
        del comp["hydrogen"]
        with pytest.raises(ValueError,
                           match="composition_vector keys"):
            _load_entry_text(tmp_path,
                             _entry_text(composition=comp))

    def test_composition_sum_not_one(self, tmp_path):
        comp = {group: 0.0 for group in CANONICAL_GROUP_ORDER}
        comp["hydrogen"] = 0.5
        with pytest.raises(ValueError, match="sums to"):
            _load_entry_text(tmp_path,
                             _entry_text(composition=comp))

    def test_composition_out_of_range(self, tmp_path):
        comp = {group: 0.0 for group in CANONICAL_GROUP_ORDER}
        comp["hydrogen"] = 1.5
        with pytest.raises(ValueError, match=r"out of \[0,1\]"):
            _load_entry_text(tmp_path,
                             _entry_text(composition=comp))

    def test_crystalline_invalid_lattice_family(self, tmp_path):
        with pytest.raises(ValueError,
                           match="lattice_family must be one"):
            _load_entry_text(tmp_path, _entry_text(
                signature={**DEF_SIG, "lattice_family": "bogus"}))

    def test_crystalline_missing_lattice_family(self, tmp_path):
        with pytest.raises(ValueError,
                           match="lattice_family must be one"):
            _load_entry_text(tmp_path, _entry_text(
                signature={"system_type": "crystalline"}))

    def test_noncrystalline_with_lattice_family(self, tmp_path):
        with pytest.raises(ValueError,
                           match="must not set lattice_family"):
            _load_entry_text(tmp_path, _entry_text(
                signature={"system_type": "amorphous",
                           "lattice_family": "tet"}),
                system_type="amorphous")

    def test_gap_negative(self, tmp_path):
        with pytest.raises(ValueError, match="gap_ev < 0"):
            _load_entry_text(tmp_path, _entry_text(
                measured={**DEF_MEAS, "gap_ev": -1.0}))

    def test_invalid_gap_kind(self, tmp_path):
        with pytest.raises(ValueError, match="invalid gap_kind"):
            _load_entry_text(tmp_path, _entry_text(
                measured={**DEF_MEAS, "gap_kind": "bogus"}))

    def test_gap_kind_metal_mismatch(self, tmp_path):
        # gap_ev=3.0 but gap_kind="none" claims a metal.
        with pytest.raises(ValueError, match="iff"):
            _load_entry_text(tmp_path, _entry_text(
                measured={**DEF_MEAS, "gap_kind": "none"}))

    def test_metal_is_consistent(self, tmp_path):
        # gap_ev=0.0 + gap_kind="none" is a valid metal.
        entry = _load_entry_text(tmp_path, _entry_text(
            measured={**DEF_MEAS, "gap_ev": 0.0,
                      "gap_kind": "none"}))
        assert entry.measured.gap_ev == 0.0

    def test_kpoint_density_nonpositive(self, tmp_path):
        with pytest.raises(ValueError,
                           match="kpoint_density must be > 0"):
            _load_entry_text(tmp_path, _entry_text(
                measured={**DEF_MEAS, "kpoint_density": 0.0}))

    def test_invalid_basis(self, tmp_path):
        with pytest.raises(ValueError, match="invalid basis"):
            _load_entry_text(tmp_path, _entry_text(
                context={**DEF_CTX, "basis": "xx"}))

    def test_empty_functional(self, tmp_path):
        with pytest.raises(ValueError,
                           match="functional must be non-empty"):
            _load_entry_text(tmp_path, _entry_text(
                context={**DEF_CTX, "functional": ""}))

    def test_cell_atom_count_nonpositive(self, tmp_path):
        with pytest.raises(ValueError,
                           match="cell_atom_count must be > 0"):
            _load_entry_text(tmp_path, _entry_text(
                context={**DEF_CTX, "cell_atom_count": 0}))

    def test_cell_volume_nonpositive(self, tmp_path):
        with pytest.raises(ValueError,
                           match="cell_volume_per_formula_unit"):
            _load_entry_text(tmp_path, _entry_text(
                context={**DEF_CTX,
                         "cell_volume_per_formula_unit": 0.0}))

    def test_flight_requires_verification(self, tmp_path):
        with pytest.raises(ValueError,
                           match=r"requires \[entry.verification\]"):
            _load_entry_text(tmp_path,
                             _entry_text(verification=None))

    def test_flight_needs_nonempty_flight_id(self, tmp_path):
        with pytest.raises(ValueError, match="non-empty flight_id"):
            _load_entry_text(tmp_path, _entry_text(
                provenance={**DEF_PROV, "flight_id": ""}))


class TestVerificationRule10:
    def test_grid_not_sorted(self, tmp_path):
        with pytest.raises(ValueError, match="not sorted"):
            _load_entry_text(tmp_path, _entry_text(
                verification={**DEF_VER,
                              "grid_values": [100.0, 50.0, 25.0]}))

    def test_converged_at_not_in_grid(self, tmp_path):
        with pytest.raises(ValueError,
                           match="converged_at not present"):
            _load_entry_text(tmp_path, _entry_text(
                verification={**DEF_VER, "converged_at": 999.0}))

    def test_converged_at_neq_kpoint_density(self, tmp_path):
        # 100.0 is in the grid but != measured.kpoint_density.
        with pytest.raises(ValueError,
                           match="converged_at != measured"):
            _load_entry_text(tmp_path, _entry_text(
                verification={**DEF_VER, "converged_at": 100.0}))

    def test_unknown_metric(self, tmp_path):
        with pytest.raises(ValueError, match="unknown metric"):
            _load_entry_text(tmp_path, _entry_text(
                verification={**DEF_VER, "metric": "free_energy"}))

    def test_confidence_out_of_range(self, tmp_path):
        with pytest.raises(ValueError,
                           match="predictor_confidence out of"):
            _load_entry_text(tmp_path, _entry_text(
                verification={**DEF_VER,
                              "predictor_confidence": 1.5}))

    def test_grid_energies_length_mismatch(self, tmp_path):
        with pytest.raises(ValueError, match="grid_energies length"):
            _load_entry_text(tmp_path, _entry_text(
                verification={**DEF_VER,
                              "grid_energies": [-1.0, -1.5]}))

    def test_grid_energies_optional(self, tmp_path):
        version = {k: v for k, v in DEF_VER.items()
                   if k != "grid_energies"}
        entry = _load_entry_text(tmp_path,
                                 _entry_text(verification=version))
        assert entry.verification.grid_energies is None


# ============================================================
#  Emitter: save_entry() / format_entry() (PSEUDOCODE 15.4)
# ============================================================

class TestEmitter:
    def test_round_trip_flight_entry(self, tmp_path):
        # load -> save -> load reproduces every field bit-exactly
        # (the %.16e format round-trips binary64).  entry_id
        # becomes the slug, which save_entry assigns.
        original = _load_entry_text(tmp_path, _entry_text())
        path = save_entry(original, str(tmp_path))
        reloaded = load_entry(path, "crystalline", {})
        assert reloaded.signature == original.signature
        assert reloaded.measured == original.measured
        assert reloaded.context == original.context
        assert reloaded.verification == original.verification
        assert reloaded.provenance == original.provenance
        assert reloaded.source == original.source
        assert reloaded.entry_id == slug_for(original)

    def test_round_trip_manual_without_verification(self, tmp_path):
        text = _entry_text(
            top={**DEF_TOP, "source": "manual"}, verification=None)
        original = _load_entry_text(tmp_path, text)
        path = save_entry(original, str(tmp_path))
        out = open(path).read()
        assert "[entry.verification]" not in out
        reloaded = load_entry(path, "crystalline", {})
        assert reloaded.verification is None
        assert reloaded.measured == original.measured

    def test_byte_deterministic(self, tmp_path):
        entry = _load_entry_text(tmp_path, _entry_text())
        assert (format_entry(entry, "crystalline-x")
                == format_entry(entry, "crystalline-x"))

    def test_dos_omitted_when_none(self, tmp_path):
        measured = {k: v for k, v in DEF_MEAS.items()
                    if k != "dos_at_fermi"}
        entry = _load_entry_text(tmp_path,
                                 _entry_text(measured=measured))
        assert entry.measured.dos_at_fermi is None
        assert "dos_at_fermi" not in format_entry(entry, "x")

    def test_float_format_and_block_alignment(self, tmp_path):
        text = format_entry(
            _load_entry_text(tmp_path, _entry_text()), "x")
        # %.16e floats; [entry.measured] aligns '=' to the widest
        # key, total_magnetization (19).
        assert ("kpoint_density".ljust(19)
                + " = 5.0000000000000000e+01") in text
        # [entry.context] aligns to cell_volume_per_formula_unit.
        assert ("basis".ljust(28) + ' = "fb"') in text
        # [entry.provenance] aligns to source_structure (16).
        assert ("flight_id".ljust(16) + ' = "seed_2026"') in text

    def test_float_arrays_multiline(self, tmp_path):
        text = format_entry(
            _load_entry_text(tmp_path, _entry_text()), "x")
        # Array opener unpadded; one element per line, trailing
        # comma; closing bracket on its own line.
        assert "grid_values = [\n" in text
        assert "    2.5000000000000000e+01,\n" in text
        assert "\n]\n" in text
        # The verification scalars align to predictor_neighbor_ids
        # (22), not to their own shorter max.
        assert ("converged_at".ljust(22)
                + " = 5.0000000000000000e+01") in text

    def test_save_entry_path_and_collision(self, tmp_path):
        entry = _load_entry_text(tmp_path, _entry_text())
        path = save_entry(entry, str(tmp_path))
        assert path.endswith(".toml")
        assert os.path.join("staging", "crystalline") in path
        assert os.path.exists(path)
        # A second save to the same slug is refused (rule 2 / 7.5).
        with pytest.raises(ValueError, match="already exists"):
            save_entry(entry, str(tmp_path))

    def test_slug_and_short_sha(self, tmp_path):
        entry = _load_entry_text(tmp_path, _entry_text())
        slug = slug_for(entry)
        assert slug.startswith("crystalline-")
        suffix = slug.split("-", 1)[1]
        assert len(suffix) == 6
        assert all(c in "0123456789abcdef" for c in suffix)
        # Deterministic, and distinct when a provenance field
        # differs.
        assert short_sha("f1", "s1", "g1") == short_sha(
            "f1", "s1", "g1")
        assert short_sha("f1", "s1", "g1") != short_sha(
            "f2", "s1", "g1")


# ============================================================
#  Predictor: predict() and its stages (PSEUDOCODE 15.5)
# ============================================================

def _comp(**weights):
    """Build a 13-d composition vector from group -> weight."""

    vector = [0.0] * len(CANONICAL_GROUP_ORDER)
    for group, weight in weights.items():
        vector[CANONICAL_GROUP_ORDER.index(group)] = weight
    return tuple(vector)


def _pentry(entry_id, *, comp=None, lattice="cubic", gap=1.0,
            spin=0.0, kpd=50.0, basis="fb", functional="gga-pbe",
            source="flight", system_type="crystalline",
            generated_at="2026-01-01"):
    """Construct an in-memory GuidanceEntry for predictor tests
    (no file I/O -- the predictor operates on loaded entries)."""

    vector = comp if comp is not None else _comp()
    if system_type == "crystalline":
        family = lattice
        onehot = tuple(1.0 if n == lattice else 0.0
                       for n in CANONICAL_LATTICE_ORDER)
    else:
        family = ""
        onehot = tuple(0.0 for _ in CANONICAL_LATTICE_ORDER)
    return GuidanceEntry(
        entry_id=entry_id, generated_at=generated_at, source=source,
        signature=Signature(system_type, vector, family, onehot),
        measured=Measured(gap, "none" if gap == 0.0 else "direct",
                          spin, 0.0, kpd, None),
        context=Context(basis, functional, 1.0e-6, 6, 100.0),
        verification=None,
        provenance=Provenance("f", "s", "c", "cur"))


def _space(entries):
    """Wrap entries in a Dataspace, partitioned by system_type."""

    by_type = {st: [] for st in VALID_SYSTEM_TYPES}
    for entry in entries:
        by_type[entry.signature.system_type].append(entry)
    return Dataspace(1, by_type, {})


class TestPredictorNonCrystalline:
    def test_returns_canonical_manual_entry(self):
        space = _space([_pentry(
            "amorphous-1", system_type="amorphous", source="manual",
            kpd=80.0)])
        query = Signature("amorphous", _comp(hydrogen=1.0), "",
                          (0.0,) * 6)
        result = predict(space, query, "fb", "gga-pbe")
        assert result.predicted_kpoint_density == 80.0
        assert result.confidence == 1.0
        assert not result.is_under_trained
        assert result.neighbor_entry_ids == ("amorphous-1",)

    def test_under_trained_without_manual_entry(self):
        # Only flight entries, no hand-seeded canonical one.
        space = _space([_pentry(
            "amorphous-1", system_type="amorphous", source="flight")])
        query = Signature("amorphous", _comp(hydrogen=1.0), "",
                          (0.0,) * 6)
        result = predict(space, query, "fb", "gga-pbe")
        assert result.is_under_trained
        assert result.confidence == 0.0


class TestSubmodelSelection:
    def test_functional_family(self):
        assert functional_family("gga-pbe") == "gga"
        assert functional_family("gga-pw91") == "gga"
        assert functional_family("lda") == "lda"

    def test_exact_submodel(self):
        pool = [_pentry("e%d" % i) for i in range(3)]   # fb,gga-pbe
        entries, under = select_submodel(pool, "fb", "gga-pbe")
        assert not under
        assert len(entries) == 3

    def test_family_fallback_picks_most_populous(self):
        # Query (mb,gga-pbe): only 1 exact, but 3 (fb,gga-pbe) in
        # the gga family -> fall back to those three.
        pool = ([_pentry("mb1", basis="mb")]
                + [_pentry("fb%d" % i, basis="fb") for i in range(3)]
                + [_pentry("lda%d" % i, functional="lda")
                   for i in range(2)])
        entries, under = select_submodel(pool, "mb", "gga-pbe")
        assert not under
        assert len(entries) == 3
        assert all(e.context.basis == "fb"
                   and e.context.functional == "gga-pbe"
                   for e in entries)

    def test_pool_fallback_when_family_thin(self):
        # Query (eb,scan): no exact, no scan family, but 4 total.
        pool = ([_pentry("a", basis="fb", functional="lda"),
                 _pentry("b", basis="mb", functional="lda"),
                 _pentry("c", basis="fb", functional="gga-pbe"),
                 _pentry("d", basis="mb", functional="gga-pbe")])
        entries, under = select_submodel(pool, "eb", "scan")
        assert not under
        assert len(entries) == 4

    def test_under_trained_when_pool_too_thin(self):
        pool = [_pentry("a"), _pentry("b")]      # 2 < k_min (3)
        entries, under = select_submodel(pool, "fb", "gga-pbe")
        assert under


class TestKnnWeights:
    def test_weights_sum_to_one_and_cap_at_k(self):
        entries = [_pentry("e%d" % i) for i in range(7)]
        # Distance increases with index; nearest is e0.
        distance = {e: float(i) for i, e in enumerate(entries)}
        pairs = knn_weights(entries, lambda e: distance[e])
        assert len(pairs) == k_neighbors           # capped at 5
        assert sum(w for _, w in pairs) == pytest.approx(1.0)
        # Nearest neighbour carries the most weight.
        weights = [w for _, w in pairs]
        assert weights[0] == max(weights)

    def test_exact_match_dominates(self):
        entries = [_pentry("near"), _pentry("far1"), _pentry("far2")]
        distance = {entries[0]: 0.0, entries[1]: 10.0,
                    entries[2]: 20.0}
        pairs = knn_weights(entries, lambda e: distance[e])
        nearest, weight = pairs[0]
        assert nearest.entry_id == "near"
        assert weight > 0.999          # 1/eps swamps the rest


class TestStages:
    def test_stage1_weighted_gap_and_full_confidence(self):
        # All neighbours share gap=2.0 -> predicted 2.0, zero
        # variance -> confidence 1.0.
        entries = [_pentry("e%d" % i, gap=2.0) for i in range(3)]
        pgap, pspin, conf, ids = stage1(
            Signature("crystalline", _comp(), "cubic",
                      tuple(1.0 if n == "cubic" else 0.0
                            for n in CANONICAL_LATTICE_ORDER)),
            entries)
        assert pgap == pytest.approx(2.0)
        assert conf == pytest.approx(1.0)
        assert len(ids) == 3

    def test_stage1_lattice_term_separates_polytypes(self):
        # Same composition, different lattice; query is cubic, so
        # the cubic entry (d1=0) dominates over the tet one.
        cubic = _pentry("cubic", lattice="cubic", gap=1.0)
        tet = _pentry("tet", lattice="tet", gap=5.0)
        query = Signature(
            "crystalline", _comp(), "cubic",
            tuple(1.0 if n == "cubic" else 0.0
                  for n in CANONICAL_LATTICE_ORDER))
        pgap, _, _, _ = stage1(query, [cubic, tet])
        assert pgap == pytest.approx(1.0, abs=1e-3)


class TestPredictEndToEnd:
    def test_confident_prediction_when_neighbors_agree(self):
        comp = _comp(transition_metal=1.0 / 3.0, chalcogen=2.0 / 3.0)
        space = _space([_pentry("e%d" % i, comp=comp, gap=2.0,
                                kpd=60.0) for i in range(4)])
        query = Signature(
            "crystalline", comp, "cubic",
            tuple(1.0 if n == "cubic" else 0.0
                  for n in CANONICAL_LATTICE_ORDER))
        result = predict(space, query, "fb", "gga-pbe")
        assert not result.is_under_trained
        assert result.predicted_kpoint_density == pytest.approx(60.0)
        assert result.predicted_gap == pytest.approx(2.0)
        assert result.confidence == pytest.approx(1.0)
        assert len(result.neighbor_entry_ids) > 0

    def test_confidence_drops_when_densities_disagree(self):
        # Same gaps (so stage 2 weights are equal) but spread-out
        # densities -> a confident gap but an uncertain density.
        comp = _comp(transition_metal=1.0 / 3.0, chalcogen=2.0 / 3.0)
        space = _space([
            _pentry("e0", comp=comp, gap=2.0, kpd=40.0),
            _pentry("e1", comp=comp, gap=2.0, kpd=50.0),
            _pentry("e2", comp=comp, gap=2.0, kpd=60.0),
            _pentry("e3", comp=comp, gap=2.0, kpd=70.0)])
        query = Signature(
            "crystalline", comp, "cubic",
            tuple(1.0 if n == "cubic" else 0.0
                  for n in CANONICAL_LATTICE_ORDER))
        result = predict(space, query, "fb", "gga-pbe")
        assert result.predicted_kpoint_density == pytest.approx(55.0)
        assert result.confidence < 1.0

    def test_under_trained_thin_crystalline_pool(self):
        space = _space([_pentry("a"), _pentry("b")])     # 2 < k_min
        query = Signature(
            "crystalline", _comp(), "cubic",
            tuple(1.0 if n == "cubic" else 0.0
                  for n in CANONICAL_LATTICE_ORDER))
        result = predict(space, query, "fb", "gga-pbe")
        assert result.is_under_trained
        assert result.confidence == 0.0
        assert result.neighbor_entry_ids == ()
