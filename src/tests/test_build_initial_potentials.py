"""test_build_initial_potentials.py -- Unit tests for the
augmented-potential-database producer (C48).

build_initial_potentials.py is the *producer* of the library /
producer / consumer split (DESIGN 5.4).  These tests cover the
parts that are pure and fully exercisable without a live Imago
toolchain.  C48.1 -- the focus here -- is the curation-manifest
reader ``load_manifest_v2`` and its nine DESIGN 5.7 validation
rules; later increments (C48.2 isolated refresh, C48.3 SCF
orchestration) add their own tests.

The reader is exercised exactly like the per-element-database
reader: hand-written TOML strings written to a temp manifest file,
no Imago runs.  conftest.py's ``SCRIPTS_DIR`` insertion lets us
import ``build_initial_potentials`` directly.
"""

import os
import tomllib
import types

import pytest

from curation_manifest import (
    load_manifest_v2,
    CurationManifest,
    ReferenceSolid,
    ReferenceEntry,
    ManifestFingerprint,
    load_structure_sources,
)
from build_initial_potentials import (
    element_path,
    is_isolated_default_for,
    build_isolated_entry,
    list_element_dirs,
    refresh_isolated_entries,
    save_databases,
    _parse_pot_file,
    _parse_coeff_file,
    make_producer_options,
    _thermsmear_for,
    make_imago_provenance,
    build_loen_units,
    harvest_fingerprints,
    curation_workspace_root,
    structure_cache_dir,
    materialize_structure,
    materialize_only,
    extract_potential,
    read_site_identity_map,
    assemble_entry_label,
    make_run_log_entry,
    make_nonconverged_log_entry,
    write_run_log,
    apply_manifest_defaults,
    build_initial_potentials,
)
import build_initial_potentials as bip
import initial_potential_db as ipdb
import mesh_climb
from kaleidoscope import CalcUnit, Flight, ReportEntry
from kaleidoscope import cluster_config
from kaleidoscope.builders.kpoint_convergence import PredictionRecord


pytestmark = pytest.mark.unit


# ============================================================
#  Manifest builders
# ============================================================

# A canonical valid single-solid manifest in the cod_id form.
# The database-wide preferred recipe (DESIGN 5.7): one declaration
#   per method, each the family's single preferred record.  A
#   [characterization] block is required (rule 2), so this leads
#   every valid manifest; tests that probe the recipe rules rebuild
#   it (e.g. via _VALID_COD_MANIFEST.replace) rather than appending a
#   second block, which TOML would reject as a duplicate table.
_CHAR_BLOCK = (
    "[characterization]\n"
    "  [[characterization.fingerprint]]\n"
    "  method = \"bispectrum\"\n"
    "  sub_spec = { twoj1 = 8, twoj2 = 8 }\n"
    "  [[characterization.fingerprint]]\n"
    "  method = \"reduce\"\n"
    "  sub_spec = { level = 2, thick = 0.5, cutoff = 5.0,"
    " tolerance = 0.05 }\n")


# A fully valid cod-sourced manifest: the required recipe followed by
#   one reference solid carrying a single default customization.  The
#   bulk of the rule tests build on this body, replacing or appending
#   to exercise one rule at a time.
_VALID_COD_MANIFEST = (
    "schema_version = 2\n\n"
    + _CHAR_BLOCK +
    "\n[[reference_solid]]\n"
    "reference_id = \"au_fcc\"\n"
    "system_type = \"crystalline\"\n"
    "basis = \"fb\"\n"
    "functional = \"wigner\"\n"
    "kpoint_integration = \"linear-tetrahedral\"\n"
    "cod_id = 9008463\n"
    "cod_revision = \"2023-04-12\"\n"
    "kpoint_spec = { density = 60.0, shift = [0.0, 0.0, 0.0] }\n"
    "scf_threshold = 1.0e-6\n\n"
    "  [[reference_solid.entry]]\n"
    "  element = \"Au\"\n"
    "  atom_site = 1\n"
    "  label = \"default_solid\"\n"
    "  default = true\n"
    "  description = \"Au in fcc bulk (Fm-3m).\"\n")


def _write(tmp_path, text, name="manifest.toml") -> str:
    """Write manifest ``text`` to ``tmp_path/name``; return path."""

    path = tmp_path / name
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    return str(path)


# ============================================================
#  Happy path
# ============================================================

class TestLoadHappyPath:
    """A fully valid manifest parses into the expected
    dataclass tree.
    """

    def test_cod_form_parses(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST)
        manifest = load_manifest_v2(path)
        assert isinstance(manifest, CurationManifest)
        assert manifest.schema_version == 2
        assert len(manifest.reference_solids) == 1
        solid = manifest.reference_solids[0]
        assert solid.reference_id == "au_fcc"
        assert solid.system_type == "crystalline"
        assert solid.basis == "fb"
        assert solid.functional == "wigner"
        assert solid.kpoint_integration == "linear-tetrahedral"
        assert solid.cod_id == 9008463
        assert solid.cod_revision == "2023-04-12"
        assert solid.structure_path is None
        assert solid.scf_threshold == pytest.approx(1e-6)
        assert solid.kpoint_spec["density"] == pytest.approx(60.0)
        assert len(solid.entries) == 1
        entry = solid.entries[0]
        assert entry.element == "Au"
        assert entry.atom_site == 1
        assert entry.label == "default_solid"
        assert entry.default is True
        assert entry.fingerprints == []

    def test_solid_less_manifest_is_valid(self, tmp_path):
        # A manifest with the required recipe but no reference solids
        # is valid: the per-solid/per-entry rules (3-9) vacuously
        # pass, and rule 7 ranges over the (empty) set of elements.
        # The [characterization] block is still required (rule 2).
        path = _write(tmp_path, "schema_version = 2\n\n" + _CHAR_BLOCK)
        manifest = load_manifest_v2(path)
        assert manifest.reference_solids == []
        assert len(manifest.characterization) == 2

    def test_structure_path_form_parses(self, tmp_path):
        # structure_path must resolve to a real file under the
        # manifest's directory.
        _write(tmp_path, "dummy structure bytes\n", name="x.skel")
        text = (
            "schema_version = 2\n\n"
            + _CHAR_BLOCK +
            "\n[[reference_solid]]\n"
            "reference_id = \"x_local\"\n"
            "system_type = \"crystalline\"\n"
            "basis = \"fb\"\n"
            "functional = \"wigner\"\n"
            "kpoint_integration = \"linear-tetrahedral\"\n"
            "structure_path = \"x.skel\"\n"
            "kpoint_spec = { density = 60.0 }\n"
            "scf_threshold = 1.0e-6\n\n"
            "  [[reference_solid.entry]]\n"
            "  element = \"Si\"\n"
            "  atom_site = 1\n"
            "  label = \"default_solid\"\n"
            "  default = true\n"
            "  description = \"local Si.\"\n")
        path = _write(tmp_path, text)
        manifest = load_manifest_v2(path)
        solid = manifest.reference_solids[0]
        assert solid.structure_path == "x.skel"
        assert solid.cod_id is None

    def test_fingerprint_declarations_parse(self, tmp_path):
        # Two per-entry fingerprint overrides on the same entry.
        # Per-entry declarations are RARE non-preferred alternates
        # (the preferred recipe lives in [characterization]), so
        # both parse with preferred False.
        text = _VALID_COD_MANIFEST + (
            "\n"
            "    [[reference_solid.entry.fingerprint]]\n"
            "    method = \"bispectrum\"\n"
            "    sub_spec = { twoj1 = 8, twoj2 = 8 }\n\n"
            "    [[reference_solid.entry.fingerprint]]\n"
            "    method = \"bispectrum\"\n"
            "    sub_spec = { twoj1 = 6, twoj2 = 4 }\n")
        path = _write(tmp_path, text)
        manifest = load_manifest_v2(path)
        fps = manifest.reference_solids[0].entries[0].fingerprints
        assert len(fps) == 2
        assert all(isinstance(f, ManifestFingerprint) for f in fps)
        assert fps[0].method == "bispectrum"
        assert fps[0].sub_spec == {"twoj1": 8, "twoj2": 8}
        # Per-entry declarations are never preferred (rule 10).
        assert fps[0].preferred is False
        assert fps[1].preferred is False

    def test_characterization_block_parses(self, tmp_path):
        # The database-wide [characterization] recipe parses into
        # manifest.characterization, one record per method, each
        # marked preferred by the reader.
        path = _write(tmp_path, _VALID_COD_MANIFEST)
        manifest = load_manifest_v2(path)
        char = manifest.characterization
        assert len(char) == 2
        assert {fp.method for fp in char} == {"bispectrum", "reduce"}
        assert all(fp.preferred for fp in char)

    def test_absent_characterization_raises(self, tmp_path):
        # A [characterization] block is required (rule 2): a manifest
        # without one is refused, so the build cannot silently
        # produce a database with no preferred descriptors.
        path = _write(
            tmp_path, _VALID_COD_MANIFEST.replace(_CHAR_BLOCK, ""))
        with pytest.raises(
                ValueError,
                match="manifest rule 2.*characterization"):
            load_manifest_v2(path)


# ============================================================
#  Missing file (require the manifest to exist)
# ============================================================

class TestMissingManifest:
    """An absent manifest is a hard error -- a missing manifest
    is not an empty curation set.
    """

    def test_missing_file_raises(self, tmp_path):
        path = str(tmp_path / "does_not_exist.toml")
        with pytest.raises(FileNotFoundError, match="not found"):
            load_manifest_v2(path)


# ============================================================
#  Validation rule firings (DESIGN 5.7)
# ============================================================

class TestRule1SchemaVersion:
    def test_wrong_version_raises(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            "schema_version = 2", "schema_version = 1"))
        with pytest.raises(ValueError, match="manifest rule 1"):
            load_manifest_v2(path)

    def test_missing_version_raises(self, tmp_path):
        body = _VALID_COD_MANIFEST.split("\n", 1)[1]  # drop line 1
        path = _write(tmp_path, body)
        with pytest.raises(ValueError, match="manifest rule 1"):
            load_manifest_v2(path)


class TestRule2RequiredSolidFields:
    def test_missing_reference_id_raises(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            'reference_id = "au_fcc"\n', ""))
        with pytest.raises(ValueError,
                           match="manifest rule 2.*reference_id"):
            load_manifest_v2(path)

    def test_missing_scf_threshold_raises(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            "scf_threshold = 1.0e-6\n", ""))
        with pytest.raises(
                ValueError,
                match="manifest rule 2.*scf_threshold"):
            load_manifest_v2(path)

    def test_missing_system_type_raises(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            'system_type = "crystalline"\n', ""))
        with pytest.raises(ValueError,
                           match="manifest rule 2.*system_type"):
            load_manifest_v2(path)

    def test_missing_basis_raises(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            'basis = "fb"\n', ""))
        with pytest.raises(ValueError,
                           match="manifest rule 2.*basis"):
            load_manifest_v2(path)

    def test_missing_kpoint_integration_raises(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            'kpoint_integration = "linear-tetrahedral"\n', ""))
        with pytest.raises(
                ValueError,
                match="manifest rule 2.*kpoint_integration"):
            load_manifest_v2(path)

    def test_invalid_system_type_raises(self, tmp_path):
        # A system_type outside the four-value domain is a hard
        # error -- the predictor switches its sub-model on it.
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            'system_type = "crystalline"\n',
            'system_type = "liquid"\n'))
        with pytest.raises(
                ValueError,
                match="manifest rule 2.*system_type.*not one"):
            load_manifest_v2(path)


class TestRule3OptionalEntryFields:
    """Rule 3 (DESIGN 5.7): entries are OPTIONAL customizations and
    every field is optional.  An absent field is filled by the
    harvest (or, for default, reads as False); none is a hard
    error."""

    def test_solid_with_no_entries_loads(self, tmp_path):
        # A reference solid may carry no customizations at all: the
        # harvest auto-discovers every environment on its own.
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            '\n  [[reference_solid.entry]]\n'
            '  element = "Au"\n'
            '  atom_site = 1\n'
            '  label = "default_solid"\n'
            '  default = true\n'
            '  description = "Au in fcc bulk (Fm-3m)."\n', ""))
        manifest = load_manifest_v2(path)
        assert manifest.reference_solids[0].entries == []

    def test_missing_default_reads_as_false(self, tmp_path):
        # An absent default is no longer an error; it reads False.
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            "  default = true\n", ""))
        manifest = load_manifest_v2(path)
        assert manifest.reference_solids[0].entries[0].default \
            is False

    def test_missing_atom_site_reads_as_none(self, tmp_path):
        # An absent atom_site is filled by the harvest; the parsed
        # customization holds atom_site is None.
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            "  atom_site = 1\n", ""))
        manifest = load_manifest_v2(path)
        assert manifest.reference_solids[0].entries[0].atom_site \
            is None

    def test_bare_entry_loads_all_optional(self, tmp_path):
        # A customization with NO fields at all is valid: every
        # field is optional, so it parses to all-None / False.
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            '  element = "Au"\n'
            '  atom_site = 1\n'
            '  label = "default_solid"\n'
            '  default = true\n'
            '  description = "Au in fcc bulk (Fm-3m)."\n', ""))
        manifest = load_manifest_v2(path)
        entry = manifest.reference_solids[0].entries[0]
        assert entry.element is None
        assert entry.atom_site is None
        assert entry.default is False
        assert entry.description is None
        assert entry.label is None


class TestRule4StructureSource:
    def test_both_sources_raises(self, tmp_path):
        # Add a structure_path alongside the cod_id.
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            'cod_id = 9008463\n',
            'cod_id = 9008463\nstructure_path = "x.skel"\n'))
        with pytest.raises(ValueError, match="manifest rule 4"):
            load_manifest_v2(path)

    def test_neither_source_raises(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            'cod_id = 9008463\ncod_revision = "2023-04-12"\n', ""))
        with pytest.raises(ValueError, match="manifest rule 4"):
            load_manifest_v2(path)

    def test_cod_id_without_revision_raises(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            'cod_revision = "2023-04-12"\n', ""))
        with pytest.raises(ValueError,
                           match="manifest rule 4.*cod_revision"):
            load_manifest_v2(path)

    def test_cod_id_not_positive_int_raises(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            "cod_id = 9008463", "cod_id = -3"))
        with pytest.raises(ValueError,
                           match="manifest rule 4.*positive"):
            load_manifest_v2(path)

    def test_structure_path_missing_file_raises(self, tmp_path):
        text = (
            "schema_version = 2\n\n"
            + _CHAR_BLOCK +
            "\n[[reference_solid]]\n"
            "reference_id = \"x_local\"\n"
            "system_type = \"crystalline\"\n"
            "basis = \"fb\"\n"
            "functional = \"wigner\"\n"
            "kpoint_integration = \"linear-tetrahedral\"\n"
            "structure_path = \"absent.skel\"\n"
            "kpoint_spec = { density = 60.0 }\n"
            "scf_threshold = 1.0e-6\n\n"
            "  [[reference_solid.entry]]\n"
            "  element = \"Si\"\n"
            "  atom_site = 1\n"
            "  label = \"default_solid\"\n"
            "  default = true\n"
            "  description = \"local Si.\"\n")
        path = _write(tmp_path, text)
        with pytest.raises(ValueError,
                           match="manifest rule 4.*missing file"):
            load_manifest_v2(path)


class TestRule5ReferenceIdUniqueness:
    def test_duplicate_reference_id_raises(self, tmp_path):
        # Two solids with the same reference_id but distinct
        # (element, label) so rule 6 does not fire first.
        text = _VALID_COD_MANIFEST + (
            "\n[[reference_solid]]\n"
            "reference_id = \"au_fcc\"\n"
            "system_type = \"crystalline\"\n"
            "basis = \"fb\"\n"
            "functional = \"wigner\"\n"
            "kpoint_integration = \"linear-tetrahedral\"\n"
            "cod_id = 1234567\n"
            "cod_revision = \"2023-01-01\"\n"
            "kpoint_spec = { density = 60.0 }\n"
            "scf_threshold = 1.0e-6\n\n"
            "  [[reference_solid.entry]]\n"
            "  element = \"Ag\"\n"
            "  atom_site = 1\n"
            "  label = \"default_solid\"\n"
            "  default = true\n"
            "  description = \"Ag bulk.\"\n")
        path = _write(tmp_path, text)
        with pytest.raises(ValueError, match="manifest rule 5"):
            load_manifest_v2(path)

    def test_non_label_safe_reference_id_raises(self, tmp_path):
        # An uppercase reference_id is rejected: it is embedded
        # verbatim in derived labels and typed into -pot (5.2.1).
        text = _VALID_COD_MANIFEST.replace(
            'reference_id = "au_fcc"', 'reference_id = "Au_FCC"')
        path = _write(tmp_path, text)
        with pytest.raises(ValueError,
                           match="manifest rule 5.*label-safe"):
            load_manifest_v2(path)


class TestRule6ElementLabelUniqueness:
    def test_duplicate_element_label_raises(self, tmp_path):
        # Second solid produces the same (Au, default_solid).
        text = _VALID_COD_MANIFEST + (
            "\n[[reference_solid]]\n"
            "reference_id = \"au_hcp\"\n"
            "system_type = \"crystalline\"\n"
            "basis = \"fb\"\n"
            "functional = \"wigner\"\n"
            "kpoint_integration = \"linear-tetrahedral\"\n"
            "cod_id = 1234567\n"
            "cod_revision = \"2023-01-01\"\n"
            "kpoint_spec = { density = 60.0 }\n"
            "scf_threshold = 1.0e-6\n\n"
            "  [[reference_solid.entry]]\n"
            "  element = \"Au\"\n"
            "  atom_site = 1\n"
            "  label = \"default_solid\"\n"
            "  default = false\n"
            "  description = \"Au hcp.\"\n")
        path = _write(tmp_path, text)
        with pytest.raises(ValueError, match="manifest rule 6"):
            load_manifest_v2(path)

    def test_optional_label_parses_as_none(self, tmp_path):
        # An entry that omits label is valid (5.2.1): the producer
        # derives the label at harvest, so the parsed entry holds
        # label is None.
        text = _VALID_COD_MANIFEST.replace(
            '  label = "default_solid"\n', "")
        path = _write(tmp_path, text)
        manifest = load_manifest_v2(path)
        entry = manifest.reference_solids[0].entries[0]
        assert entry.label is None

    def test_label_less_customizations_need_no_cross_check(
            self, tmp_path):
        # Derived labels are unique by construction (DESIGN 5.7):
        # each environment mints its own at harvest from the run
        # identity, so two label-less customizations carry no
        # cross-manifest collision and the manifest loads.
        text = _VALID_COD_MANIFEST.replace(
            '  label = "default_solid"\n', "") + (
            "\n  [[reference_solid.entry]]\n"
            "  element = \"Au\"\n"
            "  atom_site = 2\n"
            "  default = false\n"
            "  description = \"Au bulk, second site.\"\n")
        path = _write(tmp_path, text)
        manifest = load_manifest_v2(path)
        assert len(manifest.reference_solids[0].entries) == 2


class TestRule7DefaultPerElement:
    def test_zero_defaults_loads(self, tmp_path):
        # Zero default customizations is no longer an error (DESIGN
        # 5.7): the element takes its isolated baseline as the
        # default at harvest.  Load enforces only "at most one".
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            "  default = true\n", "  default = false\n"))
        manifest = load_manifest_v2(path)
        assert manifest.reference_solids[0].entries[0].default \
            is False

    def test_two_defaults_same_element_raises(self, tmp_path):
        # A second Au entry (distinct label) also marked default.
        text = _VALID_COD_MANIFEST + (
            "\n  [[reference_solid.entry]]\n"
            "  element = \"Au\"\n"
            "  atom_site = 2\n"
            "  label = \"surface\"\n"
            "  default = true\n"
            "  description = \"Au surface site.\"\n")
        path = _write(tmp_path, text)
        with pytest.raises(ValueError,
                           match="manifest rule 7.*Au has 2"):
            load_manifest_v2(path)


class TestRule8FingerprintUniqueness:
    def test_duplicate_method_subspec_raises(self, tmp_path):
        text = _VALID_COD_MANIFEST + (
            "\n"
            "    [[reference_solid.entry.fingerprint]]\n"
            "    method = \"bispectrum\"\n"
            "    sub_spec = { twoj1 = 8, twoj2 = 8 }\n\n"
            "    [[reference_solid.entry.fingerprint]]\n"
            "    method = \"bispectrum\"\n"
            "    sub_spec = { twoj1 = 8, twoj2 = 8 }\n")
        path = _write(tmp_path, text)
        with pytest.raises(ValueError, match="manifest rule 8"):
            load_manifest_v2(path)

    def test_canonical_duplicate_raises(self, tmp_path):
        # Reordered keys + int-vs-float spelling are the SAME
        # sub-spec under canonical equality, so this is a dup.
        text = _VALID_COD_MANIFEST + (
            "\n"
            "    [[reference_solid.entry.fingerprint]]\n"
            "    method = \"bispectrum\"\n"
            "    sub_spec = { twoj1 = 8, twoj2 = 8 }\n\n"
            "    [[reference_solid.entry.fingerprint]]\n"
            "    method = \"bispectrum\"\n"
            "    sub_spec = { twoj2 = 8.0, twoj1 = 8.0 }\n")
        path = _write(tmp_path, text)
        with pytest.raises(ValueError, match="manifest rule 8"):
            load_manifest_v2(path)

    def test_same_method_different_subspec_ok(self, tmp_path):
        # Same method, two distinct sub_specs on one entry -> rule 8
        # is satisfied (the database stores as many sub_specs per
        # family as the curator wants), so the file loads cleanly.
        text = _VALID_COD_MANIFEST + (
            "\n"
            "    [[reference_solid.entry.fingerprint]]\n"
            "    method = \"bispectrum\"\n"
            "    sub_spec = { twoj1 = 8, twoj2 = 8 }\n\n"
            "    [[reference_solid.entry.fingerprint]]\n"
            "    method = \"bispectrum\"\n"
            "    sub_spec = { twoj1 = 6, twoj2 = 4 }\n")
        path = _write(tmp_path, text)
        manifest = load_manifest_v2(path)
        fps = manifest.reference_solids[0].entries[0].fingerprints
        assert len(fps) == 2


class TestRule9MethodRegistered:
    """Rule 9 is enforced only when a matcher registry
    (known_methods) is supplied; with None it is skipped, which
    is how C48 loads manifests before the registry exists (C54).
    """

    _UNKNOWN = _VALID_COD_MANIFEST + (
        "\n"
        "    [[reference_solid.entry.fingerprint]]\n"
        "    method = \"nonsense\"\n"
        "    sub_spec = { twoj1 = 8, twoj2 = 8 }\n")

    def test_unknown_method_raises_with_registry(self, tmp_path):
        path = _write(tmp_path, self._UNKNOWN)
        with pytest.raises(ValueError, match="manifest rule 9"):
            load_manifest_v2(
                path, known_methods={"bispectrum", "reduce"})

    def test_unknown_method_skipped_without_registry(
            self, tmp_path):
        path = _write(tmp_path, self._UNKNOWN)
        manifest = load_manifest_v2(path)  # known_methods=None
        fps = manifest.reference_solids[0].entries[0].fingerprints
        assert fps[0].method == "nonsense"

    def test_unknown_characterization_method_raises(self, tmp_path):
        # Rule 9 covers the [characterization] block too: an unknown
        # preferred method is a hard error when a registry is given.
        # Swap the recipe's reduce record for an unknown method.
        text = _VALID_COD_MANIFEST.replace(
            'method = "reduce"', 'method = "nonsense"')
        path = _write(tmp_path, text)
        with pytest.raises(ValueError, match="manifest rule 9"):
            load_manifest_v2(
                path, known_methods={"bispectrum", "reduce"})

    def test_known_method_accepted_with_registry(self, tmp_path):
        text = _VALID_COD_MANIFEST + (
            "\n"
            "    [[reference_solid.entry.fingerprint]]\n"
            "    method = \"bispectrum\"\n"
            "    sub_spec = { twoj1 = 8, twoj2 = 8 }\n")
        path = _write(tmp_path, text)
        manifest = load_manifest_v2(
            path, known_methods={"bispectrum", "reduce"})
        assert len(
            manifest.reference_solids[0].entries[0].fingerprints
        ) == 1


# A second crystalline solid producing a different element (Ag),
# used to exercise the cross-manifest fingerprint rules.  The
# caller appends ``[[reference_solid.entry.fingerprint]]`` blocks.
_SECOND_AG_SOLID = (
    "\n[[reference_solid]]\n"
    "reference_id = \"ag_fcc\"\n"
    "system_type = \"crystalline\"\n"
    "basis = \"fb\"\n"
    "functional = \"wigner\"\n"
    "kpoint_integration = \"linear-tetrahedral\"\n"
    "cod_id = 9008464\n"
    "cod_revision = \"2023-04-12\"\n"
    "kpoint_spec = { density = 60.0, shift = [0.0, 0.0, 0.0] }\n"
    "scf_threshold = 1.0e-6\n"
    "\n"
    "  [[reference_solid.entry]]\n"
    "  element = \"Ag\"\n"
    "  atom_site = 1\n"
    "  label = \"default_solid\"\n"
    "  default = true\n"
    "  description = \"Ag fcc.\"\n")


class TestRule10ManifestPreferred:
    """Manifest rule 10 (DESIGN 5.7): the [characterization] block
    declares at most one fingerprint per method -- that single
    declaration is the family's database-wide preferred record -- and
    a per-entry declaration may not be marked preferred."""

    def test_characterization_one_per_method_ok(self, tmp_path):
        # A well-formed recipe (bispectrum + reduce, one each) loads,
        # and each characterization record is marked preferred.
        path = _write(tmp_path, _VALID_COD_MANIFEST)
        manifest = load_manifest_v2(path)
        assert len(manifest.characterization) == 2
        assert all(fp.preferred for fp in manifest.characterization)

    def test_method_declared_twice_in_characterization_raises(
            self, tmp_path):
        # The preferred recipe has exactly one home: naming a method
        # twice in [characterization] is a hard error.  Renaming the
        # recipe's reduce record to bispectrum collides the two.
        text = _VALID_COD_MANIFEST.replace(
            'method = "reduce"', 'method = "bispectrum"')
        path = _write(tmp_path, text)
        with pytest.raises(ValueError,
                           match="manifest rule 10.*twice"):
            load_manifest_v2(path)

    def test_per_entry_preferred_raises(self, tmp_path):
        # A per-entry fingerprint may not be preferred -- the
        # preferred record is fixed by [characterization].
        text = _VALID_COD_MANIFEST + (
            "\n"
            "    [[reference_solid.entry.fingerprint]]\n"
            "    method = \"bispectrum\"\n"
            "    sub_spec = { twoj1 = 8, twoj2 = 8 }\n"
            "    preferred = true\n")
        path = _write(tmp_path, text)
        with pytest.raises(
                ValueError,
                match="manifest rule 10.*may not be preferred"):
            load_manifest_v2(path)


class TestRule11ManifestPreferredSubspec:
    """Manifest rule 11 (DESIGN 5.7): the preferred sub_spec for a
    family is uniform across the whole database.  It holds
    structurally -- the preferred recipe is the single
    [characterization] declaration per method, so it cannot diverge
    between elements; there is no runtime check to fire."""

    def test_single_recipe_applies_across_elements(self, tmp_path):
        # Two solids contributing two elements (Au, Ag) share the one
        # database-wide [characterization] recipe.  There is no
        # per-element preferred record to diverge, so the manifest
        # loads and the recipe is a single bispectrum + reduce pair.
        text = (_VALID_COD_MANIFEST + _SECOND_AG_SOLID)
        path = _write(tmp_path, text)
        manifest = load_manifest_v2(path)
        assert len(manifest.reference_solids) == 2
        methods = {fp.method for fp in manifest.characterization}
        assert methods == {"bispectrum", "reduce"}


def test_harvest_fingerprints_recipe_and_override(tmp_path):
    """Every environment harvests the database-wide [characterization]
    recipe (preferred) plus any per-entry override (non-preferred);
    the preferred flag rides through onto each FingerprintRecord
    (DESIGN 5.7 / 5.6.5)."""
    dat_skl = tmp_path / "datSkl.map"
    dat_skl.write_text("DAT SKEL ELEM SPECIES TYPE\n"
                       "  2    1   Si       1    1\n"
                       "  1    2   O        1    1\n")
    result_toml = {"outputs": {"datSkl_map": str(dat_skl)}}

    # Both sub_specs use twoj2 = 4 (5 descriptor components) but
    #   differ in twoj1, so each has its own dispatched loen unit.
    recipe_spec = {"twoj1": 6, "twoj2": 4, "cutoff": 9.0}
    override_spec = {"twoj1": 4, "twoj2": 4, "cutoff": 9.0}
    rows = ("1 O  1 1 1   9 9 9 9 9   0\n"
            "2 Si 1 1 2   1 2 3 4 5   0\n")
    # Each bispectrum sub_spec has its own dispatched loen unit.
    _write_loen_descriptor(tmp_path, "au_fcc", recipe_spec, rows)
    flight = _write_loen_descriptor(
        tmp_path, "au_fcc", override_spec, rows)

    characterization = [ManifestFingerprint(
        method="bispectrum", sub_spec=recipe_spec, preferred=True)]
    spec = ReferenceEntry(
        element="Si", atom_site=1, label="t", default=True,
        description="d",
        fingerprints=[ManifestFingerprint(
            method="bispectrum", sub_spec=override_spec)])
    records = harvest_fingerprints(
        flight, _ref(), spec.atom_site, spec.fingerprints,
        result_toml, characterization)

    assert len(records) == 2
    # The [characterization] record comes first and is preferred; the
    #   per-entry override follows and is not.
    assert records[0].sub_spec == recipe_spec
    assert records[0].preferred is True
    assert records[1].sub_spec == override_spec
    assert records[1].preferred is False


# ============================================================
#  C48.2 helpers: legacy pot1/coeff1 fixtures
# ============================================================

def _write_pot(elem_dir, num_alphas, alpha_min, alpha_max,
               nuclear_z=79.0, nuclear_alpha=20.0,
               covalent_radius=1.0):
    """Write a legacy ``pot1`` file in the fixed eight-line
    layout (the same shape atomSCF and the C47 consumer emit).
    """

    with open(os.path.join(elem_dir, "pot1"), "w") as handle:
        handle.write("NUCLEAR_CHARGE__ALPHA\n")
        handle.write(f"{nuclear_z:f} {nuclear_alpha:f}\n")
        handle.write("COVALENT_RADIUS\n")
        handle.write(f"{covalent_radius:f}\n")
        handle.write("NUM_ALPHAS\n")
        handle.write(f"{num_alphas}\n")
        handle.write("ALPHAS\n")
        handle.write(f"{alpha_min:.6e} {alpha_max:.6e}\n")


def _write_coeff(elem_dir, coefficients, alphas, count=None):
    """Write a legacy ``coeff1`` file: a count line plus one
    five-column line per term (cols 3-5 are the ignored zeros).
    ``count`` defaults to the true term count; pass a different
    value to forge an inconsistent file.
    """

    if count is None:
        count = len(coefficients)
    with open(os.path.join(elem_dir, "coeff1"), "w") as handle:
        handle.write(f"   {count}\n")
        for coefficient, alpha in zip(coefficients, alphas):
            handle.write(
                f" {coefficient:.10E} {alpha:.10E}"
                f" 0.000000E+00 0.000000E+00 0.000000E+00\n")


def _make_element(pdb_root, elem, coefficients, alphas,
                  num_alphas=None, **pot_kw):
    """Create ``<pdb_root>/<elem>/`` with a pot1/coeff1 pair.

    ``num_alphas`` defaults to the term count; overriding it
    forges a pot/coeff disagreement for the consistency tests.
    Returns the element directory path.
    """

    elem_dir = os.path.join(pdb_root, elem)
    os.makedirs(elem_dir, exist_ok=True)
    if num_alphas is None:
        num_alphas = len(coefficients)
    _write_pot(elem_dir, num_alphas, min(alphas), max(alphas),
               **pot_kw)
    _write_coeff(elem_dir, coefficients, alphas)
    return elem_dir


def _empty_manifest() -> CurationManifest:
    """A manifest curating nothing (isolated baselines only)."""

    return CurationManifest(schema_version=2,
                            manifest_path="x.toml",
                            reference_solids=[])


def _manifest_curating_au() -> CurationManifest:
    """A manifest with one Au default-tagged curated entry."""

    return CurationManifest(
        schema_version=2, manifest_path="x.toml",
        reference_solids=[ReferenceSolid(
            reference_id="au_fcc",
            system_type="crystalline", basis="fb",
            functional="wigner",
            kpoint_integration="linear-tetrahedral",
            kpoint_spec={"density": 60.0},
            scf_threshold=1e-6,
            cod_id=9008463, cod_revision="2023-04-12",
            structure_path=None,
            entries=[ReferenceEntry(
                element="Au", atom_site=1,
                label="default_solid", default=True,
                description="Au bulk.")])])


# ============================================================
#  Legacy file parsers
# ============================================================

class TestParsePotFile:
    def test_parses_scalar_fields(self, tmp_path):
        elem_dir = _make_element(
            str(tmp_path), "au",
            [1.0, 2.0, 3.0], [0.15, 1.5, 1.0e8])
        pot = _parse_pot_file(os.path.join(elem_dir, "pot1"))
        assert pot.nuclear_z == pytest.approx(79.0)
        # Z is a real, not an int (Imago uses it as a real).
        assert isinstance(pot.nuclear_z, float)
        assert pot.nuclear_alpha == pytest.approx(20.0)
        assert pot.covalent_radius == pytest.approx(1.0)
        assert pot.num_gaussians == 3
        assert pot.alpha_min == pytest.approx(0.15)
        assert pot.alpha_max == pytest.approx(1.0e8)

    def test_bad_tag_raises(self, tmp_path):
        elem_dir = _make_element(
            str(tmp_path), "au", [1.0], [0.15])
        # Corrupt the first tag line.
        pot_path = os.path.join(elem_dir, "pot1")
        text = open(pot_path).read().replace(
            "NUCLEAR_CHARGE__ALPHA", "WRONG_TAG")
        with open(pot_path, "w") as handle:
            handle.write(text)
        with pytest.raises(ValueError, match="malformed pot file"):
            _parse_pot_file(pot_path)


class TestParseCoeffFile:
    def test_parses_columns_one_and_two(self, tmp_path):
        elem_dir = _make_element(
            str(tmp_path), "au",
            [1.0, -2.5, 3.0], [0.15, 1.5, 1.0e8])
        coeffs, alphas = _parse_coeff_file(
            os.path.join(elem_dir, "coeff1"))
        assert coeffs == pytest.approx([1.0, -2.5, 3.0])
        assert alphas == pytest.approx([0.15, 1.5, 1.0e8])

    def test_count_mismatch_raises(self, tmp_path):
        elem_dir = os.path.join(str(tmp_path), "au")
        os.makedirs(elem_dir)
        # Count line claims 5 but only two term lines follow.
        _write_coeff(elem_dir, [1.0, 2.0], [0.1, 0.2], count=5)
        with pytest.raises(ValueError, match="count line says 5"):
            _parse_coeff_file(os.path.join(elem_dir, "coeff1"))


# ============================================================
#  element_path / is_isolated_default_for
# ============================================================

class TestElementPath:
    def test_lowercases_element_dir(self):
        path = element_path("/root", "Au")
        assert path == os.path.join(
            "/root", "au", "s_gaussian_pot.toml")


class TestIsIsolatedDefaultFor:
    def test_true_when_manifest_empty(self):
        assert is_isolated_default_for(
            "au", _empty_manifest()) is True

    def test_false_when_manifest_curates_element(self):
        # The manifest's Au default_solid wins over the
        # baseline; comparison is case-insensitive ("Au" vs
        # the "au" directory name).
        assert is_isolated_default_for(
            "au", _manifest_curating_au()) is False

    def test_true_for_uncurated_element(self):
        # Si is not in the Au-only manifest, so its baseline is
        # that file's default.
        assert is_isolated_default_for(
            "si", _manifest_curating_au()) is True


# ============================================================
#  build_isolated_entry
# ============================================================

class TestBuildIsolatedEntry:
    def test_builds_from_pot_and_coeff(self, tmp_path):
        _make_element(str(tmp_path), "au",
                      [1.0, -2.5, 3.0], [0.15, 1.5, 1.0e8])
        entry = build_isolated_entry(
            str(tmp_path), "au", "deadbee", "2026-05-20T00:00:00Z",
            _empty_manifest())
        assert entry.label == "isolated"
        assert entry.default is True            # empty manifest
        assert "isolated Au atom" in entry.description
        assert entry.num_gaussians == 3
        assert entry.coefficients == pytest.approx(
            [1.0, -2.5, 3.0])
        assert entry.alphas == pytest.approx([0.15, 1.5, 1.0e8])
        assert entry.provenance["source"] == "atomSCF"
        assert entry.provenance["commit"] == "deadbee"
        assert entry.fingerprints == []

    def test_default_false_when_manifest_curates(self, tmp_path):
        _make_element(str(tmp_path), "au", [1.0], [0.15])
        entry = build_isolated_entry(
            str(tmp_path), "au", "c", "t",
            _manifest_curating_au())
        assert entry.default is False

    def test_pot_coeff_term_mismatch_raises(self, tmp_path):
        # pot declares 3 alphas, coeff carries only 2 terms.
        _make_element(str(tmp_path), "au",
                      [1.0, 2.0], [0.1, 0.2], num_alphas=3)
        with pytest.raises(ValueError, match="disagree on term"):
            build_isolated_entry(
                str(tmp_path), "au", "c", "t", _empty_manifest())


# ============================================================
#  list_element_dirs / refresh_isolated_entries / save
# ============================================================

class TestListElementDirs:
    def test_only_dirs_with_pot1_sorted(self, tmp_path):
        root = str(tmp_path)
        _make_element(root, "au", [1.0], [0.15])
        _make_element(root, "ag", [1.0], [0.15])
        # A sibling dir without a pot1 is skipped.
        os.makedirs(os.path.join(root, "cache"))
        assert list_element_dirs(root) == ["ag", "au"]


class TestRefreshIsolatedEntries:
    def test_creates_db_with_isolated_default(self, tmp_path):
        root = str(tmp_path)
        _make_element(root, "au",
                      [1.0, 2.0, 3.0], [0.15, 1.5, 1.0e8])
        dbs = refresh_isolated_entries(
            root, _empty_manifest(), "c", "t")
        assert set(dbs) == {"au"}
        db = dbs["au"]
        assert db.element_symbol == "Au"     # capitalized
        assert db.nuclear_z == pytest.approx(79.0)
        iso = ipdb.lookup(db, "isolated")
        assert iso.default is True

    def test_roundtrip_through_save_and_load(self, tmp_path):
        # Producing then loading the file must satisfy the
        # per-element database rules (6: isolated present;
        # 7: exactly one default).
        root = str(tmp_path)
        _make_element(root, "au",
                      [1.0, 2.0, 3.0], [0.15, 1.5, 1.0e8])
        dbs = refresh_isolated_entries(
            root, _empty_manifest(), "c", "t")
        save_databases(dbs, root)
        reloaded = ipdb.load(element_path(root, "au"))
        assert ipdb.baseline(reloaded).label == "isolated"
        assert ipdb.default_entry(reloaded).label == "isolated"

    def test_refresh_preserves_prior_harvested_entries(self, tmp_path):
        # INCREMENTAL (DESIGN 5.7): the refresh loads each element file
        # and refreshes only its isolated baseline, PRESERVING every
        # previously harvested entry; the harvest phase inserts-or-skips
        # this run's solids on top.  So a prior default_solid survives
        # the refresh -- the file grows, it is not rebuilt.
        root = str(tmp_path)
        _make_element(root, "au",
                      [1.0, 2.0, 3.0], [0.15, 1.5, 1.0e8])
        # Seed a valid v2 file: isolated (non-default) plus a
        # previously harvested default_solid (default).
        seed = ipdb.ElementDatabase(2, "Au", 79.0, 20.0, 1.0)
        seed.potentials.append(ipdb.PotentialEntry(
            "isolated", False, "old iso", 1, 0.15, 1.0e8,
            [9.0], [0.15],
            {"source": "atomSCF", "commit": "old",
             "generated_at": "old"}))
        seed.potentials.append(ipdb.PotentialEntry(
            "default_solid", True, "Au bulk", 1, 0.15, 1.0e8,
            [0.5], [0.15],
            {"source": "Imago", "commit": "old",
             "generated_at": "old", "reference_id": "au_fcc",
             "atom_site": 1, "kpoint_spec": "k",
             "scf_threshold": 1e-6, "scf_iterations": 9,
             "type_assignment": "symmetry"}))
        ipdb.save(seed, element_path(root, "au"))

        dbs = refresh_isolated_entries(
            root, _manifest_curating_au(), "new", "now")
        db = dbs["au"]
        # The prior default_solid is preserved; the isolated baseline
        #   is refreshed in place.
        labels = [e.label for e in db.potentials]
        assert "default_solid" in labels and "isolated" in labels
        # The refreshed baseline was rebuilt from current pot1/coeff1
        #   (3 terms, not the seed's 1) with the new commit.
        iso = ipdb.lookup(db, "isolated")
        assert iso.num_gaussians == 3
        assert iso.provenance["commit"] == "new"
        # The preserved default_solid keeps its old potential and
        #   provenance untouched.
        kept = ipdb.lookup(db, "default_solid")
        assert kept.coefficients == [0.5]
        assert kept.provenance["commit"] == "old"


# ============================================================
#  Producer pipeline -- incr 3c (the toolchain seam is mocked)
# ============================================================

def _ref(**overrides) -> ReferenceSolid:
    """A ReferenceSolid with sensible defaults for the helper
    tests; ``overrides`` replace individual fields."""

    base = dict(
        reference_id="au_fcc", system_type="crystalline",
        basis="fb", functional="wigner",
        kpoint_integration="linear-tetrahedral",
        kpoint_spec={"density": 60.0, "shift": [0.0, 0.0, 0.0]},
        scf_threshold=1.0e-6, cod_id=None, cod_revision=None,
        structure_path="au.skel",
        entries=[ReferenceEntry(
            element="Au", atom_site=1, label="default_solid",
            default=True, description="Au bulk.")])
    base.update(overrides)
    return ReferenceSolid(**base)


def _write_result(workspace, unit_id, calc, *, energy,
                  iterations=7, scfv="scfV.dat"):
    """Write one COMPLETED unit's status.toml and result.toml under
    the workspace (the kaleidoscope run-dir layout).  A completed
    run carries both: the completion gate reads status.toml first --
    the report entry ``collect`` builds carries the terminal status,
    checked before any result.toml is opened (DESIGN 6.2.10) -- then
    result.toml."""

    run_dir = os.path.join(workspace, "wingbeats", unit_id, *calc)
    os.makedirs(run_dir, exist_ok=True)
    with open(os.path.join(run_dir, "status.toml"), "w") as handle:
        handle.write('status = "done"\n')
    path = os.path.join(run_dir, "result.toml")
    with open(path, "w") as handle:
        handle.write(f"total_energy = {energy}\n")
        handle.write(f"scf_iterations = {iterations}\n")
        handle.write(f'outputs = {{ scfV = "{scfv}" }}\n')


def _write_failed(workspace, unit_id, calc):
    """Write one FAILED unit: a status.toml='failed' and NO
    result.toml, the shape left when a unit aborts at the
    makeinput/imago seam (DESIGN 6.2.10)."""

    run_dir = os.path.join(workspace, "wingbeats", unit_id, *calc)
    os.makedirs(run_dir, exist_ok=True)
    with open(os.path.join(run_dir, "status.toml"), "w") as handle:
        handle.write('status = "failed"\n')


# ---- pure helpers --------------------------------------------

def test_make_producer_options_emits_coded_tool_settings():
    options = make_producer_options(_ref(), "abc123")
    # Dest-keyed, coded -- each tool's own vocabulary (DESIGN
    #   6.2.10): functional->xccode, kpoint_integration->scfkpint,
    #   basis->scf_basis, scf_threshold->converg, shift->kpshift.
    assert options["scf_basis"] == "fb"
    assert options["xccode"] == 100              # wigner
    assert options["scfkpint"] == 1              # linear-tetrahedral
    assert options["converg"] == pytest.approx(1.0e-6)
    assert options["imago_commit"] == "abc123"
    assert options["kpshift"] == [0.0, 0.0, 0.0]
    # No physics-name keys leak in: the sub-model travels separately
    #   (DESIGN 6.2.8), and the swept k-density is added per grid
    #   point by the builder, not pinned in the fixed options.
    for absent in ("basis", "functional", "kpoint_integration",
                   "scf_threshold", "kpoint_shift", "kpd"):
        assert absent not in options
    # A bare integration token names no smearing width, so no
    #   thermal-smearing option leaks in: makeinput keeps its rc
    #   default.
    assert "thermsmear" not in options


def test_thermsmear_for_reads_width_or_none():
    # Bare tokens name no width -> None (makeinput keeps the rc
    #   therm_smear_main default).
    assert _thermsmear_for("gaussian") is None
    assert _thermsmear_for("linear-tetrahedral") is None
    # A smeared Gaussian names the sigma (eV) after the dash.
    assert _thermsmear_for("gaussian-0.1") == pytest.approx(0.1)
    assert _thermsmear_for("gaussian-0.25") == pytest.approx(0.25)


def test_thermsmear_for_rejects_malformed_width():
    # A dash with no number, or a non-numeric tail, is a manifest
    #   fault caught here rather than passed to makeinput blindly.
    with pytest.raises(ValueError, match="smearing width"):
        _thermsmear_for("gaussian-")
    with pytest.raises(ValueError, match="smearing width"):
        _thermsmear_for("gaussian-wide")


def test_make_producer_options_forwards_smearing_sigma():
    # A ``gaussian-<sigma>`` integration token threads the sigma to
    #   makeinput's ``thermsmear`` option (-> THERMAL_SMEARING_SIGMA)
    #   while the integration code stays Gaussian (0).
    options = make_producer_options(
        _ref(kpoint_integration="gaussian-0.1"), "abc123")
    assert options["scfkpint"] == 0
    assert options["thermsmear"] == pytest.approx(0.1)


def test_make_imago_provenance_satisfies_schema():
    prov = make_imago_provenance(
        "abc123", "2026-06-12T00:00:00Z", _ref(), 1, 7)
    assert prov["source"] == "Imago"
    assert prov["system_type"] == "crystalline"
    assert prov["scf_iterations"] == 7
    # It must pass the ipdb Imago-provenance validator unchanged.
    ipdb.require_provenance(prov, "x.toml", "default_solid")


_BISPEC_DECL = ManifestFingerprint(
    method="bispectrum",
    sub_spec={"twoj1": 4, "twoj2": 4, "cutoff": 9.0})


def test_build_loen_units_skips_python_side():
    # Python-side (reduce) declarations need no dispatched unit -- they
    #   are computed in process during the harvest -- and the default
    #   entry declares no fingerprints at all, so nothing is built.
    #   An empty [characterization] recipe adds nothing either.
    assert build_loen_units(_ref(), "au.skel", {}, []) == []


def test_build_loen_units_builds_one_run_per_subspec():
    """A bispectrum declaration yields one structure-only loen unit:
    kind "fingerprint" (so the convergence harvest skips it), the job
    overridden to ``loen``/``-scf no``, the LOEN block carried via
    ``-loeninput``, and a calc tag encoding the method and every
    sub_spec key.  Two entries that share the sub_spec collapse to a
    single unit -- one descriptor table covers every site."""
    ref = _ref(entries=[
        ReferenceEntry(element="Si", atom_site=1, label="a",
                       default=True, description="d",
                       fingerprints=[_BISPEC_DECL]),
        ReferenceEntry(element="O", atom_site=2, label="b",
                       default=False, description="d",
                       fingerprints=[_BISPEC_DECL])])
    options = {"xccode": 100, "imago_commit": "abc", "converg": 1.0e-6}
    units = build_loen_units(ref, "si.skel", options, [])

    assert len(units) == 1
    unit = units[0]
    assert unit.kind == "fingerprint"
    assert unit.id == "au_fcc"
    assert unit.structure == "si.skel"
    assert unit.options["job"] == "loen"
    assert unit.options["scf_basis"] == "no"
    assert unit.options["loeninput"] == [
        "1", "4", "4", "50", "9.0", "0.85"]
    # cutoff=9.0 formats slug-safe to "9"; keys are alphabetical.
    assert unit.calc == ("loen-bispectrum-cutoff_9-twoj1_4-twoj2_4",)
    # The base build options survive on the copy (not mutated away).
    assert unit.options["xccode"] == 100


def test_build_loen_units_covers_characterization_recipe():
    """The database-wide ``[characterization]`` recipe must build its
    own loen unit even when no entry declares a per-entry fingerprint.

    This is the seed-manifest case (Si defaults): bispectrum lives in
    ``[characterization]`` and the entries carry no overrides.  The
    harvest reads that recipe (:func:`harvest_fingerprints`), so a build
    that only looked at ``entry.fingerprints`` would dispatch no loen run
    and the harvest would fail on a missing descriptor.  One unit is
    built, tagged for the recipe's sub_spec."""
    ref = _ref(entries=[
        ReferenceEntry(element="Si", atom_site=1, label="a",
                       default=True, description="d",
                       fingerprints=[])])
    options = {"xccode": 100, "imago_commit": "abc", "converg": 1.0e-6}
    units = build_loen_units(ref, "si.skel", options, [_BISPEC_DECL])

    assert len(units) == 1
    assert units[0].kind == "fingerprint"
    assert units[0].options["job"] == "loen"
    assert units[0].calc == (
        "loen-bispectrum-cutoff_9-twoj1_4-twoj2_4",)


def test_build_loen_units_dedups_recipe_and_override():
    """A per-entry override that repeats the recipe's (method, sub_spec)
    collapses to a single loen unit -- one descriptor table serves both,
    so the calc-tag dedup spans the characterization and per-entry
    sources, not just within each."""
    ref = _ref(entries=[
        ReferenceEntry(element="Si", atom_site=1, label="a",
                       default=True, description="d",
                       fingerprints=[_BISPEC_DECL])])
    options = {"xccode": 100, "imago_commit": "abc", "converg": 1.0e-6}
    units = build_loen_units(ref, "si.skel", options, [_BISPEC_DECL])

    assert len(units) == 1


def test_sub_spec_slug_is_slug_safe_and_ordered():
    # Floats format %.6g and any forbidden character (a dot) is
    #   sanitized to "_", so the tag is always a valid run-dir slug.
    slug = bip._sub_spec_slug(
        {"twoj2": 4, "cutoff": 7.5, "twoj1": 8})
    assert slug == "cutoff_7_5-twoj1_8-twoj2_4"


def test_harvest_fingerprints_no_declarations_is_empty():
    # An entry that declares no fingerprints harvests nothing and
    #   never touches the structure (the common case so far), so the
    #   result_toml is not even read.
    entry = _ref().entries[0]
    assert entry.fingerprints == []
    # Empty recipe and no per-entry override -> nothing to harvest.
    assert harvest_fingerprints(
        None, _ref(), entry.atom_site, entry.fingerprints,
        {}, []) == []


def _fake_two_atom_structure():
    """A duck-typed structure (Si at row 1, O at row 2) exposing just
    what ReduceMatcher.compute_query reads: the 1-indexed identity
    arrays and the minimum-image distance matrix.  Stands in for the
    StructureControl the harvest would read from the run's expanded
    imago.fract-mi, so the matcher path is tested without a real
    structure read (which needs the space-group database)."""

    return types.SimpleNamespace(
        num_atoms=2, num_elements=2,
        atom_element_id=[None, 1, 2],
        atom_species_id=[None, 1, 1],
        atom_element_name=[None, "Si", "O"],
        min_dist=[[None, None, None],
                  [None, 0.0, 2.0],
                  [None, 2.0, 0.0]])


_REDUCE_DECL = ManifestFingerprint(
    method="reduce",
    sub_spec={"level": 1, "thick": 0.1, "cutoff": 5.0,
              "tolerance": 0.05})


def test_read_skeleton_to_dat_map_keys_by_skeleton(tmp_path):
    """datSkl.map (DAT# SKELETON# ELEMENT SPECIES TYPE) reads into
    {skeleton: (dat, element)}, skipping the header line."""
    dat_skl = tmp_path / "datSkl.map"
    dat_skl.write_text("DAT SKEL ELEM SPECIES TYPE\n"
                       "  2    1   Si       1    1\n"
                       "  1    2   O        1    1\n")
    result_toml = {"outputs": {"datSkl_map": str(dat_skl)}}
    assert bip.read_skeleton_to_dat_map(result_toml) == {
        1: (2, "Si"), 2: (1, "O")}


def test_harvest_fingerprints_reduce(tmp_path, monkeypatch):
    """A reduce declaration harvests one element-only shell_code
    FingerprintRecord for the named site, mapping atom_site (skeleton)
    to the structure row via datSkl.map and computing in-process via
    ReduceMatcher (DESIGN 5.7 / 5.2)."""
    dat_skl = tmp_path / "datSkl.map"
    dat_skl.write_text("DAT SKEL ELEM SPECIES TYPE\n"
                       "  1    1   Si       1    1\n"
                       "  2    2   O        1    1\n")
    result_toml = {"outputs": {"structure": "imago.fract-mi",
                               "datSkl_map": str(dat_skl)}}
    monkeypatch.setattr(
        bip, "_read_structure_with_distances",
        lambda path, cutoff: _fake_two_atom_structure())

    spec = ReferenceEntry(
        element="Si", atom_site=1, label="t", default=True,
        description="d", fingerprints=[_REDUCE_DECL])
    # An empty recipe, so only the per-entry override is harvested.
    records = harvest_fingerprints(
        None, _ref(), spec.atom_site, spec.fingerprints,
        result_toml, [])

    assert len(records) == 1
    record = records[0]
    assert record.method == "reduce"
    assert record.preferred is False     # per-entry overrides aren't
    assert record.sub_spec == _REDUCE_DECL.sub_spec
    # Si at site 1 sees a single O neighbor at 2.0 Angstrom; the
    #   stored shell_code is element-only and lowercased.
    assert record.payload == {
        "shell_code": {
            "element": "si",
            "levels": [{"distance": pytest.approx(2.0),
                        "neighbors": ["o"]}],
        }}


def _write_loen_descriptor(tmp_path, reference_id, sub_spec, rows_text):
    """Create the loen unit's run directory under a flight root and drop
    a fake ``gs_loen-fb.plot`` descriptor in it, returning the flight.
    Mirrors what kaleidoscope would have produced for the loen unit
    build_loen_units dispatched for ``(reference_id, bispectrum,
    sub_spec)``."""
    calc_tag = bip._loen_calc_tag("bispectrum", sub_spec)
    run_dir = tmp_path / "wingbeats" / reference_id / calc_tag
    run_dir.mkdir(parents=True)
    header = ("site# element species type_in_species type_flat  "
              "c0 c1 c2 c3 c4  sum\n")
    (run_dir / "gs_loen-fb.plot").write_text(header + rows_text)
    return Flight(root=str(tmp_path), units=[])


def test_harvest_fingerprints_loen_side(tmp_path):
    """A bispectrum declaration harvests the descriptor row of the loen
    unit kaleidoscope ran: the harvest reconstructs the run dir from the
    calc tag, finds the ``*loen*.plot`` descriptor, maps atom_site
    (skeleton) to its dat row via datSkl.map, and wraps that row's
    vector as a ``values`` payload (DESIGN 5.10.3 / 5.2)."""
    # datSkl.map: skeleton site 1 -> dat row 2 (Si).
    dat_skl = tmp_path / "datSkl.map"
    dat_skl.write_text("DAT SKEL ELEM SPECIES TYPE\n"
                       "  2    1   Si       1    1\n"
                       "  1    2   O        1    1\n")
    result_toml = {"outputs": {"datSkl_map": str(dat_skl)}}

    sub_spec = {"twoj1": 4, "twoj2": 4, "cutoff": 9.0}
    flight = _write_loen_descriptor(
        tmp_path, "au_fcc", sub_spec,
        "1 O  1 1 1   9 9 9 9 9   0\n"
        "2 Si 1 1 2   1 2 3 4 5   0\n")

    spec = ReferenceEntry(
        element="Si", atom_site=1, label="t", default=True,
        description="d",
        fingerprints=[ManifestFingerprint(
            method="bispectrum", sub_spec=sub_spec)])
    records = harvest_fingerprints(
        flight, _ref(), spec.atom_site, spec.fingerprints,
        result_toml, [])

    assert len(records) == 1
    assert records[0].method == "bispectrum"
    assert records[0].preferred is False    # per-entry overrides aren't
    assert records[0].sub_spec == sub_spec
    # atom_site 1 -> dat row 2 (Si) -> that row's five components.
    assert records[0].payload == {
        "values": [1.0, 2.0, 3.0, 4.0, 5.0]}


def test_harvest_fingerprints_loen_guards_numbering_desync(tmp_path):
    """When the descriptor row's self-describing element disagrees with
    datSkl.map for a site, the numbering has desynced and the loen
    harvest refuses rather than storing the wrong atom's fingerprint."""
    dat_skl = tmp_path / "datSkl.map"
    dat_skl.write_text("DAT SKEL ELEM SPECIES TYPE\n"
                       "  1    1   Si       1    1\n")   # map says Si
    result_toml = {"outputs": {"datSkl_map": str(dat_skl)}}

    sub_spec = {"twoj1": 4, "twoj2": 4, "cutoff": 9.0}
    flight = _write_loen_descriptor(                     # row says O
        tmp_path, "au_fcc", sub_spec,
        "1 O 1 1 1   1 2 3 4 5   0\n")

    spec = ReferenceEntry(
        element="Si", atom_site=1, label="t", default=True,
        description="d",
        fingerprints=[ManifestFingerprint(
            method="bispectrum", sub_spec=sub_spec)])
    with pytest.raises(ValueError):
        harvest_fingerprints(
            flight, _ref(), spec.atom_site, spec.fingerprints,
            result_toml, [])


def test_harvest_fingerprints_guards_numbering_desync(
        tmp_path, monkeypatch):
    """When the expanded structure row and datSkl.map name different
    elements for a site, the numbering has desynced and the harvest
    refuses rather than describing the wrong atom."""
    dat_skl = tmp_path / "datSkl.map"
    dat_skl.write_text("DAT SKEL ELEM SPECIES TYPE\n"
                       "  1    1   Au       1    1\n"   # map says Au
                       "  2    2   O        1    1\n")
    result_toml = {"outputs": {"structure": "imago.fract-mi",
                               "datSkl_map": str(dat_skl)}}
    monkeypatch.setattr(                            # structure says Si
        bip, "_read_structure_with_distances",
        lambda path, cutoff: _fake_two_atom_structure())
    spec = ReferenceEntry(
        element="Au", atom_site=1, label="t", default=True,
        description="d", fingerprints=[_REDUCE_DECL])
    with pytest.raises(ValueError):
        harvest_fingerprints(
            None, _ref(), spec.atom_site, spec.fingerprints,
            result_toml, [])


# ---- discover_environments (B1: one rep per distinct env) ----

def _identity(mapping):
    """A datSkl.map reader stub returning the given
    {skeleton: (element, species, type)} for any result_toml."""
    return lambda result_toml: dict(mapping)


def test_discover_environments_one_rep_per_type():
    # Four sites span two distinct (element, species, type) groups:
    #   sites 1 & 3 are Si type 1; sites 2 & 4 are Si type 2.
    #   Discovery yields one representative per group -- the lowest
    #   skeleton index -- regardless of the map's row order.
    ref = _ref(entries=[])
    identity = _identity({3: ("si", 1, 1), 1: ("si", 1, 1),
                          4: ("si", 1, 2), 2: ("si", 1, 2)})
    envs = bip.discover_environments({}, ref, identity_fn=identity)

    assert [(e.atom_site, e.type_number) for e in envs] == [
        (1, 1), (2, 2)]
    # An auto-discovered environment is not default, carries no
    #   overrides, and derives its label + a non-empty description.
    assert all(e.default is False and e.overrides == []
               for e in envs)
    assert envs[0].label == "au_fcc-si1-t1-a1"
    assert envs[0].description != ""


def test_discover_environments_layers_pinned_customization():
    # A customization pins site 4 (Si type 2) and supplies a label,
    #   default flag, description, and a fingerprint override; the
    #   discovered environment uses that site as its representative
    #   and layers the rest on, while the other group stays auto.
    override = ManifestFingerprint(
        method="reduce",
        sub_spec={"level": 1, "thick": 0.1, "cutoff": 5.0,
                  "tolerance": 0.05})
    ref = _ref(entries=[ReferenceEntry(
        element="Si", atom_site=4, label="picked", default=True,
        description="custom", fingerprints=[override])])
    identity = _identity({1: ("si", 1, 1), 2: ("si", 1, 1),
                          3: ("si", 1, 2), 4: ("si", 1, 2)})
    envs = bip.discover_environments({}, ref, identity_fn=identity)

    by_type = {e.type_number: e for e in envs}
    # The type-1 group is auto-discovered (rep = lowest site 1).
    assert by_type[1].atom_site == 1
    assert by_type[1].default is False
    # The type-2 group takes the curator's pinned site and overrides.
    assert by_type[2].atom_site == 4
    assert by_type[2].label == "picked"
    assert by_type[2].default is True
    assert by_type[2].description == "custom"
    assert by_type[2].overrides == [override]


def test_discover_environments_skips_site_less_customization():
    # A site-less customization cannot yet be matched to an
    #   environment, so it is ignored; discovery still yields the
    #   run's environments unannotated (interim limitation).
    ref = _ref(entries=[ReferenceEntry(
        element=None, atom_site=None, label="floating",
        default=True, description="x")])
    identity = _identity({1: ("si", 1, 1)})
    envs = bip.discover_environments({}, ref, identity_fn=identity)

    assert len(envs) == 1
    assert envs[0].label == "au_fcc-si1-t1-a1"
    assert envs[0].default is False


def test_discover_environments_element_mismatch_raises():
    # A customization naming an element that disagrees with its
    #   pinned site is a hard error.
    ref = _ref(entries=[ReferenceEntry(
        element="O", atom_site=1, label=None, default=False,
        description=None)])
    identity = _identity({1: ("si", 1, 1)})
    with pytest.raises(ValueError):
        bip.discover_environments({}, ref, identity_fn=identity)


def test_discover_environments_two_customizations_one_env_raises():
    # Two customizations pinning sites in the same environment is
    #   ambiguous and refused.
    ref = _ref(entries=[
        ReferenceEntry(element="Si", atom_site=1, label="a",
                       default=False, description=None),
        ReferenceEntry(element="Si", atom_site=2, label="b",
                       default=False, description=None)])
    identity = _identity({1: ("si", 1, 1), 2: ("si", 1, 1)})
    with pytest.raises(ValueError):
        bip.discover_environments({}, ref, identity_fn=identity)


def test_discover_environments_pinned_site_absent_raises():
    # A customization pinning a site the converged run does not
    #   contain is a hard error.
    ref = _ref(entries=[ReferenceEntry(
        element="Si", atom_site=9, label=None, default=False,
        description=None)])
    identity = _identity({1: ("si", 1, 1)})
    with pytest.raises(ValueError):
        bip.discover_environments({}, ref, identity_fn=identity)


# ---- insert_or_skip / find_bispectrum_duplicate (B2) ---------

_BISPEC_SUBSPEC = {"twoj1": 6, "twoj2": 4, "cutoff": 9.0}


def _entry(label, vector, *, default=False):
    """A PotentialEntry carrying one preferred bispectrum
    fingerprint with the given descriptor vector."""
    return ipdb.PotentialEntry(
        label=label, default=default, description="d",
        num_gaussians=1, alpha_min=0.15, alpha_max=0.15,
        coefficients=[1.0], alphas=[0.15],
        provenance={"source": "Imago"},
        fingerprints=[ipdb.FingerprintRecord(
            method="bispectrum", sub_spec=dict(_BISPEC_SUBSPEC),
            preferred=True, payload={"values": list(vector)})])


def _db_with(*entries):
    db = ipdb.ElementDatabase(2, "Si", 14.0, 20.0, 1.0)
    db.potentials.extend(entries)
    return db


def test_find_bispectrum_duplicate_within_floor():
    # A stored entry within the floor (0.10) of the new entry's
    #   descriptor is reported as the duplicate.
    db = _db_with(_entry("a", [0.0, 0.0, 0.0]))
    new = _entry("b", [0.0, 0.0, 0.05])          # L2 = 0.05
    dup = bip.find_bispectrum_duplicate(db, new)
    assert dup is not None and dup.label == "a"


def test_find_bispectrum_duplicate_beyond_floor_is_none():
    db = _db_with(_entry("a", [0.0, 0.0, 0.0]))
    new = _entry("b", [0.0, 0.0, 0.5])           # L2 = 0.5 > floor
    assert bip.find_bispectrum_duplicate(db, new) is None


def test_find_bispectrum_duplicate_returns_nearest():
    db = _db_with(_entry("near", [0.0, 0.0, 0.02]),
                  _entry("far", [0.0, 0.0, 0.09]))
    new = _entry("b", [0.0, 0.0, 0.0])
    assert bip.find_bispectrum_duplicate(db, new).label == "near"


def test_find_bispectrum_duplicate_ignores_keyless_entries():
    # The isolated baseline (no bispectrum) is never a match target,
    #   and a new entry without a bispectrum key never matches.
    iso = ipdb.PotentialEntry(
        "isolated", True, "iso", 1, 0.15, 0.15,
        [1.0], [0.15], {"source": "atomSCF"})
    assert bip.find_bispectrum_duplicate(
        _db_with(iso), _entry("b", [0.0, 0.0, 0.0])) is None
    assert bip.find_bispectrum_duplicate(
        _db_with(_entry("a", [0.0, 0.0, 0.0])), iso) is None


def test_insert_or_skip_appends_novel():
    db = _db_with(_entry("a", [0.0, 0.0, 0.0]))
    bip.insert_or_skip(db, _entry("b", [1.0, 0.0, 0.0]))
    assert [e.label for e in db.potentials] == ["a", "b"]


def test_insert_or_skip_skips_bispectrum_duplicate():
    db = _db_with(_entry("a", [0.0, 0.0, 0.0]))
    bip.insert_or_skip(db, _entry("b", [0.0, 0.0, 0.03]))
    # b duplicates a within the floor, so it is dropped; a stands.
    assert [e.label for e in db.potentials] == ["a"]


def test_insert_or_skip_replaces_same_label():
    # A same label (a re-harvested solid or a curator override)
    #   replaces in place, ahead of the bispectrum dedup.
    db = _db_with(_entry("a", [0.0, 0.0, 0.0]))
    bip.insert_or_skip(db, _entry("a", [5.0, 5.0, 5.0]))
    assert [e.label for e in db.potentials] == ["a"]
    assert ipdb.lookup(db, "a").fingerprints[0].payload == {
        "values": [5.0, 5.0, 5.0]}


def test_make_imago_provenance_records_symmetry_type_assignment():
    # The producer assigns types crystallographically, so every
    #   harvested entry records type_assignment="symmetry"; from it
    #   the native/witness role is derived (DESIGN 5.2.2), and both
    #   reduce and bispectrum are (exact) witnesses for these runs.
    prov = bip.make_imago_provenance("sha", "ts", _ref(), 1, 7)
    assert prov["type_assignment"] == "symmetry"


def test_materialize_structure_resolves_local_path(tmp_path):
    # The structure_path branch is a plain join under the manifest
    #   dir -- no network.
    path = materialize_structure(
        _ref(structure_path="sub/au.skel"),
        manifest_dir="/manifests", pdb_root="/data/atomicPDB")
    assert path == os.path.join("/manifests", "sub", "au.skel")


def test_materialize_structure_converts_cod_to_skl(
        tmp_path, monkeypatch):
    # The cod_id branch fetches the CIF and converts it to a skeleton
    # via cif2skl, returning the .skl path (fetch + convert are stubbed
    # so the wiring is tested without network or the ASE/binary stack).
    import cif2skl
    calls = {}

    def fake_fetch(cod_id, cod_revision, dest):
        calls["fetch"] = (cod_id, cod_revision, dest)
        with open(dest, "w") as handle:
            handle.write("# fake cif\n")

    def fake_convert(cif_path, skl_path, title=None):
        calls["convert"] = (cif_path, skl_path, title)
        with open(skl_path, "w") as handle:
            handle.write("title\nx\nend\n")
        return "227_a"

    monkeypatch.setattr(bip, "_fetch_cod_structure", fake_fetch)
    monkeypatch.setattr(cif2skl, "convert", fake_convert)
    pdb_root = str(tmp_path / "atomicPDB")
    ref = _ref(structure_path=None, cod_id=9008463,
               cod_revision="291735", reference_id="au_fcc")
    path = materialize_structure(
        ref, manifest_dir=str(tmp_path), pdb_root=pdb_root)
    assert path.endswith("au_fcc.skl")
    assert os.path.exists(path)
    assert calls["fetch"][0] == 9008463
    assert calls["convert"][0].endswith("au_fcc.cif")
    assert calls["convert"][2] == "au_fcc"          # title


def test_materialize_cod_conversion_failure_is_fatal(
        tmp_path, monkeypatch):
    # An unresolvable space group is a hard error that points the
    # curator at the structure_path escape hatch.
    import cif2skl

    monkeypatch.setattr(
        bip, "_fetch_cod_structure",
        lambda cod_id, rev, dest: open(dest, "w").write("# cif\n"))

    def boom(cif_path, skl_path, title=None):
        raise cif2skl.CifConversionError("no variant verified")

    monkeypatch.setattr(cif2skl, "convert", boom)
    ref = _ref(structure_path=None, cod_id=1, cod_revision="1",
               reference_id="bad_solid")
    with pytest.raises(RuntimeError, match="structure_path"):
        materialize_structure(
            ref, manifest_dir=str(tmp_path),
            pdb_root=str(tmp_path / "atomicPDB"))


def test_curation_workspace_root_sits_beside_databases():
    root = curation_workspace_root("/data/atomicPDB")
    assert root == os.path.join("/data", "curation", "workspace")


def test_structure_cache_dir_sits_beside_databases():
    cache = structure_cache_dir("/data/atomicPDB")
    assert cache == os.path.join(
        "/data", "atomicBDB", "cache", "structures")


# ============================================================
#  Structure pre-flight: load_structure_sources + materialize_only
# ============================================================

# A relaxed-load manifest carrying ONLY the structure sources -- no
#   kpoint_spec / scf_threshold / sub-model / entries.  The full
#   loader rejects it (rule 2); the pre-flight loader accepts it.
_SOURCES_ONLY_MANIFEST = """\
schema_version = 2

[[reference_solid]]
reference_id = "si_diamond"
cod_id = 2104737
cod_revision = "201401"

[[reference_solid]]
reference_id = "a_si"
structure_path = "a_si.skl"
"""


class TestLoadStructureSources:
    """The relaxed reader behind --materialize-only: it parses the
    structure sources from a manifest that lacks the run and harvest
    fields the full loader demands."""

    def test_parses_sources_without_run_fields(self, tmp_path):
        # A bare structure_path file must exist for nothing here --
        #   the loader defers that check to materialize time.
        path = _write(tmp_path, _SOURCES_ONLY_MANIFEST)
        sources = load_structure_sources(path)
        assert [s.reference_id for s in sources] == [
            "si_diamond", "a_si"]
        assert sources[0].cod_id == 2104737
        assert sources[0].cod_revision == "201401"
        assert sources[0].structure_path is None
        assert sources[1].structure_path == "a_si.skl"
        assert sources[1].cod_id is None

    def test_full_loader_rejects_what_relaxed_accepts(self, tmp_path):
        # The same source-only manifest the relaxed loader accepts is
        #   rejected by the full loader (missing rule-2 fields), which
        #   is exactly why the pre-flight needs its own reader.
        path = _write(tmp_path, _SOURCES_ONLY_MANIFEST)
        with pytest.raises(ValueError, match="rule 2"):
            load_manifest_v2(path)

    def test_missing_reference_id_raises(self, tmp_path):
        path = _write(tmp_path, 'schema_version = 2\n\n'
                      '[[reference_solid]]\ncod_id = 1\n'
                      'cod_revision = "1"\n')
        with pytest.raises(ValueError, match="reference_id"):
            load_structure_sources(path)

    def test_both_sources_raises(self, tmp_path):
        path = _write(tmp_path, 'schema_version = 2\n\n'
                      '[[reference_solid]]\nreference_id = "x"\n'
                      'cod_id = 1\ncod_revision = "1"\n'
                      'structure_path = "x.skl"\n')
        with pytest.raises(ValueError, match="rule 4"):
            load_structure_sources(path)

    def test_duplicate_reference_id_raises(self, tmp_path):
        path = _write(tmp_path, 'schema_version = 2\n\n'
                      '[[reference_solid]]\nreference_id = "x"\n'
                      'structure_path = "x.skl"\n\n'
                      '[[reference_solid]]\nreference_id = "x"\n'
                      'structure_path = "y.skl"\n')
        with pytest.raises(ValueError, match="duplicate"):
            load_structure_sources(path)

    def test_wrong_schema_version_raises(self, tmp_path):
        path = _write(tmp_path, 'schema_version = 1\n')
        with pytest.raises(ValueError, match="rule 1"):
            load_structure_sources(path)


class TestMaterializeOnly:
    """The pre-flight orchestration: it materializes every source and
    reports per-solid, continuing past failures."""

    def test_converts_each_source_into_cache_dir(
            self, tmp_path, monkeypatch):
        # Stub the network fetch + the converter so the wiring (and
        #   the cache_dir redirect) is tested without ASE or a binary.
        import cif2skl

        def fake_fetch(cod_id, cod_revision, dest):
            with open(dest, "w") as handle:
                handle.write("# fake cif\n")

        def fake_convert(cif_path, skl_path, title=None):
            with open(skl_path, "w") as handle:
                handle.write("title\nx\nend\n")
            return "227_a"

        monkeypatch.setattr(bip, "_fetch_cod_structure", fake_fetch)
        monkeypatch.setattr(cif2skl, "convert", fake_convert)

        # Give the structure_path solid a real file to resolve.
        (tmp_path / "a_si.skl").write_text("title\nx\nend\n")
        manifest = _write(tmp_path, _SOURCES_ONLY_MANIFEST)
        mirror = tmp_path / "mirror"

        report = materialize_only(
            manifest, str(tmp_path / "atomicPDB"),
            cache_dir=str(mirror))

        assert all(row["ok"] for row in report)
        assert {row["reference_id"] for row in report} == {
            "si_diamond", "a_si"}
        # The cod_id solid's CIF + skl land in the redirect mirror.
        assert (mirror / "si_diamond.cif").exists()
        assert (mirror / "si_diamond.skl").exists()
        # The structure_path solid resolves against the manifest dir.
        diamond = next(r for r in report
                       if r["reference_id"] == "si_diamond")
        assert diamond["source"].startswith("cod_id 2104737")

    def test_reports_failure_and_continues(
            self, tmp_path, monkeypatch):
        # One unresolvable space group must not hide the rest: the
        #   loop captures the failure and keeps going.
        import cif2skl

        monkeypatch.setattr(
            bip, "_fetch_cod_structure",
            lambda cod_id, rev, dest: open(dest, "w").write("# c\n"))

        def boom(cif_path, skl_path, title=None):
            raise cif2skl.CifConversionError("no variant verified")

        monkeypatch.setattr(cif2skl, "convert", boom)

        (tmp_path / "a_si.skl").write_text("title\nx\nend\n")
        manifest = _write(tmp_path, _SOURCES_ONLY_MANIFEST)

        report = materialize_only(
            manifest, str(tmp_path / "atomicPDB"),
            cache_dir=str(tmp_path / "mirror"))

        by_id = {row["reference_id"]: row for row in report}
        assert by_id["si_diamond"]["ok"] is False
        assert "structure_path" in by_id["si_diamond"]["message"]
        # The local structure_path solid still resolves cleanly.
        assert by_id["a_si"]["ok"] is True

    def test_cli_materialize_only_skips_dispatch(
            self, tmp_path, monkeypatch):
        # --materialize-only must short-circuit before any SCF: if it
        #   reached build_initial_potentials the stub would raise.
        import cif2skl

        monkeypatch.setattr(
            bip, "_fetch_cod_structure",
            lambda cod_id, rev, dest: open(dest, "w").write("# c\n"))
        monkeypatch.setattr(
            cif2skl, "convert",
            lambda c, s, title=None: open(s, "w").write("t\n") or "1")

        def explode(*a, **k):
            raise AssertionError("full build must not run")

        monkeypatch.setattr(bip, "build_initial_potentials", explode)

        (tmp_path / "a_si.skl").write_text("t\n")
        manifest = _write(tmp_path, _SOURCES_ONLY_MANIFEST)
        rc = bip.main([
            "--manifest", manifest,
            "--pdb-root", str(tmp_path / "atomicPDB"),
            "--materialize-only",
            "--materialize-dir", str(tmp_path / "mirror")])
        assert rc == 0

    def test_cli_materialize_dir_requires_materialize_only(
            self, tmp_path):
        # --materialize-dir without --materialize-only is a usage
        #   error (argparse exits non-zero via SystemExit).
        with pytest.raises(SystemExit):
            bip.main([
                "--manifest", str(tmp_path / "m.toml"),
                "--pdb-root", str(tmp_path / "atomicPDB"),
                "--materialize-dir", str(tmp_path / "mirror")])


def _write_scfv(path, body):
    """Write a multi-type scfV output file (the NUM_TYPES /
    TOTAL__OR__SPIN_UP / SPIN_DN layout Imago writes from fort.8)."""
    path.write_text(body)


def test_extract_potential_reads_scfv_type_block(tmp_path):
    """extract_potential selects the named site's type block from
    the multi-type scfV output (NUM_TYPES header + the
    TOTAL__OR__SPIN_UP channel), resolving the site's type via
    datSkl.map and taking cols 1-2 (coeff, alpha).  The redundant
    SPIN_DN channel is ignored."""
    scfv = tmp_path / "scfV.dat"
    _write_scfv(scfv,
                "NUM_TYPES    1\n"
                "TOTAL__OR__SPIN_UP\n"
                "   2\n"
                " 0.5 1.0 0 0 0\n"
                " 0.3 2.0 0 0 0\n"
                "SPIN_DN\n"
                "   2\n"
                " 9.9 1.0 0 0 0\n"        # must NOT be read
                " 9.9 2.0 0 0 0\n")
    datskl = tmp_path / "datSkl.map"
    datskl.write_text("DAT SKEL ELEM SPECIES TYPE\n"
                      "  1    1   Si       1    1\n")
    result_toml = {"outputs": {"scfV": str(scfv),
                               "datSkl_map": str(datskl)}}
    coeffs, alphas = extract_potential(result_toml, atom_site=1)
    assert coeffs == [0.5, 0.3]
    assert alphas == [1.0, 2.0]


def test_extract_potential_selects_the_right_type_block(tmp_path):
    """With several types in one scfV file, the site's type number
    (from datSkl.map) picks the correct block -- not just the
    first."""
    scfv = tmp_path / "scfV.dat"
    _write_scfv(scfv,
                "NUM_TYPES    2\n"
                "TOTAL__OR__SPIN_UP\n"
                "   1\n"
                " 1.1 0.5 0 0 0\n"        # type 1
                "   2\n"
                " 2.1 0.5 0 0 0\n"        # type 2
                " 2.2 1.5 0 0 0\n"
                "SPIN_DN\n"
                "   1\n"
                " 9.9 0.5 0 0 0\n"
                "   2\n"
                " 9.9 0.5 0 0 0\n"
                " 9.9 1.5 0 0 0\n")
    datskl = tmp_path / "datSkl.map"
    datskl.write_text("DAT SKEL ELEM SPECIES TYPE\n"
                      "  1    1   Si       1    1\n"
                      "  2    2   Si       2    2\n")
    result_toml = {"outputs": {"scfV": str(scfv),
                               "datSkl_map": str(datskl)}}
    coeffs, alphas = extract_potential(result_toml, atom_site=2)
    assert coeffs == [2.1, 2.2]
    assert alphas == [0.5, 1.5]


# ---- site identity (datSkl.map) + label assembly (C87) -------

def test_read_site_identity_map_keys_by_skeleton_number(tmp_path):
    # The five-column datSkl.map makeinput now writes: a header
    # line plus DAT#  SKELETON#  ELEMENT  SPECIES  TYPE.  The
    # reader keys by the skeleton number (column 2), since the
    # producer harvests by atom_site (skeleton numbering).
    dat_skl = tmp_path / "datSkl.map"
    dat_skl.write_text(
        "      DAT#  SKELETON#    ELEMENT    SPECIES       TYPE\n"
        "         1          2         si          1          1\n"
        "         2          1          o          1          2\n")
    identity = read_site_identity_map(
        {"outputs": {"datSkl_map": str(dat_skl)}})
    assert identity == {2: ("si", 1, 1), 1: ("o", 1, 2)}


def test_assemble_entry_label_builds_the_5_2_1_form():
    # <reference_id>-<element><species>-t<type>-a<site>, lowercased.
    assert assemble_entry_label(
        "si_diamond", "Si", 1, 1, 1) == "si_diamond-si1-t1-a1"
    assert assemble_entry_label(
        "forsterite", "Mg", 1, 2, 2) == "forsterite-mg1-t2-a2"


# ---- run log -------------------------------------------------

def test_write_run_log_round_trips(tmp_path):
    log_path = str(tmp_path / "curation" / "run_log.toml")
    rows = [
        # record_converged's output shape: the converged mesh, its
        #   k-density, and the flatness ladder (only mesh + density
        #   reach the run log).
        make_run_log_entry(
            _ref(),
            {"converged_mesh": [4, 4, 4],
             "converged_kpoint_density": 100,
             "grid_values": [], "grid_energies": []},
            {"scf_iterations": 7}),
        make_nonconverged_log_entry(_ref(reference_id="ag_fcc")),
    ]
    write_run_log(log_path, "abc123", "2026-06-12T00:00:00Z", rows)
    with open(log_path, "rb") as handle:
        data = tomllib.load(handle)
    assert data["imago_commit"] == "abc123"
    assert data["run"][0]["reference_id"] == "au_fcc"
    assert data["run"][0]["converged"] is True
    assert data["run"][0]["converged_kpoint_density"] == 100
    assert data["run"][0]["converged_mesh"] == [4, 4, 4]
    assert data["run"][1]["converged"] is False


# ---- the orchestrator (toolchain seam mocked) ----------------

_AU_LOCAL_MANIFEST = (
    "schema_version = 2\n\n"
    + _CHAR_BLOCK +
    """
[[reference_solid]]
reference_id = "au_fcc"
system_type = "crystalline"
basis = "fb"
functional = "wigner"
kpoint_integration = "linear-tetrahedral"
structure_path = "au.skel"
kpoint_spec = { density = 60.0, shift = [0.0, 0.0, 0.0] }
scf_threshold = 1.0e-6

  [[reference_solid.entry]]
  element = "Au"
  atom_site = 1
  label = "default_solid"
  default = true
  description = "Au in fcc bulk."
""")


def _install_climb_mocks(monkeypatch, workspace, *,
                         material="au_fcc", converged_mesh=(4, 4, 4)):
    """Install the climb-flow seams for an orchestration test.

    The producer's build phase predicts a seed and assembles a
    ``ClimbConfig`` per solid, then ``converge_by_climb`` drives the
    climb; all three need a live imago and the guidance dataspace, so
    they are mocked deterministically here (the climb's own logic is
    exercised by the converge_by_climb unit tests below).  The mocked
    ``build_climb_config`` picks ``recip_cell_volume`` so the
    converged mesh's k-density is a round number the run-log
    assertions can name: ``product([4,4,4]) / 0.64 == 100``.  The
    converged mesh's ``result.toml`` is written where Phase 3 rebuilds
    the unit and reads it back.  Returns a dict recording whether the
    guidance contribution fired."""

    def fake_predict(struct, dataspace, system_type, submodel,
                     center=None):
        seed = float(center) if center is not None else 100.0
        record = PredictionRecord(
            policy=("curator_override" if center is not None
                    else "verify_around_prediction"),
            predicted_kpoint_density=seed, confidence=0.9,
            is_under_trained=False, system_type=system_type,
            basis=submodel["basis"],
            functional=submodel["functional"],
            kpoint_integration=submodel["kpoint_integration"])
        return seed, 0.9, False, record

    monkeypatch.setattr(bip, "predict_kpoint_density", fake_predict)

    config = bip.ClimbConfig(
        classes=[0, 0, 0], recip_mag=[1.0, 1.0, 1.0],
        recip_cell_volume=0.64, mode=mesh_climb.UNIT_STEP,
        flat_needed=1, grid_width=0, start_offset=1, max_stride=8,
        cell_atom_count=2, threshold=1.0e-6, max_count=20)
    monkeypatch.setattr(bip, "build_climb_config",
                        lambda *a, **k: config)

    # A converged rung plus the >= 3 distinct rungs its stop test
    #   required (DESIGN 3.12.3), so record_converged has a ladder.
    ladder = [bip.Rung([2, 2, 2], -1.0), bip.Rung([3, 3, 3], -1.0),
              bip.Rung(list(converged_mesh), -1.0)]

    def fake_converge(materials, configs, seeds, dispatcher,
                      on_non_converged=None):
        outcomes = {m: bip.Rung(list(converged_mesh), -1.0)
                    for m in materials}
        return outcomes, {m: ladder for m in materials}

    monkeypatch.setattr(bip, "converge_by_climb", fake_converge)

    mesh_tag = ("kpt-mesh-" + "-".join(str(c) for c in converged_mesh),)
    _write_result(workspace, material, mesh_tag, energy=-1.0)

    # The guidance entry builder and store are unit-tested in
    #   test_guidance_harvest; here we only confirm the producer wires
    #   the in-memory contribution through them.
    harvested = {}
    monkeypatch.setattr(bip.guidance_harvest, "build_entry",
                        lambda *a, **k: object())
    monkeypatch.setattr(
        bip.guidance_harvest, "save_entry",
        lambda entry, root: harvested.setdefault("called", True))
    return harvested


def test_build_initial_potentials_harvests_curated_entry(
        tmp_path, monkeypatch):
    """End-to-end producer wiring with the toolchain seam mocked:
    the builder, dispatch, the scfV read, and the guidance harvest
    are all stubbed, leaving the producer's own orchestration under
    test -- build the combined flight, pick the converged point,
    assemble + save the Imago-source entry, and write the run log."""

    data_root = str(tmp_path)
    pdb_root = os.path.join(data_root, "atomicPDB")
    _make_element(pdb_root, "au", [1.0, 2.0, 3.0],
                  [0.15, 1.5, 1.0e8])
    (tmp_path / "au.skel").write_text("dummy structure\n")
    manifest_path = _write(tmp_path, _AU_LOCAL_MANIFEST)

    # The predictor, ClimbConfig assembly, and climb all need a live
    #   imago and the guidance dataspace; install those seams mocked
    #   (predict -> config -> converge, plus the converged mesh's
    #   result.toml), and record that the guidance contribution fires.
    monkeypatch.setattr(
        bip.guidance_db, "load",
        lambda root: types.SimpleNamespace(group_table={}))
    workspace = bip.curation_workspace_root(pdb_root)
    harvested = _install_climb_mocks(monkeypatch, workspace)
    # build_entry (mocked) receives the loaded structure; a stub for
    #   the loader suffices since the mock ignores it.
    monkeypatch.setattr(
        bip.guidance_harvest, "load_structure",
        lambda path: types.SimpleNamespace(num_atoms=2))

    # The loen pre-flight is the only real dispatch (the climb is
    #   mocked); the fingerprint harvest is stubbed, so nothing reads
    #   its results.
    def fake_dispatch(flight, executor=None, force=False):
        pass

    build_initial_potentials(
        manifest_path, pdb_root, data_root,
        dispatch_fn=fake_dispatch,
        prepare_fn=lambda flight, workspace: None,
        extract_fn=lambda result, site: ([0.5, 0.3], [1.0, 2.0]),
        # Site 1 is Au (cross-checks the entry's element); the
        #   fingerprint harvest (recipe + override) needs a live run,
        #   so it is stubbed away here.
        identity_fn=lambda result: {1: ("au", 1, 1)},
        fingerprint_fn=lambda *args, **kwargs: [])

    # The Au database gained the curated default_solid entry.
    database = ipdb.load(element_path(pdb_root, "au"),
                         known_methods=None)
    entry = ipdb.lookup(database, "default_solid")
    assert entry.coefficients == [0.5, 0.3]
    assert entry.num_gaussians == 2
    assert entry.provenance["source"] == "Imago"
    assert entry.provenance["reference_id"] == "au_fcc"
    # The guidance contribution fired and the run log was written.
    assert harvested.get("called") is True
    with open(os.path.join(data_root, "curation", "run_log.toml"),
              "rb") as handle:
        run_log = tomllib.load(handle)
    assert run_log["run"][0]["converged"] is True
    assert run_log["run"][0]["converged_kpoint_density"] == 100


# A label-less variant of the local manifest: the entry omits
#   ``label`` so the producer must derive it at harvest (5.2.1).
_AU_LOCAL_MANIFEST_NO_LABEL = _AU_LOCAL_MANIFEST.replace(
    '  label = "default_solid"\n', "")


def test_build_initial_potentials_derives_label_at_harvest(
        tmp_path, monkeypatch):
    """When the manifest entry omits ``label``, the producer reads
    the run's site identity (datSkl.map, via ``identity_fn``) and
    assembles the DESIGN 5.2.1 label
    ``<reference_id>-<element><species>-t<type>-a<site>``."""

    data_root = str(tmp_path)
    pdb_root = os.path.join(data_root, "atomicPDB")
    _make_element(pdb_root, "au", [1.0, 2.0, 3.0],
                  [0.15, 1.5, 1.0e8])
    (tmp_path / "au.skel").write_text("dummy structure\n")
    manifest_path = _write(tmp_path, _AU_LOCAL_MANIFEST_NO_LABEL)

    monkeypatch.setattr(
        bip.guidance_db, "load",
        lambda root: types.SimpleNamespace(group_table={}))
    workspace = bip.curation_workspace_root(pdb_root)
    _install_climb_mocks(monkeypatch, workspace)
    monkeypatch.setattr(
        bip.guidance_harvest, "load_structure",
        lambda path: types.SimpleNamespace(num_atoms=2))

    def fake_dispatch(flight, executor=None, force=False):
        pass

    # Inject the site-identity reader: site 1 is Au species 1,
    #   type 1, so the derived label is au_fcc-au1-t1-a1.
    build_initial_potentials(
        manifest_path, pdb_root, data_root,
        dispatch_fn=fake_dispatch,
        prepare_fn=lambda flight, workspace: None,
        extract_fn=lambda result, site: ([0.5, 0.3], [1.0, 2.0]),
        identity_fn=lambda result: {1: ("au", 1, 1)},
        fingerprint_fn=lambda *args, **kwargs: [])

    database = ipdb.load(element_path(pdb_root, "au"),
                         known_methods=None)
    entry = ipdb.lookup(database, "au_fcc-au1-t1-a1")
    assert entry.default is True
    assert entry.coefficients == [0.5, 0.3]
    assert entry.provenance["reference_id"] == "au_fcc"


# ================================================================
#  The [defaults] hoist (DESIGN 5.7; C104) -- shared run settings
#  live once in a top-level [defaults] block and a solid inherits
#  any it omits.  The producer folds them into each solid up front
#  (apply_manifest_defaults) so the rest of the run reads one fully
#  resolved setting per field.
# ================================================================

# Two cod-sourced solids sharing a [defaults] block: au_fcc names
#   no run settings (inherits every one), cu_fcc overrides only the
#   functional (keeping the other four inherited).  Different
#   elements keep the per-element label rule (rule 6) trivial.
_DEFAULTS_TWO_SOLID_MANIFEST = (
    "schema_version = 2\n\n"
    + _CHAR_BLOCK +
    "\n[defaults]\n"
    "basis = \"fb\"\n"
    "functional = \"wigner\"\n"
    "kpoint_integration = \"linear-tetrahedral\"\n"
    "kpoint_spec = { density = 60.0, shift = [0.0, 0.0, 0.0] }\n"
    "scf_threshold = 1.0e-6\n\n"
    "[[reference_solid]]\n"
    "reference_id = \"au_fcc\"\n"
    "system_type = \"crystalline\"\n"
    "cod_id = 9008463\n"
    "cod_revision = \"2023-04-12\"\n\n"
    "[[reference_solid]]\n"
    "reference_id = \"cu_fcc\"\n"
    "system_type = \"crystalline\"\n"
    "functional = \"pz-lda\"\n"
    "cod_id = 9008468\n"
    "cod_revision = \"2023-04-12\"\n")


def test_apply_manifest_defaults_inherits_and_overrides(tmp_path):
    """apply_manifest_defaults folds the [defaults] block into each
    solid: a solid that omits a setting inherits the default, and a
    solid that names its own value keeps it.  After the pass, every
    run-setting field is populated on every solid (no None left)."""

    path = _write(tmp_path, _DEFAULTS_TWO_SOLID_MANIFEST)
    manifest = load_manifest_v2(path)
    # Before resolution the sparse solids carry None where they
    #   omitted a setting -- au_fcc omits all five.
    au_before = manifest.reference_solids[0]
    assert au_before.basis is None
    assert au_before.functional is None

    apply_manifest_defaults(manifest)

    au_fcc, cu_fcc = manifest.reference_solids
    # au_fcc inherited every shared setting verbatim.
    assert au_fcc.basis == "fb"
    assert au_fcc.functional == "wigner"
    assert au_fcc.kpoint_integration == "linear-tetrahedral"
    assert au_fcc.kpoint_spec == {"density": 60.0,
                                  "shift": [0.0, 0.0, 0.0]}
    assert au_fcc.scf_threshold == 1.0e-6
    # cu_fcc kept its own functional but inherited the rest.
    assert cu_fcc.functional == "pz-lda"
    assert cu_fcc.basis == "fb"
    assert cu_fcc.kpoint_integration == "linear-tetrahedral"
    assert cu_fcc.scf_threshold == 1.0e-6


# A [defaults]-driven variant of the local Au manifest: the solid
#   names no run settings at all, so the harvest can only succeed if
#   the producer resolved basis / kpoint_spec / scf_threshold (used
#   by the builder submodel and pick_converged_unit) from [defaults].
_AU_DEFAULTS_MANIFEST = (
    "schema_version = 2\n\n"
    + _CHAR_BLOCK +
    "\n[defaults]\n"
    "basis = \"fb\"\n"
    "functional = \"wigner\"\n"
    "kpoint_integration = \"linear-tetrahedral\"\n"
    "kpoint_spec = { density = 60.0, shift = [0.0, 0.0, 0.0] }\n"
    "scf_threshold = 1.0e-6\n"
    """
[[reference_solid]]
reference_id = "au_fcc"
system_type = "crystalline"
structure_path = "au.skel"

  [[reference_solid.entry]]
  element = "Au"
  atom_site = 1
  label = "default_solid"
  default = true
  description = "Au in fcc bulk."
""")


def test_build_initial_potentials_resolves_defaults(
        tmp_path, monkeypatch):
    """End-to-end producer run driven by a [defaults] manifest: the
    solid names no run settings, so the harvest succeeds only because
    the producer folded the shared [defaults] into it up front.  The
    curated entry lands with the resolved basis in its provenance."""

    data_root = str(tmp_path)
    pdb_root = os.path.join(data_root, "atomicPDB")
    _make_element(pdb_root, "au", [1.0, 2.0, 3.0],
                  [0.15, 1.5, 1.0e8])
    (tmp_path / "au.skel").write_text("dummy structure\n")
    manifest_path = _write(tmp_path, _AU_DEFAULTS_MANIFEST)

    monkeypatch.setattr(
        bip.guidance_db, "load",
        lambda root: types.SimpleNamespace(group_table={}))
    workspace = bip.curation_workspace_root(pdb_root)
    _install_climb_mocks(monkeypatch, workspace)
    monkeypatch.setattr(
        bip.guidance_harvest, "load_structure",
        lambda path: types.SimpleNamespace(num_atoms=2))

    def asserting_predict(struct, dataspace, system_type, submodel,
                          center=None):
        # The submodel must carry the resolved run settings, not
        #   None: assert the [defaults] flowed through to the
        #   predictor, and the resolved kpoint_spec.density is center.
        assert submodel["basis"] == "fb"
        assert submodel["functional"] == "wigner"
        assert submodel["kpoint_integration"] == "linear-tetrahedral"
        assert center == 60.0
        record = PredictionRecord(
            policy="curator_override",
            predicted_kpoint_density=60.0, confidence=0.9,
            is_under_trained=False, system_type=system_type,
            basis=submodel["basis"],
            functional=submodel["functional"],
            kpoint_integration=submodel["kpoint_integration"])
        return 60.0, 0.9, False, record

    monkeypatch.setattr(bip, "predict_kpoint_density",
                        asserting_predict)

    def fake_dispatch(flight, executor=None, force=False):
        pass

    build_initial_potentials(
        manifest_path, pdb_root, data_root,
        dispatch_fn=fake_dispatch,
        prepare_fn=lambda flight, workspace: None,
        extract_fn=lambda result, site: ([0.5, 0.3], [1.0, 2.0]),
        identity_fn=lambda result: {1: ("au", 1, 1)},
        fingerprint_fn=lambda *args, **kwargs: [])

    # The Au entry harvested -- the run could only reach harvest
    #   because scf_threshold resolved from [defaults] and the climb
    #   converged (mocked); the resolved basis rides in provenance.
    database = ipdb.load(element_path(pdb_root, "au"),
                         known_methods=None)
    entry = ipdb.lookup(database, "default_solid")
    assert entry.coefficients == [0.5, 0.3]
    assert entry.provenance["reference_id"] == "au_fcc"


# ================================================================
#  Dispatch wiring (DESIGN 6.2.11; PSEUDOCODE 13.7) -- the producer
#  change-over from the deleted curation_executor to attaching a
#  flight Parsl Config and letting the driver auto-select.  (The
#  shared resolve_dispatch / write_resolved_dispatch helpers are
#  exercised in test_cluster_config.py, where they now live.)
# ================================================================

def test_producer_local_default_attaches_no_config(monkeypatch,
                                                   tmp_path):
    """build_initial_potentials defaults to local: resolve_dispatch
    returns no Parsl Config, and that None (with the force
    cache-bypass) is threaded into the climb's round dispatcher."""

    data_root = str(tmp_path)
    pdb_root = os.path.join(data_root, "atomicPDB")
    (tmp_path / "au.skel").write_text("dummy structure\n")
    manifest_path = _write(tmp_path, _AU_LOCAL_MANIFEST)

    monkeypatch.setattr(
        bip.guidance_db, "load",
        lambda root: types.SimpleNamespace(group_table={}))
    monkeypatch.setattr(bip, "refresh_isolated_entries",
                        lambda *a, **k: {})
    # Mock the Phase-1 predict + config; converge returns
    #   NON_CONVERGED so Phase 3 skips the harvest entirely (no
    #   element database or result.toml needed).
    monkeypatch.setattr(
        bip, "predict_kpoint_density",
        lambda struct, ds, st, sm, center=None: (
            60.0, 0.9, False, PredictionRecord(
                policy="curator_override",
                predicted_kpoint_density=60.0, confidence=0.9,
                is_under_trained=False, system_type=st,
                basis=sm["basis"], functional=sm["functional"],
                kpoint_integration=sm["kpoint_integration"])))
    monkeypatch.setattr(
        bip, "build_climb_config",
        lambda *a, **k: bip.ClimbConfig(
            classes=[0, 0, 0], recip_mag=[1.0, 1.0, 1.0],
            recip_cell_volume=1.0, mode=mesh_climb.UNIT_STEP,
            flat_needed=1, grid_width=0, start_offset=1, max_stride=8,
            cell_atom_count=2, threshold=1.0e-6, max_count=20))
    monkeypatch.setattr(
        bip, "converge_by_climb",
        lambda materials, *a, **k: (
            {m: bip.NON_CONVERGED for m in materials}, {}))
    monkeypatch.setattr(bip, "save_databases", lambda *a, **k: None)
    monkeypatch.setattr(bip, "write_run_log", lambda *a, **k: None)

    # Spy the climb dispatcher's resolved config + force.  The climb
    #   is mocked, so the dispatcher is built but never driven.
    seen = {}
    real_make = bip.make_climb_dispatcher

    def spy_make(*args, parsl_config=None, executor=None,
                 force=False, **kwargs):
        seen["parsl_config"] = parsl_config
        seen["executor"] = executor
        seen["force"] = force
        return real_make(*args, parsl_config=parsl_config,
                         executor=executor, force=force, **kwargs)

    monkeypatch.setattr(bip, "make_climb_dispatcher", spy_make)

    bip.build_initial_potentials(
        manifest_path, pdb_root, data_root,
        dispatch_fn=lambda flight, executor=None, force=False: None,
        prepare_fn=lambda flight, workspace, units=None: None,
        force=True)
    assert seen["parsl_config"] is None
    assert seen["force"] is True
    # The producer built ONE executor for the whole run (local ->
    #   LocalExecutor) and threaded that same object into the climb
    #   dispatcher, so every rung rides one pool (DESIGN 6.2.11).
    assert type(seen["executor"]).__name__ == "LocalExecutor"


# ==============================================================
#  C113: submit the orchestrator as its own batch job (6.2.11)
# ==============================================================

def test_submit_orchestrator_batch_writes_and_submits(
        tmp_path, monkeypatch):
    """submit_orchestrator_batch writes the sbatch script from the
    site's orchestrator block and returns the SLURM job id sbatch
    reports; the batch command re-runs the producer WITHOUT --submit
    (structures already materialized on the login node)."""
    site = {
        "account": "rulisp-lab",
        "worker_init": ["module load imago", "source venv/bin/act"],
        "extra_scheduler_options": [],
        "orchestrator": {"cores": 2, "memory": "8G",
                         "walltime": "24:00:00"},
    }
    monkeypatch.setattr(
        bip, "load_site_config",
        lambda profile, partition=None: site)
    monkeypatch.setattr(
        bip, "resolve_choices",
        lambda s, cli: {"dispatch": "slurm-per-job",
                        "partition": "general", "nodes": 1,
                        "walltime": "02:00:00"})
    captured = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        captured["script"] = open(cmd[1]).read()
        return types.SimpleNamespace(
            stdout="Submitted batch job 12345\n")

    monkeypatch.setattr(bip.subprocess, "run", fake_run)
    args = types.SimpleNamespace(
        profile=None, dispatch="slurm-per-job", partition=None,
        nodes=None, walltime=None)
    # The FLAG vector this run parsed -- no program name in it.
    flags = ["--manifest", "m.toml", "--dispatch", "slurm-per-job",
             "--submit"]

    job_id = bip.submit_orchestrator_batch(flags, args, str(tmp_path))

    assert job_id == "12345"
    assert captured["cmd"][0] == "sbatch"
    # The script sizes the driver from the orchestrator block and
    #   re-runs the producer with --dispatch but NOT --submit.
    assert "#SBATCH --cpus-per-task=2" in captured["script"]
    assert "--dispatch slurm-per-job" in captured["script"]
    assert "--submit" not in captured["script"]
    assert "--manifest m.toml" in captured["script"]


def test_submit_applies_the_queue_override_to_the_driver_job(
        tmp_path, monkeypatch):
    """The driver's OWN batch job is sized from the overlaid site, not
    the cluster-wide defaults.  A debug queue that caps walltime must
    cap the driver's job too: asking for 24 hours on a 30-minute queue
    is rejected or silently truncated by the scheduler, with nothing
    naming the cause (DESIGN 6.2.11 -- the loader owns every overlay).
    """
    real_load = cluster_config.load_site_config

    def fake_clusterrc():
        return {
            "partitions": ["general", "debug"],
            "worker_init": ["module load imago"],
            "account": None, "walltime": "12:00:00", "nodes": 1,
            "extra_scheduler_options": [],
            "default_topology": "slurm-per-job",
            "orchestrator": {"cores": 2, "memory": "8G",
                             "walltime": "24:00:00"},
            "queue_overrides": {"debug": {
                "walltime": "00:30:00",
                "orchestrator": {"cores": 1, "memory": "2G",
                                 "walltime": "00:25:00"},
                "extra_scheduler_options": ["#SBATCH --exclusive"]}},
        }

    monkeypatch.setattr(
        cluster_config, "_load_clusterrc_module",
        lambda: types.SimpleNamespace(
            parameters_and_defaults=fake_clusterrc))
    monkeypatch.setattr(bip, "load_site_config", real_load)

    captured = {}

    def fake_run(cmd, **kwargs):
        captured["script"] = open(cmd[1]).read()
        return types.SimpleNamespace(stdout="Submitted batch job 5\n")

    monkeypatch.setattr(bip.subprocess, "run", fake_run)
    args = types.SimpleNamespace(
        profile=None, dispatch="slurm-per-job", partition="debug",
        nodes=None, walltime=None, orchestrator_cores=None,
        orchestrator_memory=None, orchestrator_walltime=None)

    bip.submit_orchestrator_batch(
        ["--partition", "debug", "--submit"], args, str(tmp_path))
    script = captured["script"]

    # The debug queue's driver shape, not the cluster-wide one.
    assert "#SBATCH --partition=debug" in script
    assert "#SBATCH --cpus-per-task=1" in script
    assert "#SBATCH --mem=2G" in script
    assert "#SBATCH --time=00:25:00" in script
    # The queue's own scheduler directive rides along.
    assert "#SBATCH --exclusive" in script
    # None of the cluster-wide defaults leak through.
    assert "24:00:00" not in script
    assert "8G" not in script


def test_submit_ignores_the_process_argv(tmp_path, monkeypatch):
    """The inner command is built from the FLAG vector this run
    parsed, never from the process's own arguments.  A library caller
    driving main(argv) lives in a process whose sys.argv belongs to
    something else entirely (a test runner, a notebook); reading it
    would submit a batch job that re-runs that instead of the
    producer.  The script also names itself by path, so the command
    is correct however this process was launched."""
    site = {
        "account": None, "worker_init": ["module load imago"],
        "extra_scheduler_options": [],
        "orchestrator": {"cores": 1, "walltime": "01:00:00"},
    }
    monkeypatch.setattr(
        bip, "load_site_config",
        lambda profile, partition=None: site)
    monkeypatch.setattr(
        bip, "resolve_choices",
        lambda s, cli: {"dispatch": "local", "partition": "p",
                        "nodes": 1, "walltime": "02:00:00"})
    # The process argv belongs to the test runner, not the producer.
    monkeypatch.setattr(
        bip.sys, "argv", ["pytest", "-q", "--some-pytest-flag"])

    captured = {}

    def fake_run(cmd, **kwargs):
        captured["script"] = open(cmd[1]).read()
        return types.SimpleNamespace(stdout="Submitted batch job 9\n")

    monkeypatch.setattr(bip.subprocess, "run", fake_run)
    args = types.SimpleNamespace(
        profile=None, dispatch="local", partition=None, nodes=None,
        walltime=None, orchestrator_cores=None,
        orchestrator_memory=None, orchestrator_walltime=None)
    flags = ["--manifest", "seed.toml", "--dispatch", "local",
             "--submit"]

    bip.submit_orchestrator_batch(flags, args, str(tmp_path))

    command = captured["script"].strip().splitlines()[-1]
    assert "--manifest seed.toml" in command
    assert "build_initial_potentials.py" in command
    assert "--submit" not in command
    # Nothing from the surrounding process leaks into the job.
    assert "pytest" not in command
    assert "--some-pytest-flag" not in command


def test_submit_honors_orchestrator_overrides(tmp_path, monkeypatch):
    """An --orchestrator-* flag overrides the site block for this
    run, key by key, and reaches the submitted script's header
    (DESIGN 6.2.11).  Without it the settings file would have to
    grow a block per orchestrator (ARCHITECTURE 9.4)."""
    site = {
        "account": None,
        "worker_init": ["module load imago"],
        "extra_scheduler_options": [],
        "orchestrator": {"cores": 2, "memory": "8G",
                         "walltime": "24:00:00"},
    }
    monkeypatch.setattr(
        bip, "load_site_config",
        lambda profile, partition=None: site)
    monkeypatch.setattr(
        bip, "resolve_choices",
        lambda s, cli: {"dispatch": "local", "partition": "general",
                        "nodes": 1, "walltime": "02:00:00"})
    captured = {}

    def fake_run(cmd, **kwargs):
        captured["script"] = open(cmd[1]).read()
        return types.SimpleNamespace(
            stdout="Submitted batch job 777\n")

    monkeypatch.setattr(bip.subprocess, "run", fake_run)
    # Under --dispatch local the driver runs the SCFs in process, so
    #   the curator raises its memory for this run instead of
    #   editing the site file.
    args = types.SimpleNamespace(
        profile=None, dispatch="local", partition=None, nodes=None,
        walltime=None, orchestrator_cores=None,
        orchestrator_memory="64G", orchestrator_walltime=None)
    flags = ["--dispatch", "local",
             "--orchestrator-memory", "64G", "--submit"]

    assert bip.submit_orchestrator_batch(
        flags, args, str(tmp_path)) == "777"
    assert "#SBATCH --mem=64G" in captured["script"]
    # The keys the flag did not name still come from the site block.
    assert "#SBATCH --cpus-per-task=2" in captured["script"]
    assert "#SBATCH --time=24:00:00" in captured["script"]


# ==================================================================
#  The adaptive k-point mesh climb (PSEUDOCODE 4e.3 / 4e.5)
#
#  These drive climb_action and converge_by_climb with a SYNTHETIC
#  dispatcher -- a stand-in for the real climb dispatcher
#  (make_climb_dispatcher) that turns each requested mesh into a Rung
#  by evaluating a per-material energy model.  Every material here is
#  a cubic cell (one axis class, equal reciprocal magnitudes), so its
#  climb is the lockstep [n,n,n] ladder and the hand-computed energies
#  below are easy to follow.
# ==================================================================

# Cubic geometry: all three axes share one class, so counts move in
#   lockstep and a density of D lands on the [round(D**(1/3))]^3 mesh.
_CUBIC_CLASSES = [0, 0, 0]
_CUBIC_RECIP_MAG = [1.0, 1.0, 1.0]


def _cubic_config(mode, flat_needed, grid_width, start_offset,
                  max_count, max_stride=8):
    """A cubic-cell ClimbConfig with a per-atom flatness threshold of
    1.0 eV and one atom per cell, so a raw energy delta counts as
    flat when it is below 1.0 / HARTREE hartree (DESIGN 7.8).  The
    stride cap defaults to the provisional eight; the unit-step tests
    that use this helper never stride, so it is inert for them."""
    return bip.ClimbConfig(
        classes=_CUBIC_CLASSES, recip_mag=_CUBIC_RECIP_MAG,
        recip_cell_volume=1.0, mode=mode, flat_needed=flat_needed,
        grid_width=grid_width, start_offset=start_offset,
        max_stride=max_stride, cell_atom_count=1, threshold=1.0,
        max_count=max_count)


def _converging_energy(mesh):
    """A cubic energy model that keeps moving up to the [5,5,5] mesh
    then goes flat: the [2..4] steps are electron-volts apart (never
    flat), while every step from [5,5,5] up shifts a tenth of a
    milli-hartree (~0.003 eV, always flat).  So the two-sided,
    two-flat-rung climb converges at [6,6,6] -- the first rung whose
    flatness is confirmed by a flat neighbour above it."""
    n = mesh[0]
    if n <= 2:
        return 0.0
    if n == 3:
        return -1.0
    if n == 4:
        return -1.5
    return -1.6 - 0.0001 * (n - 5)


def _ceiling_energy(mesh):
    """A cubic energy model that never goes flat: every rung drops a
    full hartree (~27 eV) below the last, so the climb only ever
    stops by reaching the per-axis ceiling (non-converged)."""
    return -1.0 * mesh[0]


def _flat_everywhere_energy(mesh):
    """A cubic energy model that is flat at every rung (each step a
    tenth of a milli-hartree), so a confident parallel grid finds a
    flat interior point in its very first round."""
    return -1.6 - 0.0001 * mesh[0]


def _early_plateau_energy(mesh):
    """A cubic energy that settles by [6,6,6] but whose geometric
    stride only reads flat high up.  Every step through [4,4,4] is a
    big electron-volt drop; from [5,5,5] the plateau is a tenth of a
    milli-hartree per step.  So the bracket stride [8->16] is the
    first flat one (its span sits entirely in the plateau) while
    [4->8] is not (it straddles the steep [4->5] step) -- exactly the
    live cubic case, where the bracket lands high but the convergence
    is low ([6,6,6]).  This is what makes the refine's early-stop
    matter: the rungs above [8,8,8] never need computing."""
    n = mesh[0]
    if n <= 4:
        return -31.3 + 0.5 * (4 - n)         # -29.3, -30.3, -30.8, -31.3
    return -31.4 - 0.0001 * (n - 5)          # plateau from [5,5,5] up


class _SyntheticDispatcher:
    """A synthetic climb dispatcher for the ``converge_by_climb``
    tests.  ``send(mesh_lists)`` runs each requested mesh through the
    per-material ``energy_of[material]`` model and queues the resulting
    ``Rung``; ``next_rung()`` returns queued rungs in FIFO (send)
    order, so several materials interleave deterministically.  A mesh
    named in ``fail`` (a set of ``(material, (a, b, c))`` pairs)
    queues the ``_RUN_FAILED`` marker instead, standing in for a rung
    whose run did not complete.  An optional ``send_counter`` list
    gets one entry appended per ``send`` call, so a test can count how
    many times the climb dispatched (the opening send plus one per
    continuation rung)."""

    def __init__(self, energy_of, fail=(), send_counter=None):
        self._energy_of = energy_of
        self._fail = {(material, tuple(mesh))
                      for material, mesh in fail}
        self._queue = []            # (material, result) in flight
        self._send_counter = send_counter

    def send(self, mesh_lists):
        if self._send_counter is not None:
            self._send_counter.append(len(mesh_lists))
        for material, meshes in mesh_lists.items():
            for mesh in meshes:
                if (material, tuple(mesh)) in self._fail:
                    self._queue.append((material, bip._RUN_FAILED))
                else:
                    energy = self._energy_of[material](mesh)
                    self._queue.append(
                        (material, bip.Rung(mesh, energy)))

    def next_rung(self):
        return self._queue.pop(0)


# --------------------------------------------------------------
#  climb_action -- one material's next step
# --------------------------------------------------------------

def test_climb_action_runs_next_rung_while_energy_moves():
    """With too few rungs to judge and the top below the ceiling,
    the verdict is RUN carrying the next lockstep mesh."""
    config = _cubic_config(
        mesh_climb.UNIT_STEP, flat_needed=2, grid_width=0,
        start_offset=1, max_count=20)
    rungs = [bip.Rung([2, 2, 2], 0.0), bip.Rung([3, 3, 3], -1.0)]
    action = bip.climb_action(rungs, config)
    assert action.kind == bip._ACTION_RUN
    assert action.mesh == [4, 4, 4]


def test_climb_action_reports_ceiling_when_backstop_reached():
    """A top rung at the per-axis backstop, still not flat, yields
    the CEILING verdict."""
    config = _cubic_config(
        mesh_climb.UNIT_STEP, flat_needed=2, grid_width=0,
        start_offset=1, max_count=5)
    rungs = [bip.Rung([n, n, n], _ceiling_energy([n, n, n]))
             for n in (3, 4, 5)]
    action = bip.climb_action(rungs, config)
    assert action.kind == bip._ACTION_CEILING


def test_climb_action_reports_the_converged_rung():
    """When the per-atom energy is flat across the interior, the
    verdict is CONVERGED carrying the first flat rung itself -- the
    [5,5,5] interior, not an index into the ladder."""
    config = _cubic_config(
        mesh_climb.UNIT_STEP, flat_needed=1, grid_width=0,
        start_offset=1, max_count=20)
    rungs = [bip.Rung([n, n, n], _flat_everywhere_energy([n, n, n]))
             for n in (4, 5, 6)]
    action = bip.climb_action(rungs, config)
    assert action.kind == bip._ACTION_CONVERGED
    assert action.rung is rungs[1]           # the [5,5,5] interior


# --------------------------------------------------------------
#  _sort_by_mesh / _merge_distinct
# --------------------------------------------------------------

def test_sort_by_mesh_orders_by_point_count():
    """Rungs come back ascending by full-mesh point count, whatever
    order they arrived in."""
    rungs = [bip.Rung([5, 5, 5], -1.0), bip.Rung([2, 2, 2], -3.0),
             bip.Rung([4, 4, 4], -2.0)]
    ordered = [rung.mesh for rung in bip._sort_by_mesh(rungs)]
    assert ordered == [[2, 2, 2], [4, 4, 4], [5, 5, 5]]


def test_merge_distinct_drops_a_repeated_mesh():
    """Merging keeps the ascending ladder and drops an incoming rung
    whose mesh is already present, so no zero-delta duplicate reaches
    the stop test."""
    existing = [bip.Rung([2, 2, 2], 0.0), bip.Rung([3, 3, 3], -1.0)]
    incoming = [bip.Rung([3, 3, 3], -1.0), bip.Rung([4, 4, 4], -1.5)]
    merged = [rung.mesh for rung in
              bip._merge_distinct(existing, incoming)]
    assert merged == [[2, 2, 2], [3, 3, 3], [4, 4, 4]]


# --------------------------------------------------------------
#  converge_by_climb -- the wait-for-any loop
# --------------------------------------------------------------

def test_climb_converges_over_several_rungs():
    """A cold single-rung climb starts one rung below the predicted
    [3,3,3] seed and climbs the lockstep ladder until two flat rungs
    persist, converging at [6,6,6] and recording the whole ladder."""
    materials = ["si"]
    configs = {"si": _cubic_config(
        mesh_climb.UNIT_STEP, flat_needed=2, grid_width=0,
        start_offset=1, max_count=20)}
    seed_densities = {"si": 27.0}            # seeds the [3,3,3] mesh
    dispatcher = _SyntheticDispatcher({"si": _converging_energy})

    outcomes, rungs = bip.converge_by_climb(
        materials, configs, seed_densities, dispatcher)

    assert outcomes["si"].mesh == [6, 6, 6]
    ladder = [rung.mesh for rung in rungs["si"]]
    assert ladder == [[n, n, n] for n in range(2, 9)]


def test_grid_mode_converges_in_the_opening():
    """A confident parallel grid lays [4,4,4], [5,5,5], [6,6,6] down
    together and, finding a flat interior point, converges from that
    single opening send with no continuation."""
    materials = ["si"]
    configs = {"si": _cubic_config(
        mesh_climb.PARALLEL_GRID, flat_needed=1, grid_width=1,
        start_offset=0, max_count=20)}
    seed_densities = {"si": 125.0}           # seeds the [5,5,5] mesh
    counter = []
    dispatcher = _SyntheticDispatcher(
        {"si": _flat_everywhere_energy}, send_counter=counter)

    outcomes, _ = bip.converge_by_climb(
        materials, configs, seed_densities, dispatcher)

    assert outcomes["si"].mesh == [5, 5, 5]
    assert len(counter) == 1                 # opening send only


def test_climb_stops_at_ceiling_and_tags_non_converged():
    """A climb whose energy never goes flat stops at the per-axis
    ceiling, is reported NON_CONVERGED, and fires the mismatch
    callback for that material."""
    materials = ["metal"]
    configs = {"metal": _cubic_config(
        mesh_climb.UNIT_STEP, flat_needed=2, grid_width=0,
        start_offset=1, max_count=5)}
    seed_densities = {"metal": 27.0}
    dispatcher = _SyntheticDispatcher({"metal": _ceiling_energy})
    flagged = []

    outcomes, _ = bip.converge_by_climb(
        materials, configs, seed_densities, dispatcher,
        on_non_converged=flagged.append)

    assert outcomes["metal"] is bip.NON_CONVERGED
    assert flagged == ["metal"]


def test_materials_climb_independently():
    """Two materials climb concurrently: one converges while the other
    runs to its ceiling.  Each reaches its own verdict, and the
    mismatch callback fires only for the non-converged one."""
    materials = ["good", "bad"]
    configs = {
        "good": _cubic_config(
            mesh_climb.UNIT_STEP, flat_needed=2, grid_width=0,
            start_offset=1, max_count=20),
        "bad": _cubic_config(
            mesh_climb.UNIT_STEP, flat_needed=2, grid_width=0,
            start_offset=1, max_count=5)}
    seed_densities = {"good": 27.0, "bad": 27.0}
    dispatcher = _SyntheticDispatcher(
        {"good": _converging_energy, "bad": _ceiling_energy})
    flagged = []

    outcomes, _ = bip.converge_by_climb(
        materials, configs, seed_densities, dispatcher,
        on_non_converged=flagged.append)

    assert outcomes["good"].mesh == [6, 6, 6]
    assert outcomes["bad"] is bip.NON_CONVERGED
    assert flagged == ["bad"]


# --------------------------------------------------------------
#  The bracket-refine climb (PSEUDOCODE 4e.3; DESIGN 3.12.3)
# --------------------------------------------------------------

def test_bracket_refine_converges_at_the_same_mesh_as_unit_step():
    """The default cold search -- bracket-refine, flat_needed = 2 --
    reaches the SAME converged mesh the fine unit-step climb does on
    the same energy: [6,6,6], the first rung whose flatness persists
    over two interior rungs.  This is the regression test for the
    off-by-one refine fill: filling flat_needed + 1 rungs above the
    flat stride's bottom lets the persistence test confirm a
    convergence that sits one rung into the flat region (DESIGN
    3.12.3)."""
    materials = ["si"]
    configs = {"si": _cubic_config(
        mesh_climb.BRACKET_REFINE, flat_needed=2, grid_width=0,
        start_offset=1, max_count=20)}
    seed_densities = {"si": 27.0}            # seeds the [3,3,3] mesh
    dispatcher = _SyntheticDispatcher({"si": _converging_energy})

    outcomes, _ = bip.converge_by_climb(
        materials, configs, seed_densities, dispatcher)

    assert outcomes["si"].mesh == [6, 6, 6]


def test_bracket_refine_strides_then_fills():
    """The bracket phase strides geometrically -- one, two, four ladder
    positions -- to bracket the convergence cheaply, then the refine
    phase fills the bracket one rung at a time.  Starting one rung
    below the [3,3,3] seed, it computes [2,2,2], [3,3,3], [5,5,5],
    [9,9,9] (the strides), then fills [4,4,4], [6,6,6], [7,7,7],
    [8,8,8] (DESIGN 3.12.2 / 3.12.3)."""
    materials = ["si"]
    configs = {"si": _cubic_config(
        mesh_climb.BRACKET_REFINE, flat_needed=2, grid_width=0,
        start_offset=1, max_count=20)}
    seed_densities = {"si": 27.0}
    sent = []

    class _RecordingDispatcher(_SyntheticDispatcher):
        def send(self, mesh_lists):
            for material, meshes in mesh_lists.items():
                for mesh in meshes:
                    sent.append(list(mesh))
            super().send(mesh_lists)

    dispatcher = _RecordingDispatcher({"si": _converging_energy})
    bip.converge_by_climb(materials, configs, seed_densities,
                          dispatcher)

    assert sent == [[2, 2, 2], [3, 3, 3], [5, 5, 5], [9, 9, 9],
                    [4, 4, 4], [6, 6, 6], [7, 7, 7], [8, 8, 8]]


def test_bracket_refine_stops_filling_once_converged():
    """The refine tests the consecutive block after each fill and
    stops at the first (smallest) converged rung, so it never computes
    the wide rungs above the convergence (DESIGN 3.12.3).  On an energy
    that settles by [6,6,6] but whose stride only reads flat at
    [8->16], the bracket spans up to [11,11,11], yet the fill halts at
    [7,7,7] the moment [6,6,6] is confirmed from [4..8] -- [9,10,11]
    are never run (8 calcs, not 11)."""
    materials = ["si"]
    configs = {"si": _cubic_config(
        mesh_climb.BRACKET_REFINE, flat_needed=2, grid_width=0,
        start_offset=0, max_count=20)}
    seed_densities = {"si": 1.0}             # seeds the [1,1,1] opening
    sent = []

    class _RecordingDispatcher(_SyntheticDispatcher):
        def send(self, mesh_lists):
            for material, meshes in mesh_lists.items():
                for mesh in meshes:
                    sent.append(list(mesh))
            super().send(mesh_lists)

    dispatcher = _RecordingDispatcher({"si": _early_plateau_energy})
    outcomes, _ = bip.converge_by_climb(
        materials, configs, seed_densities, dispatcher)

    assert outcomes["si"].mesh == [6, 6, 6]
    assert sent == [[1, 1, 1], [2, 2, 2], [4, 4, 4], [8, 8, 8],
                    [16, 16, 16], [5, 5, 5], [6, 6, 6], [7, 7, 7]]
    # The rungs above the convergence were never computed.
    for wide in ([9, 9, 9], [10, 10, 10], [11, 11, 11]):
        assert wide not in sent


def test_bracket_refine_stops_at_ceiling():
    """A bracket-refine climb whose energy never goes flat strides
    toward the ceiling, fills the final up-to-the-ceiling interval,
    finds no converged rung, and stops NON_CONVERGED rather than
    looping forever (DESIGN 3.12.3, the ceiling backstop)."""
    materials = ["metal"]
    configs = {"metal": _cubic_config(
        mesh_climb.BRACKET_REFINE, flat_needed=2, grid_width=0,
        start_offset=1, max_count=5)}
    seed_densities = {"metal": 27.0}
    dispatcher = _SyntheticDispatcher({"metal": _ceiling_energy})
    flagged = []

    outcomes, _ = bip.converge_by_climb(
        materials, configs, seed_densities, dispatcher,
        on_non_converged=flagged.append)

    assert outcomes["metal"] is bip.NON_CONVERGED
    assert flagged == ["metal"]


def test_bracket_refine_resumes_after_a_false_bracket():
    """A stride whose endpoints match by coincidence (an oscillating
    energy that dips and returns) reads flat, but no interior rung of
    the filled bracket passes the two-sided test.  The search does not
    accept it: with the bracket not run up to the ceiling, it resumes
    striding from the top of the bracket (DESIGN 3.12.3)."""
    config = _cubic_config(
        mesh_climb.BRACKET_REFINE, flat_needed=2, grid_width=0,
        start_offset=1, max_count=20)
    # [2,2,2] and [4,4,4] match (both 0.0), so the stride across them
    #   read flat, but the interior [3,3,3] dips far below -- no rung
    #   is truly settled.  The bracket [2,2,2]..[4,4,4] is fully
    #   filled.
    rungs = [bip.Rung([2, 2, 2], 0.0), bip.Rung([3, 3, 3], -5.0),
             bip.Rung([4, 4, 4], 0.0)]
    state = bip.BracketRefineState(
        phase=bip._PHASE_REFINE, stride=2,
        endpoints=[[2, 2, 2], [4, 4, 4]], lo=[2, 2, 2], hi=[4, 4, 4],
        from_cap=False)

    action, next_state = bip.bracket_refine_next(rungs, state, config)

    # Resumed: a fresh bracket phase striding up from the old top.
    assert action.kind == bip._ACTION_RUN
    assert action.mesh == [5, 5, 5]
    assert next_state.phase == bip._PHASE_BRACKET
    assert next_state.endpoints == [[4, 4, 4], [5, 5, 5]]


def test_climb_next_dispatches_on_mode():
    """climb_next routes a bracket-refine material through the stateful
    state machine (advancing its search state) and a unit-step material
    through the stateless climb_action (leaving the state untouched)."""
    rungs = [bip.Rung([2, 2, 2], 0.0), bip.Rung([3, 3, 3], -1.0)]

    bracket = _cubic_config(
        mesh_climb.BRACKET_REFINE, flat_needed=2, grid_width=0,
        start_offset=1, max_count=20)
    seed_state = bip.new_bracket_refine_state([2, 2, 2])
    bracket_action, advanced = bip.climb_next(
        [bip.Rung([2, 2, 2], 0.0)], seed_state, bracket)
    assert bracket_action.kind == bip._ACTION_RUN
    assert advanced.endpoints == [[2, 2, 2], [3, 3, 3]]

    unit = _cubic_config(
        mesh_climb.UNIT_STEP, flat_needed=2, grid_width=0,
        start_offset=1, max_count=20)
    unit_action, unchanged = bip.climb_next(rungs, None, unit)
    assert unit_action.kind == bip._ACTION_RUN
    assert unit_action.mesh == [4, 4, 4]     # the next lockstep rung
    assert unchanged is None                 # stateless: passed through


def test_new_search_state_is_stateful_only_for_bracket_refine():
    """Only the bracket-refine climb carries search state; the grid and
    the unit-step climb carry an empty (None) state they never read."""
    bracket = _cubic_config(
        mesh_climb.BRACKET_REFINE, flat_needed=2, grid_width=0,
        start_offset=1, max_count=20)
    assert bip.new_search_state(bracket, [[3, 3, 3]]).phase \
        == bip._PHASE_BRACKET
    for mode in (mesh_climb.UNIT_STEP, mesh_climb.PARALLEL_GRID):
        config = _cubic_config(
            mode, flat_needed=2, grid_width=0, start_offset=1,
            max_count=20)
        assert bip.new_search_state(config, [[3, 3, 3]]) is None


def test_record_converged_records_only_the_consecutive_block():
    """The recorded flatness trace is the CONSECUTIVE block around the
    converged rung, so a sparse bracket endpoint left on the ladder is
    dropped -- keeping it could make the harvest's re-judge read a
    false early convergence across the gap (DESIGN 3.12.3 / 4e.6)."""
    config = _cubic_config(
        mesh_climb.BRACKET_REFINE, flat_needed=2, grid_width=0,
        start_offset=1, max_count=20)
    # A filled block [3,3,3]..[6,6,6] plus a stray high stride endpoint
    #   [9,9,9] with a gap below it ([7,7,7], [8,8,8] absent).
    ladder = [bip.Rung([n, n, n], -1.0 - 0.001 * n)
              for n in (3, 4, 5, 6, 9)]
    converged = ladder[2]                    # the [5,5,5] rung
    record = bip.record_converged(converged, ladder, config)

    assert record["converged_mesh"] == [5, 5, 5]
    # Point-count densities of the consecutive block only (vol = 1.0);
    #   the sparse [9,9,9] (729) is dropped.
    assert record["grid_values"] == [27, 64, 125, 216]
    assert len(record["grid_energies"]) == 4


# ==================================================================
#  The climb dispatcher (PSEUDOCODE 4e.7; DESIGN 7.7)
# ==================================================================

def _mesh_of_unit(unit):
    """Read the axial-count mesh back out of a kpt-mesh unit's calc
    tag (``("kpt-mesh-4-4-4",)`` -> ``[4, 4, 4]``)."""
    token = unit.calc[0][len("kpt-mesh-"):]
    return [int(part) for part in token.split("-")]


def _mesh_reader(energy_by_mesh, resolved_by_mesh=None):
    """A fake result reader: return each mesh unit's chosen energy
    and, unless overridden, echo the requested mesh as the resolved
    one."""
    def read_fn(workspace, unit):
        mesh = _mesh_of_unit(unit)
        resolved = mesh
        if resolved_by_mesh is not None:
            resolved = resolved_by_mesh.get(tuple(mesh), mesh)
        return {"total_energy": energy_by_mesh[tuple(mesh)],
                "kpoint_mesh": resolved}
    return read_fn


def _fake_send_off(flight, units, executor, force):
    """A no-op ``send_off``: return one ``(unit, marker)`` pair per
    unit.  The marker stands in for the future and is unused -- the
    fake ``collect_next`` below decides each unit's status -- so a
    bare None suffices."""
    return [(unit, None) for unit in units]


def _fake_collect_next(status_by_mesh):
    """A fake ``collect_next`` that returns the outstanding units in
    FIFO order, tagging each with the status chosen for its mesh
    (default ``"done"``).  A non-``"done"`` status is what the
    dispatcher turns into a ``_RUN_FAILED`` rung."""
    def collect_next(flight, outstanding):
        unit, _marker = outstanding[0]
        remaining = outstanding[1:]
        mesh = tuple(_mesh_of_unit(unit))
        status = status_by_mesh.get(mesh, "done")
        entry = ReportEntry(
            id=unit.id, calc=unit.calc, status=status, detail=None,
            wingbeat_dir="/ws", runtime_seconds=None, message=None)
        return unit, entry, remaining
    return collect_next


def _make_climb(energy_by_mesh, *, status_by_mesh=None,
                resolved_by_mesh=None):
    """Build a make_climb_dispatcher with the toolchain seam mocked:
    a no-op prepare, a fake ``send_off``, a fake ``collect_next`` that
    assigns each mesh a status, and a fake result reader."""
    return bip.make_climb_dispatcher(
        {"si": object()}, {"si": {"scf_basis": "fb"}}, "/ws",
        prepare_fn=lambda flight, workspace, units=None: None,
        send_off_fn=_fake_send_off,
        collect_next_fn=_fake_collect_next(status_by_mesh or {}),
        read_fn=_mesh_reader(energy_by_mesh, resolved_by_mesh))


def test_climb_dispatcher_reads_rungs_back():
    """send launches a material's meshes; next_rung returns each as a
    ``(material, Rung)`` pair carrying the run's total energy."""
    dispatcher = _make_climb({(4, 4, 4): -1.0, (5, 5, 5): -1.1})
    dispatcher.send({"si": [[4, 4, 4], [5, 5, 5]]})
    landed = [dispatcher.next_rung(), dispatcher.next_rung()]
    assert [(material, rung.mesh, rung.energy)
            for material, rung in landed] \
        == [("si", [4, 4, 4], -1.0), ("si", [5, 5, 5], -1.1)]


def test_climb_dispatcher_marks_a_failed_mesh():
    """A unit that did not complete comes back as
    ``(material, _RUN_FAILED)``, so the climb loop reads it as a run
    failure."""
    dispatcher = _make_climb(
        {(4, 4, 4): -1.0, (5, 5, 5): -1.1},
        status_by_mesh={(5, 5, 5): "failed"})
    dispatcher.send({"si": [[4, 4, 4], [5, 5, 5]]})
    landed = [dispatcher.next_rung(), dispatcher.next_rung()]
    assert landed[0] == ("si", bip.Rung([4, 4, 4], -1.0))
    assert landed[1] == ("si", bip._RUN_FAILED)


def test_climb_dispatcher_asserts_the_mesh_is_honoured():
    """An explicit mesh must resolve to itself; a run that reports a
    different resolved mesh fails loudly rather than mis-recording
    the rung."""
    dispatcher = _make_climb(
        {(4, 4, 4): -1.0},
        resolved_by_mesh={(4, 4, 4): [9, 9, 9]})
    dispatcher.send({"si": [[4, 4, 4]]})
    with pytest.raises(RuntimeError):
        dispatcher.next_rung()


# --------------------------------------------------------------
#  converge_by_climb -- fail-fast on a rung that will not run
# --------------------------------------------------------------

def test_climb_opening_failure_is_non_converged():
    """A material whose opening rung fails to run has no rung to stand
    on and is reported NON_CONVERGED with the mismatch callback
    fired."""
    configs = {"si": _cubic_config(
        mesh_climb.UNIT_STEP, flat_needed=2, grid_width=0,
        start_offset=1, max_count=20)}
    # The cold climb's opening is the single [2,2,2] rung (seed 27 ->
    #   [3,3,3], start_offset=1 -> one rung below); fail it.
    dispatcher = _SyntheticDispatcher(
        {"si": _converging_energy}, fail={("si", (2, 2, 2))})
    flagged = []

    outcomes, _ = bip.converge_by_climb(
        ["si"], configs, {"si": 27.0}, dispatcher,
        on_non_converged=flagged.append)

    assert outcomes["si"] is bip.NON_CONVERGED
    assert flagged == ["si"]


def test_climb_missing_continuation_rung_is_non_converged():
    """When a continuation rung fails to come back, the climb cannot
    advance, so the material stops NON_CONVERGED instead of
    re-dispatching the failing mesh forever."""
    configs = {"si": _cubic_config(
        mesh_climb.UNIT_STEP, flat_needed=2, grid_width=0,
        start_offset=1, max_count=20)}
    # The opening [2,2,2] runs; the first continuation [3,3,3] fails.
    dispatcher = _SyntheticDispatcher(
        {"si": _converging_energy}, fail={("si", (3, 3, 3))})
    flagged = []

    outcomes, _ = bip.converge_by_climb(
        ["si"], configs, {"si": 27.0}, dispatcher,
        on_non_converged=flagged.append)

    assert outcomes["si"] is bip.NON_CONVERGED
    assert flagged == ["si"]


# --------------------------------------------------------------
#  record_converged -- the climb's harvest inputs (PSEUDOCODE 4e.6)
# --------------------------------------------------------------

def test_record_converged_builds_density_mesh_and_grid():
    """record_converged turns a converged rung and its ladder into
    the density / mesh / grid harvest inputs: the density is the
    full-mesh volume density (product of counts over the reciprocal
    cell volume), the mesh is stored exact, and the grid is the
    ascending ladder's densities and raw energies."""
    config = _cubic_config(                  # recip_cell_volume = 1.0
        mesh_climb.UNIT_STEP, flat_needed=2, grid_width=0,
        start_offset=1, max_count=20)
    rungs = [bip.Rung([4, 4, 4], -1.0), bip.Rung([5, 5, 5], -1.1),
             bip.Rung([6, 6, 6], -1.1)]

    out = bip.record_converged(rungs[1], rungs, config)

    assert out["converged_mesh"] == [5, 5, 5]
    assert out["converged_kpoint_density"] == 125.0      # 5^3 / 1.0
    assert out["grid_values"] == [64.0, 125.0, 216.0]
    assert out["grid_energies"] == [-1.0, -1.1, -1.1]


# --------------------------------------------------------------
#  build_climb_config: geometry sourced from the loaded cell
# --------------------------------------------------------------

_FIXTURE_STRUCTURES = os.path.join(
    os.path.dirname(__file__), "fixtures", "structures")


@pytest.mark.integration
def test_build_climb_config_axis_classes_and_recip_mag():
    """build_climb_config recomputes a cell's axis classes and
    reciprocal-axis magnitudes from its OWN space-group operations
    and lattice, before any run (DESIGN 3.12).  A cubic cell couples
    all three axes into one class with equal |b_i|; a hexagonal cell
    couples its in-plane pair (equal |b_a| = |b_b|) and leaves the c
    axis apart -- the check that would catch the reciprocal being
    read row-wise instead of down its columns."""

    thresholds = mesh_climb.DEFAULT_POLICY_THRESHOLDS
    ref = types.SimpleNamespace(kpoint_convergence_threshold=1.0e-4)

    cubic = bip.build_climb_config(
        ref, os.path.join(_FIXTURE_STRUCTURES, "si_diamond.skl"),
        confidence=0.9, under_trained=False,
        thresholds=thresholds, max_count=20)
    # All three axes share one class; the three |b_i| are equal.
    assert len(set(cubic.classes)) == 1
    assert cubic.recip_mag[0] == pytest.approx(cubic.recip_mag[1])
    assert cubic.recip_mag[1] == pytest.approx(cubic.recip_mag[2])
    # A confident, trained prediction warrants the parallel-grid
    #   mode (mesh_climb 4e.4).
    assert cubic.mode == mesh_climb.PARALLEL_GRID
    assert cubic.threshold == 1.0e-4
    assert cubic.max_count == 20

    hexagonal = bip.build_climb_config(
        ref, os.path.join(_FIXTURE_STRUCTURES, "beo_hexagonal.skl"),
        confidence=0.9, under_trained=False,
        thresholds=thresholds, max_count=20)
    # The in-plane pair a,b share a class distinct from c, and their
    #   reciprocal magnitudes are equal while c's differs.
    assert hexagonal.classes[0] == hexagonal.classes[1]
    assert hexagonal.classes[2] != hexagonal.classes[0]
    assert hexagonal.recip_mag[0] == pytest.approx(
        hexagonal.recip_mag[1])
    assert hexagonal.recip_mag[2] != pytest.approx(
        hexagonal.recip_mag[1])


@pytest.mark.integration
def test_build_climb_config_under_trained_climbs_serially():
    """An under-trained prediction (the bootstrap regime, DESIGN 7.9)
    warrants a serial climb, not a parallel grid -- the bracket-refine
    shape by default -- and carries the stride cap through to the
    config the producer's bracket phase reads."""

    ref = types.SimpleNamespace(kpoint_convergence_threshold=1.0e-4)
    config = bip.build_climb_config(
        ref, os.path.join(_FIXTURE_STRUCTURES, "si_diamond.skl"),
        confidence=0.0, under_trained=True,
        thresholds=mesh_climb.DEFAULT_POLICY_THRESHOLDS,
        max_count=20)
    assert config.mode == mesh_climb.BRACKET_REFINE
    assert config.max_stride == \
        mesh_climb.DEFAULT_POLICY_THRESHOLDS.max_stride
