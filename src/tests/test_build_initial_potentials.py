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

import pytest

from build_initial_potentials import (
    load_manifest_v2,
    CurationManifest,
    ReferenceSolid,
    ReferenceEntry,
    ManifestFingerprint,
    element_path,
    is_isolated_default_for,
    build_isolated_entry,
    list_element_dirs,
    refresh_isolated_entries,
    save_databases,
    _parse_pot_file,
    _parse_coeff_file,
)
import initial_potential_db as ipdb


pytestmark = pytest.mark.unit


# ============================================================
#  Manifest builders
# ============================================================

# A canonical valid single-solid manifest in the cod_id form.
_VALID_COD_MANIFEST = """\
schema_version = 2

[[reference_solid]]
reference_id = "au_fcc"
system_type = "crystalline"
basis = "fb"
functional = "wigner"
kpoint_integration = "linear-tetrahedral"
cod_id = 9008463
cod_revision = "2023-04-12"
kpoint_spec = { density = 60.0, shift = [0.0, 0.0, 0.0] }
scf_threshold = 1.0e-6

  [[reference_solid.entry]]
  element = "Au"
  atom_site = 1
  label = "default_solid"
  default = true
  description = "Au in fcc bulk (Fm-3m)."
"""


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

    def test_empty_manifest_is_valid(self, tmp_path):
        # A manifest with no reference solids is valid: rules
        # 2-9 are per-solid/per-entry and vacuously pass, and
        # rule 7 ranges over the (empty) set of elements seen.
        path = _write(tmp_path, "schema_version = 2\n")
        manifest = load_manifest_v2(path)
        assert manifest.reference_solids == []

    def test_structure_path_form_parses(self, tmp_path):
        # structure_path must resolve to a real file under the
        # manifest's directory.
        _write(tmp_path, "dummy structure bytes\n", name="x.skel")
        text = (
            "schema_version = 2\n\n"
            "[[reference_solid]]\n"
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


class TestRule3RequiredEntryFields:
    def test_missing_default_raises(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            "  default = true\n", ""))
        with pytest.raises(ValueError,
                           match="manifest rule 3.*default"):
            load_manifest_v2(path)

    def test_missing_atom_site_raises(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            "  atom_site = 1\n", ""))
        with pytest.raises(ValueError,
                           match="manifest rule 3.*atom_site"):
            load_manifest_v2(path)


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
            "[[reference_solid]]\n"
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


class TestRule7DefaultPerElement:
    def test_zero_defaults_raises(self, tmp_path):
        path = _write(tmp_path, _VALID_COD_MANIFEST.replace(
            "  default = true\n", "  default = false\n"))
        with pytest.raises(ValueError,
                           match="manifest rule 7.*Au has 0"):
            load_manifest_v2(path)

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

    def test_existing_curated_entry_preserved(self, tmp_path):
        # An existing file with a curated default_solid entry
        # keeps it; only the isolated entry is rebuilt, and its
        # default flag follows the manifest (false here, since
        # the manifest curates Au's default).
        root = str(tmp_path)
        _make_element(root, "au",
                      [1.0, 2.0, 3.0], [0.15, 1.5, 1.0e8])
        # Seed a valid v2 file: isolated (non-default) plus a
        # curated default_solid (default).
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
             "scf_threshold": 1e-6, "scf_iterations": 9}))
        ipdb.save(seed, element_path(root, "au"))

        dbs = refresh_isolated_entries(
            root, _manifest_curating_au(), "new", "now")
        db = dbs["au"]
        labels = sorted(e.label for e in db.potentials)
        assert labels == ["default_solid", "isolated"]
        # The curated entry survived untouched.
        assert ipdb.lookup(db, "default_solid").default is True
        # The isolated entry was rebuilt from current pot1/coeff1
        # (3 terms, not the seed's 1) and is no longer default.
        iso = ipdb.lookup(db, "isolated")
        assert iso.num_gaussians == 3
        assert iso.default is False
        assert iso.provenance["commit"] == "new"
