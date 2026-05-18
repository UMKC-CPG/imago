"""test_initial_potential_db.py -- Unit tests for the
per-element potential database file-format library
(``src/scripts/initial_potential_db.py``).

The library is the *library* half of the producer / consumer
split documented in DESIGN 5.4: it knows the TOML file format
and nothing else.  These tests exercise that contract by
building synthetic byte strings and synthetic in-memory
databases, without ever invoking a real SCF run.

Coverage map
------------
* DESIGN 5.2 / PSEUDOCODE 11.1 -- the six load-time validation
  rules each have at least one dedicated test that fires the
  rule with the expected field/label in the error message.
* DESIGN 5.5 / PSEUDOCODE 11.2 -- the deterministic emitter has
  tests for bit-level repeatability, %.16e float formatting,
  per-block ``=`` alignment, TOML basic-string escaping, the
  canonical provenance key order (with and without the Imago
  extras), the multi-line array layout (one value per line,
  three-space indent, trailing comma), and the single-newline
  EOF rule.
* PSEUDOCODE 11.1 lookup helpers -- :func:`lookup` and
  :func:`baseline` are covered for both success and KeyError
  paths.
* Round-trip -- ``load(save(db))`` is checked to preserve every
  numerical and provenance field exactly.

Why a separate test file (not folded into one of the existing
``test_structure_control_*`` modules): ``initial_potential_db``
is a distinct script-layer module unrelated to StructureControl
or the atomic-structure manipulation surface, so it deserves
its own test surface.  conftest.py's ``SCRIPTS_DIR`` insertion
lets us import the module directly with
``from initial_potential_db import ...``.
"""

import os

import pytest

from initial_potential_db import (
    PotentialEntry,
    ElementDatabase,
    load,
    save,
    lookup,
    baseline,
    require_provenance,
)


# All tests here are pure-computation / temp-file unit tests --
# no fixtures from the element database, no Fortran binaries.
pytestmark = pytest.mark.unit


# ============================================================
#  Test helpers
# ============================================================

def _atomscf_provenance() -> dict:
    """Return a valid atomSCF-source provenance dict."""

    return {
        "source": "atomSCF",
        "commit": "abcdef1",
        "generated_at": "2026-05-18T14:00:00Z",
    }


def _imago_provenance() -> dict:
    """Return a valid Imago-source provenance dict.

    Carries the extra reference-run fields required by
    DESIGN 5.2 (and enforced by require_provenance) when
    source == "Imago".
    """

    return {
        "source": "Imago",
        "commit": "fedcba2",
        "generated_at": "2026-05-18T14:30:00Z",
        "reference_id": "COD-1011098",
        "atom_site": 1,
        "kpoint_spec": "12 12 12 0 0 0",
        "convergence_threshold": 1.0e-6,
        "scf_iterations": 28,
    }


def _isolated_entry() -> PotentialEntry:
    """Return a valid 'isolated' PotentialEntry."""

    return PotentialEntry(
        label         = "isolated",
        description   = "Isolated Au atom from atomSCF.",
        num_gaussians = 3,
        alpha_min     = 1.0e-3,
        alpha_max     = 1.0e+2,
        coefficients  = [1.0, 2.0, 3.0],
        alphas        = [1.0e-3, 1.0e-1, 1.0e+2],
        provenance    = _atomscf_provenance(),
    )


def _default_solid_entry() -> PotentialEntry:
    """Return a valid 'default_solid' PotentialEntry."""

    return PotentialEntry(
        label         = "default_solid",
        description   = "Au in fcc bulk (Fm-3m).",
        num_gaussians = 3,
        alpha_min     = 1.0e-3,
        alpha_max     = 1.0e+2,
        coefficients  = [0.1, 0.2, 0.3],
        alphas        = [1.0e-3, 1.0e-1, 1.0e+2],
        provenance    = _imago_provenance(),
    )


def _valid_db(symbol: str = "Au") -> ElementDatabase:
    """Return a fully-populated valid ElementDatabase."""

    db = ElementDatabase(
        schema_version  = 1,
        element_symbol  = symbol,
        nuclear_z       = 79,
        nuclear_alpha   = 4.0e-01,
        covalent_radius = 1.0,
    )
    db.potentials.append(_isolated_entry())
    db.potentials.append(_default_solid_entry())
    return db


def _path_for(tmp_path, elem: str) -> str:
    """Build the canonical ``<elem>/s_gaussian_pot.toml`` path
    underneath the per-test tmp_path so that rule 2's parent-
    directory check has the right neighborhood to compare
    against.
    """

    elem_dir = tmp_path / elem.lower()
    elem_dir.mkdir(parents=True, exist_ok=True)
    return str(elem_dir / "s_gaussian_pot.toml")


def _write_toml(path: str, content: str) -> None:
    """Write ``content`` to ``path`` verbatim (UTF-8, LF EOL)."""

    with open(path, "w", encoding="utf-8", newline="\n") as h:
        h.write(content)


# ============================================================
#  Reader -- happy path
# ============================================================

class TestLoadHappyPath:
    """The reader accepts a fully-valid TOML file and produces
    an ElementDatabase whose dataclass fields exactly mirror
    the file contents.
    """

    def test_round_trips_top_level_fields(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        db = load(path)
        assert db.schema_version  == 1
        assert db.element_symbol  == "Au"
        assert db.nuclear_z       == 79
        assert db.nuclear_alpha   == pytest.approx(0.4)
        assert db.covalent_radius == pytest.approx(1.0)

    def test_round_trips_entry_count_and_labels(
            self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        db = load(path)
        labels = [entry.label for entry in db.potentials]
        assert labels == ["isolated", "default_solid"]

    def test_element_symbol_case_insensitive_dir(
            self, tmp_path):
        """Rule 2 compares case-insensitively, so a directory
        named ``au`` paired with element_symbol ``"Au"`` is
        accepted -- the lowercase-element-directory convention
        of share/atomicPDB/ must not force the file content
        to be lowercase as well.
        """

        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        db = load(path)
        assert db.element_symbol == "Au"


# ============================================================
#  Reader -- validation rule firings (DESIGN 5.2)
# ============================================================

class TestRule1SchemaVersion:
    """Rule 1: schema_version must equal 1."""

    def test_wrong_schema_version_raises(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        # Substitute a different schema_version in the file.
        text = open(path).read().replace(
            "schema_version  = 1\n",
            "schema_version  = 2\n")
        _write_toml(path, text)
        with pytest.raises(
                ValueError, match="unsupported schema_version"):
            load(path)


class TestRule2ElementSymbolMatchesDir:
    """Rule 2: element_symbol matches the parent directory
    name (case-insensitive).
    """

    def test_mismatch_raises_naming_both_sides(self, tmp_path):
        # Build the file inside a directory named "ag" but
        # declare element_symbol "Au" -- a mismatch.
        path = _path_for(tmp_path, "Ag")
        db = _valid_db("Au")
        save(db, path)
        with pytest.raises(ValueError) as excinfo:
            load(path)
        msg = str(excinfo.value)
        assert "element_symbol" in msg
        # Both the value and the parent directory name appear
        # in the message so the user can fix either side.
        assert "Au" in msg
        assert "ag" in msg


class TestRule3RequiredFieldsPresent:
    """Rule 3: every required field present at top level and
    inside every [[potential]] block.
    """

    def test_missing_top_level_field_raises(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        # Drop nuclear_z from the file.
        lines = open(path).read().splitlines()
        lines = [ln for ln in lines
                 if not ln.startswith("nuclear_z")]
        _write_toml(path, "\n".join(lines) + "\n")
        with pytest.raises(
                ValueError,
                match="missing top-level field: nuclear_z"):
            load(path)

    def test_missing_entry_label_raises(self, tmp_path):
        # Hand-craft a minimal file with a [[potential]] block
        # that omits "label".  The error message names the
        # missing field explicitly.
        path = _path_for(tmp_path, "Au")
        text = (
            "schema_version  = 1\n"
            "element_symbol  = \"Au\"\n"
            "nuclear_z       = 79\n"
            "nuclear_alpha   = 0.4\n"
            "covalent_radius = 1.0\n"
            "\n"
            "[[potential]]\n"
            "description   = \"x\"\n"
            "num_gaussians = 1\n"
            "alpha_min     = 1.0e-3\n"
            "alpha_max     = 1.0e+2\n"
            "coefficients  = [1.0]\n"
            "alphas        = [1.0]\n"
            "\n"
            "[potential.provenance]\n"
            "source       = \"atomSCF\"\n"
            "commit       = \"x\"\n"
            "generated_at = \"x\"\n"
        )
        _write_toml(path, text)
        with pytest.raises(
                ValueError, match="missing required field: label"):
            load(path)

    def test_missing_entry_description_raises(self, tmp_path):
        # Save a valid db, then strip the description line from
        # the 'isolated' entry.
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        text = open(path).read().replace(
            "description   = \"Isolated Au atom from "
            "atomSCF.\"\n", "", 1)
        _write_toml(path, text)
        with pytest.raises(ValueError) as excinfo:
            load(path)
        msg = str(excinfo.value)
        assert "isolated" in msg
        assert "description" in msg

    def test_missing_provenance_field_raises(self, tmp_path):
        # The atomSCF provenance must carry source/commit/
        # generated_at.  Drop "commit".
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        text = open(path).read().replace(
            "commit       = \"abcdef1\"\n", "", 1)
        _write_toml(path, text)
        with pytest.raises(ValueError) as excinfo:
            load(path)
        msg = str(excinfo.value)
        assert "isolated" in msg
        assert "commit" in msg

    def test_missing_imago_extra_field_raises(self, tmp_path):
        # The Imago-source entry requires reference_id and
        # friends.  Drop scf_iterations from the default_solid
        # provenance and expect a clear error naming the entry.
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        text = open(path).read().replace(
            "scf_iterations        = 28\n", "", 1)
        _write_toml(path, text)
        with pytest.raises(ValueError) as excinfo:
            load(path)
        msg = str(excinfo.value)
        assert "default_solid" in msg
        assert "Imago provenance" in msg
        assert "scf_iterations" in msg

    def test_invalid_provenance_source_raises(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        text = open(path).read().replace(
            "source       = \"atomSCF\"",
            "source       = \"oracle\"")
        _write_toml(path, text)
        with pytest.raises(
                ValueError,
                match="provenance.source must be"):
            load(path)


class TestRule4LengthConsistency:
    """Rule 4: ``len(coefficients) == len(alphas) ==
    num_gaussians`` for every entry.
    """

    def test_coefficients_too_short_raises(self, tmp_path):
        # Build a db whose 'isolated' entry declares 3
        # Gaussians but only has 2 coefficients.
        path = _path_for(tmp_path, "Au")
        db = _valid_db("Au")
        db.potentials[0].coefficients = [1.0, 2.0]
        save(db, path)
        with pytest.raises(ValueError) as excinfo:
            load(path)
        msg = str(excinfo.value)
        assert "isolated" in msg
        assert "coefficients/alphas length" in msg
        assert "num_gaussians" in msg

    def test_alphas_too_short_raises(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        db = _valid_db("Au")
        db.potentials[0].alphas = [1.0, 2.0]
        save(db, path)
        with pytest.raises(ValueError) as excinfo:
            load(path)
        assert "coefficients/alphas length" in str(
            excinfo.value)


class TestRule5LabelUniqueness:
    """Rule 5: labels are unique within the file."""

    def test_duplicate_label_raises(self, tmp_path):
        # Two [[potential]] entries with the same label.
        path = _path_for(tmp_path, "Au")
        db = _valid_db("Au")
        # Replace default_solid with a second "isolated"
        # entry that is otherwise valid.
        db.potentials[1].label = "isolated"
        # Provenance can stay as Imago because the duplicate-
        # label check fires before rule 6.
        save(db, path)
        with pytest.raises(
                ValueError,
                match="duplicate \\[\\[potential\\]\\] label"):
            load(path)


class TestRule6IsolatedBaselinePresent:
    """Rule 6: every database must carry a label='isolated'
    entry so the validation harness can reach it via baseline().
    """

    def test_missing_isolated_raises(self, tmp_path):
        # Build a file that has only the 'default_solid' entry.
        path = _path_for(tmp_path, "Au")
        db = ElementDatabase(
            schema_version  = 1,
            element_symbol  = "Au",
            nuclear_z       = 79,
            nuclear_alpha   = 4.0e-01,
            covalent_radius = 1.0,
        )
        db.potentials.append(_default_solid_entry())
        save(db, path)
        with pytest.raises(
                ValueError,
                match="missing required 'isolated' baseline"):
            load(path)


# ============================================================
#  require_provenance() -- isolated unit tests
# ============================================================

class TestRequireProvenance:
    """The provenance validator is exposed as part of the
    public surface (PSEUDOCODE 11.1).  Cover both its happy
    paths and each failure mode with direct calls, not via
    load(), so the rule firings are pinned to this helper.
    """

    def test_accepts_valid_atomscf(self):
        require_provenance(
            _atomscf_provenance(), "x.toml", "isolated")

    def test_accepts_valid_imago(self):
        require_provenance(
            _imago_provenance(), "x.toml", "default_solid")

    def test_rejects_missing_base_field(self):
        prov = _atomscf_provenance()
        del prov["commit"]
        with pytest.raises(
                ValueError, match="commit"):
            require_provenance(prov, "x.toml", "isolated")

    def test_rejects_unknown_source(self):
        prov = _atomscf_provenance()
        prov["source"] = "oracle"
        with pytest.raises(
                ValueError, match="provenance.source"):
            require_provenance(prov, "x.toml", "isolated")

    def test_rejects_imago_without_extras(self):
        prov = _imago_provenance()
        del prov["reference_id"]
        with pytest.raises(
                ValueError, match="reference_id"):
            require_provenance(
                prov, "x.toml", "default_solid")


# ============================================================
#  lookup() and baseline()
# ============================================================

class TestLookup:
    """The lookup helpers do exactly what their names imply,
    and both produce informative error messages on miss.
    """

    def test_lookup_returns_matching_entry(self):
        db = _valid_db("Au")
        entry = lookup(db, "default_solid")
        assert entry.label == "default_solid"
        assert entry.description == "Au in fcc bulk (Fm-3m)."

    def test_lookup_missing_label_raises_keyerror(self):
        db = _valid_db("Au")
        with pytest.raises(KeyError, match="no_such_label"):
            lookup(db, "no_such_label")

    def test_baseline_returns_isolated_entry(self):
        db = _valid_db("Au")
        assert baseline(db).label == "isolated"


# ============================================================
#  Emitter -- bit-level determinism (DESIGN 5.5)
# ============================================================

class TestEmitterDeterminism:
    """The emitter is a pure function of the in-memory
    database: given the same input, the same bytes come out.
    This is the bit-level guarantee at the emitter level of
    the layered-reproducibility contract (DESIGN 5.5).
    """

    def test_two_saves_byte_identical(self, tmp_path):
        # Use two element subdirectories so the parent-dir
        # rule passes for both, then assert byte-equality.
        path_a = _path_for(tmp_path / "a", "Au")
        path_b = _path_for(tmp_path / "b", "Au")
        save(_valid_db("Au"), path_a)
        save(_valid_db("Au"), path_b)
        with open(path_a, "rb") as h_a, open(path_b, "rb") as h_b:
            assert h_a.read() == h_b.read()

    def test_provenance_order_independent_of_dict_order(
            self, tmp_path):
        """The emitter sorts provenance keys into a canonical
        order, so two databases that differ only in the
        insertion order of provenance keys must produce
        byte-identical output.
        """

        db_a = _valid_db("Au")
        db_b = _valid_db("Au")
        # Reverse the insertion order of the provenance dict
        # for db_b's Imago entry.
        prov = db_b.potentials[1].provenance
        reordered = {k: prov[k] for k in reversed(list(prov))}
        db_b.potentials[1].provenance = reordered

        path_a = _path_for(tmp_path / "a", "Au")
        path_b = _path_for(tmp_path / "b", "Au")
        save(db_a, path_a)
        save(db_b, path_b)
        with open(path_a, "rb") as h_a, open(path_b, "rb") as h_b:
            assert h_a.read() == h_b.read()


# ============================================================
#  Emitter -- format details (DESIGN 5.5)
# ============================================================

class TestEmitterFormat:
    """The emitted file matches DESIGN 5.5's layout: per-block
    ``=`` alignment, %.16e floats, multi-line arrays with
    three-space indent and trailing commas, blank-line
    separators, single trailing newline, and TOML basic-string
    escapes for strings that need them.
    """

    def test_top_block_alignment(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        text = open(path).read()
        # Top-block width is set by the longest top-level key
        # name, "covalent_radius" (15 chars).  "schema_version"
        # is 14 chars so it gets one space of padding before
        # the " = ".
        assert "schema_version  = 1\n" in text
        assert "covalent_radius = " in text

    def test_entry_body_alignment(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        text = open(path).read()
        # Body-block width is the max over body keys plus the
        # two array keys: max(5, 11, 13, 9, 9, 12, 6) = 13.
        # "label" (5) -> 8 spaces of padding before " = ".
        assert "label         = \"isolated\"\n" in text
        # "num_gaussians" (13) -> no padding before " = ".
        assert "num_gaussians = 3\n" in text
        # The array openers align with the rest of the body.
        assert "coefficients  = [\n" in text
        assert "alphas        = [\n" in text

    def test_provenance_alignment_grows_for_imago(
            self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        text = open(path).read()
        # The default_solid provenance carries
        # "convergence_threshold" (21 chars), which sets the
        # alignment width for that block.  "source" (6 chars)
        # then needs 15 spaces of padding.
        assert ("source                = \"Imago\"\n"
                in text)
        # The atomSCF provenance has only base keys, so the
        # alignment width is just "generated_at" (12 chars).
        assert "source       = \"atomSCF\"\n" in text

    def test_floats_use_pct16e(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        text = open(path).read()
        # 1.0e-3 -> "1.0000000000000000e-03" (17 significant
        # digits, one before the decimal point and 16 after).
        assert "1.0000000000000000e-03" in text
        # 1.0e+2 -> "1.0000000000000000e+02".
        assert "1.0000000000000000e+02" in text

    def test_multiline_array_layout(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        text = open(path).read()
        # Each array element is on its own line, indented by
        # three spaces, with a trailing comma.
        assert "   1.0000000000000000e+00,\n" in text
        assert "   2.0000000000000000e+00,\n" in text
        assert "   3.0000000000000000e+00,\n" in text
        # The closing bracket is on its own line, no indent.
        assert "\n]\n" in text

    def test_single_trailing_newline(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        with open(path, "rb") as h:
            data = h.read()
        # Last byte is a single LF; the byte before it is not
        # another LF (no blank line at EOF).
        assert data.endswith(b"\n")
        assert not data.endswith(b"\n\n")

    def test_string_escapes(self, tmp_path):
        """TOML basic strings escape backslash, double quote,
        and the named control characters.  Verify by writing
        a description that contains all three and reading
        back the literal bytes.
        """

        path = _path_for(tmp_path, "Au")
        db = _valid_db("Au")
        db.potentials[0].description = "a\"b\\c\nd\te"
        save(db, path)
        text = open(path).read()
        assert "\"a\\\"b\\\\c\\nd\\te\"" in text

    def test_unicode_passes_through(self, tmp_path):
        """Non-control Unicode characters are emitted as-is
        (the emitter only escapes control chars, backslash,
        and double quote).  Verify with a representative
        non-ASCII description.
        """

        path = _path_for(tmp_path, "Au")
        db = _valid_db("Au")
        db.potentials[0].description = "Au (fcc, lattice ≈4Å)"
        save(db, path)
        text = open(path).read()
        assert "Au (fcc, lattice ≈4Å)" in text
        # Reload to confirm tomllib accepts the bytes we wrote.
        reloaded = load(path)
        assert (lookup(reloaded, "isolated").description
                == "Au (fcc, lattice ≈4Å)")


# ============================================================
#  Round-trip (DESIGN 5.4 / PSEUDOCODE 11.1+11.2)
# ============================================================

class TestRoundTrip:
    """load(save(db)) preserves every numerical and provenance
    field exactly.  The bit-deterministic emitter + %.16e
    float format together guarantee that the round trip is
    lossless for IEEE-754 binary64 values, so this is checked
    with strict equality rather than ``approx``.
    """

    def test_round_trip_top_level(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        db = _valid_db("Au")
        save(db, path)
        re = load(path)
        assert re.schema_version  == db.schema_version
        assert re.element_symbol  == db.element_symbol
        assert re.nuclear_z       == db.nuclear_z
        assert re.nuclear_alpha   == db.nuclear_alpha
        assert re.covalent_radius == db.covalent_radius

    def test_round_trip_entry_arrays_bit_exact(
            self, tmp_path):
        """Numerical arrays survive the text round trip with
        bit-exact equality.  Use awkward IEEE-754 values
        (denormals, exponent boundaries, exact integer in
        float) to confirm that %.16e is genuinely
        round-trip-safe and that tomllib parses our emitted
        bytes back to the same binary doubles.
        """

        awkward = [
            1.0,
            0.1,
            1.0 / 3.0,
            2.2250738585072014e-308,
            1.7976931348623157e+308,
            -1.5e-7,
        ]
        path = _path_for(tmp_path, "Au")
        db = _valid_db("Au")
        db.potentials[0].num_gaussians = len(awkward)
        db.potentials[0].coefficients = list(awkward)
        db.potentials[0].alphas       = list(awkward)
        # Match the other entry to the same length so it
        # remains valid (rule 4).
        db.potentials[1].num_gaussians = len(awkward)
        db.potentials[1].coefficients = list(awkward)
        db.potentials[1].alphas       = list(awkward)

        save(db, path)
        re = load(path)
        baseline_entry = baseline(re)
        for original, recovered in zip(
                awkward, baseline_entry.coefficients):
            assert recovered == original
        for original, recovered in zip(
                awkward, baseline_entry.alphas):
            assert recovered == original

    def test_round_trip_provenance_full(self, tmp_path):
        """Provenance fields survive the round trip in full,
        including the Imago-source extras.  Cross-check both
        the type and the value so an accidental int->float
        slip in either save or load surfaces.
        """

        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        re = load(path)
        prov_iso = lookup(re, "isolated").provenance
        assert prov_iso == _atomscf_provenance()
        prov_def = lookup(re, "default_solid").provenance
        assert prov_def == _imago_provenance()
        # atom_site and scf_iterations stay as ints, not
        # floats.
        assert isinstance(prov_def["atom_site"], int)
        assert isinstance(
            prov_def["scf_iterations"], int)
        # convergence_threshold stays as a float.
        assert isinstance(
            prov_def["convergence_threshold"], float)

    def test_round_trip_idempotent_save(self, tmp_path):
        """save(load(save(db)), path2) produces byte-identical
        output to the first save.  This is the strongest
        idempotency claim the emitter makes: writing what was
        read changes nothing.
        """

        path_a = _path_for(tmp_path / "a", "Au")
        path_b = _path_for(tmp_path / "b", "Au")
        save(_valid_db("Au"), path_a)
        re = load(path_a)
        save(re, path_b)
        with open(path_a, "rb") as h_a, open(path_b, "rb") as h_b:
            assert h_a.read() == h_b.read()
