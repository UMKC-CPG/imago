"""test_initial_potential_db.py -- Unit tests for the
per-element potential database file-format library
(``src/scripts/initial_potential_db.py``).

The library is the *library* half of the producer / consumer
split documented in DESIGN 5.4: it knows the TOML file format
and nothing else.  These tests exercise that contract by
building synthetic byte strings and synthetic in-memory
databases, without ever invoking a real SCF run.

The library implements schema **version 2** (DESIGN 5.2): each
entry carries a ``default`` flag (rule 7) and an optional
``[[potential.fingerprint]]`` array (rules 8 and 9).  These
tests cover the v2 surface.

Coverage map
------------
* DESIGN 5.2 / PSEUDOCODE 11.1 -- the nine load-time validation
  rules each have at least one dedicated test that fires the
  rule with the expected field/label in the error message.
  Rule 7 (exactly one default), rule 8 ((method, sub_spec)
  uniqueness), and rule 9 (registered method) are the v2
  additions.
* DESIGN 5.5 / PSEUDOCODE 11.2 -- the deterministic emitter has
  tests for bit-level repeatability, %.16e float formatting,
  per-block ``=`` alignment, TOML basic-string escaping, the
  canonical provenance key order (with and without the Imago
  extras), the multi-line array layout (one value per line,
  three-space indent, trailing comma), the single-newline EOF
  rule, and the ``[[potential.fingerprint]]`` block layout
  (method, sub_spec inline table, multi-line payload vector).
* PSEUDOCODE 11.1 lookup helpers -- :func:`lookup`,
  :func:`baseline`, :func:`default_entry`, and
  :func:`find_fingerprint` are covered for both success and
  miss paths, including canonical sub-spec equality.
* Round-trip -- ``load(save(db))`` is checked to preserve every
  numerical, provenance, and fingerprint field exactly.

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
    FingerprintRecord,
    PotentialEntry,
    ElementDatabase,
    load,
    save,
    lookup,
    baseline,
    default_entry,
    find_fingerprint,
    canonicalize_sub_spec,
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
    """Return a valid 'isolated' PotentialEntry (not default).

    In the canonical two-entry fixture the curated bulk entry
    is the default, so the isolated baseline carries
    ``default = False``.
    """

    return PotentialEntry(
        label         = "isolated",
        default       = False,
        description   = "Isolated Au atom from atomSCF.",
        num_gaussians = 3,
        alpha_min     = 1.0e-3,
        alpha_max     = 1.0e+2,
        coefficients  = [1.0, 2.0, 3.0],
        alphas        = [1.0e-3, 1.0e-1, 1.0e+2],
        provenance    = _atomscf_provenance(),
    )


def _default_solid_entry() -> PotentialEntry:
    """Return a valid 'default_solid' PotentialEntry (default).

    Carries ``default = True`` so the canonical fixture
    satisfies rule 7 (exactly one default per file).
    """

    return PotentialEntry(
        label         = "default_solid",
        default       = True,
        description   = "Au in fcc bulk (Fm-3m).",
        num_gaussians = 3,
        alpha_min     = 1.0e-3,
        alpha_max     = 1.0e+2,
        coefficients  = [0.1, 0.2, 0.3],
        alphas        = [1.0e-3, 1.0e-1, 1.0e+2],
        provenance    = _imago_provenance(),
    )


def _valid_db(symbol: str = "Au") -> ElementDatabase:
    """Return a fully-populated valid schema-v2 ElementDatabase.

    Two entries: the atomSCF ``isolated`` baseline (not
    default) and the Imago-source ``default_solid`` entry
    (default).  Neither carries fingerprints, so the byte-level
    emitter fixtures stay stable; fingerprint behavior is
    exercised by the dedicated tests below.
    """

    db = ElementDatabase(
        schema_version  = 2,
        element_symbol  = symbol,
        nuclear_z       = 79.0,
        nuclear_alpha   = 4.0e-01,
        covalent_radius = 1.0,
    )
    db.potentials.append(_isolated_entry())
    db.potentials.append(_default_solid_entry())
    return db


def _bispectrum_fingerprint() -> FingerprintRecord:
    """Return a bispectrum FingerprintRecord for emitter and
    lookup tests.  The payload carries a real-valued ``values``
    vector of length ``2 * twoj2 + 1 = 5`` for the sub-spec
    ``{twoj1 = 6, twoj2 = 2}``.
    """

    return FingerprintRecord(
        method   = "bispectrum",
        sub_spec = {"twoj1": 6, "twoj2": 2},
        payload  = {"values": [0.1, 0.2, 0.3, 0.4, 0.5]},
    )


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
        assert db.schema_version  == 2
        assert db.element_symbol  == "Au"
        # nuclear_z is a real (Imago uses Z as a real number);
        # load coerces it to float regardless of on-disk spelling.
        assert db.nuclear_z       == pytest.approx(79.0)
        assert isinstance(db.nuclear_z, float)
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
    """Rule 1: schema_version must equal 2."""

    def test_wrong_schema_version_raises(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        # Substitute a non-v2 schema_version in the file.  A v1
        # value must be rejected too -- there is no on-disk
        # back-compatibility (DESIGN 5.2).
        text = open(path).read().replace(
            "schema_version  = 2\n",
            "schema_version  = 1\n")
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
            "schema_version  = 2\n"
            "element_symbol  = \"Au\"\n"
            "nuclear_z       = 79\n"
            "nuclear_alpha   = 0.4\n"
            "covalent_radius = 1.0\n"
            "\n"
            "[[potential]]\n"
            "default       = true\n"
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
            schema_version  = 2,
            element_symbol  = "Au",
            nuclear_z       = 79.0,
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
        assert "schema_version  = 2\n" in text
        assert "covalent_radius = " in text

    def test_entry_body_alignment(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        text = open(path).read()
        # Body-block width is the max over body keys plus the
        # two array keys: max(5, 11, 13, 9, 9, 12, 6) = 13.
        # "label" (5) -> 8 spaces of padding before " = ".
        assert "label         = \"isolated\"\n" in text
        # The v2 `default` flag slots right after the label and
        # aligns at the same width; isolated is not the default.
        assert "default       = false\n" in text
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


# ============================================================
#  Rule 7 -- exactly one default entry (DESIGN 5.2, v2)
# ============================================================

class TestRule7DefaultTag:
    """Rule 7: exactly one entry per file carries
    ``default = true``.  Zero and multiple are both hard
    errors -- selection must be declared explicitly, with no
    implicit fallback to the 'isolated' baseline.
    """

    def test_zero_defaults_raises(self, tmp_path):
        # Clear the one default flag so no entry is default.
        path = _path_for(tmp_path, "Au")
        db = _valid_db("Au")
        db.potentials[1].default = False
        save(db, path)
        with pytest.raises(ValueError) as excinfo:
            load(path)
        msg = str(excinfo.value)
        assert "exactly one" in msg
        assert "found 0" in msg

    def test_multiple_defaults_raises(self, tmp_path):
        # Flag both entries default; the count is now two.
        path = _path_for(tmp_path, "Au")
        db = _valid_db("Au")
        db.potentials[0].default = True
        save(db, path)
        with pytest.raises(ValueError) as excinfo:
            load(path)
        msg = str(excinfo.value)
        assert "exactly one" in msg
        assert "found 2" in msg

    def test_missing_default_field_raises(self, tmp_path):
        # `default` is a required per-entry field (rule 3): drop
        # it from the isolated entry and expect a missing-field
        # error that names the entry and the field.
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        text = open(path).read().replace(
            "default       = false\n", "", 1)
        _write_toml(path, text)
        with pytest.raises(ValueError) as excinfo:
            load(path)
        msg = str(excinfo.value)
        assert "isolated" in msg
        assert "default" in msg


# ============================================================
#  Rule 8 -- per-entry fingerprint uniqueness (DESIGN 5.2)
# ============================================================

class TestRule8FingerprintUniqueness:
    """Rule 8: within one entry, the pair (method, sub_spec)
    must be unique.  Two records with the same method *and* the
    same sub-spec keys-and-values are a hard error; the same
    method with a differing sub-spec is allowed.
    """

    def test_duplicate_method_subspec_raises(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        db = _valid_db("Au")
        # Two fingerprints with the identical (method, sub_spec).
        db.potentials[1].fingerprints = [
            FingerprintRecord(
                method   = "bispectrum",
                sub_spec = {"twoj1": 6, "twoj2": 2},
                payload  = {"values": [0.1, 0.2, 0.3, 0.4, 0.5]}),
            FingerprintRecord(
                method   = "bispectrum",
                sub_spec = {"twoj1": 6, "twoj2": 2},
                payload  = {"values": [0.6, 0.7, 0.8, 0.9, 1.0]}),
        ]
        save(db, path)
        # known_methods=None -> rule 9 skipped, so rule 8 is the
        # rule under test here.
        with pytest.raises(
                ValueError, match="duplicate fingerprint"):
            load(path)

    def test_same_method_different_subspec_ok(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        db = _valid_db("Au")
        # Same method, different sub-spec -> both records coexist.
        db.potentials[1].fingerprints = [
            FingerprintRecord(
                method   = "bispectrum",
                sub_spec = {"twoj1": 6, "twoj2": 2},
                payload  = {"values": [0.1, 0.2, 0.3, 0.4, 0.5]}),
            FingerprintRecord(
                method   = "bispectrum",
                sub_spec = {"twoj1": 8, "twoj2": 4},
                payload  = {"values": [0.1, 0.2, 0.3, 0.4, 0.5,
                                       0.6, 0.7, 0.8, 0.9]}),
        ]
        save(db, path)
        re = load(path)
        entry = lookup(re, "default_solid")
        assert len(entry.fingerprints) == 2


# ============================================================
#  Rule 9 -- fingerprint method must be registered (DESIGN 5.2)
# ============================================================

class TestRule9FingerprintMethodRegistered:
    """Rule 9: a fingerprint ``method`` must be a registered
    matcher name -- but only when the caller supplies the
    registry via ``known_methods``.  A ``None`` registry skips
    the rule, which is how isolated unit tests load files
    without importing makeinput.py's matcher table.
    """

    def _db_with_unknown_method(self):
        db = _valid_db("Au")
        db.potentials[1].fingerprints = [
            FingerprintRecord(
                method   = "nonsense",
                sub_spec = {"twoj1": 6, "twoj2": 2},
                payload  = {"values": [0.1, 0.2, 0.3, 0.4, 0.5]}),
        ]
        return db

    def test_unknown_method_raises_with_registry(
            self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(self._db_with_unknown_method(), path)
        with pytest.raises(ValueError) as excinfo:
            load(path, known_methods={"bispectrum", "reduce"})
        msg = str(excinfo.value)
        assert "nonsense" in msg
        assert "registered matcher" in msg

    def test_unknown_method_skipped_without_registry(
            self, tmp_path):
        # With no registry the method string is not policed, so
        # the file loads cleanly.
        path = _path_for(tmp_path, "Au")
        save(self._db_with_unknown_method(), path)
        re = load(path)
        entry = lookup(re, "default_solid")
        assert entry.fingerprints[0].method == "nonsense"

    def test_known_method_accepted_with_registry(
            self, tmp_path):
        path = _path_for(tmp_path, "Au")
        db = _valid_db("Au")
        db.potentials[1].fingerprints = [_bispectrum_fingerprint()]
        save(db, path)
        re = load(path, known_methods={"bispectrum", "reduce"})
        entry = lookup(re, "default_solid")
        assert entry.fingerprints[0].method == "bispectrum"


# ============================================================
#  default_entry()
# ============================================================

class TestDefaultEntry:
    """default_entry returns the single default-tagged entry,
    which is distinct from baseline() unless the curator marks
    the isolated entry itself as the default.
    """

    def test_returns_default_tagged_entry(self):
        db = _valid_db("Au")
        assert default_entry(db).label == "default_solid"

    def test_distinct_from_baseline_when_different(self):
        db = _valid_db("Au")
        # In the canonical fixture the bulk entry is default and
        # the isolated entry is the baseline -- different objects.
        assert default_entry(db) is not baseline(db)
        assert default_entry(db).label != baseline(db).label

    def test_same_as_baseline_when_isolated_is_default(self):
        # A curator may flag the isolated baseline as default;
        # then the two helpers return the same object.
        db = _valid_db("Au")
        db.potentials[0].default = True   # isolated
        db.potentials[1].default = False  # default_solid
        assert default_entry(db) is baseline(db)


# ============================================================
#  find_fingerprint() and canonicalize_sub_spec()
# ============================================================

class TestFindFingerprint:
    """find_fingerprint matches on method plus canonical
    sub-spec equality, and raises KeyError when no record
    matches.
    """

    def _entry_with_fp(self):
        entry = _default_solid_entry()
        entry.fingerprints = [_bispectrum_fingerprint()]
        return entry

    def test_finds_matching_record(self):
        entry = self._entry_with_fp()
        fp = find_fingerprint(
            entry, "bispectrum", {"twoj1": 6, "twoj2": 2})
        assert fp.payload["values"][0] == pytest.approx(0.1)

    def test_match_ignores_subspec_key_order(self):
        entry = self._entry_with_fp()
        # Same keys, reversed insertion order -> still a match.
        fp = find_fingerprint(
            entry, "bispectrum", {"twoj2": 2, "twoj1": 6})
        assert fp.method == "bispectrum"

    def test_match_treats_int_and_float_equal(self):
        entry = self._entry_with_fp()
        # 6 == 6.0 and 2 == 2.0 under canonical equality.
        fp = find_fingerprint(
            entry, "bispectrum", {"twoj1": 6.0, "twoj2": 2.0})
        assert fp.method == "bispectrum"

    def test_missing_method_raises_keyerror(self):
        entry = self._entry_with_fp()
        with pytest.raises(KeyError):
            find_fingerprint(
                entry, "reduce", {"twoj1": 6, "twoj2": 2})

    def test_missing_subspec_raises_keyerror(self):
        entry = self._entry_with_fp()
        with pytest.raises(KeyError):
            find_fingerprint(
                entry, "bispectrum", {"twoj1": 8, "twoj2": 8})


class TestCanonicalizeSubSpec:
    """The canonical form underpins both rule-8 uniqueness and
    find_fingerprint equality, so its key properties are pinned
    directly: key-order independence, int/float normalization,
    bool kept distinct from the integers it subclasses, and
    hashability.
    """

    def test_key_order_independent(self):
        a = canonicalize_sub_spec({"twoj1": 6, "twoj2": 2})
        b = canonicalize_sub_spec({"twoj2": 2, "twoj1": 6})
        assert a == b

    def test_int_and_float_equal(self):
        a = canonicalize_sub_spec({"twoj1": 6})
        b = canonicalize_sub_spec({"twoj1": 6.0})
        assert a == b

    def test_bool_distinct_from_int(self):
        # bool is a subclass of int; the canonical form must not
        # collapse True onto the integer 1.
        a = canonicalize_sub_spec({"flag": True})
        b = canonicalize_sub_spec({"flag": 1})
        assert a != b

    def test_string_distinct_from_number(self):
        a = canonicalize_sub_spec({"x": 1.0})
        b = canonicalize_sub_spec({"x": "1.0"})
        assert a != b

    def test_result_is_hashable(self):
        # Must be usable as a set member (rule-8 dedup relies on
        # this) and as a nested value.
        canon = canonicalize_sub_spec(
            {"twoj1": 6, "nested": {"a": 1}, "vec": [1.0, 2.0]})
        assert len({canon}) == 1


# ============================================================
#  Emitter -- fingerprint sub-blocks (DESIGN 5.5, v2)
# ============================================================

class TestEmitterFingerprintBlock:
    """The emitter writes [[potential.fingerprint]] blocks with
    method and sub_spec first, the sub_spec as a deterministic
    inline table, and a real-valued payload vector laid out as
    a multi-line array.  Entries with no fingerprints emit no
    such block.
    """

    def _db_with_fingerprint(self):
        db = _valid_db("Au")
        db.potentials[1].fingerprints = [_bispectrum_fingerprint()]
        return db

    def test_no_block_when_no_fingerprints(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(_valid_db("Au"), path)
        text = open(path).read()
        assert "[[potential.fingerprint]]" not in text

    def test_block_header_and_fields(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(self._db_with_fingerprint(), path)
        text = open(path).read()
        assert "[[potential.fingerprint]]\n" in text
        # Block width spans method/sub_spec/values; "sub_spec"
        # (8 chars) is the widest, so "method" (6) gets two
        # spaces of padding before " = ".
        assert "method   = \"bispectrum\"\n" in text
        # sub_spec is a deterministic inline table with keys in
        # sorted order and one space inside the braces.
        assert ("sub_spec = { twoj1 = 6, twoj2 = 2 }\n"
                in text)
        # The payload vector is a multi-line %.16e array.  0.1
        # is not exactly representable, so it round-trips as
        # ...0001; 0.5 is exact and ends in zeros.
        assert "values   = [\n" in text
        assert "   1.0000000000000001e-01,\n" in text
        assert "   5.0000000000000000e-01,\n" in text

    def test_round_trip_preserves_fingerprint(self, tmp_path):
        path = _path_for(tmp_path, "Au")
        save(self._db_with_fingerprint(), path)
        re = load(path)
        entry = lookup(re, "default_solid")
        assert len(entry.fingerprints) == 1
        fp = find_fingerprint(
            entry, "bispectrum", {"twoj1": 6, "twoj2": 2})
        # sub_spec ints survive as ints; payload floats survive
        # bit-exact through the %.16e round trip.
        assert fp.sub_spec == {"twoj1": 6, "twoj2": 2}
        assert fp.payload["values"] == [0.1, 0.2, 0.3, 0.4, 0.5]

    def test_fingerprint_emit_is_deterministic(self, tmp_path):
        # Two saves of the same fingerprinted database produce
        # byte-identical files, including the inline-table key
        # ordering and the payload array layout.
        path_a = _path_for(tmp_path / "a", "Au")
        path_b = _path_for(tmp_path / "b", "Au")
        save(self._db_with_fingerprint(), path_a)
        save(self._db_with_fingerprint(), path_b)
        with open(path_a, "rb") as h_a, open(path_b, "rb") as h_b:
            assert h_a.read() == h_b.read()
