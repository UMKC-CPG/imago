"""Tests for cod_fish.py -- the COD acquisition front-end (ARCHITECTURE
9.5).

These cover the offline logic: CSV parsing, advisory ranking, the
index/id resolution that lets a student pin a row by its small table
number, and the strict revision check (with the HTTP layer stubbed so no
test hits the network).
"""

import pytest

import cod_fish
from cod_fish import CodFishError


# A trimmed COD result.php CSV -- a leading comment line (skipped) then a
# header and two single-element silicon rows: a named, ambient, common
# spacegroup phase, and an unnamed high-pressure one in a different group.
_CSV = (
    "# Search results from the Crystallography Open Database\n"
    "file,a,b,c,alpha,beta,gamma,vol,celltemp,diffrtemp,"
    "cellpressure,diffrpressure,sg,sgNumber,commonname,chemname,"
    "mineral,formula,calcformula\n"
    "2104737,5.43,5.43,5.43,90,90,90,160.1,,,,,F d -3 m,227,,,"
    "Silicon,- Si -,Si\n"
    "9000001,6.62,6.62,6.62,90,90,90,290.3,,,11000000,,I a -3,206,"
    ",,,- Si -,Si\n")


class TestParseCsv:
    def test_extracts_rows_and_fields(self):
        rows = cod_fish._parse_cod_csv(_CSV)
        assert [r["id"] for r in rows] == ["2104737", "9000001"]
        assert rows[0]["sgnum"] == "227"
        assert rows[0]["name"] == "Silicon"     # mineral column
        assert rows[0]["formula"] == "- Si -"

    def test_skips_comment_lines(self):
        # The leading '#' line must not become a spurious row.
        assert len(cod_fish._parse_cod_csv(_CSV)) == 2


class TestRank:
    def test_named_ambient_consensus_outranks_oddball(self):
        rows = cod_fish._parse_cod_csv(_CSV)
        ordered = cod_fish.rank(rows)
        best, best_reasons = ordered[0]
        assert best["id"] == "2104737"          # named + ambient
        assert "ambient" in best_reasons
        assert any("named" in r for r in best_reasons)
        # The high-pressure entry is flagged non-ambient.
        worst, worst_reasons = ordered[-1]
        assert worst["id"] == "9000001"
        assert "non-ambient" in worst_reasons

    def test_ambient_when_conditions_blank(self):
        # COD leaves T/P blank for room-condition measurements.
        assert cod_fish._is_ambient(
            {"celltemp": "", "diffrtemp": "", "cellpressure": "",
             "diffrpressure": ""})

    def test_high_temperature_marks_non_ambient(self):
        assert not cod_fish._is_ambient(
            {"celltemp": "900", "diffrtemp": "", "cellpressure": "",
             "diffrpressure": ""})


class TestIndexResolution:
    def _session(self):
        return [{"id": "1111111"}, {"id": "2222222"},
                {"id": "3333333"}]

    def test_indices_map_to_ids(self):
        assert cod_fish.resolve_ids(["1", "3"], self._session()) == [
            "1111111", "3333333"]

    def test_ranges_expand(self):
        assert cod_fish.resolve_ids(["1-3"], self._session()) == [
            "1111111", "2222222", "3333333"]

    def test_large_number_is_raw_cod_id(self):
        # A seven-digit value above the row count is a raw id, not an
        # index, so raw ids still work with a session present.
        assert cod_fish.resolve_ids(["9008463"], self._session()) == [
            "9008463"]

    def test_raw_ids_without_session(self):
        assert cod_fish.resolve_ids(["9008463"], None) == ["9008463"]

    def test_non_numeric_token_errors(self):
        with pytest.raises(CodFishError, match="row index or COD id"):
            cod_fish.resolve_ids(["wat"], self._session())


class TestRevisionGuard:
    def _stub(self, monkeypatch, body):
        monkeypatch.setattr(
            cod_fish, "_http_get", lambda url, timeout=60: body)

    def test_cif_revision_parsed(self):
        assert cod_fish.cif_revision(
            b"#$Revision: 291735 $\ndata_x\n") == "291735"

    def test_fetch_verifies_matching_revision(self, monkeypatch):
        self._stub(monkeypatch, b"#$Revision: 42 $\ndata_x\n")
        # Matching revision returns the bytes without error.
        assert cod_fish.fetch_cif("123", revision="42").startswith(
            b"#$Revision")

    def test_fetch_refuses_mismatched_revision(self, monkeypatch):
        self._stub(monkeypatch, b"#$Revision: 42 $\ndata_x\n")
        with pytest.raises(CodFishError, match="does not match"):
            cod_fish.fetch_cif("123", revision="99")

    def test_fetch_without_revision_skips_check(self, monkeypatch):
        self._stub(monkeypatch, b"#$Revision: 42 $\ndata_x\n")
        assert cod_fish.fetch_cif("123") == b"#$Revision: 42 $\ndata_x\n"


class TestAutoReferenceId:
    """reference_id derived from CIF metadata: formula + H-M symbol +
    IT number + publication year (DESIGN 5.7)."""

    _DIAMOND = (b"#$Revision: 100 $\n"
                b"_chemical_formula_sum 'Si'\n"
                b"_symmetry_space_group_name_H-M 'F d -3 m :1'\n"
                b"_space_group_IT_number 227\n"
                b"_journal_year 2010\n")

    def test_builds_formula_symbol_number_year(self):
        assert cod_fish._auto_reference_id(self._DIAMOND) == \
            "si_fd-3m_227_2010"

    def test_drops_setting_and_slash_from_symbol(self):
        cif = (b"_chemical_formula_sum 'Si'\n"
               b"_symmetry_space_group_name_H-M 'P 63/m m c'\n"
               b"_space_group_IT_number 194\n_journal_year 1984\n")
        assert cod_fish._auto_reference_id(cif) == "si_p63mmc_194_1984"

    def test_reduces_formula_unit(self):
        cif = (b"_chemical_formula_sum 'Fe4 O6'\n"
               b"_symmetry_space_group_name_H-M 'R -3 c'\n"
               b"_space_group_IT_number 167\n_journal_year 1990\n")
        assert cod_fish._auto_reference_id(cif) == "fe2o3_r-3c_167_1990"

    def test_omits_year_when_absent(self):
        cif = (b"_chemical_formula_sum 'Si'\n"
               b"_symmetry_space_group_name_H-M 'I m m a'\n"
               b"_space_group_IT_number 74\n")
        assert cod_fish._auto_reference_id(cif) == "si_imma_74"

    def test_none_when_space_group_missing(self):
        cif = b"_chemical_formula_sum 'Si'\n_journal_year 2010\n"
        assert cod_fish._auto_reference_id(cif) is None


class TestManifestFragment:
    def test_fragment_carries_id_and_revision(self):
        fragment = cod_fish._manifest_fragment(
            [{"id": "9008463", "revision": "291735"}])
        assert "cod_id = 9008463" in fragment
        assert 'cod_revision = "291735"' in fragment
        assert "reference_id" in fragment

    def test_emits_schema_version_header(self):
        # A pin's output is a complete, ready-to-read sketch.
        fragment = cod_fish._manifest_fragment(
            [{"id": "1", "revision": "1",
              "reference_id": "si_imma_74_1993"}])
        assert fragment.startswith("schema_version = 2")
        assert 'reference_id = "si_imma_74_1993"' in fragment

    def test_uses_auto_reference_id_when_present(self):
        fragment = cod_fish._manifest_fragment(
            [{"id": "2104737", "revision": "201401",
              "reference_id": "si_fd-3m_227_2010"}])
        assert 'reference_id = "si_fd-3m_227_2010"' in fragment
        assert "cod_2104737" not in fragment

    def test_disambiguates_colliding_reference_ids(self):
        # Two stubs reducing to the same name: the second gets a
        #   trailing counter so reference_id stays unique (rule 5).
        fragment = cod_fish._manifest_fragment([
            {"id": "1", "revision": "1",
             "reference_id": "si_p63mmc_194"},
            {"id": "2", "revision": "2",
             "reference_id": "si_p63mmc_194"}])
        assert 'reference_id = "si_p63mmc_194"' in fragment
        assert 'reference_id = "si_p63mmc_194_2"' in fragment

    def test_falls_back_to_cod_id_without_reference_id(self):
        fragment = cod_fish._manifest_fragment(
            [{"id": "9008463", "revision": "1"}])
        assert 'reference_id = "cod_9008463"' in fragment


class TestStoichiometry:
    """Element specs carry optional counts (2Fe / Fe2); by default a
    count selects the exact reduced formula unit, while --fuzzy hands
    the count to COD's own broader match."""

    # file 2 (Fe4 O6, Z=2) reduces to Fe2O3; file 4 (Fe25 O32) and
    # file 5 (non-stoichiometric) do not -- the exact filter keeps 1,2.
    _CSV = ("# header comment\n"
            "file,formula,sgNumber\n"
            "1,- Fe2 O3 -,167\n"
            "2,- Fe4 O6 -,167\n"
            "3,- Fe3 O4 -,227\n"
            "4,- Fe25 O32 -,1\n"
            "5,- Fe0.911 O -,225\n")

    def _stub(self, monkeypatch):
        monkeypatch.setattr(
            cod_fish, "_http_get",
            lambda url, timeout=120: self._CSV.encode())

    def test_spec_count_either_side(self):
        assert cod_fish._parse_element_spec("2Fe") == (2, "Fe")
        assert cod_fish._parse_element_spec("Fe2") == (2, "Fe")
        assert cod_fish._parse_element_spec("Fe") == (None, "Fe")
        assert cod_fish._parse_element_spec("o") == (None, "O")

    def test_spec_rejects_double_count(self):
        with pytest.raises(CodFishError):
            cod_fish._parse_element_spec("2Fe2")

    def test_spec_rejects_garbage(self):
        with pytest.raises(CodFishError):
            cod_fish._parse_element_spec("42")

    def test_count_box_leading_space_for_one_letter(self):
        assert cod_fish._cod_count_box("Fe", 2) == "Fe 2"
        assert cod_fish._cod_count_box("O", 3) == " O 3"

    def test_reduced_formula_gcd_and_nonstoich(self):
        assert cod_fish._reduced_formula(
            {"Fe": 4.0, "O": 6.0}) == {"Fe": 2, "O": 3}
        assert cod_fish._reduced_formula(
            {"Fe": 0.911, "O": 1.0}) is None

    def test_default_keeps_exact_formula_unit(self, monkeypatch):
        self._stub(monkeypatch)
        ids = [r["id"] for r in cod_fish.search(["2Fe", "3O"])]
        assert ids == ["1", "2"]      # Fe2O3 and Fe4O6 (->Fe2O3)

    def test_fuzzy_skips_client_filter(self, monkeypatch):
        self._stub(monkeypatch)
        ids = [r["id"] for r in cod_fish.search(["2Fe", "3O"],
                                                fuzzy=True)]
        assert ids == ["1", "2", "3", "4", "5"]
