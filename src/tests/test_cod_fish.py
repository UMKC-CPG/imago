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


class TestManifestFragment:
    def test_fragment_carries_id_and_revision(self):
        fragment = cod_fish._manifest_fragment(
            [{"id": "9008463", "revision": "291735"}])
        assert "cod_id = 9008463" in fragment
        assert 'cod_revision = "291735"' in fragment
        assert "reference_id" in fragment
