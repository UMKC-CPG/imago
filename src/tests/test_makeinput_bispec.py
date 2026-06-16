"""test_makeinput_bispec.py -- Unit tests for the bispectrum matcher
(C55, C89; ARCHITECTURE 8.9, DESIGN 5.10, PSEUDOCODE 11.3).

The bispectrum descriptor captures the *angular* arrangement of an
atom's neighbors and is computed by the Imago Fortran engine (the
"loen" path), not in Python -- so ``BispecMatcher`` carries
``needs_loen_run = True``.  It supplies every part of the matcher
that does not itself run the engine:

* ``to_loen_input`` -- maps a ``sub_spec`` to the full ``LOEN_INPUT_DATA``
  parameter dict, with the three optional parameters defaulting to the
  values makeinput has historically hardcoded (DESIGN 5.10.5).  The
  reserved ``by_element`` key is refused until the element-aware variant
  lands (TODO C62).
* ``parse_loen_output`` -- reads a self-describing ``fort.21`` (DESIGN
  5.10.3) into one ``LoenSite`` per site: the identity columns (site#,
  element, species, type_in_species, type_flat) plus the leading
  ``twoj2 + 1`` components, skipping the header and dropping the sum.
* ``distance`` -- the Euclidean distance between two vectors.
* ``representative`` -- the element-wise mean of a species' members.
* ``build_payload`` / ``extract_query_vector`` -- the on-disk ``values``
  (de)serialization.

``compute_query`` is not implemented for bispectrum (the makegroups flow
reads ``fort.21`` instead) and is not exercised here.  conftest.py's
``SCRIPTS_DIR`` insertion lets us import the matchers module directly
without installing the package; these tests touch no Fortran binary --
the ``fort.21`` reader is fed a hand-written fixture file via ``tmp_path``.
"""

import math

import pytest

from matchers import MATCHERS, BispecMatcher


# Pure tests: no Fortran binaries, no database reads, no full build.
pytestmark = pytest.mark.unit


# ==============================================================
#  Registry + protocol contract (ARCHITECTURE 8.9)
# ==============================================================

def test_bispec_matcher_registered():
    """The registry maps the manifest ``method`` name to the class, and
    the class advertises the loen-side contract."""

    assert MATCHERS["bispectrum"] is BispecMatcher
    matcher = BispecMatcher()
    assert matcher.name == "bispectrum"
    assert matcher.needs_loen_run is True
    assert matcher.default_similarity_floor == 0.10


# ==============================================================
#  to_loen_input -- the LOEN parameter contract (DESIGN 5.10.5)
# ==============================================================

def test_to_loen_input_fills_defaults():
    """With only the required ``twoj1`` / ``twoj2``, the optional
    parameters take the documented defaults that reproduce makeinput's
    historically-hardcoded LOEN block."""

    params = BispecMatcher().to_loen_input({"twoj1": 8, "twoj2": 6})
    assert params == {
        "loenCode":     1,
        "twoj1":        8,
        "twoj2":        6,
        "max_neigh":    20,
        "cutoff":       5.0,
        "angleSqueeze": 0.85,
    }


def test_to_loen_input_passes_through_overrides():
    """Supplied optional parameters override the defaults; the snake_case
    ``angle_squeeze`` sub_spec key maps to the camelCase ``angleSqueeze``
    LOEN parameter the Fortran block reads."""

    params = BispecMatcher().to_loen_input({
        "twoj1": 8, "twoj2": 6,
        "max_neigh": 32, "cutoff": 7.5, "angle_squeeze": 0.5,
    })
    assert params["max_neigh"] == 32
    assert params["cutoff"] == 7.5
    assert params["angleSqueeze"] == 0.5


@pytest.mark.parametrize("missing", ["twoj1", "twoj2"])
def test_to_loen_input_requires_both_twoj(missing):
    """``twoj1`` and ``twoj2`` are required; dropping either is a hard
    error rather than a silent default."""

    sub_spec = {"twoj1": 8, "twoj2": 6}
    del sub_spec[missing]
    with pytest.raises(ValueError, match=missing):
        BispecMatcher().to_loen_input(sub_spec)


def test_to_loen_input_rejects_by_element():
    """The element-aware variant is not built yet (TODO C62 / D10), so
    a ``by_element`` sub_spec is refused loudly rather than silently
    producing element-blind parameters."""

    with pytest.raises(NotImplementedError, match="by_element"):
        BispecMatcher().to_loen_input(
            {"twoj1": 8, "twoj2": 6, "by_element": True})


# ==============================================================
#  parse_loen_output -- the fort.21 reader
# ==============================================================

def test_parse_loen_output_reads_per_site_records(tmp_path):
    """Each ``fort.21`` data row yields one ``LoenSite``: the identity
    columns (site#, element, species, type_in_species, type_flat) plus the
    leading ``twoj2 + 1`` components.  The header is skipped, the trailing
    sum dropped, and blank lines ignored so a trailing newline is not a
    phantom site (DESIGN 5.10.3)."""

    # twoj2 = 2 -> twoj2 + 1 = 3 components, after 5 identity columns and
    #   before the ignored sum column.
    fort21 = tmp_path / "fort.21"
    fort21.write_text(
        "  site#  element  species  type_sp  type_flat"
        "      2j_NNN   total\n"
        "      1       Si        1        1          1"
        "   1.0 2.0 3.0   6.0\n"
        "      2        O        1        1          2"
        "   0.5 0.5 0.5   1.5\n"
        "\n")               # trailing blank line: must be skipped
    sites = BispecMatcher().parse_loen_output(
        str(fort21), {"twoj1": 4, "twoj2": 2})

    assert len(sites) == 2
    first = sites[0]
    assert first.site == 1
    assert first.element == "Si"
    assert first.species == 1
    assert first.type_in_species == 1
    assert first.type_flat == 1
    assert first.vector == [1.0, 2.0, 3.0]
    assert sites[1].element == "O"
    assert sites[1].type_flat == 2
    assert sites[1].vector == [0.5, 0.5, 0.5]


def test_parse_loen_output_rejects_short_row(tmp_path):
    """A data row with fewer than 5 identity + ``twoj2 + 1`` component
    columns means the file does not match the declared ``sub_spec``; that
    is a hard error, not a truncated record."""

    fort21 = tmp_path / "fort.21"
    fort21.write_text(
        "  site#  element  species  type_sp  type_flat  2j  total\n"
        "      1       Si        1        1\n")   # missing components
    with pytest.raises(ValueError, match="columns"):
        BispecMatcher().parse_loen_output(
            str(fort21), {"twoj1": 4, "twoj2": 2})


# ==============================================================
#  distance + representative (DESIGN 5.6.5)
# ==============================================================

def test_distance_is_euclidean():
    """``distance`` is the plain L2 norm of the difference: the classic
    3-4-5 right triangle gives 5."""

    matcher = BispecMatcher()
    assert matcher.distance([0.0, 0.0, 0.0],
                            [3.0, 4.0, 0.0]) == pytest.approx(5.0)
    assert matcher.distance([1.0, 1.0], [1.0, 1.0]) == 0.0


def test_distance_rejects_unequal_length():
    """Vectors of different length cannot have come from the same
    ``sub_spec``; comparing them is a hard error."""

    with pytest.raises(ValueError, match="unequal length"):
        BispecMatcher().distance([1.0, 2.0, 3.0], [1.0, 2.0])


def test_representative_is_elementwise_mean():
    """The representative is the slot-by-slot arithmetic mean of the
    members -- the order-independent center the bispectrum species pass
    queries the database with (DESIGN 5.6.5)."""

    members = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
    assert BispecMatcher().representative(members) == [3.0, 4.0]


def test_representative_rejects_empty():
    """An empty member list has no mean; refuse rather than divide by
    zero."""

    with pytest.raises(ValueError, match="at least one member"):
        BispecMatcher().representative([])


# ==============================================================
#  Payload serialization (DESIGN 5.2 / 5.4 -- the `values` field)
# ==============================================================

def test_build_payload_and_extract_roundtrip():
    """``build_payload`` stores the vector under ``values`` and
    ``extract_query_vector`` is its inverse, handing the vector back."""

    matcher = BispecMatcher()
    vector = [0.1, -0.2, 0.3, 0.4, 0.5]
    payload = matcher.build_payload(vector)
    assert payload == {"values": [0.1, -0.2, 0.3, 0.4, 0.5]}
    assert matcher.extract_query_vector(payload) == vector


def test_distance_consistent_with_representative_center():
    """Cross-check that the mean really is the metric center the distance
    rewards: the representative of a symmetric pair is equidistant from
    both members."""

    matcher = BispecMatcher()
    a = [0.0, 0.0]
    b = [2.0, 4.0]
    center = matcher.representative([a, b])
    assert center == [1.0, 2.0]
    assert matcher.distance(center, a) == pytest.approx(
        matcher.distance(center, b))
    assert matcher.distance(center, a) == pytest.approx(math.sqrt(5.0))
