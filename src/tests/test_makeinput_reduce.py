"""test_makeinput_reduce.py -- Regression + unit tests for the reduce
grouping scheme and the Phase-2 matcher protocol (C54; ARCHITECTURE
8.9, DESIGN 5.6.4).

The reduce scheme groups atoms of one element into species by how
similar their concentric-shell neighborhoods are.  C54 refactors the
long-standing ``group_reduce`` algorithm so its two halves live behind
the uniform *matcher* protocol: the per-atom shell-configuration build
becomes ``ReduceMatcher.compute_query`` and the three equivalence
tests (level-distance tolerance, neighbor count, neighbor composition)
become ``ReduceMatcher.distance``.  Nothing a user sees may change --
``group_reduce`` must assign exactly the same species it always did.

These tests pin that contract on two axes:

* **Regression on ``group_reduce``** -- two hand-traced fixtures (one
  where two same-element atoms share an environment and group, one
  where their nearest-neighbor distances differ beyond tolerance and
  they split) lock the species assignment the algorithm has always
  produced.  They are written to pass against the pre-refactor code so
  that re-running them after the refactor proves behavior was
  preserved.
* **Unit coverage of the matcher surface** -- the registry wiring, the
  shell-code fingerprints ``compute_query`` emits, the 0-or-infinite
  ``distance`` across its element / distance / count / composition
  branches, and the ``representative`` reduction.

conftest.py's ``SCRIPTS_DIR`` insertion lets us import ``makeinput``
directly without installing the package.  The fixtures fabricate the
small slice of ``ScriptSettings`` / ``StructureControl`` state that
``group_reduce`` actually touches (a minimum-distance matrix plus the
per-atom element / species / name arrays) rather than reading a real
structure, so the tests stay pure and fast.
"""

import math
import types

import pytest

import makeinput
from matchers import MATCHERS, ReduceMatcher


# Pure tests: no Fortran binaries, no database reads, no full build.
pytestmark = pytest.mark.unit


# ==============================================================
#  Fixture builders
# ==============================================================

def _make_reduce_state(min_dist, element_id, element_name,
                       species_id, reduce_spec):
    """Fabricate the ``(settings, sc)`` pair ``group_reduce`` needs.

    Only the handful of attributes the reduce algorithm reads or
    writes are populated.  All per-atom arrays are 1-indexed (index 0
    is an unused ``None`` placeholder), matching the convention used
    throughout makeinput.

    Parameters
    ----------
    min_dist : list[list[float]]
        The symmetric, 1-indexed minimum-image distance matrix placed
        on the fake ``StructureControl`` as ``sc.min_dist``.
    element_id, element_name, species_id : list
        1-indexed per-atom element index, element name, and the
        *input* species index (the assignment reduce starts from).
    reduce_spec : dict
        One reduce specification: ``level``, ``thick``, ``cutoff``,
        ``tolerance``, ``op``.
    """

    num_atoms = len(element_id) - 1
    num_elements = max(eid for eid in element_id[1:])

    settings = types.SimpleNamespace(
        reduces=[reduce_spec],
        num_atoms=num_atoms,
        num_elements=num_elements,
        atom_element_id=list(element_id),
        atom_species_id=list(species_id),
        atom_element_name=list(element_name),
        # Written by group_reduce; sized for 1-indexed element access.
        num_species=[0] * (num_elements + 1),
        num_types=None,
    )
    sc = types.SimpleNamespace(min_dist=min_dist)
    return settings, sc


# A reduce specification shared by the regression fixtures: one shell,
# a thin 0.1 Angstrom acceptance band, a generous 5.0 Angstrom cutoff,
# and a 5% level-distance tolerance.
_ONE_SHELL = {"level": 1, "thick": 0.1, "cutoff": 5.0,
              "tolerance": 0.05, "op": "species"}


# ==============================================================
#  Regression: group_reduce species assignment (the user contract)
# ==============================================================

def test_group_reduce_groups_equivalent_environments(tmp_path,
                                                     monkeypatch):
    """Two same-element atoms with identical nearest-neighbor shells
    collapse to one species.

    Layout (distances in Angstrom): atoms 1,2 are element 1 ("Si");
    atoms 3,4 are element 2 ("O").  Atom 1's nearest neighbor is atom
    3 at 2.0 and atom 2's is atom 4 at 2.0, so the two Si atoms see an
    identical single-O shell at the same distance and must group; by
    symmetry the two O atoms group as well.  Expected final
    assignment: every atom is species 1 of its element.
    """

    monkeypatch.chdir(tmp_path)

    #            self  1    2    3    4
    min_dist = [
        [None, None, None, None, None],   # index 0 placeholder
        [None, 0.0,  3.0,  2.0,  2.6],    # atom 1
        [None, 3.0,  0.0,  2.6,  2.0],    # atom 2
        [None, 2.0,  2.6,  0.0,  3.0],    # atom 3
        [None, 2.6,  2.0,  3.0,  0.0],    # atom 4
    ]
    settings, sc = _make_reduce_state(
        min_dist,
        element_id=[None, 1, 1, 2, 2],
        element_name=[None, "Si", "Si", "O", "O"],
        species_id=[None, 1, 1, 1, 1],
        reduce_spec=dict(_ONE_SHELL),
    )

    makeinput.group_reduce(settings, sc, 0)

    assert settings.atom_species_id == [None, 1, 1, 1, 1]
    assert settings.num_species[1] == 1
    assert settings.num_species[2] == 1
    # reduce only assigns species, so every species owns exactly one
    # type after the pass (file set 0).
    assert settings.num_types[0][1] == [None, 1]
    assert settings.num_types[0][2] == [None, 1]


def test_group_reduce_splits_distinct_environments(tmp_path,
                                                   monkeypatch):
    """Two same-element atoms whose nearest-neighbor distances differ
    by more than the tolerance split into separate species.

    Same four atoms, but atom 2's nearest neighbor sits at 2.6
    Angstrom while atom 1's sits at 2.0.  The level-distance test
    allows only ``0.05 * 2.0 = 0.1`` of slack, so the 0.6 Angstrom gap
    fails and the two Si atoms cannot be the same species; the two O
    atoms split for the same reason.  Expected: species 1 and 2 within
    each element.
    """

    monkeypatch.chdir(tmp_path)

    #            self  1    2    3    4
    min_dist = [
        [None, None, None, None, None],   # index 0 placeholder
        [None, 0.0,  3.5,  2.0,  3.0],    # atom 1
        [None, 3.5,  0.0,  3.0,  2.6],    # atom 2
        [None, 2.0,  3.0,  0.0,  3.5],    # atom 3
        [None, 3.0,  2.6,  3.5,  0.0],    # atom 4
    ]
    settings, sc = _make_reduce_state(
        min_dist,
        element_id=[None, 1, 1, 2, 2],
        element_name=[None, "Si", "Si", "O", "O"],
        species_id=[None, 1, 1, 1, 1],
        reduce_spec=dict(_ONE_SHELL),
    )

    makeinput.group_reduce(settings, sc, 0)

    assert settings.atom_species_id == [None, 1, 2, 1, 2]
    assert settings.num_species[1] == 2
    assert settings.num_species[2] == 2


def test_group_reduce_rejects_non_species_op(tmp_path, monkeypatch):
    """The reduce scheme supports only species grouping; a ``types``
    op is a contract fault and raises ``MakeinputError`` rather than
    killing a worker with ``sys.exit`` (DESIGN 6.3.1)."""

    monkeypatch.chdir(tmp_path)

    spec = dict(_ONE_SHELL)
    spec["op"] = "types"
    settings, sc = _make_reduce_state(
        [[None, None, None], [None, 0.0, 2.0], [None, 2.0, 0.0]],
        element_id=[None, 1, 1],
        element_name=[None, "Si", "Si"],
        species_id=[None, 1, 1],
        reduce_spec=spec,
    )

    with pytest.raises(makeinput.MakeinputError):
        makeinput.group_reduce(settings, sc, 0)


# ==============================================================
#  Matcher protocol surface (ARCHITECTURE 8.9)
# ==============================================================

def test_reduce_matcher_registered():
    """The registry maps the manifest ``method`` name to the class,
    and the class advertises the Python-side, no-loen contract."""

    assert MATCHERS["reduce"] is ReduceMatcher
    matcher = ReduceMatcher()
    assert matcher.name == "reduce"
    assert matcher.needs_loen_run is False
    assert matcher.default_similarity_floor == 0.05


def _structure_view(min_dist, element_id, element_name, species_id):
    """A minimal duck-typed structure for the matcher: the attribute
    names match both makeinput's internal view and StructureControl."""

    num_atoms = len(element_id) - 1
    return types.SimpleNamespace(
        num_atoms=num_atoms,
        num_elements=max(eid for eid in element_id[1:]),
        atom_element_id=list(element_id),
        atom_species_id=list(species_id),
        atom_element_name=list(element_name),
        min_dist=min_dist,
    )


def test_compute_query_emits_shell_codes():
    """``compute_query`` returns one shell-code fingerprint per atom,
    each carrying its central element, the per-level distance, and the
    (element, species) multiset of that shell's neighbors."""

    min_dist = [
        [None, None, None, None, None],
        [None, 0.0,  3.0,  2.0,  2.6],
        [None, 3.0,  0.0,  2.6,  2.0],
        [None, 2.0,  2.6,  0.0,  3.0],
        [None, 2.6,  2.0,  3.0,  0.0],
    ]
    structure = _structure_view(
        min_dist,
        element_id=[None, 1, 1, 2, 2],
        element_name=[None, "Si", "Si", "O", "O"],
        species_id=[None, 1, 1, 1, 1],
    )
    fingerprints = ReduceMatcher().compute_query(structure, _ONE_SHELL)

    atom_one = fingerprints[1]
    assert atom_one.element_id == 1
    assert atom_one.element_name == "Si"
    assert atom_one.tolerance == 0.05
    # One level, holding the single O neighbor (element 2, species 1)
    # at 2.0 Angstrom.
    assert atom_one.levels[1].distance == pytest.approx(2.0)
    assert atom_one.levels[1].members == [(2, 1)]
    assert atom_one.levels[1].member_names == ["O"]


def test_distance_zero_for_equivalent_and_inf_otherwise():
    """``distance`` is 0 for shell codes that pass all three tests and
    infinite the moment any test fails -- element, level distance,
    neighbor count, or neighbor composition."""

    min_dist = [
        [None, None, None, None, None],
        [None, 0.0,  3.0,  2.0,  2.6],
        [None, 3.0,  0.0,  2.6,  2.0],
        [None, 2.0,  2.6,  0.0,  3.0],
        [None, 2.6,  2.0,  3.0,  0.0],
    ]
    structure = _structure_view(
        min_dist,
        element_id=[None, 1, 1, 2, 2],
        element_name=[None, "Si", "Si", "O", "O"],
        species_id=[None, 1, 1, 1, 1],
    )
    matcher = ReduceMatcher()
    fps = matcher.compute_query(structure, _ONE_SHELL)

    # Atoms 1 and 2 see identical single-O shells -> equivalent.
    assert matcher.distance(fps[1], fps[2]) == 0.0
    # Atom 1 (Si) versus atom 3 (O): different central element.
    assert matcher.distance(fps[1], fps[3]) == math.inf


def test_representative_returns_first_member():
    """For reduce, intra-species fingerprints agree by construction,
    so the representative is simply the first member (DESIGN 5.6.5)."""

    sentinel_first = object()
    members = [sentinel_first, object(), object()]
    assert ReduceMatcher().representative(members) is sentinel_first


# ==============================================================
#  Payload serialization (DESIGN 5.2 element-only shell_code)
# ==============================================================

def test_build_payload_is_element_only_shell_code():
    """``build_payload`` serializes a shell code into the stored,
    cross-structure form: the central element symbol plus per-level
    distance and neighbor element symbols, all lowercased, with NO
    species component (which would not transfer across structures)."""

    min_dist = [
        [None, None, None, None, None],
        [None, 0.0,  3.0,  2.0,  2.6],
        [None, 3.0,  0.0,  2.6,  2.0],
        [None, 2.0,  2.6,  0.0,  3.0],
        [None, 2.6,  2.0,  3.0,  0.0],
    ]
    structure = _structure_view(
        min_dist,
        element_id=[None, 1, 1, 2, 2],
        element_name=[None, "Si", "Si", "O", "O"],
        species_id=[None, 1, 1, 1, 1],
    )
    matcher = ReduceMatcher()
    fps = matcher.compute_query(structure, _ONE_SHELL)

    payload = matcher.build_payload(fps[1])
    assert payload == {
        "shell_code": {
            "element": "si",
            "levels": [
                {"distance": pytest.approx(2.0), "neighbors": ["o"]},
            ],
        }
    }
    # No species leaked into the transferable descriptor.
    neighbors = payload["shell_code"]["levels"][0]["neighbors"]
    assert all(isinstance(name, str) for name in neighbors)


def test_extract_query_vector_unwraps_build_payload():
    """``extract_query_vector`` is the inverse of ``build_payload`` --
    it hands back the stored ``shell_code`` dict unchanged."""

    matcher = ReduceMatcher()
    shell_code = {"element": "o",
                  "levels": [{"distance": 1.98,
                              "neighbors": ["si", "si"]}]}
    payload = {"shell_code": shell_code}
    assert matcher.extract_query_vector(payload) is shell_code
