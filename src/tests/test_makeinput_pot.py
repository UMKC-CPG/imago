"""test_makeinput_pot.py -- Unit tests for makeinput.py's
augmented-potential-database consumer helpers (C47).

makeinput.py is the *consumer* in the library / producer /
consumer split of DESIGN 5.4.  C47 wires the augmented per-element
potential database (``s_gaussian_pot.toml``) into makeinput's
potential-emission path, following the reduced (no-environment-
matcher) selection flow of PSEUDOCODE 11.3.0.

makeinput has no end-to-end test harness (a full run needs a real
structure, the contract/basis databases, and Imago binaries), so
these tests pin the two *pure* helpers that carry the new logic:

* ``_select_augmented_pot_entry`` -- the reduced entry pick:
  ``-pot LABEL`` override (fatal on a missing label), else the
  database's default-tagged entry.
* ``_write_legacy_pot_files_from_entry`` -- materializes a chosen
  ``PotentialEntry`` into the historical ``pot``/``coeff`` text
  files that the imago.dat writer and Imago's ``scfV.dat`` reader
  consume, so the augmented path rejoins the legacy code path
  byte-for-byte.

conftest.py's ``SCRIPTS_DIR`` insertion lets us import both
``makeinput`` and ``initial_potential_db`` directly.
"""

import os

import pytest

import makeinput
from makeinput import ScriptSettings
import initial_potential_db as ipdb
from initial_potential_db import PotentialEntry, ElementDatabase


# Pure-computation / temp-file unit tests -- no Fortran binaries,
# no full makeinput run.
pytestmark = pytest.mark.unit


# ============================================================
#  Fixtures
# ============================================================

def _atomscf_provenance() -> dict:
    """Return a valid atomSCF-source provenance dict."""

    return {
        "source": "atomSCF",
        "commit": "abcdef1",
        "generated_at": "2026-05-20T14:00:00Z",
    }


def _imago_provenance() -> dict:
    """Return a valid Imago-source provenance dict."""

    return {
        "source": "Imago",
        "commit": "fedcba2",
        "generated_at": "2026-05-20T14:30:00Z",
        "reference_id": "COD-1011098",
        "atom_site": 1,
        "kpoint_spec": "12 12 12 0 0 0",
        "scf_threshold": 1.0e-6,
        "scf_iterations": 28,
    }


def _au_db() -> ElementDatabase:
    """Return a small but valid schema-v2 Au database with an
    ``isolated`` baseline (not default) and a ``default_solid``
    entry (default).  Three Gaussian terms each.
    """

    db = ElementDatabase(
        schema_version  = 2,
        element_symbol  = "Au",
        nuclear_z       = 79.0,
        nuclear_alpha   = 20.0,
        covalent_radius = 1.0,
    )
    db.potentials.append(PotentialEntry(
        label         = "isolated",
        default       = False,
        description   = "Isolated Au atom from atomSCF.",
        num_gaussians = 3,
        alpha_min     = 1.5e-1,
        alpha_max     = 1.0e+8,
        coefficients    = [1.0, -2.5, 3.0],
        coefficient_std = [0.0, 0.0, 0.0],
        alphas          = [1.5e-1, 1.5e+0, 1.0e+8],
        provenance      = _atomscf_provenance(),
    ))
    db.potentials.append(PotentialEntry(
        label         = "default_solid",
        default       = True,
        description   = "Au in fcc bulk (Fm-3m).",
        num_gaussians = 3,
        alpha_min     = 1.5e-1,
        alpha_max     = 1.0e+8,
        coefficients    = [0.1, 0.2, 0.3],
        coefficient_std = [0.0, 0.0, 0.0],
        alphas          = [1.5e-1, 1.5e+0, 1.0e+8],
        provenance      = _imago_provenance(),
    ))
    return db


# ============================================================
#  _select_augmented_pot_entry (PSEUDOCODE 11.3.0 reduced pick)
# ============================================================

class TestSelectAugmentedPotEntry:
    """The reduced entry pick: -pot override (fatal on miss),
    else the default-tagged entry.
    """

    def test_no_override_returns_default_entry(self):
        db = _au_db()
        entry = makeinput._select_augmented_pot_entry(
            ipdb, db, None, "au")
        assert entry.label == "default_solid"

    def test_override_returns_named_entry(self):
        db = _au_db()
        entry = makeinput._select_augmented_pot_entry(
            ipdb, db, "isolated", "au")
        assert entry.label == "isolated"

    def test_missing_override_label_is_fatal(self):
        # A -pot label absent from the database aborts the build
        # rather than falling back -- a deliberate override must
        # not silently pick a different potential.  It now raises
        # MakeinputError instead of calling sys.exit, so the fault
        # propagates out of a kaleidoscope worker as a failed unit
        # instead of killing the worker (DESIGN 6.3.1, 6.3.5).
        db = _au_db()
        with pytest.raises(makeinput.MakeinputError) as excinfo:
            makeinput._select_augmented_pot_entry(
                ipdb, db, "no_such_label", "au")
        # The message names both the bad label and the element.
        msg = str(excinfo.value)
        assert "no_such_label" in msg
        assert "au" in msg


# ============================================================
#  _select_augmented_pot_entry precedence 2: the per-species
#  environment fingerprint match (DESIGN 5.6.5; C59)
# ============================================================

import types                                            # noqa: E402

from matchers import (ReduceMatcher, ReduceShellCode,    # noqa: E402
                      ReduceShellLevel, LoenSite)
import makegroups                                         # noqa: E402

# The descriptor sub_spec a reduce scheme stores and an entry records.
_REDUCE_SUB_SPEC = {"level": 1, "thick": 0.1,
                    "cutoff": 5.0, "tolerance": 0.05}


def _reduce_entry(label, default, neighbor, distance,
                  preferred=False):
    """A PotentialEntry carrying one reduce fingerprint: a single shell
    with the given neighbor element symbol at the given distance.

    ``preferred`` flags the record as the database's preferred reduce
    representative, which the file-dictated regime (DESIGN 5.6.5)
    consults to choose the family and sub_spec to match on."""

    record = ipdb.FingerprintRecord(
        method="reduce", sub_spec=dict(_REDUCE_SUB_SPEC),
        preferred=preferred,
        payload={"shell_code": {
            "element": "si",
            "levels": [{"distance": distance,
                        "neighbors": [neighbor]}]}})
    return PotentialEntry(
        label=label, default=default, description="d",
        num_gaussians=1, alpha_min=1.0, alpha_max=1.0,
        coefficients=[1.0], coefficient_std=[0.0], alphas=[1.0],
        provenance=_imago_provenance(), fingerprints=[record])


# The bispectrum sub_spec a loen run produces and an entry records;
# twoj2 = 2 fixes the per-vector component count at twoj2 + 1 = 3.
_BISPEC_SUB_SPEC = {"twoj1": 6, "twoj2": 2}


def _bispec_entry(label, default, values, preferred=False):
    """A PotentialEntry carrying one bispectrum fingerprint whose stored
    ``values`` vector is the given list (length twoj2 + 1 = 3)."""

    record = ipdb.FingerprintRecord(
        method="bispectrum", sub_spec=dict(_BISPEC_SUB_SPEC),
        preferred=preferred, payload={"values": list(values)})
    return PotentialEntry(
        label=label, default=default, description="d",
        num_gaussians=1, alpha_min=1.0, alpha_max=1.0,
        coefficients=[1.0], coefficient_std=[0.0], alphas=[1.0],
        provenance=_imago_provenance(), fingerprints=[record])


def _si_fingerprinted_db() -> ElementDatabase:
    """A Si database with two fingerprinted entries -- bulk-like (a Si
    neighbor at 2.35) and oxide-like (an O neighbor at 1.60) -- plus the
    oxide entry flagged default."""

    db = ElementDatabase(
        schema_version=2, element_symbol="Si", nuclear_z=14.0,
        nuclear_alpha=10.0, covalent_radius=1.1)
    db.potentials.append(_reduce_entry("si-bulk", False, "si", 2.35))
    db.potentials.append(_reduce_entry("si-oxide", True, "o", 1.60))
    return db


def _si_query(neighbor, distance):
    """A ReduceShellCode query: central Si with one neighbor of the given
    element symbol at the given distance."""

    return ReduceShellCode(
        element_id=1, element_name="Si", tolerance=0.05,
        levels=[None, ReduceShellLevel(
            distance=distance, members=[(1, 1)],
            member_names=[neighbor])])


class TestSelectAugmentedPotEntryFingerprint:
    """Precedence 2: when a reduce scheme supplies a species' query
    fingerprint, the pick matches it against the entries' stored
    fingerprints instead of defaulting (DESIGN 5.6.5)."""

    def test_match_picks_closest_entry_over_default(self):
        # A query that looks like the bulk environment selects si-bulk,
        #   even though si-oxide is the default-tagged entry.
        db = _si_fingerprinted_db()
        matcher = ReduceMatcher()
        entry = makeinput._select_augmented_pot_entry(
            ipdb, db, None, "Si",
            query=_si_query("si", 2.35), matcher=matcher,
            sub_spec=_REDUCE_SUB_SPEC, species=1)
        assert entry.label == "si-bulk"

    def test_no_match_within_floor_falls_to_default(self, capsys):
        # A query unlike either stored environment exceeds the floor on
        #   both, so the pick warns and falls back to the default entry.
        db = _si_fingerprinted_db()
        matcher = ReduceMatcher()
        entry = makeinput._select_augmented_pot_entry(
            ipdb, db, None, "Si",
            query=_si_query("au", 9.9), matcher=matcher,
            sub_spec=_REDUCE_SUB_SPEC, species=2)
        assert entry.label == "si-oxide"        # the default
        assert "similarity floor" in capsys.readouterr().out

    def test_pot_override_beats_fingerprint_match(self):
        # Precedence 1 outranks precedence 2: an explicit -pot label wins
        #   even when a fingerprint would have matched a different entry.
        db = _si_fingerprinted_db()
        matcher = ReduceMatcher()
        entry = makeinput._select_augmented_pot_entry(
            ipdb, db, "si-oxide", "Si",
            query=_si_query("si", 2.35), matcher=matcher,
            sub_spec=_REDUCE_SUB_SPEC, species=1)
        assert entry.label == "si-oxide"


class TestSpeciesQueryFingerprint:
    """``_species_query_fingerprint`` summarizes a species' atoms into one
    representative, or None when the species has no reduce fingerprint."""

    def _settings(self, atom_fingerprints):
        return types.SimpleNamespace(
            num_atoms=len(atom_fingerprints) - 1,
            atom_element_id=[None, 1, 1, 2],
            atom_species_id=[None, 1, 2, 1],
            atom_reduce_fingerprint=atom_fingerprints)

    def test_returns_first_member_of_species(self):
        seed = _si_query("si", 2.35)
        other = _si_query("si", 2.35)
        settings = self._settings([None, seed, other, _si_query("o", 1.6)])
        result = makeinput._species_query_fingerprint(
            settings, 1, 1, ReduceMatcher())
        assert result is seed                     # first member, species 1

    def test_returns_none_without_fingerprints(self):
        # No reduce scheme ran (the attribute is absent) -> None, which
        #   sends the pick to the default entry.
        settings = types.SimpleNamespace(
            num_atoms=2, atom_element_id=[None, 1, 1],
            atom_species_id=[None, 1, 1])
        assert makeinput._species_query_fingerprint(
            settings, 1, 1, ReduceMatcher()) is None


# ============================================================
#  _write_legacy_pot_files_from_entry (legacy-format emission)
# ============================================================

class TestWriteLegacyPotFiles:
    """The generated pot/coeff files match the fixed-line legacy
    format that _print_scf_pot reads positionally and that
    Imago's scfV.dat reader parses.
    """

    def _emit(self, tmp_path, entry_label="default_solid"):
        db = _au_db()
        entry = ipdb.lookup(db, entry_label)
        pot_path = str(tmp_path / "pot_aug_au")
        coeff_path = str(tmp_path / "coeff_aug_au")
        makeinput._write_legacy_pot_files_from_entry(
            db, entry, pot_path, coeff_path)
        return db, entry, pot_path, coeff_path

    def test_pot_file_eight_line_layout(self, tmp_path):
        db, entry, pot_path, _ = self._emit(tmp_path)
        with open(pot_path) as handle:
            lines = handle.read().splitlines()
        # Positions match how _print_scf_pot indexes pot_lines.
        assert lines[0] == "NUCLEAR_CHARGE__ALPHA"
        # Line 1: Z and nuclear alpha.  _print_scf_pot reads
        # float(split()[0]) as the nuclear charge.
        assert float(lines[1].split()[0]) == pytest.approx(79.0)
        assert float(lines[1].split()[1]) == pytest.approx(20.0)
        assert lines[2] == "COVALENT_RADIUS"
        assert float(lines[3]) == pytest.approx(1.0)
        assert lines[4] == "NUM_ALPHAS"
        # Line 5: the Gaussian count, read as int(split()[0]).
        assert int(lines[5].split()[0]) == entry.num_gaussians
        assert lines[6] == "ALPHAS"
        # Line 7: the alpha min/max range.
        lo, hi = lines[7].split()
        assert float(lo) == pytest.approx(entry.alpha_min)
        assert float(hi) == pytest.approx(entry.alpha_max)

    def test_coeff_file_count_and_term_lines(self, tmp_path):
        db, entry, _, coeff_path = self._emit(tmp_path)
        with open(coeff_path) as handle:
            lines = handle.read().splitlines()
        # First line is the term count (read as int by both
        # _print_scf_pot's modpot path and Imago).
        assert int(lines[0].split()[0]) == entry.num_gaussians
        term_lines = lines[1:]
        assert len(term_lines) == entry.num_gaussians
        # Imago's reader does read(8,*) over five values per term
        # and consumes only column 1 (the coefficient); the line
        # must therefore carry exactly five whitespace-separated
        # fields, with column 1 equal to the entry coefficient.
        for term_line, coeff in zip(term_lines, entry.coefficients):
            fields = term_line.split()
            assert len(fields) == 5
            assert float(fields[0]) == pytest.approx(coeff)
        # Trailing three columns are the ignored placeholders.
        for term_line in term_lines:
            fields = term_line.split()
            assert float(fields[2]) == 0.0
            assert float(fields[3]) == 0.0
            assert float(fields[4]) == 0.0

    def test_isolated_entry_emits_its_own_coefficients(
            self, tmp_path):
        # Picking a different label emits that entry's numbers.
        _, entry, _, coeff_path = self._emit(
            tmp_path, entry_label="isolated")
        with open(coeff_path) as handle:
            lines = handle.read().splitlines()
        first_term = lines[1].split()
        assert float(first_term[0]) == pytest.approx(1.0)


# ============================================================
#  _summarize_species and _subspec_cache_key
# ============================================================

class TestSummarizeSpecies:
    """``_summarize_species`` filters a per-atom descriptor array to
    one species' atoms and asks the matcher for a representative, or
    returns None when the species has no usable descriptor."""

    def _settings(self):
        # Three atoms: atoms 1,2 are element 1 species 1; atom 3 is
        # element 2 species 1.
        return types.SimpleNamespace(
            num_atoms=3,
            atom_element_id=[None, 1, 1, 2],
            atom_species_id=[None, 1, 1, 1])

    def test_returns_first_member_for_reduce(self):
        seed = _si_query("si", 2.35)
        per_atom = [None, seed, _si_query("si", 2.35),
                    _si_query("o", 1.6)]
        result = makeinput._summarize_species(
            ReduceMatcher(), per_atom, self._settings(), 1, 1)
        assert result is seed                  # first member, species 1

    def test_returns_none_when_no_member(self):
        per_atom = [None, None, None, _si_query("o", 1.6)]
        # Element 1 / species 1 atoms carry no descriptor -> None.
        result = makeinput._summarize_species(
            ReduceMatcher(), per_atom, self._settings(), 1, 1)
        assert result is None


class TestSubspecCacheKey:
    """``_subspec_cache_key`` is hashable and ignores key order so a
    whole-structure query is computed once per distinct sub_spec."""

    def test_order_independent_and_hashable(self):
        key_a = makeinput._subspec_cache_key(
            {"level": 1, "thick": 0.1, "cutoff": 5.0})
        key_b = makeinput._subspec_cache_key(
            {"cutoff": 5.0, "level": 1, "thick": 0.1})
        assert key_a == key_b
        assert len({key_a, key_b}) == 1        # usable as a dict key


# ============================================================
#  _resolve_species_query -- the two-regime query resolution
#  (DESIGN 5.6.5; the C93 fingerprint-pick decoupling)
# ============================================================

def _si_dimer_settings():
    """A two-atom Si dimer's settings: both atoms element 1 species 1,
    2.35 Angstrom apart, with the minimum-distance matrix already built
    so the reduce query can be computed without a real StructureControl.
    """

    settings = types.SimpleNamespace(
        num_atoms=2,
        num_elements=1,
        atom_element_id=[None, 1, 1],
        atom_species_id=[None, 1, 1],
        atom_element_name=[None, "Si", "Si"],
        _min_dist_made=True)
    sc = types.SimpleNamespace(min_dist=[
        None,
        [None, 0.0, 2.35],
        [None, 2.35, 0.0]])
    return settings, sc


class TestResolveSpeciesQueryUserScheme:
    """Regime A: the user grouped with -reduce, so the resolver matches
    that family at the user's sub_spec, reusing the stored per-atom
    descriptors (DESIGN 5.6.5)."""

    def test_user_scheme_uses_reduce_and_stored_query(self):
        stored = _si_query("si", 2.35)
        settings = types.SimpleNamespace(
            num_atoms=1, atom_element_id=[None, 1],
            atom_species_id=[None, 1],
            atom_reduce_fingerprint=[None, stored])
        matcher, sub_spec, query = makeinput._resolve_species_query(
            settings, None, ipdb, _si_fingerprinted_db(),
            element=1, species=1, elem_name="Si",
            user_scheme_active=True, reduce_sub_spec=_REDUCE_SUB_SPEC,
            file_reduce_descriptors={}, loen_site_cache={})
        assert matcher.name == "reduce"
        assert sub_spec is _REDUCE_SUB_SPEC
        assert query is stored                 # the stored per-atom fp


class TestResolveSpeciesQueryFileDictated:
    """Regime B: file-dictated species (no environment flag).  The
    database's preferred record selects the family and sub_spec, and the
    resolver computes the one query that family needs (DESIGN 5.6.5)."""

    def _preferred_reduce_db(self):
        # Two reduce-fingerprinted entries sharing one sub_spec; the
        # bulk entry is preferred (it sets the family/sub_spec), the
        # oxide entry is the default-tagged fallback.
        db = ElementDatabase(
            schema_version=2, element_symbol="Si", nuclear_z=14.0,
            nuclear_alpha=10.0, covalent_radius=1.1)
        db.potentials.append(_reduce_entry(
            "si-bulk", False, "si", 2.35, preferred=True))
        db.potentials.append(_reduce_entry(
            "si-oxide", True, "o", 1.60))
        return db

    def test_file_dictated_reduce_matches_crystalline_species(self):
        # The core C59 fix: a crystalline (file-dictated) Si species,
        # which no -reduce flag ever touched, still computes a reduce
        # query from geometry and matches the bulk entry instead of
        # silently taking the default.
        settings, sc = _si_dimer_settings()
        descriptor_cache = {}
        matcher, sub_spec, query = makeinput._resolve_species_query(
            settings, sc, ipdb, self._preferred_reduce_db(),
            element=1, species=1, elem_name="Si",
            user_scheme_active=False, reduce_sub_spec=None,
            file_reduce_descriptors=descriptor_cache,
            loen_site_cache={})
        assert matcher.name == "reduce"
        assert sub_spec == _REDUCE_SUB_SPEC
        # The computed query selects the bulk entry over the default.
        entry = makeinput._select_augmented_pot_entry(
            ipdb, self._preferred_reduce_db(), None, "Si",
            query=query, matcher=matcher, sub_spec=sub_spec,
            species=1)
        assert entry.label == "si-bulk"
        # The whole-structure shell codes were cached for reuse.
        assert len(descriptor_cache) == 1

    def test_no_preferred_record_yields_no_query(self):
        # An element whose database carries no fingerprints at all
        # offers no preferred record, so the resolver returns the
        # empty triple and the pick falls through to the default.
        settings, sc = _si_dimer_settings()
        result = makeinput._resolve_species_query(
            settings, sc, ipdb, _au_db(),
            element=1, species=1, elem_name="Au",
            user_scheme_active=False, reduce_sub_spec=None,
            file_reduce_descriptors={}, loen_site_cache={})
        assert result == (None, None, None)

    def test_file_dictated_bispectrum_uses_loen_seam(
            self, monkeypatch):
        # When the database's preferred record is bispectrum, the
        # resolver runs the loen seam (here stubbed) and summarizes the
        # sites of this (element, species) into a mean representative
        # that matches the bulk bispectrum entry.
        db = ElementDatabase(
            schema_version=2, element_symbol="Si", nuclear_z=14.0,
            nuclear_alpha=10.0, covalent_radius=1.1)
        db.potentials.append(_bispec_entry(
            "si-bulk-b", False, [0.1, 0.2, 0.3], preferred=True))
        db.potentials.append(_bispec_entry(
            "si-oxide-b", True, [5.0, 5.0, 5.0]))

        def _fake_loen(skeleton_path, sub_spec, **kwargs):
            return [LoenSite(
                site=1, element="Si", species=1, type_in_species=1,
                type_flat=1, vector=[0.1, 0.2, 0.3])]
        monkeypatch.setattr(
            makegroups, "loen_site_descriptors", _fake_loen)

        settings = types.SimpleNamespace(
            num_atoms=1, num_elements=1,
            atom_element_id=[None, 1], atom_species_id=[None, 1],
            atom_element_name=[None, "Si"])
        site_cache = {}
        matcher, sub_spec, query = makeinput._resolve_species_query(
            settings, None, ipdb, db, element=1, species=1,
            elem_name="Si", user_scheme_active=False,
            reduce_sub_spec=None, file_reduce_descriptors={},
            loen_site_cache=site_cache)
        assert matcher.name == "bispectrum"
        assert sub_spec == _BISPEC_SUB_SPEC
        assert query == pytest.approx([0.1, 0.2, 0.3])
        assert len(site_cache) == 1            # loen sites cached
        entry = makeinput._select_augmented_pot_entry(
            ipdb, db, None, "Si", query=query, matcher=matcher,
            sub_spec=sub_spec, species=1)
        assert entry.label == "si-bulk-b"


# ============================================================
#  -nofingerprint plumbing (DESIGN 5.6.5)
# ============================================================

class TestNoFingerprintFlag:
    """The -nofingerprint flag parses to ``no_fingerprint`` and is off
    by default, so fingerprint matching is enabled unless asked off."""

    def test_default_is_matching_enabled(self):
        settings = ScriptSettings.__new__(ScriptSettings)
        args = settings._args_from_options({})
        assert args.no_fingerprint is False

    def test_flag_sets_no_fingerprint(self):
        settings = ScriptSettings.__new__(ScriptSettings)
        args = settings._args_from_options({"no_fingerprint": True})
        assert args.no_fingerprint is True
