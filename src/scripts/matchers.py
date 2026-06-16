"""matchers.py -- The Phase-2 matcher protocol (ARCHITECTURE 8.9).

A *matcher* is a uniform interface over one local-environment descriptor
family -- the reduce concentric-shell scheme or the Fortran-side
bispectrum scheme.  This module is the single, neutral home for the
protocol so that every caller imports it *downward* and no import cycle
forms:

  - ``makeinput.py`` uses ``ReduceMatcher`` / ``ReduceStructureView`` for
    the in-process reduce grouping (``group_reduce``).
  - ``makegroups.py`` uses ``BispecMatcher`` for the sequential loen
    grouping flow (DESIGN 5.10).
  - ``build_initial_potentials.py`` and the manifest validator consult
    the ``MATCHERS`` registry.

The protocol first lived inside ``makeinput.py`` (the script that drove
all grouping).  It moved here once ``makegroups.py`` needed
``BispecMatcher``: leaving it in ``makeinput`` would have forced
``makeinput`` to import *upward* from ``makegroups`` (a cycle).  This is
also a stepping stone toward the eventual split in which ``makegroups``
owns all environment-based grouping and ``makeinput`` is a plain
input-writer.
"""

import math
from collections import Counter
from dataclasses import dataclass


# ===========================================================================
# Phase-2 matcher protocol (ARCHITECTURE 8.9, DESIGN 5.6.4)
# ===========================================================================
#
# A *matcher* is a uniform interface over one local-environment descriptor
# family.  Each matcher knows exactly one way to summarize "what an atom sees
# around it" into a comparable fingerprint, and exposes a small, fixed set of
# operations so the reduce grouping (makeinput), the makegroups bispectrum flow,
# and the initial-potential database producer (build_initial_potentials.py)
# can all drive it without caring
# which descriptor family it is.
#
# The two families Phase 2 ships with differ enormously in cost:
#
#   * ``reduce`` (ReduceMatcher) computes its fingerprint entirely in Python
#     from the structure geometry -- the concentric-shell description that the
#     long-standing reduce scheme has always used.  It needs no quantum run.
#   * ``bispectrum`` (BispecMatcher) computes its fingerprint inside the Imago
#     Fortran engine and is therefore far heavier; the makegroups
#     sequential loen flow (DESIGN 5.10) runs it.
#
# Adding a future family (e.g. SOAP) is a new Matcher subclass plus a new
# MATCHERS entry; no other code needs to change.


class Matcher:
    """Uniform interface over one local-environment descriptor family.

    Concrete subclasses (``ReduceMatcher``, ``BispecMatcher``) implement the
    operations below; this base class documents the protocol and supplies
    "not implemented" defaults so a partially-built matcher fails loudly
    rather than silently.

    Class attributes every matcher advertises:

    - ``name`` : the string written into a database entry's fingerprint
      ``method`` field and used as the ``MATCHERS`` registry key.
    - ``needs_loen_run`` : True when the descriptor is computed by the Imago
      Fortran engine (the loen path), False when it is computed in Python.
      This drives the nested-makeinput bootstrap of DESIGN 5.10.
    - ``default_similarity_floor`` : the per-matcher default distance
      threshold used when matching a query atom against stored database
      fingerprints (DESIGN 5.6.5).  A query whose nearest stored fingerprint
      is farther than this falls back to the default tag.  Users may override
      it per scheme on the command line.

    Instance operations:

    - ``compute_query(structure, sub_spec)`` : compute the per-atom
      fingerprints for ``structure`` under the parameters in ``sub_spec``.
    - ``distance(vec_a, vec_b)`` : a symmetric scalar distance between two
      fingerprints in this matcher's descriptor space.
    - ``representative(members)`` : reduce a list of member-atom fingerprints
      (the atoms of one species) into a single representative fingerprint.
    - ``to_loen_input(sub_spec)`` / ``parse_loen_output(path, sub_spec)`` :
      meaningful only when ``needs_loen_run`` is True; they translate a
      ``sub_spec`` into the Imago LOEN input block and parse the resulting
      ``fort.21`` descriptor file.
    """

    name = None
    needs_loen_run = False
    default_similarity_floor = None

    def compute_query(self, structure, sub_spec):
        """Compute one fingerprint per atom of ``structure``.

        ``structure`` exposes the atom-identity arrays and the
        minimum-distance matrix the matcher needs (``num_atoms``,
        ``atom_element_id``, ``atom_species_id``, ``atom_element_name``,
        ``min_dist``); ``StructureControl`` satisfies this directly and
        makeinput's internal ``ReduceStructureView`` mirrors it.
        """
        raise NotImplementedError(
            "compute_query is implemented by concrete Matcher subclasses")

    def distance(self, vec_a, vec_b):
        """Scalar distance between two fingerprints of this family."""
        raise NotImplementedError(
            "distance is implemented by concrete Matcher subclasses")

    def representative(self, members):
        """Reduce a species' member fingerprints to one representative."""
        raise NotImplementedError(
            "representative is implemented by concrete Matcher subclasses")

    def build_payload(self, fingerprint):
        """Serialize one fingerprint into the dict stored as a
        ``FingerprintRecord`` payload (DESIGN 5.2 / 5.4).  The producer's
        fingerprint harvest calls this; the consumer reads it back with
        ``extract_query_vector``.  Each matcher owns its payload field name
        (``values`` for bispectrum, ``shell_code`` for reduce)."""
        raise NotImplementedError(
            "build_payload is implemented by concrete Matcher subclasses")

    def extract_query_vector(self, payload):
        """Read a stored ``FingerprintRecord`` payload back into the form
        this matcher's ``distance`` expects -- the inverse of
        ``build_payload`` -- so the producer and consumer agree on the
        payload field name (DESIGN 5.4)."""
        raise NotImplementedError(
            "extract_query_vector is implemented by concrete subclasses")

    def to_loen_input(self, sub_spec):
        """Translate ``sub_spec`` into the Imago LOEN input block.

        Only loen-side matchers (``needs_loen_run`` True) implement this.
        """
        raise NotImplementedError(
            "to_loen_input is only meaningful for loen-side matchers")

    def parse_loen_output(self, path, sub_spec):
        """Parse a loen run's ``fort.21`` into per-site fingerprint vectors.

        Only loen-side matchers (``needs_loen_run`` True) implement this.
        """
        raise NotImplementedError(
            "parse_loen_output is only meaningful for loen-side matchers")


# ---------------------------------------------------------------------------
# The reduce descriptor: a per-atom shell code.
# ---------------------------------------------------------------------------

@dataclass
class ReduceShellLevel:
    """One concentric shell of the reduce descriptor for a single atom.

    Attributes
    ----------
    distance : float
        The level distance -- the distance (Angstrom) to the nearest atom
        that starts this shell.  This is the radial coordinate the reduce
        distance test compares within tolerance.
    members : list[tuple[int, int]]
        The ``(element_id, species_id)`` pair of every neighbor that falls in
        this shell, in ascending atom-index order.  Their multiset is what
        the reduce composition test compares; their count is the reduce count
        test.
    member_names : list[str]
        The element name of each member, parallel to ``members``.  Carried
        only so the ``reduceSummary`` diagnostic can list shell neighbors by
        name exactly as the historical algorithm did; it plays no part in the
        distance comparison.
    """

    distance: float
    members: list
    member_names: list


@dataclass
class ReduceShellCode:
    """The reduce fingerprint for one atom: its central element plus the
    per-level shell description (DESIGN 5.6.5, "shell_code").

    The comparison tolerance is embedded here so that ``distance`` is fully
    self-contained -- it needs no ``sub_spec`` argument, matching the matcher
    protocol's ``distance(vec_a, vec_b)`` signature (ARCHITECTURE 8.9).  All
    fingerprints in one reduce pass are computed under the same ``sub_spec``
    and therefore share the same tolerance, so reading it from the reference
    (first) operand reproduces the historical, directional reduce test where
    the tolerance band is scaled by the reduction atom's level distance.

    Attributes
    ----------
    element_id : int
        The central atom's element index.  Two atoms of different elements
        can never be the same species, so a mismatch is an infinite distance.
        This is a structure-local index, not a global atomic number, so it
        is used only for the within-structure ``distance`` test and never
        stored; the transferable form keeps the symbol instead.
    element_name : str
        The central atom's element symbol.  Carried for the stored,
        cross-structure payload (``build_payload``), where a global symbol is
        needed because the integer ``element_id`` does not transfer across
        structures.
    tolerance : float
        The fractional level-distance tolerance (e.g. 0.05 for 5%).
    levels : list[ReduceShellLevel]
        1-indexed list of shells (index 0 is an unused ``None`` placeholder),
        one entry per requested reduction level.
    """

    element_id: int
    element_name: str
    tolerance: float
    levels: list


@dataclass
class ReduceStructureView:
    """The slice of structure state ``ReduceMatcher.compute_query`` reads.

    During the makeinput species pass the authoritative per-atom species
    assignment lives on ``settings`` (it is mutated as earlier grouping flags
    run), while the minimum-distance geometry lives on the
    ``StructureControl``.  This view bundles the two so the matcher can read a
    single object whose attribute names match ``StructureControl`` -- which
    lets the producer (build_initial_potentials.py) hand a real
    ``StructureControl`` to the same ``compute_query`` unchanged.
    """

    num_atoms: int
    num_elements: int
    atom_element_id: list
    atom_species_id: list
    atom_element_name: list
    min_dist: list


class ReduceMatcher(Matcher):
    """Python-side matcher wrapping the established reduce scheme.

    The reduce fingerprint describes an atom by concentric spherical shells:
    for each requested level it records the distance to the shell and the
    element/species multiset of the neighbors in it (see ``compute_query``).
    Two atoms are "the same species" when, at every level, their shell
    distances agree within tolerance, they hold the same number of neighbors,
    and those neighbors form the same ``(element, species)`` multiset (see
    ``distance``).  This is exactly the test the historical ``group_reduce``
    applied in place; C54 simply moved it behind this class so the producer
    and the species pass reach it through one surface.
    """

    name = "reduce"
    needs_loen_run = False
    # Reduce distances are 0 (equivalent) or infinite (not), so any small
    # positive floor admits only exact matches.  0.05 is the DESIGN 5.6.5
    # starting value, tunable during the Phase-2 validation pass (C61).
    default_similarity_floor = 0.05

    def compute_query(self, structure, sub_spec):
        """Build the per-atom shell codes for ``structure``.

        Reproduces the historical reduce Phase 1 verbatim.  For each atom we
        sweep outward building concentric shells: at each level we find the
        closest still-unassigned atom (its distance defines the shell), then
        gather every still-unassigned atom whose distance falls in the band
        ``[shell_distance, shell_distance + thick]`` and within ``cutoff``.

        Parameters
        ----------
        structure : ReduceStructureView or StructureControl
            Supplies ``num_atoms``, the 1-indexed ``min_dist`` matrix, and the
            per-atom ``atom_element_id`` / ``atom_species_id`` /
            ``atom_element_name`` arrays.
        sub_spec : dict
            The reduce parameters: ``level`` (number of shells), ``thick``
            (shell acceptance band, Angstrom), ``cutoff`` (maximum neighbor
            distance, Angstrom), and ``tolerance`` (fractional level-distance
            tolerance, embedded into each returned fingerprint).

        Returns
        -------
        list[ReduceShellCode]
            1-indexed list (index 0 is ``None``) of per-atom fingerprints.
        """

        num_levels = sub_spec["level"]
        thick      = sub_spec["thick"]
        cutoff     = sub_spec["cutoff"]
        tolerance  = sub_spec["tolerance"]

        num_atoms    = structure.num_atoms
        min_dist     = structure.min_dist
        element_id   = structure.atom_element_id
        species_id   = structure.atom_species_id
        element_name = structure.atom_element_name

        fingerprints = [None] * (num_atoms + 1)

        for atom in range(1, num_atoms + 1):

            # atom_level[other] records which shell each other atom is
            # assigned to for the current atom: 0 = unassigned, -1 = the
            # atom itself (excluded from its own shells), >=1 = shell index.
            atom_level = [0] * (num_atoms + 1)
            atom_level[atom] = -1

            level_distance = [None] * (num_levels + 1)

            # Build each shell outward from the current atom.
            for level in range(1, num_levels + 1):

                # The closest still-unassigned atom starts this shell.  Ties
                # keep the earlier atom index (strict ``<``), matching the
                # historical scan order.
                closest_atom = 0
                for other in range(1, num_atoms + 1):
                    if atom_level[other] == 0:
                        if closest_atom == 0:
                            closest_atom = other
                        elif (min_dist[atom][other] <
                              min_dist[atom][closest_atom]):
                            closest_atom = other

                closest_dist = min_dist[atom][closest_atom]
                level_distance[level] = closest_dist

                # Sweep every atom whose distance falls in the shell band and
                # within the overall cutoff into this level.
                for other in range(1, num_atoms + 1):
                    d = min_dist[atom][other]
                    if (d >= closest_dist and
                            d <= closest_dist + thick and
                            d <= cutoff):
                        atom_level[other] = level

            # Collect each shell's members in ascending atom-index order,
            # recording their (element, species) pair for the comparison and
            # their element name for the diagnostic.
            levels = [None] * (num_levels + 1)
            for level in range(1, num_levels + 1):
                members = []
                member_names = []
                for other in range(1, num_atoms + 1):
                    if atom_level[other] == level:
                        members.append(
                            (element_id[other], species_id[other]))
                        member_names.append(element_name[other])
                levels[level] = ReduceShellLevel(
                    level_distance[level], members, member_names)

            fingerprints[atom] = ReduceShellCode(
                element_id[atom], element_name[atom],
                tolerance, levels)

        return fingerprints

    def distance(self, vec_a, vec_b):
        """Reduce distance: 0.0 when ``vec_b`` is equivalent to the reference
        ``vec_a``, ``math.inf`` the moment any equivalence test fails.

        The three historical tests are applied at every level; because the
        verdict is a conjunction over all levels and tests, the order in which
        they are checked does not change the result, only which test
        short-circuits.  ``vec_a`` is the reference (the species' reduction
        atom or representative): its tolerance and level distances set the
        acceptance bands, reproducing the directional historical comparison.
        """

        # Test 0: atoms of different elements can never be one species.
        if vec_a.element_id != vec_b.element_id:
            return math.inf

        num_levels = len(vec_a.levels) - 1
        for level in range(1, num_levels + 1):
            shell_a = vec_a.levels[level]
            shell_b = vec_b.levels[level]

            # Test 1: level distances within the reference's tolerance band.
            if (abs(shell_a.distance - shell_b.distance) >
                    vec_a.tolerance * shell_a.distance):
                return math.inf

            # Test 2: same number of neighbors in the shell.
            if len(shell_a.members) != len(shell_b.members):
                return math.inf

            # Test 3: identical (element, species) neighbor multiset.  With
            # the counts already equal (test 2), multiset equality is exactly
            # the historical greedy permutation match.
            if Counter(shell_a.members) != Counter(shell_b.members):
                return math.inf

        return 0.0

    def representative(self, members):
        """Return the first member's fingerprint as the species
        representative.  Intra-species reduce fingerprints agree within
        tolerance by construction, so any member is equally good (DESIGN
        5.6.5)."""

        return members[0]

    def build_payload(self, shell_code):
        """Serialize one ``ReduceShellCode`` into the element-only,
        cross-structure payload stored on disk (DESIGN 5.2): the central
        atom's element symbol plus, per level, the shell distance and the
        list of neighbor element symbols.

        The structure-local integer ids and the neighbor species are
        dropped here -- species numbering does not transfer across
        structures, so only the transferable element symbols are kept.  All
        symbols are lowercased to match the CLI element/species token
        convention.  ``levels`` is sliced from index 1 to skip the unused
        1-indexing placeholder."""

        levels = []
        for shell in shell_code.levels[1:]:
            levels.append({
                "distance": shell.distance,
                "neighbors": [name.lower()
                              for name in shell.member_names],
            })
        return {"shell_code": {
            "element": shell_code.element_name.lower(),
            "levels": levels,
        }}

    def extract_query_vector(self, payload):
        """Read the stored shell code back out of a record payload (DESIGN
        5.4).  Reduce records keep it under the ``shell_code`` key."""

        return payload["shell_code"]


@dataclass
class LoenSite:
    """One row of a loen ``fort.21`` file (DESIGN 5.10.3): a potential
    site's self-describing identity plus its bispectrum vector.

    The enriched ``fort.21`` leads each row with the site's identity so a
    consumer can map a row to its atom and type directly off the file,
    without a separate ``datSkl.map`` cross-reference.  ``vector`` holds
    the ``twoj2 + 1`` bispectrum components (the count of coupling channels
    ``j`` in the triangle range ``|j1 - j2| <= j <= j1 + j2``); the
    trailing accumulated-sum column is dropped.
    """

    site: int               # 1-based potential-site index (dat order)
    element: str            # element symbol of the site
    species: int            # per-element species index (the si1/si2 number)
    type_in_species: int    # per-element-species type index
    type_flat: int          # the single global (flat) type index
    vector: list            # the twoj2 + 1 bispectrum components


class BispecMatcher(Matcher):
    """Loen-side matcher for the bispectrum descriptor (ARCHITECTURE 8.9).

    The bispectrum captures *angular* detail about an atom's neighborhood
    -- how its neighbors are distributed over directions, not merely which
    radial shells they fall in -- and is computed by the Imago Fortran
    engine (the "loen" local-environment path), not in Python.  So
    ``needs_loen_run`` is True: producing a query vector means actually
    running the engine, which the ``makegroups`` sequential loen flow
    (DESIGN 5.10) does.

    This class supplies the parts that do *not* themselves run the engine:
    the parameter mapping (:meth:`to_loen_input`), the ``fort.21`` reader
    (:meth:`parse_loen_output`), the vector distance (:meth:`distance`),
    the species representative (:meth:`representative`), and the on-disk
    payload (de)serialization.  It does *not* implement
    :meth:`compute_query`: bispectrum vectors are not produced in-process,
    so the ``makegroups`` orchestrator runs the engine and reads the
    resulting ``fort.21`` through :meth:`parse_loen_output` instead.  The
    inherited ``compute_query`` raises ``NotImplementedError`` so any
    in-process call fails loudly.
    """

    name = "bispectrum"
    needs_loen_run = True
    # DESIGN 5.6.5 starting value for the bispectrum similarity floor.
    default_similarity_floor = 0.10

    def to_loen_input(self, sub_spec):
        """Map a bispectrum ``sub_spec`` to the LOEN parameter dict the
        ``LOEN_INPUT_DATA`` block of ``imago.dat`` expects (DESIGN 5.10.5).

        Both the producer (harvesting a fingerprint) and the ``makegroups``
        grouping flow call this, so the loen parameters they emit stay
        aligned by construction.  ``twoj1`` and ``twoj2`` are required --
        they set the angular-momentum pair and hence the output vector
        length ``twoj2 + 1`` (the count of coupling channels ``j`` in the
        triangle range ``|j1 - j2| <= j <= j1 + j2``, with ``twoj1 >=
        twoj2``).  The remaining three parameters are
        optional and default to the values ``makeinput.py`` has
        historically hardcoded, so an omitted key reproduces today's
        behavior exactly.

        Parameters
        ----------
        sub_spec : dict
            Must contain ``twoj1`` and ``twoj2`` (ints).  May contain
            ``max_neigh`` (int), ``cutoff`` (real, Bohr), and
            ``angle_squeeze`` (real).  An optional ``by_element`` key is
            reserved for the element-aware variant (TODO C62 / D10) and is
            rejected here until that lands.

        Returns
        -------
        dict
            Keys ``loenCode``, ``twoj1``, ``twoj2``, ``max_neigh``,
            ``cutoff``, ``angleSqueeze`` -- every parameter the LOEN block
            reads.
        """

        # Element-aware bispectrum changes the loen output shape and is not
        #   built yet; refuse it loudly rather than silently emitting the
        #   element-blind parameters.  C62 replaces this guard with the
        #   real `bispecByElement` handling (DESIGN 5.10.5).
        if sub_spec.get("by_element", False):
            raise NotImplementedError(
                "element-aware bispectrum (sub_spec by_element=true) is "
                "not yet implemented; see TODO C62 / D10")

        if "twoj1" not in sub_spec:
            raise ValueError(
                "BispecMatcher.to_loen_input requires sub_spec['twoj1']")
        if "twoj2" not in sub_spec:
            raise ValueError(
                "BispecMatcher.to_loen_input requires sub_spec['twoj2']")

        return {
            "loenCode":     1,
            "twoj1":        sub_spec["twoj1"],
            "twoj2":        sub_spec["twoj2"],
            "max_neigh":    sub_spec.get("max_neigh", 20),
            "cutoff":       sub_spec.get("cutoff", 5.0),
            "angleSqueeze": sub_spec.get("angle_squeeze", 0.85),
        }

    def parse_loen_output(self, path, sub_spec):
        """Read a loen ``fort.21`` file into per-site records (DESIGN
        5.10.3).

        ``fort.21`` opens with a header line, then carries one row per
        potential site of the whole structure in site-index order.  Each
        data row leads with the site's identity -- ``site#``, ``element``,
        ``species``, ``type_in_species``, ``type_flat`` -- then the
        ``twoj2 + 1`` real bispectrum components, then a trailing
        accumulated-sum column this matcher ignores.  The component count
        is ``twoj2 + 1``: the number of allowed coupling channels ``j`` in
        the triangle range ``|j1 - j2| <= j <= j1 + j2`` (the engine
        guarantees ``twoj1 >= twoj2``).

        Returning the identity alongside the vector lets the orchestrator
        (``makegroups``, DESIGN 5.10) map a row to its atom and type
        directly off the file, rather than re-deriving it through
        ``datSkl.map``.

        Parameters
        ----------
        path : str
            Path to the ``fort.21`` produced by the loen run.
        sub_spec : dict
            The same parameters used to produce the file; ``twoj2`` fixes
            the per-row component count ``twoj2 + 1``.

        Returns
        -------
        list[LoenSite]
            One record per potential site, in site-index order.
        """

        num_components = sub_spec["twoj2"] + 1
        # Five identity columns precede the components: site#, element,
        #   species, type_in_species, type_flat (DESIGN 5.10.3).
        num_identity = 5
        min_columns = num_identity + num_components

        sites = []
        with open(path) as fort21_file:
            for line_number, line in enumerate(fort21_file, start=1):
                tokens = line.split()
                if not tokens:
                    # Skip blank lines so a trailing newline does not
                    #   become a spurious zero-length site.
                    continue
                # The header row leads with "site#" (non-numeric); every
                #   data row leads with an integer site index.  Anything
                #   that does not start with an integer is the header (or
                #   stray text) and is skipped.
                try:
                    site_index = int(tokens[0])
                except ValueError:
                    continue
                if len(tokens) < min_columns:
                    raise ValueError(
                        f"{path}:{line_number}: expected at least "
                        f"{min_columns} columns (5 identity + "
                        f"{num_components} components, twoj2 + 1) but "
                        f"found {len(tokens)}")
                sites.append(LoenSite(
                    site=site_index,
                    element=tokens[1],
                    species=int(tokens[2]),
                    type_in_species=int(tokens[3]),
                    type_flat=int(tokens[4]),
                    vector=[float(token) for token
                            in tokens[5:5 + num_components]],
                ))
        return sites

    def distance(self, vector_a, vector_b):
        """Euclidean (L2) distance between two bispectrum vectors.

        Symmetric and cheap, and consistent with the element-wise mean
        used by :meth:`representative`: the mean is the point that
        minimizes the summed squared L2 distance to its members, so the
        representative sits at the center of the metric this distance
        defines.  The two vectors must share a length (the same
        ``sub_spec`` produced them); a mismatch is a hard error rather
        than a silently truncated comparison.
        """

        if len(vector_a) != len(vector_b):
            raise ValueError(
                f"bispectrum distance over vectors of unequal length "
                f"{len(vector_a)} and {len(vector_b)}; they must come "
                f"from the same sub_spec")
        return math.dist(vector_a, vector_b)

    def representative(self, members):
        """Collapse a species' member vectors into one representative: the
        element-wise arithmetic mean (DESIGN 5.6.5, ARCHITECTURE 8.9).

        Bispectrum species group atoms whose environments are *similar*,
        not identical, so the members scatter around a center; the mean is
        the order-independent summary that speaks for the whole group --
        no single atom is privileged.  All members share one length
        because they come from one bootstrap run under one ``sub_spec``,
        so the mean is well-defined slot by slot.
        """

        if not members:
            raise ValueError(
                "BispecMatcher.representative needs at least one member")
        num_members = len(members)
        num_components = len(members[0])
        return [sum(member[slot] for member in members) / num_members
                for slot in range(num_components)]

    def build_payload(self, vector):
        """Serialize one bispectrum vector into the stored record payload
        (DESIGN 5.2 / 5.4): the components live under the ``values`` key."""

        return {"values": list(vector)}

    def extract_query_vector(self, payload):
        """Read a stored bispectrum vector back out of a record payload
        (DESIGN 5.4).  Bispectrum records keep it under ``values``."""

        return payload["values"]


# The registry the species pass and the producer consult by manifest
# ``method`` name; ``initial_potential_db.load`` validates fingerprint
# methods against ``MATCHERS.keys()`` (ARCHITECTURE 8.9, rule 9).
MATCHERS = {
    "reduce":     ReduceMatcher,
    "bispectrum": BispecMatcher,
}
