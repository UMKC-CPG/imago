#!/usr/bin/env python3

"""makegroups.py -- Group the atoms of a non-crystalline structure into
species by their bispectrum local-environment fingerprints.

This is the orchestrator for the *Fortran-side* (bispectrum) matcher
described in DESIGN 5.10 and PSEUDOCODE 11.3.f.  A bispectrum
descriptor measures the angular arrangement of an atom's neighbors and
is computed inside the Imago engine's "loen" (local-environment) path,
not in Python.  Because that descriptor can only come from a *completed*
engine run, grouping by it is a short, explicit **sequence** of ordinary
program runs driven from *outside* makeinput.py -- never makeinput.py
launching a copy of itself:

  1. First makeinput.  Write a provisional imago.dat from the skeleton
     with no environment grouping (every atom keeps the species the
     skeleton already declares).  The LOEN_INPUT_DATA block of that
     imago.dat carries the chosen (twoj1, twoj2, ...) parameters, passed
     through makeinput's -loeninput flag (DESIGN 5.10.2 / 5.10.5).
  2. Run loen.  Invoke ``imago.py -loen -scf no`` against that imago.dat.
     The -scf no skips the self-consistent field cycle the bispectrum
     does not need; the run writes a self-describing fort.21 (one row per
     potential site, each row leading with its own identity -- site#,
     element, species, type_in_species, type_flat -- then the descriptor;
     DESIGN 5.10.3).
  3. Bucket and rewrite.  Read fort.21, group the atoms of each element
     by how close their fingerprints are (the shared leader-clustering of
     matchers.bucket_by_fingerprint, DESIGN 5.6.4), and rewrite the
     skeleton with explicit per-element species tags -- Si1, Si2, ... then
     O1, O2, ... each element restarting at 1 (DESIGN 5.10.4).  makeinput
     then reads those tags on a later run like any other explicit
     assignment.

Grouping applies only to **non-crystalline (P1)** systems.  A
symmetry-bearing skeleton is refused up front: rewriting its types in P1
would force it out of its space group, and the Brillouin-zone k-point
folding the band structure depends on samples the irreducible wedge
*using that space group*.  Dropping it would silently corrupt the
result, so a crystal is harvested on the witness path instead -- it is
never grouped (DESIGN 5.10.1).

The module is dual-mode:

  * Importable.  ``group_by_bispectrum(skeleton_path, sub_spec,
    similarity_floor)`` is what the producer (build_initial_potentials.py)
    calls; it returns the per-atom species assignment and leaves the
    rewritten skeleton on disk, backed up alongside the original.
  * Command line.  Running ``makegroups.py <skeleton> [options]`` does
    the same for manual use.

There is no recursion to guard against: each step is an ordinary process
run in order, and the only state passed between the two makeinput runs is
the skeleton file itself.
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from collections import OrderedDict

from matchers import MATCHERS, bucket_by_fingerprint
from structure_control import StructureControl


# The sibling scripts this orchestrator runs live next to it in the
#   source tree and, once installed, in the same bin/ directory.  We
#   resolve them relative to this file first so a checkout runs against
#   its own copies, and fall back to PATH for an installed layout.
_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))


# ===========================================================================
# Pure helpers (no engine run): the P1 guard, the species bucketing, and
# the skeleton rewrite.  These carry the grouping logic and are unit-tested
# directly, without the Fortran binary.
# ===========================================================================


def _parse_atom_tag(tag):
    """Split a skeleton atom tag into its element symbol and species
    number.

    A tag is an element symbol followed by an integer species id, e.g.
    ``Si2`` or ``o1`` (DESIGN 5.10.4; the skeleton format of
    structure_control.read_imago_skl).  A tag with no trailing digit
    means species 1, matching the engine's own convention of appending a
    ``1`` when none is written.  The element's letter case is returned as
    written so the rewrite can preserve it; callers that compare elements
    lowercase both sides.

    Parameters
    ----------
    tag : str
        The atom tag token, e.g. ``"Si2"``.

    Returns
    -------
    tuple(str, int)
        ``(element_symbol, species_number)``.
    """

    match = re.match(r'^([A-Za-z]+)(\d*)$', tag.strip())
    if not match:
        raise ValueError(
            f"cannot parse atom tag '{tag}': expected an element symbol "
            "optionally followed by an integer species id (e.g. 'Si2')")
    element = match.group(1)
    species = int(match.group(2)) if match.group(2) else 1
    return element, species


def _require_p1(structure, skeleton_path):
    """Refuse to group anything but a P1 (non-crystalline) skeleton.

    Grouping reassigns every atom a type from its fingerprint, which is
    only meaningful when the cell carries no symmetry: it must first drop
    the space group and treat the cell as P1.  Doing that to a real
    crystal would discard the space group the k-point folding relies on
    and corrupt the band structure, so a symmetry-bearing skeleton is
    rejected here rather than silently regrouped (DESIGN 5.10.1, 5.10.4).
    A crystalline reference belongs on the witness path, which leaves both
    the symmetry and the skeleton untouched.

    The skeleton is P1 when its space group resolves to number 1 (the
    mandatory ``space`` line reads ``1_a``) *and* its supercell is the
    single unit cell ``1 1 1``.

    Parameters
    ----------
    structure : StructureControl
        A structure already populated by ``read_imago_skl``.
    skeleton_path : str
        Path to the skeleton, used only for the error message.
    """

    # structure.supercell is 1-indexed (slot 0 is an unused placeholder),
    #   so the three real repeat counts are slots 1..3.
    supercell = structure.supercell[1:4]
    if structure.space_group_num != 1 or supercell != [1, 1, 1]:
        raise ValueError(
            f"{skeleton_path} is not P1 (space group "
            f"'{structure.space_group}', number "
            f"{structure.space_group_num}; supercell {supercell}); "
            "bispectrum grouping rewrites every atom's type in P1 and "
            "would drop the space group the Brillouin-zone k-point "
            "folding depends on.  Group only non-crystalline (P1, the "
            "'1_a' space line and a '1 1 1' supercell) skeletons; harvest "
            "a crystal on the witness path instead (DESIGN 5.10.1).")


def _assign_species_from_rows(rows, matcher, similarity_floor):
    """Bucket each element's atoms into species by fingerprint distance.

    The fort.21 rows arrive in engine ("dat") order, which is sorted by
    element.  The bispectrum distance is element-blind, so atoms of
    different elements must never share a species: we partition the rows
    by element and cluster each element's fingerprints on its own with the
    shared leader-clustering of ``bucket_by_fingerprint`` (DESIGN 5.6.4).

    Each bucket becomes a species, numbered from 1 within its element in
    first-appearance order (DESIGN 5.10.4).  The result keys the new
    species by the atom's identity in the *first, ungrouped* pass --
    ``(element, original_species)`` -- which is exactly what the skeleton
    atom tags carry, so the rewrite can look each atom up directly off its
    tag without a separate cross-reference file.

    Parameters
    ----------
    rows : list[LoenSite]
        The fort.21 records, as returned by
        ``BispecMatcher.parse_loen_output``.
    matcher : Matcher
        Supplies the ``distance`` / ``representative`` the bucketing uses.
    similarity_floor : float
        The largest distance at which two atoms still share a species.

    Returns
    -------
    dict
        ``{(element_lower, original_species): new_species}``, the species
        each first-pass atom is reassigned to.
    """

    # Group rows by element, preserving the order each element is first
    #   seen so the species numbering is deterministic.
    rows_by_element = OrderedDict()
    for row in rows:
        rows_by_element.setdefault(row.element.lower(), []).append(row)

    new_species_by_key = {}
    for element, element_rows in rows_by_element.items():
        vectors = [row.vector for row in element_rows]
        bucket_of = bucket_by_fingerprint(
            vectors, matcher, similarity_floor)
        for row, bucket_index in zip(element_rows, bucket_of):
            # bucket_of is 0-based; species tags start at 1 per element.
            new_species_by_key[(element, row.species)] = bucket_index + 1
    return new_species_by_key


def _retag_atom_line(line, new_species_by_key):
    """Rewrite one atom line's tag with its newly assigned species.

    Only the tag token (the first field) is replaced; everything after it
    -- the original spacing, the coordinates, an optional connection tag,
    a trailing comment, the newline -- is preserved verbatim so the
    rewrite is a faithful, minimal edit of the skeleton.  The element's
    letter case is kept as written; only the species number changes.

    Parameters
    ----------
    line : str
        The original atom line.
    new_species_by_key : dict
        ``(element_lower, original_species) -> new_species`` from
        ``_assign_species_from_rows``.

    Returns
    -------
    tuple(str, int)
        The rewritten line and the new species number it now carries.
    """

    stripped = line.lstrip()
    lead = line[:len(line) - len(stripped)]
    tag_token = stripped.split(None, 1)[0]
    # Everything from the end of the tag onward (whitespace + coordinates
    #   + any connection tag/comment + newline) is carried through as-is.
    remainder = stripped[len(tag_token):]

    element, original_species = _parse_atom_tag(tag_token)
    key = (element.lower(), original_species)
    if key not in new_species_by_key:
        raise ValueError(
            f"atom tag '{tag_token}' (element '{element}', species "
            f"{original_species}) has no bispectrum fingerprint row.  The "
            "first, ungrouped makeinput/loen pass must give every atom "
            "its own species so each produces a fingerprint row in "
            "fort.21 (DESIGN 5.10.2).")
    new_species = new_species_by_key[key]
    new_tag = element + str(new_species)
    return lead + new_tag + remainder, new_species


def _rewrite_skeleton_species(skeleton_path, new_species_by_key):
    """Rewrite a skeleton's atom tags in place with their new species.

    The skeleton carries a single ``frac N`` / ``cart N`` block followed
    by N atom lines (structure_control.read_imago_skl).  We retag each of
    those N lines and leave every other line -- title, cell, space group,
    supercell -- untouched.  The per-element species restart at 1 (the
    numbering ``_assign_species_from_rows`` produced), the convention the
    engine expects when it expands the per-element species into its flat
    type list (DESIGN 5.10.4).

    Parameters
    ----------
    skeleton_path : str
        Path to the skeleton; rewritten in place.
    new_species_by_key : dict
        ``(element_lower, original_species) -> new_species``.

    Returns
    -------
    list[int]
        ``species_of``: the new species of each atom, in skeleton (file)
        order.
    """

    with open(skeleton_path) as handle:
        lines = handle.readlines()

    rewritten = list(lines)
    species_of = []
    for index, line in enumerate(lines):
        tokens = line.split()
        # The atom list opens with a 'frac'/'cart' keyword and a count;
        #   the N atom lines follow immediately (no blank lines between),
        #   matching how read_imago_skl indexes them.
        if (tokens
                and (tokens[0].lower().startswith('frac')
                     or tokens[0].lower().startswith('cart'))
                and len(tokens) >= 2):
            num_atoms = int(tokens[1])
            for atom_offset in range(1, num_atoms + 1):
                atom_index = index + atom_offset
                rewritten[atom_index], new_species = _retag_atom_line(
                    lines[atom_index], new_species_by_key)
                species_of.append(new_species)
            # A skeleton holds exactly one atom block; stop after it.
            break

    if not species_of:
        raise ValueError(
            f"{skeleton_path} has no 'frac'/'cart' atom block to "
            "regroup; is it a valid imago skeleton?")

    with open(skeleton_path, 'w') as handle:
        handle.writelines(rewritten)
    return species_of


# ===========================================================================
# Orchestration (runs the engine): the two-process sequence that turns a
# P1 skeleton into a bispectrum-grouped one.
# ===========================================================================


def _resolve_sibling(name):
    """Find a sibling script (makeinput.py / imago.py) to run.

    Prefer the copy next to this file (a source checkout runs against its
    own scripts); fall back to PATH for an installed bin/ layout.
    """

    candidate = os.path.join(_SCRIPTS_DIR, name)
    if os.path.exists(candidate):
        return candidate
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(
        f"could not find '{name}' next to makegroups.py ({_SCRIPTS_DIR}) "
        "or on PATH; makegroups needs both makeinput.py and imago.py to "
        "run the loen sequence.")


def _loen_input_values(matcher, sub_spec):
    """Turn a bispectrum sub_spec into the six -loeninput values.

    ``to_loen_input`` is the single mapping both the producer and this
    orchestrator share (DESIGN 5.10.5), so the LOEN block makeinput writes
    matches the descriptor contract.  The values are emitted in
    LOEN-block order -- method code, the 2j1/2j2 pair, max_neigh, cutoff,
    angleSqueeze -- as strings ready for the makeinput command line.
    """

    params = matcher.to_loen_input(sub_spec)
    return [str(params['loenCode']), str(params['twoj1']),
            str(params['twoj2']), str(params['max_neigh']),
            str(params['cutoff']), str(params['angleSqueeze'])]


def _run_makeinput(work_dir, loen_values):
    """Run the first, ungrouped makeinput in the working directory.

    makeinput reads its skeleton from the fixed name ``imago.skl`` in the
    current directory (it takes no skeleton argument), which is why the
    skeleton is staged there under that name and the run is launched with
    ``cwd=work_dir``.  ``-loeninput`` fills the LOEN_INPUT_DATA block from
    the chosen parameters; no ``-pot`` is given because the bispectrum is
    purely geometric, so the potential each atom receives is irrelevant
    (DESIGN 5.10.2).
    """

    command = [sys.executable, _resolve_sibling('makeinput.py'),
               '-loeninput', *loen_values]
    subprocess.run(command, cwd=work_dir, check=True)


def _run_loen(work_dir):
    """Run the loen pass (``imago.py -loen -scf no``) on the prepared dir.

    ``-scf no`` skips the self-consistent field cycle the bispectrum does
    not need: loen reads only the structure and the LOEN block.  The run
    writes the self-describing fort.21 (DESIGN 5.10.3).
    """

    command = [sys.executable, _resolve_sibling('imago.py'),
               '-loen', '-scf', 'no']
    subprocess.run(command, cwd=work_dir, check=True)


def _find_loen_descriptor(work_dir):
    """Locate the bispectrum descriptor file the loen run produced.

    The loen pass writes its descriptor table to Fortran unit 21, which
    imago.py's output management then renames to the loen "plot"
    (1D-profile) output ``<edge>_loen<basis>.plot`` -- e.g.
    ``gs_loen-fb.plot``.  That single file is what makegroups reads: it
    is the self-describing table parse_loen_output expects (DESIGN
    5.10.3).  If it is absent the loen run did not complete normally, and
    that is a real failure to surface rather than paper over with
    fallback names -- so we raise.  (A spin-polarized run can leave more
    than one match; any is fine since the bispectrum is spin-independent
    geometry, so the first in sorted order is taken.)
    """

    matches = sorted(glob.glob(os.path.join(work_dir, '*loen*.plot')))
    if not matches:
        raise FileNotFoundError(
            f"no loen descriptor (*loen*.plot) in {work_dir}; the loen "
            "run did not complete -- check the run's loen .out log.")
    return matches[0]


def group_by_bispectrum(skeleton_path, sub_spec, similarity_floor=None,
                        work_dir=None, keep_work_dir=False):
    """Group a P1 skeleton's atoms into species by their bispectrum.

    This is the importable entry point the producer calls (DESIGN 5.10,
    PSEUDOCODE 11.3.f).  It guards that the skeleton is non-crystalline,
    runs the ungrouped makeinput + loen sequence in a scratch working
    directory, buckets the resulting fingerprints into per-element
    species, backs up the original skeleton, and rewrites it in place with
    explicit species tags.

    Parameters
    ----------
    skeleton_path : str
        Path to the P1 ``imago.skl`` to regroup.  Rewritten in place; the
        original is preserved next to it as ``<skeleton>.orig``.
    sub_spec : dict
        The bispectrum parameters: ``twoj1`` and ``twoj2`` are required,
        ``max_neigh``, ``cutoff``, and ``angle_squeeze`` default to the
        descriptor-contract values (DESIGN 5.10.5).
    similarity_floor : float, optional
        The largest fingerprint distance at which two atoms still share a
        species.  Defaults to the matcher's
        ``default_similarity_floor`` (DESIGN 5.6.5).
    work_dir : str, optional
        Scratch directory for the makeinput/loen runs.  Defaults to a
        ``makegroups_work`` directory beside the skeleton.
    keep_work_dir : bool, optional
        When True the scratch directory is left in place for inspection;
        by default it is removed once the grouping is read out.

    Returns
    -------
    list[int]
        ``species_of``: the new species of each atom, in skeleton order.
    """

    matcher = MATCHERS["bispectrum"]()
    if similarity_floor is None:
        similarity_floor = matcher.default_similarity_floor

    skeleton_path = os.path.abspath(skeleton_path)

    # 1. P1 guard: read the skeleton and refuse anything but a P1 cell.
    structure = StructureControl()
    structure.read_imago_skl(skeleton_path, use_file_species=True)
    _require_p1(structure, skeleton_path)

    # 2. Stage the skeleton in a scratch directory as imago.skl so the
    #    makeinput/loen runs never touch the user's working tree.
    if work_dir is None:
        work_dir = os.path.join(os.path.dirname(skeleton_path),
                                'makegroups_work')
    os.makedirs(work_dir, exist_ok=True)
    shutil.copy(skeleton_path, os.path.join(work_dir, 'imago.skl'))

    # 3-4. First makeinput (ungrouped, carrying the LOEN params), then the
    #      loen pass that computes the fingerprints.
    _run_makeinput(work_dir, _loen_input_values(matcher, sub_spec))
    _run_loen(work_dir)

    # 5. Read the self-describing loen descriptor into per-site records.
    rows = matcher.parse_loen_output(
        _find_loen_descriptor(work_dir), sub_spec)

    # 6. Bucket each element's fingerprints into species.
    new_species_by_key = _assign_species_from_rows(
        rows, matcher, similarity_floor)

    # 7. Back up the original skeleton, then rewrite it in place with the
    #    new per-element species tags.
    shutil.copy(skeleton_path, skeleton_path + '.orig')
    species_of = _rewrite_skeleton_species(
        skeleton_path, new_species_by_key)

    if not keep_work_dir:
        shutil.rmtree(work_dir, ignore_errors=True)
    return species_of


# ===========================================================================
# Command-line interface
# ===========================================================================


def _build_parser():
    """Build the argument parser for the ``makegroups.py`` CLI."""

    parser = argparse.ArgumentParser(
        description="Group the atoms of a non-crystalline (P1) imago "
                    "skeleton into species by their bispectrum "
                    "local-environment fingerprints.  Runs an ungrouped "
                    "makeinput plus a loen pass, buckets the atoms of "
                    "each element by fingerprint similarity, and rewrites "
                    "the skeleton in place with explicit per-element "
                    "species tags (the original is kept as "
                    "<skeleton>.orig).  A symmetry-bearing skeleton is "
                    "refused, since grouping only applies in P1.")
    parser.add_argument(
        'skeleton',
        help="Path to the P1 imago skeleton (.skl) to regroup.")
    parser.add_argument(
        '-twoj1', dest='twoj1', type=int, default=4,
        help="Larger angular-momentum parameter 2j1.  Default: 4.")
    parser.add_argument(
        '-twoj2', dest='twoj2', type=int, default=4,
        help="Smaller angular-momentum parameter 2j2; the descriptor has "
             "2j2 + 1 components.  Default: 4.")
    parser.add_argument(
        '-maxneigh', dest='maxneigh', type=int, default=50,
        help="Per-site neighbor-list cap.  Default: 50.")
    parser.add_argument(
        '-cutoff', dest='cutoff', type=float, default=9.0,
        help="Radial neighbor cutoff in Bohr.  Default: 9.0.")
    parser.add_argument(
        '-anglesqueeze', dest='anglesqueeze', type=float, default=0.85,
        help="Angular weighting (angleSqueeze).  Default: 0.85.")
    parser.add_argument(
        '-floor', dest='floor', type=float, default=None,
        help="Similarity floor: the largest fingerprint distance at "
             "which two atoms still share a species.  Default: the "
             "bispectrum matcher's built-in value (0.10).")
    parser.add_argument(
        '-workdir', dest='workdir', type=str, default=None,
        help="Scratch directory for the makeinput/loen runs.  Default: a "
             "'makegroups_work' directory beside the skeleton.")
    parser.add_argument(
        '-keepwork', dest='keepwork', action='store_true', default=False,
        help="Keep the scratch directory for inspection instead of "
             "removing it.  Default: remove it.")
    return parser


def main(argv=None):
    """Command-line entry point: group one skeleton and report the result.

    Parses the command line (or ``argv`` when given for testing), runs
    ``group_by_bispectrum``, and prints a short summary of how many
    species each element was grouped into.
    """

    args = _build_parser().parse_args(argv)
    sub_spec = {
        'twoj1':         args.twoj1,
        'twoj2':         args.twoj2,
        'max_neigh':     args.maxneigh,
        'cutoff':        args.cutoff,
        'angle_squeeze': args.anglesqueeze,
    }
    species_of = group_by_bispectrum(
        args.skeleton, sub_spec, similarity_floor=args.floor,
        work_dir=args.workdir, keep_work_dir=args.keepwork)
    num_species = len(set(species_of))
    print(f"Grouped {len(species_of)} atoms into {num_species} "
          f"species; rewrote {args.skeleton} (original kept as "
          f"{args.skeleton}.orig).")


if __name__ == '__main__':
    main()
