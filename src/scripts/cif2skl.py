#!/usr/bin/env python3

"""cif2skl.py -- Convert a CIF crystal structure into an imago.skl file,
preserving the space group.

A CIF (Crystallographic Information File) is the standard interchange
format for crystal structures -- what the Crystallography Open Database
serves and what journals publish.  Imago's own structure format is the
``imago.skl`` skeleton.  This tool bridges the two.

Crucially, an ``imago.skl`` is *not* a flat P1 list of every atom: it
stores the **asymmetric unit** plus a **space-group token**, and
``read_imago_skl``'s ``apply_space_group`` regenerates the full cell from
them.  Preserving the space group is a correctness requirement, not a
nicety -- Imago's Brillouin-zone integration samples the irreducible
wedge *using that space group*, so flattening a crystal to P1 would
silently corrupt every k-point-dependent result.  This converter
therefore recovers the CIF's space group rather than discarding it.

How it works (ARCHITECTURE 9.5):

1. ASE parses the CIF.  This is the only step that touches ``ase`` and it
   spares us a hand-written CIF parser: ASE handles CIF's flexible
   ``loop_`` blocks and symmetry expansion.  We read off the *authored
   asymmetric unit* (the raw ``_atom_site_*`` list), the cell, the
   International Tables space-group number and origin setting, and -- by
   reading the ordinary way -- the fully symmetry-expanded cell that
   serves as the verification reference.

2. The space group is resolved by **verification, not by parsing
   operation tables**.  The ``spaceDB`` settings for the CIF's IT number
   are the candidate tokens (``<IT#>_a``, ``<IT#>_b``, ...).  For each,
   the authored asymmetric unit is written with that token and run
   through ``apply_space_group``; the variant whose regenerated cell
   reproduces ASE's full expansion (atom count and positions, up to
   lattice translation and permutation, within a tolerance) is the
   answer.  Imago's own symmetry engine is the judge, so no ``spaceDB``
   operation set is ever parsed, and the origin choice falls out of the
   match because the symmetry translations are origin-dependent.

3. The skeleton is written as the asymmetric unit plus the matched token.
   If no candidate verifies it is a hard error (P1 fallback is *not*
   offered for a crystal); ``--space <token>`` forces a specific setting,
   still verified.

No symmetry library (e.g. spglib) is used: the CIF's symmetry is taken as
authored (ASE reports it) and checked against Imago's expansion, never
re-derived.  Partial / mixed site occupancy is refused for now.

Command-line usage::

    cif2skl.py input.cif                 # writes input.skl beside it
    cif2skl.py input.cif output.skl      # writes the named skeleton
    cif2skl.py input.cif --space 227_b   # force a spaceDB setting
    cif2skl.py input.cif -t "Si diamond" # set the skeleton title

Importable usage::

    from cif2skl import convert
    convert("au_fcc.cif", "au_fcc.skl")
"""

import argparse
import glob
import os
import sys
from datetime import datetime

from structure_control import StructureControl

# Fractional-coordinate tolerance for matching the regenerated cell
#   against ASE's expansion.  Loose enough for the rounding in published
#   CIF coordinates, tight enough to separate distinct origin settings.
_MATCH_TOLERANCE = 1.0e-3

# How far below 1.0 a refined site occupancy may sit and still count as
#   fully occupied.  Experimental CIFs report a full site as a value
#   very close to -- but not exactly -- 1 (e.g. 0.9999 is a rounding;
#   the site is physically full).  Anything more than this far below 1
#   is treated as genuine partial occupancy (disorder or vacancies),
#   which the skeleton format cannot express and so is refused.
_FULL_OCCUPANCY_TOLERANCE = 1.0e-2


class CifConversionError(Exception):
    """A CIF could not be converted -- unreadable, unsupported (partial
    occupancy), or its space group could not be resolved to a spaceDB
    setting.  Carries a message naming the specific failure."""


def _element_from_label(token):
    """Extract an element symbol from a CIF ``_atom_site`` label or type
    symbol (e.g. ``"Si1"`` -> ``"Si"``, ``"O2-"`` -> ``"O"``).

    The two-letter prefix is tried first so two-letter elements
    (``"Ca"``, ``"Cl"``) win over their one-letter prefixes; the
    one-letter prefix is the fallback.  Both are validated against ASE's
    element table so a stray label fails loudly rather than inventing an
    element."""

    from ase.data import chemical_symbols

    text = token.strip()
    two = text[:2].title()
    if len(text) >= 2 and two in chemical_symbols:
        return two
    one = text[:1].title()
    if one in chemical_symbols:
        return one
    raise CifConversionError(
        f"cannot read an element symbol from CIF atom token {token!r}")


def _read_cif(cif_path):
    """Parse ``cif_path`` with ASE.

    Returns a tuple of everything the resolver needs: the six cell
    parameters, the authored asymmetric unit (element symbols and
    fractional coordinates), the IT space-group number, ASE's origin
    setting, and the fully symmetry-expanded reference (lowercased
    element symbols and fractional coordinates).
    """

    try:
        from ase.io import read as ase_read
        from ase.io import cif as ase_cif
    except ImportError as exc:
        raise CifConversionError(
            "cif2skl needs the ASE package to read CIF files "
            "(`pip install ase`); the rest of the toolchain does not "
            "depend on it.") from exc

    # ----- Authored asymmetric unit, straight off the CIF tags.
    with open(cif_path) as handle:
        blocks = list(ase_cif.parse_cif(handle))
    if not blocks:
        raise CifConversionError(
            f"{cif_path}: no data block found in the CIF")
    block = blocks[0]

    # A refined occupancy within _FULL_OCCUPANCY_TOLERANCE of 1 (e.g.
    #   0.9999) is a fully-occupied site reported with rounding; only
    #   occupancy genuinely below that band is partial (disorder or
    #   vacancies), which the skeleton format cannot represent.
    occupancy = block.get("_atom_site_occupancy")
    if occupancy is not None and any(
            float(value) < 1.0 - _FULL_OCCUPANCY_TOLERANCE
            for value in occupancy):
        raise CifConversionError(
            f"{cif_path}: partial site occupancy is not supported yet; "
            f"convert a fully-occupied structure or edit the CIF")

    labels = block.get("_atom_site_label")
    types = block.get("_atom_site_type_symbol")
    frac_x = block.get("_atom_site_fract_x")
    frac_y = block.get("_atom_site_fract_y")
    frac_z = block.get("_atom_site_fract_z")
    if not labels or frac_x is None:
        raise CifConversionError(
            f"{cif_path}: CIF carries no _atom_site fractional "
            f"coordinates to read")
    asym_symbols = [
        _element_from_label((types or labels)[i])
        for i in range(len(labels))]
    asym_frac = [
        (float(frac_x[i]), float(frac_y[i]), float(frac_z[i]))
        for i in range(len(labels))]

    it_number = (block.get("_space_group_it_number")
                 or block.get("_symmetry_int_tables_number"))

    # ----- Fully expanded reference cell (ASE applies the symmetry).
    atoms = ase_read(cif_path, format="cif")
    spacegroup = atoms.info.get("spacegroup")
    setting = getattr(spacegroup, "setting", None)
    if it_number is None and spacegroup is not None:
        it_number = spacegroup.no
    cell = list(atoms.cell.cellpar())
    full_symbols = [s.lower() for s in atoms.get_chemical_symbols()]
    full_frac = atoms.get_scaled_positions()

    return (cell, asym_symbols, asym_frac,
            int(it_number) if it_number is not None else None,
            setting, full_symbols, full_frac)


def _candidate_tokens(it_number, setting, spacedb_dir):
    """List the spaceDB tokens for ``it_number`` (``<IT#>_a`` ... and a
    bare ``<IT#>`` if present), ordered so the setting ASE reported is
    tried first."""

    tokens = sorted(
        os.path.basename(path)
        for path in glob.glob(os.path.join(spacedb_dir,
                                           f"{it_number}_*")))
    bare = str(it_number)
    if os.path.exists(os.path.join(spacedb_dir, bare)):
        tokens.insert(0, bare)
    # Origin setting 2 in the CIF most often corresponds to the "_b"
    #   variant; try it first when present so the common case verifies on
    #   the first attempt (verification still decides correctness).
    if setting == 2:
        preferred = f"{it_number}_b"
        if preferred in tokens:
            tokens.remove(preferred)
            tokens.insert(0, preferred)
    return tokens


def _skeleton_lines(cell, symbols, frac, space_token, title):
    """Build the imago.skl text (asymmetric unit + space token)."""

    a, b, c, alpha, beta, gamma = cell
    lines = [
        "title\n", title + "\n", "end\n",
        "cell\n", f"{a} {b} {c} {alpha} {beta} {gamma}\n",
        f"fract {len(symbols)}\n",
    ]
    for index in range(len(symbols)):
        x, y, z = frac[index]
        lines.append(
            f"{symbols[index]} {x:.10f} {y:.10f} {z:.10f}\n")
    lines += [f"space {space_token}\n", "supercell 1 1 1\n", "full\n"]
    return lines


def _expansion_matches(full_symbols, full_frac, structure):
    """True when the expanded ``structure`` reproduces the ASE reference.

    Each reference atom must pair with a same-element atom in the
    structure whose fractional coordinates agree within tolerance, with
    wraparound (0 and 1 are the same point) and one-to-one matching."""

    if structure.num_atoms != len(full_symbols):
        return False
    imago = [
        (structure.atom_element_name[atom].lower(),
         structure.fract_abc[atom][1:4])
        for atom in range(1, structure.num_atoms + 1)]
    used = [False] * len(imago)
    for symbol, position in zip(full_symbols, full_frac):
        matched = False
        for index, (other_symbol, other_pos) in enumerate(imago):
            if used[index] or other_symbol != symbol:
                continue
            if all(_coord_close(position[axis], other_pos[axis])
                   for axis in range(3)):
                used[index] = True
                matched = True
                break
        if not matched:
            return False
    return True


def _coord_close(value_a, value_b):
    """Fractional-coordinate equality with wraparound (0 == 1)."""

    difference = abs((value_a - value_b) % 1.0)
    return min(difference, 1.0 - difference) < _MATCH_TOLERANCE


def cif_to_skeleton(cif_path, space_override=None, title=None):
    """Resolve ``cif_path``'s space group and return ``(skeleton_lines,
    space_token)`` for the asymmetric-unit skeleton.

    Raises :class:`CifConversionError` when no spaceDB setting reproduces
    the CIF's expansion (or when the forced ``space_override`` does not).
    """

    (cell, asym_symbols, asym_frac, it_number,
     setting, full_symbols, full_frac) = _read_cif(cif_path)

    if title is None:
        title = (f"Converted from {os.path.basename(cif_path)} by "
                 f"cif2skl (space group preserved).")

    spacedb_dir = StructureControl().space_group_db
    if space_override is not None:
        candidates = [space_override]
    else:
        if it_number is None:
            raise CifConversionError(
                f"{cif_path}: CIF states no space-group number; pass "
                f"--space <token> to choose a spaceDB setting")
        candidates = _candidate_tokens(it_number, setting, spacedb_dir)
        if not candidates:
            raise CifConversionError(
                f"{cif_path}: no spaceDB setting found for IT number "
                f"{it_number}; pass --space <token>")

    for token in candidates:
        lines = _skeleton_lines(
            cell, asym_symbols, asym_frac, token, title)
        structure = StructureControl()
        try:
            structure.read_imago_skl(lines=lines)
        except Exception:           # noqa: BLE001 -- a bad token is a
            continue                #   miss, not a crash; try the next
        if _expansion_matches(full_symbols, full_frac, structure):
            return lines, token

    if space_override is not None:
        raise CifConversionError(
            f"{cif_path}: the forced setting {space_override!r} does not "
            f"reproduce the CIF's {len(full_symbols)} atoms after "
            f"apply_space_group")
    raise CifConversionError(
        f"{cif_path}: could not resolve the space group -- tried "
        f"{candidates} for IT number {it_number}, none reproduced the "
        f"CIF's {len(full_symbols)} atoms.  Pass --space <token> to "
        f"force a spaceDB setting.")


def convert(cif_path, skl_path, space_override=None, title=None):
    """Convert ``cif_path`` to an imago.skl at ``skl_path`` and return
    the resolved space-group token."""

    lines, token = cif_to_skeleton(
        cif_path, space_override=space_override, title=title)
    with open(skl_path, "w", newline="\n") as handle:
        handle.writelines(lines)
    return token


def _default_skl_path(cif_path):
    """Return ``cif_path`` with its extension replaced by ``.skl``."""

    base, _ = os.path.splitext(cif_path)
    return base + ".skl"


def _build_parser():
    """Build the command-line argument parser."""

    description_text = (
        "Convert a CIF crystal structure into an imago.skl "
        "skeleton, preserving the space group.\n\n"
        "The skeleton stores the asymmetric unit plus a space-group "
        "token (not a flattened P1 cell), and the spaceDB setting is "
        "resolved by expanding each candidate and matching the CIF -- "
        "so an unresolvable space group is a hard error you can "
        "override with --space.")
    epilog_text = (
        "Examples:\n"
        "  cif2skl.py au_fcc.cif\n"
        "  cif2skl.py si.cif si_diamond.skl\n"
        "  cif2skl.py si.cif --space 227_b\n")

    parser = argparse.ArgumentParser(
        prog="cif2skl.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=description_text, epilog=epilog_text)
    parser.add_argument(
        "cif", metavar="INPUT.cif",
        help="Path to the input CIF file.  Required.")
    parser.add_argument(
        "skl", metavar="OUTPUT.skl", nargs="?", default=None,
        help="Path for the output skeleton.  Default: the input "
             "path with its extension changed to '.skl'.")
    parser.add_argument(
        "-s", "--space", default=None, metavar="TOKEN",
        help="Force a specific spaceDB setting token (e.g. '227_b') "
             "instead of resolving it automatically; the forced "
             "setting is still verified against the CIF's expansion.  "
             "Default: automatic resolution.")
    parser.add_argument(
        "-t", "--title", default=None, metavar="STR",
        help="Title recorded in the skeleton's title block.  "
             "Default: a note naming the source CIF.")
    return parser


def main(argv=None):
    """Entry point for the ``cif2skl.py`` command line."""

    args = _build_parser().parse_args(argv)
    skl_path = args.skl or _default_skl_path(args.cif)
    try:
        token = convert(args.cif, skl_path,
                        space_override=args.space, title=args.title)
    except CifConversionError as error:
        print(f"cif2skl: {error}", file=sys.stderr)
        return 1
    print(f"Wrote {skl_path}: space group {token}.")
    return 0


def record_command():
    """Append the issued command line to a file named "command" in
    the current directory, so the exact invocation can be recovered
    later.  This is a standing project convention: each run appends
    a dated block, so the file builds up a history of how the script
    was called."""

    with open("command", "a") as cmd:
        now = datetime.now()
        stamp = now.strftime("%b. %d, %Y: %H:%M:%S")
        cmd.write(f"Date: {stamp}\n")
        cmd.write("Cmnd:")
        for argument in sys.argv:
            cmd.write(f" {argument}")
        cmd.write("\n\n")


if __name__ == "__main__":
    record_command()
    sys.exit(main())
