#!/usr/bin/env python3
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

"""Expand a structure sketch into a complete curation manifest.

This is the authoring step between ``cod_fish`` and
``build_initial_potentials`` (DESIGN 5.7).  ``cod_fish`` discovers and
pins structures from the Crystallography Open Database; with
``--sketch-only`` it prints ``[[reference_solid]]`` stubs -- each
carrying only the structure's identity (``reference_id``,
``system_type``) and source (``cod_id`` + ``cod_revision``, or
``structure_path``), plus optional discovery hints.  This tool reads
that sketch (from a file, or from standard input as a pipe) and fills
in the rest of the manifest the producer needs: the database-wide
fingerprint recipe (the ``[characterization]`` block), the shared run
settings (basis, functional, k-point integration, SCF threshold) as a
single top-level ``[defaults]`` block that every solid inherits, and
any per-structure customizations (overrides on the environments the
harvest auto-discovers).

For a straightforward set of crystals ``cod_fish`` already emits a
complete, runnable manifest on its own; this tool exists for the cases
that need more -- local (non-COD) structures, a non-default recipe or
run settings, or interactively authored per-structure customizations.

``cod_fish`` and this tool and the producer all share one schema
definition, and the same shared defaults, by importing
``curation_manifest``, so nothing drifts between them.

Two modes:

* Interactive (``-i``) walks each sketched structure and prompts for
  the shared settings once and then the optional per-structure
  customizations, offering a sensible default at every step.  The
  database-wide fingerprint recipe is set automatically, so it writes
  a complete, immediately loadable manifest.  Because the prompts use
  standard input and output, interactive mode needs a sketch file
  argument and an output file (``-o``).
* Mechanical (the default) sets the shared settings once in a
  ``[defaults]`` block and the database-wide recipe once, carries each
  sketched structure through as a sparse solid, and emits the result
  with a commented copy-and-edit template for the optional
  per-structure customizations, which the curator fills in by hand.
  It reads the sketch from a file or, when none is named, from
  standard input; without ``-o`` it prints to standard output, so the
  result can be reviewed, redirected, or piped straight from
  ``cod_fish``.

The k-point density is deliberately left out of the emitted
``[defaults]`` ``kpoint_spec``: the producer predicts a starting
density and verifies it by a convergence sweep, so pinning a density
here would only override that.
"""

import argparse
import contextlib
import os
import re
import sys
import tempfile
import tomllib
from datetime import datetime

from curation_manifest import (
    CurationManifest,
    ReferenceSolid,
    ReferenceEntry,
    load_structure_sources,
    format_manifest,
    write_manifest,
    default_characterization,
)
# system_type is the one authoring vocabulary the manifest schema
#   itself validates (via curation_manifest's loader), so it is sourced
#   from its canonical owner rather than duplicated here; the canonical
#   reduce + bispectrum recipe and the [defaults] run settings likewise
#   come from curation_manifest, so nothing drifts between the two
#   authoring tools (DESIGN 5.7).
from guidance_db import VALID_SYSTEM_TYPES


# Human-readable vocabularies, used in help text and interactive
#   menus.  The producer translates these names into each tool's coded
#   form; an unrecognised name is caught later by the producer, so the
#   menus are guidance rather than a hard gate here.  (system_type is
#   not here -- it is imported from its schema owner above.)
VALID_FUNCTIONALS = (
    "wigner", "ceperley-alder", "hedin-lundqvist", "pbe")
VALID_INTEGRATIONS = ("gaussian", "linear-tetrahedral")
VALID_BASES = ("mb", "fb", "eb")


# A commented copy-and-edit template appended to mechanical-mode
#   output.  The mechanical output already carries the required
#   [characterization] recipe and every structure's run settings, so
#   it is loadable as-is: the harvest auto-discovers one
#   representative per distinct environment with no entries at all.
#   Entries are OPTIONAL customizations -- add one only to override
#   what the harvest would otherwise produce for an environment.
#   load_manifest_v2 ignores comments, so these stay inert until
#   uncommented.
ENTRY_TEMPLATE = """\
# -----------------------------------------------------------------
# Per-structure customizations -- OPTIONAL.
#
# The harvest already discovers and stores one representative per
#   distinct environment in every structure above, computing the
#   [characterization] recipe for each.  Add an
#   [[reference_solid.entry]] block (under its [[reference_solid]],
#   before the next one) ONLY to override what the harvest would
#   otherwise produce -- every field is optional:
#
# [[reference_solid.entry]]
# atom_site   = 1               # representative atom of the target
# element     = "Si"            # cross-checked at harvest
# default     = true            # this environment is the element
#                               #   default (at most one per element;
#                               #   omitted -> the isolated baseline)
# description = "describe this site"   # else auto-composed
# label       = "explicit-label"       # else derived at harvest
#
# The preferred fingerprint recipe is the database-wide
#   [characterization] block above, NOT a per-entry choice.  A
#   per-entry fingerprint is a RARE extra, non-preferred descriptor
#   harvested for this one environment (no preferred flag):
#
# [[reference_solid.entry.fingerprint]]
# method   = "bispectrum"
# sub_spec = { twoj1 = 6, twoj2 = 4 }
# -----------------------------------------------------------------
"""


def load_sketch_hints(path: str) -> dict[str, dict]:
    """Read the optional discovery hints cod_fish writes into a sketch
    -- per reference_id, the composition (``elements``) and a
    ``source_description`` -- so the interactive flow can offer them
    as the element and description defaults.  These are not manifest
    schema fields (the producer ignores them); a sketch without them
    just yields empty hints, so a hand-written sketch still works."""

    with open(path, "rb") as handle:
        raw = tomllib.load(handle)
    hints: dict[str, dict] = {}
    for ref in raw.get("reference_solid", []):
        reference_id = ref.get("reference_id")
        if reference_id is None:
            continue
        hints[reference_id] = {
            "elements": ref.get("elements") or [],
            "description": ref.get("source_description", "")}
    return hints


def shared_defaults(*, basis: str, functional: str,
                    kpoint_integration: str,
                    scf_threshold: float) -> dict:
    """Build the top-level ``[defaults]`` run-settings dict expand
    emits once (DESIGN 5.7), so each solid stays sparse and inherits
    the shared settings rather than repeating all five.

    ``kpoint_spec`` is left empty (no fixed density) so the producer
    predicts a starting density and verifies it by a convergence
    sweep; pinning one here would only override that.  ``system_type``
    is deliberately not among these five -- it is an intrinsic
    structure property, named per solid, not a shared run setting."""

    return {
        "basis": basis,
        "functional": functional,
        "kpoint_integration": kpoint_integration,
        "kpoint_spec": {},
        "scf_threshold": scf_threshold,
    }


def sparse_solid(source: ReferenceSolid, system_type: str, *,
                 entries: list[ReferenceEntry] | None = None
                 ) -> ReferenceSolid:
    """Return a reference solid that copies ``source``'s identity and
    structure source and resolves its ``system_type``, leaving every
    run setting unset (``None``) so it inherits the top-level
    ``[defaults]`` block (DESIGN 5.7).

    ``source`` is one structure stub from the sketch (as
    :func:`curation_manifest.load_structure_sources` yields it), so it
    carries a ``reference_id`` and exactly one structure source.
    ``system_type`` is the already-resolved value (from the sketch
    when it named one, else the shared default).  The
    ``source_description`` cod_fish read from the CIF rides along, so
    the harvest can compose each environment's auto-description from
    it.  ``entries`` are the per-structure customizations -- empty in
    mechanical mode, filled by the interactive walk."""

    return ReferenceSolid(
        reference_id=source.reference_id,
        system_type=system_type,
        cod_id=source.cod_id,
        cod_revision=source.cod_revision,
        structure_path=source.structure_path,
        source_description=source.source_description,
        entries=list(entries) if entries else [])


def build_mechanical(sources: list[ReferenceSolid], *, basis: str,
                     functional: str, kpoint_integration: str,
                     scf_threshold: float,
                     system_type_default: str) -> CurationManifest:
    """Set the shared defaults once in a ``[defaults]`` block and the
    database-wide ``[characterization]`` recipe once, then carry every
    sketched structure through as a sparse solid (identity, structure
    source, and system_type only), leaving the per-structure
    customizations empty for the curator.  The result is a
    structurally valid (loadable) manifest -- the required recipe
    block is present (manifest rule 2) and every solid's run settings
    resolve from ``[defaults]`` -- that simply harvests every
    environment with the default descriptions until customizations are
    added."""

    solids = [
        sparse_solid(source, source.system_type or system_type_default)
        for source in sources]
    return CurationManifest(
        schema_version=2, manifest_path="(generated)",
        characterization=default_characterization(),
        defaults=shared_defaults(
            basis=basis, functional=functional,
            kpoint_integration=kpoint_integration,
            scf_threshold=scf_threshold),
        reference_solids=solids)


def _as_bool(answer: str, default: bool) -> bool:
    """Interpret a yes/no answer; an empty answer keeps ``default``."""

    text = answer.strip().lower()
    if not text:
        return default
    return text in ("y", "yes", "true", "1")


def _ask_yes_no(ask, prompt: str, default: bool) -> bool:
    """Ask a yes/no question through the injected ``ask`` callable."""

    return _as_bool(ask(f"{prompt} (y/n)", "y" if default else "n"),
                    default)


# The longest answer any prompt here has a use for.  A curator types
#   element symbols, short descriptions and labels; anything longer is
#   a paste accident, and refusing it keeps a runaway string out of the
#   manifest rather than discovering it later during a read.
MAX_ANSWER_LENGTH = 200

# Characters that carry special meaning to a command shell.  No answer
#   this script accepts has a use for one: the fields are element
#   symbols, site numbers, short labels and short descriptions.  They
#   are refused at the prompt rather than escaped further downstream,
#   so that a value which could be misread by anything that later
#   handles a manifest never enters the file in the first place.
SHELL_METACHARACTERS = frozenset(";|&$`<>\\(){}[]")

# The shapes the two structured free-text answers may take.  Both
#   permit an empty answer, because both prompts legitimately default
#   to blank -- a label left blank is derived at harvest time, and the
#   element default is blank when the sketch carried no composition
#   hint.  Refusing blank here would leave a curator with no way past
#   the prompt.
ELEMENT_PATTERN = re.compile(r"[A-Za-z]{0,3}")
LABEL_PATTERN = re.compile(r"[A-Za-z0-9._:-]*")


def clean_answer(raw: str) -> str:
    """Normalize one typed answer before it is used for anything.

    Two things are removed.  Control characters -- anything below
    ``0x20`` plus ``DEL`` -- are stripped, because a terminal answer
    that carries an embedded newline, escape sequence or NUL is either
    a paste accident or an attempt to smuggle structure into a value,
    and no legitimate answer here contains one.  Surrounding whitespace
    goes with them.

    An over-long answer is refused outright rather than silently
    truncated, since a truncated element symbol or label would be
    wrong in a way nobody would notice until much later.

    The manifest writer already escapes what it emits (the TOML basic
    string quoting in ``curation_manifest._quote``), so this is a
    second, earlier line of defence rather than the only one: it keeps
    a value that cannot be meant from entering the program at all,
    instead of relying on every later consumer to render it safely.
    """

    cleaned = "".join(character for character in raw
                      if ord(character) >= 0x20 and character != "\x7f")
    cleaned = cleaned.strip()
    if len(cleaned) > MAX_ANSWER_LENGTH:
        raise ValueError(
            f"answer is {len(cleaned)} characters; the limit is "
            f"{MAX_ANSWER_LENGTH}")
    found = sorted(set(cleaned) & SHELL_METACHARACTERS)
    if found:
        raise ValueError(
            f"answer contains {' '.join(found)}, which cannot appear "
            f"in a manifest field; please rephrase without it")
    return cleaned


def _as_element(answer: str) -> str:
    """Return ``answer`` if it looks like an element symbol, else raise.

    One to three letters, or blank when the sketch offered no
    composition hint and the curator has nothing to fill in yet.
    """

    if not ELEMENT_PATTERN.fullmatch(answer):
        raise ValueError(
            f"{answer!r} is not an element symbol (one to three "
            f"letters)")
    return answer


def _as_label(answer: str) -> str:
    """Return ``answer`` if it is a usable label, else raise.

    Labels name an entry for later reference, so they are held to
    letters, digits and the few separators those names use.  Blank is
    allowed and means the label is derived at harvest time.
    """

    if not LABEL_PATTERN.fullmatch(answer):
        raise ValueError(
            f"{answer!r} is not a valid label (letters, digits and "
            f". _ : - only)")
    return answer


def _as_system_type(answer: str) -> str:
    """Return ``answer`` if it names a system type, else raise."""

    if answer not in VALID_SYSTEM_TYPES:
        raise ValueError(
            f"{answer!r} is not one of "
            f"{'/'.join(VALID_SYSTEM_TYPES)}")
    return answer


def _as_atom_site(answer: str) -> int:
    """Return ``answer`` as a site index, else raise.

    Sites are counted from one, so zero and negatives are refused here
    rather than becoming an out-of-range index much later.
    """

    try:
        site = int(answer)
    except ValueError:
        raise ValueError(f"{answer!r} is not a whole number") from None
    if site < 1:
        raise ValueError(f"atom_site must be 1 or greater, not {site}")
    return site


def _ask_until_valid(ask, prompt: str, default: str, convert):
    """Ask until ``convert`` accepts the answer, then return its value.

    ``convert`` turns the typed text into the value the caller wants
    and raises ``ValueError`` when it cannot.  A rejected answer is
    reported and the question is asked again, so a slip of the finger
    costs one retry rather than the whole curation session -- the
    manifest is only written after every structure has been walked
    through, so an exception here would discard all of it.
    """

    while True:
        try:
            return convert(ask(prompt, default))
        except ValueError as error:
            print(f"  -- {error}; please try again.", file=sys.stderr)


def build_interactive(sources: list[ReferenceSolid], ask, *,
                      basis: str, functional: str,
                      kpoint_integration: str, scf_threshold: float,
                      system_type_default: str,
                      hints: dict[str, dict] | None = None
                      ) -> CurationManifest:
    """Walk the curator through completing each sketched structure.

    ``ask(prompt, default)`` returns the curator's answer, or
    ``default`` when they accept it (an empty reply); injecting it
    rather than calling ``input`` directly keeps this flow testable.

    The shared method defaults are confirmed once up front (each
    pre-filled with the value passed in), and the database-wide
    ``[characterization]`` recipe is stamped automatically (the
    canonical reduce + bispectrum sub-specs, so the file satisfies
    manifest rule 2).  Then, for each structure, the curator sets its
    system type and adds optional customizations (5.2.2): each names
    an atom site and overrides what the harvest would otherwise
    produce -- the element, the default flag, the description, or the
    label.  The fingerprint recipe is NOT a per-entry choice; it
    lives in ``[characterization]``, so no per-entry fingerprints are
    authored here.  The rare non-preferred per-entry override is
    added to the emitted file by hand (the producer harvests it
    alongside the preferred recipe)."""

    print("Shared run settings (press Enter to accept the default):",
          file=sys.stderr)
    basis = ask(f"basis ({'/'.join(VALID_BASES)})", basis)
    functional = ask(
        f"functional ({'/'.join(VALID_FUNCTIONALS)})", functional)
    kpoint_integration = ask(
        f"kpoint_integration "
        f"({'/'.join(VALID_INTEGRATIONS)}, or gaussian-<sigma>)",
        kpoint_integration)
    scf_threshold = float(ask("scf_threshold", str(scf_threshold)))

    # cod_fish's per-structure discovery hints (composition +
    #   description), keyed by reference_id; empty for a hand-written
    #   sketch, in which case the element and description just prompt
    #   with blank defaults.
    hints = hints or {}

    # Elements that already hold their one default customization,
    #   tracked across the whole manifest so rule 7 (at most one
    #   default per element) holds by default.
    default_elements: set[str] = set()
    solids: list[ReferenceSolid] = []

    for source in sources:
        # Announce the structure (a printed header, NOT a prompt) so
        #   the questions that follow are clearly about this solid.
        print(f"\n--- structure {source.reference_id} ---",
              file=sys.stderr)
        hint = hints.get(source.reference_id, {})
        hint_elements = hint.get("elements") or []
        hint_description = hint.get("description", "")
        # Check the answer against the vocabulary the prompt shows.
        #   Without this the value is only validated much later, when
        #   the finished manifest is read back (guidance_db), so a typo
        #   surfaced after the whole session rather than at the prompt
        #   that caused it.
        system_type = _ask_until_valid(
            ask, f"system_type ({'/'.join(VALID_SYSTEM_TYPES)})",
            source.system_type or system_type_default,
            _as_system_type)

        entries: list[ReferenceEntry] = []
        # The first entry defaults to yes; once one is added, "add
        #   another" defaults to no, so pressing Enter ends this
        #   structure rather than looping forever.
        add_entry = _ask_yes_no(
            ask, f"add an entry for {source.reference_id}?", True)
        while add_entry:
            # Default the element to the structure's composition (the
            #   next not-yet-entered element, so a single-element
            #   solid auto-fills it) and the description to the
            #   CIF-derived hint -- both from cod_fish (DESIGN 5.7).
            element_default = (
                hint_elements[len(entries)]
                if len(entries) < len(hint_elements)
                else (hint_elements[0] if hint_elements else ""))
            element = _ask_until_valid(
                ask, "  element", element_default, _as_element)
            # A non-numeric atom_site used to raise straight out of the
            #   loop, losing every answer given so far; re-ask instead.
            atom_site = _ask_until_valid(
                ask, "  atom_site", "1", _as_atom_site)
            # Default to yes only until this element has its one
            #   default entry, so accepting defaults yields exactly
            #   one default per element (rule 7).
            is_default = _ask_yes_no(
                ask, "  default entry for this element?",
                element not in default_elements)
            if is_default:
                default_elements.add(element)
            description = ask("  description", hint_description)
            label = _ask_until_valid(
                ask, "  label (blank to derive at harvest)", "",
                _as_label)

            # No per-entry fingerprints are authored here: the
            #   preferred recipe is the database-wide
            #   [characterization] block (set once below), and a
            #   per-entry override is a rare hand-edit.
            entries.append(ReferenceEntry(
                element=element, atom_site=atom_site,
                default=is_default, description=description,
                label=label or None, fingerprints=[]))
            add_entry = _ask_yes_no(
                ask, f"add another entry for {source.reference_id}?",
                False)

        # The solid stays sparse: its run settings inherit the shared
        #   [defaults] block below (DESIGN 5.7).  system_type is the
        #   per-structure answer, and the sketch's source_description
        #   rides along via sparse_solid.
        solids.append(sparse_solid(source, system_type,
                                   entries=entries))

    return CurationManifest(
        schema_version=2, manifest_path="(generated)",
        characterization=default_characterization(),
        defaults=shared_defaults(
            basis=basis, functional=functional,
            kpoint_integration=kpoint_integration,
            scf_threshold=scf_threshold),
        reference_solids=solids)


def _input_ask(prompt: str, default: str) -> str:
    """The real prompt: show the default in brackets and return the
    curator's reply, or the default when they press Enter.

    The reply is passed through :func:`clean_answer` before it is
    returned, so nothing downstream ever sees a raw terminal line.
    Every answer this script accepts ends up as a field of the
    manifest it writes, so the value is normalized at the one point it
    enters the program rather than at each of the places it is used.

    The prompt repeats when :func:`clean_answer` refuses a reply, so
    a rejected answer costs one retry rather than raising out of the
    curation loop and discarding the session's work.  Re-asking here
    rather than in each caller means every prompt is covered, including
    the free-text ones that need no further checking of their own.
    """

    while True:
        try:
            raw = clean_answer(input(f"{prompt} [{default}]: "))
        except ValueError as error:
            print(f"  -- {error}.", file=sys.stderr)
            continue
        return raw if raw else default


@contextlib.contextmanager
def sketch_path(named: str | None):
    """Yield a filesystem path to the sketch to read.

    With a named file, yield it unchanged.  With none, capture
    standard input to a temporary file and yield that path -- so
    ``cod_fish.py pin ... --sketch-only | expand_manifest.py`` works
    as a pipe.  stdin is captured to a real file rather than streamed
    because the sketch readers (:func:`load_structure_sources`,
    :func:`load_sketch_hints`) validate and parse a path, and
    interactive mode reads the sketch twice (structures, then hints).
    The temporary file is removed on exit."""

    if named:
        yield named
        return
    data = sys.stdin.buffer.read()
    handle = tempfile.NamedTemporaryFile(
        mode="wb", suffix=".toml", delete=False)
    try:
        handle.write(data)
        handle.close()
        yield handle.name
    finally:
        os.unlink(handle.name)


def _build_parser() -> argparse.ArgumentParser:
    """Build the flat argument parser (one ``--help`` lists every
    option with its default, the convention the other Imago scripts
    follow)."""

    parser = argparse.ArgumentParser(
        prog="expand_manifest.py",
        description=(
            "Expand a structure sketch (from cod_fish) into a "
            "complete curation manifest by stamping the shared "
            "method defaults and the per-structure harvest "
            "curation."))
    parser.add_argument(
        "sketch", metavar="SKETCH", nargs="?", default=None,
        help="the sketch manifest: a schema-version-2 file carrying "
             "only the structure stubs (reference_id, system_type, "
             "and one structure source per solid).  Default: read the "
             "sketch from standard input, so `cod_fish.py pin ... "
             "--sketch-only | expand_manifest.py` works as a pipe "
             "(interactive mode still needs a file argument).")
    parser.add_argument(
        "-o", "--output", default=None, metavar="PATH",
        help="write the manifest to PATH.  Default: print to "
             "standard output (required with --interactive, whose "
             "prompts use standard output).")
    parser.add_argument(
        "-i", "--interactive", action="store_true", default=False,
        help="walk through the per-structure curation with prompts "
             "instead of emitting a fill-in template.  Default: off "
             "(mechanical template).")
    parser.add_argument(
        "--basis", default="fb",
        help=f"shared basis for every structure, one of "
             f"{'/'.join(VALID_BASES)}.  Default: fb (full basis).")
    parser.add_argument(
        "--functional", default="wigner",
        help=f"shared exchange-correlation functional, one of "
             f"{'/'.join(VALID_FUNCTIONALS)}.  Default: wigner.")
    parser.add_argument(
        "--kpoint-integration", dest="kpoint_integration",
        default="gaussian",
        help=f"shared k-point integration, one of "
             f"{'/'.join(VALID_INTEGRATIONS)} (or gaussian-<sigma> "
             f"to name a thermal smearing in eV).  Default: "
             f"gaussian.")
    parser.add_argument(
        "--scf-threshold", dest="scf_threshold", type=float,
        default=1.0e-6,
        help="shared SCF convergence threshold.  Default: 1e-6.")
    parser.add_argument(
        "--system-type", dest="system_type", default="crystalline",
        help=f"system type stamped on any structure whose sketch did "
             f"not name one, one of {'/'.join(VALID_SYSTEM_TYPES)}.  "
             f"Default: crystalline.")
    return parser


def main(argv=None) -> int:
    """Entry point for the ``expand_manifest.py`` command line."""

    args = _build_parser().parse_args(argv)

    # Interactive prompts read standard input, so the sketch cannot
    #   also be piped in: interactive mode needs a named sketch file.
    if args.interactive and not args.sketch:
        _build_parser().error(
            "--interactive needs a SKETCH file argument (its prompts "
            "read standard input, so the sketch cannot also be piped "
            "in)")

    shared = dict(
        basis=args.basis, functional=args.functional,
        kpoint_integration=args.kpoint_integration,
        scf_threshold=args.scf_threshold,
        system_type_default=args.system_type)

    # Read the sketch from the named file, or from standard input when
    #   none was given (captured to a temp path the readers can parse).
    with sketch_path(args.sketch) as path:
        sources = load_structure_sources(path)

        if args.interactive:
            if not args.output:
                _build_parser().error(
                    "--interactive requires -o/--output (its prompts "
                    "use standard output)")
            hints = load_sketch_hints(path)
            manifest = build_interactive(
                sources, _input_ask, hints=hints, **shared)
            write_manifest(manifest, args.output)
            print(f"expand_manifest: wrote {len(sources)} reference "
                  f"solids to {args.output}", file=sys.stderr)
            return 0

        manifest = build_mechanical(sources, **shared)
        text = format_manifest(manifest) + "\n" + ENTRY_TEMPLATE
        if args.output:
            with open(args.output, "w") as handle:
                handle.write(text)
        else:
            sys.stdout.write(text)
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
