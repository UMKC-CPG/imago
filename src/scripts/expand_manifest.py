#!/usr/bin/env python3
"""Expand a structure sketch into a complete curation manifest.

This is the authoring step between ``cod_fish`` and
``build_initial_potentials`` (DESIGN 5.7).  ``cod_fish`` discovers and
pins structures from the Crystallography Open Database and prints
sketch ``[[reference_solid]]`` stubs -- each carrying only the
structure's identity (``reference_id``, ``system_type``) and source
(``cod_id`` + ``cod_revision``, or ``structure_path``).  The curator
collects those stubs into a sketch file.  This tool reads that sketch
and fills in the rest of the manifest the producer needs: the shared
method defaults (basis, functional, k-point integration, SCF
threshold) and the per-structure harvest curation (which atom sites
become database entries, and the fingerprints to compute for each).

``cod_fish`` stays a pure discovery tool and never writes a manifest;
this tool owns the sketch-to-manifest expansion, and the producer only
reads the finished manifest.  Both this tool and the producer share
one schema definition by importing ``curation_manifest``.

Two modes:

* Interactive (``-i``) walks each sketched structure and prompts for
  the shared defaults once and then the per-structure entries and
  fingerprints, offering a sensible default at every step.  It writes
  a complete, immediately loadable manifest.  Because the prompts use
  standard output, interactive mode requires an output file (``-o``).
* Mechanical (the default) stamps the shared defaults onto every
  sketched structure and emits the result with a commented
  copy-and-edit template for the per-structure entries, which the
  curator fills in by hand.  Without ``-o`` it prints to standard
  output, so the result can be reviewed or redirected.

The k-point density is deliberately left out of every emitted
``kpoint_spec``: the producer predicts a starting density and verifies
it by a convergence sweep, so pinning a density here would only
override that.
"""

import argparse
import sys
import tomllib
from datetime import datetime

from curation_manifest import (
    CurationManifest,
    ReferenceSolid,
    ReferenceEntry,
    ManifestFingerprint,
    load_structure_sources,
    format_manifest,
    write_manifest,
)


# The settled canonical fingerprint sub-specs (the Si seed defaults).
#   Rule 11 fixes one sub_spec per method across the whole database,
#   so these are offered as the interactive defaults and shown in the
#   mechanical-mode template.  The reduce descriptor is a pure-Python
#   shell code; the bispectrum descriptor needs a per-structure loen
#   run but discriminates disordered environments far better.
DEFAULT_REDUCE_SUB_SPEC = {
    "level": 2, "thick": 0.5, "cutoff": 5.0, "tolerance": 0.05}
DEFAULT_BISPECTRUM_SUB_SPEC = {
    "twoj1": 8, "twoj2": 8, "cutoff": 9.0}

# Human-readable vocabularies, used in help text and interactive
#   menus.  The producer translates these names into each tool's coded
#   form; an unrecognised name is caught later by the producer, so the
#   menus are guidance rather than a hard gate here.
VALID_FUNCTIONALS = (
    "wigner", "ceperley-alder", "hedin-lundqvist", "pbe")
VALID_INTEGRATIONS = ("gaussian", "linear-tetrahedral")
VALID_BASES = ("mb", "fb", "eb")
VALID_SYSTEM_TYPES = (
    "crystalline", "amorphous", "nanostructure", "molecular")


# A commented copy-and-edit template appended to mechanical-mode
#   output.  The curator pastes one [[reference_solid.entry]] block
#   (and its fingerprints) under each [[reference_solid]] above and
#   edits the element, site, and description.  load_manifest_v2
#   ignores comments, so the mechanical output is already loadable
#   (with no entries) until these are filled in.
ENTRY_TEMPLATE = """\
# -----------------------------------------------------------------
# Per-structure curation -- fill this in.
#
# Under each [[reference_solid]] above, add one or more
#   [[reference_solid.entry]] blocks (before the next
#   [[reference_solid]]).  Each names an atom site whose converged
#   potential becomes a database entry.  Exactly one entry per
#   element across the WHOLE manifest must carry default = true.
#   Omit label to derive it at harvest from the run's site identity.
#   Copy, uncomment, and edit:
#
# [[reference_solid.entry]]
# element = "Si"
# atom_site = 1
# default = true
# description = "describe this site"
# # label = "explicit-label"   # optional override; omit to derive
#
# [[reference_solid.entry.fingerprint]]
# method = "reduce"
# sub_spec = { level = 2, thick = 0.5, cutoff = 5.0, tolerance = 0.05 }
# preferred = true
#
# [[reference_solid.entry.fingerprint]]
# method = "bispectrum"
# sub_spec = { twoj1 = 8, twoj2 = 8, cutoff = 9.0 }
# preferred = true
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


def stamp_shared_defaults(source: ReferenceSolid, *, basis: str,
                          functional: str, kpoint_integration: str,
                          scf_threshold: float,
                          system_type_default: str
                          ) -> ReferenceSolid:
    """Return a reference solid that copies ``source``'s identity and
    structure source and adds the shared method defaults.

    ``source`` is one structure stub from the sketch (as
    :func:`curation_manifest.load_structure_sources` yields it), so it
    carries a ``reference_id`` and exactly one structure source but
    placeholder run fields.  ``system_type`` is taken from the sketch
    when it named one, else the supplied default.  The ``kpoint_spec``
    is left empty (no density) so the producer predicts and verifies
    it.  Entries are carried through unchanged -- the caller fills
    them (interactively or by hand)."""

    return ReferenceSolid(
        reference_id=source.reference_id,
        system_type=source.system_type or system_type_default,
        basis=basis,
        functional=functional,
        kpoint_integration=kpoint_integration,
        kpoint_spec={},
        scf_threshold=scf_threshold,
        cod_id=source.cod_id,
        cod_revision=source.cod_revision,
        structure_path=source.structure_path,
        entries=list(source.entries))


def build_mechanical(sources: list[ReferenceSolid], *, basis: str,
                     functional: str, kpoint_integration: str,
                     scf_threshold: float,
                     system_type_default: str) -> CurationManifest:
    """Stamp the shared defaults onto every sketched structure,
    leaving the per-structure entries empty for the curator.  The
    result is a structurally valid (loadable) manifest; it simply
    harvests nothing until entries are added."""

    solids = [
        stamp_shared_defaults(
            source, basis=basis, functional=functional,
            kpoint_integration=kpoint_integration,
            scf_threshold=scf_threshold,
            system_type_default=system_type_default)
        for source in sources]
    return CurationManifest(
        schema_version=2, manifest_path="(generated)",
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
    pre-filled with the value passed in).  Then, for each structure,
    the curator sets its system type and adds entries; each entry may
    carry the canonical reduce and bispectrum fingerprints.  The
    ``preferred`` flag is managed automatically: the first entry to
    declare a given ``(element, method)`` is marked preferred, so the
    manifest satisfies the one-preferred-per-family rule (rule 10)
    without the curator tracking it."""

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

    # (element, method) pairs already marked preferred, and elements
    #   that already hold their one default entry -- both tracked
    #   across the whole manifest so rules 10 and 7 hold by default.
    preferred_seen: set[tuple[str, str]] = set()
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
        system_type = ask(
            f"system_type ({'/'.join(VALID_SYSTEM_TYPES)})",
            source.system_type or system_type_default)

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
            element = ask("  element", element_default)
            atom_site = int(ask("  atom_site", "1"))
            # Default to yes only until this element has its one
            #   default entry, so accepting defaults yields exactly
            #   one default per element (rule 7).
            is_default = _ask_yes_no(
                ask, "  default entry for this element?",
                element not in default_elements)
            if is_default:
                default_elements.add(element)
            description = ask("  description", hint_description)
            label = ask("  label (blank to derive at harvest)", "")

            fingerprints: list[ManifestFingerprint] = []
            if _ask_yes_no(
                    ask, "  attach canonical reduce + bispectrum "
                    "fingerprints?", True):
                for method, sub_spec in (
                        ("reduce", DEFAULT_REDUCE_SUB_SPEC),
                        ("bispectrum", DEFAULT_BISPECTRUM_SUB_SPEC)):
                    key = (element, method)
                    preferred = key not in preferred_seen
                    if preferred:
                        preferred_seen.add(key)
                    fingerprints.append(ManifestFingerprint(
                        method=method, sub_spec=dict(sub_spec),
                        preferred=preferred))

            entries.append(ReferenceEntry(
                element=element, atom_site=atom_site,
                default=is_default, description=description,
                label=label or None, fingerprints=fingerprints))
            add_entry = _ask_yes_no(
                ask, f"add another entry for {source.reference_id}?",
                False)

        solids.append(stamp_shared_defaults(
            ReferenceSolid(
                reference_id=source.reference_id,
                system_type=system_type, basis="", functional="",
                kpoint_integration="", kpoint_spec={},
                scf_threshold=0.0, cod_id=source.cod_id,
                cod_revision=source.cod_revision,
                structure_path=source.structure_path,
                entries=entries),
            basis=basis, functional=functional,
            kpoint_integration=kpoint_integration,
            scf_threshold=scf_threshold,
            system_type_default=system_type))

    return CurationManifest(
        schema_version=2, manifest_path="(generated)",
        reference_solids=solids)


def _input_ask(prompt: str, default: str) -> str:
    """The real prompt: show the default in brackets and return the
    curator's reply, or the default when they press Enter."""

    raw = input(f"{prompt} [{default}]: ").strip()
    return raw if raw else default


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
        "sketch", metavar="SKETCH",
        help="the sketch manifest: a schema-version-2 file carrying "
             "only the structure stubs (reference_id, system_type, "
             "and one structure source per solid).  Required.")
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
    sources = load_structure_sources(args.sketch)

    shared = dict(
        basis=args.basis, functional=args.functional,
        kpoint_integration=args.kpoint_integration,
        scf_threshold=args.scf_threshold,
        system_type_default=args.system_type)

    if args.interactive:
        if not args.output:
            _build_parser().error(
                "--interactive requires -o/--output (its prompts "
                "use standard output)")
        hints = load_sketch_hints(args.sketch)
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
