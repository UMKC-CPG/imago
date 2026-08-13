#!/usr/bin/env python3
"""Audit a release-build log against the recorded warning manifest.

The debug campaign (dev/DEBUG.md) drove the audit-build warning
count to zero, but gfortran's -Wmaybe-uninitialized family fires
only at -O2 and above, so the release build permanently shows a
set of adjudicated-benign warnings the audit build never can.  A
standing list of known-benign warnings is indistinguishable from
a standing list of unread ones: a NEW warning scrolling past
inside it would be invisible.  This tool restores the property
the zero-warning doctrine actually cares about, by diffing a
fresh build log against the committed manifest of expected
warnings so that any change -- a new warning appearing, or a
recorded one vanishing -- is loud and machine-detected.

Each warning is reduced to a group key of (source file, warning
class, variable name), the same granularity the campaign's
adjudication used.  Compiler-internal descriptor components such
as 'var.dim[1].stride' or 'var.offset' collapse to the base
variable 'var', because they all describe the same adjudicated
site.  Source line numbers are deliberately NOT part of the key:
they drift with every unrelated edit, and the adjudications in
dev/DEBUG.md are per-variable, not per-line.  The known limit of
this granularity: a genuinely new warning about an
already-listed variable in the same file is folded into its
existing group and will not be flagged.

The manifest is a committed TSV file, one row per expected group
per build variant, refreshed deliberately (with --write-manifest)
only after a change to the expected set has been reviewed and
recorded in dev/DEBUG.md.
"""

## The logs this tool audits are produced by from-clean builds of
## the two release variants, each into its OWN log (the two
## targets' compile lines interleave under a shared -j log, which
## makes variant attribution unreliable -- dev/DEBUG.md records
## the method rule):
##
##     cd build/release && make clean
##     make -j8 imagoG > uninit_real.log    2>&1
##     make -j8 imago  > uninit_complex.log 2>&1

import argparse
import re
import sys
from pathlib import Path

# The committed manifest lives beside this script so the pair
#   travels together through checkouts and reviews.
DEFAULT_MANIFEST = Path(__file__).parent / \
    "release_warning_manifest.tsv"

# A gfortran diagnostic location line: "path/file.F90:277:47:".
#   It may stand alone (with the Warning: line following a source
#   excerpt) or prefix the Warning: text on the same line.
LOCATION_PATTERN = re.compile(
    r"^(?P<path>\S+\.[fF]\d*):(?P<line>\d+):(?P<col>\d+):")

# The warning text itself, wherever it appears on a line.  The
#   trailing bracket names the warning class, e.g.
#   [-Wmaybe-uninitialized].
WARNING_PATTERN = re.compile(
    r"Warning:\s*(?P<message>.*?)\s*"
    r"(?:\[(?P<warn_class>-W[a-z0-9-]+)\])?\s*$")

# The variable or dummy-argument name is the first single-quoted
#   token of the message, when one exists.
QUOTED_NAME_PATTERN = re.compile(r"'(?P<name>[^']+)'")


def parse_build_log(log_path):
    """Extract warning group keys from one variant's build log.

    Returns (groups, unattributed, unparsed):
      groups       -- set of (file, warning class, base name)
      unattributed -- warning lines seen before any location line,
                      grouped under the file '<unknown>'
      unparsed     -- raw lines containing 'Warning:' that the
                      patterns failed to interpret at all
    Nothing is silently dropped: every 'Warning:' line lands in
    one of the three, so a format drift in a future gfortran
    shows up as noise in the report instead of a quiet miss.
    """
    groups = set()
    unattributed = []
    unparsed = []
    current_file = None

    with open(log_path, encoding="utf-8", errors="replace") as log:
        for raw_line in log:
            line = raw_line.rstrip("\n")

            # A location line updates the attribution context for
            #   any bare Warning: line that follows it.  When the
            #   warning shares the line, the same match provides
            #   the file directly.
            location = LOCATION_PATTERN.match(line)
            if location:
                current_file = Path(location.group("path")).name

            if "Warning:" not in line:
                continue

            warning = WARNING_PATTERN.search(line)
            if warning is None:
                unparsed.append(line)
                continue

            # The group's name: the quoted variable, with any
            #   descriptor component ('.dim[1].stride', '.offset')
            #   collapsed to the base variable it describes.  A
            #   message with no quoted token keeps its whole text
            #   as the name so it still forms a diffable group.
            quoted = QUOTED_NAME_PATTERN.search(
                warning.group("message"))
            if quoted:
                base_name = quoted.group("name") \
                    .split(".")[0].lower()
            else:
                base_name = warning.group("message").strip()

            warn_class = warning.group("warn_class") \
                or "<unclassified>"
            source_file = current_file if location or current_file \
                else "<unknown>"

            key = (source_file, warn_class, base_name)
            groups.add(key)
            if source_file == "<unknown>":
                unattributed.append(line)

    return groups, unattributed, unparsed


def read_manifest(manifest_path):
    """Load the manifest into {variant: set of group keys}."""
    expected = {}
    with open(manifest_path, encoding="utf-8") as manifest:
        for raw_line in manifest:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 4:
                sys.exit(f"Malformed manifest row (expected 4 "
                         f"tab-separated fields): {line!r}")
            variant, source_file, warn_class, base_name = fields
            expected.setdefault(variant, set()).add(
                (source_file, warn_class, base_name))
    return expected


def write_manifest(manifest_path, groups_by_variant):
    """Write the manifest from freshly parsed logs.

    Rows are sorted so the committed file diffs cleanly from one
    refresh to the next.
    """
    with open(manifest_path, "w", encoding="utf-8") as manifest:
        manifest.write(
            "# Expected release-build warnings, one row per "
            "(variant, file,\n"
            "#   warning class, variable) group.  Every row is "
            "an adjudicated\n"
            "#   entry recorded in dev/DEBUG.md's release-build "
            "uninitialized\n"
            "#   tranche or Phase 1 accepted-warning notes.  "
            "Refresh only via\n"
            "#   check_release_warnings.py --write-manifest, "
            "after the change\n"
            "#   in the expected set has been reviewed and "
            "recorded there.\n")
        for variant in sorted(groups_by_variant):
            for key in sorted(groups_by_variant[variant]):
                manifest.write(
                    "\t".join((variant,) + key) + "\n")


def report_variant(variant, found, expected):
    """Print the comparison for one variant; True when clean."""
    new_groups = sorted(found - expected)
    missing_groups = sorted(expected - found)

    if not new_groups and not missing_groups:
        print(f"{variant}: matches the recorded standing list "
              f"({len(expected)} warning groups).")
        return True

    print(f"{variant}: DIFFERS from the recorded standing list "
          f"({len(expected)} groups expected, "
          f"{len(found)} found).")
    for source_file, warn_class, base_name in new_groups:
        print(f"  NEW      {source_file}  {warn_class}  "
              f"{base_name}")
    for source_file, warn_class, base_name in missing_groups:
        print(f"  MISSING  {source_file}  {warn_class}  "
              f"{base_name}")
    return False


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Diff release-build warning logs against the "
                    "committed manifest of expected warnings, so "
                    "any new or vanished warning is loud.  Build "
                    "each variant from clean into its own log "
                    "(make imagoG and make imago separately) "
                    "before running this check.")
    parser.add_argument(
        "--real", metavar="LOG", default=None,
        help="build log of the real (gamma, imagoG) variant "
             "(default: not checked)")
    parser.add_argument(
        "--complex", metavar="LOG", dest="complex_log",
        default=None,
        help="build log of the complex (multi-k, imago) variant "
             "(default: not checked)")
    parser.add_argument(
        "--manifest", metavar="TSV", default=str(DEFAULT_MANIFEST),
        help=f"manifest file to check against or write "
             f"(default: {DEFAULT_MANIFEST})")
    parser.add_argument(
        "--write-manifest", action="store_true",
        help="write the manifest from the given logs instead of "
             "checking against it; requires BOTH --real and "
             "--complex so the manifest stays complete "
             "(default: check)")
    arguments = parser.parse_args(argv)

    if arguments.write_manifest:
        if not (arguments.real and arguments.complex_log):
            parser.error("--write-manifest requires both --real "
                         "and --complex logs")
    elif not (arguments.real or arguments.complex_log):
        parser.error("nothing to check: give --real and/or "
                     "--complex")

    logs_by_variant = {}
    if arguments.real:
        logs_by_variant["real"] = arguments.real
    if arguments.complex_log:
        logs_by_variant["complex"] = arguments.complex_log

    groups_by_variant = {}
    parse_noise = False
    for variant, log_path in logs_by_variant.items():
        groups, unattributed, unparsed = parse_build_log(log_path)
        groups_by_variant[variant] = groups
        for line in unattributed:
            print(f"{variant}: warning line had no preceding "
                  f"location line (kept under '<unknown>'): "
                  f"{line}")
            parse_noise = True
        for line in unparsed:
            print(f"{variant}: unrecognized warning line "
                  f"(compiler format drift?): {line}")
            parse_noise = True

    if arguments.write_manifest:
        write_manifest(arguments.manifest, groups_by_variant)
        total = sum(len(groups)
                    for groups in groups_by_variant.values())
        print(f"Wrote {total} warning-group rows to "
              f"{arguments.manifest}.")
        return 1 if parse_noise else 0

    expected_by_variant = read_manifest(arguments.manifest)
    all_clean = True
    for variant, found in groups_by_variant.items():
        expected = expected_by_variant.get(variant, set())
        if not report_variant(variant, found, expected):
            all_clean = False

    return 0 if (all_clean and not parse_noise) else 1


if __name__ == "__main__":
    sys.exit(main())
