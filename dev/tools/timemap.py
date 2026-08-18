#!/usr/bin/env python3
"""Turn imago's operation time stamps into a per-operation time map.

Every imago run writes, into its main output file (unit 20, the
``gs_scf-*.out`` / ``gs_pscf+*.out`` files), a start stamp and an
end stamp around each of the 34 named operations that
``O_TimeStamps`` (src/imago/timeStamps.f90) knows about: "Overlap
Integrals", "Secular Equation", "Make SCF Potential", and so on.
Those stamps were written as a progress indicator, but taken
together they are also a coarse profile that costs nothing to
collect and works on any run, including a multi-hour production
run that no profiler could follow.  This tool reads them back and
answers "where did the time go?" at the level of those
operations: how many times each was entered, how much wall time
it accumulated, and what share of the run that is.

This is measurement layer 2 of dev/PERFORMANCE.md ("Measurement
layers"), TODO PF2.  Its purpose is to TEST, before any code is
restructured, the two claims ARCHITECTURE 6.5 inherited from a
sibling code base and never measured on imago: that the
real-space grid loops in elecStat/exchCorr are the cheap first
target, and that the secular solve is the scaling wall.  If the
grid loops turn out to be five percent of the run, the first
target is the wrong first target, and this table is what says so.

The stamps look like this in the output file (a start block, the
operation's own printing, then an end block):

##    ***************************************************
##    ****************  Secular Equation  ***************
##    Date is: 2026/08/18 Time is: 13:25:15.298
##
##      ... whatever the operation prints ...
##
##    Date is: 2026/08/18 Time is: 13:25:15.410
##    ****************  Secular Equation  ***************
##    ***************************************************

The parser recognizes an operation START as a label line
immediately followed by a "Date is" line, and an operation END as
a "Date is" line immediately followed by the same label line.
Operations do not nest in imago today, but the parser keeps a
stack anyway so a nested pair would still be attributed
correctly rather than corrupting the table.

Two totals frame the table.  "Stamped span" is the wall time from
the first stamp to the last, i.e. the part of the run the stamps
cover.  "Unstamped" is the span minus the sum of all operations:
time spent between operations (file I/O, allocation, the parts of
the SCF loop nobody wrapped).  A large unstamped share is itself
a finding -- it means the stamps are missing something worth
naming -- and is why this tool prints it rather than hiding it.

Limits worth knowing.  Resolution is one millisecond, the
resolution of ``date_and_time``; a sub-millisecond operation
counts as zero.  The stamps are wall-clock, so on a shared node
they include whatever else was running; the baselines in
dev/PERFORMANCE.md are taken on an exclusive node for exactly
that reason.  A run that crosses midnight is handled (the date is
part of the stamp), but a stamp file from a run that died
mid-operation leaves that operation open; it is reported as
"unclosed" and excluded from the totals rather than guessed at.

Usage:
    timemap.py OUTFILE [OUTFILE ...]        one table per file
    timemap.py --csv OUTFILE                machine-readable rows
    timemap.py --calls OUTFILE              every call, in order

OUTFILE may also be a job directory, in which case the newest
``gs_*.out`` file in it is used.
"""

import argparse
import csv
import datetime
import glob
import os
import sys

## The banner and label lines are exactly 51 characters wide in the
## Fortran (format a51) and are made of asterisks around the label
## text; a label line is any 51-character line that starts and ends
## with an asterisk but is not the all-asterisk banner.
BANNER = "*" * 51
STAMP_PREFIX = "Date is: "


def is_label_line(line):
    """Return True when ``line`` is an operation label line.

    A label line is 51 characters wide, begins and ends with an
    asterisk, and carries text between the asterisks -- which is
    what separates it from the all-asterisk banner that frames it.
    """
    stripped = line.rstrip("\n")
    return (len(stripped) == 51 and stripped.startswith("*")
            and stripped.endswith("*") and stripped != BANNER)


def label_text(line):
    """Strip the asterisk padding from a label line: the op name."""
    return line.strip().strip("*").strip()


def parse_stamp(line):
    """Parse a ``Date is: YYYY/MM/DD Time is: HH:MM:SS.mmm`` line.

    Returns a datetime, or None when the line is not a stamp.
    """
    stripped = line.strip()
    if not stripped.startswith(STAMP_PREFIX):
        return None
    try:
        date_part = stripped[len(STAMP_PREFIX):len(STAMP_PREFIX) + 10]
        time_part = stripped.split("Time is:")[1].strip()
        return datetime.datetime.strptime(
            f"{date_part} {time_part}", "%Y/%m/%d %H:%M:%S.%f")
    except (ValueError, IndexError):
        return None


def parse_out_file(path):
    """Read one imago output file and return its list of calls.

    Each call is a dict with the operation label, its start and end
    datetimes, and its duration in seconds.  Unclosed operations
    (a start with no matching end, e.g. from a run that died) are
    returned separately so the caller can report them without
    letting them distort the totals.
    """
    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    calls = []
    open_stack = []          # [(label, start_datetime), ...]
    for index in range(len(lines) - 1):
        line = lines[index]
        next_line = lines[index + 1]

        # An operation START: a label line followed by a stamp.
        if is_label_line(line):
            start = parse_stamp(next_line)
            if start is not None:
                open_stack.append((label_text(line), start))
            continue

        # An operation END: a stamp followed by a label line.  Match
        #   it against the innermost open operation with that label
        #   so an (unexpected) nested pair still pairs correctly.
        end = parse_stamp(line)
        if end is not None and is_label_line(next_line):
            label = label_text(next_line)
            for depth in range(len(open_stack) - 1, -1, -1):
                if open_stack[depth][0] == label:
                    _, start = open_stack.pop(depth)
                    calls.append({
                        "label": label,
                        "start": start,
                        "end": end,
                        "seconds": (end - start).total_seconds(),
                    })
                    break

    unclosed = [label for label, _ in open_stack]
    return calls, unclosed


def summarize(calls):
    """Aggregate a call list into per-operation totals.

    Returns (rows, span_seconds, stamped_seconds) where rows is a
    list of dicts sorted by total time descending, span_seconds is
    first-stamp-to-last-stamp wall time, and stamped_seconds is the
    sum over all calls.
    """
    totals = {}
    for call in calls:
        row = totals.setdefault(call["label"], {
            "label": call["label"], "calls": 0, "seconds": 0.0,
            "max_seconds": 0.0})
        row["calls"] += 1
        row["seconds"] += call["seconds"]
        row["max_seconds"] = max(row["max_seconds"], call["seconds"])

    rows = sorted(totals.values(), key=lambda r: -r["seconds"])
    if calls:
        first = min(call["start"] for call in calls)
        last = max(call["end"] for call in calls)
        span_seconds = (last - first).total_seconds()
    else:
        span_seconds = 0.0
    stamped_seconds = sum(row["seconds"] for row in rows)
    return rows, span_seconds, stamped_seconds


def resolve_out_file(path):
    """Accept a file or a job directory; return the file to parse.

    For a directory, the newest ``gs_*.out`` inside it is chosen,
    which is the file imago.py copies the run's unit-20 output to.
    """
    if os.path.isdir(path):
        candidates = glob.glob(os.path.join(path, "gs_*.out"))
        if not candidates:
            sys.exit(f"timemap: no gs_*.out file in {path}")
        return max(candidates, key=os.path.getmtime)
    return path


def print_table(path, rows, span_seconds, stamped_seconds, unclosed):
    """Print the human-readable time map for one output file."""
    print(f"== {path}")
    if span_seconds <= 0.0:
        print("   no time stamps found")
        return
    print(f"   {'operation':<40} {'calls':>5} {'total s':>10} "
          f"{'% span':>7} {'mean s':>9} {'max s':>9}")
    for row in rows:
        share = 100.0 * row["seconds"] / span_seconds
        mean = row["seconds"] / row["calls"]
        print(f"   {row['label']:<40} {row['calls']:>5} "
              f"{row['seconds']:>10.3f} {share:>7.1f} {mean:>9.3f} "
              f"{row['max_seconds']:>9.3f}")
    unstamped = span_seconds - stamped_seconds
    print(f"   {'-- unstamped (between operations)':<40} {'':>5} "
          f"{unstamped:>10.3f} {100.0 * unstamped / span_seconds:>7.1f}")
    print(f"   {'== stamped span (first to last stamp)':<40} {'':>5} "
          f"{span_seconds:>10.3f} {100.0:>7.1f}")
    if unclosed:
        print(f"   !! unclosed operations (excluded): "
              f"{', '.join(unclosed)}")


def print_csv(path, rows, span_seconds, stamped_seconds):
    """Emit one CSV row per operation, plus the two framing rows."""
    writer = csv.writer(sys.stdout)
    for row in rows:
        writer.writerow([path, row["label"], row["calls"],
                         f"{row['seconds']:.3f}",
                         f"{100.0 * row['seconds'] / span_seconds:.2f}"
                         if span_seconds else "0"])
    writer.writerow([path, "unstamped", "",
                     f"{span_seconds - stamped_seconds:.3f}", ""])
    writer.writerow([path, "span", "", f"{span_seconds:.3f}", "100"])


def print_calls(path, calls):
    """List every call in the order it happened (for the curious)."""
    print(f"== {path}")
    for call in calls:
        print(f"   {call['start'].strftime('%H:%M:%S.%f')[:-3]}  "
              f"{call['seconds']:>10.3f}  {call['label']}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="timemap",
        description="Per-operation wall-time map from the O_TimeStamps "
                    "stamps in an imago output file (or a job "
                    "directory, whose newest gs_*.out is used).")
    parser.add_argument("paths", nargs="+", metavar="OUTFILE",
                        help="imago output file(s) or job "
                             "directory(ies)")
    parser.add_argument("--csv", action="store_true",
                        help="emit CSV rows (file, operation, calls, "
                             "seconds, percent) instead of a table. "
                             "Default: table.")
    parser.add_argument("--calls", action="store_true",
                        help="list every stamped call in order "
                             "instead of aggregating. Default: "
                             "aggregate.")
    args = parser.parse_args(argv)

    if args.csv:
        csv.writer(sys.stdout).writerow(
            ["file", "operation", "calls", "seconds", "percent"])
    for given in args.paths:
        path = resolve_out_file(given)
        calls, unclosed = parse_out_file(path)
        if args.calls:
            print_calls(path, calls)
            continue
        rows, span_seconds, stamped_seconds = summarize(calls)
        if args.csv:
            print_csv(path, rows, span_seconds, stamped_seconds)
        else:
            print_table(path, rows, span_seconds, stamped_seconds,
                        unclosed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
