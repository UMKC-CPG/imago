#!/usr/bin/env python3
"""Measure how far apart two imago scratch files are.

PSEUDOCODE 31.7 accepts the rank-k density build against a tolerance
that is MEASURED rather than chosen: the largest absolute difference
between the scratch file a recast binary writes and the one the
pre-recast binary writes on the same deck.  ``h5diff`` can only answer
"is every element within delta?", which turns finding the floor into a
bisection over deltas; this tool reads both files and reports the
floor directly.

For every dataset present in both files it prints the number of
elements, the largest absolute difference, and the largest relative
difference (absolute difference over the larger magnitude, elements
below 1e-300 in both files skipped), then the overall maximum.  Paths
whose name contains any excluded substring are skipped: eigenvectors
are excluded by default, because a degenerate subspace lets two
correct runs choose different bases and their difference measures
nothing about the density.

Datasets present in only one file, or differing in shape, are
reported as such and count as a difference.  The exit status is 0
when the overall maximum absolute difference is within --tolerance
(default 0, meaning bit-identical) and 1 otherwise, so the tool can
gate a batch script the way ``h5diff -q`` does.

Usage:
    h5maxdiff.py FILE_A FILE_B [--exclude SUBSTRING ...]
                 [--tolerance DELTA] [--top N]
"""

import argparse
import sys

import h5py
import numpy


def collect_datasets(hdf5_file):
    """Return {path: dataset} for every dataset in the file."""
    found = {}

    def visit(name, node):
        if isinstance(node, h5py.Dataset):
            found[name] = node

    hdf5_file.visititems(visit)
    return found


def compare(path_a, path_b, excludes, top):
    """Print per-dataset differences and return the overall maximum.

    Returns None when a structural difference (missing dataset or
    shape mismatch) makes a numerical maximum meaningless.
    """
    structural_difference = False
    rows = []
    with h5py.File(path_a, "r") as file_a, h5py.File(path_b, "r") as file_b:
        sets_a = collect_datasets(file_a)
        sets_b = collect_datasets(file_b)
        names = sorted(set(sets_a) | set(sets_b))
        for name in names:
            if any(substring in name for substring in excludes):
                continue
            if name not in sets_a or name not in sets_b:
                which = "A" if name in sets_a else "B"
                print(f"ONLY IN {which}: {name}")
                structural_difference = True
                continue
            data_a = sets_a[name][()]
            data_b = sets_b[name][()]
            if data_a.shape != data_b.shape:
                print(f"SHAPE {data_a.shape} vs {data_b.shape}: {name}")
                structural_difference = True
                continue
            if data_a.dtype.kind not in "fc" and data_b.dtype.kind not in "fc":
                # Integer or string data: any difference is structural.
                if not numpy.array_equal(data_a, data_b):
                    print(f"NON-NUMERIC DIFFERENCE: {name}")
                    structural_difference = True
                continue
            values_a = numpy.asarray(data_a, dtype=numpy.float64).ravel()
            values_b = numpy.asarray(data_b, dtype=numpy.float64).ravel()
            absolute = numpy.abs(values_a - values_b)
            scale = numpy.maximum(numpy.abs(values_a), numpy.abs(values_b))
            nonzero = scale > 1e-300
            relative = numpy.zeros_like(absolute)
            relative[nonzero] = absolute[nonzero] / scale[nonzero]
            rows.append((float(absolute.max()) if absolute.size else 0.0,
                         float(relative.max()) if relative.size else 0.0,
                         absolute.size, name))

    rows.sort(reverse=True)
    print(f"{'max |a-b|':>12} {'max rel':>12} {'elements':>10}  dataset")
    for max_abs, max_rel, count, name in rows[:top]:
        print(f"{max_abs:12.4e} {max_rel:12.4e} {count:10d}  {name}")
    if len(rows) > top:
        print(f"... {len(rows) - top} more datasets, all smaller")
    overall = max((row[0] for row in rows), default=0.0)
    note = " (STRUCTURAL DIFFERENCES ABOVE)" if structural_difference else ""
    print(f"OVERALL max |a-b| = {overall:.4e} over {len(rows)} datasets{note}")
    return None if structural_difference else overall


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Largest absolute difference between two HDF5 files, "
                    "dataset by dataset.")
    parser.add_argument("file_a")
    parser.add_argument("file_b")
    parser.add_argument("--exclude", action="append",
                        default=["eigenVectors"],
                        help="skip datasets whose path contains this "
                             "substring (repeatable; default: eigenVectors)")
    parser.add_argument("--tolerance", type=float, default=0.0,
                        help="largest acceptable absolute difference "
                             "(default 0.0, i.e. bit-identical)")
    parser.add_argument("--top", type=int, default=12,
                        help="datasets to list, largest first (default 12)")
    arguments = parser.parse_args(argv)
    overall = compare(arguments.file_a, arguments.file_b,
                      arguments.exclude, arguments.top)
    if overall is None or overall > arguments.tolerance:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
