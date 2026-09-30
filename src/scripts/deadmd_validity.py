#!/usr/bin/env python3
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Mohammed Belhadj Larbi

"""
deadmd_validity.py -- pass/fail check of a condensed model before it is scored.

USAGE (inside a member's lammps/ directory):
    deadmd_validity.py [-dir .] [-valence C 4 6] [-valence B 5 6 ...]
                     [-max_strain 0.30]

A model that fails is not a WORSE model, it is an INVALID one (see
dead-md_dev/DIAGNOSIS.md, ARCHITECTURE: "gate first, score second"), so
deadmd.py gives it the crash sentinel instead of a score.

Reads final.data (write_data after the final minimisation), lammps.dat (the
starting atom count), lammps.in (which map file belongs to which RXN),
the reaction map files (atoms each event deletes) and rxn_counts_evolve
(events per reaction), all written by condense.py.

HARD checks (any failure makes the model invalid):
  coordination  every atom's bond count is allowed for its element.
                Defaults: H 1, C 4, N 3, O 2, F 1, Si 4, P 3 4, S 2,
                Cl 1, Ge 4, B 3 4 5 6.  Override per element with
                -valence (e.g. a-BC:H cage carbon: -valence C 4 6).
                Catches reactions that delete their leaving atoms but never
                write their new bond (e.g. the precursorDB si-1 si-1
                template: Si left with 3 bonds).
  bookkeeping   final atoms = start atoms - sum_r deleted_r x events_r.
  strain        every bond is within +-max_strain of its Bond Coeffs r0.
REPORTED only (never fail the gate):
  fragments (connected pieces), missing / extra angles relative to the
  bond topology (linkage angles are deliberately floppy in some runs), and
  non-bonded close contacts (< 1.0 A).

Writes validity.txt; exit code 0 = valid, 1 = invalid.
"""

import argparse
import os
import re
import sys
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np

DEFAULT_VALENCE = {
    "H": {1}, "C": {4}, "N": {3}, "O": {2}, "F": {1}, "Si": {4},
    "P": {3, 4}, "S": {2}, "Cl": {1}, "Ge": {4}, "B": {3, 4, 5, 6},
}
MASS_TO_ELEMENT = {
    1.008: "H", 10.811: "B", 12.011: "C", 14.007: "N", 15.999: "O",
    18.998: "F", 28.086: "Si", 30.974: "P", 32.06: "S", 35.45: "Cl",
    72.63: "Ge",
}


def element_of(mass):
    best = min(MASS_TO_ELEMENT, key=lambda m: abs(m - mass))
    if abs(best - mass) > 0.1:
        sys.exit(f"deadmd_validity.py: no element with mass {mass}")
    return MASS_TO_ELEMENT[best]


def read_data(path):
    """Box lengths, per-atom element and position, bonds, angles, bond r0."""
    mass, typ, pos, bonds, angles, r0 = {}, {}, {}, [], set(), {}
    box, tilt, sec = [], False, None
    for line in open(path):
        if "xlo xhi" in line or "ylo yhi" in line or "zlo zhi" in line:
            lo, hi = map(float, line.split()[:2])
            box.append(hi - lo)
            continue
        if "xy xz yz" in line:
            tilt = True
        words = line.split("#")[0].split()
        if not words:
            continue
        if line[0].isalpha():
            sec = "BondCoeffs" if line.startswith("Bond Coeffs") else words[0]
            continue
        if sec == "Masses":
            mass[int(words[0])] = float(words[1])
        elif sec == "BondCoeffs":
            r0[int(words[0])] = float(words[2])
        elif sec == "Atoms":
            atom = int(words[0])
            typ[atom] = int(words[2])
            pos[atom] = np.array(list(map(float, words[4:7])))
        elif sec == "Bonds":
            bonds.append((int(words[1]), int(words[2]), int(words[3])))
        elif sec == "Angles":
            i, j, k = int(words[2]), int(words[3]), int(words[4])
            angles.add((min(i, k), j, max(i, k)))
    if tilt:
        sys.exit("deadmd_validity.py: triclinic box not supported")
    element = {a: element_of(mass[t]) for a, t in typ.items()}
    return np.array(box), element, pos, bonds, angles, r0


def expected_atoms(directory):
    """Start atoms minus the atoms deleted by every recorded event."""
    with open(os.path.join(directory, "lammps.dat")) as f:
        start = next(int(l.split()[0]) for l in f if l.strip().endswith("atoms"))
    maps = re.findall(r"react RXN(\d+) all \S+ \S+ \S+ \S+ \S+ (\S+\.map)",
                      open(os.path.join(directory, "lammps.in")).read())
    deleted = {}
    for rxn, map_file in maps:
        with open(os.path.join(directory, map_file)) as f:
            m = re.search(r"(\d+)\s+deleteIDs", f.read())
        deleted[int(rxn)] = int(m.group(1)) if m else 0
    with open(os.path.join(directory, "rxn_counts_evolve")) as f:
        counts = [int(float(x)) for x in f.read().split()[-len(deleted):]]
    return start - sum(deleted[r + 1] * n for r, n in enumerate(counts)), \
        start, deleted, counts


def gate(directory=".", valence=None, max_strain=0.30):
    """Return (valid, reasons, report lines) for one member's lammps dir."""
    allowed = {k: set(v) for k, v in DEFAULT_VALENCE.items()}
    allowed.update(valence or {})
    box, element, pos, bonds, angles, r0 = read_data(
        os.path.join(directory, "final.data"))

    def dist(a, b):
        d = pos[a] - pos[b]
        d -= box * np.round(d / box)
        return float(np.linalg.norm(d))

    reasons, report = [], []
    adj = defaultdict(set)
    for _, a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)

    # Coordination.
    bad = Counter((element[a], len(adj[a])) for a in element
                  if element[a] in allowed
                  and len(adj[a]) not in allowed[element[a]])
    if bad:
        reasons.append("coordination: " + ", ".join(
            f"{n} {e} with {k} bonds" for (e, k), n in sorted(bad.items())))
    report.append("coordination defects: " + (", ".join(
        f"{e}{k}:{n}" for (e, k), n in sorted(bad.items())) or "none"))

    # Bookkeeping.
    want, start, deleted, counts = expected_atoms(directory)
    report.append(f"bookkeeping: {len(element)} atoms, expected {want} "
                  f"(start {start}, events {counts}, deleted per event "
                  f"{[deleted[r] for r in sorted(deleted)]})")
    if len(element) != want:
        reasons.append(f"bookkeeping: {len(element)} atoms, expected {want}")

    # Bond strain against each bond type's r0.
    strain = [(dist(a, b) - r0[t]) / r0[t] for t, a, b in bonds]
    worst = max(strain, key=abs)
    n_out = sum(abs(s) > max_strain for s in strain)
    report.append(f"bond strain: worst {worst:+.1%}, {n_out} beyond "
                  f"+-{max_strain:.0%}")
    if n_out:
        reasons.append(f"strain: {n_out} bonds beyond +-{max_strain:.0%} "
                       f"(worst {worst:+.1%})")

    # Reported only: fragments, angle audit, close contacts.
    seen, frags = set(), 0
    for a in element:
        if a in seen:
            continue
        frags += 1
        stack = [a]
        while stack:
            u = stack.pop()
            if u not in seen:
                seen.add(u)
                stack.extend(adj[u])
    implied = {(min(i, k), j, max(i, k))
               for j in adj for i, k in combinations(adj[j], 2)}
    atoms = sorted(element)
    xyz = np.array([pos[a] for a in atoms])
    close = 0
    for n, a in enumerate(atoms[:-1]):
        d = xyz[n + 1:] - xyz[n]
        d -= box * np.round(d / box)
        near = np.where(np.einsum("ij,ij->i", d, d) < 1.0)[0]
        close += sum(1 for m in near if atoms[n + 1 + m] not in adj[a])
    report.append(f"fragments: {frags}; angles: {len(implied - angles)} "
                  f"missing / {len(angles - implied)} extra; non-bonded "
                  f"contacts < 1.0 A: {close}")
    return not reasons, reasons, report


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("-dir", default=".")
    ap.add_argument("-valence", nargs="+", action="append", default=[],
                    metavar="EL N", help="allowed bond counts, e.g. C 4 6")
    ap.add_argument("-max_strain", type=float, default=0.30)
    args = ap.parse_args()
    valence = {v[0].capitalize(): {int(x) for x in v[1:]}
               for v in args.valence}
    valid, reasons, report = gate(args.dir, valence, args.max_strain)
    text = "\n".join(["VALID" if valid else "INVALID"] + reasons + report)
    with open(os.path.join(args.dir, "validity.txt"), "w") as f:
        f.write(text + "\n")
    print(text)
    sys.exit(0 if valid else 1)


if __name__ == "__main__":
    main()
