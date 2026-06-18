#!/usr/bin/env python3

"""cod_fish.py -- Acquisition front-end for the Crystallography Open
Database (COD).

COD holds hundreds of thousands of experimentally determined crystal
structures, served as CIF files over plain HTTP.  This tool is a thin
front-end over COD's *own* endpoints -- it keeps no local database and
writes no SQL of its own -- with four verbs that separate browsing from
reproducible acquisition:

* ``get``    -- fetch one structure's CIF by its COD id, optionally
                verifying a pinned revision.  This is the reproducible
                path: a manifest pins ``(cod_id, cod_revision)`` and this
                fetch refuses any other revision.
* ``search`` -- find candidate structures by composition (and optionally
                author/text), printing a numbered table and saving the
                result set so later verbs can refer to rows by index.
* ``pin``    -- resolve chosen rows (by index, no eight-digit typing) to
                their current revisions and print a paste-ready manifest
                fragment -- the bridge from browsing to a pinned pull.
* ``rank``   -- advisory triage: when a composition has many entries,
                annotate and order them by interpretable signals so a
                student can spot the likely "real" phase.  It narrows the
                field; it does not decide.

The split matters for reproducibility: discovery (``search``/``rank`` and
the saved-index session) never feeds the build directly -- a manifest
only ever pins an id and a revision, which ``get`` fetches strictly.

Command-line usage::

    cod_fish.py search --elements Si              # single-element Si
    cod_fish.py search --elements Si O --max-extra 1
    cod_fish.py rank   --elements SiO2-as-elements ...
    cod_fish.py pin 2 5 7                          # by row index
    cod_fish.py get 9008463 --revision 291735 -o au.cif
"""

import argparse
import csv
import json
import os
import re
import sys
import urllib.parse
import urllib.request
from collections import Counter

COD_BASE = "https://www.crystallography.net/cod"

# The session file search writes (and pin/get read) so a chosen row can
#   be named by its small table index instead of its eight-digit id.
SESSION_FILE = "cod_fish_results.json"

# COD stamps each CIF's revision in the header as `$Revision: NNN $`.
_REVISION_RE = re.compile(r"\$Revision:\s*(\d+)\s*\$")


class CodFishError(Exception):
    """A COD request failed, or a pinned revision did not match.  Carries
    a message naming the specific failure."""


# ============================================================
#  HTTP (urllib only; no third-party request library)
# ============================================================

def _http_get(url, timeout=60):
    """GET ``url`` and return the raw bytes, or raise CodFishError.

    Strict on failure: a network outage or a COD error never returns a
    partial or substitute result."""

    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            return response.read()
    except Exception as exc:                       # noqa: BLE001
        raise CodFishError(f"COD request failed ({url}): {exc}") from exc


def cif_revision(cif_bytes):
    """Return the ``$Revision:$`` stamped in a CIF's header, or None."""

    match = _REVISION_RE.search(cif_bytes.decode("utf-8", "replace"))
    return match.group(1) if match else None


def fetch_cif(cod_id, revision=None, timeout=60):
    """Fetch one structure's CIF bytes by COD id.

    When ``revision`` is given the fetched CIF's stamped revision must
    match it exactly; a mismatch (the pinned revision has been superseded)
    is a hard error rather than a silent substitution, so a reproducible
    build never drifts from its pinned manifest."""

    data = _http_get(f"{COD_BASE}/{cod_id}.cif", timeout=timeout)
    if revision is not None:
        found = cif_revision(data)
        if found != str(revision):
            raise CodFishError(
                f"COD {cod_id}: pinned revision {revision} does not "
                f"match the served revision {found}; refusing to "
                f"substitute a different revision")
    return data


# ============================================================
#  Search (COD result.php; we only format what it returns)
# ============================================================

def search(elements, text=None, max_extra=0, limit=None, timeout=120):
    """Query COD's ``result.php`` for structures of a composition.

    ``elements`` is the list of element symbols that must be present.  By
    default the search is for *exactly* that composition: the
    distinct-element count is bounded to ``len(elements)``
    (``strictmin = strictmax = len``), so ``["Si"]`` returns single-
    element silicon phases, not every silicon-bearing compound.
    ``max_extra`` loosens the upper bound to admit up to that many further
    elements.  ``text`` is COD's free-text field (author names, etc.).

    Returns a list of candidate dicts (one per matching structure), in
    the order COD returns them.
    """

    params = [(f"el{index}", symbol)
              for index, symbol in enumerate(elements, start=1)]
    count = len(elements)
    params.append(("strictmin", str(count)))
    params.append(("strictmax", str(count + max_extra)))
    if text:
        params.append(("text", text))
    params.append(("format", "csv"))
    url = f"{COD_BASE}/result.php?" + urllib.parse.urlencode(params)
    text_body = _http_get(url, timeout=timeout).decode("utf-8", "replace")
    rows = _parse_cod_csv(text_body)
    return rows[:limit] if limit else rows


def _parse_cod_csv(text_body):
    """Parse COD's CSV result into candidate dicts, keeping the few
    columns the table and the triage need."""

    data_lines = [line for line in text_body.splitlines()
                  if not line.startswith("#")]
    candidates = []
    for record in csv.DictReader(data_lines):
        cod_id = (record.get("file") or "").strip()
        if not cod_id:
            continue
        candidates.append({
            "id": cod_id,
            "formula": (record.get("formula")
                        or record.get("calcformula") or "").strip(),
            "name": (record.get("mineral") or record.get("commonname")
                     or record.get("chemname") or "").strip(),
            "sg": (record.get("sg") or "").strip(),
            "sgnum": (record.get("sgNumber") or "").strip(),
            "a": record.get("a"), "b": record.get("b"),
            "c": record.get("c"),
            "alpha": record.get("alpha"), "beta": record.get("beta"),
            "gamma": record.get("gamma"),
            "vol": record.get("vol"),
            "celltemp": record.get("celltemp"),
            "diffrtemp": record.get("diffrtemp"),
            "cellpressure": record.get("cellpressure"),
            "diffrpressure": record.get("diffrpressure"),
        })
    return candidates


# ============================================================
#  Advisory triage (rank)
# ============================================================

def _as_float(value):
    """Parse a COD numeric cell, returning None for blank/garbage."""

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _is_ambient(candidate):
    """True when nothing marks the structure as off-ambient.

    COD usually leaves temperature/pressure blank for room-condition
    measurements, so a blank field is read as ambient.  A recorded
    temperature far from room temperature, or any recorded non-trivial
    pressure, marks it special."""

    for key in ("celltemp", "diffrtemp"):
        temperature = _as_float(candidate.get(key))
        if temperature is not None and not 250.0 <= temperature <= 320.0:
            return False
    for key in ("cellpressure", "diffrpressure"):
        pressure = _as_float(candidate.get(key))
        # COD records pressure in kPa; ~101 kPa is one atmosphere.
        if pressure is not None and pressure > 200.0:
            return False
    return True


def rank(candidates):
    """Order candidates by interpretable signals and return a list of
    ``(candidate, reasons)`` tuples, best first.

    The signals are drawn only from COD metadata -- a named mineral /
    common phase, ambient measurement conditions, agreement with the
    spacegroup the *plurality* of candidates share, and whether the cell
    volume is a gross outlier from the median.  This is triage to narrow
    "which of these is the real phase," not a validation verdict, so the
    reasons are reported alongside every row.
    """

    spacegroup_counts = Counter(
        c["sgnum"] for c in candidates if c["sgnum"])
    consensus_sg = (spacegroup_counts.most_common(1)[0][0]
                    if spacegroup_counts else None)
    volumes = [v for v in (_as_float(c["vol"]) for c in candidates)
               if v is not None]
    median_volume = (sorted(volumes)[len(volumes) // 2]
                     if volumes else None)

    scored = []
    for candidate in candidates:
        score = 0
        reasons = []
        if candidate["name"]:
            score += 2
            reasons.append(f"named:{candidate['name']}")
        if _is_ambient(candidate):
            score += 2
            reasons.append("ambient")
        else:
            score -= 2
            reasons.append("non-ambient")
        if consensus_sg and candidate["sgnum"] == consensus_sg:
            score += 1
            reasons.append(f"sg{candidate['sgnum']}(consensus)")
        volume = _as_float(candidate["vol"])
        if median_volume and volume is not None and (
                volume > 2.0 * median_volume
                or volume < 0.5 * median_volume):
            score -= 2
            reasons.append("volume-outlier")
        scored.append((score, candidate, reasons))
    scored.sort(key=lambda item: -item[0])
    return [(candidate, reasons) for _, candidate, reasons in scored]


# ============================================================
#  Session file (index <-> cod_id) and pinning
# ============================================================

def save_session(candidates, query, path=SESSION_FILE):
    """Write the search result set so ``pin``/``get`` can name a row by
    its table index."""

    with open(path, "w") as handle:
        json.dump({"query": query, "rows": candidates}, handle, indent=2)


def load_session(path=SESSION_FILE):
    """Return the saved search rows, or None when there is no session."""

    if not os.path.exists(path):
        return None
    with open(path) as handle:
        return json.load(handle).get("rows")


def _expand_index_tokens(tokens):
    """Expand index tokens like ``2`` and ranges like ``1-3`` into a flat
    integer list, leaving non-range tokens as-is for the caller."""

    expanded = []
    for token in tokens:
        if "-" in token and not token.startswith("-"):
            low, high = token.split("-", 1)
            expanded.extend(str(n) for n in range(int(low), int(high) + 1))
        else:
            expanded.append(token)
    return expanded


def resolve_ids(tokens, session_rows):
    """Map a mix of row indices and raw COD ids to COD ids.

    A small integer (<= the number of saved rows) is read as a 1-based
    index into the saved search; anything larger is taken as a raw COD id
    (COD ids are seven digits, far above any table index).  Without a
    saved session every token must be a raw id.
    """

    row_count = len(session_rows) if session_rows else 0
    ids = []
    for token in _expand_index_tokens(tokens):
        if not token.isdigit():
            raise CodFishError(f"not a row index or COD id: {token!r}")
        number = int(token)
        if session_rows and number <= row_count:
            ids.append(session_rows[number - 1]["id"])
        else:
            ids.append(token)
    return ids


def pin(tokens, session_rows=None):
    """Resolve the chosen rows to their current revisions.

    Returns a list of ``{id, revision}`` dicts.  Only the chosen few are
    fetched (to read their revisions); a broad search is never downloaded
    wholesale."""

    pinned = []
    for cod_id in resolve_ids(tokens, session_rows):
        data = fetch_cif(cod_id)
        pinned.append({"id": cod_id, "revision": cif_revision(data)})
    return pinned


# ============================================================
#  Command-line presentation
# ============================================================

def _print_table(candidates):
    """Print the numbered candidate table."""

    header = (f"{'#':>3}  {'cod_id':>8}  {'formula':<16} "
              f"{'spacegroup':<12} {'name':<16} "
              f"{'a':>7} {'b':>7} {'c':>7}")
    print(header)
    for index, c in enumerate(candidates, start=1):
        print(f"{index:>3}  {c['id']:>8}  {c['formula'][:16]:<16} "
              f"{(c['sg'] or '?')[:12]:<12} {c['name'][:16]:<16} "
              f"{(c['a'] or '?'):>7} {(c['b'] or '?'):>7} "
              f"{(c['c'] or '?'):>7}")


def _manifest_fragment(pinned):
    """Render pinned entries as a paste-ready manifest fragment."""

    lines = []
    for entry in pinned:
        lines.append("[[reference_solid]]")
        lines.append(f'reference_id = "cod_{entry["id"]}"')
        lines.append('system_type = "crystalline"')
        lines.append(f'cod_id = {entry["id"]}')
        lines.append(f'cod_revision = "{entry["revision"]}"')
        lines.append("# fill in: kpoint_spec, scf_threshold, basis, "
                     "functional, kpoint_integration, and entries")
        lines.append("")
    return "\n".join(lines)


def _cmd_get(args):
    data = fetch_cif(args.cod_id, revision=args.revision)
    out = args.output or f"{args.cod_id}.cif"
    with open(out, "wb") as handle:
        handle.write(data)
    print(f"Wrote {out} (COD {args.cod_id}, revision "
          f"{cif_revision(data)}).")
    return 0


def _cmd_search(args):
    rows = search(args.elements, text=args.text,
                  max_extra=args.max_extra, limit=args.limit)
    if not rows:
        print("No matching structures.")
        return 0
    if args.rank:
        ordered = rank(rows)
        rows = [candidate for candidate, _ in ordered]
        reasons_by_id = {c["id"]: r for c, r in ordered}
    _print_table(rows)
    if args.rank:
        print("\nTriage (advisory):")
        for index, candidate in enumerate(rows, start=1):
            print(f"  [{index}] {candidate['id']}  "
                  f"{' · '.join(reasons_by_id[candidate['id']])}")
    save_session(rows, vars(args).get("elements"))
    print(f"\nSaved {len(rows)} rows to {SESSION_FILE}; "
          f"pin them by index, e.g. `cod_fish.py pin 1 3`.")
    return 0


def _cmd_rank(args):
    args.rank = True
    return _cmd_search(args)


def _cmd_pin(args):
    pinned = pin(args.tokens, session_rows=load_session())
    print(_manifest_fragment(pinned))
    return 0


def _build_parser():
    """Build the cod_fish argument parser (four sub-commands)."""

    parser = argparse.ArgumentParser(
        description="Fetch and search crystal structures from the "
                    "Crystallography Open Database.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    getp = subparsers.add_parser(
        "get", help="Fetch one structure's CIF by COD id.")
    getp.add_argument("cod_id", help="The COD id to fetch.")
    getp.add_argument(
        "-r", "--revision", default=None,
        help="Pinned COD revision to verify; a mismatch is a hard "
             "error.  Optional; defaults to no revision check.")
    getp.add_argument(
        "-o", "--output", default=None,
        help="Output CIF path.  Optional; defaults to <cod_id>.cif.")
    getp.set_defaults(func=_cmd_get)

    def add_search_args(sub):
        sub.add_argument(
            "-e", "--elements", nargs="+", required=True, metavar="EL",
            help="Element symbols that must be present.  The search is "
                 "for exactly this composition unless --max-extra is "
                 "given.")
        sub.add_argument(
            "--max-extra", type=int, default=0, metavar="K",
            help="Allow up to K elements beyond those listed.  "
                 "Optional; defaults to 0 (exact composition).")
        sub.add_argument(
            "--text", default=None,
            help="COD free-text filter (e.g. an author name).  "
                 "Optional; defaults to none.")
        sub.add_argument(
            "--limit", type=int, default=None, metavar="N",
            help="Show at most N candidates.  Optional; defaults to "
                 "all matches.")

    searchp = subparsers.add_parser(
        "search", help="Find candidate structures by composition.")
    add_search_args(searchp)
    searchp.add_argument(
        "--rank", action="store_true", default=False,
        help="Also print advisory triage of the candidates.  "
             "Optional; defaults to off.")
    searchp.set_defaults(func=_cmd_search)

    rankp = subparsers.add_parser(
        "rank", help="Search and print advisory triage (search "
                     "--rank).")
    add_search_args(rankp)
    rankp.set_defaults(func=_cmd_rank)

    pinp = subparsers.add_parser(
        "pin", help="Resolve chosen rows (by index or id) to pinned "
                    "revisions and print a manifest fragment.")
    pinp.add_argument(
        "tokens", nargs="+", metavar="INDEX|ID",
        help="Row indices from the last search (ranges like 1-3 are "
             "allowed) or raw COD ids.")
    pinp.set_defaults(func=_cmd_pin)

    return parser


def main(argv=None):
    """Entry point for the ``cod_fish.py`` command line."""

    args = _build_parser().parse_args(argv)
    try:
        return args.func(args)
    except CodFishError as error:
        print(f"cod_fish: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
