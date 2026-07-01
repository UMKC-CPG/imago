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
                their current revisions and print a complete, runnable
                manifest (or, with ``--sketch-only``, a sketch for
                expand_manifest to finish) -- the bridge from browsing to
                a pinned pull.
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
import math
import os
import re
import sys
import urllib.parse
from datetime import datetime
import urllib.request
from collections import Counter

# The default ``pin`` output is a complete, runnable manifest: it is
#   built as a CurationManifest and rendered through the shared writer,
#   so the [defaults] run settings and the [characterization] recipe
#   come from one place (curation_manifest) and cannot drift from
#   expand_manifest or the producer (DESIGN 5.7).
from curation_manifest import (
    CurationManifest, ReferenceSolid,
    default_run_settings, default_characterization, format_manifest)

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

def search(elements, text=None, max_extra=0, limit=None, fuzzy=False,
           timeout=120):
    """Query COD's ``result.php`` for structures of a composition.

    ``elements`` is the list of element specs that must be present.  Each
    is an element symbol with an *optional* count on either side --
    ``"Fe"``, ``"2Fe"``, or ``"Fe2"`` (see :func:`_parse_element_spec`).
    By default the search is for *exactly* the listed element set: the
    distinct-element count is bounded to the number of symbols
    (``strictmin = strictmax = N``), so ``["Si"]`` returns single-element
    silicon phases, not every silicon-bearing compound.  ``max_extra``
    loosens the upper bound to admit up to that many further elements.
    ``text`` is COD's free-text field (author names, etc.).

    When any spec carries a count, how it is applied depends on
    ``fuzzy``:

    * **Default (exact).**  COD is queried by element *set* only, and the
      results are filtered here, in the client, to those whose *formula
      unit* (GCD-reduced formula) has exactly the requested count of each
      element.  This gives the clean "N per formula unit" meaning a human
      reads off a formula -- it includes Z-multiplied cells (an
      ``Fe4 O6`` entry reduces to ``Fe2 O3``) and excludes
      non-stoichiometric entries (``Fe0.911 O``).
    * **Fuzzy (``fuzzy=True``).**  The counts are handed to COD's own
      formula match (the ``"Fe 2"`` / ``" O 3"`` box form) and whatever
      COD returns is taken as-is, with no client reduction.  COD's match
      is a digit-prefix match on the stored formula string, so this is
      broader and noisier -- ``"Fe 2"`` also matches ``Fe2.645``,
      ``Fe21``, ``Fe25`` -- but it asks nothing of the formula's
      reducibility, which is useful for exploring messy or
      non-stoichiometric data.

    Returns a list of candidate dicts (one per matching structure), in
    the order COD returns them.
    """

    specs = [_parse_element_spec(token) for token in elements]
    symbols = [symbol for _, symbol in specs]
    wanted = {symbol: count for count, symbol in specs
              if count is not None}

    # Build the COD element boxes.  Default sends bare symbols (counts
    #   are applied exactly, client-side); fuzzy sends each counted
    #   element in COD's own "<sym> <count>" formula-box form.
    params = []
    for index, (count, symbol) in enumerate(specs, start=1):
        if fuzzy and count is not None:
            params.append((f"el{index}", _cod_count_box(symbol, count)))
        else:
            params.append((f"el{index}", symbol))
    params.append(("strictmin", str(len(symbols))))
    params.append(("strictmax", str(len(symbols) + max_extra)))
    if text:
        params.append(("text", text))
    params.append(("format", "csv"))
    url = f"{COD_BASE}/result.php?" + urllib.parse.urlencode(params)
    text_body = _http_get(url, timeout=timeout).decode("utf-8", "replace")
    rows = _parse_cod_csv(text_body)
    if wanted and not fuzzy:
        rows = [row for row in rows
                if _matches_stoichiometry(row, wanted)]
    return rows[:limit] if limit else rows


def _cod_count_box(symbol, count):
    """Format an element+count for COD's formula box (tips.php): the
    symbol, a space, then the count -- with a *leading* space too when
    the symbol is one letter, so COD does not prefix-match it into a
    two-letter element (``" O 3"`` not ``"O3"``)."""

    leading = " " if len(symbol) == 1 else ""
    return f"{leading}{symbol} {count}"


def _parse_element_spec(token):
    """Parse one ``--elements`` token into ``(count, symbol)``.

    The count is optional and may sit on either side of the symbol, so a
    user can write whichever reads naturally: ``"Fe"`` -> ``(None,
    "Fe")``, ``"2Fe"`` -> ``(2, "Fe")``, ``"Fe2"`` -> ``(2, "Fe")``.  A
    count on both sides (``"2Fe2"``) or any other shape is rejected.  The
    symbol is normalized to element capitalization (``"fe"`` -> ``"Fe"``).
    """

    text = token.strip()
    leading = re.fullmatch(r"(\d+)([A-Za-z]{1,2})", text)    # 2Fe
    trailing = re.fullmatch(r"([A-Za-z]{1,2})(\d+)", text)   # Fe2
    bare = re.fullmatch(r"([A-Za-z]{1,2})", text)            # Fe
    if leading:
        return int(leading.group(1)), _normalize_symbol(leading.group(2))
    if trailing:
        return int(trailing.group(2)), _normalize_symbol(
            trailing.group(1))
    if bare:
        return None, _normalize_symbol(bare.group(1))
    raise CodFishError(
        f"bad element spec {token!r}; use a symbol with an optional "
        f"count on one side, e.g. Fe, 2Fe, or Fe2")


def _normalize_symbol(symbol):
    """Capitalize an element symbol the conventional way (Fe, O, Si)."""

    return symbol[:1].upper() + symbol[1:].lower()


def _formula_counts(formula):
    """Parse a COD formula string (e.g. ``"- Fe2 O3 -"``) into a
    ``{symbol: count}`` map.  Counts may be fractional in the source
    (non-stoichiometric entries like ``"Fe0.911 O"``); they are kept as
    floats here and judged for integrality in :func:`_reduced_formula`."""

    counts = {}
    for symbol, number in re.findall(r"([A-Z][a-z]?)([0-9.]*)", formula):
        counts[symbol] = counts.get(symbol, 0.0) + (
            float(number) if number else 1.0)
    return counts


def _reduced_formula(counts):
    """Reduce a ``{symbol: count}`` map to its formula unit (divide by
    the GCD of the counts), or return None when any count is not an
    integer -- a non-stoichiometric entry has no integer formula unit and
    so cannot match an exact per-element request."""

    integral = {}
    for symbol, value in counts.items():
        if abs(value - round(value)) > 1.0e-6:
            return None
        integral[symbol] = int(round(value))
    divisor = 0
    for value in integral.values():
        divisor = math.gcd(divisor, value)
    if divisor > 1:
        integral = {symbol: value // divisor
                    for symbol, value in integral.items()}
    return integral


def _matches_stoichiometry(candidate, wanted):
    """True when the candidate's formula unit has exactly the requested
    count of each counted element.  Elements the user left bare impose no
    count constraint (only presence, already handled by the query)."""

    reduced = _reduced_formula(_formula_counts(candidate.get("formula")
                                               or ""))
    if not reduced:
        return False
    return all(reduced.get(symbol) == count
               for symbol, count in wanted.items())


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


# ============================================================
#  CIF metadata -> reference_id (DESIGN 5.7)
# ============================================================
#
# pin already downloads each chosen CIF (to read its revision), so it
# also derives a human-meaningful, label-safe reference_id from the
# CIF's own metadata: the reduced chemical formula, the space-group
# Hermann-Mauguin symbol and International Tables number, and the
# publication year -- e.g. "si_fd-3m_227_2010".  This frees the
# curator from inventing an id for every pinned structure.  The year
# also separates two phases that share a space group (4H and hcp
# silicon are both P 63/m m c); any residual collision is broken with
# a trailing counter by the fragment emitter.


def _cif_value(text, *tags):
    """Return the value of the first CIF data item matching any of
    ``tags`` (a single-line tag/value pair), stripped of quotes, or
    "" when none is present.  The exact first-token match avoids
    picking up a longer tag that merely starts with the same text."""

    for tag in tags:
        for line in text.splitlines():
            parts = line.strip().split(None, 1)
            if parts and parts[0] == tag:
                return parts[1].strip().strip("'\"") \
                    if len(parts) > 1 else ""
    return ""


def _sanitize_token(text):
    """Reduce a string to the label-safe charset (lowercase letters,
    digits, hyphen) a reference_id is restricted to (DESIGN 5.2.1)."""

    return re.sub(r"[^a-z0-9-]", "", text.lower())


def _formula_token(text):
    """Build the composition part of a reference_id from the CIF's
    ``_chemical_formula_sum`` -- the reduced formula unit, elements in
    alphabetical order, counts above one kept (Si -> "si", Fe2 O3 ->
    "fe2o3").  Returns "" when no integer formula unit is available."""

    counts = _reduced_formula(
        _formula_counts(_cif_value(text, "_chemical_formula_sum")))
    if not counts:
        return ""
    return "".join(
        symbol.lower() + (str(counts[symbol]) if counts[symbol] > 1
                          else "")
        for symbol in sorted(counts))


def _space_group_token(text):
    """Return the (Hermann-Mauguin symbol, International Tables number)
    tokens for a reference_id.  The H-M symbol is reduced to the
    label-safe charset with the origin/setting suffix dropped (so
    "F d -3 m :1" -> "fd-3m" and "P 63/m m c" -> "p63mmc")."""

    hm = _cif_value(text, "_symmetry_space_group_name_H-M",
                    "_space_group_name_H-M_alt")
    it = _cif_value(text, "_space_group_IT_number",
                    "_symmetry_Int_Tables_number")
    return _sanitize_token(hm.split(":")[0]), _sanitize_token(it)


def _auto_reference_id(cif_bytes):
    """Derive a label-safe reference_id from a CIF's metadata --
    ``<formula>_<H-M symbol>_<IT number>_<year>``, e.g.
    "si_fd-3m_227_2010".  The publication year (always present and
    numeric) both dates the entry and separates two phases that share
    a space group.  Returns None when the formula or space group
    cannot be read, so the caller falls back to a cod-id-based name."""

    text = cif_bytes.decode("utf-8", "replace")
    formula = _formula_token(text)
    hm, it = _space_group_token(text)
    if not (formula and hm and it):
        return None
    parts = [formula, hm, it]
    year = _sanitize_token(_cif_value(text, "_journal_year"))
    if year:
        parts.append(year)
    return "_".join(parts)


def _composition(cif_bytes):
    """The element symbols present in the CIF's chemical formula,
    sorted -- recorded in the sketch so the authoring tool can
    auto-fill each entry's element without re-reading the CIF."""

    text = cif_bytes.decode("utf-8", "replace")
    return sorted(_formula_counts(
        _cif_value(text, "_chemical_formula_sum")))


def _source_description(cif_bytes):
    """A human description of the pinned structure from its CIF
    metadata -- a name (the common, mineral, or systematic chemical
    name, else the formula) plus the space-group H-M symbol, its
    International Tables number, and the publication year, e.g.
    "hexagonal 4H Si-IV, P 63/m m c (194), 2018".  Recorded in the
    sketch as the authoring tool's default entry description, so the
    curator does not invent one."""

    text = cif_bytes.decode("utf-8", "replace")
    name = (_cif_value(text, "_chemical_name_common",
                       "_chemical_name_mineral",
                       "_chemical_name_systematic")
            or _cif_value(text, "_chemical_formula_sum"))
    hm = _cif_value(text, "_symmetry_space_group_name_H-M",
                    "_space_group_name_H-M_alt").split(":")[0].strip()
    it = _cif_value(text, "_space_group_IT_number",
                    "_symmetry_Int_Tables_number")
    year = _cif_value(text, "_journal_year")
    space_group = f"{hm} ({it})" if hm and it else hm
    return ", ".join(p for p in (name, space_group, year) if p)


def pin(tokens, session_rows=None):
    """Resolve the chosen rows to their current revisions, names, and
    discovery hints.

    Returns a list of ``{id, revision, reference_id, elements,
    description}`` dicts.  Only the chosen few are fetched -- to read
    each one's revision and derive, from the CIF metadata, its
    reference_id, its composition, and a human description; a broad
    search is never downloaded wholesale."""

    pinned = []
    for cod_id in resolve_ids(tokens, session_rows):
        data = fetch_cif(cod_id)
        pinned.append({
            "id": cod_id,
            "revision": cif_revision(data),
            "reference_id": _auto_reference_id(data),
            "elements": _composition(data),
            "description": _source_description(data)})
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


def _named_pins(pinned):
    """Yield ``(reference_id, entry)`` for each pinned structure.

    The ``reference_id`` is the metadata-derived name from
    :func:`_auto_reference_id`, falling back to ``cod_<id>`` when the
    CIF lacked the needed fields; should two structures reduce to the
    same name, the later ones get a trailing counter so every
    reference_id stays unique (manifest rule 5).  Both the complete
    manifest and the sketch name their solids through here, so the two
    outputs agree on every reference_id."""

    used = {}
    for entry in pinned:
        base = entry.get("reference_id") or f"cod_{entry['id']}"
        used[base] = used.get(base, 0) + 1
        name = base if used[base] == 1 else f"{base}_{used[base]}"
        yield name, entry


def _complete_manifest(pinned):
    """Render pinned entries as a complete, runnable manifest.

    This is the default ``pin`` output: a finished schema-v2 manifest
    the producer runs directly, with no hand-editing required.  It
    carries the shared ``[characterization]`` recipe and ``[defaults]``
    run settings -- both taken from :mod:`curation_manifest`, so they
    match what expand_manifest emits and what the producer resolves
    (DESIGN 5.7) -- plus one ``[[reference_solid]]`` per pinned
    structure.

    No ``[[reference_solid.entry]]`` customizations are written: the
    harvest auto-discovers one representative per distinct environment
    and falls back to each element's ``"isolated"`` baseline for its
    default (DESIGN 5.2.3), so a curator adds an entry only to override
    that.  The CIF-derived ``source_description`` rides along on each
    solid, and the whole manifest round-trips through
    :func:`curation_manifest.load_manifest_v2`."""

    solids = []
    for name, entry in _named_pins(pinned):
        solids.append(ReferenceSolid(
            reference_id=name,
            system_type="crystalline",
            # cod_id is a schema integer; pin() may carry it as text.
            cod_id=int(entry["id"]),
            cod_revision=entry["revision"],
            source_description=entry.get("description") or None))
    manifest = CurationManifest(
        schema_version=2, manifest_path="",
        characterization=default_characterization(),
        defaults=default_run_settings(),
        reference_solids=solids)
    return format_manifest(manifest)


def _manifest_fragment(pinned):
    """Render pinned entries as a sketch manifest (the --sketch-only
    output).

    Unlike the complete manifest, this is a minimal starting point for
    the authoring tool (expand_manifest): a ``schema_version`` header
    plus one ``[[reference_solid]]`` stub per pinned structure, each
    carrying the CIF-derived discovery hints (``elements`` and a
    ``source_description``) and a reminder of the run settings expand
    will stamp in.  ``cod_fish.py pin <ids> --sketch-only > sketch.toml``
    writes a file expand reads directly.  Solids are named through
    :func:`_named_pins`, so the sketch and the complete manifest agree
    on every reference_id (manifest rule 5)."""

    lines = ["schema_version = 2", ""]
    for name, entry in _named_pins(pinned):
        lines.append("[[reference_solid]]")
        lines.append(f'reference_id = "{name}"')
        lines.append('system_type = "crystalline"')
        lines.append(f'cod_id = {entry["id"]}')
        lines.append(f'cod_revision = "{entry["revision"]}"')
        # Discovery hints the authoring tool reads to auto-fill each
        #   entry's element and description.  They are not part of the
        #   manifest schema -- the producer ignores them, and the
        #   finished manifest the authoring tool writes omits them.
        if entry.get("elements"):
            joined = ", ".join(f'"{symbol}"'
                               for symbol in entry["elements"])
            lines.append(f"elements = [{joined}]")
        if entry.get("description"):
            escaped = entry["description"].replace(
                "\\", "\\\\").replace('"', '\\"')
            lines.append(f'source_description = "{escaped}"')
        lines.append("# fill in: kpoint_spec, scf_threshold, basis, "
                     "functional, kpoint_integration, and entries")
        lines.append("")
    return "\n".join(lines)


def _run_get(args):
    """Run the ``get`` verb: fetch one CIF by COD id."""

    if len(args.targets) != 1:
        raise CodFishError(
            "get needs exactly one COD id, e.g. "
            "`cod_fish.py get 9008463`")
    cod_id = args.targets[0]
    data = fetch_cif(cod_id, revision=args.revision)
    out = args.output or f"{cod_id}.cif"
    with open(out, "wb") as handle:
        handle.write(data)
    print(f"Wrote {out} (COD {cod_id}, revision "
          f"{cif_revision(data)}).")
    return 0


def _run_search(args, with_rank):
    """Run the ``search`` (or ``rank``) verb: query by composition,
    print the numbered table, and save the session."""

    if not args.elements:
        raise CodFishError(
            "search/rank needs --elements, e.g. "
            "`cod_fish.py search --elements Si`")
    rows = search(args.elements, text=args.text,
                  max_extra=args.max_extra, limit=args.limit,
                  fuzzy=args.fuzzy)
    if not rows:
        print("No matching structures.")
        return 0
    reasons_by_id = {}
    if with_rank:
        ordered = rank(rows)
        rows = [candidate for candidate, _ in ordered]
        reasons_by_id = {c["id"]: r for c, r in ordered}
    _print_table(rows)
    if with_rank:
        print("\nTriage (advisory):")
        for index, candidate in enumerate(rows, start=1):
            print(f"  [{index}] {candidate['id']}  "
                  f"{' · '.join(reasons_by_id[candidate['id']])}")
    save_session(rows, args.elements)
    print(f"\nSaved {len(rows)} rows to {SESSION_FILE}; pin them by "
          f"index, e.g. `cod_fish.py pin 1 3`.")
    return 0


def _run_pin(args):
    """Run the ``pin`` verb: resolve chosen rows to pinned revisions
    and print a manifest.  By default this is a complete, runnable
    manifest; ``--sketch-only`` prints the sketch for expand_manifest
    to finish instead."""

    if not args.targets:
        raise CodFishError(
            "pin needs at least one row index or COD id, e.g. "
            "`cod_fish.py pin 1 3`")
    pinned = pin(args.targets, session_rows=load_session())
    if args.sketch_only:
        print(_manifest_fragment(pinned))
    else:
        print(_complete_manifest(pinned))
    return 0


def _build_parser():
    """Build the cod_fish argument parser.

    A single flat parser (no sub-parsers), so one ``--help`` lists
    every option together with its default -- the convention the
    other Imago scripts follow.  The first positional selects the
    verb; each option's help notes which verb it applies to."""

    description_text = (
        "Fetch and search crystal structures from the "
        "Crystallography Open Database (COD).\n\n"
        "Verbs (the first argument):\n"
        "  get     fetch one structure's CIF by COD id (strict, "
        "revision-verified)\n"
        "  search  find candidates by composition; print a numbered "
        "table and save\n"
        "          a session so rows can later be picked by index\n"
        "  rank    like search, and also print advisory triage of "
        "the candidates\n"
        "  pin     resolve chosen rows (by index or id) to pinned "
        "revisions and\n"
        "          print a complete manifest (--sketch-only for a "
        "sketch to expand)")
    epilog_text = (
        "Examples:\n"
        "  cod_fish.py search --elements Si\n"
        "  cod_fish.py rank --elements Si O --max-extra 1\n"
        "  cod_fish.py pin 1 3 5\n"
        "  cod_fish.py get 9008463 --revision 291735 -o au.cif\n")

    parser = argparse.ArgumentParser(
        prog="cod_fish.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=description_text, epilog=epilog_text)

    # The verb to run.
    parser.add_argument(
        "command", choices=["get", "search", "rank", "pin"],
        help="Which operation to run (see the list above). "
             "Required.")
    # Positional targets, interpreted per verb.
    parser.add_argument(
        "targets", nargs="*", metavar="ID|INDEX",
        help="get: one COD id.  pin: row indices from the last "
             "search (ranges like 1-3 allowed) or raw COD ids.  "
             "Unused for search/rank.  Default: none.")
    # Composition filter (search/rank).
    parser.add_argument(
        "-e", "--elements", nargs="+", metavar="EL", default=None,
        help="search/rank: element symbols that must be present, each "
             "with an optional atom count on either side (Fe, 2Fe, or "
             "Fe2).  The search is for exactly this element set unless "
             "--max-extra is given; counts select an exact formula "
             "unit (e.g. 2Fe 3O -> Fe2O3) unless --fuzzy.  Default: "
             "none.")
    # Fuzzy count matching (search/rank).
    parser.add_argument(
        "--fuzzy", action="store_true", default=False,
        help="search/rank: apply element counts via COD's own "
             "(broader, digit-prefix) formula match instead of the "
             "exact client-side formula-unit reduction.  Default: off "
             "(exact).")
    # Extra-element slack (search/rank).
    parser.add_argument(
        "--max-extra", type=int, default=0, metavar="K",
        help="search/rank: allow up to K elements beyond those in "
             "--elements.  Default: 0 (exact composition).")
    # Free-text filter (search/rank).
    parser.add_argument(
        "--text", default=None, metavar="STR",
        help="search/rank: COD free-text filter, e.g. an author "
             "name.  Default: none.")
    # Result cap (search/rank).
    parser.add_argument(
        "--limit", type=int, default=None, metavar="N",
        help="search/rank: show at most N candidates.  Default: all "
             "matches.")
    # Advisory triage on search (the `rank` verb is the same thing).
    parser.add_argument(
        "--rank", action="store_true", default=False,
        help="search: also print advisory triage of the candidates "
             "(equivalent to the `rank` verb).  Default: off.")
    # Pinned-revision check (get).
    parser.add_argument(
        "-r", "--revision", default=None, metavar="REV",
        help="get: pinned COD revision to verify against the served "
             "CIF; a mismatch is a hard error.  Default: no revision "
             "check.")
    # Output path (get).
    parser.add_argument(
        "-o", "--output", default=None, metavar="PATH",
        help="get: output CIF path.  Default: <cod_id>.cif in the "
             "current directory.")
    # Sketch-vs-complete manifest (pin).
    parser.add_argument(
        "--sketch-only", action="store_true", default=False,
        help="pin: print a minimal sketch (structure stubs plus "
             "element and description hints) for expand_manifest to "
             "finish, instead of the default complete, runnable "
             "manifest.  Default: off (complete manifest).")
    return parser


def main(argv=None):
    """Entry point for the ``cod_fish.py`` command line."""

    args = _build_parser().parse_args(argv)
    try:
        if args.command == "get":
            return _run_get(args)
        if args.command == "pin":
            return _run_pin(args)
        # search and rank share one path; the `rank` verb and the
        #   `--rank` flag both turn the triage on.
        with_rank = args.rank or args.command == "rank"
        return _run_search(args, with_rank=with_rank)
    except CodFishError as error:
        print(f"cod_fish: {error}", file=sys.stderr)
        return 1


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
