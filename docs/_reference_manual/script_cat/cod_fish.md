---
layout: default
title: cod_fish.py
script: cod_fish.py
description: 'Searches the Crystallography Open Database (COD) by composition and downloads CIF files by COD id. Four verbs separate browsing (search, rank) from reproducible acquisition (pin, get), so a pinned id and revision can always be re-fetched exactly.'
parent: Script Catalog
---

# cod_fish

An acquisition front-end for the Crystallography Open Database (COD).
COD holds hundreds of thousands of experimentally determined crystal
structures, served as CIF files over plain HTTP. This script is a thin
front-end over COD's own endpoints: it keeps no local database and
writes no queries of its own. It offers four verbs that separate
*browsing* from *reproducible acquisition*.

| Verb     | Purpose |
|:--------:|:--------|
|`search`  | Find candidate structures by composition (and optionally author or free text). Prints a numbered table and saves the result set so later verbs can refer to rows by index. |
|`rank`    | Like `search`, and also prints advisory triage: the candidates are annotated and ordered by interpretable signals so you can spot the likely "real" phase. It narrows the field; it does not decide. |
|`pin`     | Resolve chosen rows (by index, so no eight-digit typing) to their current COD revisions and print a complete, runnable manifest. This is the bridge from browsing to a pinned pull. |
|`get`     | Fetch one structure's CIF by its COD id, optionally verifying a pinned revision. This is the reproducible path. |

**Why the split matters.** COD entries are revised over time. A
calculation is reproducible only if the structure it used can be
fetched again exactly, so discovery (`search`, `rank`, and the saved
session) never feeds a build directly. A manifest only ever pins a
`cod_id` together with a `cod_revision`, and `get` fetches that pair
strictly, refusing any other revision.

The CIF that `get` writes is converted to a skeleton with
[cif2skl.py](cif2skl.html).

## Usage

**Commandline:**

```
cod_fish.py {get,search,rank,pin} [ID|INDEX ...]
            [-e EL [EL ...]] [--fuzzy] [--max-extra K] [--text STR]
            [--limit N] [--rank]
            [-r REV] [-o PATH]
            [--sketch-only] [--save-cif [DIR]]
            [--help]
```

One flat option list serves all four verbs, so `--help` lists every
option with its default in one place. Each option's description notes
which verb it applies to.

**Examples:**

```
cod_fish.py search --elements Si                 # single-element Si
cod_fish.py search --elements Si O --max-extra 1
cod_fish.py rank   --elements 2Fe 3O             # Fe2O3, ranked
cod_fish.py pin 2 5 7                            # by row index
cod_fish.py pin 2 5 7 --save-cif cifs/ > manifest.toml
cod_fish.py get 9008463 --revision 291735 -o au.cif
```

**Variable Table:**

| Variable   | Accepted Values                     | Description |
|:----------:|:-----------------------------------:|:------------|
|ID          |integer; COD id                      | A COD structure id, as printed in the `cod_id` column of a search table or shown on the COD website. |
|INDEX       |integer or range such as `1-3`       | A row number from the most recent `search` or `rank` table. |
|EL          |Element symbol, optionally with a count: `Fe`, `2Fe`, or `Fe2` | An element that must be present. A count selects an exact formula unit. |
|K           |integer                              | Number of elements allowed beyond those listed with `--elements`. |
|STR         |String                               | COD free-text filter, for example an author name. |
|N           |integer                              | Maximum number of candidates to show. |
|REV         |COD revision number                  | The revision stamped in the CIF header as `$Revision: NNN $`. |
|PATH        |String; file path                    | Output CIF path. |
|DIR         |String; directory path               | Directory in which to keep the CIFs that `pin` fetched. |

**Option Table:**

| Option           | Verbs        | Outputs           | Defaults                     | Description |
|:----------------:|:------------:|:-----------------:|:----------------------------:|:------------|
|-e / --elements   |search, rank  |N/A                |None                          | Element symbols that must be present, each with an optional atom count on either side. The search is for exactly this element set unless `--max-extra` is given; counts select an exact formula unit (e.g. `2Fe 3O` means Fe2O3) unless `--fuzzy`. |
|--fuzzy           |search, rank  |N/A                |Off (exact)                   | Apply element counts via COD's own broader, digit-prefix formula match instead of the exact client-side formula-unit reduction. |
|--max-extra       |search, rank  |N/A                |0 (exact composition)         | Allow up to K elements beyond those in `--elements`. |
|--text            |search, rank  |N/A                |None                          | COD free-text filter, e.g. an author name. |
|--limit           |search, rank  |N/A                |All matches                   | Show at most N candidates. |
|--rank            |search        |N/A                |Off                           | Also print advisory triage of the candidates. Equivalent to using the `rank` verb. |
|-r / --revision   |get           |N/A                |No revision check             | Pinned COD revision to verify against the served CIF. A mismatch is a hard error. |
|-o / --output     |get           |CIF file           |`<cod_id>.cif` in the current directory | Output CIF path. |
|--sketch-only     |pin           |Sketch manifest on stdout | Off (complete manifest) | Print a minimal sketch (structure stubs plus element and description hints) for `expand_manifest.py` to finish, instead of a complete, runnable manifest. |
|--save-cif        |pin           |CIF files in DIR   |Off; with the flag and no DIR, the current directory | Also write each pinned structure's CIF into DIR, named `<reference_id>.cif` to match its manifest block. `pin` downloads these anyway, so this costs no extra fetch. DIR is created if absent, and an existing file of the same name is overwritten. Each path written is reported on standard error so it does not land in a redirected manifest. |
|--help            |all           |N/A                |N/A                           | Print usage information and exit. |

## Verb Details

### search

> cod_fish.py search --elements EL [EL ...] [--max-extra K] [--fuzzy] [--text STR] [--limit N]

Queries COD for structures matching the composition and prints a
numbered table with these columns:

```
  #    cod_id  formula          spacegroup   name                   a       b       c
```

The rows are then saved to a session file named
`cod_fish_results.json` in the current directory. The `pin` verb reads
this file so that you can refer to structures by their row number
instead of retyping COD ids. A new `search` overwrites the session.

Composition matching is exact by default: `--elements Si O` finds
structures containing silicon and oxygen and nothing else. Use
`--max-extra` to allow additional elements, and attach counts
(`2Fe 3O`) to restrict to one formula unit.

### rank

> cod_fish.py rank --elements EL [EL ...] ...

Takes the same options as `search`, prints the same table, and then
orders the candidates by a small set of interpretable signals drawn
only from COD metadata. The reasons are printed beside every row so
you can see why it landed where it did:

| Signal                 | Effect  | Meaning |
|:-----------------------|:-------:|:--------|
|`named:<mineral>`       | +2      | COD records a mineral or common phase name. |
|`ambient`               | +2      | Measured at ambient temperature and pressure. |
|`non-ambient`           | -2      | High-pressure, high-temperature or cryogenic measurement. |
|`sg<N>(consensus)`      | +1      | Shares the space group that the plurality of candidates have. |
|`volume-outlier`        | -2      | Cell volume is more than twice, or less than half, the median. |

This is triage to narrow "which of these is the real phase". It is
not a validation verdict. When a composition has many entries, read
the reasons and choose.

### pin

> cod_fish.py pin INDEX|ID [INDEX|ID ...] [--sketch-only] [--save-cif [DIR]]

Resolves each chosen row (a row number from the last `search`, a
range like `1-3`, or a raw COD id) to its current COD revision by
fetching the CIF, then prints a manifest on standard output.
Redirect it to a file:

```
cod_fish.py pin 1 3 5 > manifest.toml
```

By default the manifest is complete and runnable: a TOML document
with one `[[reference_solid]]` block per structure, each carrying its
`reference_id`, `cod_id` and `cod_revision`, plus default run
settings. With `--sketch-only` only a minimal sketch is printed,
containing the structure stubs and discovery hints (`elements` and a
`source_description`), for `expand_manifest.py` to finish. Both forms
name the solids identically, so a sketch and its completed manifest
always agree on every `reference_id`.

Manifests are the input to the curated reference-structure workflow
(`curation_manifest.py`, `expand_manifest.py` and
`build_initial_potentials.py`), where each pinned structure is
fetched with `get` and then run through Imago.

### get

> cod_fish.py get ID [-r REV] [-o PATH]

Fetches one CIF by COD id and writes it to PATH. With `--revision`
the revision stamped in the served CIF's header is compared to REV,
and a mismatch is a hard error rather than a warning. This is the
strict path a manifest relies on: if COD has revised the entry since
it was pinned, the fetch fails loudly instead of silently returning a
different structure.

## Session File

`search` and `rank` write `cod_fish_results.json` in the current
directory. It maps row numbers to COD ids and is read by `pin`. If no
session file exists, every argument to `pin` must be a raw COD id.
