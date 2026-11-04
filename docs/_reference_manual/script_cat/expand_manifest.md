---
layout: default
title: expand_manifest.py
script: expand_manifest.py
description: 'Authoring tool that expands a structure sketch (from cod_fish.py --sketch-only or written by hand) into a complete curation manifest for build_initial_potentials.py, filling in the database-wide fingerprint recipe and the shared basis, functional, k-point integration and SCF threshold, either mechanically or through interactive prompts.'
parent: Script Catalog
---

# expand_manifest

The authoring step between [cod_fish.py](cod_fish.html) and
[build_initial_potentials.py](build_initial_potentials.html).
`cod_fish.py pin --sketch-only` prints one `[[reference_solid]]`
stub per structure carrying only its identity (`reference_id`,
`system_type`) and source (`cod_id` plus `cod_revision`, or a local
`structure_path`), plus optional discovery hints. This tool reads
that sketch, from a file or from standard input as a pipe, and fills
in the rest of the manifest the producer needs:

* the database-wide fingerprint recipe, the `[characterization]`
  block;
* the shared run settings (basis, functional, k-point integration,
  SCF threshold) as a single top-level `[defaults]` block that every
  solid inherits;
* any per-structure customizations, which are overrides on the
  environments the harvest auto-discovers.

For a straightforward set of crystals `cod_fish` already emits a
complete, runnable manifest on its own. This tool exists for the
cases that need more: local, non-COD structures; a non-default recipe
or run settings; or interactively authored per-structure
customizations. All three tools share one schema definition by
importing `curation_manifest.py`, so nothing drifts between them.

**The k-point density is deliberately left out** of the emitted
`[defaults]`. The producer predicts a starting density and verifies
it by a convergence sweep, so pinning a density here would only
override that.

## Usage

**Commandline:**

```
expand_manifest.py [SKETCH] [-o PATH] [-i]
                   [--basis {mb,fb,eb}]
                   [--functional {wigner,ceperley-alder,hedin-lundqvist,pbe}]
                   [--kpoint-integration {gaussian,linear-tetrahedral,gaussian-<sigma>}]
                   [--scf-threshold T] [--system-type TYPE] [--help]
```

**Examples:**

```
# Pipe straight from cod_fish and review on screen
cod_fish.py pin 1 3 5 --sketch-only | expand_manifest.py

# Mechanical mode to a file, with a non-default functional
expand_manifest.py sketch.toml --functional pbe -o manifest.toml

# Interactive mode: prompts for settings and per-structure overrides
expand_manifest.py sketch.toml -i -o manifest.toml
```

**Variable Table:**

| Variable   | Accepted Values     | Description |
|:----------:|:-------------------:|:------------|
|SKETCH      |String; file path    | A schema-version-2 sketch manifest with one structure source per solid. |
|PATH        |String; file path    | Where to write the finished manifest. |
|T           |float                | SCF convergence threshold. |
|TYPE        |crystalline, or another registered system type | System type stamped on any structure whose sketch did not name one. |

**Option Table:**

| Option                 | Defaults                 | Description |
|:----------------------:|:------------------------:|:------------|
|SKETCH                  |Standard input            | The sketch to expand. Omitting it lets `cod_fish.py pin --sketch-only | expand_manifest.py` work as a pipe. |
|-o / --output           |Standard output           | Write the manifest to PATH. Required with `--interactive`, whose prompts use the terminal. |
|-i / --interactive      |Off (mechanical mode)     | Walk through the per-structure curation with prompts instead of emitting a fill-in template. |
|--basis                 |fb                        | Shared basis for every structure: `mb`, `fb` or `eb`. |
|--functional            |wigner                    | Shared exchange-correlation functional. |
|--kpoint-integration    |gaussian                  | Shared k-point integration: `gaussian`, `linear-tetrahedral`, or `gaussian-<sigma>` to name a thermal smearing width in eV. |
|--scf-threshold         |1e-6                      | Shared SCF convergence threshold. |
|--system-type           |crystalline               | System type for structures whose sketch did not name one. |
|--help                  |N/A                       | Print usage information and exit. |

## The Two Modes

**Mechanical** (the default) sets the shared settings once in a
`[defaults]` block and the database-wide recipe once, carries each
sketched structure through as a sparse solid, and emits the result
with a commented copy-and-edit template for the optional
per-structure customizations, which you fill in by hand. Without
`-o` it prints to standard output, so the result can be reviewed,
redirected, or piped straight from `cod_fish`.

**Interactive** (`-i`) walks each sketched structure and prompts for
the shared settings once, then for the optional per-structure
customizations, offering a sensible default at every step. The
fingerprint recipe is set automatically, so it writes a complete,
immediately loadable manifest. Because the prompts use standard
input and output, this mode needs a sketch file argument and an
output file.
