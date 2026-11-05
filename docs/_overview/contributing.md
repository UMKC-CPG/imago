---
layout: default
title: Contributing
nav_order: 4
---

# Contributing to Imago

Imago is open-source software developed by the Computational Physics
Group at the University of Missouri-Kansas City and licensed under
the Educational Community License 2.0. Contributions to the code
follow the design-document discipline described in `CONTRIBUTING.md`
at the root of the repository; contributions to this documentation
are welcome from anyone and are described here.

## Contributing to the documentation

There are two ways to contribute:

- **Email.** If you have questions, troubles, or solutions, email
  Dr. Rulis with them. Good questions are added to the
  documentation, and difficulties and their solutions as time
  permits.
- **Pull request.** Fork the repository, make your changes, and
  open a pull request. This is the preferred method.

Things to know before editing:

- The site is built by GitHub Pages with Jekyll from the `docs/`
  directory of the `main` branch, using the
  [Just the Docs](https://just-the-docs.com/) theme. Keep all
  documentation under `docs/`.
- Every page is a Markdown file with a YAML front matter block that
  names its `title`, its `layout` (always `default`), and its
  `nav_order` within its section. Pages of the Script Catalog also
  carry `script` and `description` fields, which the catalog index
  reads to build its list; a new catalog page needs both or it will
  not be listed.
- The site is organised into three collections: `_overview`,
  `_tutorials` and `_reference_manual`, the last holding the
  `script_cat/` and `theory/` sections. Put a new page in the
  collection it belongs to.
- Mathematics is written in LaTeX between `$$` delimiters and
  rendered by MathJax; see the Theory pages for examples.
- Verify your changes by deploying the site from your fork with
  GitHub Pages, and leave it up so that the maintainers can see the
  rendered result.
- In your pull request, list your changes, additions and removals.
- Changes to the established hierarchy are scrutinized heavily.

The best source for a script page is the script itself: every script
carries a full description in its docstring and prints every option
and default with `--help`. A page that disagrees with `--help` is
wrong.

## Contributing to the code

Imago follows a five-level document chain, from vision through
architecture, design and pseudocode to code, and a change enters at
the top and flows down. The rule with the most bite is that no file
under `src/` is edited until the pseudocode governing the change
exists. `CONTRIBUTING.md` and `CLAUDE.md` at the repository root
state the discipline, the licensing and credit rules, and the coding
style in full; read them before a first change.

Thank you for your effort in contributing.
