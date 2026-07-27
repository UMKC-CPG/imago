# Contributing to Imago

Imago is developed by the Computational Physics Group (CPG) at the
University of Missouri-Kansas City. This note covers the two things
a contributor needs to get right before writing code: how licensing
and credit work, and how the design document chain works.

## Licensing and Credit

Imago is licensed under the Educational Community License, Version
2.0 (ECL-2.0). By contributing, you agree that your contribution is
licensed under those same terms. This follows automatically from
section 5 of the license, but it is stated here so that nobody has
to go read the license to find it out.

**Copyright and authorship are two different things, and they go in
two different places.** Conflating them is the most common mistake,
and it is an easy one to make: claiming credit for your own work is
correct and encouraged, but the copyright line is not where that
claim belongs.

### Copyright: uniform, one holder

Every source file carries the same two-line header, and nothing
about it varies from file to file:

```fortran
!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis
```

For Python, the same two lines with `##` in place of `!!`, placed
after the `#!/usr/bin/env python3` shebang and before the module
docstring. The doubled comment markers (`!!` and `##`) mark these
lines as structured content so that `rewrap_prose.py` will not
reflow them.

Do **not** add a personal copyright block to a file you write or
modify. It does not increase your credit, and a file carrying a
copyright holder different from the rest of the tree has to be
tracked separately forever after.

### Authorship: per file, and in the citation metadata

Credit for who wrote what belongs in two places:

1. **An `Author:` line in the file's own header block or module
   docstring.** This is the right home for "I designed and wrote
   this." Name the people who did the work.
2. **The `authors:` list in `CITATION.cff`.** This is what actually
   propagates into other people's reference lists when they cite
   Imago, so it is the credit that reaches the published record.
   If you have made a substantial contribution, add yourself, with
   your ORCID if you have one.

If you have contributed substantially and are not in `CITATION.cff`,
that is an oversight worth raising rather than working around.

### Third-party material

If you bring in code, data, or an algorithm from an outside source,
two things must happen. Add an entry to the `NOTICE` file at the
repository root describing what was used and where it came from --
section 4(d) of the license requires that `NOTICE` travel with every
redistribution, so this is the mechanism by which upstream credit
survives downstream. And cite the source in the file itself, next to
the code that uses it, so a reader encountering the algorithm can
find the paper or book it came from without hunting.

Attribution is not a formality in scientific software. It is the
same obligation that governs a methods section.

## The Design Document Chain and Coding Style

Imago follows a five-level document chain -- vision, architecture,
design, pseudocode, then code -- and a new feature enters at the top
and flows down, each level written before the one below it exists.
The rule with the most bite: **before editing any file under `src/`,
the governing PSEUDOCODE section must already exist and must already
describe the change.** Writing the code first and back-filling the
specification defeats the purpose, because nobody reviewed the code
against it.

The chain itself lives in `dev/`, starting with `dev/VISION.md`, and
`dev/TODO.md` tracks tasks by level.

**`CLAUDE.md` at the repository root is the authoritative statement
of both the chain discipline and the coding style**, including the
80-character line limit and the documentation and naming rules. It
is written for coding agents, but it governs human contributors
equally. Read it before your first change; this file does not repeat
its rules, so where the two ever appear to disagree, `CLAUDE.md`
wins.
