---
allowed-tools: Read, Glob, Grep
description: Check consistency across the document chain and flag drift.
---

# Refine: Design Chain Consistency Check

Walk the VISION → ARCHITECTURE → DESIGN → PSEUDOCODE → Code chain,
checking that each level faithfully implements the one above it and that
no level has drifted out of sync. Follow these steps precisely.

**Step 1: Read the chain**

Read (or skim headings of) the following files in order:
- `dev/VISION.md`
- `dev/ARCHITECTURE.md`
- `dev/DESIGN.md` (section headings first; read specific sections as needed)
- `dev/PSEUDOCODE.md`
- Key source files in `src/` (guided by what DESIGN.md and PSEUDOCODE.md
  reference)

**Step 2: Check each boundary**

Work through the chain top-down, checking one boundary at a time. At each
boundary, ask: "Does the lower level faithfully implement what the upper
level claims?"

```
VISION → ARCHITECTURE
  Does the architecture serve the stated goals and principles? Are there
  architectural choices that contradict a vision principle?

ARCHITECTURE → DESIGN
  Does the design reference the correct modules and layouts? Are there
  design sections that describe structures not in the architecture, or
  architectural components with no design coverage?

DESIGN → PSEUDOCODE
  Does the pseudocode match the algorithms described in the design? Are
  there design decisions that the pseudocode silently ignores or
  contradicts?

PSEUDOCODE → Code
  Does the source code implement the pseudocode? Are there code paths
  with no pseudocode coverage, or pseudocode steps not reflected in the
  implementation?
```

Also check upward: does the code reveal something that should be
documented in a higher-level file but is not?

**Step 3: Check TODO.md alignment**

Read `dev/TODO.md`. Verify:
- Every unchecked item cites the document level it belongs to (VISION /
  ARCHITECTURE / DESIGN / PSEUDOCODE / CODE).
- No completed work is still listed as pending.
- No pending work is missing from the list.

**Step 4: Report findings**

Present each inconsistency as a numbered item in this format:

```
1. [BOUNDARY] DESIGN → PSEUDOCODE
   [FILE] PSEUDOCODE.md, §3.2
   [ISSUE] Design §8.6.3 specifies susceptibility function S = sigma *
   (1 - cos(theta)), but PSEUDOCODE.md omits the sigma term.
   [SUGGESTION] Add the full susceptibility formula to PSEUDOCODE.md §3.2.
```

If no inconsistencies are found at a boundary, state that explicitly:
"VISION → ARCHITECTURE: consistent."

**Step 5: Propose resolution**

For each finding, ask the programmer whether to:
- Fix the lower level (update the doc or code to match the level above),
  or
- Fix the upper level (update the doc to reflect what the code actually
  does), or
- Defer (add as a TODO item at the appropriate level).

Do NOT make changes without the programmer's decision on direction.
