---
allowed-tools: Read, Glob, Grep, Edit, Bash
description: >
  Audit source files for line-length violations and
  reflow to the 70-80 character target band.
---

# Lint: Line-Length Audit and Reflow

Audit source code files for line-length violations and
reflow them into the 70-80 character target band. The user
provides a target as `$ARGUMENTS`, which may be:

- A specific file (e.g., `src/imago/dos.F90`)
- A directory (e.g., `src/imago/`)
- A glob pattern (e.g., `**/*.f90`, `src/scripts/*.py`)
- Omitted — default to `src/`

## Helper scripts

Two helper scripts live in `.claude/commands/scripts/`:

- **rewrap_prose.py** — reflows comment-block paragraphs
  and Python docstring paragraphs to fill the 70-80 band.
  Supports Fortran `!` comments and Python `#` comments
  plus prose paragraphs inside multi-line triple-quoted
  docstrings.  Run with `--check` for a dry run.
  Usage: `python3 .claude/commands/scripts/rewrap_prose.py
  FILE [--check] [--lines S-E]`

  Python notes:
  - `## code` lines are treated as commented-out code and
    are never touched by the reflow engine.
  - A triple-quoted string is treated as a docstring only
    when its opening quote is the first non-whitespace
    token on its line (detection uses `tokenize`, so
    assignments such as `help_text = """..."""` are
    correctly left alone — argparse text remains a
    manual exercise).
  - Structural content inside docstrings (bullet lists,
    numbered lists, doctest blocks, section headers,
    RST directives, numpy-style `name : type` field
    lines, and code/data blocks with no sentence
    punctuation) is detected and left untouched.
  - Double-space-after-period typography is preserved
    through a reflow pass.

- **rewrap_code.py** — auto-detects over-wrapped code
  blocks (2+ lines that could fit on fewer) and compacts
  them.  Tries single-line join first, then comma-split.
  Scan mode is the default:
  `python3 .claude/commands/scripts/rewrap_code.py
  SOURCE [--check]`

  For manual control, use directive mode:
  `python3 .claude/commands/scripts/rewrap_code.py
  SOURCE --noscan DIRECTIVES [--check]`

  Directive format (one per line, # comments allowed):
    `unwrap START END` — join lines into one, strip
      gratuitous continuation parens.
    `rewrap START END` — join then re-split at commas.
    `merge-strings START END` — concatenate adjacent
      string literal contents, re-split to fill 70-80
      per line.

Follow these steps precisely.

**Step 1: Identify target files**

Resolve `$ARGUMENTS` into a list of source files. Use
Glob to expand directories and patterns. Include only
source files (`.f90`, `.F90`, `.py`, `.pl`, `.pm`). Skip
binary files, swap files, and build artifacts.

If the target resolves to more than 20 files, list the
count and ask the programmer whether to proceed or narrow
the scope.

**Step 2: Scan for violations**

Read each target file and identify lines that violate the
line-length rules:

1. **Over-length** (hard violation): lines exceeding 80
   characters. These MUST be fixed.
2. **Under-filled** (soft violation): lines shorter than
   70 characters where the content could reasonably be
   reflowed to fill closer to 70-80. This applies to:
   - Comment blocks where consecutive short comment lines
     could be merged into fewer, fuller lines.
   - String literals or argument lists that were broken
     early when more content would fit on the line.
   - Code lines split across 3-4 lines when they would
     fit on 1-2 lines (over-wrapped code).

   Do NOT flag as under-filled:
   - Lines that are naturally short (e.g., `implicit
     none`, `end if`, `return`, blank lines, single
     closing brackets/parentheses).
   - Lines where filling further would break code logic
     or readability (e.g., one-item-per-line formatting
     that aids clarity).
   - The last line of a reflowed paragraph or comment
     block (it may be a short remainder and that is
     fine).

For Fortran (`.f90`, `.F90`) and Python (`.py`) files,
also run `rewrap_prose.py FILE --check` to detect comment
blocks and (for Python) docstring paragraphs that need
reflow.  Include its output in the report.

**Step 3: Report findings**

Present a summary grouped by file:

```
Line-length audit: N file(s) scanned, M violation(s).

--- src/imago/dos.F90 ---
  Over-length (K):
    L42  (87 chars): <truncated line preview>
    L108 (83 chars): <truncated line preview>
  Under-filled (J):
    L55-L58: comment block averages 41 chars/line
    L200 (38 chars): argument list could join prev line
  Comment reflow (rewrap_prose --check):
    L55-L58: 4 -> 3 lines

--- src/scripts/condense.py ---
  Over-wrapped code (K):
    L236-239: sys.exit() split across 4 lines (fits 1)
    L301-304: merge-strings candidate (3 fragments)
```

If no violations are found in any file, say so and stop.

**Step 4: Ask for approval**

Ask the programmer:
"Shall I reflow these files to fix the violations? You
can approve all, select specific files, or skip."

Do NOT proceed without explicit approval.

**Step 5: Apply fixes (only after approval)**

Use the helper scripts first, then fall back to Edit for
anything they cannot handle.

**5a. Comment and docstring prose** — run rewrap_prose.py
on Fortran (`.f90`, `.F90`) and Python (`.py`) files:
```
python3 .claude/commands/scripts/rewrap_prose.py FILE
```
This automatically finds and reflows comment paragraphs
(both Fortran `!` and Python `#`) plus prose paragraphs
inside Python docstrings.  For Python, double-space-
after-period typography is preserved across the reflow.

**5b. Over-wrapped code lines** — run rewrap_code.py:
```
python3 .claude/commands/scripts/rewrap_code.py \
    FILE --check
```
Review the output, then apply (omit `--check`):
```
python3 .claude/commands/scripts/rewrap_code.py FILE
```
The script automatically finds over-wrapped blocks and
compacts them.  Run it repeatedly until it reports
0 changes (it is idempotent).

For blocks the scanner cannot handle automatically
(e.g., merge-strings), use directive mode:
```
python3 .claude/commands/scripts/rewrap_code.py \
    FILE --noscan /tmp/directives.txt
```

**5c. Remaining violations** — use Edit for anything the
scripts cannot handle:
- Over-length code lines that need intelligent wrapping
  at language-specific break points (Fortran `&`,
  Python implicit continuation, etc.)
- Under-filled lines that need manual judgment

Preserve meaning: never alter logic, variable names, or
executable behavior. Only whitespace and line breaks
should change.

**Step 6: Verify and summarize**

Re-scan each modified file to confirm no new violations
were introduced. Report what was changed:

```
Reflow complete:
  src/imago/dos.F90:
    rewrap_prose: 3 comment blocks reflowed
    rewrap_code: 5 directives applied (5 unwrap)
    manual edit: 2 over-length lines wrapped
  src/scripts/condense.py:
    rewrap_code: 8 directives applied (6 unwrap,
      2 merge-strings)
```
