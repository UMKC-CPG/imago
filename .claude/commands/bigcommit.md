---
allowed-tools: Read, Glob, Grep, Bash, Edit
description: Full commit -- doc propagation check, quality pass, commit message.
---

# Full Commit

Perform a thorough, documented commit. Follow each step in order.

**Step 1: Review all changes**

Run `git diff HEAD` to see all uncommitted changes. List the files and
the document level (VISION/ARCHITECTURE/DESIGN/PSEUDOCODE/CODE) for each.

**Step 2: Document propagation check**

For each changed source file or document:
- Check that dev/VISION.md, dev/ARCHITECTURE.md, and dev/DESIGN.md still
  accurately describe the system after all changes.
- Check that any completed TODO items are marked `- [x]` in dev/TODO.md.
- Check that any new design decisions are reflected in the appropriate document.
Report every inconsistency found. If any are found, pause and ask the programmer
to resolve them before proceeding.

**Step 3: Code quality check**

For each changed source file, verify:
- Lines are <= 80 characters, with 80 as the target, not a ceiling to stay well
  under -- lines shorter than necessary harm readability
- No unused imports or dead code remains
- Variable names are long and expressive, not abbreviated
- Non-obvious logic has explanatory comments
Report any violations. Fix minor issues (line length, trivial style) inline.
For substantive issues, stop and ask the programmer.

**Step 4: Confirm**

Present a brief summary of what will be committed and ask the programmer to
confirm.

**Step 5: Commit**

Stage and commit with a comprehensive message:
- First line: brief summary (imperative mood, <= 72 chars)
- Blank line
- Body: what changed, why, and which document sections are affected
- Blank line
- Include this trailer:
  `Co-Authored-By: Claude <noreply@anthropic.com>`

**Step 6: Report**

Report the commit hash and list any remaining open TODO items related to the
committed work.
