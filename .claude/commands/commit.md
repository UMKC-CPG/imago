---
allowed-tools: Bash
description: Minimal commit with a brief message. No docs or quality checks.
---

# Quick Commit

Perform a minimal commit with no ceremony. Follow these steps:

1. Run `git status` to show what has changed.
2. Run `git diff --staged` to show what is already staged (if any).
3. Ask the programmer which files to stage (or confirm staged are correct).
4. Ask for a one-line commit description if not yet provided.
5. Stage the specified files with `git add`.
6. Commit with a message that follows these rules:
   - First line: imperative mood, under 72 characters (e.g., "Fix overlap
     sign convention in gaussian.py")
   - Blank line
   - Include this trailer on the final line:
     `Co-Authored-By: Claude <noreply@anthropic.com>`
7. Report the resulting commit hash.

Do NOT update documentation, run tests, run quality checks, or take any action
beyond staging and committing.
