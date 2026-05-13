---
allowed-tools: Read, Glob, Grep
description: Focus the session: summarize TODO.md and load relevant context.
---

# Session Focus

You are starting a focused development session. Follow these steps precisely.

**Step 1: Read TODO.md**

Read `dev/TODO.md`. Identify all unchecked items (`- [ ]`) in each section: VISION,
ARCHITECTURE, DESIGN, PSEUDOCODE, CODE.

**Step 2: Summarize**

Present the pending items in this format:

```
Pending items by level:

VISION (N):
  - <item>

ARCHITECTURE (N):
  - <item>

DESIGN (N):
  - <item>

PSEUDOCODE (N):
  - <item>

CODE (N):
  - <item>
```

If a section has no pending items, write "(none)".

**Step 3: Ask**

Ask the programmer: "What would you like to work on today?"

**Step 4: Load focused context**

Based on the programmer's answer, read only what is needed:

- VISION work: read `dev/VISION.md`
- ARCHITECTURE work: read `dev/ARCHITECTURE.md`
- DESIGN work: grep `dev/DESIGN.md` for relevant section headers, then read
  only those sections (do not read the entire file)
- PSEUDOCODE work: read `dev/PSEUDOCODE.md`
- CODE work: read the relevant source files in `src/` and their corresponding
  sections in dev/DESIGN.md

**Step 5: Confirm**

Briefly state what you have loaded and confirm you are ready to begin.
