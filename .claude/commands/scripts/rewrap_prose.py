#!/usr/bin/env python3
"""rewrap_prose.py -- Reflow comment blocks and Python
docstrings to the 70-80 char band.

Algorithm (per project convention):
  1. Start each line with the correct prefix/indent.
  2. Do NOT wrap until you reach the 80-character limit.
  3. Wrap whatever word caused the overflow onto the next
     line, using the running "continuation prefix" as the
     new line's leading whitespace.

Supported languages (detected from file extension):
  .f90, .F90  -- Fortran: `!` comment blocks.
  .py         -- Python: `#` comment blocks plus prose
                 paragraphs inside triple-quoted
                 docstrings.

Language-specific conventions
-----------------------------

Fortran `!` comments:
  `! prose text`  -- prose comment; reflow candidate.
  `!code`         -- commented-out code; left alone.
  Commented-out code is recognised by the absence of a
  space immediately after `!`; prose always has a space.

Python `#` comments:
  `# prose text`  -- prose comment; reflow candidate.
  `## code`       -- commented-out code; left alone.
  Commented-out code must be marked with a leading `##`
  so the reflow engine will never touch it.

Python docstrings:
  A triple-quoted string is treated as a docstring only
  when its opening quote is the first non-whitespace on
  its line (optionally preceded by an `r`, `f`, `b`, or
  `u` string prefix).  This rule excludes strings that are
  assigned to variables, passed as arguments, or used as
  return values -- those typically hold structured content
  (argparse help, HTML templates, error messages) that is
  better managed by hand.

  Prose paragraphs inside a docstring are reflowed.
  Structural content (bullet lists, numbered lists,
  doctest blocks, section headers, RST directives, table
  rows, and numpy-style `name : type` field definitions)
  is detected and left untouched so that the structural
  meaning of the docstring is preserved.

  The first line of a docstring (the summary line that
  shares the opening triple quote) and the line that holds
  the closing triple quote are never reflowed, because
  either move would require rewriting the triple-quote
  delimiters themselves.

Usage:
    rewrap_prose.py FILE [--check] [--lines S-E]

Options:
    --check       Report blocks needing reflow; do not
                  modify the file.  Exit status is 1 when
                  any reflow would be needed, 0 otherwise.
    --lines S-E   Only process lines S through E
                  (1-based, inclusive).  Lines outside
                  this range are neither examined nor
                  touched.
"""

import sys
import os
import re
import io
import tokenize
from collections import namedtuple

MAX_COL = 80
MIN_COL = 70


# ============================================================
# Comment syntax description
# ============================================================
# Both Fortran and Python share the same paragraph-detection
# and reflow machinery; the only differences are the comment
# character and the heuristic that distinguishes commented-
# out code from prose.  A CommentSyntax object bundles the
# language-specific bits so the shared helpers can stay
# language-agnostic.

CommentSyntax = namedtuple(
    'CommentSyntax',
    ['char', 'banner_re', 'is_code'])


# Fortran convention: no space after `!` marks the line as
# commented-out code.  Prose always has a space between
# the `!` and the first word.
#
# The banner regex only matches "pure fill" separators
# such as `!!!!!!!!` or `! ==========`; the trailing
# `\s*$` anchor forbids any text on the banner line.
# This is deliberately narrower than the Python variant:
# Fortran math comment blocks in this codebase frequently
# lead with an ASCII-art divider such as
# `! ----- Case 1: e1 <= E < e2` that is mathematically
# significant rather than decorative, and must not be
# treated as a structural break.
FORTRAN_SYNTAX = CommentSyntax(
    char='!',
    banner_re=re.compile(r'^\s*!\s*[!=>*\-+]+\s*$'),
    is_code=lambda rest: bool(rest) and rest[0] != ' ',
)


# Python convention: a leading `##` marks commented-out
# code.  A single `#` (with or without trailing space)
# marks prose.  The `##` prefix is required so that the
# reflow engine will never touch a line that was
# deliberately commented out for later reactivation.
#
# The banner regex matches the same two forms as its
# Fortran counterpart -- pure fill separators such as
# `########` or `# ===========` and bannered headers
# such as `# ----- Composition -----`.
PYTHON_SYNTAX = CommentSyntax(
    char='#',
    banner_re=re.compile(r'^\s*#\s*[#=>*\-+]{3,}'),
    is_code=lambda rest: rest.startswith('#'),
)


# ============================================================
# Pure-comment-line classification (shared)
# ============================================================

def is_pure_comment(line, syntax):
    """True when the only content on LINE is a comment --
    optional leading whitespace followed immediately by the
    comment character.  Used to identify the start of a
    Pattern-B comment block (two or more consecutive pure
    comment lines sharing the same prefix)."""
    return line.lstrip().startswith(syntax.char)


def is_banner(line, syntax):
    """True for decorative separator lines such as
    '   !!!!!!!' (Fortran) or '   #######' (Python).
    Banner lines are never joined into a prose paragraph;
    they terminate any run that touches them."""
    return bool(syntax.banner_re.match(line))


def is_commented_code(line, syntax):
    """True when the comment line looks like commented-out
    source code rather than prose.  The per-language rule
    lives on the CommentSyntax object above."""
    after = line.lstrip()
    if not after.startswith(syntax.char):
        return False
    rest = after[1:]
    return syntax.is_code(rest)


def is_blank_comment(line, syntax):
    """True for comment lines with no text payload, such
    as a bare `!` or `   #   ` used as a paragraph break
    inside an otherwise continuous comment block."""
    after = line.lstrip()
    if not after.startswith(syntax.char):
        return False
    return after[1:].strip() == ''


def should_skip(line, syntax):
    """True when LINE is a banner, commented-out code, or
    a blank comment.  Such lines terminate any in-progress
    prose paragraph rather than being absorbed into it."""
    return (is_banner(line, syntax)
            or is_commented_code(line, syntax)
            or is_blank_comment(line, syntax))


# ============================================================
# Prefix parsing for comment lines
# ============================================================

def parse_comment(line, syntax):
    """Split a pure comment line into (prefix, text).

    The prefix contains the leading whitespace, the comment
    character, and any spaces between the comment character
    and the first word of the text.  Example (Python)::

        '   #   some text'  ->  ('   #   ', 'some text')

    The reflow engine uses the prefix to line up every
    output line with the original indentation."""
    stripped = line.lstrip()
    indent = len(line) - len(stripped)
    after_char = stripped[1:]
    text = after_char.lstrip(' ')
    spaces = len(after_char) - len(text)
    prefix = (' ' * indent) + syntax.char + (' ' * spaces)
    return prefix, text.rstrip()


def find_inline_marker(line, syntax):
    """Return the character index of the first `!` or `#`
    that begins an inline comment on LINE -- that is, a
    comment preceded on the same line by at least one
    non-whitespace token of code.  Returns -1 if no inline
    comment is found.

    String literals (both single and double quoted) are
    respected so that a `#` or `!` inside a string does not
    confuse the scanner.  Backslash escapes inside strings
    are also skipped; this matters for Python source, and
    is harmless for Fortran (which does not use backslash
    escapes inside string literals in the first place)."""
    in_str = False
    quote = None
    i = 0
    while i < len(line):
        ch = line[i]
        if in_str:
            if ch == '\\' and i + 1 < len(line):
                i += 2
                continue
            if ch == quote:
                in_str = False
            i += 1
            continue
        if ch in ('"', "'"):
            in_str = True
            quote = ch
            i += 1
            continue
        if ch == syntax.char:
            if line[:i].strip():
                return i
            return -1
        i += 1
    return -1


def parse_inline(line, marker_idx, syntax):
    """Split a code+comment line at MARKER_IDX into
    (prefix, text).  The prefix includes the code portion
    on the left of the comment character, the comment
    character itself, and any whitespace between the
    comment character and the first word of the comment
    payload."""
    after = line[marker_idx + 1:]
    text = after.lstrip(' ')
    spaces = len(after) - len(text)
    prefix = line[:marker_idx + 1] + (' ' * spaces)
    return prefix, text.rstrip()


# ============================================================
# Comment-block paragraph finder (Pattern A and Pattern B)
# ============================================================

def find_comment_paragraphs(lines, syntax,
                            skip_ranges=None):
    """Identify comment paragraphs that may need reflow.

    Yields (start, end, first_prefix, cont_prefix) tuples.
    Indices are 0-based; end is exclusive.  first_prefix is
    used for the first output line and cont_prefix for all
    subsequent lines; the two differ when the paragraph
    began as an inline comment, or when a hanging-indent
    comment block is detected.

    Two patterns are recognised:

      Pattern A -- a code line ending with an inline
        comment, followed by one or more pure comment
        continuation lines whose prefix matches.

      Pattern B -- a run of two or more pure comment lines
        that share the same base prefix.  A small hanging
        indent (up to two extra spaces after the comment
        character on continuation lines) is allowed so
        that wrapped bullet items survive.

    When SKIP_RANGES is supplied it must be an iterable of
    (lo, hi) tuples naming half-open line-index ranges the
    paragraph finder will ignore.  The Python driver uses
    this to avoid mis-parsing comments that live inside a
    docstring body (such comments are not really comments
    at all -- they are ordinary string content)."""
    n = len(lines)
    i = 0

    def in_skip(idx):
        if not skip_ranges:
            return False
        return any(lo <= idx < hi
                   for lo, hi in skip_ranges)

    while i < n:
        if in_skip(i):
            i += 1
            continue

        line = lines[i]

        # --- Pattern A: inline comment + continuation ---
        marker = find_inline_marker(line, syntax)
        if marker >= 0:
            first_pf, _ = parse_inline(
                line, marker, syntax)

            if (i + 1 < n
                    and not in_skip(i + 1)
                    and is_pure_comment(
                        lines[i + 1], syntax)
                    and not should_skip(
                        lines[i + 1], syntax)):
                cont_pf, _ = parse_comment(
                    lines[i + 1], syntax)
                start = i
                i += 2
                while (i < n
                       and not in_skip(i)
                       and is_pure_comment(
                           lines[i], syntax)
                       and not should_skip(
                           lines[i], syntax)):
                    pf, _ = parse_comment(
                        lines[i], syntax)
                    if pf != cont_pf:
                        break
                    i += 1
                yield (start, i, first_pf, cont_pf)
                continue
            i += 1
            continue

        # --- Pattern B: pure comment block ---
        if (is_pure_comment(line, syntax)
                and not should_skip(line, syntax)):
            first_pf, _ = parse_comment(line, syntax)
            start = i
            i += 1

            # Detect continuation prefix from line 2.
            # Accept only when the base indent is the same
            # AND the continuation has >= as many spaces
            # after the comment character as the first
            # line.  Fewer spaces means a DECREASE in the
            # hanging indent, which signals a new item
            # rather than a continuation of the current
            # one.  The hanging increase is capped at two
            # spaces to avoid accidentally joining a
            # deeply-indented sub-bullet into its parent.
            cont_pf = first_pf
            if (i < n
                    and not in_skip(i)
                    and is_pure_comment(
                        lines[i], syntax)
                    and not should_skip(
                        lines[i], syntax)):
                pf2, _ = parse_comment(lines[i], syntax)
                ind1 = len(line) - len(line.lstrip())
                ind2 = (len(lines[i])
                        - len(lines[i].lstrip()))
                if (ind1 == ind2
                        and len(pf2) >= len(first_pf)
                        and (len(pf2) - len(first_pf)
                             <= 2)):
                    cont_pf = pf2
                    i += 1

            while (i < n
                   and not in_skip(i)
                   and is_pure_comment(lines[i], syntax)
                   and not should_skip(
                       lines[i], syntax)):
                pf, _ = parse_comment(lines[i], syntax)
                if pf != cont_pf:
                    break
                i += 1

            if i - start >= 2:
                yield (start, i, first_pf, cont_pf)
            continue

        i += 1


# ============================================================
# Python docstring detection
# ============================================================

# An opening triple-quote is treated as a docstring only
# when the line's only leading content is whitespace plus
# an optional string prefix (r, f, b, u, and combinations).
# This rule excludes `x = """..."""`, `foo("""...""")`,
# `return """..."""`, and every other context where the
# triple-quoted string is an expression rather than a
# statement at module, class, or function scope.
_DOCSTRING_OPEN_RE = re.compile(
    r'^(\s*)([rRfFbBuU]{0,2})("""|\'\'\')')


# Structural markers that must never be reflowed into a
# surrounding prose paragraph.  Any line that matches one
# of these patterns is left untouched and terminates any
# paragraph it would otherwise continue.
#
# Two flavours are kept: _FIRST is consulted at the start
# of a paragraph and includes parenthesized enumerations
# like `(1)` or `(a)`, which always mean enumeration when
# they are the first token of a new paragraph.  _CONT is
# consulted for continuation lines and deliberately omits
# the parenthesized-enum pattern so that prose sentences
# may legitimately wrap onto a line whose next token is a
# mid-sentence parenthetical such as `(2) an integer that
# signifies...`.
_STRUCT_LIST_FIRST_RE = re.compile(
    r'^\s*('
    r'[-*+]\s|'            # bullet list (- / * / +)
    r'\d+[.)]\s|'          # numbered list (1. / 1))
    r'\([\w\d]+\)\s|'      # parenthesized enum (a), (1)
    r'>>>\s?|'             # doctest prompt
    r'\.\.\s|'             # RST directive
    r'\|'                  # table row
    r')')

_STRUCT_LIST_CONT_RE = re.compile(
    r'^\s*('
    r'[-*+]\s|'
    r'\d+[.)]\s|'
    r'>>>\s?|'
    r'\.\.\s|'
    r'\|'
    r')')


# Lines composed entirely of ---/===/~~~/^^^ act as
# underlines or horizontal rules beneath a section header
# in numpy-style and RST-style docstrings.
_STRUCT_UNDERLINE_RE = re.compile(
    r'^\s*[-=~^]{3,}\s*$')


# A label line looks like `Keyword:` or `Required Value:
# foo`: a short run of word-like characters (optionally
# containing spaces, dots, or underscores) followed by a
# colon and either a space-plus-content or end-of-line.
# The label run must contain at least three characters so
# that two-letter prose words (`as:`, `NB:`) are not
# mis-classified, and it is capped at 25 characters so
# that long prose sentences containing a mid-sentence
# colon do not sweep in as labels either.  Label lines
# neither start nor continue a reflowable paragraph,
# because merging them would destroy the structural
# meaning of the colon.
_LABEL_LINE_RE = re.compile(
    r'^\s*\w[\w_\. \t]{2,25}:(\s|$)')


# A dash-dash field line looks like
# `.f90, .F90  -- Fortran: comment blocks.`
# or `-i FILE  -- input file name`.  These are the
# ASCII equivalent of a numpy `name : type` row where
# the separator is ` -- ` instead of ` : `.  We require
# TWO OR MORE spaces before the `--` so that prose
# em-dashes like `return values -- those typically
# hold...` (which have a single space before `--`) are
# not treated as field separators.  Field-list items
# almost always use extra whitespace for column
# alignment, so the double-space test is a reliable
# way to tell them apart.
_DASH_FIELD_LINE_RE = re.compile(
    r'^\s*\S[^\n]{0,28}?  +--\s\S')


def is_structural_first_line(line):
    """True when LINE would be a bad choice as the first
    line of a reflowable paragraph.  Matches every kind
    of structural marker, including parenthesized
    enumerations like `(1)` or `(a)` that legitimately
    start new items when they begin a fresh paragraph."""
    if _STRUCT_LIST_FIRST_RE.match(line):
        return True
    if _STRUCT_UNDERLINE_RE.match(line):
        return True
    if _LABEL_LINE_RE.match(line):
        return True
    if _DASH_FIELD_LINE_RE.match(line):
        return True
    return False


def is_structural_continuation(line):
    """True when LINE must NOT continue a paragraph that
    is already in progress.  Parenthesized enumerations
    are deliberately accepted as continuations here so
    that prose sentences may wrap onto a line whose next
    word is a mid-sentence parenthetical of the form
    `(2) an integer that signifies...`."""
    if _STRUCT_LIST_CONT_RE.match(line):
        return True
    if _STRUCT_UNDERLINE_RE.match(line):
        return True
    if _LABEL_LINE_RE.match(line):
        return True
    if _DASH_FIELD_LINE_RE.match(line):
        return True
    return False


# Regex used by looks_like_prose() below.  A line counts as
# prose-like when it contains a "sentence terminator" --
# any of `.`, `;`, `?`, `!` -- immediately followed by
# whitespace or end-of-line.  The whitespace requirement
# excludes intra-token dots such as the `.` inside `0.85`
# or `foo.bar`, which are not sentence terminators but
# would otherwise confuse the heuristic.
_SENTENCE_TERMINATOR_RE = re.compile(r'[.;?!](\s|$)')


def looks_like_prose(block_lines):
    """True when the block of lines looks like natural-
    language prose rather than tabular data, code, or
    configuration samples.  The heuristic is simple: at
    least one line in the block must contain a sentence
    terminator.  Multi-line docstring paragraphs of real
    prose almost always contain at least one period;
    `"Example input"` blocks and reactant tables almost
    never do.  The cost of a false negative (a prose
    paragraph that happens to contain no punctuation) is
    that the user has to reflow it by hand; the cost of a
    false positive (a code block rewritten as prose) is
    that the docstring content is silently corrupted.
    We lean toward the cheaper mistake."""
    for line in block_lines:
        if _SENTENCE_TERMINATOR_RE.search(line):
            return True
    return False


def find_python_docstrings(lines):
    """Scan a Python source file for multi-line docstrings.

    Returns a list of (open_line, close_line) tuples with
    0-based, INCLUSIVE line indices.

    Detection is delegated to the standard-library
    `tokenize` module so that triple-quote state is
    tracked accurately across the whole file.  A hand-
    rolled regex scanner would mis-classify the closing
    quote of a multi-line string literal (for example the
    `description_text = \"\"\"...\"\"\"` assignment used
    for argparse help in this codebase) as a new docstring
    opening and then drag every subsequent triple-quoted
    region out of alignment.

    A string token is treated as a docstring only when
    every condition below is met:
      1. it is a triple-quoted string (optionally with an
         `r`, `f`, `b`, `u`, or combined string prefix);
      2. its opening and closing quotes live on different
         lines (a one-liner has nothing to reflow);
      3. its opening quote is the first non-whitespace
         on its starting line, i.e. the characters that
         precede the quote on that line are all
         whitespace.  This excludes strings assigned to
         variables, passed as arguments, or used as
         return values, so the argparse help text in
         condense.py and its siblings is left alone."""
    regions = []
    source = '\n'.join(lines) + '\n'
    try:
        tokens = list(
            tokenize.generate_tokens(
                io.StringIO(source).readline))
    except (tokenize.TokenizeError, IndentationError,
            SyntaxError):
        return regions

    for tok in tokens:
        if tok.type != tokenize.STRING:
            continue

        # Strip the optional string-prefix letters to
        # locate the opening quote.  Only triple-quoted
        # strings are interesting here.
        payload = tok.string
        prefix_len = 0
        while (prefix_len < len(payload)
               and payload[prefix_len] in 'rRfFbBuU'):
            prefix_len += 1
        quote_start = payload[prefix_len:prefix_len + 3]
        if quote_start not in ('"""', "'''"):
            continue

        start_line = tok.start[0] - 1
        end_line = tok.end[0] - 1
        if start_line == end_line:
            continue

        # Only accept the token as a docstring when
        # nothing other than whitespace precedes its
        # opening quote on the starting line.
        start_col = tok.start[1]
        if 0 <= start_line < len(lines):
            line_prefix = lines[start_line][:start_col]
            if line_prefix.strip():
                continue

        regions.append((start_line, end_line))

    return regions


def find_docstring_paragraphs(lines, open_line,
                              close_line):
    """Find reflowable prose paragraphs inside a docstring.

    The body of the docstring spans from open_line + 1
    through close_line - 1 (0-based); the lines that hold
    the opening or closing triple-quote delimiter are
    never reflowed because doing so would require
    rewriting the delimiter position itself.

    Yields (start, end, first_prefix, cont_prefix) tuples
    in the same format as find_comment_paragraphs.  A
    reflowable paragraph is a run of two or more
    consecutive non-blank lines where every line shares
    the same leading indent and none of the lines matches
    is_structural_doc_line().  Paragraphs are bounded on
    either side by blank lines, indent changes, structural
    lines, and the docstring delimiters themselves."""
    body_start = open_line + 1
    body_end = close_line  # exclusive

    i = body_start
    while i < body_end:
        raw = lines[i]
        stripped = raw.strip()

        if not stripped:
            i += 1
            continue

        # Reject first-line candidates that are clearly
        # structural.  These lines still advance `i` so
        # the outer loop keeps searching past them for a
        # valid paragraph start.
        if is_structural_first_line(raw):
            i += 1
            continue
        if stripped.endswith(':'):
            i += 1
            continue

        first_indent = len(raw) - len(raw.lstrip())
        start = i
        i += 1

        # Collect continuation lines that share the base
        # indent of the first line and are themselves
        # non-structural.  A line ending with a colon is
        # still rejected as a continuation, because the
        # introduced section that follows it is almost
        # always at a different indent or in a structural
        # block.
        while i < body_end:
            cont = lines[i]
            cont_stripped = cont.strip()
            if not cont_stripped:
                break
            cont_indent = (len(cont)
                           - len(cont.lstrip()))
            if cont_indent != first_indent:
                break
            if is_structural_continuation(cont):
                break
            if cont_stripped.endswith(':'):
                break
            i += 1

        if i - start >= 2:
            # Final filter: only yield when the block
            # actually looks like prose.  Example-input
            # blocks, reactant tables, and similar
            # configuration samples are rejected here
            # because merging them into a single flowing
            # paragraph would corrupt the docstring.
            block = lines[start:i]
            if looks_like_prose(block):
                indent_str = ' ' * first_indent
                yield (start, i, indent_str, indent_str)


# ============================================================
# Word extraction
# ============================================================

def extract_words_comment(block, syntax):
    """Pull every word out of a comment-block paragraph.
    The first line may be either an inline comment on a
    code line or a pure comment line; subsequent lines
    must be pure comments."""
    words = []
    for idx, line in enumerate(block):
        if idx == 0:
            marker = find_inline_marker(line, syntax)
            if marker >= 0:
                _, text = parse_inline(
                    line, marker, syntax)
            else:
                _, text = parse_comment(line, syntax)
        else:
            _, text = parse_comment(line, syntax)
        words.extend(text.split())
    return words


def extract_words_plain(block):
    """Pull every word out of a plain-text paragraph, such
    as a prose paragraph inside a Python docstring.  There
    is no comment character to strip -- only the leading
    whitespace on each line."""
    words = []
    for line in block:
        words.extend(line.split())
    return words


# ============================================================
# Reflow engine (shared)
# ============================================================

def reflow(words, first_pf, cont_pf, double_space=False):
    """Reflow WORDS into a list of lines.  The first output
    line uses FIRST_PF as its leading prefix; all
    subsequent lines use CONT_PF.  Words are packed onto
    each line until adding the next word would cross the
    MAX_COL threshold, at which point we start a new line
    with CONT_PF and continue packing.

    When DOUBLE_SPACE is true, a double space is emitted
    after any word that ends in a sentence terminator
    (`.`, `!`, or `?`) so that the old-school two-space-
    after-period typography is preserved through a reflow
    pass.  The caller decides whether to enable this on a
    per-paragraph basis by looking at the incoming text:
    condense.py's module docstring uses double spaces,
    but many Fortran source comments in the codebase do
    not, and silently flipping the convention would cause
    a noisy reformat of in-band prose."""
    result = []
    pf = first_pf
    cur = pf

    def separator(prev_word):
        # Only `.` is treated as a sentence terminator
        # here.  We deliberately exclude `!` and `?` to
        # stay consistent with uses_double_space(), which
        # also only looks for `.  ` in the source text.
        if (double_space
                and prev_word
                and prev_word[-1] == '.'):
            return '  '
        return ' '

    for idx, word in enumerate(words):
        if cur == pf:
            # First word on this line -- always add it,
            # even if that pushes the line above MAX_COL
            # (the reflow contract is "do not wrap until
            # you reach 80", not "never exceed 80"; an
            # abnormally long word is allowed through).
            cur += word
            continue
        sep = separator(words[idx - 1])
        if len(cur) + len(sep) + len(word) <= MAX_COL:
            cur += sep + word
        else:
            result.append(cur)
            pf = cont_pf
            cur = pf + word

    if cur.rstrip():
        result.append(cur)

    return result


def uses_double_space(block_lines):
    """True when the block uses the old-school double-
    space-after-period convention.  Detection is local to
    the paragraph: any in-line occurrence of `.  ` (period
    followed by two or more spaces) is taken as evidence
    that the convention is in force, so reflow will
    preserve it.  Paragraphs without any such pair use
    single-space spacing, matching what the author of the
    source wrote.

    Only `.` is checked, not `!` or `?`.  Those two marks
    coincide with the Fortran comment character and with
    a common inflection, and the `!   text` hanging-indent
    prefix of Fortran prose blocks would otherwise be
    mis-read as "exclamation followed by two spaces"."""
    for line in block_lines:
        if '.  ' in line:
            return True
    return False


def needs_reflow(lines):
    """True if the paragraph would benefit from reflow:
    either some line exceeds MAX_COL, or some non-final
    line is below MIN_COL.  The last line of a paragraph
    is allowed to be short because it may be a natural
    remainder with nothing more available to fill it."""
    for idx, line in enumerate(lines):
        ln = len(line.rstrip())
        if ln > MAX_COL:
            return True
        if idx < len(lines) - 1 and ln < MIN_COL:
            return True
    return False


# ============================================================
# Language-specific paragraph drivers
# ============================================================

Paragraph = namedtuple(
    'Paragraph',
    ['start', 'end', 'first_pf', 'cont_pf', 'words',
     'double_space'])


def find_fortran_paragraphs(lines):
    """Collect every reflowable paragraph in a Fortran
    source file.  Only `!` comment blocks are examined;
    Fortran has no equivalent of Python's triple-quoted
    docstring construct."""
    paragraphs = []
    for s, e, fp, cp in find_comment_paragraphs(
            lines, FORTRAN_SYNTAX):
        block = lines[s:e]
        words = extract_words_comment(
            block, FORTRAN_SYNTAX)
        ds = uses_double_space(block)
        paragraphs.append(
            Paragraph(s, e, fp, cp, words, ds))
    return paragraphs


def find_python_paragraphs(lines):
    """Collect every reflowable paragraph in a Python
    source file.  Both `#` comment blocks (outside any
    docstring) and prose paragraphs inside multi-line
    triple-quoted docstrings are examined.  The two lists
    are merged and sorted so that downstream processing
    can apply the changes in reverse file order without
    having to special-case either source."""
    doc_regions = find_python_docstrings(lines)

    # Build a skip list so the comment-paragraph finder
    # never tries to parse the body of a docstring as a
    # run of `#` comments.  The docstring body belongs
    # to the docstring paragraph finder below.
    skip_ranges = [(o, c + 1) for o, c in doc_regions]

    paragraphs = []

    for s, e, fp, cp in find_comment_paragraphs(
            lines, PYTHON_SYNTAX,
            skip_ranges=skip_ranges):
        block = lines[s:e]
        words = extract_words_comment(
            block, PYTHON_SYNTAX)
        ds = uses_double_space(block)
        paragraphs.append(
            Paragraph(s, e, fp, cp, words, ds))

    for open_line, close_line in doc_regions:
        for s, e, fp, cp in find_docstring_paragraphs(
                lines, open_line, close_line):
            block = lines[s:e]
            words = extract_words_plain(block)
            ds = uses_double_space(block)
            paragraphs.append(
                Paragraph(s, e, fp, cp, words, ds))

    paragraphs.sort(key=lambda p: p.start)
    return paragraphs


# ============================================================
# File processing
# ============================================================

def process_file(path, paragraph_finder,
                 check_only=False, line_range=None):
    """Find and reflow paragraphs in a source file.

    PARAGRAPH_FINDER is a callable that takes the list of
    file lines and returns a list of Paragraph records.
    Returns 0 if no changes were needed (or changes were
    applied), 1 if --check found violations that would be
    reflowed on a non-dry run."""
    with open(path, 'r') as fh:
        raw = fh.readlines()

    orig = [l.rstrip('\n').rstrip('\r') for l in raw]
    result = list(orig)

    paras = paragraph_finder(orig)

    if line_range:
        lo = line_range[0] - 1
        hi = line_range[1]
        paras = [p for p in paras
                 if p.start >= lo and p.end <= hi]

    changes = []

    # Process in reverse so line-number shifts from earlier
    # replacements do not cascade into later ones.
    for p in reversed(paras):
        block = orig[p.start:p.end]
        if not needs_reflow(block):
            continue
        if not p.words:
            continue

        new_block = reflow(
            p.words, p.first_pf, p.cont_pf,
            double_space=p.double_space)
        if new_block == [l.rstrip() for l in block]:
            continue

        changes.append(
            (p.start + 1, p.end, len(block),
             len(new_block)))

        if not check_only:
            result[p.start:p.end] = new_block

    changes.reverse()

    if check_only:
        if changes:
            for s, e, on, nn in changes:
                print(f"  L{s}-{e}: {on} -> {nn} lines")
            return 1
        print("  No reflow needed.")
        return 0

    if changes:
        with open(path, 'w') as fh:
            for line in result:
                fh.write(line + '\n')
        for s, e, on, nn in changes:
            print(f"  L{s}-{e}: {on} -> {nn} lines")
        return 0

    print("  No changes needed.")
    return 0


# ============================================================
# CLI
# ============================================================

# Extension -> paragraph-finder dispatch table.  Adding a
# new language is a matter of implementing its finder and
# registering it here.
_FINDERS = {
    '.f90': find_fortran_paragraphs,
    '.py': find_python_paragraphs,
}


def main():
    """Entry point.  Parse arguments and dispatch to the
    language-specific paragraph finder."""
    args = sys.argv[1:]
    if not args or args[0] in ('-h', '--help'):
        print(__doc__)
        sys.exit(0)

    path = args[0]
    check = '--check' in args

    lr = None
    for i, a in enumerate(args):
        if a == '--lines' and i + 1 < len(args):
            parts = args[i + 1].split('-')
            lr = (int(parts[0]), int(parts[1]))

    ext = os.path.splitext(path)[1].lower()
    finder = _FINDERS.get(ext)
    if finder is None:
        print(f"Unsupported extension: {ext}",
              file=sys.stderr)
        sys.exit(1)

    sys.exit(process_file(path, finder, check, lr))


if __name__ == '__main__':
    main()
