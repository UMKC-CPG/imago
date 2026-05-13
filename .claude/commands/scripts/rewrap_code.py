#!/usr/bin/env python3
"""rewrap_code.py -- Mechanically unwrap over-wrapped code.

Two modes of operation:

  Scan mode (default):
    rewrap_code.py SOURCE [--check]

  Scans the file for over-wrapped code blocks (3+ lines
  that could fit on fewer).  For each block, tries unwrap
  (single line) first, then rewrap (comma-split).  Skips
  blocks that cannot be compacted within 80 chars.

  Directive mode:
    rewrap_code.py SOURCE --noscan DIRECTIVES [--check]

  Reads explicit directives from a file and applies them.
  Directive format (one per line, # comments allowed):
      unwrap START END
      rewrap START END
      merge-strings START END
      detriplify START END

  unwrap: join lines into one, strip gratuitous parens.
    Fails if result exceeds 80 chars.
  rewrap: join then re-split at top-level commas.
    Falls back to unwrap when joined line fits.
  merge-strings: concatenate adjacent string literals
    and re-split to fill 70-80 per line.
  detriplify: concatenate a run of adjacent string
    literals that together encode preformatted multi-
    line text, and rewrite the enclosing call as a
    single triple-quoted string literal.  Scan mode
    triggers this operation automatically for blocks
    where the call's single argument is a run of FIVE
    OR MORE adjacent string literals and at least one
    of them contains a `\\n` escape.

  Line numbers are 1-based, inclusive.

Options:
    --noscan  Directive mode: requires a directives file.
    --check   Dry run.  Report changes without modifying
              the file.  Exit 1 if any changes needed.
"""

import sys
import os
import re

MAX_COL = 80
MIN_COL = 70


# ============================================================
# Utility helpers
# ============================================================

def has_top_level_comma(text):
    """True when text contains a comma that is not nested
    inside any brackets.  Used to distinguish tuples like
    (a, b) from parenthesized expressions like (a + b)."""
    depth = 0
    for ch in text:
        if ch in '([{':
            depth += 1
        elif ch in ')]}':
            depth -= 1
        elif ch == ',' and depth == 0:
            return True
    return False


def parens_balanced(text):
    """True when every ( [ { in text has a matching closer
    and no closer appears before its opener."""
    depth = 0
    for ch in text:
        if ch in '([{':
            depth += 1
        elif ch in ')]}':
            depth -= 1
            if depth < 0:
                return False
    return depth == 0


# ============================================================
# unwrap: join lines into one, strip gratuitous parens
# ============================================================

def strip_gratuitous_parens(joined):
    """Remove wrapper parens that exist only for continuation.

    Recognized patterns and their rewrites:
      name = (expr)    ->  name = expr
        when expr has no top-level comma (not a tuple)
        and all brackets inside expr are balanced.
      if (expr):       ->  if expr:
      elif (expr):     ->  elif expr:
      while (expr):    ->  while expr:
    Returns the (possibly rewritten) string.
    """
    # Assignment:  anything = (expr)
    m = re.match(r'^(.+=\s*)\((.+)\)\s*$', joined)
    if m:
        lhs, inner = m.group(1), m.group(2)
        if (not has_top_level_comma(inner)
                and parens_balanced(inner)):
            return lhs + inner

    # Control-flow keyword:  if/elif/while (expr):
    for kw in ('if', 'elif', 'while'):
        prefix = kw + ' ('
        if (joined.startswith(prefix)
                and joined.endswith('):')):
            inner = joined[len(prefix):-2]
            if parens_balanced(inner):
                return kw + ' ' + inner + ':'

    return joined


def do_unwrap(lines, start, end):
    """Join lines[start-1..end-1] into a single line.

    Strips gratuitous parentheses that were added solely
    to enable line continuation.

    Parameters
    ----------
    lines : list of str
        Full file contents (0-indexed).
    start, end : int
        1-based inclusive line range to join.

    Returns
    -------
    result : str
        The single joined line.
    fit : bool
        True when len(result) <= MAX_COL.
    """
    block = [lines[i].rstrip()
             for i in range(start - 1, end)]
    first = block[0]
    indent = len(first) - len(first.lstrip())
    indent_str = first[:indent]
    stripped = [ln.strip() for ln in block]

    joined = smart_join(stripped)
    joined = strip_gratuitous_parens(joined)
    result = indent_str + joined
    return result, len(result) <= MAX_COL


# ============================================================
# smart_join: shared join logic for unwrap and rewrap
# ============================================================

OPENERS = ('(', '[')
CLOSERS = (')', ']')

# Leading characters on a continuation line that must
# NOT be preceded by a separating space when the lines
# are joined.
#
# `.` is the member-access operator and must bind tightly
# to its target (`self.foo`, never `self .foo`).
# `[` is the subscript operator and binds the same way
# (`table[key]`, never `table [key]`).  Both show up on
# continuation lines when a programmer wraps an attribute
# chain or an indexing expression across several source
# lines.
#
# `,`, `:`, and `;` are separator / terminator punctuation
# and never take a leading space in Python.  Listing them
# here means that any continuation line that begins with
# one of these (a bit unusual but legal) is joined flush
# against the previous token, yielding `f(a, b)` instead
# of `f(a , b)` and `if x: y = 1` instead of
# `if x : y = 1`.
NO_SPACE_PREFIX = ('.', '[', ',', ':', ';')


def smart_join(stripped_parts):
    """Join stripped line parts with smart spacing.

    Suppresses the space between an opening bracket and
    the next token, and between the last token and a
    closing bracket.  Also suppresses the space in front
    of a continuation that begins with `.` (member
    access) or `[` (subscript), so that an expression
    like::

        self
        .molecule_to_family_map
        [mol_nm]

    unwraps to `self.molecule_to_family_map[mol_nm]` and
    not the syntactically-legal but ugly variant
    `self .molecule_to_family_map [mol_nm]`.  Handles
    both parentheses and square brackets.
    """
    joined = stripped_parts[0]
    for part in stripped_parts[1:]:
        if not part:
            continue
        if joined and joined[-1] in OPENERS:
            joined += part
        elif part[0] in CLOSERS:
            joined += part
        elif part[0] in NO_SPACE_PREFIX:
            joined += part
        else:
            joined += ' ' + part
    return joined


# ============================================================
# rewrap: join then re-split at top-level commas
# ============================================================

def find_top_level_commas(text):
    """Return indices of all commas at bracket depth 0.

    Skips commas inside nested parentheses, brackets, or
    braces, and commas inside string literals.  Returns
    a list of 0-based character positions.
    """
    positions = []
    depth = 0
    in_str = False
    quote = None
    i = 0
    while i < len(text):
        ch = text[i]
        if in_str:
            if ch == '\\' and i + 1 < len(text):
                i += 2
                continue
            if ch == quote:
                in_str = False
        else:
            if ch in ('"', "'"):
                in_str = True
                quote = ch
            elif ch in '([{':
                depth += 1
            elif ch in ')]}':
                depth -= 1
            elif ch == ',' and depth == 0:
                positions.append(i)
        i += 1
    return positions


def do_rewrap(lines, start, end):
    """Join lines then re-split at commas to fit MAX_COL.

    Algorithm:
      1. Join exactly like do_unwrap (smart join, strip
         gratuitous parens).
      2. If the joined result fits in MAX_COL, return it
         as a single line (same as unwrap).
      3. Otherwise, find the outermost open-paren on the
         first stripped token, and split the argument
         list at top-level commas.  Pack arguments
         greedily onto lines up to MAX_COL.  Continuation
         lines are indented 4 past the statement indent.

    Parameters
    ----------
    lines : list of str
        Full file contents (0-indexed).
    start, end : int
        1-based inclusive line range.

    Returns
    -------
    new_lines : list of str
        Replacement lines.
    fit : bool
        True when every output line is <= MAX_COL.
    """
    # Step 1: join (reuse do_unwrap internals).
    result, fit = do_unwrap(lines, start, end)
    if fit:
        return [result], True

    # Step 2: we need to re-split.  Work on the stripped
    # content (no leading indent).
    block = [lines[i].rstrip()
             for i in range(start - 1, end)]
    first = block[0]
    indent = len(first) - len(first.lstrip())
    indent_str = first[:indent]
    stripped = [ln.strip() for ln in block]

    joined = smart_join(stripped)
    joined = strip_gratuitous_parens(joined)

    # Find the first top-level open paren to determine
    # where the argument list starts.
    paren_pos = -1
    depth = 0
    in_str = False
    quote = None
    for i, ch in enumerate(joined):
        if in_str:
            if ch == '\\':
                continue
            if ch == quote:
                in_str = False
            continue
        if ch in ('"', "'"):
            in_str = True
            quote = ch
        elif ch == '(' and depth == 0:
            paren_pos = i
            break
        elif ch in '([{':
            depth += 1
        elif ch in ')]}':
            depth -= 1

    if paren_pos < 0:
        # No paren found -- cannot re-split intelligently.
        return [indent_str + joined], False

    # Everything up to and including the '(' is the opener.
    opener = joined[:paren_pos + 1]

    # Find the matching close paren.
    inner_start = paren_pos + 1
    depth = 1
    in_str = False
    quote = None
    close_pos = len(joined)
    for i in range(inner_start, len(joined)):
        ch = joined[i]
        if in_str:
            if ch == '\\':
                continue
            if ch == quote:
                in_str = False
            continue
        if ch in ('"', "'"):
            in_str = True
            quote = ch
        elif ch in '([{':
            depth += 1
        elif ch in ')]}':
            depth -= 1
            if depth == 0:
                close_pos = i
                break

    inner = joined[inner_start:close_pos]
    closer = joined[close_pos:]

    # Split inner at top-level commas.
    comma_positions = find_top_level_commas(inner)

    if not comma_positions:
        # No commas to split on -- nothing we can do.
        return [indent_str + joined], False

    # Build argument fragments.  Each fragment includes
    # the comma at the end (except the last one).
    # Filter out empty trailing fragments that arise
    # from trailing commas like func(a, b,).
    args = []
    prev = 0
    for cp in comma_positions:
        args.append(inner[prev:cp + 1].strip())
        prev = cp + 1
    tail = inner[prev:].strip()
    if tail:
        args.append(tail)

    # Strip trailing comma from the last arg -- it looks
    # odd jammed against the closing delimiter and serves
    # no purpose once lines are compacted.
    if args and args[-1].endswith(','):
        args[-1] = args[-1][:-1]

    # Continuation indent: 4 spaces past the statement
    # indent.
    cont_indent_str = indent_str + '    '
    cont_indent = len(cont_indent_str)

    # Strategy: try to put opener + first args on line 1,
    # then pack remaining args onto continuation lines.
    # When an argument contains a parenthesized string
    # block (e.g. help=("str1" "str2")), allow it to
    # span multiple lines by splitting at string
    # boundaries within the parens.
    output = []
    budget = MAX_COL - indent - len(opener)
    cur = indent_str + opener
    first_on_line = True

    for i, arg in enumerate(args):
        is_last = (i == len(args) - 1)
        # What goes after this arg if it's last?
        suffix = closer if is_last else ''
        token = arg + suffix if is_last else arg

        if first_on_line:
            candidate = cur + token
            if len(candidate) <= MAX_COL:
                cur = candidate
                first_on_line = False
            else:
                # Token too long for this line.  Try
                # splitting paren-strings within it.
                split = _split_paren_string_arg(
                    token, cur, cont_indent_str)
                if split is not None:
                    output.extend(split[:-1])
                    cur = split[-1]
                    first_on_line = False
                else:
                    # Opener alone, arg on next.  Try
                    # paren-string split on new line.
                    output.append(cur)
                    split = _split_paren_string_arg(
                        token, cont_indent_str,
                        cont_indent_str + '    ')
                    if split is not None:
                        output.extend(split[:-1])
                        cur = split[-1]
                    else:
                        cur = cont_indent_str + token
                    first_on_line = False
        else:
            candidate = cur + ' ' + token
            if len(candidate) <= MAX_COL:
                cur = candidate
            else:
                # Token too long for this line.  Try
                # splitting paren-strings within it
                # on the current line first.
                split = _split_paren_string_arg(
                    token, cur + ' ', cont_indent_str)
                if split is not None:
                    output.extend(split[:-1])
                    cur = split[-1]
                else:
                    # Start on a new continuation line
                    # and try paren-string split there.
                    output.append(cur)
                    split = _split_paren_string_arg(
                        token, cont_indent_str,
                        cont_indent_str + '    ')
                    if split is not None:
                        output.extend(split[:-1])
                        cur = split[-1]
                    else:
                        cur = cont_indent_str + token

    # Flush the last line.
    if not cur.endswith(closer):
        cur = cur + closer
    output.append(cur)

    ok = all(len(ln) <= MAX_COL for ln in output)
    return output, ok


def _split_paren_string_arg(token, line_prefix,
                            cont_indent_str):
    """Split an argument with parenthesized strings.

    Handles patterns like:
        name=("string1" "string2" f"string3")
    by placing the opening portion on the current line
    and continuation strings on subsequent lines indented
    4 past the continuation indent.

    Parameters
    ----------
    token : str
        The argument text, e.g.
        'help=("Name of file. " f"Default: {x}")'.
    line_prefix : str
        What's already on the current line (including
        indent and any prior tokens).
    cont_indent_str : str
        The continuation indent string.

    Returns
    -------
    list of str or None
        Output lines if splitting succeeded, else None.
    """
    # Look for  =( pattern that opens a string block.
    m = re.search(r'=\s*\(', token)
    if m is None:
        return None

    # Everything before and including the '(' .
    paren_open = m.end()
    prefix_part = token[:paren_open]

    # Find the matching ')'.
    depth = 1
    in_str = False
    quote = None
    close_idx = len(token)
    for idx in range(paren_open, len(token)):
        ch = token[idx]
        if in_str:
            if ch == '\\':
                continue
            if ch == quote:
                in_str = False
            continue
        if ch in ('"', "'"):
            in_str = True
            quote = ch
        elif ch in '([{':
            depth += 1
        elif ch in ')]}':
            depth -= 1
            if depth == 0:
                close_idx = idx
                break

    inner = token[paren_open:close_idx]
    after_close = token[close_idx:]

    # The inner content should be string literals
    # separated by whitespace (implicit concatenation).
    # Split into individual string fragments.
    fragments = _split_string_fragments(inner.strip())
    if fragments is None or len(fragments) < 2:
        return None

    # Build output lines.  Put prefix_part + first
    # fragment(s) on the current line, then continue.
    str_indent = cont_indent_str + '    '
    result_lines = []

    first_line = line_prefix + prefix_part
    frag_idx = 0

    # Pack as many fragments as fit on the first line.
    while frag_idx < len(fragments):
        candidate = first_line + fragments[frag_idx]
        # If this is the last fragment, add the closer.
        if frag_idx == len(fragments) - 1:
            candidate += after_close
        if len(candidate) <= MAX_COL:
            first_line = candidate
            if frag_idx < len(fragments) - 1:
                first_line += ' '
            frag_idx += 1
        else:
            break

    if frag_idx == 0:
        # Not even one fragment fits on the first line.
        return None

    result_lines.append(first_line.rstrip())

    # Pack remaining fragments onto continuation lines.
    while frag_idx < len(fragments):
        cur_line = str_indent
        while frag_idx < len(fragments):
            frag = fragments[frag_idx]
            candidate = cur_line + frag
            if frag_idx == len(fragments) - 1:
                candidate += after_close
            if len(candidate) <= MAX_COL:
                cur_line = candidate
                if frag_idx < len(fragments) - 1:
                    cur_line += ' '
                frag_idx += 1
            else:
                if cur_line.strip():
                    break
                # Single fragment too long; give up.
                return None
        result_lines.append(cur_line.rstrip())

    return result_lines


def _split_string_fragments(text):
    """Split implicit string concatenation into fragments.

    Given text like:
        "hello " f"world {x}" "end"
    returns:
        ['"hello "', 'f"world {x}"', '"end"']

    Returns None if the text doesn't look like string
    concatenation.
    """
    fragments = []
    i = 0
    while i < len(text):
        # Skip whitespace between fragments.
        while i < len(text) and text[i] == ' ':
            i += 1
        if i >= len(text):
            break

        # Expect a string prefix or quote.
        start = i
        # Skip string prefix letters (f, r, b, u).
        while (i < len(text)
               and text[i] in 'fFrRbBuU'):
            i += 1
        if i >= len(text) or text[i] not in ('"', "'"):
            return None

        quote = text[i]
        i += 1

        # Walk to the closing quote.
        while i < len(text):
            if text[i] == '\\' and i + 1 < len(text):
                i += 2
            elif text[i] == quote:
                i += 1
                break
            else:
                i += 1
        else:
            return None

        fragments.append(text[start:i])

    return fragments if fragments else None


# ============================================================
# merge-strings: repack string fragments to fill 70-80
# ============================================================

STRING_PREFIX_RE = re.compile(r'^([fFrRbBuU]*)(["\'])')


def extract_string_info(stripped_line):
    """Extract content from a Python string literal.

    Parameters
    ----------
    stripped_line : str
        A whitespace-stripped line that should be a string
        literal (e.g. 'f"hello {name}"').

    Returns
    -------
    (prefix, content, quote) or None
        prefix  : letter(s) before the quote (e.g. 'f').
        content : raw text between the quotes.
        quote   : the quote character (' or ").
    """
    m = STRING_PREFIX_RE.match(stripped_line)
    if not m:
        return None
    prefix, quote = m.group(1), m.group(2)
    rest = stripped_line[len(prefix) + 1:]

    # Walk to the unescaped closing quote.
    chars = []
    i = 0
    while i < len(rest):
        if rest[i] == '\\' and i + 1 < len(rest):
            chars.append(rest[i:i + 2])
            i += 2
        elif rest[i] == quote:
            return prefix, ''.join(chars), quote
        else:
            chars.append(rest[i])
            i += 1
    return None


def find_split_point(text, max_chars):
    """Find where to split text at <= max_chars characters.

    Avoids splitting inside {expr} blocks (f-string
    expressions).  Splits BEFORE a space so that the
    continuation string starts with the space character,
    matching the existing formatting convention where each
    fragment begins with a space separator.

    Returns the 0-based index where the second part begins.
    """
    if len(text) <= max_chars:
        return len(text)

    best = -1
    depth = 0
    for i, ch in enumerate(text):
        if i > max_chars:
            break
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
        elif ch == ' ' and depth == 0:
            best = i

    if best > 0:
        return best

    # No split point within budget; take the first space
    # past max_chars as a fallback.
    depth = 0
    for i in range(max_chars, len(text)):
        ch = text[i]
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
        elif ch == ' ' and depth == 0:
            return i

    return len(text)


def do_merge_strings(lines, start, end):
    """Merge adjacent string literals and repack to 70-80.

    Handles the common pattern where a function call wraps
    multiple short string-fragment lines::

        func(
            f"short fragment 1"
            f" short fragment 2"
            f" short fragment 3"
        )

    Parameters
    ----------
    lines : list of str
        Full file contents (0-indexed).
    start, end : int
        1-based inclusive line range.

    Returns
    -------
    new_lines : list of str, or None on parse failure.
    fit : bool
        True when every output line is within MAX_COL.
    """
    block = [lines[i].rstrip()
             for i in range(start - 1, end)]

    opening = None
    string_contents = []
    string_prefix = ''
    quote_char = '"'
    string_indent = None
    closing = None

    for raw_line in block:
        stripped = raw_line.strip()

        # Closing paren on its own line
        if stripped in (')', ');', '):', '),'):
            closing = stripped
            continue

        # Pure string literal line
        info = extract_string_info(stripped)
        if info is not None:
            pfx, content, qch = info
            if pfx and not string_prefix:
                string_prefix = pfx
            quote_char = qch
            string_contents.append(content)
            if string_indent is None:
                string_indent = (len(raw_line)
                                 - len(raw_line.lstrip()))
            continue

        # Opening line: ends with ( and has no string
        if stripped.endswith('('):
            opening = raw_line.rstrip()
            if string_indent is None:
                call_indent = (len(raw_line)
                               - len(raw_line.lstrip()))
                string_indent = call_indent + 4
            continue

        # Unrecognized line structure
        return None, False

    if not string_contents:
        return None, False

    merged = ''.join(string_contents)
    if string_indent is None:
        string_indent = 8

    indent_str = ' ' * string_indent
    # Per-line overhead = indent + prefix + open_quote
    #                    + close_quote
    overhead = string_indent + len(string_prefix) + 2
    max_content = MAX_COL - overhead

    # Split the merged content into lines that each fill
    # toward MAX_COL.
    content_lines = []
    pos = 0
    while pos < len(merged):
        remaining = merged[pos:]
        if len(remaining) <= max_content:
            content_lines.append(remaining)
            break
        split = find_split_point(remaining, max_content)
        content_lines.append(remaining[:split])
        pos += split

    # Assemble output lines.
    result = []
    if opening is not None:
        result.append(opening)

    for i, content in enumerate(content_lines):
        is_last = (i == len(content_lines) - 1)
        s_line = (indent_str + string_prefix
                  + quote_char + content + quote_char)
        # Try to place closing paren on the last string
        # line rather than giving it a line of its own.
        if is_last and closing is not None:
            candidate = s_line + closing
            if len(candidate) <= MAX_COL:
                s_line = candidate
                closing = None
        result.append(s_line)

    # Closing paren gets its own line only when it did
    # not fit on the last string line above.
    if closing is not None:
        close_indent = ' ' * max(0, string_indent - 4)
        result.append(close_indent + closing)

    ok = all(len(ln) <= MAX_COL for ln in result)
    return result, ok


# ============================================================
# detriplify: collapse a run of adjacent string literals
# into a single triple-quoted string
# ============================================================

# Minimum number of adjacent string literals required for
# scan mode to treat a block as a detriplify candidate.
# The threshold is high enough that short calls like
#     write("foo " "bar")
# never trigger, because the point of detriplify is to
# rescue long preformatted-text blocks that have been
# chopped across many short lines.
_DETRIPLIFY_MIN_STRINGS = 5


def do_detriplify(lines, start, end):
    """Rewrite a function call whose single argument is a
    run of adjacent string literals as one triple-quoted
    string literal.

    The block must parse as a single expression statement
    whose value is a call with exactly one positional
    argument.  That argument must itself be a plain string
    constant that Python produced by concatenating two or
    more adjacent string literals at parse time (so the
    AST collapses them into one `ast.Constant` while the
    raw token stream still shows several `STRING` tokens).

    Preconditions enforced here:
      - The AST path is exactly Expr -> Call -> Constant
        (str), with zero keyword arguments and one
        positional argument.
      - The token stream contains at least
        _DETRIPLIFY_MIN_STRINGS STRING tokens inside the
        call's argument list.
      - None of those tokens carry an `f`, `r`, or `b`
        string prefix.  An `f` prefix would cross an
        interpolation boundary; `r` and `b` change the
        semantics of escape sequences and are rejected
        defensively.  A bare `u` prefix is accepted
        because it is a no-op in Python 3.
      - The concatenated content contains at least one
        `\\n` so we only touch text that is actually
        multi-line output.
      - The concatenated content does not already contain
        `\"\"\"`, which would terminate the new literal
        prematurely.

    Returns (new_lines, True) on success, or
    (None, False) if any precondition fails."""
    import ast
    import io
    import textwrap
    import tokenize

    block = [lines[i].rstrip()
             for i in range(start - 1, end)]
    first = block[0]
    indent = len(first) - len(first.lstrip())
    indent_str = first[:indent]

    # textwrap.dedent lets ast.parse and tokenize process
    # the block without tripping over its leading indent.
    dedented = textwrap.dedent('\n'.join(block))
    try:
        tree = ast.parse(dedented)
    except SyntaxError:
        return None, False

    if len(tree.body) != 1:
        return None, False
    stmt = tree.body[0]
    if not isinstance(stmt, ast.Expr):
        return None, False
    call = stmt.value
    if not isinstance(call, ast.Call):
        return None, False
    if len(call.args) != 1 or call.keywords:
        return None, False

    arg = call.args[0]
    if not isinstance(arg, ast.Constant):
        return None, False
    content = arg.value
    if not isinstance(content, str):
        return None, False

    # Count STRING tokens and reject unsafe prefixes.
    try:
        toks = list(
            tokenize.generate_tokens(
                io.StringIO(dedented).readline))
    except tokenize.TokenizeError:
        return None, False

    n_strings = 0
    for t in toks:
        if t.type != tokenize.STRING:
            continue
        n_strings += 1
        raw = t.string
        j = 0
        while j < len(raw) and raw[j] in 'rRfFbBuU':
            j += 1
        prefix = raw[:j].lower()
        if prefix and prefix != 'u':
            return None, False

    if n_strings < _DETRIPLIFY_MIN_STRINGS:
        return None, False

    if '\n' not in content:
        return None, False
    if '"""' in content:
        return None, False

    # Reconstruct the call expression using the triple-
    # quoted literal in place of the chopped-up run.  We
    # use ast.unparse on the function expression so that
    # simple dotted calls like `lmpin.write` are round-
    # tripped accurately.
    try:
        func_source = ast.unparse(call.func)
    except Exception:
        return None, False

    opener = indent_str + func_source + '("""'
    closer = '""")'

    # Split the content at real newlines.  The opener
    # (ending in `"""`) always sits alone on its own line,
    # and the actual text always begins on the line below
    # it.  This keeps the triple-quoted block visually
    # clean regardless of whether the original content
    # started with a newline or with immediate text.  The
    # closer is attached to the final content line.
    content_lines = content.split('\n')
    result = [opener]
    for extra in content_lines:
        result.append(extra)
    result[-1] = result[-1] + closer

    return result, True


# ============================================================
# Directive parsing
# ============================================================

VALID_OPS = ('unwrap', 'rewrap', 'merge-strings',
             'detriplify')


def parse_directives(path):
    """Read a directive file.

    Returns a list of (op, start, end) tuples.  Blank
    lines and lines starting with # are skipped.
    """
    directives = []
    with open(path, 'r') as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 3:
                print(
                    f"WARNING: malformed directive on "
                    f"line {lineno}: {line!r}",
                    file=sys.stderr)
                continue
            op = parts[0]
            try:
                s, e = int(parts[1]), int(parts[2])
            except ValueError:
                print(
                    f"WARNING: bad line numbers on "
                    f"directive line {lineno}: {line!r}",
                    file=sys.stderr)
                continue
            if op not in VALID_OPS:
                print(
                    f"WARNING: unknown op '{op}' on "
                    f"directive line {lineno}",
                    file=sys.stderr)
                continue
            directives.append((op, s, e))
    return directives


# ============================================================
# Main driver
# ============================================================

def apply_directives(source_path, directives,
                     check_only):
    """Apply all directives to the source file.

    Directives are sorted by descending start line so that
    replacements do not shift the indices of later ones.
    Returns the number of changes applied (or that would
    be applied in --check mode).
    """
    with open(source_path, 'r') as fh:
        lines = [ln.rstrip('\n').rstrip('\r')
                 for ln in fh]

    sorted_dirs = sorted(directives,
                         key=lambda d: d[1],
                         reverse=True)
    changes = 0

    for op, start, end in sorted_dirs:
        if start < 1 or end > len(lines) or start > end:
            print(f"  SKIP {op} {start} {end}: out of "
                  f"range (file has {len(lines)} lines)")
            continue

        old_block = [lines[i].rstrip()
                     for i in range(start - 1, end)]

        if op == 'unwrap':
            result, fit = do_unwrap(lines, start, end)
            if not fit:
                print(
                    f"  FAIL {op} {start}-{end}: "
                    f"result is {len(result)} chars "
                    f"(> {MAX_COL})")
                continue
            if [result] == old_block:
                continue
            print(f"  {op} {start}-{end}: "
                  f"{len(old_block)} lines -> 1 "
                  f"({len(result)} chars)")
            if not check_only:
                lines[start - 1:end] = [result]
            changes += 1

        elif op == 'rewrap':
            new_lines, fit = do_rewrap(
                lines, start, end)
            if not fit:
                print(
                    f"  FAIL {op} {start}-{end}: "
                    f"could not fit within "
                    f"{MAX_COL} chars")
                continue
            if new_lines == old_block:
                continue
            print(f"  {op} {start}-{end}: "
                  f"{len(old_block)} lines -> "
                  f"{len(new_lines)}")
            if not check_only:
                lines[start - 1:end] = new_lines
            changes += 1

        elif op == 'merge-strings':
            new_lines, fit = do_merge_strings(
                lines, start, end)
            if new_lines is None:
                print(
                    f"  FAIL {op} {start}-{end}: "
                    f"could not parse string block")
                continue
            if not fit:
                over = [ln for ln in new_lines
                        if len(ln) > MAX_COL]
                print(
                    f"  FAIL {op} {start}-{end}: "
                    f"{len(over)} output line(s) "
                    f"exceed {MAX_COL} chars")
                continue
            if new_lines == old_block:
                continue
            print(f"  {op} {start}-{end}: "
                  f"{len(old_block)} lines -> "
                  f"{len(new_lines)}")
            if not check_only:
                lines[start - 1:end] = new_lines
            changes += 1

        elif op == 'detriplify':
            new_lines, ok = do_detriplify(
                lines, start, end)
            if not ok:
                print(
                    f"  FAIL {op} {start}-{end}: "
                    f"block does not match a run of "
                    f"{_DETRIPLIFY_MIN_STRINGS}+ "
                    f"adjacent string literals with "
                    f"multi-line content")
                continue
            if new_lines == old_block:
                continue
            print(f"  {op} {start}-{end}: "
                  f"{len(old_block)} lines -> "
                  f"{len(new_lines)}")
            if not check_only:
                lines[start - 1:end] = new_lines
            changes += 1

    if not check_only and changes > 0:
        with open(source_path, 'w') as fh:
            for ln in lines:
                fh.write(ln + '\n')

    return changes


def _build_triple_quote_map(lines):
    """Return a set of line indices inside triple-quoted
    strings (docstrings or multi-line string literals).

    Scans the file once, toggling a flag each time a
    triple-quote delimiter is encountered.  This is a
    rough heuristic -- it does not handle escaped triple
    quotes or mixed quote styles, but it is good enough
    for well-formed Python source.
    """
    inside = set()
    in_tq = False
    tq = None
    for i, text in enumerate(lines):
        # Check if we're inside a triple quote BEFORE
        # processing this line's delimiters.  If we
        # entered a triple quote on a previous line and
        # haven't closed it, this line is inside.
        if in_tq:
            inside.add(i)
        for delim in ('"""', "'''"):
            count = text.count(delim)
            for _ in range(count):
                if not in_tq:
                    in_tq = True
                    tq = delim
                elif delim == tq:
                    in_tq = False
                    tq = None
    return inside


def _find_detriplify_candidates(source):
    """Walk the source file's AST looking for every call
    expression whose single positional argument is a
    string constant that Python formed by concatenating a
    run of adjacent string literals at parse time.  For
    each match, we count the STRING tokens that live
    inside the call's argument span; if there are at
    least `_DETRIPLIFY_MIN_STRINGS` of them AND the
    constant value is a genuine multi-line string, the
    call is marked as a detriplify candidate.

    Returns a list of (start, end) tuples with 1-based,
    inclusive line numbers -- the same shape that the
    directive path uses, so the two paths share
    `do_detriplify` without translation."""
    import ast
    import io
    import tokenize

    try:
        tree = ast.parse(source)
    except SyntaxError:
        return []

    try:
        toks = list(
            tokenize.generate_tokens(
                io.StringIO(source).readline))
    except tokenize.TokenizeError:
        return []

    # Collect every STRING token's (start_pos, end_pos)
    # so we can later ask "how many string tokens live
    # between (lineA, colA) and (lineB, colB)?".  Tuples
    # of the form (row, col) compare lexicographically,
    # which is exactly the order we need.
    string_spans = [
        (t.start, t.end)
        for t in toks
        if t.type == tokenize.STRING
    ]

    candidates = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if len(node.args) != 1 or node.keywords:
            continue
        arg = node.args[0]
        if not isinstance(arg, ast.Constant):
            continue
        if not isinstance(arg.value, str):
            continue

        arg_start = (arg.lineno, arg.col_offset)
        arg_end = (
            arg.end_lineno,
            arg.end_col_offset,
        )

        n = 0
        for s_pos, e_pos in string_spans:
            if s_pos >= arg_start and e_pos <= arg_end:
                n += 1
        if n < _DETRIPLIFY_MIN_STRINGS:
            continue

        if '\n' not in arg.value:
            continue
        if '"""' in arg.value:
            continue

        candidates.append(
            (node.lineno, node.end_lineno))

    # Sort in file order so the caller can choose a
    # stable reverse order for in-place editing.
    candidates.sort()
    return candidates


def scan_file(source_path, check_only):
    """Scan a file for over-wrapped code blocks and fix.

    Runs two passes in order:

      1. Detriplify pass.  Uses the AST to find every
         call whose single argument is a run of five or
         more adjacent string literals encoding a real
         multi-line payload, and rewrites each one as a
         triple-quoted string.  This pass is not limited
         by the 10-line lookahead window used for
         unwrap/rewrap -- a 40-line chopped-up preformat
         block is detected just as reliably as a 7-line
         one.

      2. Over-wrapped-block pass.  For each block of 3+
         lines created by an unbalanced opening bracket
         (within a 10-line window), tries unwrap (single
         line) first, then rewrap (comma-split).  Skips
         blocks that cannot be compacted within MAX_COL.

    Returns the number of changes applied (or that would
    be applied in --check mode).
    """
    with open(source_path, 'r') as fh:
        lines = [ln.rstrip('\n').rstrip('\r')
                 for ln in fh]

    changes = 0

    # -------------------------------------------------
    # Pass 1: detriplify
    # -------------------------------------------------
    # Find call expressions whose argument is a long run
    # of adjacent string literals encoding preformatted
    # multi-line text, and convert each to a single
    # triple-quoted literal.  Process in reverse file
    # order so that line-number shifts from earlier
    # replacements do not invalidate later ranges.
    detriplify_source = '\n'.join(lines) + '\n'
    dcands = _find_detriplify_candidates(detriplify_source)
    for start, end in reversed(dcands):
        if start < 1 or end > len(lines) or start > end:
            continue
        old_block = [lines[j].rstrip()
                     for j in range(start - 1, end)]
        new_lines, ok = do_detriplify(lines, start, end)
        if not ok:
            continue
        if new_lines == old_block:
            continue
        print(f"  detriplify {start}-{end}: "
              f"{len(old_block)} lines -> "
              f"{len(new_lines)}")
        if not check_only:
            lines[start - 1:end] = new_lines
        changes += 1

    # -------------------------------------------------
    # Pass 2: unwrap / rewrap for over-wrapped blocks
    # -------------------------------------------------
    # Collect candidate blocks.  Each candidate is
    # (start, end) using 1-based inclusive line numbers.
    tq_lines = _build_triple_quote_map(lines)
    candidates = []
    n = len(lines)
    i = 0
    while i < n:
        line = lines[i]
        stripped = line.strip()

        # Skip blank lines and comment-only lines.
        if not stripped or stripped.startswith('#'):
            i += 1
            continue

        # Skip lines inside triple-quoted strings.
        if i in tq_lines:
            i += 1
            continue

        # Count bracket depth on this line.
        depth = 0
        in_str = False
        quote = None
        for ch in stripped:
            if in_str:
                if ch == '\\':
                    continue
                if ch == quote:
                    in_str = False
                continue
            if ch in ('"', "'"):
                in_str = True
                quote = ch
            elif ch in '([{':
                depth += 1
            elif ch in ')]}':
                depth -= 1

        if depth <= 0:
            i += 1
            continue

        # Line has unclosed brackets.  Look ahead up
        # to 10 lines for the matching close.
        end = -1
        scan_depth = depth
        for j in range(i + 1, min(i + 11, n)):
            jline = lines[j].strip()
            in_str = False
            quote = None
            for ch in jline:
                if in_str:
                    if ch == '\\':
                        continue
                    if ch == quote:
                        in_str = False
                    continue
                if ch in ('"', "'"):
                    in_str = True
                    quote = ch
                elif ch in '([{':
                    scan_depth += 1
                elif ch in ')]}':
                    scan_depth -= 1
            if scan_depth <= 0:
                end = j
                break

        if end >= 0 and (end - i + 1) >= 3:
            candidates.append((i + 1, end + 1))
            # Jump past this block to avoid overlapping
            # candidates.
            i = end + 1
        else:
            i += 1

    if not candidates:
        if changes == 0:
            print("  No over-wrapped blocks found.")
        return changes

    # Process candidates in reverse order so that line
    # number shifts from earlier replacements do not
    # affect later ones.  `changes` was already seeded by
    # the detriplify pass above, so we keep accumulating.
    candidates.reverse()

    for start, end in candidates:
        old_block = [lines[j].rstrip()
                     for j in range(start - 1, end)]

        # Try unwrap first (single line).
        result, fit = do_unwrap(lines, start, end)
        if fit:
            if [result] == old_block:
                continue
            print(f"  unwrap {start}-{end}: "
                  f"{len(old_block)} lines -> 1 "
                  f"({len(result)} chars)")
            if not check_only:
                lines[start - 1:end] = [result]
            changes += 1
            continue

        # Try rewrap (comma-split).
        new_lines, fit = do_rewrap(lines, start, end)
        if fit and new_lines != old_block:
            print(f"  rewrap {start}-{end}: "
                  f"{len(old_block)} lines -> "
                  f"{len(new_lines)}")
            if not check_only:
                lines[start - 1:end] = new_lines
            changes += 1
            continue

    if not check_only and changes > 0:
        with open(source_path, 'w') as fh:
            for ln in lines:
                fh.write(ln + '\n')

    return changes


def main():
    """Entry point: parse CLI arguments and dispatch."""
    args = sys.argv[1:]
    if not args or '-h' in args or '--help' in args:
        print(__doc__)
        sys.exit(0)

    source = args[0]
    check_only = '--check' in args
    noscan = '--noscan' in args

    if not os.path.isfile(source):
        print(f"Source not found: {source}",
              file=sys.stderr)
        sys.exit(1)

    if not noscan:
        # Default: scan mode.
        label = "Scanning" if check_only else "Fixing"
        print(f"{label} {source}:")
        n = scan_file(source, check_only)
        if check_only:
            print(f"\n{n} change(s) would be made.")
            sys.exit(1 if n > 0 else 0)
        else:
            print(f"\n{n} change(s) applied.")
        sys.exit(0)

    # Directive mode: need a directives file.
    remaining = [a for a in args if a != source
                 and not a.startswith('--')]
    if not remaining:
        print("Usage: rewrap_code.py SOURCE [--check]",
              file=sys.stderr)
        print("       rewrap_code.py SOURCE --noscan "
              "DIRECTIVES [--check]",
              file=sys.stderr)
        sys.exit(1)

    directive_file = remaining[0]
    if not os.path.isfile(directive_file):
        print(f"Directives not found: {directive_file}",
              file=sys.stderr)
        sys.exit(1)

    directives = parse_directives(directive_file)
    if not directives:
        print("No valid directives found.")
        sys.exit(0)

    label = "Would change" if check_only else "Applying"
    print(f"{label} {len(directives)} directive(s) "
          f"on {source}:")
    n = apply_directives(source, directives, check_only)

    if check_only:
        print(f"\n{n} change(s) would be made.")
        sys.exit(1 if n > 0 else 0)
    else:
        print(f"\n{n} change(s) applied.")


if __name__ == '__main__':
    main()
