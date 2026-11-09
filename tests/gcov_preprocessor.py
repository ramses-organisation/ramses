"""Work out which source lines no build of the test suite ever compiled.

gcov marks a line "-" both when it is genuinely not executable (a comment, a
declaration, the continuation of a statement counted on its first line) and
when the compiler never saw it because it sat in a preprocessor branch that
was switched off.  Those two cases mean completely different things for
coverage: the first is not a gap, the second is a gap no namelist can close.

Given the -D flags of every build in the run, this module walks the
#if/#ifdef/#else/#endif structure of a source file, following the #define and
#undef it contains, and returns the lines that are dead in *all* builds,
together with the condition that gated each one.
"""
import re
import warnings

COND = re.compile(r'^\s*#\s*(ifdef|ifndef|if|elif|else|endif)\b(.*)$')
# a "(" right after the name makes it a function-like macro
DEFINE = re.compile(r'^\s*#\s*(define|undef)\s+(\w+)(\()?\s*(.*)$')
IDENT = re.compile(r'\b[A-Za-z_]\w*\b')


def canonical(expr):
    """One spelling per condition, so `NVAR > NHYDRO` and `NVAR>NHYDRO` agree."""
    e = re.sub(r'\s+', '', str(expr).strip())
    e = re.sub(r'^\((.*)\)$', r'\1', e) if e.count('(') == e.count(')') and \
        e.startswith('(') and e.endswith(')') and '(' not in e[1:-1] else e
    return e


def label_of(expr):
    """Render a canonical condition the way it is written in the source."""
    e = canonical(expr)
    m = re.fullmatch(r'defined\((\w+)\)', e)
    if m:
        return f"#ifdef {m.group(1)}"
    m = re.fullmatch(r'!defined\((\w+)\)', e)
    if m:
        return f"#ifndef {m.group(1)}"
    # canonical form has no spaces; put them back around the operators only
    e = e.replace("&&", " && ").replace("||", " || ")
    return f"#if {e}"


def constant(expr):
    """Evaluate a constant expression. Python's SyntaxWarnings about
    something like `0(1)` are silenced: the caller handles the failure."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return eval(expr, {"__builtins__": {}}, {})


def macro_value(v):
    """'2+3' -> 5, '' -> 1 (defined without a value), anything else as text."""
    if not v.strip():
        return 1
    try:
        return int(constant(v))
    except Exception:
        return v


def parse_defines(defines):
    """'-DNDIM=3 -DSOLVERmhd' -> {'NDIM': 3, 'SOLVERmhd': 1}"""
    out = {}
    for tok in str(defines).split():
        if not tok.startswith("-D"):
            continue
        k, _, v = tok[2:].partition("=")
        out[k] = macro_value(v)
    return out


def eval_cond(expr, macros):
    """Evaluate a C preprocessor condition. None when it cannot be decided."""
    e = expr.split("//")[0].strip()
    if not e:
        return None
    e = re.sub(r'defined\s*\(\s*(\w+)\s*\)',
               lambda m: '1' if m.group(1) in macros else '0', e)
    e = re.sub(r'defined\s+(\w+)',
               lambda m: '1' if m.group(1) in macros else '0', e)

    def sub(m):
        n = m.group(0)
        v = macros.get(n)
        if v is None:
            return '0'                     # undefined macro is 0 in #if
        return str(v) if isinstance(v, int) else '1'
    e = IDENT.sub(sub, e)
    e = e.replace('!=', ' __NE__ ').replace('&&', ' and ').replace('||', ' or ')
    e = re.sub(r'!(?=\s*[\d(])', ' not ', e)
    e = e.replace('__NE__', '!=')
    try:
        return bool(constant(e))
    except Exception:
        return None


def negate(expr):
    """!defined(X) <-> defined(X); anything else gets wrapped."""
    e = canonical(expr)
    if e.startswith("!") and "&&" not in e and "||" not in e:
        return e[1:]
    if re.fullmatch(r'defined\(\w+\)', e):
        return "!" + e
    return f"!({e})"


def unbuilt_lines(source_lines, macro_sets):
    """
    Returns ({line number: the directive that kept it out of every build},
    [(line number, directive) for each condition that could not be evaluated
    for a build that reaches it]). An undecided condition is taken as true,
    so what it gates counts as compiled.
    """
    # live[i] is True where build i currently compiles
    n = len(macro_sets)
    # each build's macros, as the #define and #undef it compiles change them
    macros = [dict(ms) for ms in macro_sets]
    stack = []            # (live_now[], taken_before[], label)
    dead = {}
    undecided = []
    for ln, raw in enumerate(source_lines, 1):
        m = COND.match(raw)
        if not m:
            d = DEFINE.match(raw)
            if d:
                live_now = stack[-1][0] if stack else [True] * n
                for ms, on in zip(macros, live_now):
                    if not on:
                        continue
                    if d.group(1) == "undef":
                        ms.pop(d.group(2), None)
                    else:   # a function-like macro only counts as defined
                        ms[d.group(2)] = 1 if d.group(3) else macro_value(d.group(4))
            if stack and not any(stack[-1][0]):
                dead[ln] = stack[-1][2]
            continue
        kw, rest = m.group(1), m.group(2).strip()
        if kw in ("ifdef", "ifndef", "if"):
            if kw == "ifdef":
                cond = [rest.split()[0] in ms if rest else None for ms in macros]
                expr = f"defined({rest.split()[0]})" if rest else "1"
            elif kw == "ifndef":
                cond = [rest.split()[0] not in ms if rest else None for ms in macros]
                expr = f"!defined({rest.split()[0]})" if rest else "1"
            else:
                cond = [eval_cond(rest, ms) for ms in macros]
                expr = canonical(rest)
            outer = stack[-1][0] if stack else [True] * n
            if any(o and c is None for o, c in zip(outer, cond)):
                undecided.append((ln, raw.strip()))
            live = [bool(o) and (c is not False) for o, c in zip(outer, cond)]
            label = expr
            # if the enclosing block is already dead everywhere, that is the
            # real reason this code is never compiled -- keep the outer label
            if stack and not any(stack[-1][0]):
                label = stack[-1][2]
            stack.append([live, list(live), label, [expr]])
        elif kw == "elif":
            if not stack:
                continue
            _, taken, _, seen = stack[-1]
            cond = [eval_cond(rest, ms) for ms in macros]
            outer = stack[-2][0] if len(stack) > 1 else [True] * n
            if any(o and not t and c is None for o, t, c in zip(outer, taken, cond)):
                undecided.append((ln, raw.strip()))
            live = [bool(o) and not t and (c is not False)
                    for o, t, c in zip(outer, taken, cond)]
            expr = canonical(rest)
            # `#elif X` gates exactly the same code as `#if X` would
            label = expr
            if len(stack) > 1 and not any(stack[-2][0]):
                label = stack[-2][2]
            stack[-1] = [live, [t or l for t, l in zip(taken, live)],
                         label, seen + [expr]]
        elif kw == "else":
            if not stack:
                continue
            _, taken, _, seen = stack[-1]
            outer = stack[-2][0] if len(stack) > 1 else [True] * n
            live = [bool(o) and not t for o, t in zip(outer, taken)]
            # the else branch is gated by the negation of everything above it
            joined = "||".join(seen)
            expr = f"!({joined})" if len(seen) > 1 else negate(seen[0])
            label = expr
            if len(stack) > 1 and not any(stack[-2][0]):
                label = stack[-2][2]
            stack[-1] = [live, [True] * n, label, seen]
        elif kw == "endif":
            if stack:
                stack.pop()
    return dead, undecided
