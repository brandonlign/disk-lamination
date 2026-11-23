"""Exact integer polynomials as coefficient tuples (index = power of q), and a
Sturm-sequence test that a polynomial has only simple negative real roots."""
from fractions import Fraction


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)


def add(p, r):
    out = [0] * max(len(p), len(r))
    for i, c in enumerate(p):
        out[i] += c
    for i, c in enumerate(r):
        out[i] += c
    return trim(out)


def mul(p, r):
    out = [0] * (len(p) + len(r) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(r):
            out[i + j] += a * b
    return trim(out)


def normalize(p):
    """Divide by the largest power of q dividing p."""
    p = trim(p)
    k = next(i for i, c in enumerate(p) if c != 0)
    return p[k:]


def _rem(a, b):
    a = [Fraction(c) for c in a]
    while len(a) >= len(b) and any(a):
        f = a[-1] / b[-1]
        s = len(a) - len(b)
        for i, c in enumerate(b):
            a[s + i] -= f * c
        a.pop()
    return trim(a) if a else (Fraction(0),)


def _sign_changes(values):
    signs = [v > 0 for v in values if v != 0]
    return sum(1 for x, y in zip(signs, signs[1:]) if x != y)


def simple_negative_roots(p):
    """True if p (with p(0) != 0) has deg p distinct roots, all real and negative."""
    p = trim(p)
    assert p[0] != 0
    d = len(p) - 1
    if d == 0:
        return True
    seq = [tuple(Fraction(c) for c in p), tuple(Fraction(i * c) for i, c in enumerate(p) if i)]
    while len(seq[-1]) > 1 or seq[-1][0] != 0:
        r = _rem(seq[-2], seq[-1])
        if r == (0,):
            break
        seq.append(tuple(-c for c in r))
    at_minus_inf = [s[-1] * (-1) ** (len(s) - 1) for s in seq]
    at_zero = [s[0] for s in seq]
    return _sign_changes(at_minus_inf) - _sign_changes(at_zero) == d
