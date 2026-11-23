"""Check of Lemma 2.1 on random instances.

For random step lists (each step allows staying or one change in a given
direction, with a random positive integer weight), the walk enumerator is
computed by dynamic programming and compared with the determinant formula of
the proof: det(I + q A U B U^T) for equal end states, where A = diag(a_i),
B = diag(b_i), U is upper triangular of ones, and the coefficient of t for
distinct end states. Every normalized enumerator is then checked to have only
simple negative roots. Expected output: a count of instances and "ok".
"""
import random
from fractions import Fraction
from polys import add, mul, normalize, simple_negative_roots, trim


def walk_enumerator(steps, start, end):
    """steps: list of (direction, weight). A change in direction 1 is 0 -> 1."""
    state = {start: (1,), 1 - start: (0,)}
    for direction, w in steps:
        new = {0: state[0], 1: state[1]}
        src, dst = 1 - direction, direction
        moved = mul(state[src], (0, w) if direction == 1 else (w,))
        new[dst] = add(new[dst], moved)
        state = new
    return state[end]


def det_poly(M):
    """Determinant of a square matrix whose entries are polynomials (tuples)."""
    n = len(M)
    if n == 0:
        return (1,)
    total = (0,)
    for j in range(n):
        if M[0][j] == (0,):
            continue
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        term = mul(M[0][j], det_poly(minor))
        total = add(total, term if j % 2 == 0 else tuple(-c for c in term))
    return total


def gram_poly(a, b):
    """det(I + q A U B U^T) as a polynomial in q."""
    s = len(a)
    G = [[a[i] * sum(b[k] for k in range(max(i, j), s)) for j in range(s)] for i in range(s)]
    return det_poly([[add((1,) if i == j else (0,), (0, G[i][j])) for j in range(s)] for i in range(s)])


def merged_runs(steps, start, end):
    """Merge equal-direction runs and prune unusable first and last runs, as in the proof."""
    runs = []
    for d, w in steps:
        if runs and runs[-1][0] == d:
            runs[-1][1] += w
        else:
            runs.append([d, w])
    while runs and runs[0][0] == start:
        runs.pop(0)
    while runs and runs[-1][0] != end:
        runs.pop()
    return [w for _, w in runs]


def check(steps, start, end):
    enum = walk_enumerator(steps, start, end)
    if enum == (0,):
        return False
    w = merged_runs(steps, start, end)
    if start == end:
        assert len(w) % 2 == 0
        a, b = w[0::2], w[1::2]
        assert normalize(gram_poly(a, b)) == normalize(enum)
    else:
        assert len(w) % 2 == 1
        a, b = w[0::2], w[1::2]
        # The determinant is affine in the last weight t; its t-coefficient is P(t=1) - P(t=0).
        with_t = gram_poly(a, b + [1])
        without_t = gram_poly(a, b + [0])
        coeff = add(with_t, tuple(-c for c in without_t))
        assert normalize(coeff) == normalize(enum)
    assert simple_negative_roots(normalize(enum))
    return True


random.seed(1)
count = 0
for _ in range(3000):
    length = random.randint(1, 14)
    steps = [(random.randint(0, 1), random.randint(1, 5)) for _ in range(length)]
    for start in (0, 1):
        for end in (0, 1):
            count += check(steps, start, end)
print(f"Lemma 2.1: {count} random walk enumerators agree with the determinant formula and have simple negative roots: ok")
