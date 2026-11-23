"""Exhaustive check of the main theorem in polygons with 4 to n vertices (n = 9 in run_all.sh).

For every triangulation of the n-gon, every single lamination (endpoints on two
non-adjacent boundary edges) and every arc, the c-polynomial is computed twice:
by cluster mutation with the shear coordinates of the S/Z rule (geometry.py),
and by the weighted order-ideal sum of Section 3, sum_I q^(sum_i b_i x_i) divided
by its lowest power of q, where a change 0 -> 1 between consecutive crossed
diagonals is allowed when their common vertex lies to the right of the arc.
The two must agree, and each distinct polynomial is checked by Sturm's theorem
to have only simple negative roots. Expected output: one line per n.
"""
import sys
from collections import Counter
from itertools import combinations
from geometry import Seed, triangulations, crosses
from polys import normalize, simple_negative_roots


def ideal_sum(seed, order, x, y, shear):
    """sum over binary x_1..x_r obeying the fence rule of q^(sum b_i x_i), normalized."""
    if not order:
        return (1,)
    right = lambda v: 0 < (2 * v - 2 * x) % seed.N < (2 * y - 2 * x) % seed.N
    b = [shear[seed.idx[e]] for e in order]
    states = [Counter({0: 1}), Counter({b[0]: 1})]
    for i in range(1, len(order)):
        (v,) = set(order[i - 1]) & set(order[i])
        allowed_up = right(v)  # 0 -> 1 allowed if True, 1 -> 0 allowed otherwise
        new = [Counter(), Counter()]
        for s in (0, 1):
            for t in (0, 1):
                if s != t and (t == 1) != allowed_up:
                    continue
                for k, c in states[s].items():
                    new[t][k + t * b[i]] += c
        states = new
    total = states[0] + states[1]
    lo, hi = min(total), max(total)
    return normalize(tuple(total.get(k, 0) for k in range(lo, hi + 1)))


def run(n, seen):
    count = 0
    for T in triangulations(tuple(range(n))):
        seed = Seed(n, T)
        for l, r in combinations(range(n), 2):
            if (r - l) % n in (1, n - 1):
                continue
            shear = seed.shear(2 * l + 1, 2 * r + 1)
            for x, y in combinations(range(n), 2):
                if (y - x) % n in (1, n - 1):
                    continue
                c = seed.c_polynomial(shear, x, y)
                assert c == ideal_sum(seed, seed.crossed(x, y), x, y, shear), (n, T, l, r, x, y)
                if c not in seen:
                    assert simple_negative_roots(normalize(c)), c
                    seen.add(c)
                count += 1
    print(f"n={n}: {count} (triangulation, lamination, arc) triples; mutation equals the ideal sum; "
          f"{len(seen)} distinct polynomials so far, all with simple negative roots", flush=True)


if __name__ == "__main__":
    seen = set()
    for n in range(4, 1 + int(sys.argv[1] if len(sys.argv) > 1 else 8)):
        run(n, seen)
