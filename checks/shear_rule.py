"""Check of the shear formula b_i = d_{i-1} - d_i of Section 3 in polygons with 4 to 10 vertices.

For every triangulation and single lamination L, the shear coordinates from the
S/Z rule (geometry.py) are compared with the turn-bit formula: list the
diagonals 1, ..., r crossed by L in order, let d_i (1 <= i < r) be 1 when the
common vertex of diagonals i and i+1 is to the right of L, and define d_0 and d_r
in the same way using the boundary edges where L starts and ends. Diagonals not
crossed by L must have shear coordinate 0. Expected output: one line per n.
"""
import sys
from itertools import combinations
from geometry import Seed, triangulations


def turn_bits(seed, l, r):
    n, N = seed.n, seed.N
    a, b = 2 * l + 1, 2 * r + 1
    right = lambda v: 0 < (2 * v - a) % N < (b - a) % N
    start = (l, (l + 1) % n)
    end = (r, (r + 1) % n)
    # Diagonals crossed by L, in order from its starting edge.
    first = next(t for t in seed.tri if set(start) <= set(t))
    rem = {e for e in seed.diag if right(e[0]) != right(e[1])}
    order, cur = [], first
    while rem:
        (c,) = [c for c in rem if set(c) <= set(cur)]
        order.append(c)
        rem.remove(c)
        cur = next(t for t in seed.tri if set(c) <= set(t) and t != cur)
    assert set(end) <= set(cur)
    sides = [start] + order + [end]
    d = []
    for e, f in zip(sides, sides[1:]):
        (v,) = set(e) & set(f)
        d.append(int(right(v)))
    return order, d


def run(n):
    count = 0
    for T in triangulations(tuple(range(n))):
        seed = Seed(n, T)
        for l, r in combinations(range(n), 2):
            if (r - l) % n in (1, n - 1):
                continue
            shear = seed.shear(2 * l + 1, 2 * r + 1)
            order, d = turn_bits(seed, l, r)
            expected = [0] * len(seed.diag)
            for i, e in enumerate(order):
                expected[seed.idx[e]] = d[i] - d[i + 1]
            assert shear == expected, (n, T, l, r, shear, expected)
            count += 1
    print(f"n={n}: {count} (triangulation, lamination) pairs; S/Z shear coordinates equal d_(i-1) - d_i", flush=True)


if __name__ == "__main__":
    for n in range(4, 1 + int(sys.argv[1] if len(sys.argv) > 1 else 10)):
        run(n)
