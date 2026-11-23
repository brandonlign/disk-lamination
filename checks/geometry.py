"""Triangulated polygons, single laminations and cluster mutation, following the
conventions of Kang--Lee--Lim (arXiv:2508.04396v4, Section 2).

Vertices 0, ..., n-1 sit counterclockwise on a circle at angular positions
0, 2, ..., 2n-2 out of 2n; the odd position 2i+1 is the midpoint of the boundary
edge from i to i+1, where a lamination may end. Arcs and diagonals are sorted
vertex pairs.
"""
from functools import lru_cache
from itertools import combinations
from polys import add, mul


@lru_cache(None)
def triangulations(vs):
    """All triangulations of the convex polygon with vertex tuple vs, as tuples of triangles."""
    if len(vs) < 3:
        return ((),)
    out = []
    for k in range(1, len(vs) - 1):
        for left in triangulations(vs[:k + 1]):
            for right in triangulations(vs[k:]):
                out.append(left + right + ((vs[0], vs[k], vs[-1]),))
    return tuple(out)


def crosses(p1, p2, a, b, N):
    """Whether chords p1p2 and ab (distinct positions mod N) cross."""
    inside = lambda x: 0 < (x - a) % N < (b - a) % N
    return inside(p1) != inside(p2)


class Seed:
    def __init__(self, n, triangles):
        self.n, self.N = n, 2 * n
        self.tri = [tuple(sorted(t)) for t in triangles]
        sides = {e for t in self.tri for e in combinations(t, 2)}
        self.diag = sorted(e for e in sides if (e[1] - e[0]) % n not in (1, n - 1))
        self.idx = {e: i for i, e in enumerate(self.diag)}
        m = len(self.diag)
        # b_ij = 1 if diagonal j immediately follows diagonal i clockwise in a triangle.
        self.B = [[0] * m for _ in range(m)]
        for x, y, z in self.tri:
            ccw = [(x, y), (y, z), (x, z)]
            for i in range(3):
                e, f = ccw[i], ccw[(i + 1) % 3]  # f follows e counterclockwise
                if e in self.idx and f in self.idx:
                    self.B[self.idx[f]][self.idx[e]] += 1
                    self.B[self.idx[e]][self.idx[f]] -= 1

    def apexes(self, e):
        return [v for t in self.tri if set(e) <= set(t) for v in t if v not in e]

    def shear(self, la, lb):
        """Shear coordinates of the lamination with endpoints at odd positions la, lb,
        by the S/Z rule of Kang--Lee--Lim, Definition 11 and Figure 3."""
        n, N = self.n, self.N

        def meets(u, v):
            for p in (la, lb):
                if (v - u) % n == 1 and p == (2 * u + 1) % N or (u - v) % n == 1 and p == (2 * v + 1) % N:
                    return True
            return crosses(2 * u, 2 * v, la, lb, N)

        out = [0] * len(self.diag)
        for e in self.diag:
            u, v = e
            if not crosses(2 * u, 2 * v, la, lb, N):
                continue
            X = next(p for p in self.apexes(e) if 0 < (p - u) % n < (v - u) % n)
            Y = next(p for p in self.apexes(e) if p != X)
            # Counterclockwise the quadrilateral is u, X, v, Y.
            hit = {s for s in [(u, X), (X, v), (v, Y), (Y, u)] if meets(*s)}
            assert len(hit) == 2
            if hit == {(X, v), (Y, u)}:
                out[self.idx[e]] = 1
            elif hit == {(u, X), (v, Y)}:
                out[self.idx[e]] = -1
        return out

    def crossed(self, x, y):
        """Diagonals crossed by the arc xy, in order from x."""
        N = self.N
        rem = {e for e in self.diag if x not in e and y not in e and crosses(2 * e[0], 2 * e[1], 2 * x, 2 * y, N)}
        order = []
        cur = next((t for t in self.tri if x in t and any(set(c) <= set(t) for c in rem)), None)
        while rem:
            (c,) = [c for c in rem if set(c) <= set(cur)]
            order.append(c)
            rem.remove(c)
            cur = next(t for t in self.tri if set(c) <= set(t) and t != cur)
        return order

    def c_polynomial(self, shear, x, y):
        """Cluster variable of the arc xy with initial variables 1 and the lamination
        variable q, by mutating along the crossed diagonals (as a coefficient tuple)."""
        order = self.crossed(x, y)
        m = len(self.diag)
        B = [row[:] for row in self.B] + [list(shear)]
        X = [(1,)] * m + [(0, 1)]
        for c in order:
            k = self.idx[c]
            assert X[k] == (1,)  # each diagonal is still initial when it is flipped
            P = Q = (1,)
            for j in range(m + 1):
                for _ in range(max(0, B[j][k])):
                    P = mul(P, X[j])
                for _ in range(max(0, -B[j][k])):
                    Q = mul(Q, X[j])
            X[k] = add(P, Q)
            B = [[-B[a][b] if k in (a, b) else
                  B[a][b] + max(0, B[a][k]) * max(0, B[k][b]) - max(0, -B[a][k]) * max(0, -B[k][b])
                  for b in range(m)] for a in range(m + 1)]
        return X[self.idx[order[-1]]] if order else (1,)
