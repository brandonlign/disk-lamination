"""The worked examples of Kang--Lee--Lim (arXiv:2508.04396v4, Section 5.1), recomputed.

Examples 3 and 4 print a c-polynomial for an arc in an octagon whose extended
exchange matrix is printed in Examples 2 and 4. Mutating that matrix along the
crossed diagonals 1, ..., 5 (and along every run of consecutive diagonals) gives
the polynomials below; the printed ones (2q^2+2q+2 and 2q^3+4q^2+2q+1) do not
occur, and the last check shows that they occur for no single lamination of
any triangulated octagon. Example 6 concerns the 12-gon of their Figure 17; with vertices labeled
A, ..., L clockwise from the leftmost one, the triangulation is
BL, CL, DL, DK, EK, EJ, EI, FI, GI, the arc is AH and the lamination runs from
edge AB to edge GH. All roots are checked to be simple and negative.
"""
from geometry import Seed
from polys import add, mul, normalize, simple_negative_roots


def mutate_runs(B, shear):
    m = len(B)
    out = {}
    for s in range(m):
        Bx = [row[:] for row in B] + [list(shear)]
        X = [(1,)] * m + [(0, 1)]
        for k in range(s, m):
            P = Q = (1,)
            for j in range(m + 1):
                for _ in range(max(0, Bx[j][k])):
                    P = mul(P, X[j])
                for _ in range(max(0, -Bx[j][k])):
                    Q = mul(Q, X[j])
            X[k] = add(P, Q)
            Bx = [[-Bx[a][b] if k in (a, b) else
                   Bx[a][b] + max(0, Bx[a][k]) * max(0, Bx[k][b]) - max(0, -Bx[a][k]) * max(0, -Bx[k][b])
                   for b in range(m)] for a in range(m + 1)]
            out[(s + 1, k + 1)] = X[k]
    return out


def show(p):
    return " + ".join(f"{c}" if i == 0 else f"{c}q" if i == 1 else f"{c}q^{i}" for i, c in enumerate(p) if c)


examples = {
    "Example 3 (matrix of Example 2)": ([[0, -1, 0, 0, 0], [1, 0, -1, 0, 0], [0, 1, 0, 1, 0], [0, 0, -1, 0, -1], [0, 0, 0, 1, 0]],
                                        [-1, 0, 1, -1, 1], (2, 2, 2)),
    "Example 4": ([[0, -1, 0, 0, 0], [1, 0, 1, 0, 0], [0, -1, 0, 1, 0], [0, 0, -1, 0, 1], [0, 0, 0, -1, 0]],
                  [-1, 1, 0, -1, 0], (1, 2, 4, 2)),
}
for name, (B, shear, printed) in examples.items():
    runs = mutate_runs(B, shear)
    assert all(simple_negative_roots(normalize(p)) for p in runs.values())
    assert printed not in {normalize(p) for p in runs.values()}
    print(f"{name}: arc crossing 1..5 gives {show(runs[(1, 5)])}; printed {show(printed)} does not occur in this seed")

# No single lamination of any triangulated octagon gives either printed polynomial.
from itertools import combinations
from geometry import triangulations
octagon = set()
for T in triangulations(tuple(range(8))):
    sd = Seed(8, T)
    for l, r in combinations(range(8), 2):
        if (r - l) % 8 not in (1, 7):
            sh = sd.shear(2 * l + 1, 2 * r + 1)
            octagon |= {normalize(sd.c_polynomial(sh, x, y)) for x, y in combinations(range(8), 2) if (y - x) % 8 not in (1, 7)}
assert all(printed not in octagon for _, _, printed in examples.values())
print(f"All single laminations of all triangulated octagons: {len(octagon)} distinct c-polynomials, neither printed one among them")

A, B_, C, D, E, F, G, H, I, J, K, L = [(-i) % 12 for i in range(12)]  # clockwise labels
seed = Seed(12, [(A, B_, L), (B_, C, L), (C, D, L), (D, K, L), (D, E, K), (E, J, K), (E, I, J), (E, F, I), (F, G, I), (G, H, I)])
edge_mid = lambda u, v: 2 * (u if (v - u) % 12 == 1 else v) + 1
shear = seed.shear(edge_mid(A, B_), edge_mid(G, H))
c = seed.c_polynomial(shear, min(A, H), max(A, H))
assert len(seed.crossed(min(A, H), max(A, H))) == 9 and simple_negative_roots(normalize(c))
print(f"Example 6 (Figure 17): arc AH gives {show(c)}")
