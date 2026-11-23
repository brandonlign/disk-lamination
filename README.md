# Single-lamination cluster polynomials of a disk are real-rooted

Brandon Li

Paper and verification code proving Conjecture 1 of Kang, Lee and Lim
([arXiv:2508.04396](https://arxiv.org/abs/2508.04396), European J. Combin. 136
(2026) 104392). For a triangulated polygon with coefficients given by a single
lamination $L$, their $c$-polynomial $c_\gamma(q)$ of an arc $\gamma$ is the
cluster variable $x_\gamma$ with every initial cluster variable set to $1$ and
the coefficient variable set to $q$. They conjectured that it is log-concave.

**Theorem.** For every triangulation of an unpunctured disk, every single
lamination and every arc $\gamma$, the polynomial $c_\gamma(q)$ is a positive
constant or has only simple, negative real roots.

In particular $c_\gamma$ is log-concave and unimodal, and the specialization of
every cluster monomial has only negative roots. The proof writes $c_\gamma$ as
a weighted count of two-state walks along the diagonals that $\gamma$ and $L$
both cross, and shows that such a count equals $\det(I+qMM^{\mathsf T})$ for a
matrix $M$ with bidiagonal inverse. The paper also recomputes the worked
examples of Kang, Lee and Lim.

- [Paper (PDF)](paper/main.pdf)
- [LaTeX source](paper/)
- [Verification code](checks/)

## Reproduce

```sh
bash checks/run_all.sh
```

This needs Python 3 only (standard library, exact arithmetic) and takes about
fifteen seconds. It runs:

- `path_lemma.py`: Lemma 2.1, comparing walk counts with the determinant
  formula on random step lists and checking the roots.
- `shear_rule.py`: Lemma 3.1, comparing shear coordinates from the S/Z rule
  with the telescoping formula $d_{i-1}-d_i$ in all triangulated $n$-gons,
  $n\le10$.
- `disk_polygons.py`: for every triangulation, single lamination and arc of an
  $n$-gon with $n\le9$, mutation equals the ideal sum of Section 3, and every
  polynomial has simple negative roots (Sturm's theorem).
- `kll_examples.py`: the examples of Section 4.

`geometry.py` (triangulations, shear coordinates, mutation) and `polys.py`
(polynomial arithmetic, Sturm sequences) are shared helpers. The output is in
[`checks/output.txt`](checks/output.txt). The PDF builds with Tectonic
(`cd paper && tectonic main.tex`) or `pdflatex`.

## License

Code: MIT. Paper: CC BY 4.0.
