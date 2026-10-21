# Jacobian--Kummer class `0011`: exact \(\mathbb H^{10}\) Coxeter no-go

## Outcome

The subgroup class `class_0011` in the full published lattice/chamber
stabilizer of the 316-wall Jacobian--Kummer chamber does **not** produce a
hyperbolic Coxeter polytope.  It fixes a Lorentzian space of signature
\((1,10)\), so the dimension is genuinely \(10\), and the restricted chamber
has 76 exact irredundant facets.  However, 200 pairs of these genuine facets
meet with forbidden Coxeter products.

This is a metric no-go, not a statement about volume.  Vinberg's
finite-volume criterion applies to a Coxeter chamber only.  Since the exact
pair gate fails, no Vinberg volume computation is logically applicable.

No surface automorphism is assumed or required.  The involution is taken
from the full lattice/chamber stabilizer

\[
O(S_{16},D_{16})\cong (C_2)^5\rtimes S_6,
\]

of the published divisor-lattice chamber.

## Fixed slice and exact facet reduction

The representative of GAP subgroup class 11 is an involution with fixed
rank 11 in the rank-17 Néron--Severi lattice \(S_{16}\).  It preserves the
published Gram matrix, fixes the ample class, and permutes all 316 walls of
\(D_{16}\).  Its fixed wall span has exact inertia

\[
(1,10,167)
\]

on the 178 orbit-sum normals, hence projectivizes to \(\mathbb H^{10}\).
The wall action has 40 singleton orbits and 138 two-element orbits.

Exact conic reduction gives:

- 76 irredundant facets, each with an explicit strict timelike
  relative-interior witness; and
- 102 redundant orbit inequalities, each with an explicit positive rational
  conic dependence on the other orbit normals.

Thus a forbidden product below cannot be dismissed as an artifact of a
redundant wall.  The retained facet-normal squares are

| square | facets |
|---:|---:|
| \(-1\) | 32 |
| \(-2\) | 16 |
| \(-3/4\) | 16 |
| \(-4\) | 12 |

## Complete exact Coxeter-pair gate

For spacelike facet normals \(r_i,r_j\), the normalized squared product is

\[
q_{ij}=\frac{(r_i,r_j)^2}{(r_i,r_i)(r_j,r_j)}.
\]

When \(0<q_{ij}<1\), a Coxeter angle would require
\(q_{ij}=\cos^2(\pi/m)\) for an integer \(m\geq3\).  Here all products are
rational, for which the allowed values are \(1/4,1/2,3/4\); \(q=1\) is a
parallel pair and \(q>1\) an ultraparallel pair.

All \(\binom{76}{2}=2850\) pairs were classified exactly:

| relation | pairs |
|---|---:|
| orthogonal, \(m=2\) | 1,068 |
| \(m=3\) | 32 |
| \(m=4\) | 184 |
| parallel | 406 |
| ultraparallel | 960 |
| **forbidden** | **200** |

The forbidden-product histogram is:

| normalized squared product | pairs |
|---:|---:|
| \(1/9\) | 24 |
| \(1/3\) | 112 |
| \(2/3\) | 64 |

None equals \(\cos^2(\pi/m)\) for an integer \(m\).  Consequently the exact
restricted chamber is not a Coxeter chamber, and `class_0011` is eliminated
from the higher-dimensional example search.

## Reproduction

- Exact certificate:
  [`jacobian_kummer_h10_class0011_wall_certificate.json`](../data/enumerations/jacobian_kummer_h10_class0011_wall_certificate.json)
- Builder and archive checker:
  [`screen_jacobian_kummer_h10_class0011.py`](../scripts/screen_jacobian_kummer_h10_class0011.py)
- Regression tests:
  [`test_jacobian_kummer_h10_class0011.py`](../tests/test_jacobian_kummer_h10_class0011.py)

Replay with:

```bash
PYTHONPATH=src:. python3 -m scripts.screen_jacobian_kummer_h10_class0011 --check-archive
PYTHONPATH=src:. pytest -q tests/test_jacobian_kummer_h10_class0011.py
```

Primary source: I. Dolgachev and I. Shimada, *15-nodal quartic surfaces.
Part II: the automorphism group*,
[DOI 10.1007/s12215-019-00464-7](https://doi.org/10.1007/s12215-019-00464-7),
[arXiv:1908.05390](https://arxiv.org/abs/1908.05390).
