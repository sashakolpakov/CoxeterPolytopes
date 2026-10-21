# A 17-wall finite-volume Coxeter polytope in dimension 12

## Lineage correction (2026-08-06)

The polytope below is a valid finite-volume H12 Coxeter chamber, but it is
**not a new abstract chamber**.  An exact audit identifies its lattice with
entry 6 of Scharlau--Kirschmer's `List12`: the determinant and genus agree,
an explicit unimodular isometry identifies the integral Gram matrices, and
all twelve proper-face counts plus all three cusps agree.  An exhaustive
all-reflective-root Vinberg replay further proves that the 17 displayed roots
generate the full reflection group of that lattice.

The archive does not distribute its chamber roots, so no source-facet
permutation can be displayed.  The strongest proved statement is chamber
identity up to conjugacy by the full reflection group.  See the
[`exact lineage audit`](k3_rank15_order3_h12_scharlau_kirschmer_lineage.md).

## Certified polytope

The order-three symmetry

$$
(14\ 15\ 16)(17\ 18\ 19)
$$

of the complete weighted $(-2)$-curve graph for

$$
\operatorname{NS}(X)=
U\oplus E_8(-1)\oplus D_4(-1)\oplus A_1(-1)
$$

produces a finite-volume, noncompact Coxeter polytope in $\mathbb H^{12}$.
Its 17 wall orbits are

$$
1,2,3,4,5,6,7,8,9,10,11,12,13,
(14\ 15\ 16),(17\ 18\ 19),20,21.
$$

The orbit Gram has rank 13 and inertia $(1,12,4)$.  Every orbit normal is
spacelike; the two three-cycle normals have norm $-6$, and all other normals
have norm $-2$.

The exact polytope data are:

| invariant | value |
|---|---:|
| facets | 17 |
| finite vertices | 93 |
| ideal vertices / cusps | 3 |
| hyperideal extreme rays | 0 |
| elliptic rank-11 edge diagrams | 623 |

The Coxeter relations are

| relation | count |
|---|---:|
| right angle ($m=2$) | 116 |
| $m=3$ | 12 |
| $m=6$ | 2 |
| parallel | 3 |
| ultraparallel | 3 |

The three ultraparallel normalized squared products are $3,3,9$.

## Two exact volume audits

The optimized Vinberg implementation and a deliberately literal
CoxIter-Proposition-1.10 implementation were both run on the complete
17-wall diagram.

The literal implementation enumerates:

- every rank-11 spherical principal subdiagram;
- every rank-12 spherical principal subdiagram; and
- every principal subset of order 12 through 17, retaining exactly the
  componentwise parabolic diagrams of Coxeter rank 11.

It obtains 623 edges, 93 finite vertices, and 3 parabolic vertices.  Every
edge has exactly two endpoints, with endpoint-type histogram

$$
(2\text{ finite},0\text{ ideal})^{496},\qquad
(1,1)^{124},\qquad
(0,2)^3.
$$

There are no bad edges.  The literal finite and ideal endpoint sets equal
the optimized Vinberg endpoint sets, and both equal the exact extreme-ray
sets.  Every full parabolic zero-neighbourhood passes Vinberg's rank test.

An independent CoxIter run on the archived 17-node input gives

$$
(f_0,\ldots,f_{12})=
(96,623,2054,4433,6807,7680,6444,4017,1835,596,130,17,1).
$$

It again reports finite covolume, noncompactness, and three vertices at
infinity.  Its exact Euler characteristic and covolume are

$$
\chi=\frac{691}{16554295296000},\qquad
\operatorname{covol}=\pi^6\frac{691}{2688779681280000}.
$$

A third certifier, sharing neither the production fold code nor its volume
implementation, reconstructs the cone from the 21-curve Gram matrix and
checks all $2^{17}$ facet intersections.  It finds 34,733 nonempty faces and
the same full $f$-vector.  It also exhausts the weighted graph group
$S_3\times C_2$: the two inverse order-three elements give the same H12 wall
partition, whereas the two order-six elements have only nine orbits and
include an isotropic and a timelike orbit sum, so they are correctly rejected.

## Integral reflection structure

The graph action reconstructs in the archived curve basis as a unimodular
integral order-three isometry of $\operatorname{NS}(X)$.  Its characteristic
polynomial is

$$
(t-1)^{13}(t^2+t+1),
$$

so its fixed and moving ranks are 13 and 2.

For an orbit of $k$ pairwise-disjoint $(-2)$-roots, the $k$ root reflections
commute.  On the invariant subspace their product is exactly reflection in
the orbit sum, whose norm is $-2k$.  This proves integrality for both
three-cycle walls (norm $-6$), as well as for every fixed-root wall.

## Fixed lattice and literature identity

The order-three graph symmetry is the triality of the $D_4(-1)$ summand.
Its fixed sublattice in $D_4(-1)$ is $A_2(-1)$: if $c$ is the central
$D_4$ root and $a$ is the sum of the three outer roots, then

$$
\begin{pmatrix}c^2&c\cdot a\\c\cdot a&a^2\end{pmatrix}
=\begin{pmatrix}-2&3\\3&-6\end{pmatrix},
$$

and the basis $(c,c+a)$ has the $A_2(-1)$ Gram matrix.  Hence the primitive
fixed Néron--Severi lattice is

$$
U\oplus E_8(-1)\oplus A_2(-1)\oplus A_1(-1),
$$

of signature $(1,12)$ and absolute determinant $6$.  Thus this construction
is the standard $D_4\rightsquigarrow G_2$ triality folding, now used as a
complete finite-volume Coxeter construction.  Its two long walls are the
norm-$6$ orbit sums.

Broad reflectivity of this fixed lattice was already known.  Wang's second
quasi-pullback construction explicitly includes

$$
K=E_8\oplus A_1\oplus A_2
$$

among the 33 positive-definite lattices for which the pullback of the
Borcherds form $\Phi_{12}$ is reflective on $2U\oplus K(-1)$.  Its positive
Weyl vector therefore implies that $U\oplus K(-1)$ is hyperbolic reflective.
On the other hand, Wang's complete hyperbolic $2$-reflective classification
does **not** contain this rank-13 lattice: the norm-$6$ triality reflections
are essential.  See H. Wang, *The classification of 2-reflective modular
forms*, Sections 5.2 and 7, doi:10.4171/JEMS/1358, together with Borcherds'
positive-Weyl-vector criterion.

This does not coincide with the obvious full nef chambers in Roulleau's
rank-13 table: in particular, the published 17-curve case
$U\oplus E_8(-1)\oplus 3A_1(-1)$ has determinant $8$ and only $(-2)$-curve
walls.  It does, however, coincide with the determinant-six entry of the
Scharlau--Kirschmer reflective-lattice census.  Their `List12` row records
the identical complete face data and three cusps; their corresponding
`Gram12` row is explicitly integrally isometric to the negative of the fixed
lattice Gram used here.  The independent all-root audit proves that these 17
walls are the full reflection chamber, not just an unidentified finite-index
reflection subgroup.  Thus no new abstract-polytope claim remains.  The
divisor-lattice/triality derivation may still be a distinct algebro-geometric
realization, but its novelty is not established.

No realization of the graph symmetry by an automorphism of a particular
surface is asserted or needed.

## Reproducibility

- certificate:
  [`k3_rank15_order3_h12_fold.json`](../data/enumerations/k3_rank15_order3_h12_fold.json)
- CoxIter input:
  [`k3_rank15_order3_h12_fold.coxiter`](../data/enumerations/k3_rank15_order3_h12_fold.coxiter)
- independent CoxIter audit:
  [`k3_rank15_order3_h12_coxiter_audit.json`](../data/enumerations/k3_rank15_order3_h12_coxiter_audit.json)
- independent exact audit:
  [`k3_rho15_h12_cyclic_fold_audit.json`](../data/candidates/k3_rho15_h12_cyclic_fold_audit.json)
- exact Scharlau--Kirschmer lineage and all-root audit:
  [`k3_rank15_order3_h12_scharlau_kirschmer_audit.json`](../data/enumerations/k3_rank15_order3_h12_scharlau_kirschmer_audit.json)
- verifier:
  [`certify_k3_rank15_order3_fold.py`](../scripts/certify_k3_rank15_order3_fold.py)
- independent verifier:
  [`independent_certify_k3_rho15_h12_cyclic_fold.py`](../scripts/independent_certify_k3_rho15_h12_cyclic_fold.py)
- lineage/all-root verifier:
  [`audit_k3_rank15_order3_h12_scharlau_kirschmer.py`](../scripts/audit_k3_rank15_order3_h12_scharlau_kirschmer.py)
- tests:
  [`test_k3_rank15_order3_fold.py`](../tests/test_k3_rank15_order3_fold.py),
  [`test_k3_rho15_h12_cyclic_fold_audit.py`](../tests/test_k3_rho15_h12_cyclic_fold_audit.py)

Full literal replay:

```bash
PYTHONPATH=src:. python3 scripts/certify_k3_rank15_order3_fold.py \
  --check-archive
PYTHONPATH=src:. python3 \
  scripts/independent_certify_k3_rho15_h12_cyclic_fold.py \
  --verify-certificate
```

Source: S. Kondō, *Algebraic K3 surfaces with finite automorphism groups*,
Section 4.5 and Figure 5, doi:10.1017/S0027763000001653; X. Roulleau,
*An atlas of K3 surfaces with finite automorphism group*, Section 15.1;
H. Wang, *The classification of 2-reflective modular forms*, Section 5.2 and
Remark 7.4, doi:10.4171/JEMS/1358.
