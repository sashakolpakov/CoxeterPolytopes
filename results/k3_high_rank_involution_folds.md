# Two exact finite-volume K3 wall-graph folds in dimensions 13 and 14

## Result

There are two certified positive constructions, up to symmetry of the source
wall graph:

| source Néron--Severi lattice | graph involution class | fixed hyperbolic space | facets | finite vertices | cusps | Vinberg edges |
|---|---:|---:|---:|---:|---:|---:|
| $U\oplus E_8(-1)\oplus D_4(-1)\oplus A_1(-1)$ | three conjugate involutions | $\mathbb H^{13}$ | 19 | 178 | 7 | 1323 |
| $U\oplus E_8(-1)\oplus D_6(-1)$ | one involution | $\mathbb H^{14}$ | 17 | 89 | 5 | 704 |

Both chambers are **finite-volume and noncompact**.  Every listed Vinberg
edge has exactly two finite-or-ideal endpoint extensions, every full
parabolic zero-neighbourhood has rank one less than the hyperbolic
dimension, and an independent exact extreme-ray enumeration finds no
hyperideal ray.

These are genuine hyperbolic Coxeter polytopes constructed from complete K3
divisor wall configurations.  They are not, however, fixed-nef chambers of
the published automorphism actions on the corresponding general K3 surfaces.
They use exact symmetries of the weighted divisor graph which Kondō excludes
from the symplectic action; the nonsymplectic involution of the general
surface acts trivially on its Néron--Severi group.  Realization on a special
period has not been tested.  No claim of mathematical novelty is made pending
a diagram comparison with the existing Coxeter-polytope literature.

## Why a weighted divisor-graph symmetry gives a polytope

Let $r_1,\ldots,r_m$ be the complete set of $(-2)$-curve walls.  Their
classes span the nondegenerate real Néron--Severi space $V$.  If a
permutation $p$ preserves the complete Gram matrix, then

$$
r_i\longmapsto r_{p(i)}
$$

descends from the free wall module to a well-defined isometry $T$ of $V$:
the kernel of the wall map is exactly the radical of the wall Gram, and is
therefore preserved by $p$.  If $p^2=1$, its fixed space $F=V^T$ is spanned
by the orbit sums

$$
n_O=\sum_{i\in O}r_i.
$$

For $x\in F$, the pairings $(x,r_i)$ are constant on each orbit, so the
complete fixed slice of the nef chamber is exactly

$$
C\cap F=\{x\in F:(x,n_O)\geq0\text{ for every wall orbit }O\}.
$$

Thus completeness of the published curve configuration transfers to the
slice.  Exact irredundancy is still required; it is checked computationally
for every orbit wall below.

For a fixed root, the ordinary root reflection restricts integrally.  For a
two-cycle of disjoint roots $\{r,s\}$, the commuting product $s_rs_s$
restricts to the reflection in $r+s$, whose square is $-4$.  Every negative
orbit normal in the two positive folds is of one of these two types.

This argument needs a Gram symmetry, not a surface automorphism.  That is
why the resulting Coxeter polytopes are rigorous even though no
automorphism-fixed divisor-lattice realization has yet been certified for
these two particular slices.

In the two positive cases the verifier goes further than the real-span
argument: it reconstructs every wall in the archived integral curve basis
of $\operatorname{NS}(X)$.  Each graph involution becomes a unimodular
integral order-two isometry, permutes every archived wall exactly, and has
fixed/anti-fixed ranks $(14,1)$ in rank 15 and $(15,1)$ in rank 16.  Thus
these are honest **Néron--Severi lattice isometries**, not merely rational
Gram coincidences.

Their actions on the discriminant group are neither $+1$ nor $-1$.  Hence
they cannot be glued to the scalar $\pm1$ transcendental action available at
a general period.  This exactly separates the lattice construction from a
general-surface Torelli realization.  A special period could have a
non-scalar transcendental Hodge isometry, so that possibility remains open.

## Exhaustive involution inventory

The verifier enumerates the full exact weighted-graph automorphism groups by
row-colour backtracking, preserving zero as well as nonzero intersections.
It finds:

| $\rho$ | weighted graph group order | nonidentity involutions | involutions with fixed dimension $>9$ |
|---:|---:|---:|---:|
| 15 | 12 | 7 | 3 |
| 16 | 4 | 3 | 1 |
| 17 | 2 | 1 | 0 |

The three rank-15 positive involutions form one conjugacy class and hence
give isometric copies of one $19$-facet $H^{13}$ polytope.  Exact
conjugators and their induced orbit relabellings are frozen in the archive.
All remaining involutions have fixed hyperbolic dimension at most eight.

## The $19$-facet polytope in $\mathbb H^{13}$

A representative wall-orbit list is

$$
1,2,3,4,5,6,7,8,9,10,11,12,13,14,(15\ 16),17,(18\ 19),20,21.
$$

The orbit Gram has rank $14$ and inertia $(1,13,5)$.  All 19 orbit normals
are spacelike and irredundant.  Its pair labels are

| relation | count |
|---|---:|
| right angle ($m=2$) | 146 |
| $m=3$ | 14 |
| $m=4$ | 2 |
| parallel | 6 |
| ultraparallel | 3 |

The three ultraparallel normalized squared products are $2,2,9$.  Vinberg's
criterion gives endpoint histogram

$$
(2\text{ finite},0\text{ ideal})^{1005},\quad
(1,1)^{304},\quad
(0,2)^{14}.
$$

Hence all $1323$ elliptic rank-12 edge diagrams have exactly two endpoints.
The full exact Vinberg audit digest is

```text
f08875e243a45d110620af4b4d86cb33c5b417a2b286e524ecdc1d9e1015bcff
```

## The $17$-facet polytope in $\mathbb H^{14}$

The wall-orbit list is

$$
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,(16\ 18),(17\ 19).
$$

The orbit Gram has rank $15$ and inertia $(1,14,2)$.  All 17 orbit normals
are spacelike and irredundant.  Its pair labels are

| relation | count |
|---|---:|
| right angle ($m=2$) | 119 |
| $m=3$ | 14 |
| $m=4$ | 2 |
| parallel | 1 |

There are no ultraparallel pairs.  Vinberg's endpoint histogram is

$$
(2,0)^{549},\quad(1,1)^{148},\quad(0,2)^7.
$$

Thus all $704$ elliptic rank-13 edge diagrams have two endpoints.  The full
exact audit digest is

```text
8a9ae49a6c2fd562b336fc0203726fd02436332eaf25e2a74a094f9571543619
```

## Surface-automorphism boundary

Kondō determines both the weighted graph groups and the subgroup induced by
automorphisms of the general surfaces in these families.

- In rank 15 the weighted graph group is $S_3\times\mathbb Z/2$.  The central
  reflection is induced by a symplectic involution, but its fixed chamber has
  dimension only six.  The three $H^{13}$ folds are the transpositions in
  the $S_3$ factor, which Kondō excludes as symplectic actions.  The other
  factor of the general surface automorphism group acts trivially on the
  Néron--Severi lattice.

- In rank 16 the graph group is $(\mathbb Z/2)^2$.  Kondō realizes one
  reflection and excludes the second because it fixes 15 rational curves,
  contradicting the eight fixed points of a symplectic involution.  The
  excluded reflection is exactly the $H^{14}$ fold above.

- In rank 17 the unique nontrivial graph symmetry is realized
  symplectically, but its fixed chamber is only eight-dimensional.

Consequently this search yields two positive **divisor-wall-graph folds**,
but no certified new automorphism-fixed K3 example above dimension nine.
Whether a special period can realize either graph isometry is a separate
Torelli/discriminant-form problem left open here.

## Reproducibility

- full certificate:
  [`k3_high_rank_involution_folds_rho15_17.json`](../data/enumerations/k3_high_rank_involution_folds_rho15_17.json)
- exhaustive verifier:
  [`search_k3_high_rank_involution_folds.py`](../scripts/search_k3_high_rank_involution_folds.py)
- regression tests:
  [`test_k3_high_rank_involution_folds.py`](../tests/test_k3_high_rank_involution_folds.py)
- source wall data:
  [`k3_high_rank_nef_chambers.json`](../data/known_examples/k3_high_rank_nef_chambers.json)

Full replay:

```bash
PYTHONPATH=src:. python3 scripts/search_k3_high_rank_involution_folds.py \
  --check-archive
```

The frozen certificate file has SHA-256

```text
df6aa10a21faaca264fc120c2459ecc36d5e727b87c6e84411d6b0b47ef345c5
```

Sources: X. Roulleau, *An atlas of K3 surfaces with finite automorphism
group*, Sections 15.1--15.3; S. Kondō, *Algebraic K3 surfaces with finite
automorphism groups*, Sections 4.3--4.5 and Figures 3--5,
doi:10.1017/S0027763000001653.
