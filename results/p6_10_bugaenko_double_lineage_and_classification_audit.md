# `P_{6,10}`: Bugaenko-double lineage and classification audit

Frozen: 2026-08-09.

## Status

The following facts are accepted independently of any six-dimensional
classification claim.

- There is a compact arithmetic Coxeter polytope `P_{6,10}` in `H6` with ten
  facets.
- Its exact unit-normal Gram matrix has rank seven and inertia `(6,1,3)`.
  Its `f`-vector is `(31,93,125,95,42,10,1)`, and its Coxeter Euler
  characteristic is `-67/288000`.
- It is the chamber of an index-two reflection subgroup of the nine-facet
  Bugaenko `H6` group: geometrically it is the double of the Bugaenko chamber
  across one wall.
- This is a known polytope, not a 2026 example. Burcroff's Figure 5 already
  draws it and attributes it to Bugaenko.

The following claims of
[arXiv:2608.02894v1](https://arxiv.org/abs/2608.02894) are **not accepted in
this ledger**:

- uniqueness of `P_{6,10}` among compact Coxeter six-polytopes with ten
  facets;
- the resulting completion of the compact `d+4`-facet classification; and
- any dimension/facet consequence whose proof uses the assertion that every
  facet of a compact Coxeter polytope is itself Coxeter.

The issue is a proof gap in version 1, not evidence that another polytope
exists. The exact construction and compactness of the displayed survivor
remain valid.

## Exact chamber datum

Number the facets `1,...,10` and put

`c_3=1/2`, `c_4=sqrt(2)/2`, and `c_5=(1+sqrt(5))/4`.

The ordinary Coxeter edges `(i,j,m)` are

```text
(1,2,5) (2,3,3) (3,4,3) (4,5,3) (5,6,3)
(9,10,5) (2,7,4) (2,8,4) (5,9,3) (6,10,5).
```

There are three ultraparallel pairs,

```text
{6,7}, {7,8}, {8,9},
```

with exact unit-Gram weights

```text
w67 = w89 = 2*sqrt(2) + sqrt(10),
w78 = 17 + 8*sqrt(5) = w67^2 - 1.
```

Thus the exact Gram matrix `G` is reconstructed without a picture by

```text
Gii = 1;
Gij = -cos(pi/m) on an ordinary m-edge;
G67 = G89 = -w67 and G78 = -w78;
all other off-diagonal entries are zero.
```

The immutable source pin is John Mcleod's repository commit
[`7a6dcc1a7918985c19fd159c8278042e8f51495b`](https://github.com/johnamcleod/10-facets-6-d-polytopes/tree/7a6dcc1a7918985c19fd159c8278042e8f51495b).
The exact diagram and Gram matrix are in
[`paper/polytope.tex`](https://github.com/johnamcleod/10-facets-6-d-polytopes/blob/7a6dcc1a7918985c19fd159c8278042e8f51495b/paper/polytope.tex),
and its exact rank/inertia/compactness replay is
[`verify_polytope.py`](https://github.com/johnamcleod/10-facets-6-d-polytopes/blob/7a6dcc1a7918985c19fd159c8278042e8f51495b/verify_polytope.py).
That source records rank `7`, inertia `(6,1,3)`, no parabolic subdiagram,
`f=(31,93,125,95,42,10,1)`, no ideal vertices,
`chi=-67/288000`, and covolume `67*pi^3/540000`.

## Exact index-two lineage

The pinned local CoxIter corpus contains the following chain.

1. `graphs/cocompact/6-Bugaenko.coxiter` is the compact nine-wall
   Bugaenko diagram. It is label-preservingly isomorphic to
   `graphs/cocompact/6-per09-tumarkin-03.coxiter` after swapping nodes `8`
   and `9`.
2. In the index-two archive, member
   `testing_index2/6-per09-tumarkin-03_0.coxiter` is the same parent diagram,
   with the doubling wall labelled `t0`. Its only finite non-orthogonal
   incidence has label `4`, so the standard even-wall index-two operation is
   admissible.
3. Member `testing_index2/6-per09-tumarkin-03_1.coxiter` deletes `t0` and
   inserts its reflected neighbours. Its Coxeter graph, including every
   ordinary label and dotted incidence, is exactly that of `P_{6,10}`. One
   explicit isomorphism from the numbering above is

```text
1 -> s8       2 -> s4       3 -> s6       4 -> s5       5 -> s3
6 -> s1       7 -> t0_s2    8 -> s2       9 -> t0_s1   10 -> s7.
```

In the Bugaenko numbering, `t0` is wall `8`. Consequently

```text
P_{6,10} = P_Bugaenko union r_t0(P_Bugaenko),
[W_Bugaenko : W(P_{6,10})] = 2.
```

The parent has `f=(23,69,95,75,35,9,1)` and
`chi=-67/576000`; the Euler characteristic and hyperbolic volume therefore
double exactly, as the index-two lineage requires.

### Frozen local content hashes

The CoxIter corpus is pinned at commit
`b800be48d8240fbabd2d1c9542e61054584ad7fe`; its archive provenance is
recorded independently in
[`positive_example_lineage_header_controls.json`](../data/enumerations/positive_example_lineage_header_controls.json).

| object | SHA-256 |
|---|---|
| `graphs/testing_index2.tar.gz` | `3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036` |
| archive member `testing_index2/6-per09-tumarkin-03_1.coxiter` | `2d32aa1bceb38410bf2029b9acc54b5b85671cd0c3d4b72ed8cd6841fb35923c` |
| archive member `testing_index2/6-per09-tumarkin-03_0.coxiter` | `f41d2cbcc728f4e063f098cd4c1ced0221204d155d0f1987edcfff2699148750` |
| `graphs/cocompact/6-Bugaenko.coxiter` | `0ed96024d9da1861c4dc689ba65df1449c8b4e906e444af922ffef6fb50e3b47` |
| `graphs/cocompact/6-per09-tumarkin-03.coxiter` | `d7a0c92eaa30333c0482eda887ba4944983cbb71ba2765f8689e6a4501d96337` |
| `testing/tests.txt` | `9151721b589aba2a531b322706b2e42aede5f3793dbd8980b8bc6e6bca830529` |
| CoxIter `index2.cpp` | `84317a17bd52c71bda9b9ab0441383a2d8e4745dc7370bebf6c2e579196b995b` |

The local test manifest marks the Bugaenko parent and both `_0` and `_1`
members cocompact. The exact graph isomorphisms above were replayed from the
hashed text records, with edge labels included in the matching.

## Why the version-1 uniqueness proof is not imported

The decisive reduction in arXiv:2608.02894v1 assumes that every facet of a
compact Coxeter polytope is itself a compact Coxeter polytope. It uses this
in Lemmas 4.2, 4.4, and 4.5, especially to filter 387 combinatorial types down
to 54 and ultimately search only 11.

That assertion is false. Burcroff states the correct limitation immediately
before Proposition 4.6: a face is acute-angled but need not be Coxeter. The
proposition gives a sufficient Coxeter-face condition only when the elliptic
diagram defining the face has no component of type `A_n` or `D_5`; see
[Burcroff, Proposition 4.6](https://arxiv.org/abs/2201.03437). A facet is
defined by one node, hence by `A_1`, so this proposition gives no such result
for facets.

There is a small exact counterexample. Take the compact Coxeter tetrahedron
with Lanner chain `[5,3,5]` and its facet corresponding to either internal
node. If `theta` is the induced angle between the sides coming from its two
neighbouring nodes, orthogonal projection of the unit normals gives

```text
cos(theta)
  = cos(pi/5)*cos(pi/3) / (sin(pi/5)*sin(pi/3)),
cos(theta)^2 = (5 + 2*sqrt(5))/15.
```

Numerically `pi/5 < theta < pi/4`, so `theta` is not `pi/m` for any integer
`m >= 2`. This facet is therefore not Coxeter.

The exact search and CoxIter verification do validate the one surviving Gram
matrix. They do not prove uniqueness after an invalid filter has removed
unsearched types. A valid completion must either search all 387 types (with
exhaustion certificates) or replace each facet reduction by a theorem whose
hypotheses are actually verified for every eliminated type. Until then the
ledger records `P_{6,10}` as a known compact example and records uniqueness as
open, not as a theorem.

## Mixed-root algebraic-carrier opportunity

The unit Gram field `Q(sqrt(2),sqrt(5))` disguises a more useful golden
crystallographic form. Let `phi=(1+sqrt(5))/2`, set

```text
d = (1,1,1,1,1,1,2,2,1,1),
D = diag(sqrt(d_i)),
B = 2*D*G*D.
```

Then every entry of `B` lies in `Z[phi]`. Its diagonal is `2` except for
the two long roots `7,8`, whose norm is `4`. On ordinary edges its entries
are

```text
m=3: -1,       m=4: -2,       m=5: -phi,
```

and on the three dotted edges they are

```text
B67 = B89 = -4*(1+2*phi),
B78 = -4*(9+16*phi).
```

All reflection coefficients `2*Bij/Bii=Bij/d_i` belong to `Z[phi]`.
Thus `P_{6,10}` has an exact mixed norm-`2`/norm-`4` root presentation over
the golden ring. This is a concrete entry point for the polytope-to-surface
program: after quotienting the three-dimensional Gram radical, one can test
the rank-seven golden module, its rank-fourteen trace lattice, stable even
overclosures, and primitive K3 embeddings. The two long roots strongly
suggest an eigen/folded realization (orbit sums of `(-2)` classes), not a
direct full-nef surface realization. The latter is in any case excluded in
dimension six by the
[direct full-nef compactness obstruction](direct_full_nef_surface_compactness_obstruction.md).

This is an opportunity, not yet a surface carrier.

## Hidden-parent-wall warning

The doubling wall `t0` is absent from the ten-wall boundary of `P_{6,10}` but
is still a genuine root of the parent Bugaenko lattice. It bisects the double
into its two nine-wall chambers. Therefore a surface constructed from the
ambient Bugaenko divisor lattice will normally recover `t0` as an additional
root/negative-curve wall and will produce the parent chamber, not the double.

Consequently the following data are **not sufficient** for a strong carrier:

- an embedding of the visible ten-wall root span;
- a K3 period with the required Picard or eigenlattice;
- correct signature, Coxeter angles, and compactness of the visible Gram; or
- a surface realization of the nine-wall Bugaenko parent.

A strong `P_{6,10}` carrier must prove that `t0` and every other parent-lattice
root is absent from the relevant effective/eigen-nef wall stock, or else that
an actual finite action folds/exchanges the two halves so that `t0` ceases to
be a wall of the invariant cone. A complete projected-root enumeration and
hidden-cut audit is therefore mandatory. This is the exact next gate for the
mixed-root carrier lane.
