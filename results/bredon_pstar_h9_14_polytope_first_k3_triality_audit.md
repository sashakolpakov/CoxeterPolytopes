# Bredon `P_star` in `H9`: polytope-first K3/triality audit

**Status (2026-08-09): exact one-polytope carrier and return, but not a
strict class theorem.**  Starting only from the Coxeter diagram of Bredon's
exceptional 14-facet ADEG polytope, the arithmetic reconstruction has two
primitive integral root-length branches.  Exactly one has a Lorentzian
simply-laced `C3/D4` unfolding of K3 rank.  That unfolding is integrally and
equivariantly identical to the already frozen `A18_4_0_C07` K3 construction,
and its invariant nef cone returns the original polytope exactly.  The other
branch has rank 28 and three negative directions and is impossible for a K3
surface.

This is a genuine **polytope -> surface -> same polytope** computation.  It is
an `L3` exact-return calibration, not an `L4` inverse theorem: the singleton
input does not define a non-singleton class, and algebraic geometry has not
yet pushed a new conclusion back to every member of such a class.

The machine-readable certificate is
[`bredon_pstar_h9_14_k3_triality_audit.json`](../data/enumerations/bredon_pstar_h9_14_k3_triality_audit.json)
(SHA256
`65b313b2d8006556bfc46711f85721a52e993a9d402ea85b934d61499b8a6f6b`).
It pins and rechecks the three upstream K3/factory certificates rather than
silently importing their conclusions.

## 1. Published input and exact Gram reconstruction

The input is Figure 5 and Theorem 4.2 of Naomi Bredon,
[*On ADEG-polyhedra in hyperbolic spaces*](https://arxiv.org/abs/2507.05153)
(2025).  The diagram has two end vertices `L,R` joined by four paths

```text
L --3-- a_j --6-- b_j --3-- c_j --3-- R,   j=0,1,2,3.
```

In the canonical order

```text
(L,a0,b0,c0,a1,b1,c1,a2,b2,c2,a3,b3,c3,R),
```

the normalized Gram has diagonal `1`, entry `-1/2` on every `3`-edge,
entry `-sqrt(3)/2` on every `6`-edge, and zero elsewhere.  Its exact matrix is
stored in the certificate.

Contracting the `3`-edges leaves two components.  A `6`-edge must change the
primitive root square by a factor of three, so there are exactly two global
root-length orientations, interchanged by reversing long and short.

For the K3-admissible orientation,

```text
long (square 6):  L,a0,a1,a2,a3
short (square 2): b0,c0,b1,c1,b2,c2,b3,c3,R.
```

The integral Gram `B` is therefore determined without a matrix search:

- `B_ii=6` on the five long vertices and `B_ii=2` on the nine short ones;
- `B_ij=-3` on long--long `3`-edges and long--short `6`-edges;
- `B_ij=-1` on short--short `3`-edges;
- all other entries vanish.

Exact invariants are

```text
rank(B) = 10,       inertia(B) = (9,1,4),
SNF(B) = 1^6 3^4 0^4.
```

The nondegenerate divisor lattice, in divisor-sign convention, is

```text
U + A2(-1)^4,       signature (1,9),       determinant -81.
```

The fourteen primitive facet classes span this lattice with index one.  A
fixed permutation identifies `B` exactly with the Gram in the independently
constructed `A18_4_0_C07` factory certificate.

The reversed orientation has nine square-six and five square-two roots.  Its
folded Gram still has rank 10 and Lorentzian inertia, but

```text
SNF = 1^4 3^6 0^4,
candidate divisor lattice = U(3) + A2(-1)^4,
determinant = -729.
```

## 2. Literal finite volume, not a cusp-count shortcut

The exact double-description and Vinberg gates give

```text
finite vertices = 128,
ideal vertices  =   6,
hyperideal      =   0,
total vertices  = 134.
```

All extreme rays lie on one time sheet, with a strict interior vector whose
facet-pairing range is `[26,276]`.  The literal edge audit finds 671 elliptic
edge diagrams, 14 connected affine components, zero bad edge endpoints, and
zero bad parabolic neighbourhoods.  Thus `P_star` is **finite-volume and
noncompact**.  This conclusion uses the full extreme-ray/edge/parabolic
tests; the existence or number of affine subdiagrams alone would not prove
finite volume.

The replay recovers the published `f`-vector entries `134` vertices, `671`
edges, and `14` facets.  Bredon's published full vector is

```text
(134,671,1480,1909,1606,917,356,91,14).
```

No new exact volume constant is claimed.  Bredon states the volume as a
rational multiple of `zeta(5)/22295347200` without determining the rational
factor in the form needed here; the present certificate proves literal
finite covolume, not that missing factor.

## 3. The polytope-first `D4`-triality unfolding

Replace every square-six root by three mutually orthogonal square-two roots
whose sum is the long root.  A `6`-edge then becomes a `D4` star: its short
root meets all three clones once.  A long--long `3`-edge becomes three
coordinatewise simple edges, and a short--short `3`-edge remains one simple
edge.  Cycling the three clones defines `C3` triality.

For the admissible orientation this gives a 24-node simply-laced Gram `A`
with

```text
rank(A) = 18,       inertia(A) = (17,1,6),
SNF(A) = 1^14 2^4 0^6,
charpoly(A) = t^6 (t-4)^6 (t-2)^10 (t^2-4t-4).
```

The orbit-sum identity is exact:

```text
orbit_sum * A * orbit_sum^T = B.
```

After reversing signs to the divisor convention, the nondegenerate quotient
is the Nikulin lattice of type `(18,4,0)`, determinant `-16`.  The certificate
gives an explicit permutation to 24 short source walls in the archived K3
model and proves all three stronger facts:

1. their Gram is exactly `-A`;
2. they span the rank-18 Neron--Severi lattice with index one;
3. canonical clone-cycling is exactly the archived `C07` order-three action.

For the reversed root-length orientation, the same forced construction has
32 nodes, rank 28, and inertia `(25,3,4)`.  It is neither Lorentzian nor of K3
Picard rank.  Thus the surface constraints select one of the two arithmetic
orientations rather than merely recognizing a carrier found in advance.

## 4. Algebraic carrier and exact return

The surviving unfolding has a primitive K3-lattice extension with

```text
NS(X) = Nikulin (18,4,0),
T(X)  = U(2) + U(2).
```

The period/Torelli certificate produces a genuine purely non-symplectic
order-three automorphism `g`.  Its coinvariant lattice is
`-(A2 tensor D4)`, has positive minimum four, and is `(-2)`-rootless.  Hence
no hidden coinvariant roots cut the fixed cone.  All restricted long source
walls are exactly redundant, and the decisive return identity is

```text
factory cone = P_r intersect Fix(g) = Nef(X) intersect Fix(g).
```

For the normal quotient `Y=X/<g>`, pullback identifies `N^1(Y)_R` with the
fixed subspace and identifies `Nef(Y)` with the same cone.  This is a statement
about the normal quotient, not the nef cone of a resolution.

Consequently the route really starts with Bredon's polytope, constructs the
integral triality cover and K3 period, and returns the same chamber.  It is not
the weaker observation that some preselected surface happens to have a
related divisor lattice.

## 5. Why this does not yet cover all 24 ADEG polytopes

There is a useful exact criterion for this particular carrier mechanism.
Contract every connected component of the `m=3` subgraph.  A primitive
integral long/short assignment, equivalently the simply-laced `C3/D4`
unfolding above, exists exactly when the resulting `m=6` multigraph

1. has no loop, and
2. is bipartite.

Necessity follows because a `3`-edge preserves root length and a `6`-edge
switches it by a factor of three.  Conversely a bipartition supplies the two
root lengths and the clone construction above supplies the unfolding.

`P_star` passes: the contraction has two vertices joined by four parallel
`6`-edges.  But Bredon's 24-member census already contains two decisive
two-dimensional controls:

- a triangle with one `6`-edge and two `3`-edges produces an `m=6` loop after
  `m=3` contraction;
- a triangle with three `6`-edges produces an odd `m=6` cycle.

Therefore no theorem based on this integral symmetric `D4`-triality
unfolding can carry **every** ADEG polytope.  This refutes uniformity of this
route only; it is not an exhaustive census of all possible surface carriers
for all 24 rows.

## 6. Strict inverse score and lineage

| target | C | Car | Return | AG | Pushback | verdict |
|---|---|---|---|---|---|---|
| singleton `P_star` | fail: singleton | pass | pass, exact nef return | pass: K3 embedding, period and Torelli | fail | `L3`, not `L4` |
| all 24 ADEG polytopes | pass: independently defined non-singleton class | fail for this route | fail class-wide | conditional only on the `P_star` branch | fail | no inverse theorem |

The novelty split is strict:

- the abstract `P_star` is Bredon's published-new 2025 polytope;
- its arithmetic chamber/lattice lineage is already known and is identified
  in the factory ledger with Scharlau--Kirschmer `II_(9,1)(3^4)`;
- the order-three K3 deformation type is classical;
- the project-level addition is the exact marked triality reconstruction and
  the invariant-nef/normal-quotient-nef realization.

In particular this audit does **not** claim a new abstract Coxeter polytope,
a new reflection group, a new volume formula, or an AG-derived improvement
of Bredon's published `24`-row classification and `n <= 11` bound.

## Reproduction

```bash
python3 scripts/certify_bredon_pstar_h9_14_k3_triality_audit.py --check-archive
```

The certifier is
[`certify_bredon_pstar_h9_14_k3_triality_audit.py`](../scripts/certify_bredon_pstar_h9_14_k3_triality_audit.py).
