# Record 32: an exact convex multiwedge realization

## Theorem

H9/14 second-layer survivor record 32 is the nerve of an ordinary convex
simple 9-polytope with fourteen facets.  Equivalently, the explicit complex
below is the boundary of a convex simplicial 9-polytope.

This is an ordinary-polytopality theorem, not a hyperbolic Coxeter example:
no Coxeter labels or Lorentzian Gram matrix are supplied.

The graph6 string is

```text
M?????????O?O?g??
```

and the complete minimal-nonface clutter is

\[
\begin{aligned}
E={}&\{\{0,11\},\{1,12\},\{0,13\},\{2,13\}\},\\
Q_4={}&\{\{2,3,4,5\},\{3,4,5,11\}\},\\
Q_5={}&\{\{6,7,8,9,10\}\}.
\end{aligned}
\tag{1}
\]

The complete blocker has exactly 110 elements, all of order five.  Their
complements generate a pure eight-dimensional complex with face vector

\[
(14,87,318,760,1242,1400,1060,495,110).
\tag{2}
\]

Every one of its 495 ridges has degree two, the dual graph is connected, all
strict local \(g\)-vector and Macaulay tests pass, and the complex and every
face link have the mod-two homology of spheres of the required dimensions.
Convexity, however, follows from the independent construction below rather
than from these necessary incidence and homology checks.

## Inverse multiwedge

Put

\[
C=\{3,4,5\}
\]

and collapse this clone block to a carrier vertex \(c\).  Remove for the
moment the separate five-set \(A=\{6,7,8,9,10\}\).  The resulting
seven-vertex base has labels

\[
\{0,1,2,11,12,13,c\}
\]

and minimal nonfaces

\[
\{0,11\},\{0,13\},\{2,13\},\{2,c\},\{11,c\},\{1,12\}.
\tag{3}
\]

The first five form a five-cycle.  Equivalently, the allowed equatorial
cycle is

\[
11-13-c-0-2-11,
\tag{4}
\]

while \(1\) and \(12\) are the two suspension vertices.  Thus the base is

\[
\partial C_5 * \partial\Delta^1,
\]

the boundary of a convex pentagonal bipyramid.  Its face vector is
\((7,15,10)\).

Apply two simplicial wedges at \(c\), replacing it by the full clone block
\(C\).  Every base minimal nonface containing \(c\) acquires all three
clones, so

\[
\{2,c\}\longmapsto\{2,3,4,5\},\qquad
\{11,c\}\longmapsto\{3,4,5,11\}.
\tag{5}
\]

All other base missing edges remain unchanged.  Finally join with
\(\partial\Delta^4\) on \(A\).  This adds exactly the minimal nonface
\(A\) and no others.  Equations (3)--(5) therefore reproduce (1) exactly.

Simplicial wedges preserve convex polytopality, and joins of polytopal
boundaries are polytopal.  The two wedges raise the base dimension from
three to five, while the simplex-boundary join raises it from five to nine.
The vertex count is

\[
7+2+5=14.
\]

This proves the theorem.  Dually, the desired simple polytope is the product
with \(\Delta^4\) of two successive facet wedges over a pentagonal prism.

## Algebraic-geometric carrier

Choose the pentagon in (4) Delzant—for example, as the moment polygon of a
toric blow-up of \(\mathbb P^2\)—and take its product with an interval.  A
facet wedge with primitively normalized integral slack preserves the Delzant
condition, as does a product with a simplex.  Consequently record 32 also
occurs as the dual intersection complex of the reduced toric boundary on a
smooth projective toric 9-fold.  This is an SNC log-Calabi--Yau boundary
\(D\) with \(K_X+D=0\).

This supplies an algebro-geometric carrier for the incidence design.  It
does not manufacture the Coxeter angle data: those require an additional
Lorentzian divisor/intersection form and an exact Gram-matrix certificate.
Indeed the associated smooth toric 9-fold has Picard rank
`14-9=5`, so its ordinary Hodge-index form has signature `(1,4)`, not the
rank-ten signature `(1,9)` required for an `H9` wall Gram.  The toric object
is therefore an automatic **incidence** factory, while a K3-style marked
eigenlattice or another enhanced divisor structure would still have to
supply the metric.

## Independent Gale cross-check

The 110 facet complements also pass the uniform rank-four Gale
oriented-matroid gate.  The frozen certificate contains 1001 chirotope signs;
an independent verifier checks all 45,045 Grassmann--Pluecker relations and
recovers precisely the same 110 positive five-circuits.  This is consistent
auxiliary evidence only.  The wedge-and-join construction, not the oriented-
matroid gate, proves real convex realizability.

## Consequence for the frontier

Record 32 moves from the bounded-timeout class to the ordinary-polytopal
calibration class.  The 142 graph shadows are therefore partitioned as

```text
122 independently impossible + 7 ordinary-polytopal + 13 undecided.
```

The remaining undecided indices are

```text
0, 1, 3, 4, 5, 6, 7, 9, 12, 14, 17, 18, 31.
```

The compact label problem for this selected nerve has since been decided
negatively by the complete Lanner/actual-vertex certificate in
`results/h9_14_selected_nerve_coxeter_no_go.md`.  That theorem concerns this
selection only; another completion of the record-32 graph shadow, and any
noncompact ideal-face structure, remain separate questions.

## Certificate

The frozen certificate is

```text
data/enumerations/h9_14_record32_polytopality.json
```

with digest

```text
94a65380b2d56d4787a3de28d5cca5dede0540cc2ec0ff23b2a3d65b4cfbd318
```

Reproduce it with

```sh
PYTHONPATH=.:src python3 scripts/certify_h9_14_record32.py
PYTHONPATH=.:src python3 scripts/independent_h9_14_record32_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record32.py
```

The search backend is provenance only.  Both verification commands replay
from the seven displayed minimal nonfaces and the archived signs without
invoking the facet-selection or Gale solvers.
