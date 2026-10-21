# Record 17: an exact rational stellar/multiwedge realization

## Result

The selected H9/14 record-17 nerve is the boundary of a rational simplicial
convex 9-polytope on fourteen vertices.  Dually, it is the nerve of an
ordinary rational simple 9-polytope with fourteen facets.

This is an ordinary-polytopality theorem, not a hyperbolic Coxeter example.
No Coxeter labels or Lorentzian Gram matrix are supplied, and neither
compactness nor finite volume in hyperbolic space follows.

## Inverse multiwedge

Collapse the four clone pairs

\[
\{3,4\}\mapsto14,\quad \{5,6\}\mapsto15,\quad
\{7,8\}\mapsto16,\quad \{9,10\}\mapsto17.
\]

The resulting ten-vertex simplicial 4-sphere has face vector

\[
(10,40,80,80,32).
\]

It is obtained from the cyclic 5-polytope \(C(8,5)\), with moment-curve
labels in the order

\[
(1,0,11,15,17,14,16,2),
\]

by two boundary stellar subdivisions:

1. subdivide the triangle \(\{1,14,15\}\) with new vertex 13;
2. subdivide the triangle \(\{2,14,15\}\) with new vertex 12.

Gale evenness gives 20 facets initially; the two subdivisions give 26 and
then 32 facets.  Expanding the four carrier vertices by simplicial wedges
raises the dimension from five to nine and reproduces exactly the five
missing edges and eight missing 5-faces of the selected record-17 nerve.

## Exact rational check

For the eight cyclic vertices use the moment curve
\((t,t^2,t^3,t^4,t^5)\), \(t=0,\ldots,7\), in the displayed label order.
The two new vertices are

\[
v_{13}=\left(
\frac{249}{94},\frac{1053}{94},\frac{2334}{47},
\frac{21423}{94},\frac{50262}{47}\right)
\]

and

\[
v_{12}=\left(
\frac{39419}{7614},\frac{219847}{7614},\frac{658768}{3807},
\frac{8300773}{7614},\frac{27033112}{3807}\right).
\]

An exact rational determinant replay examines every five-point candidate
supporting hyperplane.  It finds exactly 32 supporting facets, all
simplicial, and exactly the abstract facets produced by the two stellar
subdivisions.  Independently, the complete minimal-nonface set of this base
is expanded through the four clone blocks and compared literally with the
record-17 clutter.

Simplicial wedges preserve rational convex polytopality, so this verifies the
full 9-dimensional realization without relying on an oriented-matroid
realizability inference.

## Algebraic-geometric carrier

Polarizing after a rational translation and applying the four rational facet
wedges gives a rational simple 9-polytope.  Its rational normal fan defines a
projective \(\mathbb Q\)-factorial toric 9-fold whose toric boundary-stratum
dual complex is this nerve.  No Delzant or smoothness assertion is made.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record17_polytopality.py`
- `scripts/certify_h9_14_record17_polytopality.py`
- `scripts/independent_h9_14_record17_polytopality_audit.py`
- `data/enumerations/h9_14_record17_polytopality.json`
- `tests/test_h9_14_record17_polytopality.py`

The certificate digest is

```text
d17332be82e6bddb0ebfe320d678238cc0fa182f2a1b03d81931217826c72e31
```

Replay with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record17_polytopality.py
python3 scripts/independent_h9_14_record17_polytopality_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record17_polytopality.py
```

