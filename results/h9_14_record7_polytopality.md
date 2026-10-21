# Record 7: an exact rational stellar/multiwedge realization

## Result

The selected H9/14 record-7 nerve is the boundary of a rational simplicial
convex 9-polytope on fourteen vertices.  Dually, it is the nerve of an
ordinary rational simple 9-polytope with fourteen facets.

This is an ordinary-polytopality theorem, not a hyperbolic Coxeter example.
No Coxeter labels or Lorentzian Gram matrix are supplied, and neither
compactness nor finite hyperbolic volume follows.

## Construction

Start with the cyclic 5-polytope (C(8,5)), using moment-curve parameters
(0,ldots,7).  Perform two boundary stellar subdivisions:

1. subdivide the edge ({0,1}) with vertex 8;
2. subdivide the triangle ({0,3,5}) with vertex 9.

The facet counts are (20\to28\to34), and the resulting ten-vertex
simplicial 5-polytope has face vector

\[
(10,41,84,85,34).
\]

Replace base carriers 3, 4, 5, and 6 by clone pairs.  These four simplicial
wedges raise the dimension from five to nine.  Under the frozen relabeling,
the missing edges become exactly

\[
\{0,12\},\quad\{0,13\},\quad\{1,13\},\quad\{2,13\},
\]

the graph shadow of record 7.  The other twelve expanded base minimal
nonfaces all have order five; there are no missing four-faces.

## Exact checks

The new point in the first stellar subdivision is displaced radially beyond
its carrier by the rational factor (1/16); the second uses (1/128).
An exhaustive exact rational determinant enumeration checks every candidate
supporting hyperplane and recovers literally the abstract 28- and 34-facet
stellar subdivisions.

The complete blocker of the expanded minimal-nonface clutter has 130
members, all of order five.  Its generated boundary has

\[
f=(14,87,320,780,1320,1554,1224,585,130).
\]

The generic exact audit verifies membership in the frozen support core,
ridge histogram `2:585`, dual connectivity, the complete face vector and
minimal-nonface clutter, every local strict (g)- and Macaulay gate,
spherical Euler characteristic for every link, and mod-two homology-sphere
status for every face link.  A standalone hostile standard-library replay
independently reconstructs the exact coordinates and supporting facets, the
multiwedge clutter, blocker, face vector, ridge degrees, and dual graph.

## Algebraic-geometric carrier

Polarizing the rational realization and applying the four rational facet
wedges gives a rational simple 9-polytope.  Its rational normal fan defines a
projective \(\mathbb Q\)-factorial toric 9-fold whose toric boundary-stratum
dual complex is this nerve.  No Delzant or smoothness assertion is made.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record7.py`
- `scripts/certify_h9_14_record7.py`
- `scripts/independent_h9_14_record7_audit.py`
- `data/enumerations/h9_14_record7_polytopality.json`
- `tests/test_h9_14_record7.py`

The discovery-only finder is
`scripts/explore_h9_14_record7_stellar.py`; the certificate and hostile
replay do not trust it.

The internal certificate digest is

```text
788c8d4a16456810e719e7c0f973149edc6e0a29bab067af28d89c47f6600759
```

Replay with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record7.py
python3 scripts/independent_h9_14_record7_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record7.py
```
