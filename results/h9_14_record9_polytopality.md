# Record 9: an exact rational stellar/multiwedge realization

## Result

The selected H9/14 record-9 nerve is the boundary of a rational simplicial
convex 9-polytope on fourteen vertices.  Dually, it is the nerve of an
ordinary rational simple 9-polytope with fourteen facets.

This is an ordinary-polytopality theorem, not a hyperbolic Coxeter example.
No Coxeter labels or Lorentzian Gram matrix are supplied, and neither
compactness nor finite hyperbolic volume follows.

## Construction

Start with the cyclic 5-polytope (C(8,5)), whose moment-curve vertices have
parameters (0,ldots,7).  Perform two boundary stellar subdivisions:

1. subdivide the edge ({2,3}) with new vertex 8;
2. subdivide the triangle ({0,5,7}) with new vertex 9.

The facet counts are (20\to28\to34), and the resulting ten-vertex
simplicial 5-polytope has face vector

\[
(10,41,84,85,34).
\]

Replace base carriers 5, 6, and 7 by clone blocks of sizes 2, 3, and 2,
respectively.  These are four simplicial wedges in total, raising the
dimension from five to nine.  The frozen singleton relabeling makes the four
missing edges exactly

\[
\{0,12\},\quad\{0,13\},\quad\{1,13\},\quad\{12,13\},
\]

which is the graph shadow of record 9.

## Exact checks

At each stellar step the new point is the carrier barycentre displaced away
from the current vertex barycentre by the rational factor (1/32).  An
exact rational determinant enumeration checks every candidate supporting
hyperplane and recovers literally the abstract 28- and 34-facet stellar
subdivisions.

After the four wedges, the complete minimal-nonface clutter consists of the
four displayed edges, the two missing four-faces

\[
\{1,3,5,6\},\qquad \{3,5,6,12\},
\]

and eight missing five-faces.  Its complete blocker has 126 members, all of
order five.  The generated boundary has

\[
f=(14,87,320,778,1308,1526,1192,567,126).
\]

The generic exact audit independently verifies membership in the frozen
ridge-support core, ridge histogram `2:567`, dual connectivity, the complete
face vector and minimal-nonface clutter, every local strict (g)- and
Macaulay gate, spherical Euler characteristic for every link, and mod-two
homology-sphere status for every face link.  A separate hostile standard-
library replay reconstructs the cyclic polytope, both stellar subdivisions,
the exact coordinate support facets, the multiwedge clutter, the blocker,
the face vector, and the ridge graph without importing the proof module.

## Algebraic-geometric carrier

Polarizing the rational realization and applying the four rational facet
wedges gives a rational simple 9-polytope.  Its rational normal fan defines a
projective \(\mathbb Q\)-factorial toric 9-fold whose toric boundary-stratum
dual complex is this nerve.  No Delzant or smoothness assertion is made.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record9.py`
- `scripts/certify_h9_14_record9.py`
- `scripts/independent_h9_14_record9_audit.py`
- `data/enumerations/h9_14_record9_polytopality.json`
- `tests/test_h9_14_record9.py`

The internal certificate digest is

```text
844b80aba9da0431fdefcdc9a2fbd93b5eb72e37f81fb0344b33a7c00848fdd2
```

Replay with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record9.py
python3 scripts/independent_h9_14_record9_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record9.py
```
