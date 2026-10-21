# Record 3: an exact rational two-edge-stellar multiwedge

## Result

The selected H9/14 record-3 nerve is the boundary of a rational simplicial
convex 9-polytope on fourteen vertices.  Dually, it is the nerve of an
ordinary rational simple 9-polytope with fourteen facets.

This is an ordinary-polytopality theorem, not a hyperbolic Coxeter example.
No Coxeter labels or Lorentzian Gram matrix are supplied, and neither
compactness nor finite hyperbolic volume follows.

## Construction

Start with the cyclic 5-polytope (C(8,5)) on moment-curve parameters
(0,ldots,7).  Stellar-subdivide the edge ({0,1}) with vertex 8, then
stellar-subdivide the edge ({4,8}) with vertex 9.  The facet counts are

\[
20\longrightarrow28\longrightarrow36,
\]

and the resulting ten-vertex simplicial 5-polytope has face vector

\[
(10,42,88,90,36).
\]

Apply simplicial wedges at base carriers 2, 3, 6, and 9.  Under the frozen
construction-to-target relabeling the two remaining missing edges are

\[
\{0,12\},\qquad\{1,13\},
\]

which is precisely the graph shadow of record 3.  The expanded clutter also
has four missing four-faces and nine missing five-faces.

## Exact checks

The two new base points are rational radial-beyond placements with factors
(1/16) and (1/64).  An exhaustive exact rational determinant census of
all candidate supporting hyperplanes recovers literally the abstract 28-
and 36-facet stellar subdivisions.

The complete blocker of the expanded clutter contains 162 sets, all of
order five.  Its generated boundary has

\[
f=(14,89,340,866,1528,1862,1504,729,162).
\]

The generic exact audit verifies frozen support-core membership, ridge
histogram `2:729`, dual connectivity, the face vector and complete minimal-
nonface clutter, all strict local (g)- and Macaulay gates, spherical Euler
characteristic for every link, and mod-two homology-sphere status for every
face link.  A standalone hostile standard-library replay independently
reconstructs the coordinates, supporting facets, multiwedge clutter,
blocker, target face vector, ridge degrees, and dual connectivity.

## Algebraic-geometric carrier

Polarizing the rational realization and applying the four rational facet
wedges gives a rational simple 9-polytope.  Its rational normal fan defines a
projective \(\mathbb Q\)-factorial toric 9-fold whose toric boundary-stratum
dual complex is this nerve.  No Delzant or smoothness assertion is made.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record3.py`
- `scripts/certify_h9_14_record3.py`
- `scripts/independent_h9_14_record3_audit.py`
- `data/enumerations/h9_14_record3_polytopality.json`
- `tests/test_h9_14_record3.py`

The internal certificate digest is

```text
500c0874643b21cc5f6c0eadaf278695d9516db1832a895afb84042fa1da4685
```

Replay with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record3.py
python3 scripts/independent_h9_14_record3_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record3.py
```
