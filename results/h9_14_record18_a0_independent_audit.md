# Record 18: corrected independent `a=0` compact audit

Status: the corrected `a=0` branch contains no compact hyperbolic Coxeter
polytope.  This result is compact-only and makes no finite-volume noncompact
claim.

The earlier assertion that the whole blocker of the eight five-sets must be
three-uniform is not used.  The sound reduction uses only the two one-centre
blocker sections: each is a nonempty complete bipartite cut, blocker duality
forces a balanced `5|5` cut, and only the 60 visible three-transversals enter
the ridge equations.  The literal blocker of all five missing edges and all
eight missing five-sets is then constructed before accepting a target.

## Exact census comparison

Two structurally different searches independently return 196 ridge-valid
literal families.  The sets of all 196 eight-row families are exactly equal.
An explicit incidence-matrix normal form, quotienting row permutations, the
eight isolate columns, and the order-eight active double-star automorphism
group, gives exactly three classes with multiplicities

```
144, 36, 16.
```

This comparison uses neither graph fingerprints nor equality of counts as a
substitute for equality of the families.

## Independent topology replay

For one independently selected representative of every class, a
standard-library-only replay reconstructs the complete face table.  Every
class has

```
f = (14, 86, 310, 735, 1204, 1372, 1052, 495, 110),
blocker sizes = {5: 110},
ridge degrees = {2: 495}.
```

The dual graph is connected; the reconstructed minimal nonfaces are exactly
the five missing edges and eight missing five-sets; every face link has the
spherical Euler characteristic and the exact mod-two homology of a sphere.
The strict local `g` and Macaulay gates also pass.

## Independent Coxeter-label replay

The replay regenerates all 420 labelled placements of the five classical
rank-five Lanner diagrams, without importing the production label gate.  It
then joins placements and tests every restriction induced by each of the 110
actual rank-nine vertices against the finite Coxeter classification.  The
three terminal sequences are

```
420 -> 96 -> 6 -> 0,
420 -> 96 -> 6 -> 0,
420 -> 96 -> 6 -> 4 -> 0.
```

Thus all three topologically valid classes have zero Coxeter labelings.
Unassigned pairs are tested at label 2; increasing a label increases the
nonnegative Coxeter adjacency matrix and cannot decrease its Perron root, so
a failed elliptic restriction cannot be repaired later.

## Artifacts

- `scripts/independent_h9_14_record18_a0_balanced_cut_audit.cpp`
- `scripts/certify_independent_h9_14_record18_a0_audit.py`
- `scripts/independent_h9_14_record18_a0_downstream_audit.py`
- `data/enumerations/h9_14_record18_a0_independent_audit.json`
- `data/enumerations/h9_14_record18_a0_downstream_independent_audit.json`
- `tests/test_h9_14_record18_a0_independent_audit.py`

The focused test suite passes: `3 passed`.
