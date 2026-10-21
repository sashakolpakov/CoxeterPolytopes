# Record 18: an exact rational stellar/multiwedge realization

## Result

**TOOT.** The selected H9/14 record-18 nerve is the boundary of a rational
simplicial convex 9-polytope on fourteen vertices.  Dually, it is the nerve
of an ordinary rational simple 9-polytope with fourteen facets.

This upgrades the previously audited abstract candidate to an ordinary-
polytopality theorem.  It is not a hyperbolic Coxeter example: no Coxeter
labels or Lorentzian Gram matrix are supplied, and neither compactness nor
finite hyperbolic volume follows.

## Inverse multiwedge

Collapse the two clone pairs

\[
\{5,6\}\mapsto5,\qquad \{7,8\}\mapsto7.
\]

The resulting twelve-vertex simplicial 6-sphere has polytope dimension seven
and face vector

\[
(12,61,175,310,341,217,62).
\]

It is obtained from the cyclic 7-polytope \(C(10,7)\), with moment-curve
labels in the order

\[
(1,0,10,5,9,2,3,4,7,11),
\]

by two boundary stellar subdivisions:

1. subdivide the tetrahedral face \(\{2,4,5,11\}\) with new vertex 12;
2. subdivide the triangular face \(\{0,5,7\}\) with new vertex 13.

Gale evenness gives 40 facets initially.  The subdivisions give 52 and then
62 facets.  Expanding carriers 5 and 7 by simplicial wedges raises the
dimension from seven to nine and reproduces exactly the five missing edges
and eight missing five-faces of the selected record-18 nerve.

## Exact rational realization

For the ten cyclic vertices use

\[
(t,t^2,\ldots,t^7),\qquad t=0,\ldots,9,
\]

in the displayed label order.  Both new vertices use the same elementary
rational beyond-placement.  If \(y\) is the barycentre of the face being
subdivided and \(p\) is the barycentre of all current vertices, put

\[
v=y+\frac1{32}(y-p).
\]

This gives

\[
v_{12}=\left(
\frac{387}{64},\frac{2649}{64},\frac{19791}{64},
\frac{782277}{320},\frac{1283031}{64},
\frac{10784145}{64},\frac{92216511}{64}\right)
\]

and

\[
v_{13}=\left(
\frac{8151}{2048},\frac{50197}{2048},\frac{366579}{2048},
\frac{14189393}{10240},\frac{22420827}{2048},
\frac{178390381}{2048},\frac{1422538659}{2048}\right).
\]

The certificate verifies both barycentric identities and examines every
seven-point candidate supporting hyperplane using exact rational
determinants.  It finds precisely the 62 simplicial facets constructed by the
two stellar subdivisions.  Independently, it derives the complete minimal-
nonface clutter of this base, expands the two clone carriers, and compares it
literally with the frozen record-18 clutter.

Boundary stellar subdivision and simplicial wedge preserve rational convex
polytopality.  Thus the exact base realization proves the full
nine-dimensional ordinary-polytopality claim without inferring real
realizability from the previously archived oriented matroid.

## Algebraic-geometric carrier

Polarizing after a rational translation and applying the two rational facet
wedges gives a rational simple 9-polytope.  Its rational normal fan defines a
projective \(\mathbb Q\)-factorial toric 9-fold whose toric boundary-stratum
dual complex is this nerve.  No Delzant or smoothness assertion is made.

## Hyperbolic boundary

The next gates remain genuinely metric:

1. find Coxeter labels compatible with all minimal nonfaces;
2. solve the exact Gram equations with Lorentzian inertia \((9,1,4)\); and
3. only for a complete wall system, apply Vinberg's full elliptic/parabolic
   endpoint-extension criterion.

Ordinary convexity, ridge degree two, and spherical links do not imply
Vinberg compactness or finite covolume.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record18_polytopality.py`
- `scripts/certify_h9_14_record18_polytopality.py`
- `scripts/independent_h9_14_record18_polytopality_audit.py`
- `data/enumerations/h9_14_record18_polytopality.json`
- `tests/test_h9_14_record18_polytopality.py`

The proof-source and certificate digests are respectively

```text
9de308c9a4dc54b84daaa2a8e4c38ca062da4dee28f352a41a2beceabaff165e
6dae56855a8380ee8cffd8b3e2572df90164f98e483e717f9323de8289d4cb90
```

Replay with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record18_polytopality.py
python3 scripts/independent_h9_14_record18_polytopality_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record18_polytopality.py
```
