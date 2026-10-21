# Exact rational ordinary-polytopal completion of H9/14 record 12

> **Subsequent compact theorem.**  The corrected whole-shadow census has 371
> residual orbits, and every one now has zero compact Coxeter label states.
> Thus neither this selected ordinary polytope nor any other completion of the
> record-12 shadow is a compact Coxeter nerve; see
> [`h9_14_record12_all_coxeter_label_gate.md`](h9_14_record12_all_coxeter_label_gate.md).

## Result

**Record 12 has an explicit ordinary convex completion.**  Its missing-edge
graph is

\[
K_{2,2}[\{0,1\},\{12,13\}]\sqcup 10K_1,
\]

and one completed nerve is the boundary of a rational simplicial convex
9-polytope with fourteen vertices.  Dually, it is the nerve of a rational
simple convex 9-polytope with fourteen facets.

This is **not** a hyperbolic Coxeter example.  No Coxeter labels or
Lorentzian Gram matrix are supplied, and ordinary convex polytopality alone
says nothing about hyperbolic compactness or finite volume.

## Construction

Begin with the cyclic simplicial 5-polytope `C(8,5)` whose vertices, in
moment-curve order, are `0,1,...,7`.  Perform two stellar subdivisions:

1. subdivide the triangular face `{0,2,4}` by a new vertex `8`;
2. subdivide the resulting triangular face `{0,1,8}` by a new vertex `9`.

Then apply one simplicial wedge at each of vertices `1,2,3,4`, adding clone
vertices `10,11,12,13`, respectively.  Finally apply the fixed relabeling
recorded in the certificate module.  The dimensions progress as

\[
5\xrightarrow{\text{two stellar subdivisions}}5
 \xrightarrow{\text{four wedges}}9,
\]

while the vertex counts progress from 8 to 10 to 14.

The two stellar subdivisions are certified by exact rational coordinates.
Vertices `0,...,7` are the moment points

\[
p(t)=(t,t^2,t^3,t^4,t^5),
\]

and the new vertices are

\[
\begin{aligned}
p_8&=(29/16,85/16,59/4,463/16,-229/4),\\
p_9&=(221/256,1279/768,169/64,-4643/768,-7671/64).
\end{aligned}
\]

They are respectively the exact radial-beyond placements with parameters
`1/8` and `1/32`.  Exhaustive rational determinants recover exactly the 26
facets after the first subdivision and the 34 facets after the second.

## Completed nerve

The minimal nonfaces are the four missing edges

```text
{0,12}, {1,12}, {0,13}, {1,13},
```

the two missing four-faces

```text
{2,3,8,12}, {5,7,10,13},
```

and the eight missing five-faces

```text
{1,2,4,7,9},   {1,2,5,7,10},  {0,3,5,8,10},
{3,5,7,8,10},  {0,3,6,8,11},  {0,4,6,9,11},
{2,4,6,9,11},  {4,6,9,11,13}.
```

Their complete blocker is pure of rank five and has 126 elements.  These are
the complements of the 126 facets.  The full face vector is

\[
(14,87,320,778,1308,1526,1192,567,126).
\]

Thus this is the `a=2` Dehn--Sommerville schedule.  The independent global
audit verifies:

- all 126 facet complements lie in the frozen record-12 support core;
- every ridge has degree two (567 ridges total);
- the dual graph is connected;
- the minimal-nonface order bound is five;
- every strict local `g` and Macaulay condition used by the production model;
- spherical Euler characteristic for every face link; and
- the exact mod-two homology of a sphere for every face link.

The hostile replay imports neither the construction module nor the
certificate driver.  It independently reconstructs the cyclic facets,
stellar operations, rational hulls, wedge-expanded clutter, complete
blocker, ridge graph, full face table, and generic face/link audit.

## Algebro-geometric interpretation and scope

The normal fan of the dual rational simple polytope defines a projective
Q-factorial toric 9-fold.  Its toric boundary-stratum dual complex is this
nerve.  No Delzant or smoothness claim is made.

The result moves record 12 out of the theorem-undecided class by producing
one explicit ordinary-polytopal completion.  It does not classify all
completions of the record-12 graph shadow.  A separate Coxeter-label/Gram
analysis would be needed even for this selected nerve.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record12.py`
- `scripts/certify_h9_14_record12.py`
- `scripts/independent_h9_14_record12_audit.py`
- `data/enumerations/h9_14_record12_polytopality.json`
- `tests/test_h9_14_record12.py`

The certificate digest is

```text
cd98e07ccb15d04654c951e3271239be88d842b21c7bfa5f3c7f4c48d5875ab1
```

Reproduce it with

```sh
PYTHONPATH=src:. python3 scripts/certify_h9_14_record12.py --verify-only
PYTHONPATH=src:. python3 scripts/independent_h9_14_record12_audit.py
PYTHONPATH=src:. pytest -q tests/test_h9_14_record12.py
```
