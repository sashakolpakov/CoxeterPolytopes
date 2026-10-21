# Exact rational ordinary-polytopal completion of H9/14 record 0

## Result

**Record 0 has an explicit ordinary convex completion.**  Its missing-edge
graph consists of

\[
\{0,13\},\qquad\{1,13\},
\]

and eleven isolated vertices.  One completed nerve is the boundary of a
rational simplicial convex 9-polytope with fourteen vertices.  Dually, it is
the nerve of a rational simple convex 9-polytope with fourteen facets.

This is not a hyperbolic Coxeter example: no Coxeter labeling or Lorentzian
Gram matrix is supplied, and ordinary polytopality implies neither
hyperbolic compactness nor finite volume.

## Construction

Start with `C(8,5)` in moment-curve order `0,1,...,7`.  Stellar-subdivide
edge `{0,1}` by vertex `8`, then the resulting edge `{6,8}` by vertex `9`.
Apply one simplicial wedge at each of vertices `0,1,6,8`, adding clone
vertices `10,11,12,13`, and use the frozen relabeling.

For `0 <= t <= 7`, use the moment points

\[
p(t)=(t,t^2,t^3,t^4,t^5),
\]

and the two exact rational radial-beyond points

\[
\begin{aligned}
p_8={}&(5/16,-9/16,-179/32,-36,-7235/32),\\
p_9={}&(4849/1536,82289/4608,326777/3072,
        45875/72,11707433/3072).
\end{aligned}
\]

Both use parameter `1/16`.  Exhaustive exact determinants recover 20 cyclic
facets, 28 facets after the first stellar subdivision, and 34 after the
second.

## Completed nerve

Besides the two edges, its minimal nonfaces are six four-faces

```text
{0,3,5,10}, {0,5,8,12}, {1,4,6,11},
{2,3,9,10}, {4,5,7,13}, {6,8,11,12},
```

and seven five-faces

```text
{0,2,5,7,9},    {0,3,6,10,11},  {1,2,4,7,9},
{1,3,6,10,11},  {1,4,7,8,12},   {2,4,5,7,9},
{4,5,7,8,12}.
```

The complete blocker is pure of rank five and contains 168 facet
complements.  The face vector is

\[
(14,89,341,874,1555,1911,1554,756,168),
\]

so this is the `a=6` Dehn--Sommerville schedule.  The exact global audit
verifies inclusion in the frozen support core, degree two for all 756
ridges, connected dual graph, the Lanner-order bound and all strict local
`g`/Macaulay gates, spherical link Euler characteristics, and exact mod-two
sphere homology for every face link.

The source-pinned hostile replay imports neither the record-0 construction
nor its certificate driver.  It independently reconstructs the rational
hulls, expanded clutter, complete blocker, ridge graph, packed face table,
and full face/link audit.

## Algebro-geometric interpretation and scope

The normal fan of the dual rational simple polytope defines a projective
Q-factorial toric 9-fold whose toric boundary-stratum dual complex is this
nerve.  No Delzant or smoothness assertion is made.

The result gives one ordinary-polytopal completion of record 0.  It does not
classify all completions of that graph shadow and does not furnish a Coxeter
realization.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record0.py`
- `scripts/certify_h9_14_record0.py`
- `scripts/independent_h9_14_record0_audit.py`
- `data/enumerations/h9_14_record0_polytopality.json`
- `tests/test_h9_14_record0.py`

The certificate digest is

```text
9b2ab118ebe84e31ba72a2c8b57cf24bf0be12a9314a29356ac457be26da51a4
```

Reproduce it with

```sh
PYTHONPATH=src:. python3 scripts/certify_h9_14_record0.py --verify-only
PYTHONPATH=src:. python3 scripts/independent_h9_14_record0_audit.py
PYTHONPATH=src:. pytest -q tests/test_h9_14_record0.py
```

