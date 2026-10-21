# Exact rational ordinary-polytopal completion of H9/14 record 4

## Result

**Record 4 has an explicit ordinary convex completion.**  Its missing-edge
graph is the three-edge path

\[
\{0,12\},\quad\{0,13\},\quad\{1,13\},
\]

together with nine isolated vertices.  One completed nerve is the boundary
of a rational simplicial convex 9-polytope with fourteen vertices.  Dually,
it is the nerve of a rational simple convex 9-polytope with fourteen facets.

This is **not** a hyperbolic Coxeter example.  There are no Coxeter labels or
Lorentzian Gram matrix, and ordinary polytopality does not imply hyperbolic
compactness or finite volume.

## Construction

Start with `C(8,5)` in moment-curve order `0,1,...,7`.  Stellar-subdivide
edge `{0,1}` by vertex `8`, then edge `{0,7}` by vertex `9`.  Apply one
simplicial wedge at each of vertices `3,4,5,6`, adding clone vertices
`10,11,12,13`, and use the frozen target relabeling.

The exact rational realization uses the moment points

\[
p(t)=(t,t^2,t^3,t^4,t^5),\qquad 0\leq t\leq7,
\]

and

\[
\begin{aligned}
p_8={}&(5/16,-9/16,-179/32,-36,-7235/32),\\
p_9={}&(2705/768,57745/2304,271585/1536,
        358073/288,13407505/1536).
\end{aligned}
\]

Both are exact radial-beyond placements with parameter `1/16`.  Exhaustive
rational supporting-hyperplane determinants recover 20 facets for `C(8,5)`,
28 after the first subdivision, and 36 after the second.

## Completed nerve

Its minimal nonfaces are the three edges above, the two four-faces

```text
{1,2,4,9}, {1,2,5,10},
```

and the ten five-faces

```text
{1,3,5,8,10},   {2,4,6,9,11},   {2,4,7,9,12},
{2,5,7,10,12},  {3,5,7,8,10},   {3,5,8,10,13},
{3,6,7,8,11},   {3,6,8,11,13},  {4,6,7,9,11},
{4,6,9,11,13}.
```

The complete blocker is pure of rank five with 146 members.  The full face
vector is

\[
(14,88,330,823,1424,1708,1364,657,146),
\]

the `a=2` Dehn--Sommerville schedule.  The independent global audit verifies
that all facets lie in the frozen support core, all 657 ridges have degree
two, the dual graph is connected, the Lanner-order and all strict local
`g`/Macaulay gates hold, and every face link has both the correct Euler
characteristic and exact mod-two sphere homology.

The hostile replay imports neither the record-4 construction nor certificate
module.  It independently reconstructs the rational hulls, expanded clutter,
complete blocker, ridge graph, packed face table, and full face/link audit.

## Algebro-geometric interpretation and scope

The normal fan of the dual rational simple polytope defines a projective
Q-factorial toric 9-fold whose toric boundary-stratum dual complex is this
nerve.  No Delzant or smoothness statement is made.

This moves record 4 out of the theorem-undecided class by giving one explicit
ordinary-polytopal completion.  It does not classify every completion of the
same graph shadow, and it does not furnish a Coxeter realization.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record4.py`
- `scripts/certify_h9_14_record4.py`
- `scripts/independent_h9_14_record4_audit.py`
- `data/enumerations/h9_14_record4_polytopality.json`
- `tests/test_h9_14_record4.py`

The certificate digest is

```text
0625af6609333638803df22868791179f3a659c47819bd762a07f1bd069cb1fc
```

Reproduce it with

```sh
PYTHONPATH=src:. python3 scripts/certify_h9_14_record4.py --verify-only
PYTHONPATH=src:. python3 scripts/independent_h9_14_record4_audit.py
PYTHONPATH=src:. pytest -q tests/test_h9_14_record4.py
```

