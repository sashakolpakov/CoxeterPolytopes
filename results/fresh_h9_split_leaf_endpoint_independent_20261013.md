# Independent verification of the fixed-bank exclusion

The exact endpoint and cover argument checks out: **no finite-volume
Coxeter completion can retain the sixteen seed facets and add facet
mirrors only from this frozen bank of 3,034 distinct mirrors**. This
does not exclude new mirrors, changes to the seed facets, or other H9
constructions.

The [independent audit script](../scripts/fresh_h9_split_leaf_endpoint_independent_20261013.py)
and [result](fresh_h9_split_leaf_endpoint_independent_20261013.json)
verify the argument without a SAT solver, numerical normals, or either
construction producer. The arithmetic comes from the separate exact
Fraction normal audit. Source file hashes agree with the final endpoint,
degree-graph, and cover certificates.

## Endpoint signs checked with reversed extensions

Each of the two intersections of eight mutually orthogonal seed walls
is a hyperbolic line. Its two ideal endpoints are expressed as
`P+sqrt(d)R` and `P-sqrt(d)R`, where `d=45+74*phi`. The audit directly
checks `P^2<0`, `R^2>0`, `(P,R)=0`, `P^2+d*R^2=0`, future orientation,
exactly eight zero seed pairings, and strict negativity on the other
eight seed walls.

It then recomputes **all 12,160 cap-endpoint pairing signs**. The
constructor places the cap radical in the inner quadratic extension;
this audit instead places the endpoint radical inside and the cap
radical outside. Exact comparisons of squares in the inner field handle
dependent and zero radicals. All four nonnegative lists agree exactly:
each has fourteen entries, including two entries with zero pairing.

These lists give necessary conditions for finite volume. If every
selected cap is strictly negative at an endpoint, move slightly farther
out on the same line, replacing `sqrt(d)` by `sqrt(d)+epsilon` with
`epsilon>0`. The resulting vector is spacelike. The eight incident
seed pairings remain zero, and all other finitely many inequalities
remain strict for sufficiently small epsilon. A further sufficiently
small perturbation toward the common strict center makes the eight
remaining inequalities strict while preserving positive norm. This
gives an open infinite-volume escape. Therefore every finite-volume
selection must include a cap with **nonnegative**, not necessarily
positive, pairing at each endpoint.

## Exhaustive covers and a common exact escape

Using the separately certified
[degree compatibility supergraph](fresh_h9_split_leaf_3034_degree_graph_20261013.md),
the audit enumerates all `14^4=38,416` choices of one cap from each
necessary list. Exactly 36 distinct choices are compatible even in this
conservative supergraph: the pairs `(i,j)` with `0<=i<=5` and
`6<=j<=11`.

Every compatible selection satisfying all four conditions must contain
one of those pairs. It therefore lies in that pair's common-neighbor
set, including the pair itself. All 36 such sets have 28 entries, and
their union contains 332 entries. The audit reconstructs every set and
checks exact agreement with the
[cover certificate](fresh_h9_split_leaf_endpoint_cover_scout_20261013.json).

In the sixteen seed-normal coordinates, the single vector

```
w = [8,8,8,1,-3,-3,-3,-3,-3,-3,-3,5,5,5,5,5]
```

has

```
(w,w) = (7242-2800*phi)/11 > 0,
(w,S) = -116-232*phi < 0.
```

Here `S` is the verified common strict timelike center. The audit checks
that `w` pairs **strictly negatively with all sixteen seed walls and
all 332 cap entries**. Thus every remaining selection has a strict
future spacelike escape and infinite volume. Selections that fail an
endpoint condition already have the preceding endpoint escape. This
exhausts the literal bank without an UNSAT claim about the earlier SAT
run.

The compatibility graph retains every possible finite Coxeter label
allowed by the cap fields, every potentially parallel or dotted pair,
and unresolved intervals. Its exact rejection proof is a separate
prerequisite, previously replayed in full. The cover calculation uses
that graph rather than the earlier numerical graph. The separate normal
and cutter audits also verify cap norms, the Lorentz signature, and all
literal covector identities.

```sh
python scripts/fresh_h9_split_leaf_endpoint_independent_20261013.py --verify
```

The endpoint producer's replay and the cover/witness producer's replay
also passed independently during this review. The conclusion is a
construction obstruction for this fixed dictionary, not a new polytope
and not a general H9 nonexistence theorem.
