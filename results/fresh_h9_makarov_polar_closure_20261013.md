# Exact exclusion of a polar bank, allowing old-wall deletion

The [exact verifier](../scripts/fresh_h9_makarov_polar_closure_20261013_verify.py)
and [certificate](fresh_h9_makarov_polar_closure_20261013_verify.json) pass.
No Coxeter-compatible subset of the **1,564 literal polar recipes**, together
with any subset of the preceding 27 walls, can give finite volume.
This permits deleting any or all old walls. It is a statement about this
finite bank, not an exhaustive theorem about every polar of the chamber,
arbitrary new walls, or all H9 polytopes.

## Why this was a different construction attempt

Each proposal truncates a spacelike projective corner by its polar
hyperplane, giving nine right-angle pairings automatically. The previous
finite-edge construction prescribed finite endpoint diagrams instead.
Seven polar proposals pass the numerical angle gates against all 27 old
walls, but none cuts the original strict escape. The broader experiment
therefore allowed old-wall deletion and used all 1,564 proposed polars.

The full bank passes a numerical enclosure scout when pairwise Coxeter
compatibility is ignored. A compatible covering search nevertheless fails.
The exact result below establishes that this incompatibility is decisive;
it does not infer failure from a time limit or from the numerical graph.

## Exact recipe definitions and witnesses

The [frozen recipe file](fresh_h9_makarov_polar_closure_20261013_recipes.json)
specifies nine original wall indices for each proposed polar. Their exact
normals come from the frozen 27-wall chamber. The checker proves rank nine
by interval Gaussian elimination with certified nonzero pivots. It defines
the unique null line, proves that line spacelike, chooses the orientation
with negative pairing against the common center, and normalizes to square
two. Thus the exact normal is defined algebraically; the numerical corner
coordinate is not used to decide an exact sign.

All intervals use integer arithmetic with outward rounding at scale 10^70.
Twenty-one saved integer vectors are proved spacelike, in the correct
center half-space, and strict for all 27 original wall inequalities. They
are therefore retained after **every** deletion of old walls.

For each witness and every polar recipe, its full pairing interval is
saved. A polar is included among a witness's possible cutters whenever the
interval does not certify a strictly negative pairing. This is a
conservative superset: every polar that could remove the witness is kept.
All polars outside that set leave the witness strict.

## Exact incompatibility and exhaustive covering test

Both normalized polar normals have square two. Positive mutual pairings
are forbidden. Other pairs are rejected only when their enclosing interval
is disjoint from all permissible Coxeter values -2 cos(pi/m) and the
parallel/dotted range. Rigorous Machin/Taylor intervals cover m=2,...,120;
the entire interval from -2 up to the label-120 bound is retained. Thus
every larger label and every unresolved comparison remain allowed.

The checker exhausts the necessary covering problem: choose mutually
compatible polar recipes that hit every witness's possible-cutter set.
At each node it chooses an uncovered witness and branches on each remaining
possible cutter. Selecting one removes only candidates with a certified
forbidden pair, plus itself. After a failed branch, excluding that candidate
from later branches is valid because every covering solution must choose
some cutter of the selected witness. There is no heuristic pruning or
unproved pair rejection. A node limit raises an error; it cannot produce
a successful exclusion certificate.

The conservative search exhausts in **35 nodes**, using 31,308 certified
pair conflicts. Its full tree and every conflict interval are saved.
Hence any compatible subset misses at least one witness's possible-cutter
set. That integer spacelike vector is strict for all selected polar walls
and for every retained old wall. Together with the common strict timelike
center, it gives an open patch at infinity and infinite volume.

Only the archived finite bank is excluded. The proposal stage used a
numerical projective-corner census with unresolved exterior incidences;
the exact result makes no completeness claim for that census. Walls outside
these nine-wall polar definitions, or another seed, remain possible.

```sh
python3 scripts/fresh_h9_makarov_polar_closure_20261013_verify.py --verify
```

The numerical SAT and DFS pilots are retained separately as proposal
history. They are not dependencies of the exact geometric sign checks.

The [proof review and separate covering replay](fresh_h9_makarov_polar_closure_20261013_review.md)
pass. The latter independently rechecks the saved sign decisions, all
31,308 conflict intervals, and the compatibility cover. It uses the saved
intervals rather than independently reconstructing their geometric values.
