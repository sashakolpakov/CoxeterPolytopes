# Prepared fixed28-wall and bounded edge-walk diagnostic

19 October 2026. **Prepared source, not executed by this lane.** I
wrote and fully read the new source and parsed its AST without
importing it. No scientific body, edge walk, numerical job, Lean
process or build ran here. The coordinator owns execution and its
separate ledgers; this preparation note makes no execution claim.

Source: `scripts/check_h19_twenty_eight_wall_seed_20261019.py`,
324 lines,16,418 bytes, SHA256
`f9ccdb1dbb5fe6a810bf456be56b3f15d1b83be69fd0d080931713eae3738600`.
Proposed output: `results/h19_twenty_eight_wall_seed_20261019.json`.
The source is stable for the coordinator's reviewed bounded job.

## Fixed algebra and inherited inputs

The source pins the native prefix,26/27-wall certificates,27-wall
checker, executed omega3 search certificate and source, pure field
helpers, and the new completion report. Saved source maps are
rehashed. No imported producer or checker body is called.

It reconstructs all28 roots, checks integral coordinates, norms and
full-lattice reflection ratios, and checks the new h covector exactly.
All378 pairs are classified and compared with the expected totals:
36 angle3,26 angle5,190 right,126 strictly ultraparallel. Its27 new
h products must also equal the saved search ledger.

The old strict center is verified against all28 walls. Each acute
projection is checked directly as a timelike point strict against
every unselected wall, giving28 actual facet witnesses. Four finite
points p,p_s,p_q,p_h are checked against all28 inequalities, with
exactly nine active normals and positive-definite rank-nine Gram.
They are distinct and off both cap hyperplanes. The two complete
41-vertex H8 facets are inherited; those82 vertices plus the four
new checked points give the lower bound86, not a vertex census.

The nine original p directions have exact first positive wall hits
recomputed against all28 walls. Each must be a feasible timelike
point with exactly nine positive-definite incident normals. This
checks nine finite initial edges without claiming global boundedness.

## Bounded finite-vertex walk

The default is --walk-limit64; --walk-limit0 disables the walk.
Every point is normalized to the affine slice x=1. At a finite vertex,
the active covectors restricted to the other nine coordinates form
an invertible9-by-9 matrix. The source uses a NEW row-pivoted exact
Gauss-Jordan solver, verifies its residual, and never applies the
unpivoted inverse helper to that general active matrix.

For each active wall i it solves dx=0,N*d=-e_i, checks all nine
derivative equations, and computes the least positive slack/slope
hit over EVERY28 inequality. Every proposed endpoint is checked
against all walls. Finite endpoints must have exactly nine active
facets with positive-definite Gram; their active set is the exact
vertex key and their affine point is also re-solved independently.

It stops at the first feasible spacelike or ideal endpoint, a no-hit
affine ray,64 discovered finite vertices, or an empty queue. A
spacelike/ideal witness retains its complete28-product ledger. A
no-hit ray records its nonzero spacelike direction and all nonpositive
slopes; the separately checked bounded-slice premises would make
that an inconsistency to investigate, not a compactness certificate.
Finite traversal records only compact incidences and coordinates.

No stop is automatically promoted to finite volume. In particular
the status always retains global_finite_volume_proved=false, including
queue exhaustion. A capped walk with remaining work is partial.
The source retains at most512KiB compact walk records and caps the
final indented certificate below1MiB before exclusive creation.

## Independent read of the coordinator's closure criterion

I independently read
`results/h19_affine_slice_edge_walk_completeness_20261019.md`, SHA256
`c9cfc2806534bc8ef9175aef90e05b616ce163b5a7aa0219fe9679ba5bd37b90`.
**Written PASS, no flaw found.** This is separate from any walk outcome.

The literal old0/r inequalities force x>=0 and y=0 when x=0.
The eight negative E8 simples force Cz>=0. The retained old9 has
strictly positive spatial coefficients, so its inequality at x=y=0
forces Cz=0 and z=0. Thus the closed cone has no nonzero x=0 vector.
The normalized-unbounded-sequence argument then proves that its
x=1 slice is a bounded Euclidean polytope. The script checks the
literal coefficient premises of this argument.

At a finite vertex, acuteness and a strict center rule out any
positive linear dependence of incident normals. Splitting a general
dependence into positive parts in its positive-definite orthogonal
space proves independence. Thus finite vertices are simple, their
active affine matrices are invertible, and all nine -e_i directions
are their genuine Euclidean edges.

If a later exact walk explores all nine edges at every discovered
finite vertex and truly empties the queue with only future timelike
endpoints, its visited graph is a closed component of the FULL
Euclidean vertex graph. The explicitly external connectedness theorem
then gives all vertices. The explicitly external convex-hull theorem
and convexity of the future timelike cone put the whole compact
slice uniformly inside that cone. Hyperboloid renormalization would
therefore give a compact polytope. No partial cap or timeout meets
this criterion. Its hypotheses must be checked before any promotion.

## Required execution controls and scope

There is no inner lock. The coordinator must use the parent-lock
bounded wrapper, required40GiB/6GiB start gates and the ongoing
resource/owner-death guards. The child has a45-second alarm; the
wrapper must retain59-second/1GiB caps. The intended first command,
only after coordinator review and scheduling, is

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
  scripts/check_h19_twenty_eight_wall_seed_20261019.py --write --walk-limit64
```

Replay must use the identical walk limit and --verify in its own
separately gated process. Generation refuses to replace a certificate.
A failure or timeout is not a completed diagnostic. No root census,
all-vertex completeness, finite-volume example, compact H9,
nonreflectivity or dimension record is asserted by source preparation.
