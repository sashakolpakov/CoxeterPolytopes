# Edge-driven continuation of the H9 H3 face

This bounded numerical continuation finds new finite-contact chains involving
earlier new walls, but no compact H3 closure. It is a discovery diagnostic,
not an exact cap certificate or a continuous-family exclusion. The exact
55-wall chamber and the frozen two-contact package are unchanged.

The [source](../scripts/fresh_nonarithmetic_h3_edge_completion_20261016.py)
and [saved output](fresh_nonarithmetic_h3_edge_completion_20261016.json)
start from the actual `S6={0,1,2,3,5,6}` face recorded in the
[two-contact data](fresh_nonarithmetic_h3_two_contact_20261016.json).
Its 20 projective extreme rays comprise 13 finite vertices and seven
exterior rays. The source reconstructs 28 finite supporting edges by
requiring a positive-definite two-side residual Gram and exactly two
extreme rays on their common face. Fifteen edges already have two finite
ends; nine have one exterior end and four have two. Intersecting their
projective segments with the null quadric gives **17 open ideal edge ends**.
These are true edge-end targets, rather than radial projections of the
exterior vertices from an unrelated center.

For each open edge, a proposed third wall must give a positive-definite
three-normal residual Gram. This ensures a finite endpoint. Its two
prescribed products leave a norm conic. Open all-other-dotted intervals
are partitioned at the current edge-end pairing roots before sampling.
Every new sample must retain the center and pass all current native
Coxeter contact gates.

The new continuation mechanism solves an additional finite contact to a
third wall whenever at least one of the three anchors is already a new
wall. These isolated equations are solved separately from the open
all-other-dotted intervals. For example, the run finds chains with
anchors `[30,44,new55]`, labels `[2,3,2]`, followed by
`[8,44,new56]`, labels `[2,2,2]`. Such new-new endpoint equations were absent
from the original three-old-anchor bank.

The bounds are six added walls, beam width six, at most 40 clips per state,
and prescribed labels `2,3,4,5,6,8,10`. The beam retains diversity in which
initial edge ends have been cut. It checks 31 intermediate states, 10,846
finite-endpoint conic branches, 29,060 continuous samples, and 530,642
third-contact equations. Numerical candidate deduplication leaves 34,867
candidate occurrences across those states. Near-affine projected triples
are conservatively rejected as possible ideal endpoints; these floating
tests are not independently certified affine obstructions.

The best retained six-cap chain has **23 finite vertices, nine exterior
rays and 19 open edge ends**. It cuts six of the original 17 ideal ends;
replacement open edges remain. Its maximum squared Klein radius decreases
from approximately `1.00116362` to `1.00040511`, but that is still exterior.
Other retained branches achieve smaller maximum radius while leaving more
open ends. These numbers explain why smaller exterior radius or additional
finite vertices alone do not establish closure.

No exact wall is promoted from this run. The next construction domain is
the six nonzero A2 profiles of an actual facet of this H3 face. They change
the source projection and are outside every right-S6 family tested here.

Generation and a fresh saved replay pass. Syntax and source manifests are
checked separately. The source imports the frozen numerical conic and
clipping helpers; its finite-label bound, parameter samples, beam cutoff,
near-zero thresholds and floating clipping preclude any universal failure
claim.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_nonarithmetic_h3_edge_completion_20261016.py --verify
```
