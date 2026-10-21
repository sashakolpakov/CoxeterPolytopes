# Research priority: bounds near the known examples

The user's requested scale is within approximately two to five dimensions
of the largest known examples: an unrestricted finite-volume ceiling around
23–26, and a compact ceiling around 10–13. The proposed arithmetic compact
bound 21 and unrestricted finite-volume bound 163 remain research candidates.
Neither a further one-dimensional improvement nor a closed finite search bank
fulfills this larger objective.

The [exact scalar checker](../scripts/fresh_dimension_gap_method_barriers_20261016.py)
and [data](fresh_dimension_gap_method_barriers_20261016.json) identify two
specific barriers in the current unrestricted estimates. These are restrictions
of particular arguments, not nonexistence theorems or claims about the best
possible dimension bounds.

The current [163 candidate](fresh_nikulin_two_budget_frontier_20261015_dimension163.md)
uses affine rank N=n-1, total target (2N-2)/9, and a balanced endpoint-pair
lower bound equal to one half the retained energy minus two. Its loss-four
estimate gives (N-37)/9, reaching the strict narrowness threshold 14 at N=163.
Even the hypothetical removal of every unit of loss would give only

    (N-19)/9.

That reaches 14 at N=145, corresponding to hyperbolic dimension 146. Thus
optimizing only the nonnegative loss while keeping this target, normalization,
and narrowness constant cannot bring the proof near dimension 26. The balanced
height ratio one remains an allowed case. Better targets or additional
geometric constraints would be a substantive change, outside this barrier.

There is a stronger existing warning: the
[rank-162 local configurations](fresh_nikulin_rank162_local_obstruction_20261016.md)
satisfy the current Gram, endpoint-budget and narrow-pair tests, but have
explicit infinite-volume escapes. A useful new theorem must use information
that those local configurations omit, such as simultaneous completion of
other faces, additional actual framing walls, or a justified change of
the selected narrow face. Arbitrary linear combinations of roots cannot be
treated as actual walls.

The older Prokhorov scalar estimate has a separate barrier. With the inherited
mean-facet and triangular-bipyramid proportion bounds, its total-weight lower
estimate per three-face is

    (n^2-70n+216)/(4(n-4)(n-3))   for even n,
    (n^2-68n+147)/(4(n-3)(n-2))   for odd n.

These are negative through dimension 66. Consequently, decreasing nonnegative
edge weights alone cannot give a contradiction in the twenties using those
unchanged scalar bounds. A better treatment of exceptional three-faces or
additional incidence relations is required. This limitation is about the
specified inequalities, not every refinement of Prokhorov's method.

For compact arithmetic polytopes the independently recorded
[small-component barriers](fresh_sqrt5_discount_structural_floor_20261016.md)
show that the existing rank-at-least-six discount family cannot even exclude
dimension 17. Reaching the low teens requires modifying weights on small
components or imposing relations between different faces. Larger discounts
also enlarge the required quadrilateral domains; the current baseline-one
attachment queues alone would not certify them.

The active construction searches in dimensions 9 and 20 remain relevant:
they can produce new examples, and their exact failed completion mechanisms
can suggest global constraints missing from the dimension arguments. A local
obstruction becomes a dimension theorem only after an exhaustive reduction
showing that every hypothetical polytope encounters it.

```sh
python3 scripts/fresh_dimension_gap_method_barriers_20261016.py --verify
```
