# Intrinsic facet incidence and the combined H20 cap bank

The new finite numerical bank cannot close the target **compact H3 face**
`S16+151`. It combines 50,314 roots from the one-parameter families, solved
partners, earlier three-contact bank and jointly solved clusters. Requiring
each selected wall to acquire enough finite neighbors reduces the bank to
1,227 possible roots. A complete finite-graph search excludes every compatible
cover of the 61 saved boundary samples that also satisfies those necessary
degree conditions. A separate binary-program implementation agrees.

This is a finite floating-point diagnostic, not an exact exclusion of the
continuous families or a nonexistence theorem for finite-volume H20 chambers.
No new finite-volume example or construction dimension is established.

## The geometric reason for the screen

Let a finite acute hyperbolic chamber have a strict interior timelike vector
`H`, with outward spacelike facet normals `e_i`. For any positive definite
principal Gram matrix `B`, its off-diagonal entries are nonpositive, so
`B^{-1}` is entrywise nonnegative. Write `h_i=(H,e_i)<0` on that subset.
The orthogonal projection

`H' = H - sum_i (B^{-1} h)_i e_i`

is timelike, lies on the selected hyperplanes and is strictly inside every
other wall: the added root coefficients are nonnegative and every remaining
root pairing is nonpositive. The path of these projections stays in the
same time component. Thus a positive definite selected Gram gives an actual
face, including after any further acute walls containing the common center
are added. It cannot disappear merely by being called redundant.

Applied to `S16+151+r`, this shows that every candidate with positive
projected norm is an actual polygon facet of the target H3 face. If that
face is compact, each such polygon needs at least three finite neighboring
facets. Removing a candidate that cannot acquire three neighbors is therefore
a necessary-condition reduction. The removal can be iterated. During a
clique search, the same test applies after restricting to roots compatible
with the already selected roots. A selected wall with fewer than three
possible finite neighbors makes that branch impossible.

The test is necessary only; degree three does not establish a compact
polygon, much less closure of the full chamber. In particular, a future
search allowing ideal vertices needs the corresponding ideal-incidence
conditions and is not excluded by this compact-face screen.

## Projection, rather than ambient angle, determines incidence

Let `X_j` be an old normal projected into `S16`'s orthogonal complement,
with Lorentz metric `M`, and let `u=X_151`. To project once more into
`151`'s orthogonal complement, use

`qbar_j = X_j.M.X_j - (X_j.M.u)^2/(u.M.u)`

and

`pairbar_rs = r.M.s - (r.M.u)(s.M.u)/(u.M.u)`.

A pair gives a finite edge when its projected two-by-two Gram is positive
definite. The native old squared norm must **not** replace `X_j.M.X_j`.
In the saved numerical geometry the only positive old projected norms are
for walls **95,142,145,150,152,156,157**. In particular, old anchors31 and132
do not meet this H3 in hyperbolic space.

Among the 41,984 continuous-family samples, 21,695 have zero old finite
neighbors, 13,594 have one, and only6,695 have two. Most therefore require
two or three new finite contacts. All40 sampled-bank covers and all16
jointly solved clusters have a selected facet of degree below three.
This identifies a structural defect, independent of how small their
remaining exterior radius becomes.

## Numerical results and independent check

The [original all-sample cover](fresh_h20_continuous_sample_joint_cover_20261016.json)
used 41,984 samples plus219 solved partners, with no coverage-ranking
truncation. Forty compatible covers were checked by full three-dimensional
halfspace intersection. The best maximum squared Klein radius was
1.0002072692459414, still exterior. At iteration41 the complete search on
61 boundary samples exhausted947,013 states and found no cover.
Its [source](../scripts/fresh_h20_continuous_sample_joint_cover_20261016.py)
and fresh saved replay pass.

The [incidence screen](fresh_h20_intrinsic_facet_incidence_20261016.json)
adds the8,009 earlier three-contact roots and102 roots from the16 coordinated
clusters. Among50,314 rows it finds4,382 numerical finite intrinsic edges.
Iterated degree removal leaves1,227 rows:559 original continuous samples,
58 solved partners,523 three-contact roots and87 coordinated roots.
All61 boundary samples still have cutters in that remaining core.
The [source](../scripts/fresh_h20_intrinsic_facet_incidence_20261016.py)
and fresh saved replay pass.

The [combined cover search](fresh_h20_incidence_core_joint_cover_20261016.json)
requires coverage, pairwise Coxeter compatibility and degree at least three
for **every selected new facet**, including any extra walls added solely
to supply finite neighbors. It finishes in52 states with complete
infeasibility; neither the node guard nor the iteration guard is reached.
The [source](../scripts/fresh_h20_incidence_core_joint_cover_20261016.py)
and fresh saved replay pass.

The [separate binary-program review](fresh_h20_incidence_core_milp_review_20261016.json)
imports neither incidence producer nor either clique-search implementation.
It reconstructs the1,227-root core directly from the saved root inventories
and checks61 sample constraints,309,203 incompatible pairs and865 conditional
degree constraints. Its310,129-row program is infeasible. This separately
checks the remaining core, **not** the preceding all-bank removal calculation
or the exact signs of the numerical geometric inputs.
The [review source](../scripts/fresh_h20_incidence_core_milp_review_20261016.py)
and fresh saved replay pass.

## Next construction domain

Two prescribed finite old-side contacts, together with the151 contact,
leave only one continuous parameter per new wall. The incidence census
shows why coordinating those families is restrictive: many walls still
need two or three new finite contacts.

The next search retains151 and **one** genuine finite old-side contact,
leaving a two-dimensional norm quadric. A second old pairing can be
continuously dotted instead of prescribed to a finite Coxeter angle.
At fixed dotted pairing the existing conic interval calculation applies;
varying that pairing supplies another parameter for joint angle recovery.
Intervals should be split at every saved boundary-sample pairing root,
rather than represented only by a coarse parameter grid. These additional
families are outside the exclusion above.

The [old-edge frontier](fresh_h20_open_edge_frontier_20261016.json), computed
by its [separate source](../scripts/fresh_h20_open_edge_frontier_20261016.py),
gives a more directed starting point. There are seven finite vertices and
13 finite edges among the seven old facets. Eight edges already have both
finite endpoints. The other five are:

| Open edge | Existing third facet | Endpoint candidates in the remaining degree core |
|---|---|---:|
| 95–142 | 156 | 170 |
| 95–157 | 156 | 2 |
| 142–152 | 150 | 44 |
| 145–152 | 150 | 2 |
| 145–157 | 156 | 19 |

The endpoint test requires a positive definite projected three-root Gram,
not just three pairwise finite contacts. Every finite edge of a compact
H3 polyhedron needs two finite endpoints. Closing these five old edges,
then generating new walls at unfinished edges involving the new walls,
is therefore a concrete continuation strategy. New conic frames may use
one or zero old facets, outside the previous fixed old-anchor inventories.
The frontier generation and fresh saved replay pass. It remains numerical;
no completeness assertion for this continuous edge-driven search is made.

```sh
python3 scripts/fresh_h20_continuous_sample_joint_cover_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_intrinsic_facet_incidence_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_incidence_core_joint_cover_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_incidence_core_milp_review_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_open_edge_frontier_20261016.py --verify
```
