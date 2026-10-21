# Distance nine, complete local packing, and raw39 corners

The unrestricted finite-volume dimension ceiling is now **108**, with a
complete internally checked computer-assisted proof. Its compact-facet
ambient consequence is **79**. Compact arithmetic remains at **21**; the
noncompact rational arithmetic case is separate. No new finite-volume
example has been found.

## Unrestricted finite volume

The [distance-nine proof](fresh_selective_distance9_dimension108_20261017.md)
applies to every finite-sided, finite-volume real hyperbolic Coxeter
polytope, without an arithmetic, compactness, simple-cusp, affine-factor-rank
or facet-count hypothesis. A corner is paid when the marked sides have
distance at most nine or lie in a connected component of order at most 17.

The high-order classification exhausts 1,262 ports, 4,442 ordered disjoint
joins and 1,332 multiple-cross-edge subsets. All disjoint cases fail a
necessary positive-definiteness or Lorentzian-index test. Among 312,316
candidate identifications, 20,718 overlaps are consistent. Inverse and
actual-edge tests reject 19,769, leaving 949 marked contexts in 54 oriented
core types. All additional cross contacts are excluded: 60,831 would give
a short marked path and 888 fail an actual-edge positivity condition.

The [small-core lemma](fresh_distance9_low_core_20261017.md) exhausts all
weighted graphs of orders four and five after proving finite label bounds.
It finds 32 and 14 canonical types. There are 104 low-first joining-label
cases, all excluded; every multiple-edge choice is also covered. Six
surviving low-first overlaps, and their reversals, have isolated preferred
marked components of order at most 12. This supplement is essential to
the complete domain; a high-order catalogue alone does not suffice.

In the new 69-context subset, some outside attachments have residual
square one. Globally bounding every attachment would be invalid. The
new argument instead bounds the component
containing the marked pair. Projection inequalities leave at most one
relevant extension and at most one compatible outside contact. A singular
residual Gram matrix isolates the component; in the remaining cases its
inverse either forbids a tail or forces a single terminal root. The latter
root is orthogonal to every further supporting root, so it cannot connect
indirectly to an unbounded outside component. The
[independent local review](fresh_distance9_supporting_extensions_review_20261017.md)
reconstructs all 69 cases and their 24 relevant extension records.

For ordinary positive-residual extensions, the inverse-diagonal inequality
for a connected positive-definite norm-two Coxeter diagram of order m is
at least m/(m+1). It bounds entire indirectly attached tails, not just
neighbors of the core. Every preferred marked component has order at
most 17, and the new 69-case subset has order at most 16.

The finite payment is at most 9n; ideal payment per link vertex is at most
sum min(9,r_i/2). Cap-18 concentration reduces the new dimensions to 756
exact profiles. All margins are positive, with minimum 45873/2476099.
The completed 116 theorem supplies the higher-dimensional tail. The
ambient-108 boundary profile has margin -4659/260642, so these scalar
weights do not give 107. The compact-facet ambient consequence is 79.

The [full independent review](fresh_selective_distance9_dimension108_20261017_review.md)
reconstructs the high-order classifications and extension calculations
with Fraction arithmetic, checks 4,371 saved nonnegative witnesses, and
recomputes the scalar bounds. Primary generation and fresh saved replay,
and independent review generation and fresh saved replay, all pass.
The small-core proof has a separate
mathematical/source review; its full graph enumeration shares the stated
arithmetic dependencies. These are internal checks, not external refereeing
or Lean formalization. The [163 proof-status note](dimension163_proof_status_20261017.md)
remains applicable to the earlier completed proof. No publication priority
or current-best literature claim is made.

## Compact arithmetic: new masses and new actual corners

The previous enlarged finite model has an exact optimum between
12.50758862113165 and 12.507594662166. The following calculations keep
its new face masses fixed and increase vertex allowances as necessary.
No global optimizer was run for this checkpoint. Local packing LPs supplied
certificates that are replayed with exact arithmetic.

The [complete local packing reconstruction](fresh_mixed_actual_corner_full_model_20261017_local_cones.md)
covers all 2,753 positive simplex features on 760 forests through 5,173
literal block orbits. Its 4,363 exact mixture terms cover the entire
represented nonnegative local-packing family. There are 675 exact local
optima and 85 covering brackets, with no resource-open forest. The 132
new violated charges raise 205 vertex allowances and repair the witness
to 12.405107765824441. All 432,021 exported model columns are checked.

The [joint-port reconstruction](fresh_mixed_actual_corner_full_model_20261017_joint_ports.md)
combines these simplex features with 473 positive actual quadrilateral
rows. Its 3,060,201 selected inequalities contain 183 violations, requiring
33 allowance increases. The repaired ratio is 12.38638766981714.

The [raw39 actual-corner calculation](fresh_raw39_actual_corner_compatibility_20261017.md)
adds all 31 positive classes of the selected canonical square core with
unambiguous complete provenance. It reconstructs 124 corners and 1,246
literal orientations on 88 forests. Exact fixed-field and incidence gates
exclude 21,047 of 21,115 pairs. The 129 resulting capacities contain four
violations and raise three allowances. The final exact fixed-witness ratio is

    247898037929594052000000000000000 /
    20020725665111814911449997327399
    = 12.382070564084598... > 12.

Monotone allowance increases preserve all earlier satisfied constraints.
This remains an obstruction in the stated partial dimension-13 weighting
model over Q(sqrt(5)), not a new arithmetic dimension theorem. Mixed
raw39/raw23/raw25 compatibility and broader face, forest and field coverage
remain open. Root mathematical/source review found no gap in the new
raw39 calculation; it did not independently repeat its census.

## Construction: original single-X incidence closed

The [continuous high-R certificate](fresh_zero_arc_cap_high_r_20261017.md)
covers all ten previously open symbolic high-label families and 150 X
recipes. It includes both norm branches, every endpoint label, and
unbounded continuous A/X parameters. The saved exact interval trees have
18,240 nodes, 9,180 leaves and maximum depth 25; 6,211 leaves are certified
gaps between consecutive Coxeter angles. A finite helper search supplies
gap witnesses, not a cutoff on the candidate labels.

A native Coxeter-angle obstruction strengthens the lower A parameter.
Two algebraically regular charts control the limiting parallel frames
without division by the vanishing frame determinant. Necessary intrinsic
nonintersection with old walls 30, 44 and 63 provides the correct projected
threshold -sqrt(2*qX). The all-parameter contradiction then excludes a
required finite A-X corner.

Together with the frozen low-R certificate, this closes the original
specified single-zero-arc-X incidence. It does not exclude the distinct
two-connector disk, different side incidences or retained-wall deformations.
The interval arithmetic and finite-corner classifier are shared with the
predecessor; the new charts and source received separate narrow mathematical
reviews, not an independent interval implementation.

The control remains 63 walls, 199 finite vertices, ten target-face open
ends and certified infinite volume. The changed seed loses eleven old
finite vertices and six compact two-faces. No changed chamber was promoted.

## Archive and next work

The [incremental auditor](../scripts/audit_distance9_and_raw39_frontier_20261017.py)
and [saved audit](distance9_and_raw39_frontier_20261017_audit.json) preserve
and rehash the [previous archive](distance10_and_mixed_corner_frontier_20261017_audit.json).
They check new source bindings, primary/review correspondence, exact scalar
and fixed-mass compositions, literal packing mixtures, scoped construction
trees, report links and Python syntax. The archive check is not an additional
independent replay of the complete geometry.

Distance eight requires a fresh analysis of order-three cores, including
unbounded black-to-black labels. The arithmetic lane should test mixed
actual-corner families before repeating the same finite optimization.
The construction lane should change the now-excluded single-X incidence.
The research objective remains open.
