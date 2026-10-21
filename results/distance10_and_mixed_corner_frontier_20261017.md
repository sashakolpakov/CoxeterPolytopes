# Distance 10, mixed corner capacities, and continuous cap parameters

The unrestricted finite-volume dimension bound is now **116**, with a
complete internally checked computer-assisted proof. The compact-facet
ambient bound is **87**. The compact-arithmetic ceiling remains **21**;
the noncompact rational arithmetic case is separate. No new finite-volume
example or construction dimension has been established.

## Unrestricted dimension 126 to 116

The [new proof](fresh_selective_distance10_dimension116_20261017.md)
applies to every finite-sided, finite-volume real hyperbolic Coxeter
polytope. It has no arithmetic, compactness, simple-cusp, affine-factor-rank,
or facet-count hypothesis. It pays corner pairs at distance at most 10,
or in connected components of order at most 17.

The complete high-order classification has 384 marked contexts and 527
candidate ports, of which 161 pass the inverse test. Among 1,478 ordered
disjoint joins, 363 have a nonnegative supporting-white obstruction.
Every possible cross bond must be simple. All 35 single-join inertia
types have two negative directions; all 120 multiple-edge subsets have
an actual-edge obstruction. Multiple possible cross edges are explicitly
covered, rather than assuming the uniqueness used at distance 11.

An essential small-order supplement covers radius-three cores of order
five. Proper positive-semidefinite subdiagrams prove that their labels
are at most 5 and that only three tree supports are possible. All 243
recipes are checked exactly. Four survive the preliminary filters, and
their 24 ordered marked joins all have signature (13,2,0). The finite
label range follows from a proof, not a search cutoff.

Of 22,881 candidate overlaps, 4,862 are consistent. Inverse tests exclude
4,296, and actual-edge positive-vector certificates exclude 254 more.
The remaining 312 contexts have 1,248 positive-semidefinite adjacent
corners. Every additional cross contact is excluded by a short marked
path or an actual-edge obstruction.

The [supporting-extension lemma](fresh_distance10_extension_20261017.md)
handles the new positive-residual cases. At most one supporting white
can contact a nonsingular core. Its residual square is 0, 1/2, or 2/3.
Projection onto all indirectly attached supporting whites bounds the
entire tail by zero, one, or two roots, respectively. Every preferred
marked component consequently has order at most 17. This argument covers
arbitrarily large ambient diagrams. The inherited cusp-factor incidence
theorem transfers the payment to nonsimple ideal corners.

Finite cost is at most 10n; ideal cost per link vertex is at most
sum min(10,r_i/2). Cap-20 concentration reduces the new dimensions to
1,015 profiles. Every exact curvature margin is positive, with minimum
145234/4084101. The completed 126 theorem supplies the higher-dimensional
tail, giving 116. The profile (20,20,20,20,20,15) has margin
-673055/16336404 at ambient dimension 116, so the same scalar weights
do not give 115. The compact-facet calculation gives ambient dimension 87.

Primary generation and fresh saved replay pass. The
[separate review](fresh_selective_distance10_dimension116_20261017_review.md)
reconstructs the rational calculations independently and uses independent
Q(sqrt(5)) tree arithmetic for the small-order supplement. It checks all
905 saved nonnegative witnesses, all preferred corners, and every scalar
profile, alongside a written geometric and completeness audit. Its final
generation and fresh saved replay also pass. These are internal agent
reviews, not external refereeing or proof-assistant formalization.

The [163 proof-status note](dimension163_proof_status_20261017.md) remains
applicable: that earlier result was already a complete internally checked
proof. The gap to examples was not a missing mathematical step. No
publication-priority or current-best literature claim is made here.

## Compact arithmetic: mixed capacities and one new solve

The [mixed raw23/raw25 calculation](fresh_mixed_raw23_raw25_corners_20261017.md)
covers 112 actual classes through 1,010 corners, 12,682 literal
orientations and 379 vertex forests. It tests all 535,065 literal-star
pairs, retaining 18,338 and excluding the rest by exact necessary
incidence, acuteness and fixed-field arithmetic conditions. Of the 1,251
exported capacity inequalities, 245 mix the two source families.

With the old face masses fixed, 67 inequalities violate its vertex
allowances. Raising 65 allowances repairs that witness from 12.38203 to
**11.863463132596562**, below the target 12. All previously satisfied
inequalities remain satisfied by monotonicity. This removes that fixed
witness; it does not establish feasibility or a dimension theorem.

One [enlarged finite-model solve](fresh_mixed_actual_corner_full_model_20261017.md)
preserves all 82,806 face obligations, 6,558 vertex budgets, 430,666
inherited variables, 96,176 zero-budget signed columns and 52 signed
vertex-cost columns. It represents all 500 raw25 and 1,251 mixed exports
by 1,355 distinct appended columns. Exact primal and dual certificates
bound the optimum by

    661061434478917472 / 52852828351105989
    <= opt <= 6253797331083 / 500000000000,

or **12.50758862113165 <= opt <= 12.507594662166**. Thus the enlarged
finite model still obstructs the target. The new face masses differ;
the earlier complete local-packing and joint-family certificates cannot
be reused without reconstruction. That reconstruction is the next
arithmetic task. Thirteen newly positive pure raw23/raw25 rows lie outside
the selected cohort, but no claim is made that these omissions alone
explain the obstruction.

The geometric capacity gates use the fixed ground field Q(sqrt(5)).
Complete ambient face, forest, field and charge coverage remains
necessary for a dimension theorem. The compact-arithmetic bound stays 21.

## Construction: a continuous obstruction and another incidence

The [continuous single-X gate](fresh_zero_arc_cap_continuous_20261017.md)
excludes the specified zero-boundary-arc bridge over every low R endpoint
recipe, with no cutoff on the continuous A or X parameters or the finite
X contact labels. Both norm branches of all 96 R recipes are covered.
The 94 branches not excluded immediately give 1,246 X cases and 1,288
interval boxes, resolved by 1,266 monotonicity leaves and one angle-gap
leaf. The certified A-X product is everywhere below -2, contradicting
its required finite corner.

Crucially, the gate uses intrinsic nonintersection with walls 30 and 63,
whose threshold is -sqrt(2*qX). It includes finite native contacts that
become intrinsically nonintersecting after projection. It does not assume
the stronger native dotted-contact condition. The ten symbolic high R
families and other incidences remain outside the theorem. The exact
interval engine is shared with earlier packages; the new calculation has
a narrow source and mathematical review, not an independent arithmetic
implementation.

A [two-connector scout](fresh_zero_arc_cap_two_connectors_20261017.md)
changes the old8 chain to A-X-Y-R. The six-cap disk has 11 vertices,
25 edges and 15 triangles. Four sampled X seeds and 94 low R roots give
13,888 Y norm equations and 27,776 branches. No coherent prefix survives.
This is a bounded numerical search, not an all-parameter exclusion.

The certified control remains 63 walls, 199 finite vertices, ten
target-face open ends and infinite volume. The changed-seed deletion
costs eleven old finite vertices and six compact two-faces; no recovery
of those costs or new chamber is claimed.

## Archive and next work

The [incremental auditor](../scripts/audit_distance10_and_mixed_corner_frontier_20261017.py)
and [saved audit](distance10_and_mixed_corner_frontier_20261017_audit.json)
rehash every frozen predecessor file before extending the archive. They
check the new source bindings, independent scalar composition, matching
primary/review data, mixed fixed-mass loads, finite-model bounds, scoped
construction counts, report links and Python syntax. Frozen historical
reports retain their original wording and hashes.

A further distance reduction must reconsider smaller cores and all
supporting extensions. The arithmetic lane should separate the complete
implicit families at the new masses before another optimization. The
construction lane retains high R parameters, more connectors, changed
incidences and explicit retained-wall deformations as separate domains.
The research objective remains open.
