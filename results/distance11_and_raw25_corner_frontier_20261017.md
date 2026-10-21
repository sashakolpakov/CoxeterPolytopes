# Distance 11, raw-25 corners, and changed cap adjacency

The unrestricted finite-volume dimension bound is now **126**, with a
complete internally checked computer-assisted proof. Primary generation
and fresh replay pass, as do a separate mathematical review and independent
Fraction implementation. The compact-facet ambient bound is **95**.
Compact arithmetic remains **21**. No new finite-volume example has been
established.

## Unrestricted dimension 136 to 126

The [new argument](fresh_selective_distance11_dimension126_20261017.md)
pays corner pairs at distance at most 11, together with pairs in connected
components of order at most 17. It gives dimension at most **126** for
finite-sided, finite-volume real hyperbolic Coxeter polytopes, without
arithmeticity, compactness, simple-cusp, or facet-count assumptions. Its
compact-facet consequence gives ambient dimension at most **95**.

The exact classification tests 200 distant ports in all 384 marked
minimal-core contexts, leaving 55 ports. Among 425 disjoint joins, 119
have an explicit nonnegative vector of nonpositive square in the supporting
white diagram. The other 306 each have a unique possible cross edge.
Positive definiteness forces it to be simple; their full Grams then have
two negative directions, excluding them all.

The overlap census has 2,137 candidate identifications and 1,010 consistent
ones. Inverse tests exclude 856; exact positive-vector certificates in
the support of an actual edge exclude another 40. The remaining 114
contexts comprise 60 shared-vertex, 27 shared-edge, nine shared-three-node-path,
and 18 oriented mixed L22/L25 cases. All unspecified cross contacts are
checked: 8,280 create a short path and 45 contradict positive definiteness.

The shared-vertex radical relation and shared-edge inverse-Gram argument
extend to the larger marked domains. The new
[mixed overlap lemma](fresh_distance11_extension_review_20261017_mixed.md)
isolates its core through an actual radical relation. The
[three-node overlap lemma](fresh_distance11_extension_review_20261017.md)
admits two possible extra roots, each lying in the core span. Their mutual
product is -2, so at most one can be supporting. All other supporting
whites are orthogonal to that extension. The marked component consequently
has order at most 16. Actual cusp-factor incidence gives the same corner
payment at nonsimple ideal vertices.

The costs are 11n at finite vertices and sum min(11,r_i/2) per cusp-link
vertex. Cap-22 concentration reduces dimensions 127 through 136 to 1,095
profiles. Every exact margin is positive; the smallest is 40603/559682.
The completed 136 proof excludes higher dimensions. At ambient 126, the
profile (22,22,22,22,22,15) has negative margin, so the same scalar weights
do not prove 125.

The [separate review](fresh_selective_distance11_dimension126_20261017_review.md)
reconstructs the critical arithmetic with Python integers and fractions
and reviews the geometric implications. Internal verification is distinct
from external refereeing and proof-assistant formalization. The earlier
[163 proof-status note](dimension163_proof_status_20261017.md) continues
to apply: that argument was complete and internally checked; the large
gap to examples was not an unproved step in it.

## Compact arithmetic: a new corner cohort

The [raw-25 calculation](fresh_raw25_actual_corner_compatibility_20261017.md)
reconstructs 55 positive face classes, including 36 newly admitted rows
with unique saved provenances. Their 440 finite prism corners give 6,260
literal orientations in 296 rank-13 vertex forests. Every external-wall
product at a shared finite vertex is determined by the full contact
covectors and the common inward orientation.

Of 149,510 literal-star pairs, 141,850 are excluded by shared-face
incidence, acuteness, the Coxeter angle alphabet, or the fixed-ground-field
arithmetic conditions. The arithmetic gates here are specifically for
Q(sqrt(5)); they are not assertions for unrestricted Coxeter polytopes or
arbitrary arithmetic ground fields.

Explicit proper colorings and capacity-one clique certificates give 500
necessary incidence inequalities. Twenty-eight violate the current fixed
witness. Raising 25 vertex allowances, with every face mass unchanged,
repairs it from 12.43082950 to

    67685455440946980000000000000000 /
    5466425626097138727892456269673
    = 12.38203170968089...

The value is still above the target 12. No optimizer was run, and this is
neither an optimum of the enlarged model nor a new arithmetic dimension
bound. All inherited inequalities remain satisfied by monotonicity. The
remaining merged rows have not acquired constraints from unexamined
provenances.

Compact arithmetic remains bounded by 21 in the repository's completed
internal proof. The noncompact rational case and its known sharp bound
are separate. The unrestricted minimum-cusp-factor-rank-9 class remains
bounded by 29.

## Construction: removing a forbidden adjacency

The [first five-cap scout](fresh_broken_cap_adjacency_20261017.md)
subdivides the previously forbidden A-R edge by a new cap X. Its disk
has 10 vertices, 22 edges, and 13 triangles. The old A-R finite-corner
obstruction no longer applies. A bounded search over 22 samples of the
continuous A curve, 94 earlier R seeds, and 50,072 X norm equations finds
no coherent bridge. This is not an all-parameter exclusion.

The [event continuation](fresh_broken_cap_adjacency_events_20261017.md)
finds 57 coherent X roots using contacts with old walls 8, 30, and 63.
Allowing the R contact with wall 44 to vary gives 4,476 norm equations,
but no required R bridge. Most otherwise coherent roots have a dotted
contact with 44, so they cannot provide the needed finite corners.

A [different five-cap placement](fresh_broken_cap_adjacency_split_patch_20261017.md)
makes X meet isolated old wall 8 and cover the old edge 30-63. Two
quadrilateral disks, [7,8,X,30] and [8,44,63,X], are then closed by separate
cap pairs. This removes the requirement that one cap cover both old edges
at 44. The bounded search tests 41 X seeds and 5,665 norm equations, with
no complete pair.

The [exact successor gate](fresh_broken_cap_adjacency_gate_20261017.md)
excludes that last X placement over all seven A2 profiles and every finite
label. Its 148 X recipes and 344 A-X cases are covered by 691 certified
interval leaves. It uses the three mandatory finite corners involving
S6+7+8+A, S6+8+A+X, and S6+30+63+X. The source reuses the earlier interval
engine with the new literal frame; it is not an independent arithmetic
implementation. The original zero-boundary-arc X with dotted contacts to
30 and 63 remains outside this exact exclusion.

The construction control is unchanged: 63 walls, 199 finite vertices,
10 target-face open ends, and certified infinite volume. No new
finite-volume example or construction dimension has been established.

## Verification and next work

The [incremental auditor](../scripts/audit_distance11_and_raw25_corner_frontier_20261017.py)
and [saved audit](distance11_and_raw25_corner_frontier_20261017_audit.json)
rehash the frozen predecessor before adding this checkpoint. They verify
the new scalar composition, the fixed-witness repair, and the construction
interval-tree counts. In particular, they match all 200 ports, 1,010
overlap classifications, 114 preferred corners, and 1,095 scalar profiles
between the primary and independent outputs. All 500 raw-25 loads and
the 344 construction interval trees are checked. The reports distinguish
exact exclusions from bounded numerical searches. Historical files are
preserved.

The next arithmetic step is compatibility between the raw-23 and raw-25
families at the same actual finite vertex; their separate capacities do
not yet impose all mixed constraints. The next distance reduction needs
a fresh classification of supporting extensions, because longer shared
paths can admit a positive orthogonal residual. For construction, the
zero-boundary-arc X remains a real open parameter domain; the latest
exact exclusion applies to the alternative X placement only.
