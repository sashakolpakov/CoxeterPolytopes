# Independent review of the selective distance-10 proof

**FINAL / PASS.** The mathematical/source review, independent exact generation, and fresh saved replay all pass. No gap remains identified in the stated unrestricted dimension-116 and compact-facet dimension-87 conclusions.

The [independent checker](../scripts/fresh_selective_distance10_dimension116_20261017_review.py) and its [certificate](fresh_selective_distance10_dimension116_20261017_review.json) review the [primary proof](fresh_selective_distance10_dimension116_20261017.md), [source](../scripts/fresh_selective_distance10_dimension116_20261017.py), and [data](fresh_selective_distance10_dimension116_20261017.json). The checker imports no producer arithmetic, SymPy, numerical linear algebra, or optimizer. Its rational calculations use Python `Fraction`; the small-order supplement uses separately implemented exact arithmetic in Q(sqrt(5)).

## Independent finite reconstruction

Starting from the frozen [26-type high-rank catalogue](fresh_high_rank_quasilanner_census_20261017.md), the checker independently reconstructs all 384 marked contexts and the complete relevant ports, joins, partial identifications, and corner decompositions:

| Check | Exact result |
|---|---:|
| High-rank candidate / inverse-admissible ports | 527 / 161 |
| Ordered disjoint joins / white obstructions | 1,478 / 363 |
| Multiple-edge subsets, all rejected | 120 |
| Distinct single-join inertia types, all with two negative directions | 35 |
| Candidate / consistent overlaps | 22,881 / 4,862 |
| Inverse / actual-edge rejections | 4,296 / 254 |
| Surviving overlap contexts | 312 |
| Adjacent corner decompositions | 1,248 |
| Extra cross-edge shortcut / actual-edge rejections | 21,393 / 168 |
| Saved nonnegative producer witnesses independently evaluated | 905 |
| Reduced cusp profiles | 1,015 |

Each rational Gram is constructed afresh from the literal diagrams. For diagrams with label 4, the checker recomputes a positive diagonal sqrt(2) congruence on the whole joined graph, rather than joining independently chosen root-length gauges. A separate rational elimination supplies signatures, inverse diagonals, null vectors, and positive or nonnegative-square tests. Every saved nonnegative witness is evaluated against the independently reconstructed principal Gram, with its gauge and square checked exactly.

The overlap enumeration is injective and checks all prescribed contacts on the common nodes. The independent preferred-corner records agree with the producer, including contact domains, finite/affine components, residual norm budgets, candidate-pair Cauchy inequalities, and all extra cross-edge checks. The four adjacent corners are tested even when their two marked roots lie in different components; only a connected preferred corner is used for payment.

## Essential small-order supplement

Distance 10 requires an order-5/radius-3 disjoint supplement that was unnecessary at distance 11. A small core can join a radius-7 port at distance 11, so discarding all cores of order at most 5 would leave a gap.

The final primary proof supplies the complete all-label reduction: the three white nodes form a chain, the two black nodes have three possible tree attachments, and the only extra cycle is a black-end triangle ruled out in a proper four-node subdiagram. Every bond lies in a connected proper four-node positive-semidefinite diagram. A label at least 6 produces an indefinite connected triple, or an affine equality whose connected fourth-node extension is indefinite. Thus labels 3, 4, 5 exhaust this domain. This is a classification argument, not an assumed numerical cutoff.

The independent checker computes every principal tree determinant by the matching-polynomial recurrence, which depends only on squared bond magnitudes 1, 2, and (3+sqrt(5))/2. All 243 recipes agree with the producer: 209 have an indefinite proper principal subdiagram, 14 are positive semidefinite, four have a nonspherical black deletion, 12 fail the port test, and four survive. Their normalized far-port inverse pairs are (4/5,3/4), (3/4,4/5), (1,4/5), and (4/5,1).

Each surviving low core has one possible joining edge to L25/node 9, in three high marked contexts and either orientation. A second, independent leaf-Schur algorithm computes inertia (13,2,0) for all four joined trees, excluding all 24 ordered cases. Its zero-leaf rule uses a two-node indefinite pivot and requires no irrational matrix entries. Exact signs in Q(sqrt(5)) are decided by rational square comparison.

## Geometric and completeness audit

The actual-edge tests require positive definiteness of W plus a single black root, including edges with ideal endpoints. This makes the inverse cutoff strict at 2. Cross-core endpoints are white under the assumed absence of a short marked path. The three-white-node test forces every possible disjoint cross bond to be simple. All subsets of the remaining possible edges are then exhausted: a single edge has two negative directions, and every multiple-edge subset has a nonnegative actual-edge obstruction. Consequently additional acute bonds cannot repair the excluded cases.

For overlaps, every unspecified bond between exclusive core parts either gives the required short path or has a nonnegative actual-edge obstruction. The latter remains valid when further nonpositive contacts are added. This fixes the core Gram before its extension arguments are used.

A radical of an indefinite realized core is an actual linear relation: the corresponding zero-square vector lies perpendicular to a timelike vector and therefore vanishes. The independently checked preferred contact domain excludes zero-radical coordinates. Any nonzero acute contact vector must meet both signs of the relation, yielding the checked short path.

For nonsingular cores, supporting whites lie in a Euclidean orthogonal complement after projection. A new root's core projection and residual satisfy the exact norm budget. The positive inverse end blocks exclude multiple contacts and every larger Coxeter label. Candidate-pair residual Cauchy inequalities either contradict acuteness or force a short marked path, so at most one white can contact the core.

The [separate supporting-extension lemma](fresh_distance10_extension_20261017.md) controls the entire indirectly attached component. The other supporting whites form a positive-definite principal block E. Projection cost is at least 1/2 per contacted component, at least 2/3 for a connected component of order at least two, and at least 3/4 for order at least three. Residual squares 0, 1/2, and 2/3 therefore allow at most zero, one, and two additional roots. Arbitrarily long indirect tails are excluded analytically; no rank cutoff is used. All 312 selected marked components remain of order at most 17.

At ideal corners, the inherited factor-fiber theorem places the whole varying affine factor in W plus the adjacent pair. Extra cusp walls outside that face star cannot enlarge the paid marked component. This is essential for nonsimple cusps and does not assume that an induced face is itself a Coxeter polytope.

## Curvature and scope

The weight pays distance at most 10 or connected component order at most 17. The inherited triangle argument gives two short paid corners; the new classification gives one paid quadrilateral corner. The finite cost is at most 10n, and the normalized ideal cost is at most sum min(10,r_i/2).

The checker independently evaluates every cap-20 reduced profile for link dimensions 116 through 125, using both integer convolution and the inclusion-exclusion absolute-moment formula from the [shape-curvature framework](fresh_cusp_link_shape_curvature_20261017.md). All 1,015 margins are positive. The minimum is 145234/4084101; the tested boundary profile at ambient dimension 116 has margin -673055/16336404. The latter is a limitation of this scalar payment, not a realized polytope. The compact-facet parity polynomials have coefficients (12,176,4) and (96,180,4).

The [completed dimension-126 proof](fresh_selective_distance11_dimension126_20261017.md) supplies the excluded higher-dimensional tail. With the written geometric and all-rank curvature inputs, these checks establish the stated dimension-116 and compact-facet dimension-87 consequences. This is a complete internal mathematical/source and independent exact-arithmetic review of the new step. It does not independently reconstruct published finite classifications or every older theorem. It is not external refereeing, proof-assistant formalization, a new construction, or a priority claim.

Final generation 21363 and fresh saved replay 80114 both terminated successfully. The review binds the frozen primary source, certificate, and proof, and preserves the completed 126 predecessor. All 11 direct review bindings, all 23 producer bindings, source AST, and nine local report links pass. The review source and data are frozen; no producer files were edited. The independent low-order implementation and all 905 saved-witness evaluations are included in the fresh replay.
