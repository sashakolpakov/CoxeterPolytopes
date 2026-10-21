This certificate identifies the new quadrilateral branches that must be checked before using the proposed adjacent-3 discount. It screens every one of the 8,298 reviewed minimal marked unions and every accepted first-white extension in the complete contact census. It does not assert a new dimension bound or a complete attachment closure.

Write (B) for the sum of the four existing baseline weights: a connected marked pair has weight one at graph distance at most two, with the additional H4 weight for its 5-edge or its end pair at distance three. The proposed change subtracts δ only when the two black corner roots are joined by a 3-edge in an entire spherical component of rank at least six. Here the distance-two discount is zero and (0\leq\delta<2/3). Every additional compensation weight is nonnegative. The retained facet root is white or an external spectator; it is never one of the black pair.

An existing path in a finite Coxeter forest cannot change when white vertices are added. An H4 component is already an entire finite component. Thus the baseline cannot decrease under finite extension. The number of discounted corners is at most (B), so (B\geq3) remains safe: (B-\delta B\geq1), with strict inequality at (B=3). At (B=2), a new obligation requires two discounted corners. At (B=1), it requires its unique paid corner to be discounted. A baseline-zero core has no adjacent black 3-edge; adding white roots cannot create such an edge. Consequently the new discount does not affect that old branch, although the current compensation coefficients must still pay its stopping states.

The extension screen retains small A, B and D components which may grow to rank six. It does not assume that current rank below six is safe. H3 can only grow to H4; H4 and F4 have no larger connected spherical extension, so none can acquire the discount. The complete finite-type classification justifies this distinction. All relevant types are retained, including E6, E7 and E8; there is no A-tail-only restriction.

The exact screen leaves 1,216 potentially affected minimal classes. Their inherited contact statuses are:

| Status | Classes |
|---|---:|
| Finite, with accepted first-white contacts | 789 |
| Symbolic-dotted contact domain | 178 |
| Sealed singular finite core | 141 |
| Sealed nonsingular core, already rank-excluded | 93 |
| Already excluded by an affine principal subdiagram | 15 |

Sealing fixes the complete marked components, so only **35** of the 141 sealed singular classes actually contain the discounted configuration: 30 have baseline one and one eligible corner, and five have baseline two and two eligible corners. Their literal rows and all four corner forests are exported in `immediate_sealed_weight_obligations`. Every allowed retained-facet mark and common spectator context must be tested before treating these rows as paid.

Among accepted first extensions of the 1,216 classes, 4,497 finite and 2,351 symbolic-dotted literal extensions remain potentially affected. Another 1,021 finite and 1,230 dotted extensions become permanently safe by the baseline argument. These are literal counts, not isomorphism-class counts. The surviving rows are exported in `outgoing_roots` for a fresh closure queue. A rank-21 corner permits at most 23 vertices in its full quadrilateral star, but that rank bound does not by itself complete either queue. In particular, the previously completed baseline-zero finite and dotted queues must not be presented as coverage of these new branches.

The [script](../scripts/fresh_sqrt5_adjacent_discount_frontier_20261016.py) and [result](fresh_sqrt5_adjacent_discount_frontier_20261016.json) bind the [complete minimal-union census](fresh_sqrt5_all_quad_unions_20261015.md), its marked quotient and the [complete one-white contact census](fresh_sqrt5_all_quad_contacts_20261015.md), including their independent review data. Graph decomposition uses the frozen separate graph-review helper. This package replays graph transitions and provenance; it imports the prior arithmetic decisions rather than duplicating their exact field calculations. Generation and fresh saved replay passed. No numerical candidate payment is promoted to an exact conclusion here.
