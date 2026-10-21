# M6 correction: two endpoint mirrors, one simple wall

20 October 2026. Independent paper and full source-delta review PASS,
without execution, for the corrected inputs:

* preparation SHA256
  `2afd545c3c3166cd07381949b972fb012280b18ed673a090025b6aac30537dd2`;
* `scripts/check_golden_dyadic_m6_actual_edge_20261020.py`, SHA256
  `9e48ab80a337f43925a1be7ee6761b02da5eaa0c60df6dfb649bc640ec5e6949`.

The unchanged historical pre-execution review
`review_golden_dyadic_m6_edge_preexecution_batch_guard_20261020.md`,
SHA256 `10db918c67a420e6c5b4d2f248d75a86df39613ae3c9461e274bf053f093e1ab`,
was too strong about the number of accepted mirrors. Its uniqueness
assertion is explicitly superseded here. The bounded execution caught
that incorrect expected count; no successful old certificate is being
inferred from that review.

The additional mirror is genuine and must remain in the census:

    g=2h+β1,
    g²=4+2+4(−1)=2,
    (g.β1,...,g.β5)=(0,−1,0,0,0).

Thus it is the q=2,k=2 D6 profile with R=4. Its time and spatial sum
are exactly 2aH and 2ZH; its edge collision parameter is tH, not an
earlier one. Its coordinates are the ones now asserted by the code.

It is not an actual simple opposite wall. The initial full chamber
lies in h<=0 because h is an integral arrangement mirror, without
assuming h already a facet. It also lies in the actual β1<=0 halfspace.
Consequently g=0 on the chamber forces h=β1=0. These normals are
independent, with two-by-two Gram determinant 1, so this is at most
a codimension-two stratum. A first simple endpoint wall must be one
of the complete accepted set, and g is excluded by this geometric
positive-cone argument. Therefore h is the unique new SIMPLE wall.
The extra g is a positive root in the same B6 endpoint chamber; it
does not replace that chamber by a D6 germ.

I checked the full corrected code. The box, six R loops, nine profiles,
square-root routine, glue and actual-old-wall filters are unchanged.
The code now requires exactly two accepted records with the exact
profile set {(1,1),(2,2)}, asserts both full vectors, and verifies the
linear identity, norm, doubled time/Z and exact tie. Both records remain
in the certificate, with h separately selected as simple and the positive
coefficients recorded. No candidate was removed by a narrower filter.
The earlier reversal/Gram/future and known-double inputs are unchanged.

No new flaw was found. This is a corrected PRE-EXECUTION source verdict;
corrected generation and an independent saved replay are still separate
requirements. The published compact-example, actual-wall, sealed-face
and full-germ geometric inputs retain their previously stated scope.
