# Star case0: complete own-contact profiles at last internal label4 or5

**This entire contact family has no compact H9 Coxeter realization.** The first three blocks `(3,4)`, `(6,7)`, `(9,10)` retain both leaves with own contacts `(3,3)`. Their internal labels may be arbitrary. The last block `(12,13)` has internal label4 or5, and its two leaves11,17 have arbitrary finite own labels. Every unlisted meeting label and dotted distance remains free.

The [exact checker](../scripts/fresh_h9_star_lastblock_20261014.py) and [certificate](fresh_h9_star_lastblock_20261014.json) enumerate290 ordered own-profile pairs, certify168 fresh contact bounds and compute77,952 residual determinants. All32 zero determinants have exact nonpositive actual-vertex witnesses. This excludes the [saved51-visit state](fresh_h9_star_zero_20261014_opposite_search_batch_003.json), while leaving the full incidence open.

The source binds the literal378-vertex incidence and the earlier independently reviewed [all-own-three](fresh_h9_star_zero_20261014_own_three_review.md) and [opposite-single](fresh_h9_star_zero_20261014_opposite_single_review.md) families. It does not modify those frozen certificates.

## The continuous high-own-label reduction

The actual four-root supports force all cross-block contacts right: any nonright crossing contains the nonpositive relaxed path4–3–4. Thus the eight block roots have a positive direct-sum Gram. The reviewed two-dense gate excludes a label at least6 in any of the first three blocks, using only their fixed own contacts. These three internal labels reduce to4 or5.

Suppose one last anchor, denoted x, has an own label at least6 to a root a. Let its positive cosine be α, its other own cosine be β≥0, and let c be the last internal cosine. Put
\[
D=(1-c^2)^{-1}\ge2,\qquad
X=D(\alpha^2+\beta^2+2c\alpha\beta).
\]
Here α²≥3/4 and X≥3D/4. An actual finite support retaining a and all six foreign block roots has own cost at least3/4. Any foreign nonright contact contributes at least1/2, so x is right to every foreign block root.

The actual triangle consisting of x, a and leaf14 contains the label≥6 bond. Either additional nonright bond produces a nonpositive three-node path, since `1−3/4−1/4=0`. Consequently x is right to14 and14 is right to a. The saved certificate checks all four literal choices of x and a.

The retained-own-root budget for14 permits at most one foreign contact, necessarily π/3. Its own first-block projection has square at least `1+sqrt2/2`. If14 has no last-block contact, its residual is negative and orthogonal to the negative residual of x, immediately giving two negative directions.

Otherwise its sole last-block contact is π/3 at the other root. Write Q for its projection square outside the last block; thus Q≥`1+sqrt2/2`, and the total projection square is `Q+D/4`. The residual two-by-two determinant is exactly
\[
\begin{aligned}
&(X-1)(Q-1+D/4)-(Dc\alpha+D\beta)^2/4\\
&\qquad=(X-1)(Q-1)+D(\alpha^2-1)/4\\
&\qquad\ge (3D/4-1)\sqrt2/2-D/16
\ge\sqrt2/4-1/8>0.
\end{aligned}
\]
The identity uses `D(1−c²)=1`. Both residual diagonal entries are negative, so this is a negative definite plane, impossible in Lorentz signature. This proof allows an arbitrary second own contact β; it does not assume an opposite-single profile. All four last own labels therefore lie in `{2,3,4,5}`.

## The finite census

The two own missing triples must be Lannér, so their reciprocal-label sums are strictly less than1. This leaves11 ordered profiles for last internal4 and13 for last internal5, hence `11²+13²=290` ordered pairs.

If the two profiles have no common nonright root, they are opposite-single profiles already excluded by the frozen all-label theorem. Otherwise an actual three-root triangle forces the two anchors mutually right: testing a nonright mutual contact gives nonpositive quadratic form on the all-ones vector. At either own root, the necessary triangle determinant is `1−cos²(a)−cos²(b)>0`; nonpositive cases are recorded and removed. The pair `(3,3)/(3,3)` is inherited from the all-own-three theorem.

For any remaining anchor profile other than `(3,3)`, one own label is at least4. The retained-root budget gives own cost at least1/2, so every foreign contact is right. A `(3,3)` profile permits either no foreign contact or one π/3 contact to one of six foreign roots. Leaf14 has the same seven conservative possibilities. Both retained-root orientations are checked, so no increase on the other own root is overlooked.

For each remaining pair, two fresh exact actual-vertex witnesses exclude a label4 between14 and the two anchors. The witnesses use the minimum first-three internal labels4 and the literal enumerated own profiles; every unassigned contact is right. Increasing an internal cosine or any omitted contact cannot restore positive definiteness of a matrix with nonpositive off-diagonal entries. Thus the two remaining labels are2 or3. There are168 such witnesses, each saved with exact leading LDL pivots.

Projection off the positive eight-root span leaves a three-by-three residual Gram for11,17,14. The complete16 choices of internal labels4/5 and all permitted projections/remaining labels give77,952 cases. A nonzero residual determinant makes the eleven-normal principal Gram have rank11, impossible in the ten-dimensional ambient space of H9. There are32 zero cases; each is rejected by a saved nonpositive actual-vertex LDL pivot. No zero or degenerate residual case is omitted.

All arithmetic is exact in `Q(sqrt2,sqrt5)`. The graph classifier is used only to locate candidate bad vertices: every rejection used by the certificate has a directly recomputed exact field LDL witness. No numerical fitting or finite bound on the original unrestricted own labels is assumed. The theorem addresses last internal labels4 and5 only; other last-internal regimes and changed first-three own profiles remain outside this certificate.

## Resumed search and preserved frontier

The [gate binder](../scripts/fresh_h9_star_lastblock_20261014_gates.py) adds the two proven last-internal values to the154 frozen conditions. The [new search](../scripts/fresh_h9_star_lastblock_20261014_search.py) resumes the exact51-visit cursor without repeating its explored prefix. It validates the saved paths and next-child cursors and preserves every pending sibling.

[Batch004](fresh_h9_star_lastblock_20261014_search_batch_004.json) uses52 additional visits and stops at a full mandatory-label/bridge survivor, for103 cumulative visits and nine pending frames. Its first three blocks are B2 with all own contacts3. The last block has an unrestricted finite label `m>=7`, and both leaves11,17 have own profiles `(2,3)`: both are right to12 and have π/3 contact to13. This is a shared-root high-label regime, outside the opposite-single and last-internal4/5 theorems. Category7 remains the entire unbounded finite-label tail, not a chosen angle or finite cutoff.

The saved state is not a full metric candidate: omitted meeting labels, dotted distances, rank and Lorentz signature still need completion. The incidence remains open. A useful next step is a continuous Schur analysis of this shared-root high block; its inverse diagonal is `1/(1−cos²(pi/m))`, and the projection budgets must be re-established for the decreased own contacts.

Primary generation and saved replay pass. The gate binder and batch004 also have saved replay commands:

```sh
python scripts/fresh_h9_star_lastblock_20261014.py --verify
python scripts/fresh_h9_star_lastblock_20261014_gates.py --verify
python scripts/fresh_h9_star_lastblock_20261014_search.py --batch 4 --verify
```

The primary imports the frozen exact algebraic-field helper; its exact arithmetic is not described as an independent implementation. The separately archived review states its own reconstruction scope.

The [independent review](fresh_h9_star_lastblock_20261014_review.md) passes generation and saved replay. Its standalone Fraction implementation imports no producer or project arithmetic: it reconstructs all290 profile domains, six budgets, four high-own triangles,24 cross guards,168 fresh contact bounds,77,952 determinants and all32 zero-case vertex failures, and verifies the continuous high-own identity and lower bound. This is a separate arithmetic implementation, not merely a rerun of the primary engine.
