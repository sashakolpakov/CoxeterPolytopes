# Star case0: opposite single contacts on the last block

**The specified opposite-single family has no compact H9 realization when the last internal label is at least4.** The first three blocks retain all six leaves with own contacts (3,3). The last two leaves, 11 and17, are right to roots12 and13 respectively and have arbitrary finite contact labels at least4 to roots13 and12. All four internal labels are at least4; other finite labels and dotted distances remain free. The checked swap of roots12 and13 gives the other orientation.

The [exact producer](../scripts/fresh_h9_star_zero_20261014_opposite_single.py) and [certificate](fresh_h9_star_zero_20261014_opposite_single.json) close this continuous family with **832 positive polynomial rectangles and64 exact actual-vertex failures**. The scope includes the saved [second search survivor](fresh_h9_star_zero_20261014_family_search_batch_002.json), but not the full star incidence. The last-internal-label3 case is treated by a [separate argument](fresh_h9_star_zero_20261014_opposite_three.md).

Write the positive eight-root block as
\[
S=(3,4;6,7;9,10;12,13).
\]
The actual four-root supports force every cross-block contact right: a nonright crossing would contain the relaxed affine path4–3–4, and increasing any bond cannot restore positive definiteness. The eight roots therefore have a positive direct-sum Gram. A literal application of the reviewed two-dense or strong-contact theorem excludes any internal label at least6. For the first three blocks the own (3,3) missing triples already force internal label at least4. The final internal≥4 hypothesis is explicit here; it does not follow merely from its opposite-single missing triples.

There remain16 internal-label choices in {4,5}. If the last internal label is4, its own missing triple (2,4,m) requires m≥5. If it is5, (2,5,m) requires m≥4. Thus, for the two nonright own cosines α andβ, the complete continuous rectangles are
\[
(\alpha,\beta)\in[\cos(\pi/5),1]^2\quad(n=4),\qquad
(\alpha,\beta)\in[\cos(\pi/4),1]^2\quad(n=5).
\]
Including the upper endpoint1 is a relaxation and avoids any upper label bound.

Actual spherical supports retain one own root and all six foreign roots for each anchor. Its own projection cost is at least1/2; any foreign nonright contact contributes at least1/2 because every foreign inverse diagonal is at least2. Hence both anchors are right to all foreign roots. Their mutual label is2 or3: the actual three-root support with root13 has the nonpositive4–4 path if the mutual label were at least4.

For the third normal14, the analogous budget leaves at most one foreign π/3 contact, giving seven possible projections. Four fresh exact vertex witnesses, one for each third–anchor pair and last-internal choice, exclude label4 at the minimum own cosines. They remain nonpositive as either own cosine or any first-three internal label increases. Thus its two anchor contacts also have labels2 or3. These witnesses are recomputed for the new right contacts; no monotonicity from a previously nonright contact is assumed.

Let c be the last internal cosine, \(D=(1-c^2)^{-1}\), and \(C=cD\). In the residual Gram of anchors11,17 and third14, the two anchor diagonals and their mutual entry are
\[
1-D\alpha^2,\qquad1-D\beta^2,\qquad-\delta-C\alpha\beta.
\]
The other two entries are \(-u-w_a\alpha\) and \(-v-w_b\beta\), where δ,u,v are each0 or1/2. The third diagonal is1 minus its exact projection cost. The saved certificate records every resulting bivariate determinant coefficient.

The16 block choices, seven third projections and eight contact choices give896 cases. In832 cases all nine Bernstein coefficients on the stated closed rectangle are strictly positive: **7,488 positive coefficients** in total. The remaining64 cases have an exact nonpositive leading pivot in an actual finite-vertex Gram at the minimum own labels. Since this normalized Gram has nonpositive off-diagonal entries, increasing the own cosines preserves failure of positive definiteness. None of the64 cases is silently discarded on a numerical test.

A positive residual determinant makes the eleven-normal principal Gram have rank11, contradicting the ten-dimensional Lorentzian ambient space of H9. The alternative vertex witnesses contradict compactness directly. These alternatives exhaust the stated family. Generation and a saved replay pass; the separate review records the independent arithmetic scope.

The [independent review](fresh_h9_star_zero_20261014_opposite_single_review.md) passes generation and saved replay; it records the separately reconstructed arithmetic and literal-domain checks.
