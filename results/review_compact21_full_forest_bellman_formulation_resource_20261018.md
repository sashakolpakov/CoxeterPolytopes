# Independent mathematical review of the complete forest epigraph

The resource/research agent read the entire root note
`compact21_full_forest_bellman_formulation_20261018.md` and finds its abstract
extended formulation correct, subject to the stated literal additivity premise.
This review is separate from implementation execution and any numerical result.

The fixed entire A1 count and capped other three hub counts determine all
component coefficients. Normalized crosses assign their vertex sum to the
destination: A1, I2(5) and H3 contributions contribute one rooted bonus per
destination root. For A4, each of its two source-root orbits contributes one
half of its rooted bonus at each destination root. A4-to-A4 bonuses contribute
one eighth of the sixteen ordered source-root pairs per A4 component. The
implementation must compare this complete coefficient convention with the
frozen literal forest support; it does not suffice to compare numerical costs
at one sparse vector.

The untracked list contains rank-two and rank-three items, so all residual
ranks except one are reachable. Acyclic Bellman inequalities in increasing
rank bound every multiset, including repeated components, and their minimal
solution is the exact maximum at each rank. Sharing a state family between
different exact hub tuples is valid because their component coefficients are
identical after the specified caps, while their fixed hub costs remain separate.
U(0) is the constant zero; nonnegative potentials do not restrict the intended
maxima because all old component weights are nonnegative.

The exceptional ordinary budgets occur in contexts with zero other hubs, so
each removed ordinary final budget corresponds to one exact hub tuple. A
reserved off-target item excludes the removed target and permits an arbitrary
remaining multiset through its Bellman state. Pure H4/E8 rank-20 solutions are
5H4, 3H4+E8 and H4+2E8; pure D5/E8 rank-21 has only D5+2E8. The explicit rows
therefore cover precisely the cases not covered by off-target items. Bellman
states at these terminal ranks may remain without imposing an obsolete final
upper bound. The four positive source adjustments require their stronger
explicit rows; retaining their old ordinary budgets is redundant and harmless
because t is nonnegative.

No unproved incidence completeness beyond the separately reviewed D5 and H4
packets is introduced by this algebraic formulation. It proves equivalence to
their conservative full finite-forest model, not to geometric realizability of
every abstract forest. Native-face completeness and the final dimension
composition remain distinct theoretical obligations.

Implementation plan: assert all 510 old support forest rows, including both
D5 states and all strengthened source rows, have exactly reconstructed
coefficients and right-hand sides. Then discard numerical auxiliary values,
recompute the minimal Bellman states from the rounded rational bonus vector,
and verify every exact hub budget, off-target budget and explicit exception.
A separate full separator and independent coverage review can subsequently
cross-check the resulting candidate. No solved bound or Lean verification is
asserted by this review.
