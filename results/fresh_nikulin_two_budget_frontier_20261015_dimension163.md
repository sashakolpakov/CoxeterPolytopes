# Reviewed research proof candidate: dimension at most 163

**Candidate theorem.** Every finite-sided, finite-volume real hyperbolic
Coxeter polytope has dimension **n≤163**, without arithmeticity, compactness,
a simple-cusp hypothesis, or a facet-count restriction.

The new [ordered-budget proof](fresh_nikulin_two_budget_frontier_20261015.md)
and [independent review](fresh_nikulin_two_budget_frontier_20261015_review.md)
sharpen the frozen [164 candidate](fresh_nikulin_defect_frontier_20261015_dimension164.md).
This [case binder](../scripts/fresh_nikulin_two_budget_frontier_20261015_dimension163.py)
and its [certificate](fresh_nikulin_two_budget_frontier_20261015_dimension163.json)
check the local outputs, complete geometric case split and recursive source
hashes. They preserve every preceding artifact. Historical pending-review
wording in the frozen primary report is superseded by the independent review
bound here. External validation, novelty, and current-best literature status
remain unestablished.

If the chosen narrow facet has an ideal vertex, its affine factors have
total rank N=n−1. Use the inherited actual framing walls, common height ratio
M≥1, prescribed omissions, and separate Schur budgets C₀≤2 and C₁≤2. For
M≤21/20 the new all-rank continuous envelope proves total target loss at
most four. When unit charge fails, select any failing factor. An ordinary
E8 factor leaves endpoint budgets at most (10/9,4/3), whose residual loss
is at most20/9; its own target is16/9. An E7-type factor leaves (1,3/2)
or (3/2,1), whose residual loss is at most22/9; its own target is14/9.
Both sums are four, and additional failures are included in the residual
product. With no failing factor the two total budgets already give four.

The independently checked envelope covers all 91 finite models and all
3,630 continuous support domains, plus analytic A/B/C/D tails at every
remaining rank. It uses no finite-angle cutoff. The two same-row E7
possibilities are excluded by separate exact convex certificates, and the
28-record envelope is packed with both endpoint budgets in their stated
order. The review uses independently reconstructed physical Grams,
Fraction convex-certificate checks and a separate packing algorithm.

The targets sum to (2N−2)/9, including the sole prescribed avoid-node loss.
Thus the near-range Lorentz pairing is at least M(N−37)/9. For M≥21/20 the
unchanged reviewed four-thirds estimate, with the same omissions, gives at
least M(N−43)/9. At N=163+t these lower bounds are respectively at least
14+t/9 and14+7t/60. They contradict the strict narrow-facet inequality
−(x,y)<14 for norm-two walls for every t≥0. This branch excludes n≥164.

If the chosen facet has no ideal vertex, finite-sidedness and finite volume
make it compact. The frozen
[compact-facet proof](fresh_nikulin_narrow_20261015_compact_facet.md), with
its [hostile review](fresh_dimension225_hostile_20261015.md), gives n≤143.
These exhaustive cases yield the candidate n≤163. The published narrowness
theorem and the reviewed incidence, common-height and compact-facet arguments
remain explicit inherited inputs; the local arithmetic certificate alone is
not a proof of those geometric statements.

Both scalar estimates at N=162 are below14, so this package claims no
dimension162 bound. Its loss-four local witness concerns the relaxed charge
domain only, not existence of an ambient hyperbolic polytope.

Generation and saved replay **PASS**. The binder checks all recursively
referenced manifests, the reviewed threshold at N=163, and the preserved
compact-facet143 branch.

```sh
python3 scripts/fresh_nikulin_two_budget_frontier_20261015_dimension163.py
python3 scripts/fresh_nikulin_two_budget_frontier_20261015_dimension163.py --verify
```
