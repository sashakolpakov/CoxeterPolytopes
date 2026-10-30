# Mixed right-angled corner package: proof map and reproduction

For current reading, see the [current sharp right-angled manuscript](../arxiv/sharp_dimension_bounds_8_and_4/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

The supporting argument
collects the universal weighted-face, low-facet rigidity, complete-corner
and opposition-moment arguments. None assumes arithmeticity or a restricted
four-facet spectrum. The dimension ceiling remains 12.

|Claim|Written proof and verification boundary|
|---|---|
|Finite/cusp ratios and forced P4 in dimensions 11/12|[Weighted packet](../results/right_angled_mixed_weighted_triangles_research_20261018.md); exact DS/coefficient/model checks and [root review](../results/review_right_angled_mixed_weighted_triangles_20261018.md)|
|P4 uniquely attains facet-plus-cusp sum 15; P5 uniquely has 16 facets|[Rigidity packet](../results/right_angled_mixed_low_facet_rigidity_research_20261018.md); all 4368 P5 active-wall choices and [independent review](../results/right_angled_mixed_low_facet_rigidity_independent_review_20261018.md)|
|Complete five-corner implies P5 or at least 20 facets; otherwise 11L<=V|[Complete-corner proof](../results/right_angled_mixed_complete_five_corner_research_20261018.md); actual cusp-width and unique-predecessor arguments, exact double attaining 20 facets and L/V=2/22|
|No closed complete-corner opposition component for n>=9|[Opposition proof](../results/right_angled_mixed_opposition_moments_research_20261018.md); invariant h and finite neighbor-sum contradiction, reviewed written geometry|
|Universal n=12 bound 2B<184.372f|Same opposition packet; exact centered moments, positive reciprocal remainder and edge-ratio summation|
|All linear clique objectives reduce to 77 multipartite vectors|Same packet; written symmetrization proof with arbitrary signed coefficients; exact vector generation|
|Current scalar inequalities still have a numerical solution|Exact DS, incidence, polygon, moment-cap and integer-scaled Macaulay replay; the assignment is not a geometric polytope and its saturation violates complete-corner propagation|

The published inputs are the usual Lorentz normal and finite/cusp link
facts, right-angle inheritance, truncation and Dehn--Sommerville, the
eight-dimensional 82-facet lower bound recorded in Alexandrov, and the
four-dimensional 10-facet and facet-plus-cusp 15 bounds of
Potyagailo--Vinberg. Their equality classifications, actual-face closure,
cusp-width restrictions, opposition invariant and moment argument are
internal written proofs, not relabelled external assumptions.

Run from the repository root:

```sh
python3 -B scripts/check_mixed_corner_package_20261018.py --build
```

This replays all five small standard-library checkers sequentially, then
builds the single paper. It requires the common resource gate and exclusive
computation lock. It retains the PDF and a compact execution ledger in
[mixed_corner_package_audit_20261018.json](mixed_corner_package_audit_20261018.json)
and removes temporary TeX files even on failure. No Lean process runs.

The initial weighted/rigidity TeX received a full mathematical review by
the right-angled lane. The complete-corner and opposition additions were
derived jointly and checked in both lanes; their arithmetic has separate
small exact replays. This is internal mathematical verification. Full
Lean verification of the new geometric arguments remains pending.

Final serial package execution passed in chunks `e60253`/`bb25f0`: all
five checkers and the PDF build passed, with zero final-pass warnings.
Temporary TeX files were removed. A draft paragraph overfull box was
corrected before that final build. The completed audit reports 198.86 GiB
free disk and 8.46 GiB available memory.
