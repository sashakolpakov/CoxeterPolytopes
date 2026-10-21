# Actual H3 count across all existing ordinary cross payments

The [source](../scripts/fresh_sqrt5_h3_all_actual_count_split_20261016.py) and
[exact certificate](fresh_sqrt5_h3_all_actual_count_split_20261016.json) separate
actual H3 count in all six existing A3/A4/D4/D5 source-root orbits. This adds
392 coefficients to 267 existing count-one families, giving 32,446 variables.
Each root pair still pays β/(3c). The separate central-root β/c family is
unchanged. The two earlier A3/A4 multiplicity barriers are explicitly broken.

All 18,035 inherited marked forests were recounted literally and checked
against their DP component sums. The face constraints were rebuilt from all
996 triangle states, 37 paid-quadrilateral stopping cases, and 1,979 inherited
baseline-zero quadrilaterals: 102,182 distinct profiles in total. There are
1,912,875 triangle, 133,484 paid-quad, and 227,408 inherited-quad common contexts.
The source also checks 559 literal pair-count/marked-root normalization cases.

The fixed-budget-35 finite-cut problem remains infeasible. Exact weak duality
proves

\[
T\ge\frac{205603142035642}{5413000027065}\approx37.9832146698.
\]

This is a certified lower bound, not an exact optimum. The numerical finite-cut
minimum is approximately 37.9878083808, with adjacent discount approximately
0.2248069526. The exact certificate has 326 vertex, 197 triangle, 1,019 quad,
and one distance-two-zero terms. No full maximizing-forest oracle was needed
because these finite literal cuts already exclude target 35. Generation and
fresh-process saved replay both pass.

The dominant new coupling is clarified by the
[four-triangle H4 barrier](fresh_sqrt5_h4_triangle_compensation_barrier_20261016.md):
H4 compensation needed in one complementary-component context is charged to
three H4 copies beside E8. Further H3 count refinement cannot remove that
particular barrier.

This excludes a payment family, not a geometric polytope. Complete newly
affected paid-quadrilateral extension coverage remains open. The frozen
compact arithmetic dimension-21 candidate is unaffected by this unsuccessful
attempt to exclude dimension 21.
