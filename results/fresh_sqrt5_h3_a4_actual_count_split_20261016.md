# Splitting the H3 count for A4 endpoints

The [source](../scripts/fresh_sqrt5_h3_a4_actual_count_split_20261016.py) and
[certificate](fresh_sqrt5_h3_a4_actual_count_split_20261016.json) give independent
coefficients at actual H3 counts 2–5 to the existing H3×A4-end payment orbit.
Count 1 keeps the original coefficient, and each individual pair still pays
β/(3c). This adds 67 coefficients, for 32,121 variables total. It changes no
other bonus family and retains the E6 edge-orbit exemption.

The previous c=5 face/c=1 vertex barrier is explicitly broken: quad 7661 now
uses a different coefficient from the marked 2A1+H3+4A4 vertex. All 18,035
inherited forest cuts were rebuilt by literal pair enumeration and checked
against independently assembled DP component sums. All 996 triangle states,
37 paid-quad stopping rows, and inherited baseline-zero quads were recomputed.

The fixed-35 finite-cut LP is nevertheless infeasible. Its exact dual proves

\[
T\ge1195729577/31999936\approx37.36662401450.
\]

This is a lower bound, not an optimum. The numerical finite-cut minimum is
approximately 38.19079549855, with adjacent discount approximately
0.17771813260. The exact dual has 288 vertex, 144 triangle, 952 quad, and one
distance-two-zero terms. No full maximizing-forest oracle was needed because
these valid finite cuts already obstruct target 35. All dual rows are checked
exactly; generation and fresh saved replay pass.

The next dominant multiplicity defect is the H3×A3-end payment, isolated in
the [A3 barrier](fresh_sqrt5_h3_a3_multiplicity_barrier_20261016.md). A separate
source/accounting read checked the actual-count lookup, unchanged denominator,
facet losses, module hooks, and literal/DP reconstruction. This was not an
independent full LP review.

The result excludes this payment model at 35, not a geometric polytope.
Complete newly affected paid-quadrilateral extension coverage remains open.
