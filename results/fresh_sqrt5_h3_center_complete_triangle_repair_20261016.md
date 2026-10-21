# The repaired uniform-discount family is obstructed

Adding the complete 996-state affected-triangle catalogue and 37 paid-quad
stopping rows makes the current payment family infeasible at target 35. The
[source](../scripts/fresh_sqrt5_h3_center_complete_triangle_repair_20261016.py)
and [exact dual](fresh_sqrt5_h3_center_complete_triangle_repair_20261016.json)
give the rigorous restricted-family lower bound

\[
T\ge \frac{1128458321210371}{28991000115964}
\approx 38.92443574546.
\]

This is a lower bound, not a certified optimum. The numerical finite-cut
relaxation is approximately 38.92561739167. It has both discounts zero. There
is no new dimension bound and no infeasibility claim about geometric polytopes.

The new model replaces the old 196 selected triangles, including their one
conjugate-inadmissible row, with all 996 states of the
[complete triangle producer](fresh_sqrt5_triangle_extension_complete_20261016_rank22.json).
It imposes all 1,912,875 triangle context/mark cases, all 133,484 context/mark
cases for the 37 new paid quads, and the inherited 1,979 baseline-zero quad rows
with 227,408 common contexts. After deduplication there are 102,182 profiles.
Explicit inequalities fix the distance-two discount to zero and bound the
adjacent discount by 2/3; these restrictions participate in exact weak duality.
All 18,035 inherited vertex cuts are literal valid marked rank-21 forests.

The dual has 298 vertex terms, 1,027 quad terms, one distance-two-zero bound,
and **no triangle terms**. The dominant new row is paid quadrilateral 4266,
whose profile is exactly 1 minus the adjacent discount, with no available
bonus coefficients. It forces that uniform discount to zero. The
[small exact barrier](fresh_sqrt5_uniform_adjacent_discount_barrier_20261016.md)
isolates the simpler consequence T≥38 and identifies the E6 root-pair orbit
responsible. Thus the failure is not evidence that adding the triangle
catalogue alone creates the obstruction.

The largest vertex dual masses use marked-A1 forests 2A1+H4+E8+E7 and
A1+D5+E8+E7, followed by forests with two E7 or two D7 components. Every vertex
term is recounted by literal root pairs in the exact replay. Every face term
is reconstructed from its common external counts and facet mark; no corner
gets an independently chosen spectator context.

Generation and fresh saved replay pass. This remains a primary restricted
payment-model certificate; complete affected paid-quadrilateral extension
coverage is open. That missing coverage cannot repair the current lower bound,
but it is mandatory before any feasible successor rule can imply a theorem.
