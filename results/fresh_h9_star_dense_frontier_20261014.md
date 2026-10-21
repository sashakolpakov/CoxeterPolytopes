The case 0 exclusion now allows ten previously fixed own contacts to be arbitrary finite labels at least 3. Facet 14 still has the two own labels exactly 3, and all four internal block labels are explicitly at least 4. This is a larger excluded label family on the same star incidence; it does not exclude the whole incidence.

The [binder](../scripts/fresh_h9_star_dense_frontier_20261014.py) and [certificate](fresh_h9_star_dense_frontier_20261014.json) append one predicate to the previous 159, giving 160 certified predicates. The new certificate binds the separate [five-family dependency review](fresh_h9_star_dense_frontier_20261014_dependency_review.md) and its exact data. The five inherited arithmetic proofs already have independent reviews, covering 155,680 residual cases in total; this extension does not repeat that arithmetic.

The common eleven-normal rank core is

\[
\{3,4,6,7,9,10,12,13\}\cup\{11,14,17\}.
\]

The exact core contacts are \((3,14)=(4,14)=3\). The widened contacts are

\[
(2,3),(2,4),(5,6),(5,7),(6,15),(7,15),
(8,9),(8,10),(9,16),(10,16).
\]

Each has an endpoint outside the core. Thus all residual Gram entries, determinant formulas and polynomial coefficients in the five proofs are literally unchanged. This is not a relaxation of the entire Gram matrix, which would not preserve its rank.

Outside the core, the changed contacts enter only lower-bound reductions and tests that an actual vertex Gram is not positive definite. Increasing a finite Coxeter label decreases its normalized off-diagonal Gram entry. For a symmetric matrix with nonpositive off-diagonal entries, replacing a nonpositive quadratic witness by its componentwise absolute value can only decrease its quadratic value. Further decreasing those entries therefore preserves the failure of positive definiteness. This also handles saved LDL witnesses whose original coordinates have mixed signs.

The independent dependency checker inspects all five proofs, including the inherited all-own3 and opposite-single branches. It verifies17 bound-only high-block reductions, six bound-only logical reductions, the earlier 30 actual-vertex guards and all stored literal witness supports. The four internal lower bounds are explicit hypotheses; they are not inferred from a missing triple after its own contacts have changed. The first three internal labels reduce as before, and the complete last-internal-at-least4 partition remains valid for arbitrary last-anchor own profiles.

The [resumable search](../scripts/fresh_h9_star_dense_frontier_20261014_search.py) continues the frozen [186-visit cursor](fresh_h9_star_lastblock_complete_20261014_search_batch_006.json). It validates every stored path, assignment, domain and pending sibling cursor before applying the new predicate. Category 7 continues to mean every finite integer label at least 7; dotted distances remain unrestricted. A full mandatory-label survivor still needs omitted meeting-contact and Lorentzian metric completion.

The [saved batch007](fresh_h9_star_dense_frontier_20261014_search_batch_007.json) stops at the next mandatory-label survivor after 229 new visits, **415 cumulative visits**, with nine pending frames. The widened predicate pruned 84 states. Other rejection counts are 304 failed spherical-vertex child trials, 51 orthogonal negative-block states, and 78 states covered by older analytic predicates. Child-trial counts are distinct from newly entered DFS nodes.

The next state has internal block labels \((4,4,4,5)\). Every own profile is \((3,3)\), except facet 16, whose profile is \((2,5)\): \((9,16)=2\) and \((10,16)=5\). This is genuinely outside the new theorem, since one of the ten widened contacts is right. It is preserved for later metric/domain analysis; no metric existence or whole-incidence conclusion is asserted. The earlier conditional last-internal-3 cross-right guards remain a separate available tool, but that is not this saved state's internal label.

Validation: binder generation and saved replay passed. The batch uses the unchanged previously reviewed resumable engine, with only the new source tag, gate path and required gate status changed. Its initial state is the saved 186-visit cursor, not an empty restart. The new dependency review is recursively bound through the gate. The final saved batch replay also passed, reproducing all counts, rejection data, the survivor and every pending cursor. No jobs remain in this prefix.
