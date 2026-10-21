# Nonpath tree cuts and cube wedges in dimension nine

The bounded construction family contains **33 distinct simple rational 9-polytope incidences with 18 facets**, obtained from 240 recipes. Four incidence classes are excluded by exact necessary Coxeter constraints. Six have necessary **partial-label** survivors; the other 23 hit explicit search limits. **No compact hyperbolic Coxeter polytope or feasible hyperbolic Gram matrix is claimed.**

The complete recipe and outcome index is [the summary certificate](fresh_h9_tree_cube_wedges_20261014_summary.json). The [geometry certificate](fresh_h9_tree_cube_wedges_20261014.json) contains every rational normal, offset, vertex, cut intersection and explicit facet isomorphism.

## Exact constructive family

Start with the five-cube `-1 <= x_i <= 1`. Use either the star with edges `01,02,03,04`, or the T tree with edges `01,02,03,34`. These are the two nonpath trees on five vertices. For every permutation of its four edges, cut the ridge on the corresponding two original positive facets. If `r=x_i+x_j`, choose the cut level `2-epsilon`, where `epsilon` is half the smallest positive value of `2-r` at a current vertex outside the ridge. This removes precisely the vertices on that ridge, and its crossing edges determine all new vertices by rational interpolation.

Choose one coordinate whose negative facet is left unwedged. Wedge each of the other four negative facets once. The negative facets still meet at the untouched all-negative cube vertex, so the four simultaneous wedge formulas are valid. A wedged negative inequality becomes `-x_i +/- z_j <= 1`. The old cube bounds imply `|x_i| <= 1` and `|z_j| <= 2`; the origin remains strict. Thus these are bounded, full-dimensional rational polytopes, with actual facet and vertex ranks checked exactly.

There are `2 * 24 * 5 = 240` recipes. Canonicalization of the minimal-nonface incidence graphs proposes 33 representatives: five star and 28 T incidences. Every identification is checked as an explicit bijection of both the complete vertex lists and the complete minimal-nonface clutters. None is isomorphic to the frozen [path incidence](fresh_h9_coupled_cube_wedges_20261014.md). Each new minimal-nonface hypergraph is connected, so none is a nontrivial product.

The representative geometries have 378–488 vertices, with 14,988 vertices in total. Missing faces have sizes two or three. The order of the small cuts matters; the certificate includes all orders, not just tree isomorphism types.

## Necessary all-label filter

The [bridge solver](../scripts/fresh_h9_tree_cube_wedges_20261014_bridge_csp.py) adapts the frozen path solver without changing its predecessor. Wedge-twin pairs and the bridge clauses are derived separately for each geometry.

Every missing triple must be Lanner. Its labels are represented by `2,3,4,5,6,>=7`; the last symbol is a quantified category, not a substituted Gram entry. Reciprocal-angle inequalities and spherical vertex classification are constant on this category. Every missing pair keeps an arbitrary strictly dotted distance.

All unspecified meeting contacts are initially relaxed to right angles. At an actual vertex, if even the smallest nonright label 3 destroys positive definiteness, that contact must remain right in every completion. This follows from monotonicity of the largest eigenvalue of the nonnegative weighted adjacency matrix. If two disjoint induced missing-face blocks have every cross contact forced right, their two negative directions are orthogonal, which is incompatible with Lorentz signature.

The search assigns the mandatory triangles and applies these bridge tests. Passing this filter does **not** assign every meeting contact or solve the global rank/signature equations.

| Outcome | Incidence classes | Original recipes |
|---|---:|---:|
| Exact all-label exclusion | 4 | 96 |
| First necessary partial-label survivor | 6 | 52 |
| Explicit node limit; unresolved | 23 | 92 |

The excluded classes are star cases 1–4: the unwedged negative facet is a leaf. They exhaust at 16,155, 19,716, 19,738 and 23,610 nodes. This excludes those complete incidence classes for all omitted finite labels and all strict dotted distances.

The six partial survivors are:

| Case | Tree | Ordered cuts | Unwedged coordinate | Vertices | Nodes to first survivor |
|---:|---|---|---:|---:|---:|
| 0 | star | 01,02,03,04 | 0 | 378 | 9 |
| 8 | T | 01,02,03,34 | 3 | 462 | 839 |
| 17 | T | 01,03,02,34 | 3 | 438 | 895 |
| 22 | T | 01,34,03,02 | 3 | 462 | 895 |
| 25 | T | 03,01,02,34 | 1 | 432 | 1,023 |
| 26 | T | 03,01,02,34 | 2 | 432 | 967 |

Each state is saved in `fresh_h9_tree_cube_wedges_20261014_case_XX.json`. The remaining 23 cases have a 5,000-node limit and retain their first unexplored path. There are 198,847 nodes across the saved searches. No exhausted claim is made for a limited case, nor for an incidence merely because one partial assignment later fails a metric test.

Case 0 supplies the first concrete metric target. Its mandatory graph consists of four diamonds: two Lanner triangles share a label-4 twin edge, with their other four bonds labeled 3. The optional meeting contacts and all 15 missing-pair distances remain free. The separate `fresh_h9_star_cube_metric_20261014` lane tests this **first assignment only**; it does not replace the whole-incidence scope recorded here.

## Verification and scope

- The geometry producer and its saved replay pass. The separate [FLINT primal replay](../scripts/fresh_h9_tree_cube_wedges_20261014_primal.py), with [certificate](fresh_h9_tree_cube_wedges_20261014_primal.json), imports no producer arithmetic. It checks all vertex signs and exact active sets, rank nine, all 18 facet ranks, the explicit bounded box, and complete minimal nonfaces. Every edge at every listed simple vertex has its second endpoint listed; connectivity of the graph of a bounded polytope proves the vertex list complete.
- The separate [incidence review](fresh_h9_tree_cube_wedges_20261014_incidence_review.md) reconstructs all 240 stellar/wedge recipes without producer imports and checks all 105,184 vertex incidences under the saved isomorphisms. It does not independently recanonicalize the 33 representatives.
- The [label replay](../scripts/fresh_h9_tree_cube_wedges_20261014_label_replay.py), with [certificate](fresh_h9_tree_cube_wedges_20261014_label_replay.json), reuses DFS traversal but replaces the Dynkin classifier with independent Fraction Schur elimination over `Q(sqrt(5))`, cycle witnesses, and an exact large-label three-node obstruction. It reproduces the counts and branch hashes for all ten completed searches. Limited cases are not promoted to exclusions.
- The [summary verifier](../scripts/fresh_h9_tree_cube_wedges_20261014_summary.py) binds all 33 outcomes and both arithmetic replays.

This is a new bounded incidence construction and necessary-label search, distinct from the earlier path, product-of-simplices and 14/15-facet surgery families. The unresolved work is metric completion of surviving assignments or further necessary search in the explicitly limited cases, not a claimed finite-volume example.
