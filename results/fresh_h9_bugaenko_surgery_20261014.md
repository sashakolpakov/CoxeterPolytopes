# Bugaenko H8 wedge and two face truncations

No new Coxeter polytope was obtained. The completed incidence census gives **2,953 distinct admissible 14-facet simple H9 types** in this surgery family. Exact label joins exclude eleven of them. Further bounded searches produced several necessary labelings, but the tested metric proposals fail exact rank or signature conditions. The other 2,942 incidence types remain unresolved by this checkpoint.

This family was not present in the earlier polygon, 120-cell, higher-base, or targeted stellar-subdivision archives. It starts with the actual 11-facet compact Bugaenko H8 polytope, wedges one facet, and then truncates two faces, each of codimension two through five. The second face may use the first new facet. The two truncations are ordered; final isomorphic nerves are merged.

## Complete polytopal incidence census

The [known source](../data/known_examples/bugaenko_compact_h8.json) has 41 vertices and minimal nonfaces

```text
{9,10}, {0,1,2,3,4}, {0,1,2,3,9},
{4,5,6,7,8}, {5,6,7,8,10}.
```

Its combinatorial symmetry blocks have sizes `4,1,4,1,1`, arranged as a pentagon multiwedge. There are three wedge-carrier orbits. The first truncation is enumerated by counts in these blocks. For the second, each block is split by membership in the first truncation face; the new facet is handled separately. Stabilizers of the first face give safe reductions. Additional symmetries after a cut can cause repetition, but cannot remove a missing recipe. Colored bipartite incidence graphs are finally canonicalized with nauty.

| Wedge block | First-face orbits | Ordered two-cut recipe orbits | Passing missing-face order |
|---|---:|---:|---:|
| First block of size four | 77 | 19,876 | 1,285 |
| Central singleton | 52 | 15,790 | 208 |
| Outer singleton | 107 | 37,847 | 4,493 |
| Total | 236 | 73,513 | 5,986 |

The other 67,527 recipes retain a minimal nonface of order greater than five and cannot be compact hyperbolic Coxeter polytopes. The 5,986 passing operation orbits yield 2,953 canonical nerves. These have 432,736 vertex occurrences in total, all with nine incident facets.

Polytopality is constructive, not inferred from a sphere test. In an affine realization `a_i x≤b_i`, wedge the chosen inequality into `a_f x+t≤b_f` and `a_f x−t≤b_f`. Truncate a face F by

```text
sum(i in F) a_i x ≤ sum(i in F) b_i − epsilon.
```

Take epsilon positive and smaller than the least positive sum of F-slacks at a vertex outside F. Repeating this construction realizes every saved recipe. These are Euclidean/Klein convex realizations of the incidence; they do not provide Coxeter metrics.

The [producer and complete catalog](../scripts/fresh_h9_bugaenko_surgery_20261014_incidence.py), [JSON](fresh_h9_bugaenko_surgery_20261014_incidence.json), and [separate incidence checker](../scripts/fresh_h9_bugaenko_surgery_20261014_incidence_independent.py) are archived. The independent checker reconstructs every recipe from vertex incidences and recovers every minimal nonface from the **entire downward face complex**, without importing the producer's clutter formula. It checks all 432,736 vertices, every ridge degree, and every canonical relabeling. Its [saved result](fresh_h9_bugaenko_surgery_20261014_incidence_independent.json) passes. The independent implementation does not duplicate nauty or the operation-orbit enumeration; those have a fresh producer replay.

## Exact label gates and bounded metric work

The [first 25 complete-table joins](fresh_h9_bugaenko_surgery_20261014_labels_0_25.json) exclude cases

```text
1, 4, 5, 6, 8, 9, 13, 16, 17, 18, 24.
```

Each exclusion reaches an empty join before the 5,000-state limit. Thirteen other joins hit that explicit limit. Case 0 instead completes with 864 necessary partial labelings; only one of these is fully rejected by the subsequent independent meeting-label/rank replay. All fourteen nonempty cases remain unresolved. The tables impose Lannér labels on every minimal nonface and spherical diagrams at every actual vertex. Unknown meeting pairs are relaxed to right angles, which is a necessary relaxation by monotonicity of the nonnegative cosine matrix. Finite labels at least seven share one category in these tests: they can occur only as isolated rank-two spherical components, while the rank-three Lannér condition is constant over that category. Parallel pairs are not conflated with this category.

A [bounded DFS](fresh_h9_bugaenko_surgery_20261014_dfs_0_5.json) found one necessary partial labeling on each of cases 0, 2, and 3. Completing **all** still-unassigned meeting-pair categories, not simply holding them right, gives respectively 74, 8, and 135 spherical-compatible extensions. Every one of these 217 extensions has a fixed 11-by-11 Gram minor with nonzero determinant modulo 241. The minor avoids every dotted distance and every variable label at least seven. Thus its rank is at least eleven, impossible in the ten-dimensional Lorentz space of H9.

The modular argument uses the ring generated by `sqrt(2), sqrt(3), sqrt(5)` with powers of two inverted. All three radicands have square roots modulo 241. A nonzero image of a Gram determinant proves that determinant is nonzero in characteristic zero; it is not a numerical rank test.

The [separate finite-category and rank implementation](../scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py) re-enumerates these 217 extensions using its own finite-diagram classifier and computes determinants with SymPy's finite-field matrix arithmetic. It also replays all 2,358 archived minor-witness records from a [larger bounded search](fresh_h9_bugaenko_surgery_20261014_rank_search.json), accounting for 30,162 rejected full category extensions. That search uses a 3,000-node limit and a 20,000 raw free-extension limit. Its eight limited cases remain unresolved.

The larger search's one category surviving its modular-minor filter, on case 318, has two orthogonal Lannér triangles `{5,6,7}` and `{4,8,9}`, each with labels `(2,3,m)` for finite `m≥7`. They give at least two negative directions, so this particular category is also impossible. This excludes the saved labeling, not the whole case-318 incidence.

The [independent rank certificate](fresh_h9_bugaenko_surgery_20261014_rank_independent.json) contains these exact conclusions. The [numerical metric scout](fresh_h9_bugaenko_surgery_20261014_metric_scout.json) is preserved separately: it tested the first 217 extensions by a rank-one negative Schur complement, three initial values, dotted parameters at most 30, and a continuous relaxation for labels at least seven. Best residuals were about 0.0675, 0.3041, and 0.4806. Those numerical failures alone prove nothing; the exact minors supply the stated exclusions.

No full metric, compactness certificate, or finite-volume candidate survived. This checkpoint supplies a complete new incidence catalog and reusable exact rank gates, not an exclusion of the entire surgery family.

## Replay

```sh
python scripts/fresh_h9_bugaenko_surgery_20261014_incidence.py --verify
python scripts/fresh_h9_bugaenko_surgery_20261014_incidence_independent.py --verify
python scripts/fresh_h9_bugaenko_surgery_20261014_labels.py --end 25 --max-states 5000 --verify
python scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py --verify
```

The independent rank replay validates saved exact certificates and the three full meeting-label completions. It does not independently repeat the larger DFS traversal or the numerical optimization. No search process remains active after the final replays.
