# Fifteen-facet construction by repairing the saved Gale nerves

No compact H9 example was found. This round constructs **45,679 distinct realizable fifteen-facet nerves** by one additional face truncation of the 770 saved fourteen-facet Gale parents. Exact local and global necessary conditions exclude 44,345 of them; **1,334 remain open without a completed global label search**. The [frozen ledger](fresh_h9_fifteen_facet_20261014_checkpoint.json) lists every open case and binds the certificates below. This family is not the class of all fifteen-facet polytopes.

The constructive input is the independently realized [770-parent Gale catalogue](fresh_h9_gale_frontier_20261014.md). A sufficiently small rational truncation of any codimension-two through five face gives a simple H9 combinatorial candidate with fifteen facets. Its dual operation is stellar subdivision at the corresponding two- through five-node face. This supplies an actual polytopal realization, without assuming that the original parent already admits a Coxeter metric.

If a cut face F is not contained in an old obstruction support U, the induced nerve on U is unchanged. A potentially repairing cut must therefore lie in the intersection of all retained obstruction supports. The [incidence census](fresh_h9_fifteen_facet_20261014_incidence.json) checks this criterion before constructing candidates; every new exclusion is then tested again on the resulting fifteen-node nerve.

| Stage | Exact count |
|---|---:|
| Literal codimension-two through five faces of 770 parents |2,124,486|
| Cuts leaving a certified old induced obstruction unchanged |2,049,175|
| Potentially repairing literal cuts |75,311|
| Operation orbits under explicitly checked parent permutations |45,888|
| Distinct resulting fifteen-facet nerves |45,679|
| Recomputed inherited local/two-Lannér gates exclude |39,931|
| Two new exact induced-pattern gates exclude, additionally |4,240|
| Completed global label searches exclude, additionally |174|
| Remaining cases |1,334|

For 490 parents, every permitted cut already fails the unchanged-obstruction test. The remaining 280 contribute the eligible cuts. All minimal nonfaces in the constructed nerves have order at most five. Symmetry reduction uses literal generators checked on both old minimal nonfaces and maximal faces; it does not require an assumed full automorphism group.

The [new local-pattern trials](fresh_h9_fifteen_facet_20261014_patterns.json) tested 30 selected patterns with a 10,000-state limit. Two failures are complete: an eight-node pattern with missing-face masks `[213,234,244,248]`, and a nine-node pattern with masks `[15,55,115,496]`. Their 3,641 and 687 applications overlap in 88 nerves. The 19 local label survivors and nine node-limited trials give no exclusion.

The global searches first tested 50 relatively simple nerves, then 100 cases chosen across different parents, and finally [24 metric-priority cases](fresh_h9_fifteen_facet_20261014_metric_priority.json) with more dotted pairs and distinct missing-face count signatures. All 174 searches terminate completely, using 2,847 states in total. None reaches a complete necessary label assignment, so no numerical metric fit was used or reported as a negative certificate. Their respective limits were 3,000, 2,000 and 5,000 states per case; no global case hit its limit. The [first 150 results](fresh_h9_fifteen_facet_20261014_search.json) preserve detailed propagation counters and sample witnesses.

The fifteen-facet search allows missing triangles and every finite integer Coxeter label. Labels 2 through 6 are literal; category 7 represents all finite m≥7. At a spherical vertex such a large label can occur only in an isolated I2(m) component. In a Lannér triangle, replacing m≥7 by 7 preserves the exact spherical/hyperbolic category: the other two labels `(2,2)` are always spherical, while every other pair gives reciprocal sum at most `1/2+1/3+1/7<1`. Larger Lannér tables have no such labels. Fixed eleven-minor tests skip category 7. Thus no fourteen-facet-specific missing-triangle or maximum-label restriction is transferred here. The successful exclusions use spherical/Lannér propagation and incompatible disjoint Lannér blocks.

The [separate incidence replay](fresh_h9_fifteen_facet_20261014_incidence_independent.json) rebuilds every candidate from set-valued maximal faces, recovers its entire minimal-nonface clutter from the downward face complex, and checks 9,817,894 vertex incidences. It checks every source cut, every operation orbit, all 209 duplicate-recipe isomorphisms and every literal local-gate application. Both new local failures are rerun with a separate join traversal and finite-diagram classifier. The old rational Gale realizations, 44 previously proved local obstructions, two-block theorem and classical Lannér tables remain explicit inherited dependencies.

The [150-case arithmetic replay](fresh_h9_fifteen_facet_20261014_search_independent.json) and [24-case replay](fresh_h9_fifteen_facet_20261014_metric_priority_independent.json) use a separately implemented finite-tree classifier and SymPy finite-field determinants. They reuse the DFS propagation, canonical-state traversal and classical Lannér tables; this is separate arithmetic, not a second independently written full search. The [peer source/proof review](fresh_h9_fifteen_facet_20261014_review.md) checks the fifteen-node port and operation-coverage argument without claiming a duplicate census.

Replay commands from the repository root:

```sh
PYTHONPATH=src:. python scripts/fresh_h9_fifteen_facet_20261014_incidence_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_fifteen_facet_20261014_search_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_fifteen_facet_20261014_metric_priority_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_fifteen_facet_20261014_checkpoint.py --verify
```

The next constructive frontier is precisely the 1,334 case IDs in the ledger. Their metric feasibility is unresolved. There is no claim that every truncation of every fourteen-facet polytope, every fifteen-facet H9 incidence type, or every compact H9 construction is excluded.
