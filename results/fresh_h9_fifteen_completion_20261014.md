# Completion of the saved fifteen-facet truncation frontier

**All 45,679 distinct nerves in the saved one-truncation family are excluded.** The 1,334 formerly open cases now have complete exact label obstructions, and all 1,334 pass a separate arithmetic replay. No compact H9 example was found. The [final ledger](fresh_h9_fifteen_completion_20261014_certificate.json) has no unsearched, node-limited or surviving case. Earlier sources and certificates remain unchanged.

This run revalidates the exact 1,334 open case IDs in the [previous ledger](fresh_h9_fifteen_facet_20261014_checkpoint.json). It completes their necessary-label searches with a limit of 20,000 states per case. Cases with more dotted pairs and missing triangles are tested first, to prioritize potential metric flexibility. The driver stops dispatching at the first full-category survivor; only such a survivor is eligible for the prepared Lorentzian metric solver. Numerical failure is never an exclusion criterion.

| Completed stage | Count |
|---|---:|
| Previously open nerves searched | 1,334 |
| Exact complete exclusions | 1,334 |
| Node-limited or unsearched cases | 0 |
| Complete necessary label survivors | 0 |
| Search states in this pass | 15,263 |
| Largest individual search | 109 states |
| Separate arithmetic replays passed | 1,334 |
| Distinct one-cut nerves excluded, including earlier gates | 45,679 |

The complete global searches across both checkpoints number 1,508 and use 18,110 states. No metric fit was launched, because no complete necessary labeling survived. In this new pass the exact pruning records contain 9,076 failures to extend a Lannér table, 3,949 forced orthogonal-Lannér obstructions, 519 empty relation joins and 12 nonspherical-vertex failures. These are search-event counts, not disjoint nerve counts. Every branch closes before any modular-minor test is reached: the summed number of such tested minor codes is zero.

The [primary records](fresh_h9_fifteen_completion_20261014.json) bind each immutable case file, the unchanged fifteen-node engine, the incidence source and the relevant relation tables. Every case is resumable separately. The [separate arithmetic replay](fresh_h9_fifteen_completion_20261014_independent.json) reconstructs its maximal faces independently, substitutes a separate finite-tree classifier and SymPy arithmetic for the finite-field determinants, and compares the entire result. It reuses the DFS traversal, propagation and classical Lannér tables. This is a separate arithmetic replay, not a second independently written complete search. All cases pass, as do generation and saved-result replay of the [final certificate checker](../scripts/fresh_h9_fifteen_completion_20261014_certificate.py). The two complete jobs have terminated; no background search remains.

The label domain is all finite integer Coxeter labels. Labels 2 through 6 are literal. Category 7 represents every finite label at least 7: such an edge can belong to a spherical diagram only as an isolated rank-two component, and its contribution to a Lannér triangle has the same category throughout that range. Indeed, the other two labels `(2,2)` always give a spherical triangle; every other pair has reciprocal sum at most `1/2+1/3+1/7<1`. The complete Lannér tables of orders four and five have no such large labels. Fixed eleven-minor tests skip entries in category 7. No maximum-label or missing-triangle restriction particular to fourteen facets is used.

The inherited incidence census contains 2,124,486 literal codimension-two through five cuts of the 770 independently realized fourteen-facet Gale parents. Of these, 2,049,175 leave an old induced obstruction unchanged. The remaining 75,311 cuts reduce, by explicitly verified symmetries and isomorphisms, to 45,679 distinct fifteen-facet nerves. The [previous report](fresh_h9_fifteen_facet_20261014.md) already excludes 44,345 of those nerves and documents the separate set-valued incidence replay. The present searches address precisely its remaining 1,334 IDs, without sampling or altering the earlier family.

A successful completion of this ledger also covers every nontrivial standard single-face truncation of these parents. For a cut of codimension at least six, stellar subdivision removes the cut support while retaining every proper subset. That support is therefore a minimal nonface of order at least six. In a compact Coxeter realization such a minimal nonface must be Lannér, whose order is at most five. This extension concerns a single sufficiently small face truncation. It does not apply to several subsequent cuts that destroy the resulting missing face, or to an arbitrary hyperplane cut changing a larger part of the incidence structure.

The next genuinely different constructive family is **direct uniform rank-five Gale configurations on fifteen points**, rather than another cut of these parents. A concrete exact realization test uses fifteen-by-five integer matrices with every five-row minor nonzero. Positive six-element circuits supply the complementary nine-vertex facets; a positive circuit sum supplies a strictly positive kernel vector. Row rescaling and a rational basis of the transpose kernel then construct a nine-dimensional primal realization. Facet coverage checks that every label is an actual vertex. Minimal nonfaces of order at most five and the existing literal Lannér/finite-vertex/rank gates can then be applied.

To force this next family outside all standard one-face truncations, one can require that no vertex link has a join factor equal to a simplex boundary. A new stellar vertex created by truncating a support `F` necessarily has link `boundary(F) * link(F)`. In the minimal-nonface hypergraph of that link, the boundary factor appears as a component consisting of the single missing face `F`. Excluding every such factor at every vertex is a sufficient family-selection filter; it is not itself a Coxeter obstruction. No direct fifteen-by-five Gale search or claim that a metric survivor exists is made here. The completed one-cut family does not classify all fifteen-facet compact-H9 incidence types.

The inexpensive final ledger replay and an optional fresh arithmetic replay are:

```sh
PYTHONPATH=src:. python scripts/fresh_h9_fifteen_completion_20261014_certificate.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_fifteen_completion_20261014_replay.py --jobs 2
```

The second command reruns every completed case without using the resumable audit cache. The archived streaming replay already performs this arithmetic once per case; it is not necessary to repeat it merely to check the final hashes and counts.
