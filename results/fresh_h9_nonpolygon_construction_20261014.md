# Compact H9 construction from nonpolygon bases

No compact H9 polytope was found. Three new **globally polytopal, nonproduct** incidence constructions were tested before metric fitting. Every tested type fails an exact necessary Coxeter-label join. These results concern the displayed families; they do not exclude general compact H9 polytopes or all multiwedges.

| Constructive family | H9 facets | Exact tested scope | Outcome |
|---|---:|---|---|
| Six wedges over a simple dodecahedron | 18 | All 12,376 positive weight compositions; 82 symmetry types survive the missing-face bound | All 82 label joins empty |
| Six edge truncations of a cube, then wedges on its six original facets | 18 | One lexicographic truncation order for each of 30 six-original-edge subset orbits | All 30 label joins empty |
| A vertex-truncated dodecahedral prism, then five specified wedges | 20 | Three target-vertex symmetry classes | All three label joins empty |

No state cutoff was reached. The largest intermediate label set contains 1,332 states. No numerical Gram fit, optimization result, or arithmetic-lattice assumption is used in these exclusions.

## Why these are actual and different polytopes

For an ordinary wedge, replace a facet inequality `a·x≤b` by `a·x+z≤b` and `a·x−z≤b`. The new coordinate satisfies `|z|≤b−a·x`; boundedness and simplicity are preserved, while dimension and facet count each increase by one. Repetition gives a multiwedge.

If a simple base has facet blocks `B_i` of positive sizes `w_i`, each base minimal nonface `M` becomes the union of its blocks. A base vertex with incident facet set `F` gives all vertices obtained by keeping every element of blocks in `F` and omitting exactly one element from every other block. All resulting vertex sets in this report contain nine facets.

These constructions are not the previously excluded `Δ2³×I³` models or polygon multiwedges. The [separate incidence audit](fresh_h9_nonpolygon_construction_20261014_incidence_audit.json) recovers each original block from identical minimal-nonface membership. The recovered base dimensions are three, three, and four, respectively. Every minimal-nonface hypergraph is connected, excluding any nontrivial product decomposition. The thirty edge-truncation recipes need not be pairwise inequivalent; their stated scope is a recipe list, not an unproved isomorphism census.

## Dodecahedral bases: all 82 weight types

The [first producer](../scripts/fresh_h9_nonpolygon_construction_20261014.py) starts with the twelve exact icosahedral vertices

```text
(0, ±1, ±φ), (±φ, 0, ±1), (±1, ±φ, 0),  φ²=φ+1.
```

Their common squared norm is `2+φ`, and adjacent vertices have dot product `φ`. All twenty triangular supporting planes are checked exactly, with strict inequalities at every other vertex. The vertices are extreme by their norm supports; twenty triangular facets attain the three-dimensional bound `2v−4`, proving boundary completeness. Central symmetry puts the origin in the interior. The polar is therefore an actual simple dodecahedron with twelve facets.

The complete 120-element symmetry group is enumerated and checked for distinctness, closure, and preservation of the triangular boundary. A triangle flag determines an automorphism by propagation across adjacent triangles, giving the matching upper bound `20×6=120`.

Six wedges require `Σw_i=18`. Every minimal nonface of a compact Coxeter polytope carries a Lannér diagram, of order at most five. For this flag base the minimal nonfaces are its 36 disjoint facet pairs, so their two weights must sum to at most five. Of all `C(17,11)=12,376` positive compositions, 6,644 pass. They form exactly 82 orbits under the complete base symmetry group. All are processed in the [certificate](fresh_h9_nonpolygon_construction_20261014.json). Their vertex counts range from 218 to 600.

All 82 types fail after a small number of mandatory Lannér placements. This exclusion covers every six-wedge weight assignment on the dodecahedral base, not just equal weights or a selected metric ansatz.

## Edge truncations as a combinatorial repair

The [second construction](../scripts/fresh_h9_nonpolygon_construction_20261014_edgecut.py) starts from an octahedron, dual to the cube. It chooses six distinct original edges; there are 924 subsets and 30 orbits under the 48-element octahedral symmetry group. One lexicographic order is used for each representative. Other orders can produce different types and are not classified here.

Each edge stellar subdivision is realized by a rational point on the radial line through its midpoint. Its scale is chosen strictly beyond the two incident facets and strictly below every other facet. All new triangular supporting facets are checked exactly after every step. Dualizing gives six actual edge truncations of the cube. The resulting simple base has twelve facets; wedging once on each of the six original facets produces an eighteen-facet H9 polytope.

This was a constructive attempt to remove finite-vertex obstructions rather than alter metric parameters on an old incidence type. Edge truncations preserve the relevant missing-face order bound, whereas truncating a vertex incident to three doubled blocks would create a forbidden missing six-set. Nevertheless, all thirty selected repairs fail after three mandatory Lannér-four placements. Their [archive](fresh_h9_nonpolygon_construction_20261014_edgecut.json) includes the exact rational convex realizations and all resulting vertex sets.

## A four-dimensional base with fewer large missing faces

The [third construction](../scripts/fresh_h9_nonpolygon_construction_20261014_fourbase.py) deliberately avoids the repeated Lannér-four pattern above. Its simple base is a dodecahedron times an interval. Choose three mutually meeting dodecahedral facets and both interval facets for the five wedges. Before truncation this would retain a product structure.

To break that product, truncate a prism vertex whose three dodecahedral facets are disjoint from the chosen triple. In the dual, the prism is the convex hull of the icosahedron in `R³×{0}` and the two interval poles. The new vertex is `101/100` times the barycenter of the selected tetrahedral facet. Exact arithmetic proves that it lies beyond precisely that facet and below all thirty-nine others. This is an actual stacking operation, so the dual base is a simple four-polytope with fifteen facets. The five wedges give twenty H9 facets.

There are ten eligible opposite dodecahedral vertices; the weighted-triangle stabilizer gives three classes, with the interval reversal identifying the two ends. Every resulting polytope has 508 vertices. Its minimal nonfaces include only **one order-four and one order-five set**, alongside smaller missing faces. This reduces the initial metric constraints substantially compared with the earlier branches.

Even these two mandatory diagrams are incompatible with the complete actual-vertex restrictions. All three joins become empty at the second placement. The [certificate](fresh_h9_nonpolygon_construction_20261014_fourbase.json) records the old supporting normals, the new stacking point and exact inequalities, the resulting boundary, and the full H9 incidences.

## Exact metric gate and independent incidence check

Each missing four- or five-set is assigned every one of the complete classical Lannér placements: 108 labelled order-four diagrams and 420 labelled order-five diagrams. The previously audited [classical table implementation](../scripts/independent_h9_14_record30_a1_labels.py) supplies these finite lists. The join requires agreement on shared pairs and tests every actual nine-facet vertex restriction.

Unassigned meeting pairs are relaxed to right angles. This is safe for exclusion: a normalized acute Gram is `I−A` with `A` entrywise nonnegative. Restoring any omitted angle can only increase `A`, so Perron–Frobenius monotonicity prevents a failed positive-definiteness test from being repaired. Missing triples, disjoint-pair distances, and all unassigned labels are left unconstrained. Empty joins therefore exclude every compact Coxeter metric on the tested type, not merely labels bounded by five. Finite tests are exact diagram classifications, and every saved join ends before its state cap.

The [separate stdlib audit](../scripts/fresh_h9_nonpolygon_construction_20261014_incidence_audit.py) imports no producer or diagram-classification code. For all 115 records it enumerates every nine-subset directly from the minimal nonfaces and checks equality with the saved vertex list, verifies minimality of every missing face and degree two of every ridge, reconstructs the wedge blocks and base dimension, and checks the nonproduct criterion. This direct check covers 5,949,320 nine-subsets and passes both generation and saved-result replay. It does not constitute a second implementation of the label joins. All three primary saved-result replays also pass.

## Reproduction and remaining scope

```sh
python scripts/fresh_h9_nonpolygon_construction_20261014.py --limit 1000 --max-states 5000 --verify
python scripts/fresh_h9_nonpolygon_construction_20261014_edgecut.py --verify
python scripts/fresh_h9_nonpolygon_construction_20261014_fourbase.py --verify
python scripts/fresh_h9_nonpolygon_construction_20261014_incidence_audit.py --verify
```

The argument depends on the complete classical Lannér tables and finite-diagram classification already bound by the source manifests. It establishes no new universal dimension bound. Other nonpolygon bases, different truncation orders and weight patterns, and constructions that avoid these mandatory-diagram conflicts remain open. No metric candidate from this checkpoint reaches a signature or compactness test.
