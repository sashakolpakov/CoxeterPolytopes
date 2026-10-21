# Higher-dimensional bases for compact H9 constructions

No new Coxeter polytope was found. This round constructs and excludes three fresh, explicitly polytopal families: three cyclic-base multiwedge families, fourteen clustered five-wedge recipes on the simple 120-cell, and a different pentagon five-wedge recipe on that cell. These are necessary-angle obstructions to the specified compact H9 realizations, not an exclusion of all higher-base constructions or of noncompact polytopes.

The work changes the base dimension and incidence structure from the previously excluded polygon and product families. A search of the repository found no prior search of these exact recipes. This is a statement about repository coverage, not literature novelty.

| Construction | Complete scope checked | Result |
| --- | --- | --- |
| Polar of C(10,4), raised to dimension nine | All 2,002 positive weights summing to 15 | Every nerve has a missing face of order above five |
| Polar of C(11,5), raised to dimension nine | All 1,001 positive weights summing to 15 | Same obstruction |
| Polar of C(12,6), raised to dimension nine | All 364 positive weights summing to 15 | Same obstruction |
| Simple 120-cell, five distinct wedges | Four chosen dual vertices form a tetrahedron; all fourteen orbits of the fifth vertex | Every necessary compact Coxeter labeling is excluded |
| Simple 120-cell, five distinct wedges | The five vertices in the link of the old edge (0,33) | Every necessary compact Coxeter labeling is excluded |

## Convex realization and complete incidence

The cyclic bases use rational moment-curve coordinates. Every supporting facet and every minimal nonface is enumerated exactly. A separate checker obtains the same facets from Gale evenness and enumerates all weight compositions independently. A compact hyperbolic Coxeter polytope cannot have a minimal nonface of order above five, so these three families stop before metric fitting. See the [cyclic certificate](fresh_h9_highbase_construction_20261014_cyclic.json) and [producer](../scripts/fresh_h9_highbase_construction_20261014_cyclic.py).

The second construction begins with 120 explicit points over Q(phi), where phi squared equals phi plus one. Their exact supporting tetrahedra form the 600-cell: 720 edges, 1,200 triangles, and 600 tetrahedra. The clique census is 120, 720, 1,200, 600, then zero, so its boundary is flag and has precisely 6,420 minimal nonface pairs. Supporting inequalities, strict inequalities at every other vertex, connected facet adjacency, and two tetrahedra per ridge certify the whole convex boundary. Central symmetry places the origin inside.

Taking the polar gives the simple 120-cell. If its inequalities are `V_i · x <= 1`, a wedge at facet i replaces that inequality by `V_i · x + z_i <= 1` and `V_i · x - z_i <= 1`. Five distinct wedges therefore give an explicit bounded simple 9-polytope with 125 facets. This establishes actual polytopality before any Coxeter labels are chosen.

The first recipe doubles the four vertices of a fixed tetrahedron and one further vertex. The 24 exact symmetries of its tetrahedron stabilizer divide the 116 other vertices into fourteen orbits. All fourteen are tested; this does not enumerate arbitrary five-vertex subsets of the 600-cell. The second recipe doubles the induced pentagon `(34,39,40,41,42)`, the link of edge `(0,33)`. It contains no tetrahedron, so is outside the first family.

For every recipe the expanded minimal nonfaces have orders two, three, or four. The complete actual vertex lists contain 17,908, 17,882, 17,850, and eleven times 17,810 vertices for the clustered family, and 17,800 for the pentagon. The total is **267,350 actual vertices across fifteen recipes**. Each nerve's minimal-nonface hypergraph is connected, so none is a product nerve.

Completeness does not depend on sampling nine-subsets of 125 facets. For any face, the base blocks included in full form a base clique. Its size is at most four. Since the five wedges contribute five extra possible entries, a face has at most nine entries. Equality forces a base tetrahedron and exactly one omitted member in each doubled block outside it. The independent checker reconstructs every vertex using this identity and checks all ridge incidences.

## Exact angle joins

A minimal nonface of order three or four must carry a Lannér diagram. The producer uses the repository's complete 108 labeled Lannér-four placements and its finite Coxeter diagram classifier. It adds a complete order-three alphabet: labels 2 through 6, and a category representing every finite integer label at least 7. There are 175 labeled triangle categories with reciprocal sum below one.

This last category is a safe exact compression for these tests. With a label at least 7, the other labels `(2,2)` give a spherical triangle for every finite value; every other pair has reciprocal sum at most `1/2 + 1/3 + 1/7 = 41/42`. A finite connected diagram with such an edge cannot have a third node. The category does **not** include parallel or dotted pairs. Missing pairs are left unconstrained; they never belong to an actual vertex.

The joins require shared meeting pairs to have consistent categories. At each actual vertex, unassigned meeting pairs are temporarily made right angles. This is a necessary relaxation: decreasing the nonnegative adjacency entries cannot destroy positive definiteness of a spherical Coxeter Gram matrix. Every surviving restriction is tested against the full finite-diagram classification. An empty state set therefore excludes every compact Coxeter labeling of that recipe without any floating-point angle tests, metric optimization, or volume heuristic.

Thirteen clustered recipes reach zero states in the [primary certificate](fresh_h9_highbase_construction_20261014_cell.json). The primary path for case 0 reaches its explicit 20,000-state limit; that path alone is inconclusive and is preserved as such. A [separate reordered join](fresh_h9_highbase_construction_20261014_refined.json), on the same incidence and with the same cap, gives **108 → 556 → 40 → 80 → 0**, with no truncation. The obstruction uses only six facets.

The [pentagon certificate](fresh_h9_highbase_construction_20261014_pentagon.json) reaches **108 → 8 → 0** after three mandatory Lannér-four placements, again without truncation. No full Gram candidate remains to test for Lorentzian signature or compact volume.

## Validation and limits

All five saved-result replays pass. The [independent stdlib checker](../scripts/fresh_h9_highbase_construction_20261014_independent.py) and its [certificate](fresh_h9_highbase_construction_20261014_independent.json) independently reconstruct the cyclic facets and weights, the 600-cell geometry and symmetry orbits, and every expanded vertex and ridge. It imports no producer or label-join code. It does **not** implement a second Lannér catalogue or a second label join; those are covered by the producer's deterministic exact replays. Source and artifact hashes bind the dependencies.

A [separate source and proof review](fresh_h9_highbase_construction_20261014_review.md) checks the label compression, necessary finite-vertex relaxation, and truncated-path scope. Its own Fraction checker reproduces all 175 triangle categories and binds the untruncated zero-state certificates; it does not duplicate the entire label join.

Run from the repository root:

```sh
python scripts/fresh_h9_highbase_construction_20261014_cyclic.py --verify
python scripts/fresh_h9_highbase_construction_20261014_cell.py --limit 100 --max-states 20000 --verify
python scripts/fresh_h9_highbase_construction_20261014_refined.py --verify
python scripts/fresh_h9_highbase_construction_20261014_pentagon.py --verify
python scripts/fresh_h9_highbase_construction_20261014_independent.py --verify
```

The complete claim is limited to the displayed three cyclic bases and fifteen 120-cell recipes. Other five-wedge subsets, repeated wedges, other bases, and constructions with ideal vertices remain outside this certificate. The results eliminate specific incidence choices; they do not lower the general possibility of a new compact H9 polytope to a proved impossibility.
