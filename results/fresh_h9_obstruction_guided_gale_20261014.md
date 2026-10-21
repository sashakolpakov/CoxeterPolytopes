# Obstruction-guided rank-five Gale search

This bounded continuation found **74 new realizable fifteen-facet nerves**, all outside the standard single-face-truncation family and outside the prior 1,063 saved canonical keys. Every new nerve is exactly excluded: 70 by a literal inherited local pattern and four by the two-Lannér-block gate. No full ambient label assignment survived to metric solving. The result concerns these saved proposals, not all fifteen-facet polytopes.

The main outputs are the [guided-walk ledger](fresh_h9_obstruction_guided_gale_20261014.json), [positive-circuit ledger](fresh_h9_obstruction_guided_gale_20261014_circuits.json), [separate incidence audit](fresh_h9_obstruction_guided_gale_20261014_independent.json), and [final certificate](fresh_h9_obstruction_guided_gale_20261014_certificate.json). Both search state files are terminal; no job remains running. All preceding archives were preserved.

The search objective changed materially. The [new support collector](../scripts/fresh_h9_obstruction_guided_gale_20261014_supports.py) scans unions of two order-three-through-five missing faces and records every recognized literal obstruction support on those carriers. It is a finite inherited dictionary, not a complete oracle for all diagram obstructions. Each proposed mutation targets rows in a failing support; the score counts these supports, with penalties for large missing faces and simplex-boundary link factors. Temporary invalid states are allowed, but only fully checked admissible nerves are saved.

The first batch used four different prior seeds, cases `908,638,748,661`, and exactly 3,000 proposals per seed. For 10,772 proposals, a rank-one two-row update moved just beyond an exact five-minor sign boundary; 1,228 proposals were longer targeted jumps. These are proposal heuristics, not mathematical exclusion rules. The lowest attained admissible support count fell from the selected seed's eight to four.

| First-batch outcome | Count |
|---|---:|
| Order-greater-than-five minimal nonface | 10,187 |
| Simplex-boundary factor in a vertex link | 315 |
| A label is not an actual vertex | 4 |
| Eligible proposals, including repeats | 1,494 |
| Total proposals | 12,000 |

The eligible proposals comprised 122 old known keys, 1,300 repetitions within the new batch, and **72 new saved nerves**. All 72 were exactly obstructed, and their separate integer realization/application replay passed.

The second batch addressed a specific weakness of minor-sign walking: many sign changes do not alter the actual facet structure. It instead forced a positive dependence among five selected Gale rows, while compensating another row to preserve the zero sum, and perturbed to both sides. Only boundaries whose complementary ten rows split into two groups of five were proposed. The targeted old missing face had to lie on one side of this split. This constructs a potential change of a ten-vertex facet triangulation; it does not guarantee that the rest of the nerve keeps small missing faces.

Of 2,000 bounded carrier attempts, 648 met this balanced positive-circuit condition. Their two perturbation sides gave 1,296 exact proposals: 1,293 had a large missing face, one had a truncation link, and **two were new admissible nerves**. Both were locally obstructed. The final independent certificate reconstructs the literal five-row positive dependence, the complementary sign partition, the perturbation, and both resulting Gale configurations.

A separate [exact flip diagnostic](fresh_h9_obstruction_guided_gale_20261014_flip_diagnostic.json) explains the same bottleneck combinatorially. For the twenty selected low-obstruction nerves, it exhausts the valid 5↔5 bistellar recipes using an exact star equality. There are 625 recipes: 543 introduce an order-greater-than-five missing face, and all 82 remaining recipes retain a literal inherited obstruction. These newly formed complexes are simplicial spheres by the star replacement; no Euclidean polytopality is claimed for them. They are **not** included in the 74 realized nerves, and this is not a census of all flip sizes or all parents.

Every saved realized nerve has all 3,003 five-row Gale minors nonzero. Positive six-row circuits supply a strictly positive kernel vector and the exact nine-vertex facets of a rational simplicial nine-polytope. Every label is an actual vertex, every ridge has two incident facets, and all minimal nonfaces have order at most five. No vertex link has an isolated minimal-nonface component, so none can be the stellar subdivision of a face of any **simple fourteen-facet parent**. Literal minimal-nonface sets and saved colored-clutter `pynauty` canonical keys are distinct; there is no independent reimplementation of graph-isomorphism canonicalization.

The arithmetic review reuses the separately implemented integer Bareiss/full-face/link-ideal checker, not the proposal geometry engine. Across the 74 new realized cases it checks **222,222 five-minors, 370,370 six-circuits, 1,110 vertex links, and 1,226 archived guidance witnesses**. Each local obstruction is applied only after its entire induced nerve is identified explicitly. Its earlier complete label/determinant proof remains a dependency; this round does not claim a fresh re-enumeration of those libraries. The inherited proofs include every finite integer label, including the rigorously quantified category for labels at least seven. No fourteen-facet-specific restriction is imposed.

The first 72 cases completed their separate arithmetic replay. The final command freshly rechecks the two positive-circuit cases and binds the complete first audit, all saved case hashes, the terminal states, and the flip diagnostic. Generation and saved-result replay passed:

```sh
PYTHONPATH=src:. python scripts/fresh_h9_obstruction_guided_gale_20261014_certificate.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_obstruction_guided_gale_20261014_flip_diagnostic.py --verify
```

No full-label search, metric fit, or node-limited case was reached. The remaining research frontier is unsampled incidence, rather than an unresolved label assignment in this package. The diagnostic cautions against simply extending the same 5↔5 walk: a different nonuniform Gale arrangement or a controlled sequence through larger missing faces would need a new realization and obstruction audit. No compact H9 example or global fifteen-facet impossibility theorem follows from this checkpoint.
