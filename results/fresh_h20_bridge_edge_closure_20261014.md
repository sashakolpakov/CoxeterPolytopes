# Two nonright caps close a bridge-face edge

The [exact certificate](fresh_h20_bridge_edge_closure_20261014.json) retains the [96-wall chamber](fresh_h20_e7_bridge_face_20261014.md) and adds two compatible caps. They close its specified open edge at **two distinct finite vertices**, both of type `A10+2A4+2A1`. All 98 listed walls are actual facets. The old rational strict escape survives, so the chamber still has infinite volume.

There are 193 new pairs: 34 right, four π/3, and 155 strictly dotted. Each cap has π/3 contacts with old walls 61 and 95, is right to the other 17 roots of the previous edge support, and is dotted to every other wall, including its mate. This is a bounded two-cap construction, not a further sampled cap bank.

## Exact symmetric construction

Write `S` for the 18-root support `S61` of the preceding construction and `r` for its norm-two wall 95. The dual weight `ω` in `span(S)` at wall 61 satisfies

`ω·N61=1`, `ω·Ni=0` for the other support roots, and `ω²=2/3`.

Let `X` be the saved future timelike point on the old edge, perpendicular to `S` and `r`, and put `A=−X²>0`. With `A0` the projection of old wall 0 off `S`, define

`W = X²(A0+r/2) − (A0·X)X`.

The checker verifies `W²>0` and `W⊥X,S,r`. It chooses a positive integer `m` by exact tests and defines

`Δ=A(m²W²−5/6)`,

`R±=A(−ω−r/2±mW)+√Δ X`.

The accepted choice is `m=1`. Here `Δ>0`, the new normals have square `2A²`, and dividing each by the positive scalar `A` gives norm two. The scaled coordinates avoid unnecessary field inverses. They use the inherited exact field and one additional positive square root; the JSON records the encoding.

The norm and contact identities follow directly from the orthogonal decomposition:

`(−ω−r/2)²=7/6`, `R±²=2A²`,

`R±·N61=R±·r=−A`,

`R+·R−=A²(2−2m²W²)`.

Every remaining new pairing is checked exactly with outward rational sign enclosures. The inherited strict center is strict at both caps, so it is a common strict future timelike point for all 98 inequalities.

## Actual finite endpoints and edge

The two endpoint representatives are

`P±=mW²X±√Δ W`.

They satisfy

`P±²=−5AW²/6<0`, `P±·X=−AmW²<0`,

`R±·P±=0`, `R∓·P±=−2Am√Δ W²<0`.

Thus both are future timelike and strictly inside the opposite cap. The checker verifies every old inequality: exactly the shared 19 roots `S∪{95}` vanish, and all other old walls are strict. Adding the corresponding new root gives a rational normalized 20-root Gram with final Schur pivot

`2−2/3−1/2=5/6>0`.

These are actual finite vertices, rather than merely compatible diagrams. Their strict opposite-cap pairings prove distinctness. The segment between them is feasible and has exactly the shared positive-definite rank-19 support active in its interior, proving a bounded finite edge.

All previously known finite vertices were strict at wall 95, so neither endpoint duplicates them. The [acute-extension preservation lemma](acute_extension_preservation_20261014.md) retains at least 116 old finite vertices, giving **at least 118 known finite vertices** in the new chamber. It also makes all 98 walls actual facets. The three old proper bounded cusps remain strict at both new caps; those six signs are explicitly checked. No complete global vertex or edge census is claimed.

The separate [independent review and small census](fresh_h20_bridge_edge_closure_20261014_review.md) prove that these are exactly the two finite vertices on the union of walls 95, 96 and 97. Wall 95 had no finite vertex before this addition: its possible old neighbors contain two disjoint dotted pairs, which no 19-subset can avoid. Each of the two new walls has exactly the 19 shared support roots as finite-angle neighbors, leaving only its certified endpoint. This local completeness statement does not enumerate the other old vertices.

## Global status and verification

The previous rational witness remains strictly negative at both new walls and at all 96 old walls. Its square is unchanged:

`20181174579467302741 / 10^18 > 0`.

It is future oriented and spacelike, so the enlarged chamber still has an open ideal-boundary region and infinite volume. Closing this one edge has not completed the hyperbolic-plane face or the ambient chamber.

The independent review identifies a further untouched component: the alternative support `S75`, obtained by replacing wall 61 with wall 75 in `S`, still bounds an entire `H2` plane. All three walls added after the 95-wall base have strictly future timelike projections onto that plane. Its inverse support Gram gives a useful necessary condition: any later acute Coxeter wall that induces a boundary geodesic there must be right to wall 75. Other support contacts remain possible; this is a constraint on the retained face, not an exclusion of all later completions.

The [producer/verifier](../scripts/fresh_h20_bridge_edge_closure_20261014.py) binds the frozen 96-wall construction, its 95-wall predecessor, and the independent whole-H2 geometry certificate. All acceptance tests are exact; no numerical angle or norm test determines acceptance.

```
PYTHONPATH=src:. python scripts/fresh_h20_bridge_edge_closure_20261014.py --verify
```

Generation and saved-result replay passed. The [independent new-cap review](fresh_h20_bridge_edge_closure_20261014_review.md) and its small census checker also passed generation and saved replay. Separate sparse Fraction arithmetic in the 128-coordinate radical extension reconstructs all 193 new pairs, both physical endpoints and their 196 wall inequalities, actual facet points, cusp signs and the escape. The predecessor's full pair and vertex/link certificates are inherited rather than repeated as independent work here.
