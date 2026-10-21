# Independent review of the symmetric bridge-edge closure

**PASS: independent generation and saved replay.** The separate [checker](../scripts/fresh_h20_bridge_edge_closure_20261014_review.py) and [certificate](fresh_h20_bridge_edge_closure_20261014_review.json) reconstruct the two new normals, all 193 new pairings and the two physical endpoint vectors from the frozen [construction data](fresh_h20_bridge_edge_closure_20261014.json). The [construction report](fresh_h20_bridge_edge_closure_20261014.md) correctly retains infinite-volume scope.

The review imports no producer arithmetic. It uses the earlier independent sparse Fraction multiplication in the 16-element multiquadratic basis, with three successive positive quadratic extensions and outward rational sign bounds. All sign decisions resolve at 192 bits. The saved artifact binds sources and data, allowing report links to be added without changing those bindings.

Let `S` be the PD18 right support of the old norm-two cap `r`, let `ω` be the support fundamental weight at wall61, and let `X` be a future timelike point on the old edge. Thus `ω²=2/3`, `X·S=X·r=0`, and `A=-X²>0`. If `W` is perpendicular to `X,S,r`, with `B=W²>0`, put

```
Δ = A(m²B − 5/6),
R± = A(−ω − r/2 ± mW) + sqrt(Δ) X,
P± = mB X ± sqrt(Δ) W.
```

For `m>0` and `Δ>0`, direct orthogonal expansion gives

```
R±² = 2A²,
R±·r = −A,
R±·wall61 = −A,
R+·R− = A²(2 − 2m²B),
P±² = −5AB/6,
R±·P± = 0,
R∓·P± = −2Am sqrt(Δ) B < 0,
P±·X = −AmB < 0.
```

These identities certify the normalization, future endpoint orientation, and strict separation of the two endpoints once all remaining wall inequalities are checked. The 20 incident roots have positive-definite Gram: after removing `S`, the normalized remaining two-root Gram is `[[2,-1],[-1,4/3]]`, with final Schur pivot `5/6`.

Two distinct feasible finite vertices sharing the old independent PD19 support `S+r` bound its actual finite edge. Along that one-dimensional hyperbolic line, each cap supplies the opposite endpoint. Linear feasibility gives the joining segment; the exact endpoint incidences and strict opposite-cap signs identify it as the full edge.

The induced angles in the old H2 support face need not themselves be Coxeter angles. Here their squared cosine with the old geodesic is `3/8`. The ambient angles are π/3; the nonzero contact to the support is essential to this distinction.

## Independent arithmetic and local census

The checker verifies all 98 normal squares, the 193 new pairings (`34 right + 4 pi3 + 155 dotted`), the retained strict future center, both normalized PD20 endpoint Grams, all 196 endpoint inequalities, and distinctness via the opposite-cap signs. It separately constructs actual new facet points without field division and verifies every other inequality strictly. The normalized endpoint diagrams have path components of sizes `[1,1,4,4,10]`, hence type `A10+2A4+2A1`; a separate literal graph check confirmed the paths.

The [small census checker](../scripts/fresh_h20_bridge_edge_closure_20261014_review_census.py) and [certificate](fresh_h20_bridge_edge_closure_20261014_review_census.json) give a complete finite-vertex statement on the three walls added after the 95-wall base. Wall95 has precisely the 20 old finite-angle neighbors `S+[0,28]`. Among them, the disjoint pairs `(61,28)` and `(87,0)` are strictly dotted. The first pairing is exactly `-11/2`; both inequalities are independently reconstructed. Every 19-subset retains at least one dotted pair, so wall95 had no finite vertex in the 96-wall chamber.

Each new wall has exactly 19 finite-angle neighbors, namely `S+[95]`; its mate and all remaining walls are strictly dotted. Its unique possible finite vertex is the certified endpoint. Therefore **exactly two finite vertices occur on the union of walls95,96,97 in the 98-wall chamber**. This is not a global old-vertex census. The lower bound of 116 previously certified vertices is inherited and, by the reviewed acute-extension lemma, yields at least 118 known vertices after the addition.

All three old bounded cusps have strictly negative pairings with both new walls. The independent checker also rechecks the retained rational escape against all 98 walls. Its positive square remains

`20181174579467302741 / 10^18`.

Thus the old edge is closed by two finite vertices while global infinite volume persists.

## A second entire H2 plane survives

The exact review additionally projects walls95,96,97 off the alternative support `S75`, obtained by replacing61 with75 in `S`. Each has strictly negative projected square and future orientation, checked using the projection of the enlarged strict center. The [earlier independent whole-face certificate](fresh_h20_bridge_faces_20261014_review.md) already proves that every nonsupport old95 projection is future timelike. Consequently **the entire S75 hyperbolic plane is still present in the 98-wall chamber**, with no induced boundary from any of these three added walls.

This gives a concrete remaining geometric obligation, independently of the retained strict escape. It does not rule out later walls having spacelike projection onto that plane.

There is also a precise necessary contact condition for such a later wall. The saved inverse support Gram has diagonal entry `5` at wall75 and is entrywise nonnegative. For a cap of square `q>0`, any nonright Coxeter contact to this norm-two wall has squared pairing at least `q/2`. Its projection cost is therefore at least `5q/2`, giving strictly negative residual square. **A cap that induces a boundary geodesic on this retained H2 face must be right to wall75.** This statement permits other support contacts and is conditional on retaining the support; it is not a prohibition on completing the chamber by a different family.

## Scope and reproduction

The full old96 pair, vertex and bounded-link censuses are inherited dependencies. New pairs, actual points, the escape and three extra whole-face projections are independently rebuilt. No complete98 vertex census, finite volume or compactness is asserted.

```
python3 scripts/fresh_h20_bridge_edge_closure_20261014_review.py --verify
python3 scripts/fresh_h20_bridge_edge_closure_20261014_review_census.py --verify
```

No issue was found in the symmetric construction, normalization, endpoint geometry or volume statement.
