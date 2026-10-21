# Two coupled finite-edge caps on the completed E7 chamber

**Exact positive construction, still infinite volume.** Two mutually orthogonal norm-two caps extend the frozen [88-wall bridge](fresh_h20_e7_causal_residuals_20261014.md) to 90 actual facets. They preserve all 44 specified old finite vertices and both filled cusp links, close 31 outgoing edges, and create exactly 40 finite vertices incident to the new caps. Among these 40 and the old 44 vertices there are 145 finite edges with a new endpoint. The inherited strict escape survives, and a new parallel pair has an explicitly unbounded cusp link.

The [primary checker](../scripts/fresh_h20_e7_finite_edgewalk_20261014_exact.py) and [certificate](fresh_h20_e7_finite_edgewalk_20261014_exact.json) use self-contained exact recipes; they do not depend on the numerical scout. The [independent review](fresh_h20_e7_finite_edgewalk_20261014_review.md) reconstructs the geometry with separate sparse squarefree arithmetic and independently enumerates the finite supports. The [frontier checker](../scripts/fresh_h20_e7_finite_edgewalk_20261014_frontier.py) and [certificate](fresh_h20_e7_finite_edgewalk_20261014_frontier.json) give the precise remaining finite-edge and cusp status.

## Construction

Use the numbering of the completed 87-wall E7 chamber; its E7 roots are `[12,13,14,15,18,31,65]`. The frozen bridge is wall 87, and the two new walls are 88 and 89. Let `ω13` be the E7 fundamental weight dual to wall 13, so `ω13²=12`. Let `V4` and `V10` be the saved old finite points, with squares `−8/5` and `−32/15`.

For each point, remove wall 13 from its 20 incident walls. Let `p1` be the projection onto the first resulting positive 19-dimensional span with products −1 at walls 0 and 14 and zero at the other 17. Define `p2` analogously at `V10`, with products −1 at walls 10 and 14. Their exact squares are `47/30` and `22/15`; both have product `2/3` with wall 13. The new normals are

```text
r88 = p1 − (2/3)ω13 + (7/4)V4,
r89 = p2 − (2/3)ω13 + (3/2)V10.
```

Orthogonality of the three summands proves their norm-two identities directly, and exact physical pairings give `r88·r89=0`. These caps have nonright contacts on the shared edge walls, so they lie outside the earlier all-right transverse family. Their coefficients remain in `Q(√2,√79,√11)`; only the retained bridge uses the further quadratic extension.

All 177 new unordered pairs pass: 39 right, five π/3, one parallel and 132 dotted. Both caps are dotted to the bridge. A common strict timelike center remains. Since every off-diagonal pairing is nonpositive, projecting that center onto each wall gives a timelike point strict on every other wall, proving all 90 facets actual.

## Actual endpoints and complete new-cap vertex census

On an old edge releasing E7 wall `j`, use its exact dual `ωj` and point `V`. For a proposed cap `r`, the intersection is

```text
X = V − tωj,    t = (r·V)/(r·ωj).
```

The checker requires `t>0`, future timelike `X`, zero products with exactly the shared 19 walls and the new cap, and strict products with every other wall of the 90-wall chamber. The old vertex is strict at both new caps. Consequently the open segment from `V` to `X` is strict on every nonincident inequality, proving this is the first intersection and an actual finite edge endpoint. There are 19 such endpoints on wall 88 and 12 on wall 89. Their incident sets are distinct, and their new-cap incidence distinguishes them from every old vertex.

The full finite-neighbor graphs of the two new caps contain 3,080 possible rank-20 cliques involving at least one new cap. Exact rational Gram tests reject 3,040 and accept 40. Their future timelike projections are checked against all 90 inequalities. Thus the new-cap finite-vertex census is complete: the 31 adjacent endpoints and nine vertices incident to both new caps. This does not assert completeness of the inherited 88-wall finite-vertex census.

The 145 finite edges listed among the 84 verified vertices comprise 31 old–new edges and 114 new–new edges. The complete finite-endpoint status of the 800 directed germs at the 40 new vertices is also determined: 259 have another finite endpoint and 541 do not. Any germ retaining a new cap is covered by the complete new-cap vertex census; removing the sole new cap at one of the other 31 vertices gives its certified edge back to the specified old vertex. Ideal endpoints of the remaining 541 germs are not fully classified.

## Remaining volume defects

The unique new parallel pair is `(9,88)`. Its future null sum is feasible for all 90 walls and incident to exactly 20. The full horospherical normal span has rank 18, short of the required 19. The frontier certificate additionally constructs a positive-norm vector perpendicular to every incident normal and to the null ray. It gives a genuine lineality direction in the horospherical link, proving that link is unbounded. The two previously completed cusp links remain unchanged because both new walls are strict there.

The frozen bridge's rational strict escape remains strict at both caps. Its square is

```text
259854917001622621 / 62500000000000000 > 0.
```

Together with the common strict timelike center this gives a strict open set of ideal directions and proves infinite volume. The finite-edge improvement is therefore a local construction, not a finite-volume completion.

## Search scope and verification

The [exploratory scout](../scripts/fresh_h20_e7_finite_edgewalk_20261014_scout.py) used the 308 directions obtained by releasing one E7 wall at each of the 44 known vertices. It enumerated 15,198 finite-label profiles with exact inverse-Gram Schur signs, plus 272 symbolic `I2(m)` tails. Finite spherical diagrams are forests of maximum degree three, so a new node has at most three nonright contacts, in distinct old components; the exact Schur bound tests the remaining choices. Labels at least seven require a single old isolated A1 component. The numerical search sampled finite omitted-wall labels and several tail labels; it is not an all-parameter or all-label endpoint exclusion. It found 76,278 numerical recipes, none cutting the then-saved escape. The full numerical archive is deterministic gzip, [scout.json.gz](fresh_h20_e7_finite_edgewalk_20261014_scout.json.gz), and is not needed to replay the exact construction.

Exact construction and frontier generation and saved replays passed. The independent review is recorded in the linked review certificate; source/data manifests bind the final frozen inputs. No predecessor artifact was changed.
