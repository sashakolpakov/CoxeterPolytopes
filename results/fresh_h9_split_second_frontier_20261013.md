# Five exact face closures beyond the 74-wall chamber

Five compatible new caps retain every wall of the [74-wall split-seed checkpoint](fresh_h9_split_endpoint_walk_exact_20261013.md). The resulting **79-facet Coxeter chamber** has sixteen certified finite vertices in the two reached components. Ten previously open edges close, producing **two compact triangular faces and three compact quadrilateral faces**. The exact frontier still has 106 open ends, and twelve strict spacelike witnesses prove infinite volume. No compact or finite-volume example is claimed.

The starting case is pair `walk2720,walk6746`, case 0 of the preceding ten-cap family. This continuation follows outgoing edges from the newly created finite vertex at `walk2720`; it does not repeat the previous ten caps or delete an existing facet.

## Construction and coherent selection

Each of the eight open outgoing edges has an exact positive-definite rank-eight wall span. The proposal stage prescribes up to three non-right attachments with labels 3 or 5, requires the nine-wall endpoint diagram positive definite, and uses an additional label 2, 3 or 5 to an anchor in the existing vertex field. It finds 270 numerical proposals compatible with all 74 walls. Five mutually compatible proposals close ten distinct old open rays, two per cap.

The exact equations are solved over the degree-at-most-eight field of `walk2720`. Nine linear pairing equations leave `p=u+t v`; the norm-two equation is quadratic over that field. Thus the new normals have degree at most sixteen over the rationals. Every numerical branch proposal is checked by exact norm, center, pairing and actual-vertex certificates.

| New cap | Non-right couplings within its eight-wall edge | Outside anchor | New compact face |
|---|---|---|---|
| `second33` | `wall22`: label 3 | `seed0`: label 3 | triangle |
| `second45` | `seed5`: label 5; `seed7`: label 3 | `seed4`: label 5 | triangle |
| `second4` | `seed4,walk2720`: label 3 | `seed2`: label 2 | quadrilateral |
| `second136` | `seed4,walk2720`: label 3 | `seed2`: label 5 | quadrilateral |
| `second5` | `seed3`: label 3 | `seed2`: label 2 | quadrilateral |

The other prescribed edge couplings are right angles. The five caps also pass every mutual pairing gate; `second33` and `second45` meet at `pi/10`. The complete chamber has 3,081 exact pair checks and an empty connected affine-subdiagram census after 1,579 exact frontier checks. A common strict timelike center certifies full dimension. Acuteness then proves each selected inequality supports a genuine facet by projecting that center onto the corresponding wall.

The [cap archive](../data/enumerations/fresh_h9_split_second_frontier_20261013_caps.json) saves the complete tower fields, line equations, norm quadratics, branches, all retained pairings and an explicit actual finite-vertex point for each cap. The [geometry archive](../data/enumerations/fresh_h9_split_second_frontier_20261013_geometry.json) saves their joint 79-wall chamber and its escapes.

## Exact vertex and edge frontier

The [exact frontier certificate](../data/enumerations/fresh_h9_split_second_frontier_20261013_frontier.json) upgrades the numerical active-set proposals to exact finite vertices and exact endpoint classifications.

For a proposed nine-wall vertex, its physical Gram `M` is positive definite and `M^-1` is entrywise nonnegative. Projecting the common strict timelike center onto those nine mirrors therefore gives a future timelike point. It has exactly those nine zero inequalities; every other wall remains strict because all cross pairings are nonpositive. This proves an actual simple finite vertex without constructing an unnecessarily large compositum for its coordinates.

For every resulting rank-eight spherical edge, the checker tests every remaining wall as a possible ninth wall. A dotted pairing with an edge wall rules out a positive-semidefinite extension immediately. All other extensions receive an exact inertia check over `Q(phi,sqrt(2+phi))`, which contains every finite Gram entry here, including label 10. A finite endpoint requires a positive-definite extension; an ideal endpoint with an active extra wall requires a positive-semidefinite singular extension. The exhaustive check finds no such singular extension.

Consequently an edge with only its starting finite endpoint reaches the ideal boundary with every non-edge wall strict. These classifications are exact and do not depend on floating intersection times.

| Reached frontier | Starting 74 walls | Final 79 walls |
|---|---:|---:|
| Actual finite vertices | 6 | 16 |
| Edges joining two finite vertices | 4 | 19 |
| Edges with one finite endpoint and an open ideal end | 46 | 106 |
| Connected finite-vertex components | 2 | 2 |
| Independent finite-graph cycles | 0 | 5 |
| Exact ninth-wall tests after dotted exclusions | 54 | 216 |

The final five cycles each lie in the intersection of seven independent spherical walls, hence on an actual two-dimensional face. Their three or four finite edges close the convex face completely. The archive gives the seven common walls and ordered vertex wall sets for each polygon. These are compact faces of the H9 chamber; no claim is made that the induced polygon angles are themselves Coxeter angles.

The frontier census is complete for the two reached finite-vertex components. It does not exclude additional disconnected components. The sixteen vertices are thus a certified lower bound for the whole chamber, not a global vertex census.

## Remaining gap and verification

The first saved strict future spacelike escape has norm

`(11764468488123570049 - 6234856568308221881 phi) / 16000000000000000000 > 0`.

Its unperturbed proposal has support `[0,4,6,9]` in the independent basis `[0,1,2,3,4,5,6,7,8,12]`. It is strict on all 79 walls. The exact open-edge census gives another geometric explanation for the remaining infinite volume: closing selected old ends has exposed further outgoing edges.

The original strict escape of the 74-wall starting chamber also remains strict after all five additions; the replay checks this directly. The face closures therefore improve the resolved local frontier without closing that already identified gap.

A separate bounded scout at `walk6746` yields 371 numerical no-removal proposals. None of those, nor any of the first 270 proposals, cuts the remote large escape supplied by the parallel search. Their largest numerical pairings there are approximately `-1.30926` and `-1.94982`, respectively. This is a numerical diagnostic for the stated scouts, not a universal obstruction or an exact theorem about all their continuous cap families. No additional caps from that opposite-side pass are claimed certified here.

The [producer-free replay](../scripts/fresh_h9_split_second_frontier_20261013_verify.py) and saved-result comparison both pass. They recompute all five caps, all 79 normal/Gram/norm/center identities, every chamber pair, the complete affine census, ten explicit finite-vertex point checks, the exact vertex/edge frontier and all twelve strict escapes. Mutations of a cap coefficient, an escape to a timelike point, and an open edge to a claimed closed edge are rejected.

A separate source review found no issue in the positive-Gram projection proof, the exhaustive ninth-wall endpoint criterion, the reached-component scope, or the compact two-face argument. That review is mathematical and source-level; it is not represented as a second run of the producer.

The [independent geometry replay](fresh_h9_split_second_frontier_independent_geometry_20261013.json) also passes using a separate implementation of the tower arithmetic and frozen source extraction. It checks all 79 normal/Gram/norm/strict-center identities, the first strict spacelike escape, and the `second33` actual finite vertex with all 79 inequalities and its physical positive rank-nine Gram. Its scope is deliberately narrower than the complete angle, affine and edge census above.

```sh
PYTHONPATH=src:. python scripts/fresh_h9_split_second_frontier_20261013_verify.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_split_second_frontier_20261013_frontier.py --verify
```

The archive-only `selected()` loader in [the exact base module](../scripts/fresh_h9_split_second_frontier_20261013_base.py) returns `G,I,h,profiles,geometry` for the frozen 79-wall chamber. Its [physical source snapshot](../data/enumerations/fresh_h9_split_second_frontier_20261013_source.json) makes the numerical scout reproducible in the same coordinate frame.

This is a constructive finite-angle edge and face continuation. It is not a complete Allcock lattice algorithm: no first-root theorem, exhaustive residual-coset enumeration, or unrestricted metric branching is asserted. The five face closures are exact progress; the chamber remains infinite-volume.
