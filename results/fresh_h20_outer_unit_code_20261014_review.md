# Independent review of the outer-code layers

**PASS.** The compressed family formulas, new pair counts, retained endpoint geometry, old escape cut, and new strict escape are correct. The construction has **7,545,383 actual facets and infinite volume**. The new families also satisfy a structural obstruction: they are all dotted to the original central facet, which has no finite vertices but has actual `pi/5` ridges. Thus arbitrary further additions dotted to that facet cannot yield finite volume.

The [separate checker](../scripts/fresh_h20_outer_unit_code_20261014_review.py) imports no producer or project arithmetic. It uses standard-library `Fraction` arithmetic in the basis `1,sqrt(5)`, independently reconstructs the binary parity-check code and its holes, and recomputes all compressed new-new pair counts with coordinate generating polynomials. Its [saved certificate](fresh_h20_outer_unit_code_20261014_review.json) binds the frozen source, data, inherited local audit and geometric obstruction note. Generation and saved replay pass. The [primary producer](../scripts/fresh_h20_outer_unit_code_20261014.py) also passes a fresh saved replay.

## Counts and all-family relations

The zero masks form a 32-word even code with minimum distance four. The 192 odd holes are exactly the masks at distance three from that code. Their weight distributions and the number of nonzero sign coordinates give the following distinct normal counts.

| Family | Norm squared | Actual facets |
|---|---:|---:|
| Code | 2 | 1,139,200 |
| Holes | 1 | 5,357,568 |
| Small half-unit centers | 1/2 | 524,288 |
| Retained binary orbit | 1 | 524,288 |
| Retained weighted core | mixed | 39 |

The old full vertex-cap bank is the code's zero-word subfamily and is counted once. With matching outer signs, each common nonzero unit coordinate contributes `2+2t^4` to the squared-distance generating polynomial and each single nonzero coordinate contributes `2t`. This gives the complete multiplicities independently of the producer's sign-count formula. The six pair-class tables agree exactly with the primary archive.

An outer sign change adds at least eight to the squared distance. The one equality case is explicitly retained: **code centers with identical unit coordinates and only the heavy sign reversed are parallel**, not dotted. All other outer-sign changes are strictly dotted. The independent checker includes the resulting 1,139,200 ordered extra parallel pairs before dividing same-family counts by two. Small-family cross counts use distance `5/2+2r`; these reproduce the right, `pi/4`, parallel and dotted classifications. There are no unclassified positive or forbidden-angle pairs.

New-seed relations are right, `pi/4`, parallel or dotted as recorded; every new-central relation is strictly dotted. Exact distance lower bounds show that every new-binary relation is strictly dotted, for both binary full-unit triples and all sign choices. This verifies **3,681,055,408,128** such relations by bounds rather than literal enumeration. The complete retained binary old-old pair proof is inherited; it is not repeated by this review.

## Geometry and escape

The common point `(0,100,1)` is strictly timelike and strict on every wall. A finite collection of distinct spacelike normals with a common strict center and nonpositive cross pairings has every wall as an actual facet: projecting the center onto that wall remains timelike and strict on all the others. This proves the actual-facet claim without a new vertex census.

All 524,288 original box-vertex endpoints remain finite. Their norm is `-2`; their incident walls are the nineteen corresponding seed walls and their own norm-two code cap. This is a positive-definite orthogonal rank-20 incidence. The other code caps have positive squared distance at least four, while the hole and small-cap norm squares are below two, so every other new inequality is strict there. The retained binary and central inequalities are also strict. This is preservation of those specified endpoints, not an exhaustive corner census.

The all-zero-unit codeword cuts the previously saved escape exactly. The new rational-coordinate escape has norm `1/4`, minimum normalized seed gap `1/100`, and the following minimum ball powers:

| Family | Minimum power |
|---|---|
| Code | `4813/10000 + phi/1250` |
| Holes | `5013/10000 + phi/1250` |
| Small | `4713/10000 + phi/1250` |
| Binary, relaxed sign minimum | `26313/10000 + 2401phi/1250` |
| Central | `102413/10000 + 9801phi/1250` |

Every entry is strictly greater than `1/4`. Therefore every cap pairing is negative. Its pairing with the common center is negative as well. This is a strict future spacelike escape for the entire literal chamber, proving infinite volume.

## Completion obstruction

The [dotted-only completion proof](fresh_h20_weighted_binary_20261014_dotted_completion_review.md) applies to finite-sided, full-dimensional, acute-compatible Coxeter completions retaining the stated facets. It uses a strict interior point of the **proposed completion** to orient projected dotted normals. Each is future timelike in the retained facet plane, so its inequality is strict on that entire hyperbolic plane; such additions do not alter the facet or create finite vertices on it.

A finite-volume Coxeter ridge of finite label outside `{2,3,4,6}` cannot meet an ideal vertex. On a bounded horospherical Coxeter link, the product of its two reflections would have irrational linear trace, contradicting preservation of the Euclidean translation lattice. Thus such a ridge must have a finite vertex. This is a statement about original polytope ridges, not a claim that every facet inherits Coxeter angles.

The central facet of the 39-wall core already supplies the contradiction. On its plane `D=2,U=1`, a finite vertex would need all nineteen seed coordinates to be `+1` or `-1`, but the resulting norm is `16+8phi>0`. It has no finite vertices. For any positive golden seed wall, `(e_g,2,1)` is an actual common ridge point of norm `phi-3<0`, strict on all other core seed walls; the ridge angle is `pi/5`. Every binary and every new code, hole or small wall is dotted to this central facet. The independent checker verifies these literal norm identities.

The older binary-facet version is also valid. As a direct actuality check, projecting the center onto a reference binary/golden pair gives

```text
p = h + (167+173phi/2)r + (169/2+173phi/4)w,
p^2 = -141089/8 - 149517phi/16 < 0.
```

Both coefficients are positive, so all other inequalities remain strict. The independent checker verifies the orthogonalities and norm, while the absence of finite vertices on binary facets is inherited from the [separate local audit](fresh_h20_weighted_binary_20261014_local_review.md).

Neither obstruction excludes replacing a retained facet or adding walls that meet it. No new complete cusp or finite-vertex census is claimed. The constructive requirement is now explicit: new cap layers dotted to both retained obstructed facets cannot finish the chamber, regardless of their size.
