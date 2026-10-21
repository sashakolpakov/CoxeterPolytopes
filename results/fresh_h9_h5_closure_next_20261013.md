# Coupled H5 endpoints close a compact hexagonal face

Two exactly certified caps retain the previous 123 walls and close an actual compact H2 hexagon inside the distinguished H5 face. The resulting 125-wall chamber still has an exact strict spacelike escape. This is a local face closure, not a finite-volume H9 example or a compact H5 face.

The frozen [certificate](../data/enumerations/fresh_h9_h5_closure_next_20261013_verify.json) is replayed by the [verifier](../scripts/fresh_h9_h5_closure_next_20261013_verify.py), using the discrete [recipes](../data/enumerations/fresh_h9_h5_closure_next_20261013_recipes.json). The [previous checkpoint](fresh_h9_h5_closure_20261013.md) supplies the 123 retained walls and their verified relations. A separate [160-digit reconstruction review](fresh_h9_h5_closure_next_20261013_review.md) checks the new formulas and strict inequalities with the same interval engine; it is not an independent implementation of all inherited roots.

## The two metric constraints

Both new normals have squared norm two and are right to the four walls defining the H5 face: `seed2`, `compact1`, `compact2`, `compact3`.

| New wall | Additional right contacts | Additional angle pi/3 contacts |
| --- | --- | --- |
| `endpoint0` | `crosspair0`, `roof6` | `nextH4`, `wall17`, `secondAdaptive0` |
| `linked0` | `crosspair0`, `nextH4`, `roof5`, `roof6` | `endpoint0` |

All other retained-wall relations are strictly dotted. These are 247 new pair checks: 14 right, four pi/3, and 229 dotted. The first cap cuts the previous strict escape; its covector pairing is approximately -0.452755. The second cap is chosen to connect endpoints, rather than to cut the next escape.

For each cap, nine prescribed pairings determine an orthogonal particular solution and a one-dimensional complement. The selected negative square-root branch has norm two by the projection identity. Invertibility, positive radicand, orientation, and every strict dotted inequality have rational interval certificates. Exact inherited finite Gram entries are used directly. The complement norm is evaluated as `c^2 - cp^T M^-1 cp`, avoiding cancellation from squaring its coordinate expression. Finite-angle equalities follow from the defining linear equations, not from overlapping intervals.

The complete finite-neighbor subset check considers 512 subsets for each cap. Totals are 320 subsets excluded by a dotted pair, 672 positive definite, 32 indefinite, and zero positive affine diagrams. Any newly introduced affine diagram would occur in this census. Together with the inherited audit, the new chamber remains free of positive affine subdiagrams.

## The compact hexagon and the remaining gap

The compact H2 face is perpendicular to

`seed2, compact1, compact2, compact3, crosspair0, nextH4, roof6`.

Its six boundary sides, in cyclic order, are

`seed6 — wall17 — endpoint0 — linked0 — roof5 — roofPolar — seed6`.

The checker tests every other wall for incidence with this face and certifies exactly these six facets. Each adjacent pair gives an actual finite vertex: its physical rank-nine Gram is positive definite, and every nonincident inequality is strict. Every side has two finite endpoints and the complete boundary graph is one six-cycle. Thus the convex H2 face is compact. No claim is made that its intrinsic angles are Coxeter angles.

The H5 census tests all 65,780 five-facet subsets: 65,738 contain a dotted pair, 29 are positive definite, and 13 are indefinite. The finite vertices form one connected component. The progression is:

| Checkpoint | H5 facets | Finite vertices | Finite edges | Open edges |
| --- | ---: | ---: | ---: | ---: |
| Previous 123-wall chamber | 24 | 26 | 42 | 46 |
| New 125-wall chamber | 26 | 29 | 46 | 53 |

The larger open-edge count is why this checkpoint stops at the coherent two-cap face closure. It does not show an improvement toward global compactness.

The saved future spacelike vector is strict against all 125 walls. Its exact norm is

`(844076039855942916265101992852977 - 521662187518674164182600834515677 phi) / 160000000000000000000000000000`,

approximately 0.0555612819184. Projection into the distinguished H5 slice remains spacelike, with norm approximately 0.0555611868315, and is strict against all its nonfixed walls. Both the full H9 chamber and the H5 face therefore have infinite volume. The next numerical gap is supported by `nextH4`, `compact0`, `roof0`, `roof1`, `closure0`; it is preserved as a target, not as a completeness assertion.

## Replay and scope

```sh
python scripts/fresh_h9_h5_closure_next_20261013_verify.py --verify
```

The frozen replay passes. It imports no numerical scout, anchors predecessor certificates by hashes, and checks the cap identities, all new pairs, positive-affine census, actual finite vertices, complete hexagon boundary, and both escapes. Numerical candidate generation and selection were bounded; no exhausted-family or global nonexistence claim follows from them. The new compact hexagon is a geometric face inside an infinite-volume H9 Coxeter chamber, not the requested new compact H9 polytope.
