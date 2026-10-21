# Exact finite-vertex additions to the split-seed chamber

Ten targeted caps preserve all 72 walls of the selected split-seed chamber and cut its saved leading escape. Their exact compatibility graph is `K6,4`, giving 24 maximal two-cap selections. Each selection has **74 genuine Coxeter facets and six certified actual finite vertices**, with no positive affine principal subdiagram. Nevertheless, one exact strict spacelike point satisfies **all 82 inequalities** from the original chamber and all ten caps. Every selection from this finite family therefore still has infinite volume. No compact or finite-volume example is claimed.

This is a construction with coherent finite-angle attachments, not an exhaustive implementation of Allcock's algorithm. It neither exhausts the edge-cap families nor closes the broader changed-seed search.

## Exact construction

The original sixteen seed normals have rank ten. The independent basis used here is `[0,1,2,3,4,5,6,7,8,12]`; its exact Gram has inertia `(9,1)`. The common center is the sum of the sixteen seed normals, with norm `-64-128phi < 0`. Coordinates are saved as normal coefficients `c` and covectors `p=Gc`. With `h=-G(center)`, the chamber inequalities in the dual affine section are `c·y >= 0`, and every retained wall has `c·h > 0`.

Each chosen rank-eight spherical edge block and its outside anchor lie in the same quartic field. Nine exact pairing equations leave a line `p=u+t v`. Substitution into `pᵀG⁻¹p=2` gives a quadratic over that quartic field. The numerical scout selects a proposed real branch; exact arithmetic then verifies the branch, norm, all sixteen seed pairings, common-center sign, target cut, and every relation to the 72 retained walls. The resulting normals lie in extensions of degree at most eight over the rationals.

For example, scout cap `2747` uses the edge

`[seed0,seed1,seed3,seed4,seed5,seed6,seed7,wall22]`

and anchor `wall23`. The last two walls come from entries 6 and 10 of the endpoint-family archive and share radicand `-4+169phi/16`. The cap meets `seed3,seed6,seed7` at `pi/3`, the anchor at `pi/5`, and the other five edge walls orthogonally. The incident rank-nine vertex diagram is `D4 + I2(5) + 3 A1`. Cap `2735` exchanges the first new `pi/3` coupling from `seed3` to `seed1`.

Both first caps individually give exact 73-wall extensions. Their mutual pairing is positive, so they cannot be retained together. The [first two-cap archive](../data/enumerations/fresh_h9_split_endpoint_walk_exact_20261013.json) and its [geometry](../data/enumerations/fresh_h9_split_endpoint_walk_exact_20261013_geometry.json) preserve that initial check.

## Joint configurations

The exact ten-cap family consists of scout IDs

`2720,2721,2722,2723,2735,2747,6746,6747,6748,6758`.

All 24 pairs between the first six and last four caps are dotted. All 21 within-part pairs have positive pairing. Thus `K6,4` is the exact compatibility graph, without a missing-finite-label qualification. Every maximal compatible subset has two caps, and all 24 were tested as actual chambers. No existing wall is removed.

For each chamber, the all-pair Coxeter gate passes, the complete connected affine-subdiagram census is empty, and each of its two new vertices is verified by a future timelike point with exactly nine zero inequalities, strict inequalities at all other walls, and a positive-definite physical incident Gram. The four inherited finite vertices are independently checked against **all ten caps at once**, so each 74-wall chamber has at least six distinct actual finite vertices. This is a certified lower bound, not a complete exact vertex count.

Every wall is a genuine facet: orthogonally projecting the common strict timelike interior onto one wall remains timelike and makes every other inequality strict, since all selected normal pairings are nonpositive.

The [family archive](../data/enumerations/fresh_h9_split_endpoint_walk_exact_20261013_family.json) saves the exact line equations, quadratic branches, target cuts and finite-vertex points. The [joint geometry archive](../data/enumerations/fresh_h9_split_endpoint_walk_exact_20261013_family_geometry.json) saves all 24 selected chambers, pair relations, complete affine censuses and their strict escapes.

## Remaining end and replay

The [checkpoint](../data/enumerations/fresh_h9_split_endpoint_walk_exact_20261013_checkpoint.json) saves a common future strict spacelike escape with norm

`(11764468488123570049 - 6234856568308221881 phi) / 16000000000000000000 > 0`.

All 82 wall inequalities are strictly satisfied. Its unperturbed proposal has support `[0,4,7,9]` in the ten-element basis above, corresponding to seed IDs `[0,4,7,12]`. This escape excludes completion by any subset of these ten caps even if pair compatibility is ignored. Further progress requires other caps or coherent replacements; the certificate is not a global obstruction for the seed or for all possible edge extensions.

The [producer-free verifier](../scripts/fresh_h9_split_endpoint_walk_exact_20261013_verify.py) passes:

- ten exact norm-two cap constructions and all retained pairings;
- exhaustive enumeration of the ten-node compatibility graph's 24 maximal cliques;
- 64,824 chamber pair checks, with 3,300 distinct pairs recomputed;
- 36,792 exact affine-frontier checks;
- 62 actual finite-vertex checks, including the four inherited vertices against all 82 inequalities;
- 288 chamber escapes and the one common escape.

The checker also rejects a forged cap coefficient and a timelike replacement for an escape. Input and output hashes anchor the literal source archives. A separate source-level review by the H20 research agent found no issue in the coefficient/covector convention, common-center signs, branch-validation logic, or actual-vertex proof; that review did not rerun the producer or claim a second arithmetic implementation.

A separate [existing-wall edge-frontier diagnostic](fresh_h9_split_edge_frontier_20261013.json) reconstructs the same normals and rechecks the four old finite vertices and the new `2747` vertex exactly. Its numerical traversal finds that fifth vertex from the old finite-vertex component. The diagnostic searches the displayed wall set only; it is not a search for the first mirror of a full reflection group.

```sh
PYTHONPATH=src:. python scripts/fresh_h9_split_endpoint_walk_exact_20261013_verify.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_split_endpoint_walk_exact_20261013_verify.py --family --verify
```

The [archive-only loader](../scripts/fresh_h9_split_endpoint_walk_exact_20261013_base.py) reconstructs the original 72 normals without a numerical search. The [small-family producer](../scripts/fresh_h9_split_endpoint_walk_exact_20261013_family.py) and [checkpoint producer](../scripts/fresh_h9_split_endpoint_walk_exact_20261013_checkpoint.py) reproduce the bounded experiment. Finite vertices and absence of affine diagrams are necessary progress toward compactness; the exact escape remains decisive here.
