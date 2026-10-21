# A remote finite-vertex component, followed by a complete symmetry-family obstruction

No finite-volume H9 example was found. A new six-cap construction retains the previous96 walls, cuts their saved remote escape and creates three actual finite vertices. The resulting102-facet chamber has no positive affine circuit, but another exact strict spacelike escape remains. The larger symmetry proposal is now closed: a single exact point survives all102 walls and all13,824 images of the six source caps under the complete weighted16-seed symmetry group. Choosing different compatible subsets or deleting old walls within that finite dictionary cannot give finite volume.

## The exact102-wall construction

The [previous96-wall checkpoint](fresh_h9_face_bridge_20261013.md) contains a compact H4 prism. Its saved escape survived the earlier17,872 finite-anchor proposals. Instead of extending that same finite component, this construction transports six existing caps by the seed isometry

```text
U: coordinates1 and2 exchange, coordinates5 and6 exchange.
```

The sources are `compact0, compact1, compact2, compact3, facePrism, faceGap`. The [recipe](../data/enumerations/fresh_h9_prism_completion_20261013_caps.json) binds the preceding source by hash. The [producer-free exact checker](../scripts/fresh_h9_prism_completion_20261013_verify.py) rechecks the seed Gram isometry and source construction, then certifies the transported normals. Their norms remain two. Their finite pair equalities follow from the isometry: `wall17` and `nextH4` are fixed because each is perpendicular to seeds1,2,5,6. Every other unprescribed pairing is rigorously strictly dotted, checked using outward-rounded integer intervals.

The [exact certificate](../data/enumerations/fresh_h9_prism_completion_20261013_verify.json) records:

- All591 new pairs:38 right,4 with angle pi/3,549 dotted.
- All992 possible new positive-semidefinite subsets: every one is positive definite. There are no new affine circuits; the inherited absence is replayed.
- Three actual finite vertices, with diagrams `A2+7A1`, `A3+6A1`, and `A3+A2+4A1`.
- Their complete finite component: two finite edges and23 open rays. Every possible ninth wall is checked, including exclusion of ideal endpoints.
- Strict common timelike interior and nonpositive pairings, which certify that all102 inequalities are true facets by orthogonal projection onto each wall.

The transported `faceGap` cuts the old strict escape. This is a separate finite component: each new wall has at most five finite old neighbors, so it cannot close an inherited rank-eight spherical edge. The construction therefore does not join the old prism component to the new one.

The whole chamber remains infinite-volume. Its next strict future spacelike escape has square

```text
(4764512044018601985931865
 -2484892311721629081203273 phi) / 16000000000000000000000000 > 0.
```

The numerical proposal has support `[0,1,3,4,5,8,9]`; its active old walls are saved in the recipe. Its rational perturbation is checked strictly against all102 walls. A numerical vertex census supplied this proposal; it is not used as a compactness proof.

## The full symmetry direction is now excluded for these six sources

The coordinate swap above belongs to a larger group. Write the16 seed normals as `U0,...,U7,V0,...,V7`, with two blocks of four indices. In the10-coordinate basis `[U0,...,U7,V0,V4]`, all seed vectors are integral:

```text
Vj = V0 + U0 - Uj       for j=1,2,3,
Vj = V4 + U4 - Uj       for j=5,6,7.
```

The two layers are orthogonal8-cliques. Their cross Gram has three distinct weights: matched indices, different indices in the same block, and indices in different blocks. The unique matching forces a common index permutation on both layers; the other weights recover the4+4 partition. Thus the full weighted Gram automorphism group consists of independent permutations inside the two blocks, a block swap, and a layer swap. Its order is `24² * 2 * 2 = 2304`.

The [orbit checker](../scripts/fresh_h9_prism_completion_20261013_orbit_obstruction.py) verifies this pattern exactly, enumerates every permutation, checks its action on all16 integer seed vectors, and proves Gram invariance. It then checks every one of the `2304 * 6 = 13824` literal source-cap images. No compatibility filtering or numerical deduplication occurs in this proof.

The [common-escape certificate](../data/enumerations/fresh_h9_prism_completion_20261013_orbit_obstruction.json) and its [rational point recipe](../data/enumerations/fresh_h9_prism_completion_20261013_orbit_obstruction_point.json) give square

```text
(2578425889049542695642049
 -1361826337575782759254129 phi) / 16000000000000000000000000
 ≈ 0.02343403675 > 0.
```

Every102 old inequality and all13,824 image inequalities are strict at this future spacelike point. There is also an explicit common strict timelike center: the sum of the16 seed normals, with integer coordinates `(4,0,0,0,4,0,0,0,4,4)` and square `-64-128phi`. Every group element fixes this center. Its strict six-source inequalities consequently transfer to every image, and its102 old inequalities are checked directly. The escape's pairing with this center is `-412311744841/250000000000`.

The segment from that strict timelike center to the strict future spacelike point crosses the light cone while every inequality remains strict. An open portion of the sphere at infinity is therefore exposed, proving infinite volume. Deleting inequalities only enlarges this intersection. Hence **no selection from the102 old walls and these13,824 literal images can yield a finite-volume chamber**, even with arbitrary old-wall deletions or replacements using this bank.

This is a finite-family theorem. It says nothing about different source normals, a continuous deformation of these sources, other isometries, or unrestricted H9 completion.

## Further bounded tests and limits

The following are numerical diagnostics, separately scoped from the exact results above:

- The [new-frontier scout](../data/enumerations/fresh_h9_prism_completion_20261013_frontier_scout.json) tests the23 new open edges with at most two nonright edge relations and one finite old anchor. Labels are2,3,4,5,6,8,10,12. No cutter of the102-wall escape retained the16 seed walls with at most two optional conflicts.
- The [continuous endpoint scout](../data/enumerations/fresh_h9_prism_completion_20261013_continuous_scout.json) examines both full norm-two hyperbola branches for1,657 spherical profiles using labels3 and5. Intersecting the target-cut and center conditions with every other pairing strictly dotted produced no feasible interval. Floating interval endpoints do not constitute an exact universal obstruction.
- The [full-group image scout](../data/enumerations/fresh_h9_prism_completion_20261013_full_orbit_scout.json) found replacement candidates obstructed only by `wall22`, so optional-wall replacement was considered explicitly. The exact common-escape theorem includes those candidates and all omitted images as well.
- The [whole-orbit numerical volume gate](../data/enumerations/fresh_h9_prism_completion_20261013_full_orbit_geometry.json) proposed the common escape, using all images before compatibility screening. Only the subsequent exact replay supports the infinite-volume conclusion.

A future constructive step needs normals outside the six-source symmetry family. Enlarging its compatible subset, completing its orbit, or removing old blockers cannot finish the volume. The102-wall checkpoint remains an exactly verified intermediate chamber, not a new finite-volume Coxeter polytope.

## Replay

```sh
python scripts/fresh_h9_prism_completion_20261013_verify.py --verify
python scripts/fresh_h9_prism_completion_20261013_orbit_obstruction.py --verify
python scripts/fresh_h9_prism_completion_20261013_regressions.py
```

The [rejection tests](fresh_h9_prism_completion_20261013_regressions.json) reject a forged source hash and the already-cut96-wall target mislabeled as the remaining escape. The checkers do not import numerical scouts or LP solvers. Exact algebraic identities prove finite pair equalities; strict signs use the inherited rigorous integer interval arithmetic. The [independent source and proof review](fresh_h9_prism_completion_20261013_source_review.md) found no issue; it is distinct from these full replays.
