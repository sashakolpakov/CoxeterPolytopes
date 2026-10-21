# Polar extensions, a second arithmetic-host obstruction, and compact contact search

No new finite-volume Coxeter polytope or dimension is established. This phase gives a verified 25-facet infinite-volume H9 chamber, excludes a different arithmetic H9 lattice completely, and advances the compact incidence search through a new unbounded-label obstruction. Every previous frozen certificate is preserved.

## A new arithmetic host is now ruled out

The lattice `L55 = E8 ⊥ [[−2φ,−1],[−1,4]]` has determinant of absolute field norm 55 and is not similar over `Q(√5)` to the previously excluded discriminant-11 host. Its complete primitive reflection norms, up to unit rescaling, are `2`, `2(2+φ)`, `2(4−φ)` and `2(7+φ)`. The zero sector is E8, with 240 roots and eight simple walls; its binary coordinate vector is not a reflector.

The [height-20 pilot](fresh_h9_beta55_pilot_20261014.md) gives twelve actual walls and is independently reproduced with separate arithmetic. The [height-100 search](fresh_h9_beta55_vinberg_h100_20261014.md) gives 44 walls at 132 control levels. Its complete cone has 129 finite vertices and no cusps. The [independent enumeration](fresh_h9_beta55_vinberg_h100_20261014_review.md) visits 61,144 nodes in reversed coordinate order and reproduces every accepted root, compared with the producer's 90,454 nodes.

An exact integral isometry maps complete finite corner 3 to corner 42, with positive vertex scale `2+2φ`, determinant `−1` and trace `5+8φ>10`. The [symmetry certificate](fresh_h9_beta55_corner_symmetries_20261014.md) and [independent review](fresh_h9_beta55_corner_symmetries_20261014_review.md) verify the matrix, integral inverse, outward-normal images and both basepoint decompositions. The finite-corner preservation lemma shows that no later simple wall can cut or touch either corner. Thus the map stabilizes the full reflection chamber; its infinite order proves infinite volume.

This excludes every finite-volume Coxeter polytope whose wall reflections preserve this literal lattice. It does not exclude merely commensurable lattices, other arithmetic hosts or arbitrary geometric walls. No higher-height run is needed here.

## Six polar caps extend the nonarithmetic chamber

The [first complete polar census](fresh_h9_quartet_polar_20261014.md) tests all 285 exterior vertices of the nineteen-wall chamber. Four compatible polar walls can be added together. The [first extension](fresh_h9_quartet_polar_extension_20261014.md) has 23 actual facets, 69 finite vertices and 213 finite edges. Its 413 exterior rays expose two further compatible polar walls.

The [second extension](fresh_h9_quartet_polar_second_20261014.md) adds those two walls and reaches:

| quantity | value |
|---|---:|
| actual facets | 25 |
| finite vertices | 81 |
| ideal vertices | 0 |
| finite edges | 247 |
| open finite-corner directions | 235 |
| hyperideal extreme rays | 513 |

Every preceding finite vertex survives. All 300 wall pairs have exact Coxeter relations. The old strict spacelike escape survives every added wall, so the resulting chamber has infinite volume.

The [independent review](fresh_h9_quartet_polar_extension_20261014_review.md) uses a separate algebraic-number implementation and independently coded complete-face cone clipping. It reproduces all three full cone censuses, all 1,211 polar-candidate tests, the final facet witnesses and finite-edge graph, and the strict escape. The complete finite-angle alphabet in this quartic field is `2,3,4,5,6,10`; it is derived from the field's mixed signature and unique quadratic subfield, without a finite-label cutoff.

Four further polar candidates, indices 62,95,96,128 in the 25-wall cone, are mutually compatible but all miss the escape. They are saved as the explicit next frontier and have not been added. Repeating these additions alone has not resolved the volume obstruction; nonpolar walls and changes to retained walls remain open.

## Compact-H9 contact search

The [new dense-third obstruction](fresh_h9_star_third34_20261014.md) covers every third own-label pair at least 3 and 4, in either order, with two single-π/3 anchors. It allows shared or opposite anchor orientations, arbitrary noncore contacts, block labels at least four and anchor-block label at least seven. Six continuous polynomial identities have 66 strictly positive coefficients. The [independent review](fresh_h9_star_third34_20261014_review.md) checks all 30 literal face guards, four budgets, six polynomials and eight orientations.

Combining this with the preceding `(3,3)` theorems covers every dense third profile. Twelve previously reviewed literal core placements give 48 new predicates, raising the bound search to 301 predicates. [Batch 012](fresh_h9_star_third34_20261014_search_batch_012.json) adds 5,404 visits and reaches 26,241 cumulative visits, with nine pending frames. The next assignment has third own profile `(2,5)`, outside the dense-third theorem. This is a necessary contact assignment, not a metric realization; the full incidence remains open.

The following [single-third theorem](fresh_h9_star_single_third_20261014.md) covers the changed own profile `(2,k)` for every finite `k≥5`, with owner internal label at least four, three other internals at least seven, and shared single-π/3 anchors. Its [independent review](fresh_h9_star_single_third_20261014_review.md) reconstructs 48 literal guards, three budgets and eleven continuous polynomials with 78 positive coefficients. Three projected roots would have nonzero Gram determinant in a two-dimensional complement, giving the contradiction.

Twelve literal placements add another 48 predicates. The search now uses **349 reviewed predicates**. [Batch 013](fresh_h9_star_single_third_20261014_search_batch_013.json) adds 710 visits, reaching **26,951 cumulative visits**, with nine pending frames. Its next state retains owner label four and third profile `(2,5)`, but all three high blocks now have **opposite** single-π/3 anchor profiles. The shared-anchor theorem does not cover it. The complete incidence and its unrestricted metric realization remain open; the cursor and siblings are preserved.

## Two scoped H20 construction obstructions

The [direct B2 coupling experiment](fresh_h20_b2_h18_coupling_20261014.md) starts from the complete forty-wall H18 source and produces a coherent 42-facet H20 chamber. More strongly, one exact strict spacelike escape satisfies every wall in the 74-normal union of all allowable inducing profiles. Hence every direct B2 coupling in this fixed-source family has infinite volume, regardless of the profile assignment. Extra noninducing walls and changed sources remain outside the theorem.

The [independent review](fresh_h20_b2_h18_coupling_20261014_review.md) reconstructs the source coordinates, all profile norms and escape signs, all 861 example pairings, the 42 facet witnesses and the actual source face with separate exact arithmetic. It does not repeat the inherited H18 finite-volume source proof.

The separate [direct A2 test](fresh_h20_a2_h18_coupling_20261014.md) covers all seven inducing profiles, including single π/5 contacts. The connected source π/3 component, together with an orthogonal pair in it, forces 24 source facets to stay ordinary. A common strict escape satisfies every wall in the remaining 138-normal relaxed bank. The [independent review](fresh_h20_a2_h18_coupling_20261014_review.md) expands that entire bank in a separate multiquadratic field and checks all norms and strict signs. Thus the direct A2 family also has infinite volume. Its additional-wall and changed-source limitations are the same as for the B2 theorem.

## Verification and next directions

The [archive audit](polar_extensions_and_discriminant55_obstruction_2026_10_15_audit.json) binds this phase's source hashes, reports and local links and checks the frozen [preceding archive](full_norm_vinberg_and_star_contacts_2026_10_14_audit.json). Archive integrity and mathematical proof have separate scopes; each report states its own completeness and review limits.

Do not restart reflection searches in either excluded literal arithmetic host. Preserve the compact search's saved stack and siblings. The geometric search needs a wall construction that addresses the surviving escape, rather than a facet-count increase alone. Nonarithmetic H20 remains the genuinely missing finite-volume dimension target; compact H9 remains the target for raising the compact dimension record. The [literature restrictions](arithmetic_dimension_bounds_and_effective_high_dimension_search.md) and [Vinberg/Allcock distinction](vinberg_allcock_search_strategy_20261013.md) remain in force.
