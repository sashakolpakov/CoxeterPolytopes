# Inherited cusp footprint of the 860-facet H22 chamber

Exactly **303** of the **1,095,990** inherited H21 cusps have bounded links in the frozen 860-facet chamber. The remaining **1,095,687** links are unbounded. Any finite-facet, finite-volume Coxeter completion retaining all 860 actual facets requires at least **2,191,374 additional facets**, hence at least **2,192,234 facets in total**.

This is an all-real geometric lower bound for retaining this particular core. It is not a nonexistence theorem, a lower bound for arbitrary H22 polytopes, or a claim that repairing all inherited cusps would suffice for finite volume. Deleting or replacing retained facets can escape its premises.

The [exact checker](../scripts/fresh_h22_inherited_cusp_footprint_20261013.py) and [certificate](fresh_h22_inherited_cusp_footprint_20261013.json) use the independently replayed [complete H21 corner census](fresh_h21_source_corner_walk_20261013.md) and the [860-facet construction](fresh_h22_multiquadratic_20261013.md). They do not enumerate a million literal rays. Instead, the checker verifies the 303 actual null projections and uses Lorentzian geometry to exclude every other inherited-cusp incidence.

The [independent source/proof review and fresh replay](fresh_h22_inherited_cusp_footprint_20261013_review.md) pass, including the finite-facet and retained-core restrictions in the lower-bound theorem.

## Exact footprint

Write the raised form as `Q ⊕ <2>` and a cap normal as `r=(u,b)`. Its inherited-cusp incidences depend only on its old projection u. The certificate reconstructs all projections using the [separate squarefree-radical arithmetic](../scripts/fresh_h22_multiquadratic_20261013_review.py), with no producer interval or projection code imported.

| Retained normals | Number | Old projection |
| --- | ---: | --- |
| Original vertical H21 walls | 210 | Old spacelike facet normals |
| Paired rational parent caps | 602 | 301 distinct future null rays |
| Added types 4 and 5, with reflected mates | 4 | Two further distinct future null rays |
| All other added types, with reflected mates | 44 | Strictly future timelike |
| Pure vertical caps | 0 | — |

For each of the 303 null directions the checker verifies its primitive rational projective representative, null norm, future orientation, all 210 old inequalities, and rank-21 full active set. It independently canonicalizes that active set under the certified 241,920-element group and matches an actual cusp in the complete census. Thus neither duplicate directions nor merely feasible nonvertex points inflate the count.

The bounded directions are distributed as follows; the other 15 cusp orbits have no bounded literal cusp in this chamber.

| Norm-two cusp type | Full orbit size | Bounded literal cusps |
| --- | ---: | ---: |
| A5⁴ | 280 | 280 |
| D4⁵ | 42 | 21 |
| D12 + E8 | 241,920 | 2 |

The root-type names come from the separate [quotient-root classification](fresh_h22_extension_feasibility_20261013.md); the footprint computation uses exact active sets and orbit indices, not the names.

At each listed cusp δ, the two cap normals have the same old projection `lambda*delta`, opposite nonzero heights, equal positive norm q, and mutual pairing `-q`. Their positive-kernel relation has weights `(1,1)`. Both are orthogonal to every old incident wall. The old affine components span the 20-dimensional old Euclidean quotient; this additional affine A1 spans the new coordinate direction. The combined rank is 21, and its positive-kernel identities leave no nonzero recession direction. The actual horospherical link is therefore bounded.

All 22 remaining positive-height cap types have strictly negative old-projection norm, verified exactly; their reflected copies have the same projection. Such projections pair strictly negatively with every future inherited null ray, so none is incident to any other inherited cusp.

## Why every missing cusp needs two new facets

Let a proposed completion retain all 860 actual facets. Its new facet normals must pair nonpositively with the 210 retained vertical normals. Thus every old projection u belongs to the closed source chamber cone C.

The source center h is a positive multiple of the sum of the 210 normals, and these normals span the ambient space. Therefore `(h,u)<0` for every nonzero u in C. The slice `(h,u)=-1` is compact. Its complete vertex census consists entirely of future timelike or null rays, so C consists of future causal vectors and zero. This conclusion is geometric and imposes no arithmetic or field restriction on a proposed new cap.

For a future inherited null ray δ, Lorentzian Cauchy gives

`(u,delta) <= 0`,

with equality for nonzero u only when `u=lambda*delta`, `lambda>0`. Hence a new facet incident to an inherited cusp has either that specific null projection or zero old projection. A nonzero null projection is incident to exactly one projective inherited cusp.

The zero-projection exception cannot occur here. A pure vertical normal `(0,b)` has positive pairing `2*b*c` with any retained cap of the same height sign c. The frozen chamber retains both height signs, so every nonzero pure vertical normal conflicts with a retained facet.

For a cap incident to δ, therefore, `r=(lambda*delta,b)` with `b!=0`. On the Euclidean quotient at `(delta,0)`, its normal points purely in the added coordinate direction. Both signs of b are necessary to bound that direction.

One can see necessity directly, without presuming a cusp-link classification for the proposed completion. If one height sign is missing, choose that unblocked sign of a sufficiently small displacement `(delta,0)+epsilon*e_last`. This vector is spacelike, with squared norm `2*epsilon²`; all incident inequalities are nonpositive. Every nonincident inequality remains strict for sufficiently small epsilon because the completion has finitely many facets. A sufficiently small perturbation toward its strict timelike interior makes all inequalities strict while preserving positive norm. The resulting strict spacelike escape proves infinite volume.

Thus each of the 1,095,687 missing inherited cusps needs at least two new facets, one of each height sign. A new facet can supply this incidence for only one inherited cusp. Consequently

`additional facets >= 2 * 1,095,687 = 2,191,374`,

and retaining the existing 860 gives the stated total lower bound. The argument allows arbitrary real cap coordinates and arbitrary Coxeter labels. It says nothing against changing the retained core.

## Replay and dependencies

A separate [source and theorem review with a fresh saved-result replay](fresh_h22_inherited_cusp_footprint_20261013_review.md) passes. It checks the projection classification, distinct cusp assignments, bounded links, and the scope of the lower bound.

The certificate binds the source data, complete corner census and its independent audit, frozen 860-facet certificate and its independent radical replay, symmetry archive, and checker code by SHA-256. The finite computation certifies the exact projection classification, all 303 distinct cusp rays and affine-A1 link completions, and the counts. The geometric argument above turns those checks into the lower-bound theorem.

```sh
PYTHONPATH=src:. python scripts/fresh_h22_inherited_cusp_footprint_20261013.py --verify
```
