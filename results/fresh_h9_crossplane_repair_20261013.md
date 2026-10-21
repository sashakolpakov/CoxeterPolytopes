# H9: metric exchanges and an orthogonal pair across two remote gaps

The new exact checkpoint has **90 actual Coxeter facets**, retaining the previous 88 walls and adding an orthogonal pair. Both new caps cut both saved symmetric gaps. It has no positive affine subdiagram, but an exact strict future spacelike escape proves **infinite volume**. No compact or finite-volume example is claimed.

This is a different construction from the nine previously excluded positive normal planes: finite-angle equations are exchanged among old walls, and a reflection equation couples two caps in different directions. The pair replaces the separate 89th-wall experiment; it does not retain that extra wall.

## Exact construction

The source is the [88-wall bridge](fresh_h9_remote_next_20261013.md). All normals have square two. Put `D = seed6 − seed7`; its exact square is four. Prescribe the following pairings for the first new normal `r`:

| Old wall | Dihedral label | Pairing |
|---|---:|---:|
| seed1, seed2, wall17 | 2 | 0 |
| seed4, seed6 | 3 | −1 |
| seed3, seed8, nextH4 | 5 | −φ |

Also prescribe `(r,D)=2`. These nine independent linear conditions leave one coordinate. If `A` is the ordered list of these eight old normals followed by `D`, `M=Gram(A)`, and `p` is the prescribed pairing vector, set

`u = Σ (M⁻¹p)ᵢ Aᵢ`,

`v = c − proj_span(A)(c)`,

`r = u − sqrt((2−u²)/v²) v`,

where `c` is the existing strict timelike center. The verifier certifies that the inverse pivots are nonzero, `v²>0`, and `2−u²>0`. Orthogonal projection gives the prescribed pairings and square two identically; equality is not inferred from a small numerical residual.

Define `s = r−D`, the reflection of `r` in `D⊥`. Then `s²=2` and `(r,s)=0` exactly. The finite relations of `s` are the same list with seed6 replaced by seed7. The unchanged listed walls are exactly perpendicular to `D`; this includes an exact check of the equal dotted seed8 pairings.

Each cap has eight finite relations and 80 strictly dotted relations to the old walls. All 160 dotted inequalities, both strict center signs, and both caps' cuts of both old strict spacelike witnesses are certified with outward rational intervals. The old walls' strict inequalities at each target are checked again. The complete resulting Coxeter pair table has 4,005 entries.

All 90 inequalities are genuine facets. Project the strict timelike center onto any wall. Its square remains negative, and all other inequalities stay strict because every wall pairing is nonpositive. Thus the projected point lies in the relative interior of that facet.

The [exact recipe](../data/enumerations/fresh_h9_crossplane_repair_20261013_pair_caps.json) and [producer-free replay](../data/enumerations/fresh_h9_crossplane_repair_20261013_pair_verify.json) include the complete formulas, inequalities and escape.

## Compactness and volume gates

A positive-semidefinite subdiagram containing either new cap can use only its finite old neighbors. Enumerating the two single-cap cases and their common-neighbor case requires 640 subsets. Exact elimination over `Q(φ)` finds 138 positive-definite subsets, 108 indefinite subsets and 394 subsets excluded by a dotted pair. No positive affine subset and no rank-nine positive-definite subset occur. The old 88-wall absence of affine circuits is inherited from its archived census.

Consequently neither new wall creates a finite vertex. Absence of affine circuits does **not** prove that a compact completion exists.

The saved escape is a rational dual vector `y`; the physical point is `−G⁻¹y`. Its square is

`(4451686641999449761613929 − 2281192286026352419181229 φ) / 16000000000000000000000000 > 0`,

approximately `0.0475400`. Its pairing with the future center is negative, and all 90 dual wall inequalities are strictly positive. This gives an open spacelike patch in the projective cone and proves infinite volume. The numerical point used to propose it has support `[1,3,4,7,8,9]`; the exact saved point is deliberately perturbed to make every inequality strict.

As a numerical diagnostic only, full affine Qhull enumeration lowered the largest vertex square from approximately `0.0513527` to `0.0475425`. It found 65,406 projective vertices, of which 65,384 were spacelike. These floating counts are not used in the exact volume proof.

## Search and discarded alternatives

The bounded scouts are preserved for reproducibility, not as family-wide impossibility theorems:

* [Rank-seven two-anchor scout](../data/enumerations/fresh_h9_crossplane_repair_20261013_scout.json): all 36 relevant seven-wall subsets, stated labels `{2,3,4,5,6,8,10,12}`, and 11,079,414 numerical norm branches; no accepted cap.
* [Full nine-wall metric scout](../data/enumerations/fresh_h9_crossplane_repair_20261013_face_scout.json): 1,795,500 profiles with three through six nonright entries; no accepted cap.
* Adaptive equation exchange first found a profile whose new wall joined three mutually orthogonal seeds by labels `3,4,3`. Its Schur complement is `2−(1+2+1)/2=0`, forcing a positive affine circuit. It was rejected before exact cap construction.
* The revised [compact-aware exchange](../data/enumerations/fresh_h9_crossplane_repair_20261013_compact_exchange_scout.json) produced an [exact 89-wall alternative](../data/enumerations/fresh_h9_crossplane_repair_20261013_verify.json). Its 512-subset census has no affine circuit, but it leaves a symmetry-related strict escape. Its reflected mate has a non-Coxeter mutual angle, motivating the exact orthogonal coupling above.
* The [coupled scout](../data/enumerations/fresh_h9_crossplane_repair_20261013_pair2_scout.json) proposed 13 records. After deduplication and possible third reflections, [eight pair/triple configurations](../data/enumerations/fresh_h9_crossplane_repair_20261013_geometry_scout.json) received numerical whole-cone checks. All remained infinite in those diagnostics. Only the reported pair was exactified; no exact theorem is claimed for the other proposals.

These computations do not exclude other caps, coherent replacements, or compact completions of the original seed. Further progress needs to reach the surviving gap and create new finite incidences; the present pair repairs two specific ends but creates no new finite corner.

## Replay

```sh
python scripts/fresh_h9_crossplane_repair_20261013_verify.py --verify
python scripts/fresh_h9_crossplane_repair_20261013_pair_verify.py --verify
python scripts/fresh_h9_crossplane_repair_20261013_regressions.py
```

The first two use no numerical search or SciPy. Exact norm/equality proofs use the construction identities; strict signs use integer outward intervals at 80 decimal places. [Regression checks](fresh_h9_crossplane_repair_20261013_regressions.json) reject a timelike replacement escape, an old escape cut by the new walls, and the wrong norm branch. This arithmetic is shared with the previously audited source chamber; it is not presented as a second independent field implementation.

A [separate source and algebra review](fresh_h9_crossplane_repair_20261013_pair_review.md) also passes. Its 160-digit execution nests all 465 changed interval enclosures inside the saved 80-digit enclosures, with exact metadata equality.
