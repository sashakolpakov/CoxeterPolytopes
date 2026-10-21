# An exact cutter on a previously unbounded H2 bridge face

The [exact certificate](fresh_h20_e7_bridge_face_20261014_exact.json) adds one norm-two wall to the frozen [95-wall chamber](fresh_h20_e7_edge_cusp_repair_20261014.md), retaining all its walls. The resulting **96 actual facets** cut the previous strict escape and introduce an actual rank-19 edge. A new exact rational strict escape remains, so this chamber still has infinite volume. This is a geometric construction, with no arithmetic-lattice assumption.

The new wall has 18 right contacts, one π/3 contact with old wall 0, one π/6 contact with old wall 28, and 75 dotted pairs. Its right support is

```
S61 = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,17,31,61,87].
```

The support has positive-definite rational Gram of type `A10+A4+A2+2A1`. The separate [whole-face certificate](fresh_h20_bridge_faces_20261014_review.md) proves that its old codimension-18 face is the entire hyperbolic plane: every nonsupport old normal projects to a strictly future timelike vector. Thus this addition supplies a geodesic boundary on a face that previously had none. It does not complete that plane to a finite-area polygon.

## Exact construction and geometry

Let `P` be orthogonal projection off the support, and put `A0=P(N0)` and `A1=P(N28)`. Their two-dimensional Gram has negative determinant. The checker solves for the unique `p` in their span with

`p·N0 = −1`, `p·N28 = −√3`.

It projects the old strict center onto the one-dimensional complement of the support and both anchors, obtaining `z` with `z²>0`. With `K=(2−p²)/z²>0`, it sets

`r = p + √K z`.

The stored positive branch is the unique one passing every old acuteness inequality and cutting the old escape. The norm identity `r²=2`, the two anchor equalities, all 18 support equalities, and every dotted inequality are checked exactly. The field is the inherited multiquadratic field with its bridge square root, followed by the positive square root of the saved `K`; the certificate explicitly records the coordinate convention. No numerical angle tolerance is used for acceptance.

Projecting the old strict center first off the support and then onto `r⊥` gives an explicit future timelike point with exactly `S61∪{95}` active. All other 77 inequalities are strict. This proves an actual edge of the enlarged chamber. The checker also constructs a common strict future timelike center and a relative interior point of the new facet.

Exact pairwise acuteness and the common strict center imply that all 96 listed walls are actual facets, by the [acute-extension preservation lemma](acute_extension_preservation_20261014.md). The same lemma strictly preserves every old finite vertex. In particular, the at least 116 distinct finite vertices certified by the predecessor and its [repaired-cusp walk](fresh_h20_e7_repaired_cusp_walk_20261014.md) remain. All three previously bounded cusps are strictly retained; their new-wall pairings are also checked explicitly. This package does not enumerate new finite vertices or claim a complete vertex or edge census.

## Renewed infinite-volume witness

The previous saved escape has positive pairing with `r`. To find a replacement, the checker uses one of the previously certified open cusp-edge limits `Z`, at which the new wall is strict. If `U` is minus the sum of the 19 physical duals of its active positive-definite support, then `Z·U=0`, `U²>0`, and every active pairing of `U` is −1. For a suitable positive `t` and small positive rational `ε`,

`X = Z + εU + ε²tH`

is strictly feasible, future oriented, and spacelike. Decimal evaluation is used only to propose a rational rounding; every final sign and norm is rechecked exactly. The saved rational witness satisfies

`X² = 20181174579467302741 / 10^18 > 0`.

All 96 wall pairings are strictly negative. Together with the strict timelike center, this gives an open ideal-boundary region and proves infinite volume. Cutting the earlier witness has not closed the chamber globally.

## Bounded numerical search scope

The separate [scout](fresh_h20_e7_bridge_face_20261014_scout.json) tested two specified PD18 bridge supports and the finite anchor-label list `2,3,4,5,6,8,10,12,16,24,30,60`, together with parallel contacts. It found 668 numerical Coxeter proposals, of which 32 cut the old escape. The exact construction above certifies one selected proposal independently of those numerical acceptance tests.

The [independent numerical polygon diagnostic](fresh_h20_e7_bridge_polygon_20261014.md) recomputes all 120,438 same-face pairs. It finds 40,265 compatible pairs, all dotted, so the bounded pool does not provide a polygonal completion. This is a numerical diagnostic of that pool, not an all-label or all-parameter exclusion. Caps with nonzero contacts to the support, and other support choices, remain outside it. There is no claim that the present 96-wall chamber is a finite-volume polytope.

## Reproduction

The [producer and verifier](../scripts/fresh_h20_e7_bridge_face_20261014_exact.py) binds its exact predecessor data and arithmetic dependencies.

```
PYTHONPATH=src:. python scripts/fresh_h20_e7_bridge_face_20261014_exact.py --verify
```

Final generation and saved-result replay passed. The [independent new-cap review](fresh_h20_e7_bridge_face_20261014_exact_review.md) also passed generation and saved replay. Its separate sparse Fraction arithmetic and outward rational sign enclosures reconstruct all new norm, pair, point and escape data without importing the producer's arithmetic. The whole-H2 premise has its own independent arithmetic certificate linked above. Predecessor vertex counts and bounded-link censuses remain explicitly inherited.
