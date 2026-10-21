# Escape-cutting caps and a permanent cusp obstruction

Three exact mutually dotted caps cut the renewed escape of the frozen 49-wall chamber, but their 52-wall chamber cannot be completed to finite volume by acute-compatible additions retaining all its walls. A positive affine G2 circuit is orthogonal to a retained pi/5 ridge. This obstruction is independent of the complete cone computation.

The [cap certificate](fresh_nonarithmetic_escape49_20261016_cluster.json) and [checker](../scripts/fresh_nonarithmetic_escape49_20261016_cluster.py) use the positive-definite support `S=[0,1,2,3,5,6,7,30]` and anchor33. Their nonright contact profiles are respectively `{30:3,7:6,33:6}`, `{3:3,7:5,33:6}`, and `{5:6,6:4,33:6}`; every omitted contact in `S+{33}` is right, and every other old contact is strictly dotted. The three new mutual pairs are strictly dotted. The inherited center is strict and all 52 facets are actual.

The formulas are exact algebraic recipes, not rounded normal vectors. Project the anchor33 into the Lorentz two-plane perpendicular to S, obtaining `a`, and choose its positive spacelike orthogonal direction `k`. For the prescribed contact vector `c`, put `pS=NS*GS^-1*c`, `p=pS+alpha*a`, with alpha chosen for the anchor contact. The cap is `p+sqrt((2-p²)/k²)*k`. The checker verifies the required exact orthogonal identities and certifies every unspecified strict inequality with outward rational intervals.

The [complete cone](fresh_nonarithmetic_escape49_20261016_cone.json), reconstructed by the [clipping checker](../scripts/fresh_nonarithmetic_escape49_20261016_cone.py), has **178 finite vertices, one ideal ray and 2,278 exterior rays**. All 178 old finite vertices survive. The caps cut 5, 21 and 10 old extreme directions respectively, with 22 distinct old directions in their union. A renewed strict spacelike escape is certified explicitly. Thus the chamber itself has infinite volume, regardless of the stronger completion obstruction below. This is a complete producer cone replay, not an independent second cone implementation.

## Five retained walls already prevent completion

The [small exact obstruction](fresh_nonarithmetic_escape49_20261016_cusp_obstruction.json) and [checker](../scripts/fresh_nonarithmetic_escape49_20261016_cusp_obstruction.py) verify the normalized Gram on walls `[30,33,49]`:

```
[ 2       0       -1       ]
[ 0       2       -sqrt(3) ]
[-1      -sqrt(3)  2       ].
```

Its positive kernel is `(1,sqrt(3),2)`. Hence `nu=n30+sqrt(3)*n33+2*r49` is nonzero future null, with `nu.H<0`. It is orthogonal to precisely `[0,1,2,3,5,6,30,33,49]` in this chamber and strictly negative on its other walls. The two retained normals0,1 are orthogonal to the affine circuit and have Gram `[[2,-phi],[-phi,2]]`. The certificate also constructs a strict timelike point on their actual ridge, with every other wall inequality strict.

For any proposed acute-compatible completion, every added normal has nonpositive pairing with nu, because nu is a positive sum of retained outward normals. Using the proposed completion's own strict center shows that nu remains in its ideal boundary; the original center need not survive. Finite volume would force a bounded Euclidean Coxeter link there. Projection of a strict link point preserves the old0,1 facets and their pi/5 ridge. Their reflection product acts on the eight-dimensional Euclidean link with trace `5+phi`, which is irrational, whereas a cocompact Euclidean reflection group preserves a full translation lattice and has integral linear traces. This is the previously established [fixed-core cusp argument](fresh_h9_fixed_core_cusp_obstruction_20261013.md), also explained in the [forbidden ridge-label lemma](fresh_h20_weighted_binary_20261014_dotted_completion_review.md).

Thus retaining just walls `[0,1,30,33,49]` rules out finite volume after any finite collection of acute-compatible additions. Deleting or replacing an anchor is outside the theorem. The untouched49 chamber is not excluded by this witness. The next construction starts again from that frozen49 chamber and filters affine circuits of every rank, rather than appending to this obstructed52 chamber.

## Discovery and negative screens

The [bounded old-bank rescore](fresh_nonarithmetic_escape49_20261016_rescore.json) checks 225,226 numerical proposals from nine saved banks against all49 walls and the renewed escape; none is a compatible cutter. This is a bounded numerical result only.

The [polar screen](fresh_nonarithmetic_escape49_20261016_polar.json) examines all 2,199 exterior-ray directions. Its separate [exact all-label gate](fresh_nonarithmetic_escape49_20261016_polar_gate.json), produced by [this checker](../scripts/fresh_nonarithmetic_escape49_20261016_polar_gate.py), proves 2,194 directions strictly noncutting and separates the other five from consecutive allowed cosine values. There are no unresolved polar cutters. The conclusion concerns these complete-cone polar directions, not arbitrary normals.

The [targeted numerical scout](fresh_nonarithmetic_escape49_20261016_target_scout.json) uses the four positive-definite eight-root subsets of the escape's active set `[0,1,2,3,5,6,7,30,33]`. It allows at most two nonright support contacts from labels3,4,5,6,10, and one outside anchor from labels2,3,4,5,6,7,8,10,12. Both isolated conic branches are tested. There are 2,964 support profiles, 256,690 real isolated branches, and 1,846 numerical cutter proposals. Numerical proposals alone carry no exact acceptance or all-label exclusion claim. Only the three explicit algebraic recipes above are promoted in this package.

All seven scripts (`rescore`, `polar`, `polar_gate`, `target_scout`, `cluster`, `cone`, `cusp_obstruction`) passed generation and saved `--verify` replay. The cone and cusp obstruction are distinct certificates: complete clipping measures this chamber, while the five-wall theorem excludes every coherent finite-volume extension retaining its anchors. Independent review is recorded separately by the parent checkpoint.
