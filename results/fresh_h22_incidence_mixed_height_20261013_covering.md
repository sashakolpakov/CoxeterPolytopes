# Covering the old chamber: an exact reformulation and a missing cusp

The mixed 1,336-facet construction has an explicit **untouched old cusp**, in addition to its previously saved remote escape. An independent stdlib checker verifies that no new wall is incident there and that the link has an unbounded vertical direction. It also supplies a simpler strict escape of square **8,386** against the entire 2,296-wall relaxed union.

[Checker](../scripts/fresh_h22_incidence_mixed_height_20261013_covering.py) · [certificate](../data/enumerations/fresh_h22_incidence_mixed_height_20261013_covering.json) · [construction report](fresh_h22_incidence_mixed_height_20261013.md).

A separate [source/mathematical review and fresh replays](fresh_h22_incidence_covering_20261013_review.md) passed for the covering criterion and the literal missed-cusp certificates.

The obstruction from timelike caps missing old cusps was already used in the repository's [earlier cusp-orbit construction](fresh_h22_cusp_orbits_20261013.md), which repaired 15 sampled symmetry orbits totaling 1,065,708 old rays. That census was explicitly nonexhaustive. The contribution here is the literal missed-ray certificate for the new branches and the covering reformulation; the cusp obligation itself is not new.

## The covering criterion

Retain the 210 vertical walls of the finite-volume old H21 chamber P. Assume a finite family of caps with a strict common interior point, positive upper and lower height bounds, and old projections in the future causal cone. For an old point x normalized by x²=−1, an upper cap `(u,a)`, a>0, imposes

`z ≤ H_r(x)=−(u,x)/a`.

For lower caps use the same function with |a| and the bound `z≥−H_r(x)`. If there are several caps on a side, take the minimum of their functions. The retained old cone contains only future causal vectors, since its original chamber has finite volume.

Under these hypotheses, the H22 chamber has finite volume **if and only if** both sides satisfy `min_r H_r(x)≤1` for every x in the entire old chamber P. If a minimum exceeds one, choose z between one and that minimum; after a small perturbation toward a strict interior point this is a strict spacelike escape. Conversely, both inequalities place the full polyhedral cone inside the future light cone. A projective slice is then a finite polytope with only finite or ideal vertices, which has finite hyperbolic volume.

The word *entire* is essential. Coverage of a compact truncation or a numerical sample is insufficient. If all old timelike points are covered, the required old-null limits follow by continuity and homogeneity: for an old null ray δ, each side's minimum must be zero. Thus no separate cusp assumption is needed for the literal all-point equivalence. In a practical computation that covers only a truncated region, the omitted cusp neighborhoods need separate proofs.

If u is future timelike and a²>q=r², write `c=u/√(a²−q)`. The inequality H_r(x)≤1 is exactly the closed hyperbolic ball

`d(c,x) ≤ arcosh(a/√(a²−q))`.

For a future null u it is a closed horoball. These are intersections with P; their centers need not be interior points of P for the equivalence itself. In this branch they are feasible old points.

For the mixed branch, put `ρ=(√6+√3)/2`. Each side has 322 horoballs `−(δ,x)≤√2/ρ`, 240 balls of radius `arcosh√3`, and one ball centered at `u/√20` with radius `arcosh√(11/10)`. The two height signs have the same spatial covering sets. This reformulation is a useful search objective, not a claim of a new theorem or of covering already achieved.

## A missed cusp that finite balls cannot repair

The certificate checks the future null vector

`δ=(20,−47,−32,−63,−93,−76,−58,−40,−20,80,52,106,160,135,111,87,64,−20,−40,−29,−16,42)`.

It satisfies all 210 old inequalities. Its 24 incident old normals have physical span rank 21 and Gram rank 20. All added caps are strictly negative at `(δ,0)`:

* Maximum pairing with the 322 null projections: −2.
* Maximum pairing with twice the 720 timelike projections: −4.
* Pairing with the new projection u: −111.

Therefore the added-coordinate vector is a positive-norm recession direction of the link. The link has quotient rank 20 in the required 21-dimensional horosphere and leaves a full vertical direction free. More directly, `(δ,ε)` is spacelike for every nonzero ε, preserves all old inequalities, and remains strict on all caps for sufficiently small ε.

The explicit vector `(100δ+h,100)` makes every inequality strict. Its square is 8,386 and its center pairing is −814. All signs are independently verified using integer arithmetic and squared radical bounds.

This also excludes repairing this particular cusp using any **finite** collection of additional caps whose old projections are future timelike: all remain strict at δ. Null projections not proportional to δ remain strict as well. Geometrically, finitely many finite-radius balls cannot cover the end of this cusp.

There is a stronger necessary condition when every old wall is retained as an acute-compatible facet. Any added normal `(u,a)` has `(u,w)≤0` for all old walls, so u lies in the old chamber cone. That cone is future causal. If u≠0, incidence at `(δ,0)` therefore forces u to be null and positively proportional to δ. To bound both vertical directions at δ, one needs such incident caps on both height sides. This is necessary, not sufficient for compatibility or a bounded link.

The exceptional old projection u=0 gives a pure vertical wall, which can bound one height side globally. Both opposite pure vertical halfspaces would collapse the interior. In the current mixed and π7 branches, even a single pure vertical addition is incompatible with a retained pair of opposite-height caps: its pairing with one member is positive. Thus that exception does not repair these fixed cores. A useful continuation needs new null projections based at the missing δ, or replacement of retained walls.

## The same missed ray in the π/7 branch

The [separate Fraction interval checker](../scripts/fresh_h22_incidence_prism_20261013_uncapped.py) and [output](fresh_h22_incidence_prism_20261013_uncapped.json) verify strict negativity at this δ for every one of the 23 new π/7 caps. Reflection in the added-coordinate plane leaves those pairings unchanged. Consequently the same old cusp remains untouched in both the 835-wall branch and its 837-wall extension by two reflected mates, even though that extension correctly repairs two other specified cusps.

There is also a direct strict escape for all these walls: `(100δ+h,40)` in the `diag(Q,2)` coordinate convention has square **1,586** and center pairing −814. The checker proves strictness for all 322 unscaled null projections, a superset of the 301 parent projections, and for **both height signs of all 23** π/7 caps. No statement is made here about unrelated later caps.

## Explicit next proposal, not inserted

The paired normals `(3δ,±1)` in `diag(Q,2)` have norm two and mutual pairing −2. Their physical-coordinate form is `(3δ,±√2)`. They are natural candidates for the missing vertical affine A1 at this cusp. Against the mixed branch, old nonzero pairings are integers at most −3, and the worst √3-cap cross pairing is `−6+√6<−2`; the scaled-null and √22-cap pairings are also dotted by the saved projection maxima. Their positive-height member cuts the π7 strict escape by exactly `−24+80=56`.

These are **proposed next normals**, not part of any frozen chamber in this report. Their relations to all 23 π7 caps and the later 839-wall extension have not been audited here, and no resulting volume claim is made. Rather than continuing one-cusp-at-a-time growth, a more useful next step is a complete old-chamber cusp-orbit census or walk, followed by a coherent covering test that includes every required old cusp direction.

```sh
python scripts/fresh_h22_incidence_mixed_height_20261013_covering.py --verify
python scripts/fresh_h22_incidence_prism_20261013_uncapped.py --verify
```
