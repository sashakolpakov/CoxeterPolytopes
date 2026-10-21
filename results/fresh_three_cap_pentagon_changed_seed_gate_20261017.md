# All-label mandatory-frame gate for the changed pentagon

FINAL/FROZEN. Exact generation and fresh saved replay PASS. No new chamber is constructed. The underlying exact control remains the infinite-volume 63-wall, 199-finite-vertex chamber with ten S6 open ends.

Delete walls 12, 47, 52, 64 and 65 and consider the exterior disk with boundary [7,8,44,63,30]. The [preceding scout](fresh_three_cap_pentagon_changed_seed_20261017.md) records its deletion cost and boundary geometry. The retained chord 30–44 and finite triangle [30,44,63] lie outside this disk. The gate excludes one-/two-center disk closures and three-center central-triangle closures with positive boundary-arc lengths, while retaining all other native walls. It does not exclude arbitrary three-facet disks, centers owning no old boundary edge, different old-side intersections, further deletions or deformations, or general finite-volume completion.

Every new actual S6 facet is right to H4 and has one of seven finite A2 profiles. At a prescribed finite endpoint, an edge with label m>=6 is an isolated I2(m) component. Thus the complete finite rank-at-most-five catalogue with symbolic label 6 gives all low recipes and all high-label families. The new [source](../scripts/fresh_three_cap_pentagon_changed_seed_gate_20261017.py) uses the frozen exact physical interval arithmetic and literal parent contacts; it does not use numerical roots or boundary scores for acceptance.

Three mandatory frames [30,7,8], [7,8,44] and [63,30,7] have nine low recipes and two high families, all with zero A2 profile. Both real norm branches are checked for each recipe. All eighteen low branches fail a strict retained-wall gate. Both branches in each high family fail on the entire closed interval 2cos(pi/m) in [sqrt(3),2], by an obtuse product or a strict consecutive-angle gap.

Every cyclic four-contact frame is also impossible. The five frames have 282 low recipes and 41 high families. All low recipes have a rigorously nonzero norm residual. Among the high families, twelve have negative discriminant; twenty-five have two exact interval roots, of which 46 are outside the complete high-cosine interval and four are in strict consecutive-angle gaps. Two two-variable families have norm strictly below two throughout their entire cosine square. One family has only the positive norm root c=2, the parallel limit: it cannot be a finite endpoint and is not mislabeled as a finite Coxeter root.

The remaining two-variable family is on [8,44,63,30], with contacts [2,m,2,l], m,l>=6. Root 44 is exactly orthogonal to the other three frame roots. Its norm equation is

\[
\frac{c_{44}^2}{2}+a c_{30}^2=2,
\qquad a=(G_{\rm frame}^{-1})_{30,30}>0,
\qquad c_i=2\cos(\pi/m_i).
\]

Since 3<=c30^2<4, exact bounds place c44 strictly between the label-15 and label-18 thresholds. Thus m is only 16 or 17. For m=16, the forced c30 lies strictly between the label-9 and label-10 thresholds. For m=17, it lies strictly between the label-6 and label-7 thresholds. Neither is a Coxeter angle. This is a complete finite reduction of the unbounded label family, not a sampled label cutoff.

For one or two centers, some cap covers at least three consecutive old boundary edges, hence has one of the excluded four-contact frames. For three centers forming a central triangle, positive arc lengths summing to five are 1+1+3 or 1+2+2. The first again contains an excluded four-contact frame. Each of the five rotations of the second contains one of the excluded three-contact frames. No choice of mutual angles or an additional norm branch can repair a nonexistent mandatory cap.

The [certificate](fresh_three_cap_pentagon_changed_seed_gate_20261017.json) records all exact frame Grams, elimination pivots, norm polynomials, real-root intervals and strict cosine comparisons, together with all 1+2+2 rotations. Shared arithmetic and the inherited finite classification are explicit dependencies; this is not an independent second arithmetic implementation.

The obstruction identifies a genuine change needed for this boundary: a larger disk, a center with different boundary incidence, or another retained-wall change. It gives no reason to enlarge the same label grid. No further seed is searched in this package.

Validation: generation session 4807 and fresh saved replay 45185 both terminated successfully. Syntax, all direct source/data hashes and local report links pass. The preceding numerical scout separately passed generation 72457 and fresh replay 40667.
