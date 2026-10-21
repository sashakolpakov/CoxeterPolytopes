# Fixed-profile endpoint belts at the E8 partner edges

The exact theorem in this package concerns a specified endpoint family, not all possible closures of the retained H3 face. At the three old edges 159/160, 159/161 and 160/161, take any of the respectively 57, 45 and 57 fixed finite endpoint profiles. Require every unprescribed old162 contact except wall158 to be strictly dotted; allow **every** Coxeter contact with158, including parallel and dotted contacts. Every Coxeter-compatible pair of such caps on different old edges is necessarily dotted. Consequently these caps cannot supply the finite adjacencies required by the proposed three-cap belt.

The elementary geometric argument below strengthens this to a **finite-batch no-go within precisely this family**: no finite coherent Coxeter extension of the168-wall chamber, using only these caps, has finite volume. The three unbounded singleton151 endpoint families, additional finite contacts to other old walls, and different intermediate facets remain open. The arithmetic exclusion already holds using the old162 subset, so the six later cusp walls cannot invalidate it. No new wall or finite-volume construction is claimed here.

The underlying exact chamber is the [162-wall E8 reconstruction](fresh_h20_e8_exact_cluster_20261017.md), with the [six-wall cusp completion](fresh_h20_e8_f4_cusps_20261017.md). The complete fixed-profile input is the frozen [endpoint census](fresh_h20_e8_partner_endpoints_20261017.md).

## Exact all-label separation

The [exact source](../scripts/fresh_h20_e8_endpoint_belt_20261017_exact_gate.py) and [saved certificate](fresh_h20_e8_endpoint_belt_20261017_exact_gate.json) use the literal 21-root core Gram and rational outward intervals from the frozen exact162 engine. Generation and fresh saved replay have passed, and this package is final/frozen. The root agent's separate source/proof read found no gap in the branch, conic, tail or cosine arguments. That read did not duplicate the arithmetic, and this remains a primary exact certificate.

Let W be one old edge's rank19 spherical support. For a prescribed norm-two cap, its pairings with W fix the projection p. Put h=2-p²>0. Let

\[
 u=n_{158}-\operatorname{proj}_W n_{158},\qquad
 \sigma=u^2>0,\qquad e=u/\sqrt\sigma.
\]

The old corner W+158 is finite. Let V be its future unit timelike vertex vector. Write x=−(r,n158) and o=(p,n158). Every possible cap with these pairings is

\[
 r=p-\frac{x+o}{\sqrt\sigma}e
   \mathbin{\pm}\sqrt{\frac{(x+o)^2}{\sigma}-h}\,V.
\]

Only the plus branch can be coherent with the retained walls. Indeed the projection onto the whole spherical corner W+158 has nonpositive coefficients in its simple roots, since all prescribed pairings are nonpositive and the inverse spherical Gram is entrywise nonnegative. It therefore pairs positively with any common strict interior point. The minus branch adds another positive term, and a zero radical cannot repair that contradiction. This argument allows the common interior point to move; it imposes no fixed-center restriction.

On the plus branch set

\[
 a=\frac{\sqrt h}{2}(V-e),\quad
 b=-\frac{\sqrt h}{2}(V+e),\quad
 r(k)=p+ak+b/k,\quad k>1.
\]

Then x=\(\sqrt{\sigma h}(k+k^{-1})/2-o\), strictly increasing with k. The finite values at labels2 through6 are 0,1,√2,φ,√3. The script checks every such value against all required old162 dotted gates. The accepted counts are258,221,267, respectively.

All remaining Coxeter possibilities have x≥9/5: the exact interval check proves2cos(π/7)>9/5, and this also includes x=2 and all x>2. For each profile let k0 correspond to x=9/5 and substitute k=k0(1+t), t≥0. Multiplication by1+t gives a quadratic vector polynomial with coefficients

\[
 p+ak_0+b/k_0,\qquad p+2ak_0,\qquad ak_0.
\]

For two tail profiles, every coefficient of
\((r\cdot s+2)(1+t)(1+v)\) is strictly negative. For a discrete profile against a tail, the analogous three coefficients are strictly negative. Thus all such pairs are dotted. The remaining discrete/discrete pairs are either dotted or lie in certified strict gaps between consecutive finite Coxeter values. In particular the48 gaps occur between edge families159/160 and160/161. They are not declared non-Coxeter by numerical tolerance.

The certificate covers8,379 tail/tail profile pairs,38,913 discrete/tail pairs,39,831 tail/discrete pairs and184,911 discrete/discrete pairs. The largest coefficient upper bounds in the three continuous cases are approximately−11.9541,−2.7244 and−0.18422, respectively; the saved rational bounds, rather than these decimal summaries, certify strictness. There are49 rejected discrete roots and746 retained ones.

The cosine comparisons use rational outward intervals: Machin's formula for π with alternating arctangent remainders, followed by the cosine series with an explicit remainder. A floating-point estimate only chooses which consecutive labels to try; both strict bounds are subsequently proved by rational interval arithmetic. No label cutoff is used in the exact theorem.

## Finite-batch consequence within the stated family

Here a coherent extension means finitely many distinct spacelike wall normals with all pairings nonpositive and a common future timelike strict interior point. Each added normal is rescaled to square2. It has one of the159 fixed endpoint profiles, and every unprescribed old162 contact other than158 is strictly dotted. This consequence does not allow other new intermediate facets.

Two elementary acute-projection facts provide the actual-face passage. If a spherical subset S has Gram G_S and H is a strict interior point, the coefficients of its projection are G_S^{-1}((H,n_i))_i<0. Therefore H minus that projection is future timelike, orthogonal to S, and strictly negative on every other wall. Thus each spherical rank20 subset is an actual finite vertex. Also, a new acute wall is strictly negative at every old finite vertex: writing it as its projection onto the old spherical corner plus a multiple of the future unit vertex, a nonpositive vertex coefficient would make its product with every common strict center positive. Finite vertices of an acute chamber are simple. One direct reason is that a dependent positive-semidefinite acute Gram has a nonnegative kernel vector on a connected singular component, contradicting the negative product of all its normals with a strict center.

For an old edge support W, the old corner W+158 and every proposed W+r are spherical rank20 subsets. They are actual finite vertices in any coherent extension. At most one distinct cap can have this same W: two would give three finite vertices on the same convex one-dimensional face, while coincident vertices contradict simplicity.

Caps on different old edges cannot meet, by the exact pair theorem. This remains true intrinsically in T-perpendicular space: orthogonal projection subtracts the nonnegative Stieltjes correction from their mutual product and decreases their positive squared norms. The same observation excludes every unprescribed strictly dotted old wall from being a side of the new T+r face. The six cusp-completion walls miss the entire T-space by the frozen height gate.

Wall158 contributes no side either. Its old T-face is the compact triangle with the three old finite vertices. These vertices remain strictly inside every cap, so their whole positive cone, and hence the compact triangle, remains strictly inside. A finite ambient angle with158 therefore does not create an intrinsic intersection with that old triangle.

Each added T+r facet consequently has only its two prescribed old-partner sides, meeting at the actual finite corner W+r. Its hyperbolic polygon is a wedge, with an open ideal arc. Walls outside T∪{r} cannot vanish identically on that facet, by finite-vertex simplicity. Their finitely many additional zeros on the ideal arc can therefore be avoided. At a remaining arc point every wall other than the intentional facet support is strict. A small spacelike displacement within the facet remains feasible. Perturbing this vector toward the common strict interior point makes all wall inequalities strict while retaining positive square and the future affine chart. This is a strict spacelike escape vector for the full ambient chamber, so its volume is infinite.

The empty batch retains the already certified168-wall escape. Thus every finite batch in this precisely defined family fails to close the chamber. The root agent separately read this geometric argument and found no gap; no additional arithmetic reconstruction was performed for that read.

## Earlier finite and numerical diagnostics

These auxiliary outputs motivated the exact theorem and retain their narrower scopes:

- [Finite-pattern source](../scripts/fresh_h20_e8_endpoint_belt_20261017.py) and [profile data](fresh_h20_e8_endpoint_belt_20261017_profiles.json):50,390 spherical corner patterns from all159 fixed profiles and mutual labels2,3,4,5,6,8,10. Other mutual I2 labels are not claimed exhausted by this preliminary census. The intended seven-facet belt has V=10, E=15 and F=7, so Euler counting does not obstruct it.
- [Numerical conic screen](../scripts/fresh_h20_e8_endpoint_belt_20261017_geometry.py) and [data](fresh_h20_e8_endpoint_belt_20261017_geometry.json): all5,057 mutual-angle cases fail coefficient tests when all unprescribed old contacts are dotted. These coefficients and conic endpoints are ordinary floating-point evidence; the exact theorem above supplies a separate proof.
- [Finite158 source](../scripts/fresh_h20_e8_endpoint_belt_20261017_finite158.py) and [data](fresh_h20_e8_endpoint_belt_20261017_finite158.json):10,127 numerical caps at158 labels2..64 or parallel, checking all168 ambient gates. Their cross-edge pairs consist of33,998,058 dotted pairs and48 non-Coxeter gaps, with no finite/parallel pair. The numerical conic domain retains the inherited chosen-center restriction; it is not the basis for the exact coherence argument.
- [Singleton151 source](../scripts/fresh_h20_e8_endpoint_belt_20261017_param151.py) and [data](fresh_h20_e8_endpoint_belt_20261017_param151.json):5,951 further numerical caps with151 labels3..64 and158 labels2..32 or parallel. Across the new cross-edge pairs,51,982,158 are dotted and72 are non-Coxeter gaps; none is finite or parallel. The unbounded151 family remains open.

Each of these four auxiliary generators has passed generation and a fresh saved replay. The separate frozen [simultaneous old-facet event diagnostic](fresh_h20_e8_endpoint_events_20261017.md) permits finite contacts to every retained wall and tests actual intrinsic intersections with the seven earlier H3 facets. It finds no accepted cap in its stated bounded domain, but is not an all-label theorem.

The obstruction therefore concerns finite batches of fixed-profile endpoint caps with dotted other old contacts. A future construction must change that contact pattern, use the remaining singleton151 families, or introduce intermediate facets with a different incidence pattern.

A separate [all-label singleton-pair calculation](fresh_h20_singleton151_pair_gap_20261017.md) excludes the recurrent singleton151/profile37 pair for every integer label, by a monotone algebraic formula. That small result does not complete the remaining singleton151 family and is not needed for the finite-batch theorem above.
