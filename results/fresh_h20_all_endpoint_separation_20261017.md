# All finite endpoint profiles in the dotted-other family

This package extends the [fixed-profile theorem](fresh_h20_e8_endpoint_belt_20261017.md) to the three unbounded singleton151 families. Under the same contact restrictions, **no finite batch of finite-endpoint caps can produce finite volume**. The restriction is that each cap completes one of the old edges159/160,159/161,160/161 to a spherical rank20 corner, while every unprescribed old162 contact except158 is strictly dotted. Every Coxeter contact with158 is allowed. All168 walls are retained, all added walls must be mutually Coxeter-compatible, and a common future timelike strict interior point is required.

This does not cover parabolic endpoint profiles, additional finite contacts with other old walls, or other intermediate facets. In particular an ideal endpoint is not silently included in the spherical deletion census. No new wall or finite-volume construction is claimed.

The [source](../scripts/fresh_h20_all_endpoint_separation_20261017.py) and [saved certificate](fresh_h20_all_endpoint_separation_20261017.json) pass exact generation and fresh saved replay and are final/frozen. They inherit the frozen fixed/fixed arithmetic, reconstruct the additional interval vectors from the exact21-root core Gram, and check the new low-label and continuous high-label cases. The root agent separately read the mathematical argument and source and found no gap. That read did not repeat the arithmetic: this package uses the original interval engine and remains a primary certificate.

## Complete finite-profile reduction

The literal old supports W have types A12+A4+A2+A1, B12+A4+A2+A1 and A12+A4+A2+A1. In each case the isolated A1 is wall151, which belongs to the fixed rank17 support T. The [complete finite endpoint census](fresh_h20_e8_partner_endpoints_20261017.md) gives57,45,57 fixed profiles and exactly one remaining family per edge: a single contact151 of label m≥3. The latter has square cost2cos²(π/m), strictly below2 for every finite m.

The labels3,4,5 give nine additional fixed profiles. Their possible158 contacts split into labels2..6 and the full continuous tail x158≥9/5, exactly as in the frozen theorem. Every low-label discrete root is checked against the required old162 dotted gates. Tail/tail and discrete/tail products are certified by coefficientwise strict negativity after the conic substitution; discrete/discrete exceptions are certified obtuse or in strict gaps between consecutive Coxeter values. This proves that compatible caps from these low profiles and the old fixed profiles cannot meet inside the T-face.

## Why high singleton labels need only orthogonality checks

Let r have a singleton151 contact of label m≥6. If another cap s actually meets r inside the hyperbolic T-face, the Gram of T+r+s is positive definite. The I2(m) subdiagram on151,r must be isolated in that finite diagram. Indeed, if s has any nonright bond to either of its nodes, the normalized three-node determinant is at most

\[
 1-\cos^2(\pi/6)-\cos^2(\pi/3)=0;
\]

an additional third bond only decreases it. Therefore an actual neighbor must satisfy s·151=s·r=0. Two singleton caps cannot be neighbors, and a fixed profile with a nonright151 contact cannot be a neighbor. Only the35,26,35 fixed profiles right to151 remain.

For this argument a hyperbolic intersection, rather than only an ambient finite angle, is essential. Projection onto T subtracts a nonnegative Stieltjes correction and may turn a finite ambient pair into intrinsically disjoint hyperplanes. The test here uses positive definiteness of T+r+s as the necessary actual-intersection condition.

## A continuous high-label reduction to m=6

Fix an old edge W and its actual finite corner W+158. Let V be the future unit timelike vertex vector. Let p151 and p158 be the projections with respective prescribed products−1 at151 or158, and zero at the other corner roots. With

\[
 d=2\cos(\pi/m),\qquad x=-(r,n_{158}),\qquad z=x/d\ge0,
\]

the coherent norm-two branch is

\[
 \frac r d=p_{151}+z p_{158}
   +\sqrt{Q(z)-2/d^2}\,V,
 \qquad Q(z)=(p_{151}+z p_{158})^2.
\]

Every coherent added wall s is strictly negative on the old finite vertex V. To see this without fixing the interior point, project s onto the spherical corner W+158. Its coefficients are nonpositive. A zero or past-pointing vertex coefficient would make s pair nonnegatively with every common strict center, a contradiction. Thus s·V<0.

For m≥6, d≥√3. Holding z fixed increases the radical as d increases, and therefore decreases (r/d)·s. The largest possible pairing is attained at d=√3. This comparison does not require the comparison root itself to satisfy the other old gates: it is a relaxation used only to bound the pairing.

It consequently suffices to prove strict negativity for the m6 conic and **every real x158≥0**. Its conic parameter starts at the exact x158=0 endpoint. For each fixed profile right to151 on a different old edge, the source checks all coefficients of its pairing with this conic: nine coefficients against the continuous158 tail, or three against an accepted discrete158 root. Every coefficient is strictly negative. Hence r·s<0 throughout the full high-label domain, ruling out the required orthogonality. No high-label grid or limit extrapolation is used.

The exact certificate checks192 high-tail and874 high-discrete cases. The low singleton extension contributes981 tail/tail,4,803 discrete/tail,4,605 tail/discrete and22,537 discrete/discrete cases. All continuous coefficients are strictly negative; the largest upper bounds are below−3.1194 in the high orthogonality checks and below−1.1802 in the low dottedness checks. There are44 accepted low discrete roots, one exact old-wall rejection, and13 strict consecutive-angle gaps. The saved rational bounds are authoritative; these decimals only summarize their margins.

## Actual facets and the finite-batch conclusion

The actual-face argument from the fixed-profile report now applies to the complete finite endpoint domain. Acute projection from a common strict center makes every spherical W+r an actual finite vertex. The old W+158 vertex remains strict for every cap, so at most one cap can occur on each old edge. The old compact158 triangle remains strictly inside all caps. All other unprescribed old162 walls are intrinsically disjoint by the dotted-contact condition, and the six cusp walls miss the T-space by the frozen height gate.

Across distinct old edges, the fixed and low-label pair checks forbid actual neighbors; the high-label reduction forbids the orthogonality necessary for an actual neighbor. Thus every added T+r facet has only its two prescribed old-partner sides. It is a hyperbolic wedge with an open ideal arc. Choose a point on that arc away from the finitely many zeros of walls outside T∪{r}, perturb slightly spacelike within the facet, and then toward the common strict interior point. This gives a strict spacelike vector in the full chamber cone, proving infinite volume. The empty batch retains the old certified escape.

The next distinct endpoint possibility is parabolic: projection cost exactly2 and a new root of the form p+λν through an old edge's forward ideal ray. Such profiles require their own affine classification and compatibility analysis; they are not ruled out here.
