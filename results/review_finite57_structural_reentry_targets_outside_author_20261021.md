# Outside-author review: finite57 structural reentry targets

21 October 2026. **PASS for the three analytic deductions and the stated
conditional research reduction. No new dimension ceiling is proved.**

Reviewed source: `results/finite57_structural_reentry_targets_20261021.md`,
SHA256 `822953841f6b551e590496b6d617afb67c5f5eceba37ad1fd9566c26d9ef53f7`.
Independent source/manual mathematics only. No enumeration, optimizer,
scientific, Lean or TeX job; no author source or manuscript edits.
The current TeX hash and all five report hashes at source lines139--160
match their actual files. I read the relevant shape, triangle and
concentration proofs as well as the displayed diagnostic boundary.

## 1. Local link sharpness: correct

Source lines25--41 give an actual local construction for every Cartesian
product of simplices. Prescribing heights0,...,r_i on one simplex is
possible because its vertices are affinely independent. At its vertex k,
exactly k of that factor's edges point downward. In the product, height
j=sum k_i therefore equals indegree j; the number of vertices at that
height is exactly the coefficient h_j of product(1+...+t^r_i).

Every edge has distinct endpoint heights before perturbation. A small
generic perturbation consequently preserves all edge orientations and
all index values. With a cut near d/2, vertices of index below d/2
are below it and vertices above d/2 are above it. If d is even, the
middle-index vertices can move to either side: both assignments contribute
d/2. The perturbation can avoid every vertex and still retain all
other side assignments. Each uncut edge is counted exactly once, giving
E0=sum h_j min(j,d-j)=G, not merely an asymptotic approximation.

This establishes sharpness from the local product shape and arbitrary
linear section alone. It does not assert simultaneous realization of
these sections at every vertex of one hyperbolic chamber, nor geometric
realization of a costly formal forest or negative-margin cusp. The
source correctly leaves global coupling open.

## 2. Functional interpolation: correct for every real t>=5

Source lines43--52 concern the existing nonnegative payment rule w.
That nonnegativity, explicit in the current manuscript, is relevant
for faces with at least five sides. For s3 or4 the old face payment is
at least5-s; for s>=5, its nonnegative payment is also at least5-s.
Thus the total new payment on an s-gon is at least

    (t/5)(5-s)+s(t-5)/5=t-s.

The new individual weights are nonnegative since t>=5. At a vertex
there is one corner per link edge, so the added constant contributes
exactly (t-5)E/5. Independently, the incidence/index argument gives
kappa_t=tG/2-E. Subtraction yields identically

    kappa_t-C_t=(t/5)(kappa_5-C).

This applies to finite and ideal vertices, with their own E and G.
Since t/5>0, every local obstruction keeps its sign. Neither integer
t nor a special cusp rank is required. The source correctly restricts
this impossibility statement to the specified automatic interpolation;
it does not rule out a new face inequality with additional geometry.

## 3. Uniform affine-A triangle exclusion: correct

Source lines56--76 use the actual varying affine A cycle F of order
m>=10 and marks a,b at cyclic distance at least5. At the other two
corners of the actual triangle, (F-a)+c and (F-b)+c are principal
positive semidefinite diagrams, whether those corners are finite or
ideal. Normalize their roots to square2. The deleted cycle is the
positive definite A_(m-1), with inverse diagonal k(m-k)/m and all
inverse entries nonnegative. Every nonzero actual contact with c has
magnitude at least1. Consequently a contact at a single position with
inverse diagonal>2 already violates the Schur requirement that the
total principal cost be at most2; other contacts cannot cancel it.

For 3<=k<=m-3, the diagonal is at least3(m-3)/m>=21/10>2.
Thus deletion of a restricts all possible contacts to its four
nonmarked radius-two neighbours; deletion of b gives the analogous
restriction. These neighbourhoods are disjoint at cyclic distance>=5,
including m10's antipodal pair. The opposite deletion also forbids a
contact with either mark itself, since its diagonal is at least
5(m-5)/m>=5/2. Hence c is orthogonal to the entire affine factor F.

The positive affine kernel is a nonzero vector on the actual ideal
corner ray. Orthogonality to F gives c zero height there. But the third
side of a genuine triangle does not contain that corner and has strictly
positive height, a contradiction. Arbitrary supporting roots and other
cusp factors do not affect this principal obstruction. This generalizes
the inspected order10/11 arguments without requiring a new core census.
It supplies no quadrilateral payment by itself.

## Conditional target and status boundary: correct

The initial diagnostic claims are bound to the already completed rank57
report, not claimed as a new run: finite formal maximum372 versus364,
the base-charge witness2D17+D23, and41 negative cusp profiles, including
the worst (10,10,10,10,16). The source does not assert realizability or
apply a compact farthest-vertex theorem to arbitrary noncompact corners.

The proposed finite7/ideal4 rule, all-family factor envelope and global
finite budget<=6n in lines78--114 remain explicitly prospective. Pair
counts alone are not used as a coverage theorem; protected triangle
roles, distinct payments, separate positive edge budgets and full-face
compensations remain part of the unproved quadrilateral obligation.

The further conditional concentration reduction in lines116--127 is
valid. For c(r)=min(r/2,4), transferring one rank from a<=b<=7 to b
preserves cost, including b7 becoming8. Repetition leaves at most one
rank below8. For 9<=a<=b, transfer until a8 preserves both plateau
costs4. Convex order, with the other summands retained, does not increase
normalized shape curvature. The increasing bounded integer sum of rank
squares ensures termination. Hence the worst normalized margin reduces
to (ell,8,...,8,b), ell<=7 and b>=8, with zero omitted. This reduction
does not evaluate any of those still-unproved margins.

Finally beta54=324=6*54; beta55-6*55=15/2,
beta56-6*56=14, and beta57-6*57=22. If one simultaneous geometric
rule supplied the stated finite budget and nonnegative ideal margins
in all four ambient dimensions54--57, the positive global-extrema
correction would contradict the payment upper bound, including at the
finite coefficient equality in dimension54. Combined with the accepted
higher tail this would imply ceiling53. It is a correct sufficient
conditional argument, with cusp total ranks53--56, and is not an achieved
bound or an authorization for computation.

No blocker was found. The three proved analytic deductions are distinct
from the three remaining obligations: quadrilateral coverage, the full
finite budget, and the terminal ideal margins. No prospective envelope,
margin or lower ceiling receives verification credit from this review.
