# Large affine cycles: full outside-author mathematical review

21 October 2026. **PASS at the stated actual-face input boundary.**
Reviewed source:
`finite57_uniform_large_cycle_four_rule_and_boundary_profile_20261021.md`,
SHA-256 `e2a637b6a96ec650fe161afa5216a131a9493a5653cb160f38fb5b6473e297b8`.
This is an outside-author manual/source audit. I derived the profile,
determinant, whole-corner and boundary arguments before reading any review;
no existing PASS supplied mathematical evidence. No scientific, symbolic,
numerical, Lean or TeX job, authored-source edit or bound promotion occurred.

## 1. Actual geometry and complete native profiles

The inherited framework explicitly allows nonsimple cusps. For an actual
two-face, the supporting n-2 roots are elliptic, every edge adds one
root to a positive-definite set, and an ideal corner restores one WHOLE
varying affine factor F. All other selected roots form an orthogonal
positive-definite forest Z. This is a selected actual face star, not a
claim that its n selected roots exhaust every incident cusp facet.

The four projected side normals span a Lorentz three-space, so the
complete star has negative index1 and corank1. The positive cycle kernel
is an actual nonzero cusp ray. Each opposite side misses that ideal
corner, giving h,k>0. The contact rows on Z may overlap and have arbitrary
rank and native magnitudes; neither h=k nor a shared extension budget
is used. The two edge-positive budgets remain separate throughout.

I checked the port list directly with the norm-two inverse of A_(m-1).
At m>=26, the third-node diagonal is at least69/26>2. The endpoint
diagonal is at least25/26, allowing native1 or sqrt2 but excluding phi;
the second-node diagonal is at least24/13, allowing only1. Positivity
of the inverse excludes all larger contacts. Two unit endpoint contacts
have cost exactly2. Endpoint/second-node pairs at opposite ends already
have cost3-1/m>2; pairs on the same end and two second nodes cost more.
Three native contacts have diagonal contribution at least3(m-1)/m>2.
Thus the five listed types and their two arc orientations are exhaustive.
The omitted opposite-mark contact is not bounded above by this test.

An adjacent black bond gives a DIFFERENT unchanged distance1 unit.
Otherwise all three such products vanish. A single short-arc port gives
an actual b-to-c or a-to-d path of length5 or4. I checked that its output
is a WHOLE finite component after Z is restored: Z can attach only through
the new black root; the remote old endpoint remains a simple leaf whose
neighbor has degree2. No classical affine type has that configuration,
and its contained A_(m-1) path is too large for every exceptional affine
type. All other selected components are positive-definite parts of Z.
Consequently that adjacent corner is finite. Its component has rank at
least26 and its marked distance<=5 unit survives every current deep,
protected-leaf and contextual change. This does not confuse a finite
subdiagram with a whole actual component.

After those outputs are removed, the moments and norm constants in the
five-row table are correct. In particular the heights are the sums of
cycle contacts, so their reciprocal bounds are1,1,1/sqrt2,1,1/2,
respectively. The isolated lower bound uses h>0 plus the native gap;
it does not import the stronger, inapplicable all-simple-cusp exclusions.

## 2. Determinant, spectators and all-rank margin

I derived the Schur identity with A=1/5,B=1/Q from the inverse entries.
The residual a,b block is (A+B)[[1,-1],[-1,1]]. For each arc the
cross term between opposite-end port rows is l_j*m_j divided by its
length. Every retained port lies within two steps of its end; P5 and
Q>=21 leave no overlap exception. A spectator subtraction changes the
normalized determinant by exactly

`(A+B)*(x/h-y/k)^t G_Z^-1*(x/h-y/k)>=0`.

This sign is important: arbitrary spectators strengthen the contradiction.
No individual nonempty-coordinate lower bound is added, and the proof
does not replace the separate strict edge bounds by a sum budget.

The source's scalar F is the normalized necessary corank equation.
For each retained row |u(l1-l2)|<=2, hence the total difference has
absolute value<=4 and its square is bounded by4 times the sum of the
two absolute differences. Substituting the native reciprocal-height
bounds gives the five H upper bounds in lines110–114. All are at most
2(A+B). The only nontrivial comparisons reduce to
A-B-4AB>=4/35>0; the second-node comparison is twice that same
inequality. Therefore, for BOTH profiles simultaneously,

`F>=1-4(A+B)>=1-4(1/5+1/21)=1/105>0`.

This contradicts the actual corank equation for every Q>=21, every pair
of native heights, and every positive-definite Z. Both-end3 profiles
are included directly, even though their new cycle corner also loses
its distance5 entry. No zeroed corner is paid recursively by another.

The all-order triangle lemma also checks directly. Deleting either mark
permits contacts only in its radius-two end neighborhoods; at cyclic
distance>=5 these neighborhoods are disjoint. Each omitted mark is
forbidden by the other deletion's inverse diagonal, at least5/2.
The third side would be orthogonal to the whole positive affine kernel,
contradicting its positive height. Spectators do not repair an excluded
positive-semidefinite principal diagram.

Giving this independent theorem first priority on any face with a changed
actual cycle5 corner is sufficient. Such a quadrilateral receives an
unchanged full unit and such a triangle does not occur. If a face has
no changed actual corner, all its weights and prior coverage remain
unchanged. No earlier route needs to spend a newly zero contribution.
At order m>=26, exactly m unordered distance5 pairs are removed, so the
A-factor total is4m and its fiber charge decreases by V. This changes
neither the other affine families nor a whole-family ideal envelope.

## 3. Exact cycle25 boundary and actual realization

For the stated native3 profile at P5,Q20 I recomputed the four-black
Schur matrix. Eliminating its bottom2I leaves -(1/4)ones(2,2),
of rank1 and negative index1. With W=A4+A19 positive of rank23,
the full27-root Gram has inertia(25,1,1) and rank26. Every W+one
black edge is strictly positive definite. The four complete corner
diagrams are precisely cycle25; A24+A1 with the A24 mark at zero-based
position4; A4+A19+2A1 with marks on distinct A1 factors; and the reflected
A24+A1 corner. Thus this is not a determinant-zero matrix with the
wrong local rank or corner signatures.

After W projection, a^2=b^2=1/4,c^2=d^2=2, a.b=-1/4,
a.c=b.d=-1, and the other black products vanish. U=a+b is null,
V=c+d has square4, and U.V=-2. The point H=U+V/4 has square-3/4
and products -1/4,-1/4,-1/2,-1/2 with a,b,c,d. The adjacent b,c;
c,d; d,a intersections are actual finite right-angle corners; U is
the actual ideal a,b corner and has negative products with c,d.
These halfplanes give an ideal Lambert quadrilateral. A sufficiently
small correction with strictly negative products on W yields a strict
ambient chamber point while retaining every side inequality. Thus the
local acute Lorentz realization is genuine.

I checked the proposed cusp padding literally. For orthogonal norm-two
z_j, beta_j=nu-z_j has norm2, z_j.beta_j=-2, beta_j.c=beta_j.d=-1,
and all other added products are zero. The full old cusp is a product
of cycle25 and the added affine A1 factors, all with the same actual
null ray nu; its Euclidean link is locally bounded. On the quadrilateral
z_j=0, so beta_j restricts to nu and is strictly negative at every
hyperbolic face point and the three finite corners. It meets only the
old ideal corner and neither clips a side nor changes a finite diagram.

I read the CURRENT TeX ordinary and contextual tables, not predecessor
summaries. A1/A24 position4 has no ordinary or contextual payment.
Two distinct marked A1 components also have none, even with unmarked
A4,A19 and arbitrary extra whole A1 padding. No same-component or
protected-leaf value applies to either distinct-component marked pair.
Therefore all four contributions become zero if only the cycle5 unit
is removed. The equality 1-4(1/5+1/20)=0 is exact. This proves the
claimed boundary of these local assumptions and unchanged finite weights;
it does NOT assert a finite-volume completion of the whole chamber.

## Source binding and limits

The inspected current manuscript is `manuscripts/finite_volume57.tex`,
SHA-256 `c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.
The six named inherited sources were read for their indicated mathematical
inputs: actual factor/fiber geometry, Schur identity, native port list,
whole-corner/priority arguments, and the analytic triangle lemma. No old
scientific receipt was treated as a fresh mathematical execution.

No mathematical correction is requested. The result is the A-family
distance-four transfer for whole varying cycles of order>=26, plus the
local order25 boundary. The order10..25, B/C/D-family, finite7/4-budget,
negative-margin-cusp and any ceiling53 obligations remain outside it.
