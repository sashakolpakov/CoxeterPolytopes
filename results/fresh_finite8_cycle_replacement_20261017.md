# Finite-eight and marked-leaf-nine replacements

**Completed internally checked local theorem; verification scope is recorded below.** Let Q be an actual quadrilateral of a finite-sided finite-volume hyperbolic Coxeter polytope. No arithmetic or simple-cusp hypothesis is imposed. Use the following corner payments:

* A finite corner pays one if its marked roots are connected at distance at most8, or at distance9 with one mark a leaf of its **whole actual connected component**.
* An ideal corner pays one at marked distance at most6, or at distance7,8,9 with a mark within6 steps of a feature of its varying affine factor. Features are vertices of degree at least3 and endpoints of bonds of label at least4.
* The already proved half-coupon between a whole A1 component and the unique middle of a whole A15 component is retained.

Then an ideal corner in an affine cycle of order at least15 at marked distance7,8 or9 is replaceable by these payments at the other corners. The only equality case uses two of the stated half-coupons. The same replacement holds for a classical affine tree at those distances if both marks are farther than6 from its features.

Adding the existing unit payments on whole finite components of order at most17 and whole affine factors of order at most14 preserves this statement. This package repairs a specific finite-distance-nine dependency of the [mixed ideal-six proof](fresh_ideal_distance6_mixed_replacement_20261017.md). It does **not** establish the complete distance-eight high-core overlap theorem, cover every preferred small-factor pattern, prove a total finite budget, or improve the ambient dimension bound.

The [source](../scripts/fresh_finite8_cycle_replacement_20261017.py) and [certificate](fresh_finite8_cycle_replacement_20261017.json) contain740 exact cycle rectangle values,66 single-profile finite diagrams, three isolated exceptional analyses and192 tree cases.

## Actual star, profiles and spectators

Write the actual face star as H=W+a+b+c+d, with opposite black pairs {a,c},{b,d}. The white diagram W and each actual edge diagram W+one black are positive definite. At the chosen ideal corner W+a+b, the [factor-fiber lemma](fresh_cusp_factor_payment_framework_20261017.md) gives the varying irreducible affine factor F containing a,b and a positive-definite diagram Z orthogonal to F. This remains true at nonsimple cusps.

Normalize roots to square2. A direct adjacent black bond gives a unit payment, so in the unpaid analysis b.c=a.d=c.d=0. Let nu be the positive affine kernel representing the actual chosen cusp. The opposite sides do not contain this cusp, so h=-c.nu and k=-d.nu are strictly positive. Every nonzero native contact magnitude is at least1.

For a cycle of order m, let the two old marked arcs have lengths P,Q, with P the marked distance and Q>=P. Here P is7,8 or9 and m>=15. Deleting a leaves the long A_(m-1) path. The complete one-root positive-semidefinite extension list is isolated; endpoint3 or endpoint4 on either end; second-node3 on either end; or both-end3. Finite exceptional extensions are too small. The [weak-profile theorem](fresh_distance8_ideal_corner_replacement_20261017.md) supplies this classification, exact heights, Schur identity and its all-rank interpretation. Nothing here assumes absence of small affine diagrams elsewhere in the polytope.

For a single attachment at position j=1 or2 on an arc of length L, the adjacent marked distance is L-j+1. The old endpoint opposite this attachment is a simple leaf whose neighbor has degree2. Both remain fixed under all spectator contacts, because Z is orthogonal to the old cycle. The adjacent connected component has at least m>=15 nodes. No classical affine A/B/C/D graph has such a fixed simple leaf with a degree2 neighbor; exceptional affine graphs are too small. By the actual factor-fiber theorem, a connected component containing both marks of an ideal corner would have to be affine. Therefore this adjacent corner is finite, including after all supporting whites have been restored.

Distances at most8 now pay immediately. The new issue is precisely distance9.

## Which distance-nine marks remain leaves?

An endpoint4 profile is a finite B_m component with the new mark c at its doubled terminal. A second-node3 profile is a finite D_m component with c a short fork leaf. Every spectator is orthogonal to every old path root, so a spectator touching this component can touch only c.

The exact inverse diagonal at c is

    B_m doubled terminal:  m/2,
    D_m short fork leaf:   m/4.

These follow by one Schur complement of A_(m-1): its endpoint inverse diagonal is (m-1)/m and its second-node diagonal is 2(m-2)/m. The native attachment squares are2 and1 respectively. Since m>=15, both displayed values exceed2. Any spectator with nonzero contact magnitude at least1 would give negative Schur norm in the principal extension by that root, contradicting positive semidefiniteness of the actual corner. Hence every spectator is orthogonal to c. The new mark remains a leaf of the whole actual component, and the distance-nine leaf rule pays.

An endpoint3 profile instead gives A_m, with inverse diagonal m/(m+1). It may extend. If c has no spectator neighbor it remains a whole-component leaf and again pays. Thus an unpaid endpoint3 profile at distance9 **must** have a spectator neighbor. This necessary contact supplies a useful determinant correction; no finite list of spectator forests is needed.

## An unequal-height spectator bound

Let G be the positive-definite norm-two Gram matrix of Z, and let x,y>=0 be the native contact-magnitude vectors from c,d to Z. If they share a nonzero coordinate, c-e-d is an actual white-interior path of length2, which pays at that adjacent corner. Otherwise their supports are disjoint.

Suppose c is the unpaid endpoint3 profile just described, and choose one of its spectator neighbors i. Then x_i>=1 and y_i=0. For z=x/h-y/k, Cauchy in the positive-definite Gram metric gives

    z^T G^(-1) z >= (z^T e_i)^2/(e_i^T G e_i)
                    = x_i^2/(2h^2) >= 1/(2h^2).

This uses no equality between h and k and no lower bound on other coordinates of z. The analogous bound holds for d. If both sides require neighbors, use **one fixed side only** throughout the height rectangle; two coordinate bounds are not added.

Write A=1/P, B=1/Q, u=1/h, v=1/k and use the inherited profile moments l,m and norm constants c0,d0. Put L=u*l1+v*m1 and M=u*l2+v*m2. Before the spectators, the normalized determinant is

    F=1-2(A*L+B*M)-(A+B)(c0*u^2+d0*v^2)-A*B*(L-M)^2.

Spectator elimination adds (A+B) z^T G^(-1)z. The new bound therefore replaces c0=1 by c0=1/2 for one required-contact endpoint3 profile, or replaces d0 analogously. The modified constants remain nonnegative. Consequently this lower bound is still nonincreasing in A,B and concave in u,v. Its minimum over the permitted height rectangle occurs at a vertex, exactly as in the inherited determinant proof. A strictly positive lower bound contradicts the actual star's corank1.

## Complete cycle partition

For each P, split Q into its permitted values up to10 and the single monotone range Q>=11. No endpoint-only sample is being extrapolated: once Q>=11 the surviving profile list is constant and the reciprocal-arc monotonicity applies. The retained lists follow only from the finite-eight and leaf-nine deductions above. In particular, an endpoint3 profile on an arc of length9 uses the fixed-side spectator correction.

The source evaluates all740 rectangle vertices for185 ordered profile pairs across the11 regions. Of these pairs,182 have a strict gap; the only three others are isolated/isolated, treated in the next section. The strict-pair minima are:

| P | Q or lower bound | Strict-pair minimum |
|---:|---:|---:|
|7|8|11/56|
|7|9|1/9|
|7|10|1/35|
|7|>=11|5/77|
|8|8|1/4|
|8|9|1/18|
|8|10|1/10|
|8|>=11|1/11|
|9|9|1/9|
|9|10|2/15|
|9|>=11|1/9|

The checks use exact Q(sqrt(2)) arithmetic and also verify that the modified-profile formula equals the original formula plus its claimed correction. The66 explicit single-profile diagrams verify positive definiteness, marked distances, the fixed opposite leaf and the inverse diagonals. Their role is corroboration; arbitrary rank is covered by the classification, inverse formulas and monotone tails.

## The three isolated cases

For isolated profiles, h,k are native opposite-contact magnitudes. They are1 or at least sqrt(2). With alpha=1/P+1/Q, if either is at least sqrt(2), the determinant is at least1-3*alpha>0 in all three exceptional regions. Hence h=k=1 before applying the unscaled spectator-gap argument.

The zero determinant then requires

    D=(x-y)^T G^(-1)(x-y)=4-1/alpha.

The three required values are4/15 at(P,Q)=(7,8),1/16 at(7,9), and0 at(8,8). The [unscaled gap proof](fresh_ideal_distance6_mixed_replacement_20261017.md) gives D>1/2 for two nonempty disjoint supports, D>=1/2 for exactly one, and D=0 if both vanish. Thus the first two are impossible without an already paid common-neighbor path.

For(8,8), D=0 implies x=y. A nonzero common coordinate pays by a length-two path; otherwise every spectator is orthogonal to c,d as well as F. This is exactly the sealed cycle16 equality core of the [half-coupon theorem](fresh_cycle16_half_coupon_transfer_20261017.md): two opposite simple leaves on the16-cycle, with no further supporting contacts. Its two adjacent finite corners have whole A1/A15 components, with the A15 mark its unique middle. Their two half-payments sum to one. The geometric equality-core and sealing proof is inherited; this package makes no new finite-volume construction claim.

## The tree branch already uses only finite paths of length at most7

The classical tree argument of the mixed ideal-six theorem needs no finite-distance-nine payment. Outside feature radius6, the outer finite components have type B of rank at least8 or D of rank at least9, so its complete eight-profile list remains valid.

At P=7, the only cost-four profile, inner second-node3, gives a path of length6 and pays. At P=8 it gives a path of length7. If that corner is finite it pays under finite-eight. If ideal, the old attachment node remains a degree3 feature and the new mark is one step away, so the ideal feature rule pays. Spectators cannot remove its three existing edges. All other profiles give positive gaps at least1/7 and1/4 respectively. At P=9 all ordered profiles already give gap at least1/9. The192 exact cases here recheck those bounds and record the finite-eight reason explicitly.

This proves the stated tree transfer with the cheaper finite rule. Together with the cycle result it removes that specific local obstruction to a future finite-eight theorem. It does not turn the old distance-nine preferred-corner census into a complete distance-eight census: the latter has additional overlap configurations and separate small-factor obligations.

## Verification and review scope

Generation and a fresh separate-process saved replay both terminated exit0/PASS for740 cycle rectangle values,66 finite profile checks, three exact isolated cases and192 tree cases. Syntax, all11 direct source hashes and all7 local report links passed. A separate full mathematical/source review found no gap in the fixed-leaf finite-corner argument, inverse-diagonal sealing, unequal-height fixed-side Cauchy correction, constant-profile monotone tails, isolated cases or tree payment transfer. That review did not independently reproduce the740 rectangle values. An additional conceptual review accepted the unequal-height coordinate Cauchy bound and the decision to use only one fixed side. Source, data and report are final/frozen. External refereeing and formalization are not claimed.
