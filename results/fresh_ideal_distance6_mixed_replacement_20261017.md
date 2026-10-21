# Mixed finite-nine and ideal-six quadrilateral replacements

**Completed local proof; final verification status below.** Use unit finite payments for marked distance at most9. At ideal corners use unit payments for marked distance at most6, or for distance7,8,9 when an endpoint is within6 graph steps of a feature. A feature is a vertex of degree at least3 or an endpoint of a bond of label at least4.

The following statements hold for actual quadrilaterals in finite-sided finite-volume hyperbolic Coxeter polytopes, without an arithmetic or simple-cusp assumption.

* An ideal marked pair at distance7 in an affine cycle of order at least15 has another corner paid by finite distance9 or ideal distance6.
* At cycle order14, the same conclusion holds unless the face has two sealed finite corners of types A1/A13 and A2/A13, with the A13 mark its unique middle. A half-payment at each of these whole-component pairs pays the face.
* For a classical affine tree, an ideal marked pair at distance7,8 or9 whose endpoints are farther than6 from all features has another corner paid by the rule above.

The [source](../scripts/fresh_ideal_distance6_mixed_replacement_20261017.py) and [certificate](fresh_ideal_distance6_mixed_replacement_20261017.json) also give the all-rank costs of this basic ideal rule. They do **not** establish a new dimension bound or remove the separate preferred-component obligations in small factors. In particular the remaining order15 D-factor payment obstruction is not resolved here. Triangle payments and the total finite budget are separate tasks.

## Actual face setup and inherited classification

Write the face star as H=W+a+b+c+d, with opposite pairs {a,c},{b,d}. If a corner has a nonzero adjacent black product it already pays, so in the unpaid analysis all adjacent black products are zero. At a chosen ideal corner, W+a+b is a varying affine factor F plus a positive-definite spectator diagram Z orthogonal to F. The [actual factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md) gives this description even at nonsimple cusps.

The positive affine kernel represents the actual ideal point. The opposite sides c,d do not contain it, so their heights h,k are strictly positive. We retain the complete finite/affine one-root profile lists and normalized Schur determinant from the frozen [weak-profile proof](fresh_distance8_ideal_corner_replacement_20261017.md) and [mixed finite-nine/ideal-seven proof](fresh_ideal_distance7_mixed_replacement_20261017.md). These classification inputs do not assume a global exclusion of small affine subdiagrams.

Eliminating the spectator block changes the normalized determinant by

    alpha * (x/h-y/k)^T G_Z^(-1) (x/h-y/k) >= 0,

where alpha>0 and x,y are nonnegative native contact-magnitude vectors of c,d with Z. Thus a strictly positive spectator-free determinant contradicts the actual face star's corank1. A negative relaxed determinant alone proves nothing; the small-cycle analysis below handles its equality condition exactly.

## Cycle factors of order at least17

Choose the shorter old marked arc as P=7. The other arc has Q>=10. The complete long-A profiles are isolated, endpoint3 on either arc, endpoint4 on either arc, second-node3 on either arc, and both-end3.

A single-attachment profile gives an adjacent component containing the old A_(m-1) path and the new root. The old endpoint opposite the attachment stays a simple leaf with a degree2 neighbor, since spectators are orthogonal to the old factor. No large affine A/B/C/D diagram has this configuration, and exceptional affine diagrams are too small. Therefore this adjacent marked component and its actual corner are finite. This is the fixed-leaf argument already proved in the [cycle16 transfer](fresh_cycle16_half_coupon_transfer_20261017.md); it applies here to old paths of rank at least13.

For an attachment at position j=1 or2 on an arc of length L, the finite marked distance is L-j+1. Under unpaidness it must exceed9, forcing L>=9+j. Thus all single profiles on the short arc7 are excluded. On the long arc an endpoint profile requires Q>=10 and a second-node profile Q>=11.

Only five profiles remain. The normalized determinant is nonincreasing in the reciprocal arc lengths and concave in the reciprocal heights, as in the frozen formal identity. For every ordered profile pair, use Q>=11 if a second-node profile occurs and Q>=10 otherwise. The checker evaluates all100 height-rectangle vertices, with exact minimum1/35. Positive heights and arbitrary positive-definite spectators therefore cannot produce a zero determinant. This proves the cycle statement for every order at least17, not just for tested orders.

## The exact spectator gap at orders14,15,16

Here P=7 and Q=7,8,9. Every single profile already yields a finite marked distance at most9. The checker verifies the18 literal positive-definite single-profile diagrams, their marked distances and their fixed opposite leaf. Only isolated and both-end3 profiles remain. The12 ordered profile pairs give48 exact rectangle values. Every pair involving both-end3 has a strictly positive minimum. Hence only isolated/isolated can remain.

For isolated profiles the heights are exactly the native opposite contact magnitudes. They are1 or at least sqrt(2). If either height is at least sqrt(2), the normalized determinant has lower bound

    1-3*alpha > 0,       alpha=1/7+1/Q <=2/7.

Consequently h=k=1. This step precedes the contact gap: no lower bound on coordinates of a divided contact vector is being assumed. The normalized determinant equation now requires

    D = (x-y)^T G_Z^(-1)(x-y) = 4-1/alpha.

Its required values are

| Cycle order | Q | Required D |
|---:|---:|---:|
|14|7|1/2|
|15|8|4/15|
|16|9|1/16|

Any common nonzero coordinate of x,y gives a path c-e-d of length2, already paid. Otherwise their supports are disjoint. Each nonzero native contact magnitude is at least1.

If both vectors are nonzero, choose i in the support of x and j in the support of y. Cauchy in the positive-definite Gram metric gives

    D >= (x_i+y_j)^2 / (4-2G_ij) > 1/2.

Indeed x_i+y_j>=2, while positive definiteness of the two-node principal block gives -G_ij<2. This argument allows any finite labels and arbitrary other spectator components.

If just one vector is nonzero, say x, Cauchy with the coordinate vector e_i gives D>=x_i^2/2>=1/2. Equality forces x_i=1 and x=(1/2)G_Z e_i. Nonnegativity of x and nonpositivity of the off-diagonal Gram entries then force i to be isolated in Z and x to have its single unit contact there. The other vector is zero. If both vectors vanish, D=0, which does not equal any required value.

Thus orders15 and16 are impossible in the unpaid analysis. Order14 forces precisely one simple contact to one isolated A1 spectator, on one side only. This proves the first statement and reduces the second to a completely specified equality core. The strict bound and its equality case are universal positive-definite Gram arguments, not finite enumerations of spectator forests.

## The sealed order14 equality core

Number the cycle0,...,13, take a=0,b=7, and attach c=14 to0, d=15 to7, and the single spectator e=16 to14, all by simple bonds. The reflected configuration is equivalent. Every other contact among these17 roots is zero. Equality in the spectator gap has already shown that all remaining supporting spectators are orthogonal to every root of this core.

The two adjacent finite corners have whole marked components A13/A2 and A13/A1. In both A13 paths the old marked root is6 edges from either endpoint, so it is the unique middle. The third finite corner has the marks in A2 and A1, with unmarked old-path components. The two specified half-coupons contribute1. They refer to whole actual finite components; orthogonality of all remaining spectators prevents enlargement or hidden joining of these components.

The source checks the literal core signature(15,1,1), its radical, the supporting positive-definite white block, every corner and both coupons. It also verifies that the equality configuration is an actual local quadrilateral: a strict common center has norm -204, each wall has a strict timelike facet witness, and the complete face cone has three finite rays and one ideal ray.

For scope, the full chamber is infinite volume. Its complete nonnegative dual cone, subject to one radical balance equation, has58 rays:27 finite,5 ideal and26 exterior. A positive displacement of an exterior ray toward the strict center gives an explicitly checked strictly feasible future null vector. This local equality example is not promoted as a finite-volume construction or a novelty claim. The transfer theorem is conditional on an actual quadrilateral and does not require that this entire small chamber be finite volume.

## Classical trees outside feature radius6

Both marks are on the simple spine. At a doubled end, distance greater than6 from both endpoints of the doubled bond leaves an outer B component of rank at least8. At a fork end, distance greater than6 from the fork center leaves an outer D component of rank at least9. These are literal node counts, not distances that count a doubled bond as an ordinary extra edge.

The same eight complete classical extension profiles as in the preceding corridor proof apply at these thresholds. A B8 component contains a doubled bond, so it cannot acquire an exceptional E extension; F4 and its affine analogue are too small. A D9 component already has more nodes than any finite E diagram or any proper subdiagram of affine E8. A replacement meeting both components produces an even larger diagram. Consequently the large B/D extension classification has no newly available exceptional profile at the reduced radius.

Each profile has cost g=2*l/h+c0/h^2 at most4. The unique cost4 profile is an inner second-node simple attachment; every other profile has cost at most3. With marked spine distance P, the normalized determinant is

    F = 1-(g_c+g_d)/P.

For P=7, an inner second-node attachment produces an adjacent path of length6, so is already paid. Removing it leaves F>=1/7. For P=8, it produces a path of length7. In the core the new marked root is a simple fork leaf, at distance1 from the degree3 attachment node; these three old/new incident edges remain in that actual corner. Spectators may add neighbors to the new root, but the old attachment node remains a feature. Therefore this corner is finite-distance9 paid or ideal-distance7 within the new feature neighborhood. Removing the profile leaves F>=1/4. For P=9, no removal is needed: F>=1/9.

The source checks192 ordered profile/distance cases:30 immediate payments and162 positive gaps. Spectators only increase the determinant. This proves the tree statement in every rank satisfying the actual distance-from-features condition.

## All-rank costs of the basic ideal rule

For order m>=18, the exact counts are

    A: 6m,
    C: 6m-21+min(m-7,16)+min(m-8,16)+min(m-9,16),
    B: 6m-14+min(m-9,15)+min(m-10,15)+min(m-11,15),
    D: 6m-7+min(m-11,14)+min(m-12,14)+min(m-13,14).

The A cycle has exactly m pairs at each distance1,...,6 and no features. For C, the path contributes6m-21 short pairs. At each distance7,8,9, the allowed start positions lie in two endpoint intervals of length8, so their union has min(m-k,16) positions.

For B, the simple spine has m-2 nodes. Short pairs total6m-20, including the one extra pair between the two fork leaves. At each of the three long distances there are two fork-leaf pairs, plus two spine start intervals of lengths7 and8. This gives the constant6m-14 and the three minima with cap15. For D, short pairs total6m-19, including the two fork-leaf sibling pairs. At each long distance there are four leaf-spine pairs, and two spine intervals of length7. For m>=18 no opposite-fork leaf pair has distance at most9. This gives6m-7 and the cap14 minima. The interval-union formulas cover their overlap exactly.

Each expression has successive linear slopes9,8,7,6. Its normalized maximum is85/12 at C order24,36/5 at B order25, and95/13 at D order26. The checker records all12 finite/infinite linear branches and proves the bound on each from its endpoint or positive infinite-tail slack. The332 literal graph counts at orders18,...,100 corroborate the formulas; the proof of the infinite tail is the spine count and linear-branch analysis, not extrapolation.

These are costs of the basic rule only. The inherited preferred-component images may require additional pairs at small orders, and the order14 coupon needs a total finite-budget proof if used globally. This package makes neither inference automatically.

## Verification scope

Generation and fresh saved replay both passed100 large-cycle rectangle values,48 small-cycle rectangle values,18 single-profile finite diagrams,192 tree cases,58 complete equality-core rays and332 affine graph counts. The12 all-rank linear branches and exact spectator-gap equality are also checked or proved as described. Syntax, all11 direct source hashes and all6 local report links passed. A separate full mathematical/source review found no gap in the height-one gate, unscaled spectator gap and equality, sealed finite components, fixed-leaf argument, persistent P8 feature, or all-rank count formulas. It did not independently rerun the100/48 rectangle values or192 tree-profile cases. An additional narrow conceptual review checked the unscaled spectator gap and the required ordering of the height reduction, without reproducing the cycle determinant substitutions. Source, data and report are final/frozen. External refereeing and proof-assistant formalization are not claimed.
