# Local finite-eight and ideal-five quadrilateral transfers

**Completed internally checked local theorem; verification scope is below.** Consider an actual quadrilateral of a finite-sided finite-volume hyperbolic Coxeter polytope. Keep the following payments:

* Finite unit: marked distance at most8, one whole marked connected component of order at most17, or distance9 with a mark a leaf of its whole connected component.
* Ideal unit: marked distance at most5, or distance6,7,8,9 with a mark within5 steps of a feature in the varying affine factor. A feature is a vertex of degree at least3 or an endpoint of a bond of label at least4.
* Finite half-coupons: whole A1 against whole A15 marked at its unique middle; and whole A1 against whole A17 marked at distances5 and11 from its two endpoints.

An ideal cycle corner of order at least15 at marked distance6,7,8 or9 is replaceable using these payments at the other corners. The same holds for an ideal classical-tree corner at these distances when both marks are farther than5 from every feature. The half-coupons occur only in explicitly sealed cycle equality cores. There is no simplicity or arithmetic assumption, and supporting spectator diagrams are arbitrary.

This is a local theorem. The [source](../scripts/fresh_finite8_ideal5_local_replacement_20261017.py) and [certificate](fresh_finite8_ideal5_local_replacement_20261017.json) do not establish a complete high-core overlap classification, remove all small-factor obligations, prove the total finite budget or assert a new ambient dimension bound. In particular this statement alone does not replace the frozen dimension80 proof.

## Inherited actual-star facts and the new issue

Use H=W+a+b+c+d with opposite pairs {a,c},{b,d}, all roots of square2. The white block and W plus any one black are positive definite. At the chosen ideal corner the [actual factor-fiber lemma](fresh_cusp_factor_payment_framework_20261017.md) gives a varying irreducible affine factor F containing a,b and a positive-definite spectator block Z orthogonal to F. The positive affine kernel is the actual cusp ray. The opposite sides have strictly positive heights h,k because they do not contain that ray.

The frozen [finite-eight/leaf-nine theorem](fresh_finite8_cycle_replacement_20261017.md) proves the cycle statement at marked distances7,8,9. In its cycle proof the only alternative ideal payments are direct adjacent bonds or shared spectator paths of length2. All its longer single-profile alternatives are finite, and its sole equality uses the sealed A15-middle coupons. Thus its **cycle** conclusion already survives the present reduction from ideal-six to ideal-five. The tree branch requires a new radius-five argument and is handled separately below.

It remains to analyze a cycle with shorter marked arc P=6 and other arc Q>=9, since its order is at least15. We retain the complete eight one-root profiles and exact normalized determinant from the [weak-profile theorem](fresh_distance8_ideal_corner_replacement_20261017.md). A single profile is finite even after spectators are restored: the fixed opposite simple leaf and its fixed degree2 neighbor exclude all sufficiently large affine completions, as proved in the predecessor.

Its adjacent marked distance is L-j+1. Distances at most8 pay. At distance9, endpoint4 and second-node3 have inverse diagonals m/2 and m/4 at their new marked leaf, both greater than2; no spectator can contact that root, so the whole-component leaf-nine rule pays. Endpoint3 either remains a leaf or has a spectator neighbor. In the latter case, absent a shared neighbor of c,d, the predecessor's unequal-height Cauchy estimate adds at least alpha/(2h^2), or alpha/(2k^2) for one fixed selected side, where alpha=1/6+1/Q. This replaces that endpoint profile's norm constant1 by1/2 and preserves the concavity and reciprocal-arc monotonicity. No two coordinate lower bounds are added.

## The complete P=6 partition

The retained profiles are constant on each of Q=9,10,11,12 and Q>=13. For the last range the determinant lower bound is nonincreasing in reciprocal Q, so Q=13 controls the entire infinite tail. The source evaluates400 exact height-rectangle vertices for100 ordered profile pairs. Ninety-three pairs have strictly positive minima. Seven pairs need the discrete-height argument below.

The only exceptional profiles are isolated and endpoint3 on the long arc. For an isolated profile h is1 or at least sqrt(2); for endpoint3, h=x+1 is1 or at least2. This includes every finite Coxeter label and parallel or dotted opposite products. For each of the seven pairs, the source evaluates both smaller height rectangles obtained by imposing either higher-height alternative. All14 resulting minima are strictly positive. Consequently every exceptional case has h=k=1 **before** the unscaled spectator gap is used.

The normalized zero-determinant equation now requires D=(x-y)^T G_Z^(-1)(x-y) to have the following value:

| Q | Ordered profiles | Required D | Resolution |
|---:|---|---:|---|
|9|isolated / isolated|2/5|excluded by gap|
|10|isolated / isolated|1/4|excluded by gap|
|10|isolated / endpoint3, either order|1/16|excluded by gap|
|10|endpoint3 / endpoint3|0|sealed finite A11|
|11|isolated / isolated|2/17|excluded by gap|
|12|isolated / isolated|0|two A1/A17 half-coupons|

If x,y share a nonzero spectator coordinate, their length-two marked path already pays. Otherwise their supports are disjoint. The [unscaled spectator-gap theorem](fresh_ideal_distance6_mixed_replacement_20261017.md) gives D>1/2 for two nonempty disjoint supports and D>=1/2 for exactly one nonempty support. Thus all five positive required values, each strictly below1/2, are impossible. If D=0, positive definiteness gives x=y; after the shared-neighbor alternative has been removed, both vanish. Every remaining supporting spectator is therefore orthogonal to F,c,d, sealing the complete equality core.

This is an all-label argument. The seven exceptional rectangles are not treated by inferring a positive gap from a negative relaxed minimum; their native height possibilities and exact spectator equality are resolved explicitly.

## Two sealed cycle cores

For Q=10, number the16-cycle0,...,15, with a=0,b=6. At the endpoint3/endpoint3 equality, c=16 attaches simply to15 and d=17 attaches simply to7; their native opposite contacts are zero. All other new contacts vanish. At the opposite corner W+c+d, the marked component is the path on7,...,15,c,d, of order11. Every other spectator is orthogonal to it. Hence this is a whole finite A11 component, paid by the retained small17 rule. In fact the two adjacent finite marked components also have order16 and pay. The source checks the core signature(16,1,1), positive-definite white block and all four corners.

For Q=12, number the18-cycle0,...,17 with a=0,b=6. The equality attaches c=18 simply to0 and d=19 simply to6, with no other contacts. Its two adjacent finite corners have whole components A1 and A17. In the A17 path the old mark is five steps from one endpoint and eleven from the other. Both half-coupons are therefore of exactly the new stated type and sum to one. The third finite corner has its marks in separate A1 components. The source checks signature(18,1,1), the white block, all four corner signatures, both full component decompositions and the two endpoint distances.

These literal cores describe a necessary equality inside an actual quadrilateral. No claim that either complete small chamber has finite volume is needed or made. The whole-component coupon interpretation follows from the proved spectator equality, not merely from a truncated diagram.

## The radius-five tree corridor

Outside distance5 from all features, both marks lie on the undecorated simple spine. The outer components left after deleting a mark have type B of rank at least7 or D of rank at least8. The inner component containing the other mark is larger by at least the six-step marked separation. The classical eight-profile list of the weak-profile theorem therefore still applies, apart from one potential finite/affine classification addition at the smallest outer rank.

The outer D8 component can acquire an affine E8 extension by prolonging a short fork leaf by a simple bond. This entire affine subdiagram consists of the outer whites and c: it lies in the actual **edge** diagram W+c, which is positive definite. It is therefore impossible. The source checks its positive kernel explicitly, with D8 chain0--1--2--3--4--5, fork leaves6,7 and new root8 attached to6: the coefficients are(1,2,3,4,5,6,4,3,2). Finite exceptional extensions cannot have more than8 nodes, and the doubled B component cannot sit inside a simply-laced E diagram. Larger inner components and profiles joining both old components have no newly available exceptional extension. Thus no arbitrary rank cutoff or omitted D8 case enters the profile reduction.

Use the same eight scalar profiles, with maximum single-side cost g=2*l/h+c0/h^2. Its maximum4 occurs only at inner second-node3; all other costs are at most3, with equality only at inner endpoint3 with h=1. The normalized determinant before spectators is

    F=1-(g_c+g_d)/P.

At P=6, inner second-node3 gives an actual adjacent path of length5 and pays. At P=7 or8 it gives a path of length6 or7. A finite such corner pays under finite-eight. At an ideal such corner, the old attachment node retains three incident edges and is a feature one step from the new marked root, so the radius-five ideal rule pays. At P=9, every profile pair already has positive gap at least1/9.

After these immediate payments, all P=7,8 cases have positive gaps. At P=6 the sole zero case is two inner endpoint3 profiles. Put u=1/h and v=1/k. Here

    F=[(3-2u-u^2)+(3-2v-v^2)]/6,

with 0<u,v<=1. Equality forces u=v=1 and zero spectator correction. A shared spectator contact pays by a path of length2; otherwise all spectators are orthogonal to F,c,d. The middle white path A5 with the two new endpoints c,d is then a whole finite A7 marked component. Its marked distance is6, so it pays even by finite-eight, without using the small17 allowance. The two outer finite blocks remain separate.

The source records256 ordered tree cases:45 immediate payments,210 strictly positive gaps, and this one sealed A7 equality. This completes the stated local tree transfer at P=6,7,8,9 in every rank satisfying the feature-distance condition.

## Verification scope

Generation and a fresh separate-process saved replay both terminated exit0/PASS for400 cycle rectangle values, all14 high-height gates for the seven exceptional ordered pairs, two literal cycle equality cores,256 tree cases and the affine E8 positive-kernel exclusion. Syntax, all14 direct input/self hashes and all6 report links passed. Two separate full mathematical/source reads found no gap in the height-first cycle analysis, fixed-side correction, whole-component A11/A17 sealing, radius-five D8 edge exclusion, persistent features or sealed A7 equality. Neither review independently reran the400/14/256 arithmetic, re-enumerated the inherited profile catalogue, or asserted a global bound. The finite/affine classification, actual factor-fiber theorem and formal determinant remain explicitly inherited proof inputs. Source, data and report are final/frozen; external refereeing and formalization are not claimed.
