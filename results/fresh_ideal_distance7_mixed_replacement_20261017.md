# Mixed finite-nine and ideal-seven local payments

**Completed internally reviewed local results.** The following statements concern actual faces of finite-sided finite-volume hyperbolic Coxeter polytopes, with no arithmeticity or simple-cusp assumption.

1. Every triangle has two distinct adjacent black pairs joined by white-interior paths of length at most7.
2. Suppose an ideal corner of an actual quadrilateral has marked distance8 or9 in a varying affine factor of order at least17. Either one mark is within7 graph steps of a feature, or some corner is paid by finite distance at most9 or ideal distance at most7. A feature is a vertex of degree at least3 or an endpoint of a bond of label at least4.

The [source](../scripts/fresh_ideal_distance7_mixed_replacement_20261017.py) and [certificate](fresh_ideal_distance7_mixed_replacement_20261017.json) prove the finite exact parts below. These are local payment statements, not an ambient dimension bound. In particular they do not replace the old whole-component payments on factors of order16 or17 without an additional argument.

There is also an exact order16 obstruction to the naive local rule: an actual quadrilateral in an infinite-volume H16 Coxeter chamber receives no payment from finite distance9 and ideal distance7 with whole components only through order15. Its complete cone and strict escape are certified below. It is **not** a finite-volume counterexample to statement2 or to a proposed global theorem.

## Two length-seven paths at every triangle

The triangle star has the supporting positive-definite white diagram W and three black roots. It is indefinite, while deleting any black gives the positive-semidefinite star at an actual corner. Consequently a minimal indefinite principal subdiagram L contains all three black roots. The inherited Lanner/quasi-Lanner classification bounds its order by10, as in the [completed dimension108 proof](fresh_selective_distance9_dimension108_20261017.md).

If L has order at most9, take a minimal spanning tree connecting its three marks. Its three arms have total length at most8. If the branch point is one of the marks, the two nonempty arms each have length at most7. Otherwise all three arms are positive; adjoining the shortest arm to either other arm gives two paths of length at most7. These paths contain no third mark internally, so their interiors consist of whites. This argument also covers every lower-order minimal core, independently of labels.

Only the three order10 diagrams in the [complete frozen catalogue](fresh_high_rank_quasilanner_census_20261017.md) remain. The checker tests all360 black triples in these three literal graphs. For each pair it deletes the third black before finding a shortest path. Every triple has at least two such paths of length at most7. The largest second-smallest length is exactly7. These paths lie in the corresponding actual corner diagrams. At an ideal corner the factor-fiber theorem puts both marks and the path in the varying affine factor. Thus two triangle payments survive lowering ideal short distance to7.

## Quadrilateral setup and cycle factors

Use H=W+a+b+c+d with opposite pairs {a,c},{b,d}, and square-two normals. Assume that no corner receives finite-distance9 or ideal-distance7 payment. In particular there is no adjacent black bond: such a bond would have distance1 at its actual corner. Hence b.c=a.d=c.d=0.

At the selected ideal corner, W+a+b is its varying affine factor F plus a positive-definite spectator forest Z orthogonal to F. The positive affine kernel gives the actual future ideal point nu. The opposite sides c,d do not contain it, so h=-c.nu and k=-d.nu are positive. These facts hold at nonsimple cusps by the [actual factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md).

We use the complete eight long-A profiles, weak height bounds and four-root determinant from the frozen [ideal distance-nine replacement](fresh_distance8_ideal_corner_replacement_20261017.md). No small-affine exclusion is restored here. For a cycle, write its two old marked arcs as P,Q. Their minimum is8 or9, and P+Q≥17. It is safe to relax this to P,Q≥8 and P+Q≥17.

A single endpoint3/endpoint4/second-node3 profile for c has an important extra consequence. In the adjacent corner W+b+c, the old cycle minus a is a path of rank at least16. The old endpoint opposite the attachment remains a simple leaf, and its neighbor remains degree2: both are old F roots, spectators are orthogonal to them, and the single attachment is at the other end. This fixed leaf and neighbor cannot occur in a connected affine A/B/C/D diagram. Exceptional affine diagrams are too small to contain the old path. Thus the full component containing b,c is finite, even after arbitrary spectators are restored. The actual factor-fiber theorem then prevents this corner from being ideal: at an ideal face corner its two marks would have to lie in the varying affine component.

Let the attachment be on an arc of length P, at position j=1 for an endpoint or j=2 for a second node. The finite marked distance from c to b is P−j+1. If it were at most9 this corner would already pay. Therefore an unpaid single profile forces

    P >= 9+j:  P>=10 at an endpoint,  P>=11 at the second node.

This handles attachments on a long arc as well as on a short one; it does not incorrectly assume that every single-profile corner has distance at most9. The same reasoning applies to d. Isolated and both-end3 profiles impose no such additional lower bound.

Use the weak heights h>=1 for isolated, endpoint3 and second-node3; h>=sqrt(2) for endpoint4; and h>=2 for both-end3. For each ordered profile pair, impose its arc lower bounds. If both lower bounds remain8, the integer condition P+Q>=17 splits the domain into P>=8,Q>=9 or P>=9,Q>=8. Otherwise the lower bounds already sum to at least17.

With A=1/P, B=1/Q, u=1/h, v=1/k, and the inherited moment data, the normalized determinant is

    F=1-2*(A*L+B*M)
        -(A+B)*(c0*u²+d0*v²)-A*B*(L-M)².

It is nonincreasing in A,B and jointly concave in u,v. The checker evaluates the four height-rectangle vertices for every resulting reciprocal-arc box. All272 exact values are at least1/18. Therefore F>0 throughout the domain, contradicting the full face star's corank1.

As before, eliminating arbitrary spectator whites changes this normalized determinant by

    alpha*(p/h² - 2q/(h*k) + r/k²) >= 0

for a positive-semidefinite [[p,q],[q,r]] and alpha>0. It cannot remove the contradiction. Thus the cycle part of statement2 is complete for every factor order at least17.

## Classical trees away from features

In these orders only affine A/B/C/D types occur. For a B/C/D factor suppose both marks lie farther than7 from every feature. They are on its simple spine. At a doubled end, the nearer doubled-bond endpoint is a feature, so the outer finite B component has rank at least9. At a fork end, the fork center is a feature, so the outer finite D component has rank at least10. These are literal node counts, not distances that accidentally count a doubled edge as a simple extension.

The eight complete large-B/D extension profiles and their all-rank Schur constants are those proved in the [preceding replacement lemma](fresh_distance8_ideal_corner_replacement_20261017.md). They include both one-component affine alternatives, endpoint4 and second-node3. Each root's normalized scalar cost g=2*l/h+c0/h² is at most4; the only profile attaining this bound is inner second-node3. The next-largest maximum is3 for inner endpoint3. The determinant after eliminating the factor whites is

    F=1-(g_c+g_d)/P,

where P is the old marked spine distance.

For P=9, F>=1/9 directly. For P=8, an inner second-node3 attachment gives a path from the new root to the opposite old mark of length P−1=7. That actual corner is already paid, whether finite or ideal. Excluding that profile for either root leaves g_c,g_d<=3, so F>=1/4. The checker records all128 ordered profile/distance cases:15 immediate length-seven paths and113 strictly positive gaps. Spectators again only increase the determinant. This proves the tree part, and therefore statement2 with feature radius7.

## Exact local obstruction at cycle order16

Take a simple cycle on0,...,15, add a simple leaf16 at0 and a simple leaf17 at8, and set every other contact to zero. Let G be its norm-two Gram. Its exact signature is(16,1,1), so the quotient by its radical is a Lorentz space for H16. The checker gives a vector h satisfying G h=-1 and h²=-204. Thus there is a strict common interior. For every root i, h+e_i/2 is a timelike witness on its wall and strictly inside every other wall, proving that all18 walls are actual facets.

Set W={1,...,7,9,...,15}, with black order(0,8,16,17). W is positive definite. Its four projected side normals have Gram

    [[ 1/4,-1/4,-1, 0],
     [-1/4, 1/4, 0,-1],
     [-1,    0,  2, 0],
     [ 0,   -1,  0, 2]].

The complete face cone has four rays: three finite vertices and one ideal vertex. Hence this is an actual quadrilateral. At the ideal corner the marked roots0,8 are antipodal at distance8 in the16-node affine A cycle, which has no features. At each other corner the marked roots lie in different finite components: A15+A1, A15+A1, and A7+A7+A1+A1. Thus finite-distance9 does not pay any of them. The ideal-short7/whole-factor15 rule pays no corner of this actual quadrilateral.

The full chamber is nevertheless **infinite volume**. In dual contact coordinates y=-Gx, the cone is y>=0 subject to one radical relation. Its extreme rays are the two zero-coordinate axes and the64 balanced pairs of positive and negative radical coordinates. The exact census gives35 finite,7 ideal and24 exterior rays. This is a complete cone census because a nonnegative vector in one homogeneous balance equation decomposes into those axes and balanced pairs.

For a particularly explicit escape, let v be the ray with y=e4. Then v²=1 and h.v=-9. With

    t=(sqrt(285)-9)/204>0,

the vector v+t*h is null, pairs negatively with h and has strictly negative product with all18 wall normals. It gives an open ideal escape region. No finite-volume realization or completion is claimed. This example shows exactly why the local mixed replacement theorem cannot simply include factor order16 without additional payment or global geometry.

## Verification scope

Generation and fresh saved replay both terminated successfully, checking the272 cycle rectangles,128 tree cases,360 marked triples and66 cone rays. Syntax, all11 direct source hashes and all7 local report links also passed. A separate internal mathematical/source review read this complete source/report and the inherited replacement proof, finding no gap in the triangle argument, cycle finite-corner and arc restrictions, rectangle argument, radius-seven tree reduction, spectator correction, or order16 cone and escape. That review did not independently reimplement the arithmetic or rederive the inherited classification tables. The source binds the frozen complete high-core classification and the existing determinant/factor-fiber proofs; it makes no new all-label enumeration or ambient dimension claim. External refereeing and proof-assistant formalization remain separate. Source, data and report are final/frozen.
