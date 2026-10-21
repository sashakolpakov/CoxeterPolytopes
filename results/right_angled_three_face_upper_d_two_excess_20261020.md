# The upper signed three-face defect is at most twice the excess

20 October2026. Paper lemma proposed by the root coordinator and independently checked here. No scientific job, optimization, checker, Lean or build ran. This is an ordinary incidence inequality for actual right-angled three-faces, not a new dimension bound.

For an actual finite-volume right-angled three-polyhedron write F,s,p for facets, cusps and finite vertices, nu=F+s, and A_j for finite vertices with exactly j ideal-ending edges. Put d=A2-3A0. Use the actual P3 recognition, bad-triangle attachment, one-cusp ten-facet exclusion and nu<=10 classification proved in `manuscripts/right_angled_dimension10.tex`: the only small types are P3,Q8,R7, all with d=0. The local three-positive-excess proof also shows A3>0 in a non-P3 face implies F+2s>=15.

**Theorem.** Every non-P3 three-face satisfies

    d <= 2(nu-10).

Consequently, for the manuscript's ordinary sums,

    D <= 2E.

The known lower bound D>=-30E remains independent. Neither inequality concerns a maximum-height subset.

## 1. Elementary reduction

Let e_FI,e_II count finite--ideal and ideal--ideal edges. Euler and degrees give

    p+4s=T=2nu-4,
    e_FI=A1+2A2+3A3=4s-2e_II.

Since p=A0+A1+A2+A3,

    d<=p,
    d=e_FI-p-2A0-2A3<=4s-p.

Thus

    d <= min(T-4s,8s-T).                        (1)

For nu>=14, the continuous maximum of the right side is T/3, which is at most2(nu-10). For nu12 its integer-s maximum is4, attained at s3 or4; for nu13 it is6, attained at s4. These are precisely2(nu-10). The small types already settle nu<=10.

At nu11, s0 violates polygon area and s1 is excluded by the one-cusp lemma. For s2, (1) gives d<=-2; for s>=4 it gives d<=2. Only s3 remains. Its counts are F8,p6 and total polygon excess

    sum_polygons(u+2v-5)=F+2s-12=2.

Since F+2s=14<15, the three-excess lemma gives A3=0. The parity of e_FI=12-2e_II then shows d is even. Equation(1) bounds it by6. Only d4 or d6 need exclusion.

## 2. Excluding d6

The identities force A2=6, A0=A1=A3=0 and e_II=0. Every A2 corner determines one polygon in which its finite vertex has ideal neighbors on both sides. Bad triangles require an ideal--ideal edge, so none occur.

A nonbad polygon contributes at most twice its excess to this corner count. For the total excess2 here this can also be checked directly. Excess1 types have u+2v=6, and only the alternating two-finite/two-ideal quadrilateral contributes, with two corners. An excess2 polygon has u+2v=7; its only potentially contributing allowable type is a three-finite/two-ideal pentagon, with at most one such corner. A one-finite/three-ideal quadrilateral has at least two ideal--ideal edges and is excluded. Thus at most4 A2 corners are possible, contradicting6.

## 3. Excluding d4

The identities force A2=4,A1=2,A0=A3=0 and e_II=1. There is at most one bad triangular facet: two would share the sole ideal--ideal edge and force the entire polyhedron to be P3 by the attachment theorem.

Every A2 corner again contributes to exactly one polygon. The excess1 possibilities contribute either0 or2, with2 only for an alternating quadrilateral. A single excess2 polygon contributes at most1; the one-finite/three-ideal quadrilateral is impossible because it requires at least two ideal--ideal edges. Together with at most one bad triangle this single-excess2 case cannot account for4 corners.

The only remaining excess partition is two excess1 polygons. Their combined A2 count is even and at most4. Since a bad triangle would add1, accounting for exactly4 forces no bad triangle and two alternating quadrilaterals. These use all the excess and have no ideal--ideal edge. Both facets at the sole ideal--ideal edge must therefore have zero excess. A zero-excess polygon containing such an edge can only be a bad triangle. The attachment theorem then gives P3, a contradiction.

Thus d<=2 at nu11 as well, completing the theorem. The correction s2 bound-2, rather than-6 in a preliminary message, does not affect the conclusion. All uses are actual facets, polygon corners and edges of one three-face; no planarity or geometric type is silently imposed by a count pattern.
