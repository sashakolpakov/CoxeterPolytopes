# Nonincident cusp carriers: exact multiplicity obstructions

The known finite-volume H17 pyramid has **1,048 unpaid quadrilateral faces** under the fine corner rule, of which987 are compact. Every one has the two affine E8 factors of the same nonincident cusp in its original face star. This disproves two natural proposed bounds: charging at most one unpaid face per edge of the carrier cusp link, or at most six per link vertex per small factor. It does not disprove a stronger bound with a larger constant, or one that uses additional geometry.

Separately, an explicit family of finite-volume H3 Coxeter chambers has arbitrarily many compact quadrilateral faces carried by a fixed-combinatorial-type cusp. Those faces are fully paid. Thus a useful global estimate must depend on the **payment deficit**, not merely on the presence of a nonincident carrier.

The [checker](../scripts/fresh_nonincident_carrier_incidence_20261017.py) and [data](fresh_nonincident_carrier_incidence_20261017.json) preserve the distinction. No dimension bound is asserted.

## Exact chamber and complete two-face ledger

Use the19-root shared-L25 Gram from the frozen [selective-distance13 proof](fresh_selective_distance13_dimension146_20261017.md) and its [independent reconstruction](fresh_selective_distance13_dimension146_20261017_review.md). Its nodes split into

* L={0,...,8};
* R={10,...,18};
* the shared node9.

The Gram radical is

`(-6,-3,-4,-2,-5,-4,-3,-2,-1,0,6,3,4,2,5,4,3,2,1)`.

In covector coordinates p_i=x·n_i the chamber cone is p_i<=0 with z·p=0. Its complete extreme rays are the single-coordinate ray supported on9 and the81 two-coordinate rays supported on one node of L and one of R. A rational18-by18 inverse rechecks all squares and the strict common timelike center. There are80 finite vertices and two ideal vertices: the apex supported on9, and the base vertex supported on{3,13}. The frozen chamber proof supplies the proper cusp links and actual-facet checks; the new calculation independently recovers all82 rays and their signs.

The projective combinatorics is a pyramid over Delta8×Delta8. A two-face allows four covector coordinates, at least one from each of L,R; its supporting whites are the other15 roots. This is exhaustive: the restricted cone has dimension3 exactly in these cases, since its positive and negative z coordinates impose one independent relation. The actual sides are the coordinate inequalities whose zero sets contain an edge, not every allowed coordinate. The checker also recognizes every15-root support as a finite forest.

There are3,456 two-faces:

| Family | Triangles | Quadrilaterals |
|---|---:|---:|
| In the base Delta8×Delta8 |1,512|1,296|
| Containing the apex |648|0|

The two affine E8 factors at the apex are exactly L and R. Their positive kernel coefficients are the absolute radical coordinates on the corresponding factor. Each positive root combination has covector−e_9, so both determine the same actual cusp. Its Euclidean link is Delta8×Delta8, with81 vertices and648 edges. The other cusp has one affine D16 factor and contributes no small factor of rank at most8.

For each E8 factor separately, the complete incidence ledger is:

| Two-face relation to the carrier apex | Count |
|---|---:|
| Incident triangles |324|
| Nonincident triangles |756|
| Nonincident quadrilaterals |1,296|
| Unpaid nonincident quadrilaterals |1,048|
| Compact unpaid nonincident quadrilaterals |987|

Every quadrilateral has all19 normals in its original wall star. Both E8 carriers are therefore present, and neither is incident to that quadrilateral.

## Literal quadrilateral payments

At a finite corner the checker uses Vinberg's actual rule: weight2 in a connected component of order at most4; weight1 in orders5 through7; in a larger classical component, weight1 for adjacency or membership in a connected end segment of order at most6 invariant under interchange of fork leaves. Otherwise the weight is0. The definition is the one in [Vinberg, §6.1–2, printed pp.80–81](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt), also bound in the [frozen finite-component checker](fresh_compact_dimension_bound_20261015.md). It is not replaced by a short-distance finite rule.

Every connected marked component at a finite quadrilateral corner is A or D. Thus no literal E8 excluded-pair list is needed. At the other ideal vertex, the varying factor is affine D16 and the checker uses the frozen [distance-three and invariant-end-segment rule](fresh_cusp_spectator_role_payments_20261017.md). The small-factor weight2 modification from the [nonincident-defect lemma](fresh_small_cusp_defects_20261017.md) does not alter any quadrilateral corner here: the small-factor cusp is not on these faces.

| Total quadrilateral payment | All quads | Compact quads |
|---|---:|---:|
|0|1,048|987|
|1|243|240|
|2|5|5|

The requested compact quadrilateral has opposite marked pairs{1,2} and{11,12}. Its four omitted pairs and connected marked components are

| Omitted roots | Marked roots | Component | Payment |
|---|---|---|---:|
|1,11|2,12|A17|0|
|1,12|2,11|A16|0|
|2,11|1,12|A16|0|
|2,12|1,11|A15|0|

These are actual corners of the exact compact quadrilateral, not merely spherical subsets of a proposed star.

The complete triangle incidence ledger is included, but **individual triangle fine weights are not claimed**. Some triangle corners use E8; the inherited finite checker records the number of its four omitted pairs without labeling them. The primary text was reopened, but fresh image retrieval did not yield a readable diagram. No excluded pair is guessed from distance. This limitation has no effect on any quadrilateral count or counterexample above.

## What the unpaid example rules out

Write Q_bad for the number of unpaid quadrilaterals and sum only over cusps which carry a small factor. This chamber gives

`Q_bad=1048 > 648 = sum f1(link_v)`.

Also, with s_small(v) the number of small factors and V_v the link-vertex count,

`Q_bad=1048 > 6·2·81=972`,

and even the compact unpaid count987 exceeds972. Consequently an inequality

`Q_bad <= c sum s_small(v)V_v`

would require c>=524/81 in general, and c>=329/54 if restricted to compact unpaid quadrilaterals. This remains true if each defective face is assigned to just one of its carriers: all available small carriers lie at the same cusp, and the right side already includes both factors. These are necessary tests on a proposed multiplicity theorem, not an upper bound established by this example.

## Unbounded carrier incidence with all faces paid

For every integer N>=1, take the upper-halfspace region

`−1/2 <= x <= 1/2`,

`−1/sqrt2 <= y <= (2N−1)/sqrt2`,

outside the N unit hemispheres centered at(0,sqrt2·k), k=0,...,N−1. In the Lorentz metric diag(−1,1,1,1), norm-one outward normals are

`X±=(1/2,±1,0,1/2)`,

`Y−=(1/sqrt2,0,−1,1/sqrt2)`,

`Y+=((2N−1)/sqrt2,0,1,(2N−1)/sqrt2)`,

`B_k=(k²,0,sqrt2·k,k²−1)`.

The pairings give:

* X±·B_k=−1/2, hence angle pi/3;
* B_j·B_k=1−(j−k)², hence adjacent caps meet at pi/2 and all others are dotted;
* the end cap meets its end y-wall at pi/4; other cap/end-wall pairs are dotted;
* the vertical x/y pairs are right, and each opposite vertical pair is parallel.

Every cap is an actual facet. Its visible patch lies over the rectangle x∈[−1/2,1/2], y∈[(2k−1)/sqrt2,(2k+1)/sqrt2], clipped only by the two global end walls when appropriate. It is a compact quadrilateral. The bottom is the upper envelope of these hemispheres. The nearest center has horizontal squared distance at most1/4+1/2=3/4, so the bottom height is at least1/2 everywhere. The finite horizontal rectangle and the integral of z^(−3) from this positive lower height to infinity prove finite volume. Only the point at infinity is ideal; its link is a bounded rectangle.

All finite vertices are

`(x,y,z)=(±1/2,(2j−1)/sqrt2,1/2)`, j=0,...,N.

Thus there are2N+2 finite vertices and one proper ideal vertex. Endpoint finite corners are B3, interior finite corners A3. Each cap's four planar angles therefore has fine weight2, total8. The x-wall pair is the same affine A1 factor in every cap's star and its cusp is not on any cap. The other cusp factor is also A1, independent of N, and the link has four vertices.

This is an all-N symbolic argument; the six exact samples N=1,2,3,7,16,64 only check the displayed identities. It proves that no bound on **all** such carrier incidences can depend solely on cusp factor ranks or the combinatorics of its link. It does not refute a bound on unpaid incidences: these N faces are fully paid.

## Verification and scope

Generation and a fresh saved replay pass:82 exact rays,3,456 two-face incidences, all1,296 quadrilateral payments, and six exact strip checks. Syntax, ten direct source hashes and all local report links pass. A separate internal symbolic-geometry review of the strip confirms its angles, complete vertex description, positive bottom height and sole proper rectangular cusp; it does not duplicate the H17 payment census. The known H17 chamber is not a new construction, and the strip family is used only as an elementary multiplicity obstruction; no novelty claim is made. The earlier dimension proofs and their frozen files are unchanged. Source, data and report are final/frozen.
