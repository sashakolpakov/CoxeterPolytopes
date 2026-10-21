# Quadrilateral payment at a large affine factor with spectators

**Final internally reviewed ideal-corner payment lemma.** Let P be a finite-sided finite-volume hyperbolic Coxeter polytope. Suppose its diagram has no connected affine subdiagram of order at most9. If an actual quadrilateral has an ideal corner whose varying irreducible affine factor has order m>=18, that quadrilateral receives payment at least1 under the finite Vinberg rule and the affine rule described below. Other irreducible factors at the cusp may be present. This statement does not assert payment for quadrilaterals with four finite corners, or an ambient dimension bound.

The affine rule pays a marked pair at graph distance at most3 or in a fork-invariant end segment of order at most6. Its cost is at most three times the factor order in every rank, by the [role-payment package](fresh_cusp_spectator_role_payments_20261017.md). The finite rule is the original published Vinberg rule. The new [checker](../scripts/fresh_large_affine_factor_quad_payment_20261017.py) and [data](fresh_large_affine_factor_quad_payment_20261017.json) verify the additional scalar and inverse-entry calculations below. Complete extension lists and unbounded parameters are handled in the proof, not by a rank or label cutoff.

## The original face star and the spectator correction

Write the four corners as W+a+b, W+b+c, W+a+d and W+c+d. The supporting white diagram W has order n-2 and is positive definite. Each W+one black supports an actual edge and is positive definite, even if an endpoint is ideal. The full star has rank n+1 and corank1, by projecting its four side normals to the three-dimensional Lorentz space perpendicular to W.

At the chosen ideal corner, W+a+b consists of the varying affine factor F, containing both a and b, together with a positive-definite spectator forest Z orthogonal to F. This is the [actual factor-fiber description](fresh_cusp_factor_payment_framework_20261017.md); it does not assume the ideal vertex is simple. The full cusp may have other incident normals not in the face star. Every normal has square2.

The positive kernel of F gives a nonzero vector nu representing this actual ideal vertex. Since c and d are opposite sides of the quadrilateral, neither contains that vertex. Therefore h=-c.nu and k=-d.nu are strictly positive. This argument replaces orthogonality to the entire simple-cusp normal span; a single factor need not span that null hyperplane.

First eliminate only the white roots of F, leaving the spectators. The frozen cycle and tree formulas give a four-root Schur determinant. Eliminating Z subsequently subtracts a positive-semidefinite matrix [[p,q],[q,r]] from the c,d block. The normalized determinant changes by

    alpha*(p/h^2 - 2q/(h*k) + r/k^2) >= 0,

where alpha>0 is the old marked-pair residual diagonal. The heights h,k are unchanged. This is the exact [PSD spectator identity](fresh_cusp_factor_payment_framework_20261017.md). Consequently every strictly positive scalar contradiction remains valid with any number of spectator roots and arbitrary allowed contacts to them. A finite-component payment requires a separate check that its marked roles remain paid; those checks appear below.

If any of b.c, a.d or c.d is negative, the corresponding adjacent marked pair is paid. In the remaining analysis all three products are zero. The opposite contacts x=-a.c and y=-b.d are arbitrary Coxeter contacts: zero, a finite cosine magnitude, parallel2, or dotted. No field restriction is imposed.

## Cycle factors

For an affine A factor, deleting a leaves the connected A_(m-1) path. Its complete one-root finite/affine extensions are exactly the five classes in the [frozen cycle proof](fresh_simple_cusp_affine_a_quadrilateral_20261017.md). Spectators do not enlarge this principal-diagram list: every connected component of a positive-semidefinite Coxeter diagram is finite or affine.

All height exclusions in that proof use only affine diagrams of orders3,4 and9. They therefore hold under the present hypothesis. The old cycle has m>=18, so its unused arc supplies the same induced affine-E8 leaf witness. Both arc lengths are at least4 at an unpaid corner. The all-rank reciprocal-arc argument gives a strictly positive determinant before spectator elimination, and the PSD identity preserves that contradiction. Thus the complete cycle branch transfers.

## A bounded affine witness for every excluded simple leaf

Consider a classical affine tree with main-path marks a,b, distance P>=4, and outer ends B_r or D_r. D2 denotes the two fork leaves. The factor order is m=r+s+P+1. A new simple leaf at a, its first inner neighbor, or its second inner neighbor has available opposite all-3 path lengths

    L=r-1+j,  R=m-3-L,  j=0,1,2.

These lengths stop before a doubled terminal bond and may end at one chosen fork leaf. If L,R>=2, their sum m-3>=15 ensures that one is at least5. The leaf together with path arms of lengths2 and5 induces affine E8, of order9. Extra fork leaves are omitted and introduce no additional edge among the selected vertices.

Only the following near-end exceptions have L<2. At j=0 they are B1, B2 and D2; the whole short end, a, the new leaf and the first inner neighbor give respectively an affine four-node doubled star, affine B4, or affine D4. At j=1 the only exception is B1, giving affine B4 on five nodes. There is no exception at j=2. The right path has length at least4-j. Thus every simple-leaf exclusion needed here uses an explicit induced affine diagram of order at most9, even when the outer end is arbitrarily long. This replaces the older proof's use of a possibly long affine B/D witness.

## Complete main-path replacements and their bounds

Deleting a separates its outer B/D end from the inner B/D component containing b. A connected one-root extension contains whole old components. The original finite/affine classification in the [short-end proof](fresh_simple_cusp_affine_tree_short_ends_20261017.md) still gives every join touching both components. The one addition is a D2 join using one outer leaf by a bond4 and the inner long port by a bond3, leaving the second A1 separate. It has outer projection data q=1,t=1/sqrt(2), hence g<=2/(1+1/sqrt(2))<14/9. Other extra doubled joins would put a doubled bond in the interior of a large classical component and are excluded by the finite/affine classification.

An outer-only extension is still strictly spherical: outer+c is a subdiagram of the actual edge W+c. Thus q<2, and the arbitrary-contact Cauchy--Schwarz argument from the short-end proof applies unchanged. Its small F4/E8 exclusions have order at most9. The D2 arbitrary-contact bound is unchanged. These give g<=14/9, with no contact-label cutoff.

The inner component may now acquire its own affine extension while the outer end remains separate. For a large B/D inner component, the complete additional profiles are a bond4 at its long endpoint and a bond3 at its next node. Together with isolated and ordinary endpoint3, their residual data are:

| Profile | c0 | l | h | Bound on g=2l/h+c0/h^2 |
|---|---:|---:|---|---:|
| isolated |2|0|x|1|
| inner endpoint3 |1|1|x+1|1|
| inner endpoint4 |0|sqrt(2)|x+sqrt(2)|2sqrt(2)/(1+sqrt(2))|
| inner next-node3 |0|2|x+1|4/(1+sqrt(2))|

For isolated x=1 and endpoint3 x=0, the preceding leaf argument applies. Isolated x=sqrt(2) gives a doubled star (or affine C2 at the doubled end). Endpoint3 x=1 gives affine A2. Endpoint4 x=0 gives the four-node affine doubled star. Next-node3 x=0 is the interior leaf just treated, and x=1 gives an affine four-cycle. Thus the stated lower height bounds hold in every rank. The case x=0 for an isolated root would give h=0, already excluded by the actual quadrilateral geometry.

For all ordinary profiles, therefore,

    g <= g_* = 4(sqrt(2)-1).

The weighted-tree Green identity is unchanged. The contact supports lie at the first or second surviving path node; at P=4 the two second nodes can coincide, and the same Green entry l*m/P remains valid. Hence

    det(S)/(h^2*k^2) = 1-(g_c+g_d)/P
                     >= 3-2sqrt(2) > 0.

The complete small-inner exceptions remain D6→E7 and D7→E8. This is the same finite/affine deletion classification: exceptional affine extensions have order at most9 and are excluded, while the newly allowed classical affine profiles are the two rows just added. Because m>=18, at most one inner component is small; the other has rank at least9.

## Finite exceptional roles and the extra E8 scalar rows

For D6→E7, the added d is the leaf at the end of the length-two arm. Its inverse diagonal is2. Spectators can touch that E7 only at d; the [role lemma](fresh_cusp_spectator_role_payments_20261017.md) makes every nonzero extension affine E7, which is forbidden here. Thus no spectator touches it, and the original finite E7 payment persists.

For D7→E8, d is again the end of the length-two arm, now with inverse diagonal4. A nonzero spectator contact costs at least4, exceeding square2. Thus every spectator is orthogonal to d and to this entire E8 component. In the D2/P5 case the marked branching vertex is paid by Vinberg's published special-vertex lemma. This argument uses the actual d position, not an arbitrary embedding of E7 or E8.

For the main-path D3/P4 exception, retain the old notation with outer A3 roots L,0,L', a the next main node, and d attached to L'. The frozen formulas give alpha=1/4, A=1/2, D=5/4 and v=1/k<=2/3. All old scalar contradictions survive the PSD correction. There are two new inner profiles and the formerly directly paid outer-only profiles to check.

The complete outer-only finite extensions of A3 are an endpoint bond3,4 or5, or a middle bond3. An endpoint5 cannot coexist with d in the positive-semidefinite opposite corner: its connected five-node diagram contains a bond5 and is neither finite nor affine. Endpoint4 at the same leaf as d gives the forbidden affine doubled star. The remaining cases have the following exact data and height bounds:

| c profile | C | E | upper bound for u=1/h |
|---|---:|---:|---|
| endpoint3, same leaf as d |5/4|3/4|1/(sqrt(2)+1/2)|
| endpoint3, other leaf |5/4|1/4|1/(sqrt(2)+1/2)|
| middle3 |1|1/2|1/(1+sqrt(2))|
| endpoint4, other leaf |1/2|sqrt(2)/4|1/(1+sqrt(2)/2)|

For endpoint3, x=0 gives the leaf affine E8 and x=1 an affine four-cycle. For middle3, x=0 gives affine D4 and x=1 affine A2. For endpoint4, x=0 gives affine B4 with one doubled terminal arm. All these witnesses have order at most9. Substitution into

    F=(1-Bu-Av)^2-alpha*(C*u^2+D*v^2+2E*u*v)

with B=0 gives a strictly positive minimum at the displayed positive upper bounds and v=2/3. The checker calculates all four minima exactly.

For the two new inner profiles E=0 and (B,C)=(sqrt(2)/4,1/2) or (1/2,1). Both have u<=1/(1+sqrt(2)). Their polynomials are respectively

    1-u/sqrt(2)-v-v^2/16+sqrt(2)*u*v/4,
    1-u-v-v^2/16+u*v/2.

Each is affine in u and concave in v, so the four rectangle corners suffice; all eight values are strictly positive. Spectator elimination only increases F. This exhausts the additions in the main-path E8 exception, rather than assuming a small finite component remains small after attaching spectators.

In the fork-leaf P6/E8 exception the old scalar cases are unchanged. The opposite D8 sibling pair permits only an affine-E8 spectator extension by the role lemma, so here it permits none. Its finite end-segment payment persists. The doubled-mate profile still contains a forbidden four-node affine doubled star. No new terminal profile is available beyond the four complete finite/affine profiles already in the frozen tree proof.

## Terminal marks and conclusion

At a doubled terminal the old bound is g<=5/2; at a fork leaf it is g<=4. Their complete four-profile lists and all exclusions use only affine diagrams of order at most9. Replacing the opposite ordinary main-path bound by g_* still gives

    1-(5/2+g_*)/6 > 0,
    1-(4+g_*)/6 > 0.

The exceptional opposite E7/E8 roles were checked above. For two terminal marks, R=m-1>=17 and the old total bound is at most8. The supports are separated and the same Green identity applies. All positive inequalities persist with arbitrary spectators.

Every unpaid large classical ideal-corner case is therefore either paid by a separately preserved finite component or contradicts the actual corank-one determinant. This proves the stated ideal-corner branch.

The no-small-affine hypothesis follows, in particular, if every irreducible cusp factor has rank at least9. Requiring every factor to have rank at least17 would also ensure that every ideal corner falls in the present large-factor branch. That does **not** yet establish a dimension theorem: the preliminary U,V decomposition in the frozen all-finite-corner transfer excluded potentially large affine private diagrams by simplicity, beyond its later local carriers of order at most17. That remaining use is not removed by this report.

## Verification status

Generation passes the generic margins, six new E8 profile calculations, the E8 inverse entry and two explicit affine-star kernels. The all-rank Green identity and old exceptional rows are inherited from frozen proofs, with the changed hypotheses checked above. A fresh saved replay also passes. Syntax, all seven direct source hashes and local links pass. Separate narrow root and construction mathematical/source reviews pass the actual factor geometry, bounded leaf witnesses, expanded profile accounting, E7/E8 roles, six new scalar rows and spectator correction. They do not claim an independent re-enumeration of every inherited finite/affine extension list or a duplicate arithmetic replay. The source, data and report are final and frozen. No previously frozen file is changed.
