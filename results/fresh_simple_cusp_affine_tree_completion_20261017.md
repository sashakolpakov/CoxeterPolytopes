# All classical affine-tree quadrilateral corners: proposed completion

**Final internally reviewed lemma: generation and fresh saved replay PASS.** Scope: an all-simple finite-volume hyperbolic Coxeter polytope of dimension n>=18, with the weights in the [frozen interface](fresh_simple_cusp_vinberg_weight_frontier_20261017.md). The theorem proved below is that every quadrilateral having an affine B/C/D corner receives its required payment. The all-finite-corner transfer is a separate obligation; no dimension bound is asserted here.

Normalize every root to square2. Denote the old ideal corner by W+a+b, and the other corners by W+b+c, W+a+d, W+c+d. Assume the old corner is unpaid. Then a,b have graph distance at least4 and are not together in either prescribed six-node end segment. If any of b.c, a.d, c.d is nonzero, the corresponding corner is paid immediately. We therefore impose all three products zero. All one-root extensions below must be spherical, or must be the entire connected affine diagram of order n: a smaller affine component is forbidden by the [global simple-cusp lemma](fresh_simple_cusp_affine_exclusion_20261017.md).

## Cusp normalization and the common Schur identity

In an affine tree give cusp coefficient1 to every main-path node and fork center, coefficient1/sqrt(2) to a terminal node on a bond4, and coefficient1/2 to a fork leaf. Write the null vector as nu=sum nu_i n_i. Replace the marked roots in the Schur calculation by nu_a*a and nu_b*b, leaving c,d at norm2. The scaled old Gram is a graph Laplacian: conductance1 on a main-path or terminal bond4 edge, and conductance1/2 on a fork edge. Its effective resistance between a,b is R. Eliminating W gives the marked block R^-1[[1,-1],[-1,1]].

For c write h=-c.nu>0, and similarly k for d. Positivity follows directly: all old contacts are nonpositive and every cusp coefficient is positive; equality would make c orthogonal to all n cusp normals, hence proportional to their nonzero null direction, contrary to c^2=2. If its contacts are localized at the a end, its residual square and cross contact to the scaled b have the form C=c0+l^2/R and B=l/R. For d use D=d0+m^2/R and A=m/R. Provided the contact supports do not cross along the surviving path, the projection cross term is E=lm/R. Consequently the [formal four-root determinant identity](fresh_simple_cusp_affine_a_quadrilateral_20261017.md) gives

`det Schur/(h^2 k^2) = 1-(g_c+g_d)/R`, where `g_c=2l/h+c0/h^2`.

After eliminating the n−2 independent white roots, all four residual normals lie in a three-dimensional Lorentz space. The residual a,b,c determinant is −R^-1 h^2<0, so their span has rank3 and the four-root Schur matrix has exactly corank1. In particular its determinant is zero. A strict positive lower bound is a contradiction. Here is an explicit all-rank derivation. In the scaled white coordinates the Dirichlet inverse is the Green matrix of the weighted tree grounded at a,b. Put resistance coordinate z along their unique path. For two separated contact supports, with left coordinate z and right coordinate z', the Green entry is z(R-z')/R. Multiplying by the cusp coefficients in the native contact vectors gives E=lm/R. A left diagonal entry equals its one-ended tree value minus z^2/R. A branch off the path adds its common branch resistance to that one-ended value. Thus ordinary first and second path nodes give (c0,l)=(1,1),(0,2); a bond4 to the first gives (0,sqrt(2)). At a fork-leaf mark, its mate lies on a resistance2 branch attached at z=2 and has cusp coefficient1/2, yielding (1,1), or (0,sqrt(2)) for a bond4. These identities give every displayed terminal row. Outer components grounded directly at a have z=0, yielding the arbitrary outer cost q and zero cross moment. The two exceptional E8 supports overlap these local regions and are computed separately. No finite-rank cutoff is used.

There are three exhaustive positions for a marked root: a main-path node, including a fork center; a terminal bond4 node; or a fork leaf. Two terminal marks at the same fork are already paid. The cases below exhaust the remaining possibilities, up to exchanging a,b.

## Main-path marks, including fork centers

Let P be their graph distance. The outer components beyond a have type B_r, r>=1, or D_r, r>=2; D2 here means two orthogonal A1 leaves, with a their common fork center. There are P-1 intervening white path nodes and an analogous outer end beyond b. The total order is n=r+s+P+1 and R=P.

The [short-end proof](fresh_simple_cusp_affine_tree_short_ends_20261017.md) applies to B_r and D_r, r>=3. Its outer inverse formulas and arbitrary-contact argument do not require an upper bound on r. Add D2 with J=I/2, old attachment v=(1,1), and t=Jv=(1/2,1/2). For a nonnegative outer contact vector w let q=w^T Jw<2 and t0=v^T Jw. The strict inequality follows because W+c supports an actual edge and is positive definite, so its subdiagram consisting of the outer block and c is positive definite. If the opposite native magnitude x=-a.c is at least1, q>=t0^2 and t0>=1/2 give `(2-q)/(x+t0)^2<=1`. If x=0 and at least two outer contacts are nonzero, t0>=1 gives q+t0^2>=2. A single outer contact of magnitude1 creates a proper affine E8 with arms1,2,5; magnitude sqrt(2) creates a proper affine B4 (one doubled terminal arm at a fork). Thus any remaining single contact has magnitude d>=phi>3/2, and its contribution is `8/d^2-2<14/9`.

When c joins the inner component to D2, classification gives either one leaf by a bond3 (the other A1 remains separate), or both leaves by bonds3. The inner contact is always a bond3 at its long endpoint. The respective outer data (q,t0) are (1/2,1/2) and (1,1), with h>=1+t0. Their contributions `2/h+(1-q)/h^2` are at most14/9 and1. A doubled bond to a single leaf would create a proper affine component beside the unused A1; a bond5 cannot occur in a large finite classical component. All other B/D joins are the complete short-end list and contribute at most13/9. Isolated and inner-only profiles contribute at most1 by the same proper affine B/D, A2 and small-star exclusions as in that proof. Hence every ordinary main-path profile has g<=14/9.

The only additional finite one-root extensions when an inner component is small are D6 to E7 and D7 to E8. B_j with j>=5 and D_j with j>=8 only admit the long-end bond3 finite extension. Here a small inner rank is r+P<=8; since P>=4 and n>=18 the opposite outer component has rank>=9. D6 to E7 is directly paid, since every pair in a finite component of order<=7 is paid. D7 means either r=2,P=5 or r=3,P=4. In the first case the marked a is the E8 branching vertex, so Vinberg's finite special-vertex lemma pays it ([original paper, §6.3, printed p81](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt)). The second case is treated below. Thus in every nonexceptional case

`det Schur/(h^2 k^2)>=1-28/(9P)>=2/9>0`.

### The main-path D3/P4 to E8 exception

Use an old D-type end whose fork center is0, fork leaves are L,L', and a is the next main-path node1. Let b be main-path node5. The new d has its single bond3 to L', so W+a+d has an E8 component. The arbitrary opposite contact is y=-b.d and k=y+1/2. At y=0 the old affine tree plus d contains a proper affine E8, so y>=1 and v=1/k<=2/3. Exact elimination gives alpha=1/4, A=1/2, D=5/4.

An outer-only c makes c,d belong to the same spherical component of W+c+d of order at most5; that corner is paid. The remaining c profiles and normalized determinant bounds are:

| c profile | reciprocal bound | normalized determinant | lower bound |
|---|---|---|---|
| isolated | u^2<=1/2 | 1-v-v^2/16-u^2/2 |1/18|
| inner-only bond3 | u<=1/2 |1-u/2-v-(2u-v)^2/16 |11/144|
| join through fork center | u<=1/2 |1-u/2-v-v^2/16 |1/18|
| join through one fork leaf | u<=2/5 |1-u/2-v-(u+v)^2/16 |14/225|

For the leaf join, x=0 creates an affine A4 cycle, so h=x+3/2>=5/2. The displayed expression uses the larger cross projection E=3/4 when c,d use the same leaf; the other leaf has E=1/4 and improves it. The inner-only and leaf-join quadratics are concave, so their minima on the stated rectangles occur at corners. Every bound is strictly positive.

## A terminal bond4 mark

Let a be that terminal node. An unpaid old pair has graph distance P>=6 from a to a main-path b; R=P and nu_a=1/sqrt(2). Deleting a leaves a large finite B or D diagram. Its complete finite-or-connected-affine one-node extensions are: isolated; bond3 to its long endpoint; bond4 to that endpoint; or bond3 to its next node. The last two are full affine extensions. This list follows by deleting the new node in the classical finite/affine diagrams: finite B/D permits only an added long-end bond3 at rank>=9; affine C/B adds a long-end bond4, and affine B/D adds a second-node bond3. No exceptional finite or affine diagram has the required large order. The corresponding (c0,l,h) are

`(2,0,x/sqrt(2))`, `(1,1,1+x/sqrt(2))`, `(0,sqrt(2),(x+2)/sqrt(2))`, `(0,2,1+x/sqrt(2))`.

In the first case x=0 is impossible by orthogonality to the full cusp, x=1 induces affine F4, and x=sqrt(2) induces affine C2. The conservative bound x>=sqrt(2) gives g<=2. In the second case x=0 induces the affine star(4,3,3), so x>=1 and again g<2. The third gives g=4/(x+2)<=2. In the fourth, x=0 gives a proper affine B4 with the doubled terminal at the end of its length-two arm. Thus x>=1 and g<=8-4sqrt(2)<5/2. The uniform bound is therefore5/2. The opposite main-path profile has g<=14/9: its near-terminal inner component is B_P of rank>=6, with no exceptional finite extension. Hence the determinant is at least `1-(5/2+14/9)/P>=35/108>0`.

## A fork-leaf mark

Let a be one leaf, with mate L and center0. An unpaid pair with a main-path b has graph distance P>=5; R=P+1 and nu_a=1/2. Deleting a leaves a large finite B or D diagram. Complete finite-or-connected-affine one-node extensions are: isolated; a bond3 to L (finite); a bond4 to L (full affine); or a bond3 to0 (full affine). The same complete large B/D deletion list applies. Their (c0,l,h) are

`(2,0,x/2)`, `(1,1,(x+1)/2)`, `(0,sqrt(2),(x+sqrt(2))/2)`, `(0,2,1+x/2)`.

For an isolated c, x=1 gives affine E8, x=sqrt(2) affine B4, and x=0 is impossible. Hence x>=phi>sqrt(2), giving g<=4. For the mate-leaf profile, x=0 gives affine E8 and x=1 gives affine A3; thus h>=(1+sqrt(2))/2 and g<=8-4sqrt(2)<4. The mate-leaf bond4 profile gives g=4sqrt(2)/(x+sqrt(2))<=4 even without excluding x=0. The center profile gives g<=4 for all x>=0. The opposite inner component is D_(P+1). Its D6 to E7 extension is paid. Apart from D7 to E8 below, the opposite profile has g<=14/9. Thus `1-(4+14/9)/R>=2/27>0`.

### The fork-leaf P6 to E8 exception

Here b is main-path node5, P=6, and d has a bond3 to the other fork leaf L. The same affine E8 exclusion gives k=y+1/2>=3/2. Literal elimination gives alpha=1/7,A=5/14,D=8/7. If c attaches by a bond3 to L, the opposite corner has a finite D8 component with c,d as its two fork leaves (the far outer end is separated by the deleted b). They lie in the same fork-invariant end segment, so the finite weight rule pays it. If c attaches by a bond4 to L, the opposite corner contains the proper affine star(4,3,3), impossible. For the center profile c→0, the opposite component is finite E8 with arms1,2,4; it must not be called affine. Here x=0 supplies a proper affine D4 centered0 with leaves a,c,L and main node1, and x=1 supplies affine A2(a,0,c). Therefore u=1/(1+x/2)<=2-sqrt(2). Its exact normalized determinant is `1-4u/7-5v/7-v^2/28`, at least `(36sqrt(2)-40)/63>0`. For isolated c put u=2/x<4/3 using x>=phi>3/2, and v=1/k<=2/3. Then

`det Schur/(h^2 k^2)=1-5v/7-v^2/28-2u^2/7>0`,

because the expression is zero at u=4/3,v=2/3 and strictly decreases in both positive variables.

## Two terminal marks

They must be at opposite ends. The effective resistance is R=n-1, for a pair of doubled terminals, a doubled terminal and a fork leaf, or two fork leaves. Their respective total g bounds are5,13/2,8. Since n>=18, all are strictly smaller than R. The local supports are separated, so the same factorized Schur formula applies without an E8 overlap.

## Status

The [exact checker](../scripts/fresh_simple_cusp_affine_tree_completion_20261017.py) and [saved data](fresh_simple_cusp_affine_tree_completion_20261017.json) bind the complete proof inputs. Generation passes 1,191 unpaid marked-pair contexts at orders18,19,23, 82,542 generic profile-pair formula probes, 14 literal exceptional E8 calculations, the scalar margins, and the additional affine B4/D4 positive-kernel witnesses. The finite checks verify explicit formulas; the Laplacian/Green derivation and finite/affine classification supply their all-rank scope. Generation and a fresh saved replay both terminate with PASS. AST and all 11 direct source hashes pass. Separate internal mathematical/source reviews by the root and construction agents pass the corrected four-profile terminal lists, cusp normalization, complete residual classification and exceptional Schur formulas. Neither review claims to duplicate the 82,542 producer probes; the construction agent separately recomputed the terminal and exceptional formulas. The separate [terminal and E8 audit](fresh_simple_cusp_affine_tree_terminal_audit_20261017.md) also passes generation and a fresh saved replay, including 48 independent physical normalization checks and symbolic exceptional determinants. That audit binds the final producer source/data hashes. This report and its producer source/data are now frozen.

The two additional terminal affine replacements and the finite D8/E8 descriptions in the fork-leaf exception are included explicitly. No dimension conclusion is asserted: the all-finite-corner transfer and final composition remain separate obligations.
