# Transferring all nine long pairs of an affine D14 factor

**Completed local proof; verification status is recorded below.** Let an actual quadrilateral in a finite-sided finite-volume hyperbolic Coxeter polytope have an ideal corner whose varying affine factor is D14, of order15. Choose one marked root from each fork triple, where a fork triple consists of a fork center and its two short leaves. There are nine such pairs, at distances10,11,12.

Each such quadrilateral has another unit-paid corner under the retained rules

* finite: marked distance at most9, or the marks in one whole connected component of order at most17;
* ideal: marked distance at most6, or the varying factor of order at most14;

or it has two sealed finite corners, each eligible for a half-payment between a whole small component of rank at most3 and a whole D12 or D14 component, with the latter mark its fork center or a short fork leaf. Consequently these half-payments permit removal of all nine specified D14 ideal pairs. Additional retained payments only strengthen this conclusion.

This is a universal face-local theorem, not an enumeration of the63 inherited preferred-core contexts. It permits arbitrary supporting white roots and nonsimple ideal vertices. No total finite-vertex budget or ambient dimension claim is made here. The [checker](../scripts/fresh_affine_d14_long_pair_transfer_20261017.py) and [certificate](fresh_affine_d14_long_pair_transfer_20261017.json) verify the literal formulas and all equality configurations.

## Actual face setup

Write the star as H=W+a+b+c+d, with opposite pairs {a,c} and {b,d}. Normalize roots to square2. The white diagram W and every actual edge diagram W plus one black root are positive definite. Each adjacent corner is positive semidefinite, and the full face star has corank1 after retaining its actual Lorentz span.

At the selected ideal corner W+a+b, the [factor-fiber lemma](fresh_cusp_factor_payment_framework_20261017.md) gives the varying affine factor F containing a,b, plus a positive-definite spectator diagram Z orthogonal to F. The positive affine kernel of F represents that actual ideal point. The two opposite sides do not contain the point, so their heights h,k are strictly positive. This does not assume that the full ambient ideal vertex is simple.

If any adjacent black product is negative, that marked distance1 corner is already paid. Hence assume b.c=a.d=c.d=0. A common nonzero spectator contact of c,d similarly gives a paid path c-e-d of length2. No restriction is imposed on the opposite products a.c or b.d. Each is either zero or has magnitude at least1; if its magnitude is greater than1, it is at least sqrt(2). Parallel and dotted products are included.

Number F by the simple spine0--1--...--10, left leaves11,12 attached to0, and right leaves13,14 attached to10. Normalize the cusp coefficients to1 at all spine nodes and1/2 at the four leaves. The resulting scaled Gram is a graph Laplacian: spine conductance1 and fork-edge conductance1/2. For a in {0,11,12} and b in {10,13,14}, the resistance is

    R=10+2*(a is a leaf)+2*(b is a leaf).

Eliminating F minus a,b leaves R^(-1)[[1,-1],[-1,1]] on the scaled marked roots. For a replacement near the a end, write its residual square as c0+l²/R and its cross contact with the scaled b as -l/R. The complete profiles below have separated left/right supports, so the projection cross term is lm/R. The [weighted-tree Schur identity](fresh_simple_cusp_affine_tree_completion_20261017.md) therefore gives

    det(S)/(h² k²) = 1-(g_c+g_d)/R,
    g_c=2l/h+c0/h².

Only its normalization, finite/affine classification and Green identities are used; none of that earlier report's global small-affine exclusions is imported. This source recomputes the end data for all nine present marked pairs and the separated projection term for every ordered profile pair, and rechecks the formal determinant identity.

Eliminating Z subtracts a positive-semidefinite matrix from the c,d block and adds

    (1/R)*(x/h-y/k)^T G_Z^(-1)(x/h-y/k)

to the normalized determinant, where x,y are the nonnegative native contact vectors of c,d with Z. The heights are unchanged. Thus a strictly positive spectator-free determinant is already impossible in an actual quadrilateral.

## Complete fork-leaf profiles

Suppose a is a short fork leaf, with mate L and center0. Deleting a leaves finite D14, with long endpoint L and second node0. Its complete one-root positive-semidefinite extensions at this rank are: an isolated root; a simple bond to L; a bond4 to L; or a simple bond to0. Finite D can extend only at its long end; the last two possibilities are the connected affine B/D extensions. A bond5 or an exceptional finite/affine extension cannot contain the old D14. Multiple contacts would create a cycle in a graph already containing the remote fork, which is not a finite/affine type. This is the same complete large-D extension list used in the [weak-profile proof](fresh_distance8_ideal_corner_replacement_20261017.md).

Put x=-a.c. The exact data are:

| Contacts of c in F minus a | c0 | l | h | Cost bound |
|---|---:|---:|---|---|
|none|2|0|x/2|8 at x=1; at most4 otherwise|
|L by3|1|1|(x+1)/2|8 at x=0; at most3 otherwise|
|L by4|0|sqrt(2)|(x+sqrt(2))/2|at most4|
|0 by3|0|2|1+x/2|at most4|

In the isolated profile h>0 implies x>0. The native gap gives either x=1 or x>=sqrt(2). In the mate3 profile, either x=0 or x>=1. The costs decrease with x, and the displayed bounds follow at these exact endpoints. Thus

    leaf end: g<=8, and every nonextreme profile has g<=4.

The two extreme profiles both say that c is a simple new leaf attached to one of the two old fork leaves: either the selected leaf a or its mate L. In both cases h=1/2. No small affine diagram is being forbidden to obtain these bounds.

## Complete fork-center profiles, including all finite labels

Suppose a=0. Deleting a leaves a large finite D12 block B on nodes1,...,10,13,14, together with two isolated A1 leaves11,12. Its long endpoint is1 and its second node is2.

First consider the principal diagram B+c. By the complete large-D extension list, c is isolated from B, is joined to1 by3 or4, or is joined to2 by3. The latter two give connected affine components of order13. A positive affine kernel forces every other acute root of the actual corner to be orthogonal to that component. It is therefore the entire varying factor containing the marks b,c, and is already paid by the retained whole-factor cutoff14. The source checks both literal positive kernels.

If c is joined to1 by3, the component is finite D13. Adding its possible contacts with the two old A1 leaves gives precisely: neither leaf; one leaf by3; one leaf by4; or both leaves by3. This is again a direct finite/affine classification: the large component excludes H/F/E, while an extra doubled terminal gives an affine B component whose positive kernel forbids another contact. The both-simple case is affine D14 of order15 and is retained in the analysis. The contact x=-a.c remains arbitrary because a is absent from this adjacent corner.

The four unpaid large-block profiles have data

| Additional outer contacts | c0 | l | h | Maximum cost |
|---|---:|---:|---|---:|
|none|1|1|x+1|3|
|one leaf by3|1/2|1|x+3/2|14/9|
|one leaf by4|0|1|x+1+sqrt(2)/2|4-2sqrt(2)|
|both leaves by3|0|1|x+2|1|

It remains to allow arbitrary contacts with the two A1 leaves while c is orthogonal to B. Write their magnitudes as u,v>=0. Exact elimination gives

    c0=2-(u²+v²)/2, l=0, h=x+(u+v)/2.

The actual edge W+c is positive definite, so (u²+v²)/2<2. This strict edge condition is needed; the corner condition alone would not suffice. The following bounds require no label cutoff:

* With no outer contact, h=x>0, hence x>=1 and g<=2.
* With one outer contact and x>=1, c0<=3/2 and h>=3/2, hence g<=2/3.
* With one outer contact and x=0, g=8/u²-2. It is6 only at u=1; otherwise u>=sqrt(2) gives g<=2.
* With two nonzero contacts, u,v>=1, so h>=1 and c0<=1, giving g<=1.

Consequently

    center end: g<=6, and every nonextreme unpaid profile has g<=3.

The unique extreme is a simple new leaf c attached to either old outer fork leaf, with x=0. Again h=1/2. The case of an arbitrarily large opposite product and every admissible finite outer label has been covered.

## Strict nonextreme gaps and the discrete equality

For center/center, center/leaf and leaf/leaf ends, the resistances are10,12,14. If at least one end is nonextreme, the cost sums are at most9,11,12, respectively. Thus the normalized determinant is at least1/10,1/12,1/7. Spectators cannot cancel these positive gaps.

Both ends must therefore be extreme. The old factor now has a simple new leaf c on one old left fork leaf and a simple new leaf d on one old right fork leaf. Both heights have been proved to equal1/2, and the total cost equals R+2. The determinant equation, including arbitrary spectators, is

    0=-2/R+(4/R)*D,
    D=(x-y)^T G_Z^(-1)(x-y).

Hence D=1/2. The height reduction precedes use of any unscaled contact gap.

If x,y have a common nonzero coordinate, the already-paid c-e-d path resolves the face. Otherwise their supports are disjoint. The [positive-definite spectator gap](fresh_ideal_distance6_mixed_replacement_20261017.md) states:

* If both supports are nonempty, Cauchy with e_i-e_j gives D >= (x_i+y_j)²/(4-2G_ij)>1/2.
* If exactly one is nonempty, D>=1/2, with equality only for one unit contact to an isolated A1 node of Z.
* If both vanish, D=0.

For completeness, in the second case equality in Cauchy with e_i forces x_i=1 and x=(1/2)G_Z e_i. Nonnegative x and nonpositive off-diagonal Gram entries force i to be isolated and every other coordinate zero. Thus the equality conclusion concerns the entire spectator diagram, including arbitrarily long indirect tails.

Exactly one extra isolated supporting node e therefore attaches simply to c or d. All other spectators are orthogonal to the whole18-root equality core. This is precisely the zero-residual extra wall allowed by the smaller17-root Lorentz core; it is not silently excluded by a strict inverse-diagonal test.

## Whole finite components and the two half-payments

The equality core is completely specified: the15-node old F, the two simple leaves c,d attached to old fork leaves, and one simple e attached to c or d. The checker treats all nine mark pairs, all four choices of the two attachment leaves and both choices of the extra-spectator side, giving72 cases. In every case the full core has signature(16,1,1), its supporting white block is positive definite, its original corner is affine F plus an isolated A1, and the other three corners are finite.

Deleting a selected fork center leaves a D12 component containing the opposite mark, while the new local marked component is A2 or A3. Deleting a selected fork leaf leaves D14 at the other end and a local A1 or A2 if c attaches to that selected leaf. The large mark is always its fork center or a short fork leaf. These give the two half-coupon corners.

If a selected fork leaf instead has the new attachment on its mate, an adjacent finite corner has both marks in one connected component of order at most17 and already contributes a unit. The literal cases give32 configurations with exactly two halves and no old unit;32 with an old unit and an additional half; and8 with three old units. Every configuration receives at least1.

These are whole actual finite components. The spectator equality has proved orthogonality of every remaining supporting root to the equality core; no undisclosed supporting tail can enlarge or join the coupon components. At a finite ambient vertex the star already contains all incident facets. At the original nonsimple ideal vertex, the factor-fiber theorem is what identifies F as the full varying factor; no compact-style triangle transitivity is assumed.

The nine long pairs may therefore all be removed from the affine-D14 local rule after adding the stated half-coupons. This reduces that factor's all-pair cost from105 to96 when all other pairs are retained. The cost of the new finite coupons and the resulting global comparison are separate obligations.

## Verification scope

Generation passes96 exact end-data calculations,256 separated-support cross identities, two positive affine13 kernels,14 scalar endpoint evaluations, four strict end-type margins and all72 equality configurations with216 finite corners. The determinant and spectator correction are checked as formal identities. The profile classification, continuous label bounds and spectator equality are analytic arguments; the finite evaluations are not a label cutoff or an enumeration of ambient chambers.

Final generation22597 and a fresh separate-process saved replay both terminated exit0/PASS. The affine component uses the standard rank label affineD14 for its15 nodes; finite D labels also denote rank. Syntax, all14 direct source hashes and all6 local links passed. A separate full mathematical/source read found no gap in the complete center/leaf extension lists, all-label outer-contact bounds, strict native gaps, the height1/2 reduction before the spectator gap, or whole-component sealing. That review did not independently repeat the96/256 formula calculations or72 equality configurations. Source, data and report are final/frozen. No external refereeing or proof-assistant formalization is claimed.
