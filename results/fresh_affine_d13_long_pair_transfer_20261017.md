# Universal long-pair transfer in the14-node affine D13 factor

**Completed internally checked local theorem; verification scope below.** Let an actual quadrilateral of a finite-sided finite-volume hyperbolic Coxeter polytope have an ideal corner whose varying affine factor is D13, of order14. Take one mark from each fork triple (the fork center and its two short leaves). There are nine pairs, at distances9,10 and11.

Each such quadrilateral is paid elsewhere by finite distance at most8 or a whole finite marked component of order at most17; ideal distance at most5 or a whole varying affine factor of order at most13; or two sealed finite half-coupons. The latter join a whole A1,A2,A4 or A5 component to a whole D11 or D13 component marked at its fork center or a short fork leaf. The [finite-eight coupon budget](fresh_finite8_a17_d11_coupon_budget_20261017.md) includes these types, with the small side enlarged to rank at most5.

This permits removal of the nine specified ideal pairs, reducing an otherwise all-pair D13 factor cost from91 to82. It is a universal local theorem, not a claim that all other pairs in every small factor are paid, or that all finite-eight overlap cases have been covered. No new ambient bound or finite-volume construction is asserted.

The [source](../scripts/fresh_affine_d13_long_pair_transfer_20261017.py) and [certificate](fresh_affine_d13_long_pair_transfer_20261017.json) recompute the14-node profile formulas, the newly non-strict cases and all72 extreme equality configurations. They preserve the frozen [D14 predecessor](fresh_affine_d14_long_pair_transfer_20261017.md).

## Actual star and normalization

Write H=W+a+b+c+d, with opposite pairs {a,c},{b,d}. All roots have square2. W and each actual edge diagram W plus one black are positive definite. Each corner is positive semidefinite. The full actual face star has corank1 in its Lorentz realization.

At the selected ideal corner, the [factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md) gives the varying affine factor F containing a,b, plus a positive-definite spectator diagram Z orthogonal to F. The positive affine kernel represents the actual cusp, and the opposite sides have strictly positive heights. No simple-cusp or arithmetic assumption is needed. A nonzero adjacent black product already pays at distance1. A shared spectator neighbor of c,d gives a paid path of length2. These alternatives can be removed from the unpaid analysis.

Number the simple spine0,...,9, with left fork leaves10,11 and right leaves12,13. The affine kernel nu has coefficient1 on the spine and1/2 on the leaves. The marks range over a in{0,10,11} and b in{9,12,13}. The effective resistance is

    R=9+2*[a is a leaf]+2*[b is a leaf].

The inherited all-rank Green identities apply at these ranks. The source also directly eliminates F minus the two marks, checking its normalized old2x2 block,96 individual end formulas and256 separated contact cross terms. For each end write its height h, moment l and norm constant c0. Its scalar cost is g=2*l/h+c0/h^2. The normalized determinant before spectators is

    F0=1-(g_c+g_d)/R.

Writing x,y for the native nonnegative contacts to Z, spectator elimination adds

    (x/h-y/k)^T G_Z^(-1)(x/h-y/k)/R >=0.

This is the same formally checked Schur identity as the predecessor; no small-affine exclusion from an earlier simple-cusp proof is imported.

## Complete end profiles

Deleting a short marked leaf leaves a large finite D13 diagram. Its complete one-root positive-semidefinite extensions are isolated, simple attachment at its long endpoint (the old mate), doubled attachment there, and simple attachment at the next node (the old fork center). An exceptional finite/affine extension cannot contain this D13. A multi-contact cycle is incompatible with the surviving remote fork in any finite/affine classification. For native opposite magnitude t, the exact table is unchanged from the predecessor:

| Profile | c0 | l | h | Cost bound |
|---|---:|---:|---:|---:|
|isolated|2|0|t/2|8, equality only t=1|
|mate3|1|1|(t+1)/2|8, equality only t=0|
|mate4|0|sqrt(2)|(t+sqrt(2))/2|4|
|center3|0|2|1+t/2|4|

The isolated height must be positive. Its non-extreme native magnitude is at least sqrt(2), giving cost at most4. The non-extreme mate3 has t>=1 and cost at most3. Thus the two extreme leaf profiles have height1/2 and are precisely a new simple leaf attached to either old fork leaf; every other leaf profile has cost at most4.

Deleting a marked fork center leaves D11 plus two isolated A1 roots. The large block has a long endpoint and its second node. A new root can be isolated from the large block, attach simply or doubly at the endpoint, or attach simply at the second node. The doubled-end and second-node alternatives produce a12-node affine component containing the new root and the far mark. The source verifies their positive kernels. Further acute contacts to that component vanish by its kernel, so the actual varying factor has order12 and is paid by the small13 rule.

If the new root attaches simply at the long endpoint, its further connections to the two old A1 roots are exactly: neither; one by3; one by4; or both by3. The large finite block excludes H/F/E alternatives, while an affine component's positive kernel forbids an extra acute contact. All these are the same complete finite/affine lists used in the D14 proof, with the large ranks now11 and13, still beyond every exceptional finite diagram.

For completeness, arbitrary contacts with only the two outer A1 roots are not replaced by a finite label search. If their magnitudes are u,v and the native opposite magnitude is t, then

    c0=2-(u^2+v^2)/2,   l=0,   h=t+(u+v)/2.

The actual edge diagram gives u^2+v^2<4. With no outer contact, positive height gives cost at most2. With one outer contact and t>0, the cost is at most2/3. With t=0 and a single magnitude u, the cost is8/u^2-2: its maximum6 occurs only at u=1, and otherwise the native gap gives cost at most2. With two nonzero outer contacts, h>=1 and c0<=1, so the cost is at most1. Thus the unique extreme center profile is a new simple leaf attached to either old outer leaf, native opposite contact zero, and height1/2.

All non-extreme center profiles have cost at most3. Equality3 occurs only for the simple inner-end attachment, with native opposite contact zero and height1. The other joined profiles have bounds14/9,4-2sqrt(2),1; an isolated center profile has bound2. These formulas cover every finite label and parallel or dotted native opposite contact. The source retains the14 exact scalar formula/derivative checks and recomputes the present pair margins.

## The new non-extreme zero cases

The resistance and bounds are:

| End types | R | Extreme sum | Maximum with a non-extreme end |
|---|---:|---:|---:|
|center / center|9|12|9|
|center / leaf, either order|11|14|11|
|leaf / leaf|13|16|12|

A strict inequality in the last column gives F0>0 and contradicts actual corank1 after the nonnegative spectator correction. Equality in the first two rows requires a center inner-end3 profile of cost3 at height1, and an extreme opposite end. Every other non-extreme combination has a strict gap.

For these zero cases, the spectator correction must vanish. Positive definiteness gives x/h=y/k. If either vector is nonzero they have a common support coordinate and the path c-e-d pays. Otherwise both vanish: every spectator is orthogonal to F,c,d.

At the adjacent corner obtained by deleting that marked center, the new simple inner-end root extends D11 to D12. Both marks lie in this whole finite D12 component. Orthogonality of every remaining spectator prevents enlargement or joining. It is paid by the small17 finite rule. The source checks all12 oriented configurations, their full signatures(14,1,1), positive-definite white blocks, all corners and the specified whole D12 unit payment. There is no assumption that a zero relaxed determinant alone is impossible.

## Extreme heights and the three-quarter spectator equality

It remains that both ends are extreme. Their heights are **first proved to be h=k=1/2** by the native profile analysis. Their cost sum is R+3. Thus the zero-determinant equation becomes

    -3/R + 4D/R=0,    D=(x-y)^T G_Z^(-1)(x-y),

and requires D=3/4. This is an unscaled native-contact quantity only after the height reduction.

A shared spectator neighbor is already paid. If both disjoint supports are nonempty, choose i in one and j in the other. If i,j are not adjacent, Gram-metric Cauchy gives D>=(x_i+y_j)^2/4>=1, impossible. If they are adjacent, the actual white-interior path c-i-j-d has length3 and pays under finite-eight or ideal-five. Consequently an unpaid extreme case has exactly one nonempty spectator support.

The frozen [strict spectator spectrum](fresh_small_spectator_difference_gap_20261017.md) applies: Z+c and Z+d are positive definite principal subdiagrams of actual edge diagrams. Its one-support D=3/4 equality forces one unit contact at an endpoint of a **whole simple A3 component**, with every other spectator contact zero. Hence all other supporting spectators are orthogonal to the entire equality core. The possible two-support A3 equality from the general spectrum was not silently discarded; it has the already paid length-three path just described.

## All72 sealed equality configurations

The core now consists of the14-node affine factor, two simple leaves c,d attached to old fork leaves, and a simple A3 chain attached at one endpoint to c or d. There are nine mark pairs, four choices of the old attachment leaves and two choices of the spectator side, giving72 configurations.

Each19-root Gram has signature(17,1,1); its15-root white block is positive definite. The original corner is affine F with finite spectators, and the other three corners are finite. The source reconstructs all216 finite corner diagrams and classifies their whole marked components. Their half-coupon types are exactly:

| Small component | Large component | Occurrences |
|---|---|---:|
|A5|D11|24|
|A2|D11|24|
|A4|D13|24|
|A1|D13|24|

Every large mark is a fork center or short fork leaf. The A5 case is why the safe budget must permit small rank5: deleting an old fork center retains an outer leaf, the new root and the three spectator-chain roots. When a marked fork leaf is deleted, the corresponding small component is A4 instead. These are whole finite components, not truncated fragments, because the spectator equality has sealed the entire core.

Exactly32 configurations have two halves and no unit;32 have a small17 unit and an additional half;8 have three small17 units. Thus every face receives at least one unit in total. The finite budget already admits all four types. This concludes the universal local transfer, independently of whether any particular core by itself defines a finite-volume chamber.

## Verification scope

Generation and a fresh separate-process saved replay both terminated exit0/PASS for96 end identities,256 cross identities, two affine12 kernels,14 native-cost rows, all12 non-extreme zero cases and72 extreme cores with216 finite corners. Syntax, all14 direct source hashes and all6 report links passed. The source reuses the frozen predecessor's profile formulas and cost functions but recomputes the present14-node Grams and all new cases; it is a new producer, not a claim of source-independent arithmetic review. A separate full mathematical/source read found no gap in the complete D11/D13 profile reuse, all-label outer bound, normalized zero-case sealing, height-first three-quarter equality, short-path two-support alternative or whole-component coupon roles. That read did not independently repeat the96/256/72 calculations or review a global composition. The all-label classification and spectator equality remain explicit analytic proof steps. Source, data and report are final/frozen. No complete finite-eight overlap theorem, ambient dimension improvement, external refereeing or formalization is claimed.
