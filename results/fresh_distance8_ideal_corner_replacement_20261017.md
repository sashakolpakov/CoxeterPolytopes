# Replacing an ideal distance-nine quadrilateral payment

**Final, internally checked local theorem.** Let Q be an actual quadrilateral face of a finite-sided finite-volume hyperbolic Coxeter polytope, with no simplicity assumption at ideal vertices. Suppose an ideal corner has its two marked roots at graph distance9 in the varying irreducible affine factor F, of order greater than17. Then either a mark lies within14 steps of a feature of F, or an adjacent black pair of Q is joined by a path of length at most8 whose internal vertices are supporting white roots. A feature means a vertex of degree at least3 or an endpoint of a bond with label at least4.

There is no global exclusion of small affine diagrams. The theorem is compatible with arbitrarily many spectator factors and supporting whites. It is a local replacement lemma, not a new ambient dimension bound.

The [checker](../scripts/fresh_distance8_ideal_corner_replacement_20261017.py) and [certificate](fresh_distance8_ideal_corner_replacement_20261017.json) verify the exact weak-profile inequalities. Their all-rank classification and geometric interpretation are proved below. The corresponding affine graph rule has the separately proved [cost at most8m on a factor of order m](fresh_affine_distance8_feature_payment_20261017.md).

## Actual corner and spectator reduction

Write the face star as H=W+a+b+c+d, with opposite pairs {a,c} and {b,d}. All roots have square2. The white diagram W is positive definite; W plus any one black supports an actual edge and is also positive definite. Each adjacent corner is positive semidefinite. The full star has rank n+1 and corank1: its four side normals project to a spanning three-dimensional Lorentz space perpendicular to W.

At the chosen ideal corner W+a+b, the [factor-fiber lemma](fresh_cusp_factor_payment_framework_20261017.md) gives the varying affine component F containing a,b, plus a positive-definite spectator forest Z orthogonal to F. This statement concerns the actual face star, not the entire set of facets at a nonsimple cusp.

Assume for contradiction that there is no path of length at most8 from {a,c} to {b,d} with white interior. In particular b.c=a.d=c.d=0. Let nu be the positive affine kernel combination of F. It is a nonzero future null vector representing this actual ideal corner. Neither opposite side c nor d contains that corner, so

    h=-c.nu>0,  k=-d.nu>0.

The opposite contacts x=-a.c and y=-b.d may be zero, any finite Coxeter magnitude, parallel2 or dotted. Every nonzero magnitude is at least1. No upper label or distance cutoff is imposed on x,y.

Eliminate the whites of F first. The frozen [four-root determinant identity](fresh_simple_cusp_affine_a_quadrilateral_20261017.md) applies. Eliminating Z afterwards subtracts a positive-semidefinite matrix [[p,q],[q,r]] from the c,d block. By the [spectator identity](fresh_cusp_factor_payment_framework_20261017.md), this changes the normalized determinant by

    alpha*(p/h² - 2q/(h*k) + r/k²) >= 0,

where alpha>0. The heights are unchanged because Z is orthogonal to F. Thus a strictly positive determinant before the spectator subtraction contradicts the required corank1 afterwards as well.

## Affine cycles: complete weak profiles

Suppose F is affine A. Its two a-to-b arcs have edge lengths P,Q≥9, since a,b have graph distance9. Deleting a leaves A_(P+Q−1), of rank at least17. Its complete one-root positive-semidefinite extensions give the following c profiles: isolated; one simple bond at either endpoint; one doubled bond at either endpoint; one simple bond at the second node from either endpoint; or simple bonds at both endpoints. In ranks at least17, finite exceptional diagrams cannot occur. The connected finite extensions are A, B and D, and the only connected affine extension of the entire path is its simple cycle. Disconnected spectators do not enlarge this list, because the principal old path plus c is itself positive semidefinite. This is the same full finite/affine classification used in the [previous cycle proof](fresh_simple_cusp_affine_a_quadrilateral_20261017.md); here none of its small-affine exclusions is used.

For an attachment to an arc, let j=1 at an endpoint and j=2 at its next node. Let l=(l1,l2) be attachment magnitude times j on the two arcs. The norm constants, heights and weak bounds are:

| Profile | c0 | l | h | Minimum h |
|---|---:|---|---|---:|
|isolated|2|(0,0)|x|1|
|endpoint3 on arc i|1|unit vector i|x+1|1|
|endpoint4 on arc i|0|sqrt(2) times unit vector i|x+sqrt(2)|sqrt(2)|
|second-node3 on arc i|0|2 times unit vector i|x+1|1|
|both endpoints3|0|(1,1)|x+2|2|

The isolated minimum follows from h>0 and the native Coxeter gap. All other minima use only x≥0. Define analogous d0,m,k for d, and u=1/h, v=1/k.

There is one additional necessary length inequality. If c and d both attach on the same arc, at positions j_c,j_d from its opposite ends, their path through that arc has length P−j_c−j_d+2. It has white interior because j_c,j_d≤2 and P≥9. The no-short-path assumption therefore implies

    P >= max(9,7+j_c+j_d).

Apply the corresponding inequality to Q. If one attachment is absent, the original bound9 suffices. This refinement is essential for two second-node attachments: their common arc must have length at least11.

Put A=1/P, B=1/Q, L=u*l1+v*m1 and M=u*l2+v*m2. The normalized four-root determinant is exactly

    F=1-2*(A*L+B*M)
        -(A+B)*(c0*u²+d0*v²)-A*B*(L-M)².

For fixed u,v≥0, F is nonincreasing in A and B. For fixed A,B≥0 it is jointly concave in u,v: its quadratic terms are negative diagonal squares and one negative square. It therefore attains its minimum on the permitted height rectangle at a rectangle vertex. For each of the64 ordered profile pairs, the checker uses the largest permitted reciprocal arc lengths, then checks all four rectangle vertices in Q(sqrt(2)). All256 values are at least5/81. Hence F≥5/81>0 for every actual arc length and every allowed height. This is an arbitrary-rank, continuous-height inequality, not extrapolation from finitely many lengths.

Consequently an ideal affine-A corner of marked distance9 cannot occur in an unpaid quadrilateral. Additional spectator contacts only increase this positive determinant.

## Feature-free classical tree corridor

In order greater than17, the other affine factors are the classical B/C/D trees. Suppose both marks are farther than14 from every feature. They lie on the undecorated simple-bond spine, at distance P=9. Beyond each mark, the outer component is a finite B or D diagram attached at its long endpoint and has rank at least9. Indeed a fork feature is at least15 simple steps away, and a doubled bond's nearer endpoint is also at least15 steps away; retaining the entire end decoration only increases the outer order. Deleting either mark therefore leaves two large finite B/D components.

A complete one-root positive-semidefinite extension of these two old components has eight possible profiles: isolated; an endpoint3, endpoint4 or second-node3 attachment to either component; or simple attachments to both long endpoints. For a single component of rank at least9, its finite one-root extension only prolongs the long endpoint by3. Its affine one-root extensions add4 at that endpoint or3 at its second node, whenever the old type permits that profile. We safely allow both affine alternatives for either old type. If both old components are joined, their two remote decorations force the corresponding affine C/B/D diagram, joined by simple endpoint bonds. No other connection fits the finite/affine classification; an additional branch or internal doubled bond would violate it. In particular this argument explicitly retains the one-component affine extensions that were excluded in the earlier simple-cusp corridor proof.

Normalize nu to coefficient1 along the spine and at all relevant ports. For c, let the outer component be the one beyond a, and the inner component the one containing b. After eliminating the white roots of F, its scalar data and weak cost are:

| Profile | l | c0 | Minimum h | Maximum 2*l/h+c0/h² |
|---|---:|---:|---:|---:|
|isolated|0|2|1|2|
|outer endpoint3|0|1|1|1|
|outer endpoint4|0|0|sqrt(2)|0|
|outer second-node3|0|0|1|0|
|inner endpoint3|1|1|1|3|
|inner endpoint4|sqrt(2)|0|sqrt(2)|2|
|inner second-node3|2|0|1|4|
|both endpoints3|1|0|2|1|

These constants follow from exact all-rank inverse entries. The long endpoint and next node of either finite B_r or D_r have inverse diagonals1 and2. The [corridor inverse-column argument](fresh_simple_cusp_affine_tree_corridor_20261017.md) gives the endpoint column. For the second-node column, put1 at the long endpoint and2 on every subsequent spine node, with value sqrt(2) at the doubled B terminal or1 at each D fork leaf. Direct multiplication by the native norm-two Gram gives the unit vector at the second node, proving its inverse diagonal2 in every rank. On the middle A_(P−1) path, the corresponding entries are (P−1)/P and2(P−2)/P. The cross term for the separated inner contacts is l*m/P; their positions are at most2 from opposite ends and P=9, so no overlap exception occurs.

Thus B=l/P, C=c0+l²/P, A=m/P, D=d0+m²/P and E=l*m/P in the four-root identity. Cancellation gives

    det(S)/(h²*k²)
      = 1-[2*l/h+c0/h²+2*m/k+d0/k²]/P.

Each single-root cost is at most4. With P=9, the determinant is at least1/9>0. The checker verifies all eight cost bounds and64 ordered pairs exactly. The same estimate holds for every P≥9 under the same outer-rank hypothesis. Once again spectator elimination preserves the contradiction.

## Consequence and verification scope

The two branches prove the stated local theorem. In the tree case, failure of the feature alternative forces the corridor contradiction. In the cycle case there are no features, and the weak-profile contradiction applies directly. Every resulting short path selects an actual adjacent corner of Q; if that corner is ideal, the factor-fiber lemma places both marked roots and the path in its varying affine factor.

Generation and a fresh separate-process saved replay both terminated exit0 with `PASS_UNIVERSAL_IDEAL_DISTANCE9_CORNER_REPLACEMENT`. All16 direct source bindings and source syntax passed. The parent agent read the complete source and proof without finding a gap; that was a mathematical/source review rather than an independent arithmetic replay.

The [separate independent review](fresh_distance8_ideal_corner_replacement_20261017_review.md) reconstructed the algebra without importing producer arithmetic. It matched all256 cycle values, eight tree costs and64 pair gaps, and checked the full formal determinant, PSD correction and four symbolic all-rank inverse identities. Its full mathematical/source read also checked the actual factor-fiber geometry, positive heights, complete enlarged extension lists, arc-length refinement and long-tail reduction. That report records its own final bound replay and hashes.

The source rechecks the formal determinant identity and binds the actual-corner, spectator and finite/affine classification inputs. It does not enumerate new high-core diagrams, assert that a finite numerical rank range exhausts them, or provide an ambient dimension theorem. External refereeing and formalization remain separate.
