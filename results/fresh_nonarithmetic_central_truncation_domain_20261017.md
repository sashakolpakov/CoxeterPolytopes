# All-label obstruction to a three/five-degree central truncation patch

For every one of the twenty specified literal ring recipes, no norm-two
Coxeter wall inducing an actual S6 facet can have finite contact with all
four ring walls and make three consecutive ring corners finite. This
excludes the first center in the new three/five-degree two-center disk,
before imposing the second center, the remaining ambient wall gates, or a
common strict center. All finite contact labels and arbitrary real
coefficient fields are covered.

This is an exact obstruction for these twenty ring recipes. It is not an
obstruction to changing a ring parameter, a different patch, or walls
missing S6. The exact ring recipes were selected by an earlier numerical
search; their complete ambient native gates are not promoted by this
result. No new chamber or finite-volume example is asserted.

The [checker](../scripts/fresh_nonarithmetic_central_truncation_domain_20261017.py)
and [exact data](fresh_nonarithmetic_central_truncation_domain_20261017.json)
use the frozen [central-ring construction engine](fresh_nonarithmetic_central_ring_contact_gate_20261017.md).
This is a fresh computation of the enlarged high-label domain. The
low-label exclusions are inherited only after literal equality of every
saved recipe, projected Gram matrix, inverse pivot and ring-edge label.

## The changed topology

Write the old boundary cycle as a,b,c,d=[7,12,52,30] and the four ring
walls in cyclic order L,P,R,Q. Retain the eight peripheral triangles

`abL, bcL, cdP, daQ, cLP, dPR, dRQ, aQL`.

Use the central triangles

`TPR, TRQ, TQL, SLP, SLT, SPT`.

Here T meets all four ring walls, while S truncates the L-P-T corner.
The central degrees are five and three, rather than the four/four degrees
of the [previous two-center test](fresh_nonarithmetic_two_central_caps_20261017.md).
Every vertex link and edge incidence is checked: V=10,E=23,F=14 and
Euler characteristic one. The incidence counts are a=8,b=11,x=4,y=4,z=6.
Under the inherited genuine-boundary, simple, no-ideal disk hypotheses,
this would reduce intrinsic open ends by four. This combinatorial identity
alone does not establish geometric realization.

All four choices of the truncated consecutive ring corner are included.
The remaining three corners form a three-edge path on the four ring
walls. Each ring wall still meets T at one of these finite corners, so all
four T-ring contacts must be finite native Coxeter angles.

The previous one-center gate required all four corners to be finite and
therefore did not already establish the present obstruction.

## Complete all-label reduction

The retained support is S6=H4+A2. An actual new S6 facet is right to H4
and has one of the seven A2 contact profiles: zero, or a single label
3,4,5 at either A2 node. This is the inherited
[finite-extension classification](fresh_nonarithmetic_a2_escape_partition_20261017.md).
The four projected ring normals span the residual Lorentz four-space.
For prescribed native contacts, subtracting the A2 inverse-Gram products
gives the residual contacts c, and the unique norm test is

`c^T G^(-1) c = 2 - profile^T A2^(-1) profile`.

When all four labels lie in {2,3,4,5}, the frozen exact gate rejects all
7×4^4 assignments per ring: 35,840 in total. This norm test did not use
any finite-corner hypothesis, so it applies unchanged here.

Now suppose a T-ring bond has label m>=6. Its magnitude is at least
sqrt(3). In a finite norm-two Gram matrix, another nonright bond incident
to either endpoint would give a three-root determinant at most

`8 - 2(3+1) = 0`.

Thus the high bond is isolated within each required finite corner.
In particular T and that ring wall have zero A2 profiles; that ring wall
is right to its neighbors on the retained three-edge path; and T is right
to those neighbors. The checker constructs every eligible high vertex
and every possible high-vertex subset. An independent set in this path
has at most two vertices, so there are no cases with three or four high
contacts. No upper bound on a high label is imposed.

## Exactly one high contact

For each eligible high vertex, all other contacts lie in {2,3,4,5}, with
path-neighbor contacts forced to be right. There are 704 such three-contact
frames across the 80 ring/corner choices. Their frame Grams are
nonsingular with a one-dimensional positive orthogonal kernel. The saved
particular solution and kernel therefore give the complete norm-two
solution set. Every frame has exactly two real branches, for 1,408 total.

The factored projection formula computes the remaining native product
without multiplying expanded residual coordinates. If h is its negative,
then a finite label m>=6 requires sqrt(3)<=h<2. Exact outward rational
intervals give:

- 1,160 branches below sqrt(3);
- 207 branches at or above the parallel threshold two;
- 41 branches strictly inside the range.

Every one of the last 41 satisfies a strict consecutive-angle bracket

`2 cos(pi/m) < h < 2 cos(pi/(m+1))`.

The observed lower bracket labels are 6,7,8,9,10,11,12,14,15,29. Floating
point proposes the neighboring integer; acceptance uses only the strict
rational interval inequalities. Cosine thresholds use the frozen rational
Machin/alternating-Taylor bounds, with exact radicals for the special
labels. Consequently there is no untested label between the bracket
endpoints. The smallest certified gap is greater than 0.00002015.

The data also saves auxiliary three-corner principal minors: thirteen
branches have all three corners certified finite. Some remaining auxiliary
minors enclose zero and are deliberately left unclassified. None is used
to reject a branch: all 41 are excluded by the independent strict cosine
gap, so no zero-incidence decision is needed for the theorem.

## Exactly two high contacts

There are six eligible high-pair configurations. The other two T-ring
contacts are forced right, and T has zero white profile. Write the high
magnitudes as u,v. If a,b,c are the appropriate entries of the inverse
four-ring Gram, the norm condition is

`a u^2 + 2 b u v + c v^2 = 2`.

It suffices to consider the larger closed real square
`(u,v) in [sqrt(3),2]^2`, which includes every pair of finite high labels.
Both partial derivatives have a fixed strict sign on this entire square
in each case. The exact minimum is therefore bounded at the corresponding
corner. The six certified lower bounds are respectively greater than
3.4924,70.2111,10.4298,13.9672,11.8400,682.1901. Each is strictly greater
than two. Thus even arbitrary real high bonds in the allowed range fail;
no label enumeration or field assumption enters this step.

## Control, deletion costs and verification

The frozen control remains the
[64-wall chamber with 198 finite vertices](fresh_nonarithmetic_replace47_close12_20261017.md)
and its [complete intrinsic boundary](fresh_nonarithmetic_replace47_close12_boundary_20261017.md):
S6 has 25 finite vertices, seven exterior rays and fifteen open ends. Its
strict escape still proves infinite volume.

The selected ring family changes12 and deletes47,53. Its inherited loss
ledger contains the four old finite supports S6+[7,8,12],S6+[7,12,47],
S6+[8,12,52],S6+[12,47,52]. Wall53 had no old finite vertex. The prior
compact8-heptagon and12-quadrilateral are touched and are not presumed
preserved. Since the required first central wall is excluded, no changed
ambient finite-vertex census or complete new chamber boundary is claimed.

Generation and a fresh saved replay pass for all twenty rings, eighty
omitted-corner cases, 1,408 one-high branches and six two-high rectangles.
The checker asserts literal agreement with the inherited low-contact
certificate, source hashes, complete branch counts, every strict angular
gap, every monotonicity sign and every norm margin. The interval engine
and exact ring-construction algorithms are shared with the frozen parent;
this is not an independent second arithmetic implementation.
