# An affine-A corner of a quadrilateral is paid in every rank

Let P be a finite-volume hyperbolic Coxeter n-polytope, with every vertex
simple and n>=14. If an actual quadrilateral has a corner whose full
n-node diagram is affine A_(n-1), then the
[proposed mixed weight](fresh_simple_cusp_vinberg_weight_frontier_20261017.md)
pays that quadrilateral: either the two marked roots at this corner have
cycle distance at most3, or the marked roots at the opposite corner are
adjacent. This resolves the complete affine-A-cycle branch, with no bound
on labels or arc lengths. Other affine corner types remain outside this
result; **no new ambient dimension bound is asserted**.

The [checker](../scripts/fresh_simple_cusp_affine_a_quadrilateral_20261017.py)
and [data](fresh_simple_cusp_affine_a_quadrilateral_20261017.json) use exact
Fraction arithmetic in Q(sqrt(2)), a formal determinant expansion, and256
rectangle-corner checks. The argument below explains why these finite
checks certify all ranks rather than a bounded search.

## Geometry and complete replacement profiles

Use roots of norm2. Write a,b for the two marked roots of the affine
corner, c opposite a, and d opposite b in the quadrilateral. If a,b have
cycle distance<=3, the ideal weight already pays. Otherwise the two cycle
arcs have edge lengths P,Q>=4, P+Q=n. Deleting a,b leaves the white diagram
W=A_(P-1)+A_(Q-1). The roots a,b attach to the first and last ends,
respectively, of both paths. In particular a.b=0.

The actual adjacent corner W+b+c is either spherical or, if ideal, a
connected affine diagram on n nodes. The old diagram W+b is the single
long path A_(n-1). Consequently its complete one-root extension list is:

| Contacts of c with W+b | Result |
|---|---|
| none | A_(n-1)+A1 |
| one bond3 at either endpoint | A_n |
| one bond4 at either endpoint | B_n |
| one bond3 at either next-to-end node | D_n |
| bonds3 at both endpoints | affine A_(n-1) |

This is a classification argument in every n>=14. A new root touching
the connected old path produces an irreducible diagram of order n.
Exceptional finite and affine types are too small. Among the classical
finite and affine diagrams, deleting one node to leave the entire path
A_(n-1) gives precisely the displayed possibilities. In particular the
other three large affine types do not have such a deletion. The finite
and affine classification is the standard primary input recalled in
[Vinberg, section2](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

Here the endpoints and next-to-end nodes lie near a, not b, because both
arc lengths are at least4. Thus c.b=0 in every profile. The analogous
list for d uses the last ends near b and gives d.a=0. The isolated A1
profile is explicitly retained; its height comes from its omitted-root
contact with a or b.

If c.d<0, the opposite corner has adjacent marked roots, so it is paid
by either the finite or ideal rule. It remains to contradict c.d=0.
Write x=-c.a>=0 and y=-d.b>=0. Each is an arbitrary Coxeter contact:
zero, 2cos(pi/m) for m>=3, parallel2, or dotted>2. No arithmetic field
restriction is imposed on x,y.

## Small affine diagrams force uniform height bounds

The [frozen global exclusion](fresh_simple_cusp_affine_exclusion_20261017.md)
for this all-simple class forbids every connected affine subdiagram of
order below n. We use only affine A2,A3,B3,E8, of orders3,4,4,9.
The checker supplies a strictly positive kernel and positive-definite
one-node deletions for each of these four literal diagrams.

A simple leaf attached anywhere on the old n-cycle creates a forbidden
affine E8: take the leaf arm of length1 and the two cycle arms of lengths2
and5. These nine vertices induce precisely the tree with arms(1,2,5).
The unused complementary cycle arc has length n-7>=7, so its two selected
ends have no extra edge. A positive kernel has center6 and successive arm
coefficients (3), (4,2), and (5,4,3,2,1).

Let w be c's vector of nonnegative contact magnitudes with W, and put
h=x+sum(w). The complete low-contact exclusions are:

| Profile | Excluded x | Reason | Sufficient remaining bound for h |
|---|---|---|---|
| isolated |0,1,sqrt(2)| orthogonality to the full cusp; leaf E8; four-node B3 star | h>sqrt(2) |
| endpoint3 |0,1| leaf E8; A2 triangle | h>=1+sqrt(2) |
| next-to-end3 |0,1| leaf E8; A3 cycle | h>=1+sqrt(2) |
| endpoint4 |0| B3 star at the old endpoint | h>=1+sqrt(2) |
| both endpoints3 |0,1| A3 cycle; A2 triangle | h>=2+sqrt(2) |

For isolated x=sqrt(2), the B3 star has center a, two simple cycle
neighbors and the double-bond new leaf. For endpoint4 and x=0, its center
is the first white root, with neighbors a, the next white root and c.
The A2 and A3 diagrams in the other rows are literal triangles or
four-cycles on the two new/marked roots and their stated white neighbors.
The Coxeter alphabet has no value strictly between0 and1 or strictly
between1 and sqrt(2). After excluding isolated x=1,sqrt(2), the next
finite value exceeds sqrt(2); parallel and dotted contacts are larger.

For completeness, isolated x=0 cannot be a new spacelike root: it would
be orthogonal to every normal at the full simple cusp. Those n normals
are linearly independent, span the null hyperplane nu-perp, and have
positive null sum nu. Their orthogonal complement is R nu. Thus a root
orthogonal to all of them would be null. This also proves h>0. All other
profiles have h>0 directly. The same arguments apply to d and its height
k=y+sum(z).

## Exact Schur identity

Set a0=1/P, b0=1/Q, alpha=a0+b0. The norm2 A_(P-1) inverse is

    J_ij=min(i,j)*(P-max(i,j))/P,  1<=i,j<=P-1.

After eliminating W, the residual a,b Gram is
alpha[[1,-1],[-1,1]]. For c define

    B=v_b^t J w,  C=2-w^t Jw,

and for d define A=v_a^t Jz, D=2-z^t Jz. Finally E=w^t Jz.
Here v_a,v_b are the nonnegative old endpoint contact vectors and J is
the inverse Gram of all W. Since the first and last inverse rows sum
entrywise to1, the residual a.c product is B-h and the residual b.d
product is A-k. Under the remaining assumption c.d=0, the full residual
matrix in order a,b,c,d is

    S = [[ alpha, -alpha, B-h,   -A  ],
         [-alpha,  alpha, -B,    A-k ],
         [ B-h,    -B,     C,    -E  ],
         [-A,      A-k,   -E,     D  ]].

The actual quadrilateral star has n+2 roots in Lorentz dimension n+1,
rank n+1, negative index1 and corank1. In particular det S=0. The formal
24-permutation determinant check gives

    det S=(hk-kB-hA)^2-alpha*(C*k^2+D*h^2+2E*h*k).             (1)

Only this necessary determinant condition will be needed. No Gram
matrix with the wrong rank or signature is promoted to an actual face.

For each profile introduce its two endpoint moments l=(l1,l2) and c0:

| Profile | l | c0 |
|---|---|---:|
| isolated |(0,0)|2|
| endpoint3 on arc i |unit vector e_i|1|
| endpoint4 on arc i |sqrt(2)e_i|0|
| next-to-end3 on arc i |2e_i|0|
| both endpoints3 |(1,1)|0|

Then B=a0*l1+b0*l2 and C=c0+a0*l1^2+b0*l2^2. For the d-profile
write m,d0; then A=a0*m1+b0*m2, D=d0+a0*m1^2+b0*m2^2 and
E=a0*l1*m1+b0*l2*m2. The latter formula includes the boundary P=4:
the two next-to-end positions can coincide, and the same inverse
formula still gives4/P.

Put u=1/h>0, v=1/k>0 and L=u*l+v*m. Dividing(1) by h^2*k^2 and
substituting the displayed expressions yields

    F=1-2(a0*L1+b0*L2)
        -(a0+b0)*(c0*u^2+d0*v^2)-a0*b0*(L1-L2)^2 = 0.      (2)

## The all-rank strict inequality

Use the conservative closed bounds

    u<=1/sqrt(2)                    for isolated,
    u<=sqrt(2)-1                    for endpoint3/4 or next-to-end3,
    u<=1-sqrt(2)/2                  for both endpoints3,

and likewise for v. For fixed a0,b0>0, F is jointly concave in(u,v):
its quadratic part is the sum of a negative diagonal quadratic and a
negative square. Therefore its minimum on each closed height rectangle
is attained at a rectangle corner.

At the enlarged reciprocal-arc values a0=b0=1/4, the exact checker tests
all8*8*4=256 corners in Q(sqrt(2)). Every value is nonnegative. Exactly
three values are zero: both isolated profiles at their upper bounds,
and the two choices of both next-to-end3 profiles on the same arc at
their upper bounds. These are the only finite algebraic checks needed;
concavity covers every real height in every rectangle.

Finally, for u,v>0 both partial derivatives with respect to a0,b0 are
strictly negative. For example

    dF/da0=-2*L1-(c0*u^2+d0*v^2)-b0*(L1-L2)^2 < 0.

All terms are nonpositive. If L1=0 and the middle term vanishes, some
nonzero moment remains in L2 and the last term is negative; if both
moments vanish, both profiles are isolated and the middle term is
negative. The other derivative is identical with the two arcs exchanged.

Actual P,Q>=4 and P+Q>=14 give a0,b0<=1/4 with at least one strict.
Hence F(actual a0,b0;u,v)>F(1/4,1/4;u,v)>=0, contradicting(2).
Thus c.d cannot vanish. Its negative native product makes the opposite
marked pair adjacent, proving the asserted quadrilateral payment.

## Verification and scope

Generation and a fresh saved replay pass all256 exact rectangle checks,
all four positive-kernel affine witnesses and their one-deletion tests,
and the formal determinant identity. Syntax, all eight direct source
hashes and all four local report links pass. An independent internal
mathematical/source read finds no gap in the full extension list,
low-contact exclusions, height bounds, boundary P=4 cross term,
concavity or strict reciprocal-arc monotonicity. That read is distinct
from an independent replay of the256 arithmetic values.
The finite/affine classification and the inherited actual-cusp affine
exclusion are explicit proof inputs. No unbounded rank or label is
replaced by a numerical cutoff. The original mixed weight interface and
all earlier61/127/163 results remain frozen and unchanged.
