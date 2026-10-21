# Nonisolated D19/D20 quad continuations: finite opposite products excluded

19 October 2026. **New written research proof, awaiting independent
review.** No calculation program, enumeration, optimizer, Lean process
or TeX build was run. Earlier reports and frozen sources are unchanged.

This note completes one branch of the remaining fork/whisker problem.
It proves that the three possible nonisolated continuations of an
unpaid adjacent corner cannot have a finite ambient product with the
opposite old whisker. The dotted and parallel ambient cases remain;
there is no new refund or dimension bound here.

## 1. The actual pattern and its finite continuation inventory

Let n be 19 or 20. Let a finite vertex of a compact arithmetic
simplest-type hyperbolic Coxeter polytope have whole diagram D_n.
Write a for its fork, b,e for its two short whiskers, and

    a--l1--l2-- ... --lL,     L=n-3,

for its long arm. Every old bond has label 3. Put N=n-2, so L=N-1.
All roots below are actual outward wall normals, normalized to square 1.

Consider an actual quadrilateral whose successive side roots are
a,b,c,d, with opposite pairs (a,c),(b,d). Its white support is

    W={e,l1,...,lL},     W=A1+A_L.

Assume the old Vinberg weight at the adjacent corner (a,d) is zero.
This is a necessary condition for a quadrilateral to require its entire
unit of payment from the selected old fork/whisker angle (a,b).
No assertion that every quadrilateral is unpaid is made.

The corner W+a+d is finite. Before adjoining d, W+a is the path

    e--a--l1-- ... --lL

of type A_(n-1). Finite classification in these ranks permits d to be
isolated, to attach at an endpoint with bond 3 or 4, or to attach at a
next-to-end node with bond 3. Two contacts would create a cycle.
Larger exceptional finite diagrams cannot contain this A18/A19 path.
The near-end choices have positive weight at (a,d): attachment at a
makes the pair adjacent, and attachment at e puts it in an end segment
of order 3. Thus a zero-weight pair leaves precisely:

|Type|Contact of d with W+a|Bond|
|---|---|---:|
|I|None|None|
|S|lL|3|
|B|lL|4|
|D|l_(L-1)|3|

In particular a.d=e.d=0 in all four cases. These are a complete
inventory of *possible* zero-weight continuations. We need only the
necessity of this list, rather than a converse assertion about weights.
The original rule features used are adjacency positivity and positivity
for a pair in a classical end segment of order at most 6; see the
[reviewed rank-four note, Section 1](compact29_rank4_surviving_end_refund_20261018.md)
and Vinberg's Section 6, printed pp.80--81.

Type I is the isolated-d branch handled in the
[separate D19/D20 classification](arithmetic_dn_fork_whisker_quad_attachment_classification_20261019.md).
Here we consider S,B,D and put

    b.d=-beta.

The acute-wall condition gives beta>=0. In fact beta>0. Project off
the positive-definite white space W. The root b is unchanged and has
square 1; the root d projects to a spacelike vector of square u>0.
The projected sides b,d are opposite sides of the acute compact
quadrilateral. A positive-definite two-side Gram matrix would give an
actual common corner by the two-normal M-matrix projection argument.
A positive null combination in the semidefinite case would give a
feasible ideal ray. Thus

    beta^2>u>0.

For S,B,D respectively, direct A_L inversion gives

    u=(N+1)/(2N),     1/N,     2/N.

We do not infer beta>=1 in the ambient Gram matrix. The purpose of
this note is to exclude the finite possibilities

    beta=cos(pi/m),     integer m>=3.

## 2. Exact determinant costs and arithmetic normalization

Let G be the old unit-root D_n Gram matrix. Its relevant inverse
entries, obtainable from the standard D_n fundamental weights, are

    (G^-1)_bb=n/2,
    (G^-1)_(lL,lL)=2,       (G^-1)_(b,lL)=1,
    (G^-1)_(l_(L-1),l_(L-1))=4,
    (G^-1)_(b,l_(L-1))=2.

For the actual new root d let q=d_D G^-1 d_D^T. The determinant ratio
of the old roots plus d is 1-q. Its three physical formulas are

    S: q=n*beta^2/2 + beta + 1/2,
    B: q=n*beta^2/2 + sqrt(2)*beta + 1,
    D: q=n*beta^2/2 + 2*beta + 1.                 (1)

The old D_n roots span a positive n-space in Lorentz signature (n,1).
Their finite vertex lies strictly inside the new wall: a compact simple
vertex already has its n incident walls. Hence d has a nonzero
component in their timelike perpendicular line, and q>1. In particular
the n+1 old-plus-new roots are independent.

At each nonidentity embedding of the arithmetic ground field K,
the corresponding independent Gram block is positive definite.
Consequently its determinant ratio satisfies

    sigma(q)<1.                                  (2)

As usual these are determinant-ratio statements after positive
diagonal normalization of the arithmetic roots. They do not assume
that every unit-normal Gram entry belongs to K.

There is, however, a closed cycle consisting of the b.d edge, the
old path from b to the indicated long-arm node, and the new contact
from that node to d. In S and D, all contacts on that cycle other than
b.d have label 3. Its norm-two cyclic product is, up to sign, 2*beta.
Thus

    B0:=2*beta is an algebraic integer in K.        (3)

In type B exactly one other cycle bond has label 4. The analogous
conclusion is

    B1:=2*sqrt(2)*beta is an algebraic integer in K. (4)

The q formulas in (1) can therefore be conjugated as polynomials
in B0 or B1. No choice of the sign of a square root at an embedding
is being imposed.

## 3. Type S: a simple far-end attachment

Write B0=2*beta. From (1)--(2), every nonphysical conjugate satisfies

    (-2-2*sqrt(n+1))/n < sigma(B0)
        < (-2+2*sqrt(n+1))/n.

For n=19,20 this implies |sigma(B0)|<3/5. Indeed
2*(1+sqrt(n+1))/n<3/5 follows by squaring from n>160/9.
If beta is finite, then 1<=B0<2 physically. If its degree were at
least 3, its nonzero integer norm would have absolute value less than

    2*(3/5)^2=18/25<1.

Hence [Q(B0):Q]<=2. The degree formula for 2*cos(pi/m), m>=3,
now leaves only m=3,4,5,6.

For m=4 or 6 the other conjugate beta' is at most -sqrt(2)/2.
Using beta'^2>=1/2 and beta'>-1 gives

    q' > n/4-1/2>1,

contrary to (2). For m=5 the other conjugate is (1-sqrt(5))/4, giving

    q'=[3n+12-(n+4)*sqrt(5)]/16>1.

The last inequality is equivalent to
3n-4>(n+4)*sqrt(5), and its squared difference is
4*(n^2-16n-16)>0 for n=19,20 (already for n>=17).

For m=3, B0=1 and q is rational and greater than 1. Any nonidentity
embedding of K would fix q, contradicting (2); hence K=Q. A compact
simplest-type chamber in these dimensions cannot have rational ground
field, by rational isotropy of an indefinite rational quadratic form
of dimension at least 5. This excludes the last case.

## 4. Type D: a simple far next-to-end attachment

Again write B0=2*beta. Now (1)--(2) give

    sigma(B0)*(n*sigma(B0)/8+1)<0,

and therefore

    -8/n<sigma(B0)<0.

For finite beta, 1<=B0<2. Any degree at least 2 would give a nonzero
integer norm of absolute value less than

    2*(8/n)<=16/19<1.

Thus B0 is a rational integer, necessarily 1, so beta=1/2 and m=3.
The rational q>1 again forces K=Q and contradicts compactness in
these dimensions. Type D has no finite beta.

## 5. Type B: a far-end attachment of label 4

Put B1=2*sqrt(2)*beta. Formula (1) becomes

    q=n*B1^2/16+B1/2+1.

Consequently every nonphysical conjugate satisfies

    -8/n<sigma(B1)<0.                            (5)

For finite beta, sqrt(2)<=B1<2*sqrt(2)<3. Degree at least 3 would
give a nonzero integer norm of absolute value less than

    3*(8/19)^2=192/361<1.

In degree 2 its norm is a negative integer with absolute value below
24/19<2, hence is -1. Its other conjugate is -1/B1, so its integer
trace is B1-1/B1, lying strictly between 0 and 3. It is therefore
1 or 2.

Trace 1 gives B1=(1+sqrt(5))/2. Since sqrt(2)<B1<2, the value
beta=B1/(2*sqrt(2)) lies strictly between 1/2 and cos(pi/4).
There is no cos(pi/m), integer m>=3, in that open interval.

Trace 2 gives B1=1+sqrt(2), and hence

    beta=(2+sqrt(2))/4.

This lies strictly between cos(pi/5)=(1+sqrt(5))/4 and
cos(pi/6)=sqrt(3)/2. The first comparison follows from
1+sqrt(2)>sqrt(5); the second from sqrt(2)<3/2. Again no integral
finite Coxeter label exists in that interval.

In degree 1 the only integer in [sqrt(2),2*sqrt(2)) is 2. It gives
beta=1/sqrt(2) and rational q=n/4+2>1. As before, (2) forces K=Q,
which is incompatible with compactness. This completes type B.

## 6. Exact consequence and the remaining branch

**Local conclusion.** At n=19 or 20, an actual fork/whisker
quadrilateral whose adjacent corner (a,d) has zero original Vinberg
weight has only the four d-continuations in Section 1. The isolated
continuation is ruled out under the additional unpaid-corner
orthogonalities c.b=c.d=0 by the separately reviewed classification.
For each of the three nonisolated continuations, a finite ambient
product b.d is impossible, independently of every contact of c.
Thus every surviving nonisolated case has

    b.d<=-1.

This includes the parallel ambient value; it has not been removed.
The projected opposite sides are still strictly ultraparallel because
u<1 in all three cases.

The next useful geometric input is that the actual opposite corner
W+c+d is finite. It constrains the *joint* contacts of c and d, beyond
the separate W+c inventory. The full quad Schur relation may then
bound the remaining dotted parameter. Neither step is proved here.
In particular this note does not yet justify reducing a fork/whisker
weight by one, does not verify all triangle or quadrilateral consumers,
and does not change the compact arithmetic ceiling 20.

The inherited arithmetic inputs are Vinberg cyclic integrality,
positive conjugate forms for arithmetic simplest type, and rational
isotropy in dimension at least 5. No high-dimensional field list,
ground-field degree bound, finite-angle cutoff or executable evidence
is used in this new branch.
