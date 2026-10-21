# The remaining far-simple D19/D20 quadrilateral cannot be wholly unpaid

19 October 2026. **New written research proof, awaiting independent
review.** No enumeration, numerical calculation, proof program, Lean
process or build was run. All earlier reports are unchanged.

## 1. Precise statement

Let n be 19 or 20, N=n-2 and L=n-3=N-1. Let an actual finite vertex
of a compact arithmetic simplest-type Coxeter polytope have whole
diagram D_n, with fork a, whiskers b,e and simple long arm

    a--l1-- ... --lL.

Suppose an actual quadrilateral has successive side roots a,b,c,d,
white support W={e,l1,...,lL}, and

    d.lL=-1/2,
    d perpendicular to every other old root except b,
    c.b=c.d=0,
    b.d=-beta,     beta>=1.

Then the opposite corner (c,d) has positive original Vinberg weight,
or the configuration is impossible. In particular it cannot be a
quadrilateral whose three corners other than (a,b) all have zero
original weight. The positive-corner alternative is essential; we
do not claim that all quadrilaterals with the displayed contacts are
geometrically impossible.

Together with the preceding isolated, finite-product, far-fork and
far-double notes, this would show that every actual quadrilateral
at the old D19/D20 fork/whisker angle has another corner of positive
original weight. All new constituent notes await independent review.
Even this statement is not a simultaneous refund theorem: the other
positive corner might itself be selected by a proposed modification.

## 2. Complete joint finite inventory

The opposite corner W+c+d is finite. Before c is added it contains
A1(e) and the long path

    l1-- ... --lL--d

of type A_(L+1). Also c.d=0. Finite classification therefore gives
exactly the following possibilities for the contacts of c with W:

|Long-arm contact of c|Contact with e|Permitted bonds|
|---|---|---|
|None|None or e|Any finite e bond|
|l1|None or e|Long bond 3; e bond 3 or 4 if present|
|l1|None|Long bond 4|
|l2|None|Long bond 3|
|lL|None|Long bond 3|

Indeed a one-root finite extension of a path of this rank is a path
extension at an endpoint, or a D-fork extension at a next-to-end
vertex. The far endpoint is d, whose contact with c is prohibited.
A second contact creates a cycle. Joining the isolated e can only
lengthen the near-end path; it cannot turn an internal 4-bond or arms
of lengths 1,2,and at least 15 into a finite diagram. The exceptional
finite diagrams have too few nodes.

In the last row W+c+d is D_(n-1)+A1(e), with fork lL and its two
short whiskers c,d. The marked pair (c,d) lies in an end segment of
order 3, hence has positive original Vinberg weight. This is the
positive-corner alternative in the statement, and disposes of that
row when testing a wholly unpaid quadrilateral.

## 3. Common Schur bound for the other rows

Use square-one roots and put x=-a.c. Thus x is zero, a finite cosine
cos(pi/m), m>=3, or at least 1. Every nonzero x is at least 1/2.
After projecting off W, use the notation of the
[far-fork positive-root formula](arithmetic_dn_far_fork_quad_obstruction_20261019.md):

    a_bar^2=s, b^2=1, c_bar^2=t, d_bar^2=u,
    a_bar.c_bar=-X, a_bar.d_bar=-r, c_bar.d_bar=-v,
    a.b=-1/2, b.c=0, b.d=-beta.

Here the d parameters are

    s=(N+2)/(4*N),     u=(N+1)/(2*N),     r=1/(2*N).

Formula (1) of that note and its monotonicity in beta give

    X <= v*(1/2+r)/(1-u)
          +sqrt[(t+v^2/(1-u))*Z],
    Z=s+(r+r^2+u/4)/(1-u)=(N+3)/(2*(N-1))<=5/8.       (1)

The next three rows of the inventory can be removed geometrically.

### A simple contact with l1

Let c.e=-gamma, where gamma is 0,1/2 or sqrt(2)/2. Then

    t=(N+1)/(2*N)-gamma^2,
    v=1/(2*N),
    X=x+(N-1)/(2*N)+gamma/2.

The first term in (1) is

    F=(N+1)/(2*N*(N-1))<1/28,

and t+v^2/(1-u)=N/(2*(N-1))-gamma^2.
For gamma=0, (1) is less than

    1/28+sqrt[(17/32)*(5/8)] <1/28+3/5<2/3.

Any nonzero x gives X>=33/34, a contradiction. Thus x=0,
and the induced diagram {a,l1,b,e,c,l2} is affine D5, with coefficient
2 on its two adjacent forks a,l1 and coefficient 1 on its four
leaves. Compactness excludes it.

If gamma>=1/2, (1) is less than

    1/28+sqrt[(9/32)*(5/8)]
       =1/28+3*sqrt(5)/16<1/2.

But X>1/2 already from its displayed value, even if x=0. This excludes
both remaining e labels in this row.

### A 4-bond contact with l1

Here c.e=0 and

    t=1/N,     v=sqrt(2)/(2*N),
    X=x+sqrt(2)*(N-1)/(2*N).

The first term in (1) is sqrt(2)*F<3/56, and
t+v^2/(1-u)=1/(N-1)<=1/16. Hence

    X<3/56+sqrt(5/128)<3/56+1/5<1/3.

Its displayed value is greater than (N-1)/(2*N)>=8/17>1/3.
This row is impossible.

### A simple contact with l2

Here c.e=0 and

    t=2/N,     v=1/N,     X=x+(N-2)/(2*N).

The first term in (1) is 2*F<1/14, and
t+v^2/(1-u)=2/(N-1)<=1/8. Therefore

    X<1/14+sqrt(5/64)<1/14+9/32<2/5.

But X>2/5 since (N-2)/(2*N)>=15/34>2/5. This excludes this row too.

## 4. No long-arm contact: two arithmetic exceptions

It remains that c has no long-arm contact. Put c.e=-gamma, with
gamma zero or any finite cosine. Then

    t=1-gamma^2,     v=0,     X=x+gamma/2.

Inequality (1) reduces to

    (x+gamma/2)^2/(1-gamma^2)<=Z<=5/8.                (2)

If gamma=0, the only native nonzero x permitted by (2) are 1/2 and
sqrt(2)/2: the next value cos(pi/5) has square greater than 5/8.
The value x=0 makes c perpendicular to the whole positive old D_n
space and is impossible. The value x=1/2 gives the affine D4 star
with center a and leaves b,e,l1,c. Thus only

    x=sqrt(2)/2,     gamma=0                         (3)

remains in this subcase.

If gamma>=1/2 and x>=1/2, the left side of (2) is at least 3/4,
which is impossible. Hence x=0. Inequality (2) then gives

    gamma^2<=2*(N+3)/(3*N+5)<=5/7<3/4.

Only labels 3,4,5 are possible for gamma. Label 3 gives the affine
E8 arms (1,2,5) at a; label 4 gives the affine B4 on
{a,b,l1,e,c}. The explicit positive-null witnesses are those in the
far-fork note and do not involve d. The only remaining subcase is

    x=0,     gamma=cos(pi/5).                       (4)

We now exclude (3) and (4) arithmetically, without a ground-field list.

## 5. Exception (3): the opposite product is rational nonintegral

The projected c square is 1 and X=sqrt(2)/2, with v=0.
The exact quad determinant relation gives

    (N-2)*beta^2-2*beta-N=0.

Its positive solution is beta=N/(N-2); its other solution is -1.
The old simple path from b to lL together with d is a closed cycle,
so 2*beta is an algebraic integer in the arithmetic ground field K.
But

    2*beta=2+4/(N-2)

is rational and not an integer for N=17,18. This is impossible.
No assertion that the 4-bond unit cosine itself belongs to K is used.

## 6. Exception (4): the quadratic discriminant cannot be a field square

Put

    u0=gamma^2=(3+sqrt(5))/8,
    B=2*beta.

The two-cycle at e,c puts u0 in K. The simple long cycle at b,d
puts B in O_K. The determinant costs against the old D_n space are

    q_c=n*u0/2,
    q_d=n*beta^2/2+beta+1/2,
    q_cd=gamma*((n-2)*beta+1)/2.

The old roots and c,d number n+2 in dimension n+1. Their two-root
Schur complement has determinant zero, so

    (q_c-1)*(q_d-1)=q_cd^2.

After expansion this becomes

    A*B^2-4*(1-u0)*B+C=0,                          (5)
    A=(2*n-2)*u0-n,
    C=4-2*(n+1)*u0.

Here A>0 and C<0. Hence (5) has exactly one positive real root.
At B=8/3, nine times its left side is

    (110*n-50)*u0-64*n-60>0.

Indeed u0>13/20 gives a lower bound (15*n-185)/2>0 at n=19,20.
Consequently the physical B satisfies

    2<=B<8/3.                                      (6)

We next justify a ground-field conclusion rather than assuming it.
The old D_n plus d block is physically independent and Lorentzian,
with determinant ratio 1-q_d<0. At every nonphysical embedding the
corresponding independent block is positive definite. Thus q_d'<1.
As in the finite-product note this gives

    |B'|<3/5                                       (7)

for every nonphysical conjugate. If B had degree at least 3, its
nonzero integer norm would have absolute value less than

    (8/3)*(3/5)^2=24/25<1,

by (6)--(7), an impossibility. Also K=Q(B): a nonidentity embedding
of the totally real K fixing B would fix q_d>1, contrary to conjugate
positive definiteness. Every proper subfield has such an embedding.
Since u0 in K generates Q(sqrt(5)), it follows that

    K=Q(B)=Q(sqrt(5)).

The discriminant of (5), using 16*u0^2=12*u0-1, is

    Delta=4*n*(n-10)*u0-n^2+16*n+16.

If B belonged to K, this discriminant would be a square in K: from
(5), (2*A*B-4*(1-u0))^2=Delta. Its norm to Q would therefore be
a rational square. The two actual values are

|n|Delta|Norm(Delta)|
|---:|---|---:|
|19|(431+171*sqrt(5))/2|9889|
|20|236+100*sqrt(5)|5696|

But 99^2<9889<100^2 and 75^2<5696<76^2. Neither integer is a rational
square. This contradiction excludes exception (4).

The nonphysical determinant inequalities in this section are invariant
under the positive diagonal normalizations of arithmetic roots.
They involve B, u0 and squared or closed-cycle products in K, not an
unsupported conjugation of individual unit-normal entries.

## 7. Exact scope of the result

The joint finite inventory in Section 2 is exhaustive. Its far-leaf
row supplies a retained positive original weight at (c,d); every
other row has now been excluded for beta>=1. The finite beta branch
of the same d continuation is separately excluded in the
[finite-product note](arithmetic_dn_nonisolated_quad_finite_opposite_obstruction_20261019.md).

Together with the other new local notes, this is evidence for the
following statement awaiting composition review: an actual
quadrilateral at a full D19/D20 fork/whisker corner cannot have all
three other original corner weights zero. It is not yet a proof
that both fork/whisker weights may be deleted simultaneously at
every D19/D20 vertex. A positive angle supplied here must remain
paid under any proposed global rule; that separate consumer check
has not been performed in this note. The compact arithmetic ceiling
20 is unchanged.
