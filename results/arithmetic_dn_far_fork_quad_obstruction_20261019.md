# Closing the far-fork quadrilateral continuation in D19/D20

19 October 2026. **New written mathematical argument, awaiting independent
review.** No computation, proof process, optimizer or build was run.
Earlier reports are unchanged. This is a local result, not a new global
weight rule or dimension bound.

## 1. The branch and its conclusion

Let n be 19 or 20, put N=n-2 and L=n-3=N-1, and let an actual finite
vertex of a compact hyperbolic Coxeter polytope have whole diagram D_n.
Its fork is a, its short whiskers are b,e, and its long arm is

    a--l1--l2-- ... --lL.

All old edges have label 3. Suppose an actual quadrilateral has cyclic
side roots a,b,c,d and white support W={e,l1,...,lL}. Assume

    c.b=c.d=0,
    d.l_(L-1)=-1/2,
    d is perpendicular to every other old root except possibly b.

The additional products c.a, c.e and c.li are initially unrestricted.
The finite opposite corner W+c+d will constrain them. Set b.d=-beta.

**Geometric claim.** This configuration is impossible when beta>=1.
No arithmetic input enters this geometric claim.

Combining it with the new
[finite-product obstruction, type D](arithmetic_dn_nonisolated_quad_finite_opposite_obstruction_20261019.md)
would remove the entire far-fork branch in the compact arithmetic
D19/D20 problem. That combination, like its constituent new notes,
awaits independent review. Neither constituent proves that every
quadrilateral is in this branch.

The condition c.b=c.d=0 is necessary if both corresponding original
Vinberg corner weights are zero, since adjacent marked roots always
have positive weight. It is an explicit hypothesis here.

## 2. Joint finite classification at the opposite corner

The white support with d is

    W+d = A1(e) + D_(L+1).

The D factor has fork l_(L-1), whiskers lL,d, and long endpoint l1.
Its rank is at least 17. The new root c is perpendicular to d, and
W+c+d is finite, being the actual opposite corner of the quadrilateral.

A finite connected diagram containing this whole D_(L+1) can only
be a larger D diagram: a path cannot contain its fork, and the finite
exceptional diagrams have too few nodes. A new attachment at a fork
leaf would create arms of lengths 1,2,and at least 14, which is not
finite. An internal attachment would create a second fork or a vertex
of degree four. Two attachments would create a cycle. Consequently
the only possible contact of c with the D factor is a simple bond to
its long endpoint l1.

If c also contacts e, the resulting diagram must still be of type D,
so c.e=-1/2. If c has no contact with the D factor, its e contact can
be any finite dihedral bond, or zero. Thus there are exactly two kinds
of joint continuation to consider:

1. c has no long-arm contact and c.e=-gamma, where gamma is zero or
   cos(pi/ell), ell>=3;
2. c.l1=-1/2, all other long-arm contacts vanish, and gamma is 0 or 1/2.

This classification includes every outside contact because the old
D_n is the whole rank-n diagram, with no further old factors.

## 3. A decreasing root formula for an actual quadrilateral

The following formula is also useful for the other nonisolated cases.
Project off the positive-definite white support W. In order a,b,c,d
the projected Gram matrix is

    [ s    -1/2   -X    -r ]
    [-1/2   1      0   -beta]
    [-X     0      t    -v ]
    [-r   -beta   -v     u ].

Here s,t,u>0 and X,r,v>=0. The opposite finite corner W+c+d gives
tu-v^2>0. The roots are actual normals in dimension n+1, and there
are n+2 of them, so this matrix has determinant zero.

Write D=beta^2-u. We are considering beta>=1 and below u=2/N<1,
so D>0. Eliminating c, or expanding the determinant, gives

    D*X^2-2*v*(beta/2+r)*X
      = t*[s*D+r*beta+r^2+u/4]+v^2*(s-1/4).

In this application s=(N+2)/(4N)>1/4, so the right side is positive.
The nonnegative root X must therefore be the positive root. Simplifying
its discriminant yields the exact expression

    X = v*(beta/2+r)/D
        + sqrt[(t+v^2/D)*(s+(r*beta+r^2+u/4)/D)].       (1)

For fixed nonnegative r,v and positive s,t,u, every displayed factor
is nonincreasing as beta increases through [1,infinity), provided
u<1. Indeed D increases,

    d/d beta [(beta/2+r)/D]
        =(-beta^2/2-2*r*beta-u/2)/D^2 < 0,

and, with C=r^2+u/4>0,

    d/d beta [(r*beta+C)/D]
        =(-r*beta^2-r*u-2*C*beta)/D^2 < 0.

The term t+v^2/D is nonincreasing too. Thus the right side of (1)
is at most its value at beta=1. This is an exact inequality, with no
numerical optimization or parameter cutoff.

For the far-fork d contact, direct inversion of A1+A_L gives

    s=(N+2)/(4N),     u=2/N,     r=1/N.

Let x=-a.c in the ambient Gram matrix. Acuteness and the Coxeter
property give

    x=0, or x=cos(pi/m) with m>=3, or x>=1.

In particular any nonzero x is at least 1/2.

## 4. No long-arm contact of c

In the first continuation of Section 2, white projection gives

    t=1-gamma^2,     v=0,     X=x+gamma/2.

Putting beta=1 in (1) therefore gives

    X^2/t <= Z,     Z=(N+6)/(4*(N-2))<=23/60<1/2.      (2)

If gamma=0 and x is nonzero, (2) forces x<cos(pi/4), so x=1/2.
But then the actual roots {a,b,e,l1,c} form the affine D4 star:
center a and four simple leaves b,e,l1,c. This is impossible in a
compact chamber. If gamma=x=0, then c is perpendicular to the whole
old positive D_n space. Its perpendicular is a timelike line, which
cannot contain a unit spacelike normal. This also is impossible.

Now suppose gamma>=1/2. If x>=1/2 then

    (x+gamma/2)^2/(1-gamma^2)>=3/4,

contrary to (2). Hence x=0. Equation (2) now implies

    gamma^2 <= (N+6)/(2*N+4) <=23/38.

But cos^2(pi/5)=(3+sqrt(5))/8>23/38, so the only possible finite
labels are gamma=1/2 or gamma=sqrt(2)/2.

For gamma=1/2, the induced root diagram on

    {a,b,e,c,l1,l2,l3,l4,l5}

is affine E8, with arms of lengths 1,2,5 at a. For gamma=sqrt(2)/2,
the induced diagram on {a,b,l1,e,c} has simple bonds a-b,a-l1,a-e
and a terminal 4-bond e-c. It is affine B4. Both contradict
compactness, with explicit positive-null certificates in Section 6.
This closes the first continuation.

## 5. A simple contact of c with l1

The second continuation has gamma=0 or 1/2, and

    t=(N+1)/(2*N)-gamma^2,
    v=1/N,
    X=x+(N-1)/(2*N)+gamma/2.

At beta=1 the two square-root factors in (1) reduce to

    t+v^2/(1-u) = (N-1)/(2*(N-2))-gamma^2,
    s+(r+r^2+u/4)/(1-u)=Z.

The first term in (1) is (N+2)/(2*N*(N-2))<1/20. Also

    (N-1)/(2*(N-2))<=8/15,
    Z<=23/60,
    (8/15)*(23/60)<1/4.

It follows that X<11/20. If x were nonzero, however, it would give

    X>=1/2+(N-1)/(2*N)>=33/34,

a contradiction. Thus x=0.

When gamma=0, the actual roots {a,l1,b,e,c,l2} form an affine D5
diagram with adjacent forks a,l1 and simple leaves b,e at a and c,l2
at l1. When gamma=1/2, the roots {a,e,c,l1} form an induced
simple 4-cycle, affine A3. Again both are impossible by compactness.
This closes the second continuation and proves the geometric claim.

## 6. Explicit positive-null witnesses and compactness

All coefficients here are for square-one roots. A coefficient vector
v on an induced root set is a positive null vector precisely when
each coefficient equals the sum of its neighbor coefficients times
the corresponding contact cosine.

* Affine D4: coefficient 2 at a, coefficient 1 at each of
  b,e,l1,c.
* Affine E8 with arms 1,2,5: coefficient 6 at a; coefficient 3
  at b; coefficients 4,2 along e,c; coefficients 5,4,3,2,1
  along l1,...,l5.
* Affine B4: coefficients 2 at a,e; coefficients 1 at b,l1;
  coefficient sqrt(2) at c. Its only nonsimple bond is e-c,
  of cosine sqrt(2)/2.
* Affine D5: coefficients 2 at a,l1, coefficients 1 at
  b,e,c,l2.
* Affine A3: coefficient 1 at a,e,c,l1 on the simple 4-cycle.

These witnesses are verified directly by the indicated neighbor
equations. They do not require a computational affine-diagram table.
All unlisted pairings vanish by the whole old D_n diagram and the
complete new contacts in Section 2, together with x=0 where used.

For completeness, let z be the sum of these actual normals with their
positive coefficients. Its pairing with an interior future timelike
point is negative, so z is nonzero. The kernel equation gives z^2=0.
Every selected wall is orthogonal to z, and acuteness makes every
other outward wall nonpositive against z. Choosing the future ray
therefore gives a feasible ideal point, contrary to compactness.

## 7. Consequence for the refund programme

The whole far-next-to-end d continuation is now closed at n=19,20,
under the stated orthogonalities: beta>=1 is excluded by this
geometric argument, and finite beta is excluded arithmetically in
the companion note. The earlier isolated-d continuation is separately
closed under the same c.b=c.d=0 hypotheses.

The two far-end d continuations, of labels 3 and 4, still have dotted
or parallel ambient b.d branches. A unit refund requires those too,
as well as a simultaneous proof for every selected triangle and
quadrilateral. No such refund or new compact arithmetic ceiling is
claimed in this note.
