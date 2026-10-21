# A silver companion and the far-end opposed-fork exception list

19 October 2026. **New written local model and conditional classification,
awaiting independent review.** No computation or build was run. The
golden survivor note and all previous sources are unchanged.

## 1. The second model

Use the complete root pattern, notation and white projection of the
[golden opposed-fork note](arithmetic_dn_golden_opposed_fork_quad_survivor_20261019.md),
but now put

    x=beta=(1+sqrt(2))/2.

Thus both opposite products a.c=-x and b.d=-beta are dotted ambient
products. All other products are unchanged: the old vertex is whole
D_n, c.e=c.lL=c.d=-1/2, and a.d=b.c=0. Here n is 19 or 20,
N=n-2, L=n-3 and s=(N+2)/(4*N).

Again x*(beta-1)=1/4. The projected Gram blocks are therefore

    P_+=[-x,-1/2;-1/2,1-beta],
    P_-=[x+2*s,-1/2;-1/2,1+beta].

The first has rank one and one negative eigenvalue; the second is
positive definite. The full physical matrix has signature
(n positive,1 negative,1 zero). Its four adjacent corners are the
same D_n,A_(n-1)+A1,D_n,A_(n-1)+A1 diagrams. The centered Klein
parallelogram proof in the golden note gives a genuine compact
hyperbolic quadrilateral in the white perpendicular space.

At the nonidentity embedding of Q(sqrt(2)), x'=beta'=(1-sqrt(2))/2.
The block P_+ is positive semidefinite of rank one. The first
diagonal of P_- is x'+2*s>0, and its determinant is

    [7*N+6-(5*N+2)*sqrt(2)]/(4*N)>0.

For N=17 this numerator is 125-87*sqrt(2), whose comparison of
squares is 125^2-2*87^2=487>0. For N=18 it is 132-92*sqrt(2),
with 132^2-2*92^2=496>0. Thus the full conjugate Gram matrix is
positive semidefinite of rank n+1.

For norm-two roots every entry is integral over Q(sqrt(2)): the
two exceptional entries are both -(1+sqrt(2)). The same full-support
kernel relation holds,

    a+c+e+l1+...+lL-2*x*(b+d)=0.

Consequently every proper principal conjugate submatrix is positive
definite, so no proper principal physical submatrix is singular.
In particular an affine-subdiagram test does not remove the model.

These are genuine local Lorentz and arithmetic data, not a completed
compact n-polytope. If the model occurs in an actual arithmetic
polytope, its ground field is exactly Q(sqrt(2)), by the independent
old-D_n-plus-c determinant argument in the golden note.

## 2. Which opposed-fork branch is being classified

Suppose an actual quadrilateral at a whole D19/D20 fork/whisker
corner has its opposite corner also a whole D_n fork/whisker corner,
and both adjacent original corner weights are zero. The latter
forces a.d=b.c=0. Write c.d=-1/2 and c for the new fork.

Since d has degree one in this new whole D_n and d is perpendicular
to a, the zero-weight continuation list forces d perpendicular to
the entire white support W. Any nonisolated d from that list would
already have one white neighbor and would have degree two after
its contact with c, so could be neither a short leaf nor a fork.

The new fork c must contact e and one endpoint of the long A_L,
each simply. There are two distinct rooted possibilities: l1 and lL.
This section addresses the FAR choice lL. The projected c square
is s, but its white cross projection with a is s, not 1-s.
The correct formula is therefore

    X=x+s,     x=-a.c,
    x*(beta-1)=1/4,     beta=-b.d>1.               (1)

The near choice has a different additive constant and cannot stand
in for this branch.

For this classification only, retain the published high-dimensional
arithmetic restriction that quadratic ground fields at these
dimensions belong to

    Q(sqrt(2)), Q(sqrt(3)), Q(sqrt(5)), Q(sqrt(6)).

As in the previously reviewed isolated-quad note, this input is used
to reject specific quadratic fields, not to assume all fields have
degree two. Higher degrees are excluded below by their own arguments.

## 3. Finite a.c: only the golden profile survives

The value x=0 is impossible by (1). Independently, the roots
a,e,c,lL,...,l1 form an induced simple affine cycle, with positive
kernel coefficients all 1. The value x=1/2 gives the induced
simple affine triangle a,e,c and is impossible by compactness.

Now let x=cos(pi/m), m>=3. The cycle a-e-c-a puts 2*x in the
ground field as an algebraic integer. The old D_n plus c block is
physically independent. Its projection cost is

    q_c=1-s+2*N*(x+s)^2.

At every nonphysical embedding it is less than 1, so

    |x'+s|<sqrt(s/(2*N))<1/10.                    (2)

The elementary cosine-moment lemma in Section 5 of the
[isolated-d classification](arithmetic_dn_fork_whisker_quad_attachment_classification_20261019.md)
applies with kappa=s and r^2=s/(2*N). Its positive test coefficient
is (3*N-10)/(16*N), so the degree of cos(pi/m) is at most 3.
The complete resulting list is m=3,4,5,6,7,9.

Apart from m=5, every nontrivial value in that list violates (2).
For m=4 or 6 use their negative conjugate, at most -sqrt(2)/2,
and s<3/10. For m=7 a positive conjugate makes x'+s>s>1/10.
For m=9 its cos(7*pi/9) conjugate is less than -sqrt(2)/2 and
again violates (2). The m=3 value was excluded above.

Thus only x=cos(pi/5) remains. Equation (1) gives

    beta=(3+sqrt(5))/4,

exactly the golden profile already documented. Its admissibility
does not use the quadratic-field restriction.

## 4. Dotted or parallel a.c: only the silver profile survives

Suppose x>=1. Set z=4*beta^2. It is an algebraic integer in the
ground field K by the length-two cyclic product. By (1),

    4<z=4*(1+1/(4*x))^2<=25/4.                  (3)

The old D_n plus d block is physically independent, with determinant
ratio 1-n*z/8<0. Nonphysical positive definiteness gives

    0<z'<8/n<=8/19=:v.                          (4)

Moreover K=Q(z): a nonidentity embedding fixing z would preserve
the physical negative determinant ratio. This is the usual totally
real extension argument, rather than a ground-field degree assumption.

We classify (3)--(4) directly. Degree at least 4 would give a
nonzero positive integer norm less than

    (25/4)*v^3=3200/6859<1.

In degree 3 the norm is a positive integer below
(25/4)*v^2=400/361<2, hence is 1. Write the two small conjugates
as r,t. Then z>1/v^2=361/64. Its integer trace T satisfies

    T=z+r+t>=z+2/sqrt(z)>361/64+4/5>6,
    T<25/4+2*v<8,

so T=7. The integer second coefficient S has bounds

    S=z*(r+t)+r*t >=2*sqrt(z)+1/z>19/4+4/25>4,
    S<(25/4)*(2*v)+v^2=1964/361<6.

Thus S=5, and z would have minimal polynomial

    T^3-7*T^2+5*T-1.

Its discriminant is -44, so it has only one real root. This contradicts
total reality. Degree 3 is impossible.

In degree 2 the positive integer norm is below (25/4)*v=50/19<3,
so it is 1 or 2. The trace is greater than 4 and less than
25/4+v<7, hence is 5 or 6. Trace 5 gives fields Q(sqrt(21)) or
Q(sqrt(17)); trace 6 gives fields Q(sqrt(2)) or Q(sqrt(7)).
The published quadratic-field restriction excludes all but Q(sqrt(2)),
with norm 1 and

    z=3+2*sqrt(2).

In degree 1, z is 5 or 6 and K=Q; rational isotropy excludes compact
simplest-type arithmetic polytopes in these dimensions.

Finally z=3+2*sqrt(2) gives

    beta=(1+sqrt(2))/2,
    x=1/(4*(beta-1))=(1+sqrt(2))/2,

the silver profile in Section 1. Its norm-two cycle products are
indeed integral, so it cannot be rejected by the test that removes
the other quadratic candidates.

## 5. Capacity and a possible direction for a balanced rule

The same exact capacity proof as in the golden note applies to
the silver profile. The two new candidate walls for golden/silver
faces at the two old fork/whisker angles have equal projection costs;
their contact rows differ by (e_b-e_e)/2. Their forced scalar product
is 3/4>0. Hence at most one such profile can occur at a given old
vertex. Golden and silver profiles cannot occur in one arithmetic
polytope at all, because they force different actual ground fields.

These are capacity facts for the two exact profiles, not a complete
consumer theorem. If all other newly unpaid face patterns are
excluded, one could try a physical activation rule which retains
the unit on an exceptional quadrilateral and refunds the other
fork/whisker unit at that vertex. The capacity theorem would then
give at least one unit of saving per full D_n vertex. This would
still need every triangle, every other quadrilateral and every
whole-forest charge audited before being a valid rule. Even one
guaranteed saving would leave bare D20 cost 37, short of target 35.

The current result is the precise local exception list for the FAR
opposite-selected branch under the stated published field input:
golden Q(sqrt(5)) and silver Q(sqrt(2)). It does not establish the
global realizability of either profile, a universal refund, or a
new dimension ceiling.
