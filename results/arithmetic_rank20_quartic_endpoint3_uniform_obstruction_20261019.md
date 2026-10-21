# Quartic exclusion of both endpoint3-only dotted rows

19 October 2026. Root-agent hand proof, independently checked and
written out by batch_guard. No enumeration, calibration, scientific
job, Lean, TeX or build was run. Earlier proof and worker sources
are preserved unchanged.

## 1. Statement, fields and arithmetic interfaces

Let K be K15,K16,K20 or K24, where Km=Q(theta) and
theta=2cos(2pi/m), with its largest real root as physical embedding.
In the isolated-d, spectator-orthogonal actual-quadrilateral setup,
suppose x=-a.c>=1 is an EXPLICIT hypothesis. The endpoint3-only
contact is impossible at either endpoint, for every N>=8 for which
the exact row interfaces below hold.

For the rank20 frontier these are N=8,...,14. This statement does
not replace or extend the geometric inventory for c meeting the
spectator, nor assert that arbitrary ambient opposite walls have
product<=-1. It concerns endpoint3 alone, not endpoint3+e3.

The interfaces, derived from actual cycles, conjugate positivity and
corank, are those of
`arithmetic_rank20_isolated_long_dotted_height_bounds_20261019.md`,
SHA256 `9fc0583124cba17afed393bb569f83835e8f2c4be8ce5cfa64839f0369254516`:

    A=2x in O_K,       A_0>=2,
    Z=NA+h,           h=1 or N-1,
    C=[2Z^2-(N+2)(N+1)]/N in O_K,
    delta in O_K,     C*delta=4(N+1),
    Z_i^2<N+1 for every nonphysical embedding.

All four fields are totally real of degree4. Their presentations are

|field|minimal polynomial of theta|
|---|---|
|K15|X^4-X^3-4X^2+4X+1|
|K16|X^4-4X^2+2|
|K20|X^4-5X^2+5|
|K24|X^4-4X^2+1|

The cyclotomic presentations and degrees are named field inputs.
For K16,K20,K24 we use O_K=Z[theta], the standard
maximal-real-cyclotomic integer-ring theorem. For an integral element
B we may consequently write

    B=p+b*theta+c*theta^2+q*theta^3,
    p,b,c,q in Z.                                  (1)

The K15 argument does NOT need a quartic integral-basis assertion
or a classification of units. It uses only its quadratic subfield
F=Q(sqrt5) and O_F=Z[phi], phi=(1+sqrt5)/2. This is the ordinary
quadratic integer-ring input. The subfield can also be seen explicitly:

    sqrt5=2theta^3-6theta+1.

Squaring this expression and reducing by the displayed K15 polynomial
gives5. The physical sign is positive, equivalently from the
triple-angle identity at theta=2cos(2pi/15). Thus [K15:F]=2;
write tau for its nonidentity relative automorphism. Relative traces
and norms of integral elements belong to O_F, and the norm tower
identity applies. These standard number-field facts are the only
additional algebraic interfaces used below.

## 2. Uniform physical and conjugate bounds

At each nonphysical embedding C_i<-(N+1). Taking norms in
C*delta=4(N+1) gives

    0<C_0<256(N+1),
    Z_0^2<(N+1)(257N+2)/2=:H_N.                    (2)

For the FAR role h=1,

    (12N+1)^2-H_N=N(31N-211)/2>0

for N>=8. Therefore

    2<=A_0<12,       -1/2<A_i<1/4.                 (3)

The conjugate interval follows from |NA_i+1|<sqrt(N+1), using
sqrt(N+1)<=N/2-1 and sqrt(N+1)<=N/4+1 for N>=8.

For the NEAR role set U=A+1, so Z=NU-1. Since

    ((25/2)N-1)^2-H_N=N(111N-618)/4>0,

we obtain

    3<=U_0<25/2,     -1/4<U_i<1/2.                (4)

The same square-root comparisons prove the conjugate interval.
They follow by squaring from N(N-8)>=0.

In(3), no A_i is zero, since embeddings are injective and A is
nonzero. If any nonphysical A_i were positive, its absolute value
would be less than1/4, giving

    abs(Norm(A))<12*(1/4)*(1/2)^2=3/4,

impossible for a nonzero algebraic integer. All THREE A_i are
therefore negative. The norm is a negative integer with absolute
value less than12/8=3/2; hence

    Norm(A)=-1,       A_0>8.                       (5)

Write b_i=-A_i. Then 0<b_i<1/2, and with S=sum b_i,

    Tr(A)=A_0-S,       13/2<Tr(A)<12.               (6)

Similarly, in(4) any negative nonphysical U_i would imply

    abs(Norm(U))<(25/2)*(1/4)*(1/2)^2=25/32<1.

Thus all THREE U_i are positive. Their product with U_0 is a
positive integer less than25/16. Consequently

    Norm(U)=1,       U_0>8,
    8<Tr(U)=U_0+sum U_i<14.                        (7)

These norm arguments use the full integer ring but no power basis.

## 3. K16 and K24: far role

Their polynomials have the form X^4-4X^2+Q, with Q=2 or1.
Put D=4-Q, so D=2 or3 respectively. Their roots are +/-t,+/-u,
where t is physical and

    t^2=2+sqrt D,       u^2=2-sqrt D.

The trace in(1) is T=4p+8c, a multiple of4. By(6), T=8,
so p=2-2c. The two small-root values are both negative and
greater than-1/2. Their sum is

    -1<A(u)+A(-u)=4-2c*sqrt D<0.

Thus

    2/sqrt D<c<5/(2sqrt D).

For D=2 this is sqrt2<c<5/(2sqrt2); for D=3 it is
2/sqrt3<c<5/(2sqrt3). Both intervals are strictly between1 and2.
Neither contains an integer c. Both fields are excluded in the far role.

## 4. K20: far role

Its roots are +/-t,+/-u, with

    t^2=(5+sqrt5)/2,       u^2=(5-sqrt5)/2.

Here T=4p+10c is even, so(6) leaves T=8 or10. The small pair gives

    -1<A(u)+A(-u)=T/2-c*sqrt5<0.

If T=10, then5<c*sqrt5<6, which contains no integer c.
If T=8, then4<c*sqrt5<5, forcing c=2 and p=-3.
The large pair is then

    A(t)+A(-t)=4+2sqrt5<17/2,

using sqrt5<9/4. Since A(-t)>-1/2, we get A_0<9.
Equation(6) now gives 0<S=A_0-8<1. By the arithmetic-geometric
mean inequality,

    abs(Norm(A))=A_0*product b_i
        <=A_0*(S/3)^3<9/27<1,

contradicting(5).

## 5. K15: far role

Set V=Norm_(K/F)(A)=A*tau(A). Its physical value is negative,
while its other F-conjugate is the product of two negative values
of A. From(3),

    -6<V_0<0,       0<V'<1/4,
    Norm_F(V)=Norm_K(A)=-1.

Write V=(m+n*sqrt5)/2 with m,n integers of the same parity.
Its trace m lies in(-6,1/4), so m is one of -5,...,0, and

    m^2-5n^2=-4.

The values (m^2+4)/5 for these six m are

    29/5, 4, 13/5, 8/5, 1, 4/5.

Only4 and1 are integer squares. The required physical sign forces
(m,n)=(-4,-2) or(-1,-1). The latter gives V'=(sqrt5-1)/2>1/4,
so necessarily

    V_0=-2-sqrt5,       V'=sqrt5-2.

Now let R=Tr_(K/F)(A)=A+tau(A). Its physical value is positive
(indeed greater than15/2) and less than12. Its other conjugate is
the sum of the two remaining negative A-values. Their absolute
values are less than1/2 and their product is V'=sqrt5-2. Thus

    -1<R'<-19/20.                                  (8)

For the upper bound, AM-GM gives -R'>=2sqrt(V'). Since
sqrt5>223/100, one has 4V'>23/25>361/400, which proves the
strict comparison with19/20.

Write R=m+n*phi, with m,n integers. Since

    0<n*sqrt5=R_0-R'<13,

we have 1<=n<=5. The elementary bounds
223/100<sqrt5<9/4 give

    -5/8<phi'<-123/200.

Combining with(8), the necessary interval for m is one of:

|n|strict rational outer interval for m|
|---|---|
|1|(-77/200,-13/40)|
|2|(23/100,3/10)|
|3|(169/200,37/40)|
|4|(73/50,31/20)|
|5|(83/40,87/40)|

None contains an integer. This excludes the far K15 role without
a quartic integral basis or a unit classification.

## 6. K16 and K24: near role

Apply(1) to U=A+1. Its trace is a multiple of4; by(7), T=12.
Thus p=3-2c and the two small positive conjugates give

    0<U(u)+U(-u)=6-2c*sqrt D<1.

For D=3 this requires5<2c*sqrt3<6, which contains no integer c.
For D=2 it forces c=2,p=-1. The large pair is then

    U(t)+U(-t)=6+4sqrt2.

Since 0<U(-t)<1/2 and sqrt2>11/8, we get U_0>11.
Its trace12 therefore gives 0<S=sum U_i=12-U_0<1. AM-GM yields

    Norm(U)=U_0*product U_i
        <=U_0*(S/3)^3<(25/2)/27<1,

contradicting(7).

## 7. K20: near role

By(7), the even trace T=4p+10c is10 or12. The two small
positive conjugates satisfy

    0<T/2-c*sqrt5<1.

For T=12 this requires5<c*sqrt5<6 and no integer c exists.
For T=10 it forces c=2; but then p=(10-20)/4=-5/2 is
not an integer, contrary to the power-basis input.

## 8. K15: near role

Set V=Norm_(K/F)(U). Its two real values are positive and satisfy

    0<V_0<25/4,       0<V'<1/4,
    Norm_F(V)=Norm_K(U)=1.

In fact V_0=1/V'>4. Write V=(m+n*sqrt5)/2 as before. Its
integer trace satisfies 4<m<13/2, so m is5 or6. The norm equation
m^2-5n^2=4 would then give n^2=21/5 or32/5, neither an integer.
This excludes the final near role, again using only O_F rather
than a classification of quartic units.

## 9. Consequence and boundary

Both endpoint3-only dotted rows are excluded in K15,K16,K20,K24
under the explicit isolated-d, spectator-orthogonal actual-quad
interfaces, uniformly for N>=8. The field and integer-ring inputs
are separated above from the actual geometry and elementary norm
and coefficient arguments.

This does not cover endpoint3+e3 in every quartic field/N, nor c
meeting the spectator, nor the incomplete triangle and simultaneous
consumer composition. It does not install a global refund, change
a rank20 allowance, or lower a compact arithmetic dimension ceiling.
