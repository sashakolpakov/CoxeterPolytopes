# A common norm obstruction for three dotted long-contact types

19 October 2026. Root-agent argument, independently checked and
written out by batch_guard. Paper only: no scientific job,
enumeration, Lean, TeX, build or calibration was run.

This note replaces the need to enumerate next-to-end3, endpoint4,
and endpoint3+e4 dotted long-contact types. It leaves the prepared
endpoint3-only worker unchanged. It does
not establish all simultaneous half-refund consumers or a new bound.

## 1. Common algebraic obstruction

Let K be any totally real number field of degree d, N>=8 an integer,
and V,delta elements of its full integer ring O_K. Distinguish one
physical embedding, denoted0. Neither of the following systems is
possible:

**Far system:**

    V_0>=2,
    -4/N<V_i<0 for every i nonphysical,
    C=NV^2+4V-2,       C*delta=8.

**Near system:**

    V_0>=3,
    0<V_i<4/N for every i nonphysical,
    C=NV^2-4V-2,       C*delta=8.

There is no field-list, degree, monogenicity or integral-basis
assumption. Ordinary integrality of V and delta is the input.

In either system C is integral, positive physically, and every
nonphysical conjugate satisfies C_i<-2. Indeed, the nonconstant
part of C_i is respectively

    N*V_i*(V_i+4/N),    N*V_i*(V_i-4/N),

which is strictly negative in the stated interval. Since C*delta=8,
the nonzero integral norm of delta yields

    abs(Norm(C))<=8^d.

If d>=2, the strict bounds abs(C_i)>2 imply

    C_0<8^d/2^(d-1)=2*4^d.                         (1)

For d=1, integrality directly gives C_0<=8. In the far system
C_0>=4N+6>8, and in the near system C_0>=9N-14>8. Thus
degree1 is excluded separately; no vacuous strict-product argument
is used for it.

### Far system

For d>=2, C_0>NV_0^2 and abs(V_i)<4/N. Therefore

    abs(Norm(V))^2
       =V_0^2*product_(i>0) V_i^2
       <[2*4^d/N]*(16/N^2)^(d-1)
       =(8/N)^(2d-1)
       <=1.

But V is a nonzero algebraic integer, so its norm is a nonzero
integer. This is a contradiction, including at the endpoint N=8.

### Near system, N=8

Here C=2D with

    D=4V^2-2V-1,       D*delta=4.

D belongs to O_K and D=2(2V^2-V)-1, so its principal ideal is
coprime to(2). Since (D) divides(4), D must be a unit. Physically
D_0>=29>1. At each nonphysical embedding 0<V_i<1/2, and

    D_i=4V_i*(V_i-1/2)-1<-1.

Consequently abs(Norm(D))>1, contradicting that D is a unit.
This argument uses ordinary ideal divisibility in O_K and does
not presume unramified2 or a particular power basis.

### Near system, N>=9

All conjugates of V are positive. Norm(V) is therefore a positive
integer, and the strict nonphysical bounds give

    V_0>(N/4)^(d-1)>=(9/4)^(d-1)=k_d.

For u>=k_d and N>=9, the function Nu^2-4u-2 increases with u,
because 2Nu-4>0. Hence

    C_0>9k_d^2-4k_d-2=:g_d.

For d=2,

    g_2=553/16>32=2*4^2.

Since k_(d+1)=(9/4)k_d, direct subtraction gives

    g_(d+1)-4g_d=(153/16)k_d^2+7k_d+6>0.

Induction proves g_d>2*4^d for every d>=2. This contradicts(1).
Both algebraic systems are now excluded in every degree.

## 2. Application to the actual long-contact rows

Use the exact actual-quadrilateral setup and arithmetic interfaces
of `arithmetic_rank20_isolated_long_dotted_height_bounds_20261019.md`,
SHA256 `9fc0583124cba17afed393bb569f83835e8f2c4be8ce5cfa64839f0369254516`.
In particular, d is isolated from the old white graph, c is
orthogonal to the spectator, and x=-a.c>=1 is an EXPLICIT hypothesis.
No assertion that arbitrary ambient opposite facet walls must have
product<=-1 is introduced.

Write N=r-2. In the rank20 problem 8<=N<=14. The algebraic proof
itself works for any N>=8 whenever the same actual-quad interfaces
have been established. The unit-normal finite contacts and actual
cycles provide the indicated integral variable A. Set
delta=4*beta^2-4, where beta=-b.d>1; cyclic-product integrality
gives delta in O_K. Corank supplies C*delta=8. The nonphysical
proper-block positivity supplies Z_i^2<4.

### Next-to-end simple contact3

Here A=2x is integral and A_0>=2. The two roles are

    far:  Z=NA+2,
          C=(Z^2-2N-4)/N=NA^2+4A-2;

    near: Z=NA+N-2,
          U=A+1,
          Z=NU-2,
          C=(Z^2-2N-4)/N=NU^2-4U-2.

In the far case Z_i^2<4 gives -4/N<A_i<0, so take V=A in
the far system. In the near case it gives 0<U_i<4/N, while
U_0>=3, so take V=U in the near system. Both are impossible.

### Endpoint double contact4

The integral ACTUAL CYCLE variable is A=2sqrt2*x, so A_0>=2sqrt2.
It is not necessary, and would be incorrect in general, to impose
sqrt2 in K merely from this notation. Its product A is the field
element supplied by the cycle. The two roles are

    far:  Z=NA+2,
          C=(Z^2-2N-4)/N=NA^2+4A-2;

    near: Z=NA+2N-2,
          U=A+2,
          Z=NU-2,
          C=(Z^2-2N-4)/N=NU^2-4U-2.

The far system applies with V=A, since A_0>=2sqrt2>2. The near
system applies with V=U, since U_0>=2sqrt2+2>3. The same strict
conjugate intervals follow from Z_i^2<4. No square root is adjoined
to K at any step of the norm proof.

### Endpoint simple contact3 together with e4

Here A=2x is integral from the actual simple-long cycle. The actual
a,e,c triangle has cyclic product, up to sign, A*sqrt2, so this
product belongs to K and is integral. Since A is nonzero, division
gives the ACTUAL positive physical sqrt2 in K. It is integral by
its monic polynomial X^2-2. No separate integer-ring theorem or
choice of signs at other embeddings is needed.

The original row is

    Z=2NA+2h+N*sqrt2,
    C=(Z^2-2N-4)/N,
    C*delta=8,
    Z_i^2<4 for every i nonphysical,

with h=1 or N-1. Set V=2A+sqrt2, which is integral and satisfies
V_0>=4+sqrt2>2. For the far role h=1 one has

    Z=NV+2,
    C=NV^2+4V-2,
    -4/N<V_i<0 for every i nonphysical.

The far system applies. For the near role set U=V+2; then

    Z=NU-2,
    C=NU^2-4U-2,
    0<U_i<4/N for every i nonphysical,
    U_0>=6+sqrt2>3.

The near system applies. In both cases the signs of the conjugate
sqrt2 values are already encoded by the actual field embeddings;
replacing them all by a positive square root would be incorrect.

Thus all THREE dotted long-contact types above are impossible over
every totally real ground field in the stated isolated-d,
spectator-orthogonal actual-quad setup. The distinct endpoint3-only
and endpoint3+e3 types, spectator-contact rows, and full triangle
and simultaneous-consumer composition are not asserted settled by
this theorem. In particular no global refund or dimension ceiling
has been changed here.
