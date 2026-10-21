# A universal endpoint3-only norm obstruction

19 October 2026. Root-agent hand proof, independently checked and
written out by batch_guard. No enumeration, scientific job,
calibration, Lean, TeX or build was run. Earlier field-specific
proofs and prepared worker sources remain unchanged historical
evidence; no worker certificate is claimed.

## 1. Algebraic theorem and actual-quad scope

Let K be a totally real number field of degree d, with one physical
embedding. Let N>=8, h=1 or N-1, and suppose A,delta belong to its
full integer ring and satisfy

    A_0>=2,
    Z=NA+h,
    C=[2Z^2-(N+2)(N+1)]/N in O_K,
    C*delta=4(N+1),
    Z_i^2<N+1 at every nonphysical embedding.       (1)

Then the NEAR role h=N-1 is impossible for every K. The FAR role
h=1 is impossible if d=1 or d>=3. If a far solution exists in
degree2, necessarily

    K=Q(sqrt2),       A_0=1+sqrt2.

This last necessary exception does not assert that a solution of
all the original geometric and arithmetic conditions exists.

These are the exact endpoint3-only interfaces of
`arithmetic_rank20_isolated_long_dotted_height_bounds_20261019.md`,
SHA256 `9fc0583124cba17afed393bb569f83835e8f2c4be8ce5cfa64839f0369254516`.
In that application A=2x is supplied by an actual cycle, x=-a.c>=1
is an EXPLICIT hypothesis, d is isolated from the old white graph,
and c is orthogonal to the spectator. This note does not infer
dottedness of arbitrary ambient opposite walls. In the rank20
frontier N=8,...,14; the algebraic result itself is valid for
all N>=8 whenever(1) has been established.

No explicit integral basis, field list, unit classification, or
primitivity of A is assumed. We use only ordinary integral norms,
traces and characteristic polynomials over the FULL degree d.
If A lies in a proper subfield, its embeddings here are still
counted with their full multiplicities.

## 2. Common norm and physical bounds

Physically Z_0>=2N+1, since h>=1. Hence

    C_0>=7N+5>4(N+1)>0.

For d=1, C and delta are positive integers with C*delta=4(N+1),
contradicting this inequality. Henceforth d>=2.

Every nonphysical C_i is less than -(N+1). Taking norms in(1),
the nonzero integer norm of delta gives

    C_0<4^d*(N+1),
    Z_0^2<(N+1)*[N(4^d+1)+2]/2=:B_d,N.            (2)

Put H=3*2^(d-2). A direct expansion yields

    (NH+1)^2-B_d,N
      =(N/2)*[(4^d/8-1)N+4H-4^d-3].

The bracket increases with N>=8 and at N=8 equals4H-11>0.
Thus the far role has

    2<=A_0<H.

For the near role set U=A+1, so Z=NU-1. Similarly,

    [N(H+1/4)-1]^2-B_d,N
      =(N/2)*[(4^d/8+H-7/8)N-4H-4^d-4].

Again the bracket increases with N>=8 and at8 equals4H-11>0.
Consequently

    3<=U_0<H+1/4.                                  (3)

The two elementary comparisons

    sqrt(N+1)<=N/2-1,
    sqrt(N+1)<=N/4+1

hold for N>=8 by squaring. Applied to the nonphysical Z bounds,
they imply

    far:  -1/2<A_i<1/4,
    near: -1/4<U_i<1/2.                            (4)

## 3. Far role: sign and unit norm

No conjugate of the nonzero A can vanish. If any nonphysical A_i
were positive, (4) would give

    abs(Norm(A))<H*(1/4)*(1/2)^(d-2)=3/4,

contradicting the nonzero integer norm. Thus every nonphysical
A_i is negative. Write b_i=-A_i, so 0<b_i<1/2. Then

    abs(Norm(A))<H*(1/2)^(d-1)=3/2,

and therefore abs(Norm(A))=1. For each nonphysical i,

    b_i=1/[A_0*product_(j nonphysical,j!=i) b_j]
       >=2^(d-2)/A_0
       >1/3.                                      (5)

The first comparison is strict for d>=3; equality there is possible
when d=2 because the product is empty. The final strict inequality
follows in every case from A_0<H.

### Far role, d>=4

Set V=2A+1, which is integral. By(5),

    0<V_i<1/3,       0<V_0<2H+1.

Hence

    0<Norm(V)<[3*2^(d-1)+1]/3^(d-1)<=25/27<1.

The last bound is attained as an outer bound at d=4 and decreases
with d, since it is 3(2/3)^(d-1)+(1/3)^(d-1). This contradicts
integrality of the nonzero norm.

### Far role, d=3

Here H=6. The unit norm is Norm(A)=1, and b_1,b_2<1/2 imply
A_0>4. The integral trace T=A_0-b_1-b_2 lies strictly between3
and6, so T is4 or5. The characteristic polynomial over K/Q is

    f(X)=X^3-TX^2+SX-1,

where S is an integer. Since the norm is1,

    S=Tr(1/A)=1/A_0-1/b_1-1/b_2<-15/4,

so S<=-4. All three roots exceed -1/2, so f(-1/2)<0, giving

    S>-(2T+9)/4.

For T=4 or5 this forces S=-4. The two polynomials then have
discriminants

|T|S|discriminant|
|---|---|---|
|4|-4|-59|
|5|-4|-231|

A real-rooted polynomial has nonnegative discriminant, including
the possibility of repeated real roots. Both cases are impossible.

### Far role, d=2

Now H=3, the other conjugate is in(-1/2,0), and Norm(A)=-1.
The integral trace lies in(3/2,3), so it equals2. Its characteristic
polynomial is X^2-2X-1. The positive physical root is1+sqrt2,
and its irrationality forces the degree2 field to be Q(sqrt2).
This proves only the stated necessary exception.

## 4. Near role: positivity and unit norm

No conjugate of the nonzero integral U can vanish. If any
nonphysical U_i were negative, (3),(4) would imply

    abs(Norm(U))<(H+1/4)*(1/4)*(1/2)^(d-2)
       =3/4+2^(-d-2)<1.

Thus every nonphysical U_i is positive. The positive integral norm
satisfies

    0<Norm(U)<(H+1/4)*(1/2)^(d-1)
       =3/2+2^(-d-1)<2.

Consequently Norm(U)=1.

### Near role, d=2

We have 3<=U_0<13/4 and 0<U_1<1/2. Therefore its integral
trace lies in(3,15/4), an interval with no integer.

### Near role, d=3

Here U_0>4, because the other two positive conjugates are less
than1/2 and their product with U_0 is1. Its trace is an integer
T with4<T<29/4, so T is5,6 or7. Write its integral
characteristic polynomial as

    f(X)=X^3-TX^2+SX-1.

We have S=Tr(1/U)>4, so S>=5. As the two nonphysical roots
are below1/2 and the physical root exceeds4, f(1/2)<0. Hence

    S<(2T+7)/4.

This excludes T=5,6. The remaining case is T=7,S=5, whose
discriminant is -44, again incompatible with real roots.

For explicit arithmetic, the discriminant of X^3-TX^2+SX-1 is
T^2S^2-4S^3-4T^3-27+18TS; this gives all three negative
values above. No enumeration or polynomial-root approximation
enters these contradictions.

### Near role, d>=4

Put k=2^(d-2), so H=3k. From Norm(U)=1 and the bounds
0<U_j<1/2, each nonphysical conjugate satisfies

    U_i>k/U_0>4k/(12k+1).

The integral element V=2U-1 therefore has

    abs(V_i)<(4k+1)/(12k+1),
    0<V_0<6k-1/2.

It is nonzero physically, so its norm is a nonzero integer.
For d=4, k=4 and the explicit bound is

    abs(Norm(V))<(47/2)*(17/49)^3
       =230911/235298<1,

a contradiction. For d>=5 the ratio (4k+1)/(12k+1) decreases
with k and is at most17/49<7/20. Also V_0<3*2^(d-1). Hence

    abs(Norm(V))<3*(7/10)^(d-1)<=7203/10000<1,

the final contradiction. This proves the near exclusion in all
degrees, together with the separate d=1 argument.

## 5. Consequence and remaining boundary

The two endpoint3-only dotted rows are excluded in every totally
real field of degree>=3, and the near row is excluded in every
degree. The only degree2 far possibility permitted by this argument
has actual ground field Q(sqrt2). This does not create an exception
in a higher-degree field merely containing Q(sqrt2): the proof
uses all embeddings of the actual K throughout.

In particular both rows are excluded in the eight current fields
Q(sqrt3),Q(sqrt6),K9,K11,K15,K16,K20,K24, without any of their
explicit integer-ring presentations. The previous field-specific
endpoint3-only arguments remain valid but are superseded for this
purpose by the present proof.

This theorem does not include endpoint3+e3, c meeting the spectator,
or the incomplete triangle and simultaneous-consumer composition.
No new global refund, rank20 allowance, or compact arithmetic
dimension ceiling is asserted.
