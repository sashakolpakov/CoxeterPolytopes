# The cubic B16 short cap has a unique polynomial

21 October 2026. Complete manual supporting candidate by root, pending
outside-author mathematical review and acceptance. No worker, enumeration
program, symbolic import, test, Lean or TeX process ran. The proof below
uses integer inequalities and three groups of possible traces; it does
not assume an unexecuted field census.

## Statement and inherited interface

Assume the actual compact arithmetic seventeen-dimensional situation of

`results/arithmetic17_b16_actual_i_short_opposite_isolation_and_degree_restriction_20261021.md`,

SHA-256
`51910bc233e16870d2322e66c9f2c2766076d757749088e7743277b8ddbed248`.

That accepted written result proves that the short B16 opposite is
isolated within B16 and dotted. With its unit spatial expression
`gamma P16+beta e`, set B=4gamma squared. The actual ground field
is K=Q(B); its degree is two or three. Its necessary bounds are

    B in O_K, 4<B<=(441/16)(1-beta^2),
    0<B_sigma<1/4 at every nonphysical embedding,
    beta=0 or cos(pi/m), m>=3.

**Proposition.** If the actual field has degree three, then beta=0
and B is the unique root in(24,25) of

    F(X)=X^3-25X^2+10X-1.

In particular the cubic branch has one actual field, of discriminant473,
and the short cap is also perpendicular to the old isolated I wall.
For this generator O_K=Z[B], N(B)=1 and N(4B-1)=3.

These are necessary conditions, not an exclusion or an existence
statement. The quadratic branch remains untouched. The B16+I row
and the global arithmetic ceiling seventeen remain open and unchanged.

## 1. Norm and trace leave a short integer problem

Write the three conjugates as a,u,v, with a the physical B. Then

    4<a<=441/16, 0<u,v<1/4.

The positive algebraic-integer norm is a positive rational integer and
satisfies

    0<N(B)=auv<(441/16)(1/4)^2=441/256<2.

Consequently auv=1. It follows from uv<1/16 that a>16.
The monic minimal polynomial therefore has the form

    P(X)=X^3-S X^2+T X-1,
    S=a+u+v, T=a(u+v)+uv,

where S,T are positive rational integers. Since

    16<S<441/16+1/2=449/16,

the complete trace range is 17<=S<=28. There are no missing
higher-degree or repeated-root cases: this proposition explicitly
assumes the degree-three branch of the inherited theorem, and number
fields in characteristic zero are separable.

## 2. The interval(0,1/4) determines both coefficients

Because u,v<1/4<a,

    P(1/4)=(1/4-a)(1/4-u)(1/4-v)<0.

Multiplication by64 gives `16T-4S-63<0`. Since 4T-S is an
integer, this is exactly

    4T-S<=15.                                      (2)

At the root u, strict u<1/4 also gives

    0=P(u)<-(S-1/4)u^2+T u-1.

The quadratic on the right is strictly positive at u and has negative
leading coefficient. Its discriminant must be strictly positive:
`T^2-4S+1>0`. The integer T squared minus4S is therefore nonnegative:

    T^2>=4S.                                      (3)

Equations(2)--(3) and 17<=S<=28 are now complete finite bounds.
They have precisely one solution:

| Trace range | Upper bound from(2) | Comparison with(3) |
|---|---:|---|
|17<=S<=20|T<=8|T squared<=64<68<=4S, impossible|
|21<=S<=24|T<=9|T squared<=81<84<=4S, impossible|
|25<=S<=28|T<=10|100>=T squared>=4S>=100; hence S25,T10|

This proves P=F. In particular the cubic necessary-condition example
in the inherited source is the ONLY possible cubic polynomial, rather
than one member of an unspecified list. No search computation is used.

## 3. Physical root and the old I contact

The values

    F(0)=-1, F(1/5)=1/125, F(1/4)=-3/64,
    F(24)=-337, F(25)=249

put two roots in(0,1/4) and one in(24,25). These are all the roots
of this cubic. The physical a>16 is therefore the last one, and a>24.

If beta were nonzero, its native finite contact would give beta>=1/2.
The inherited common short-edge upper would then imply

    a<=(441/16)(3/4)=1323/64<24,

contradicting a>24. Hence beta=0. This is orthogonality to the old
I wall; it makes no assertion about the short cap's contact with the
opposite wall across the I edge or with other original facets.

## 4. Exact field data and the remaining obligation

The rational-root test gives F(1)=-15 and F(-1)=-37, so F is
irreducible. Its polynomial discriminant is

    25^2*10^2-4*10^3-4*25^3-27+18*25*10
      =62500-4000-62500-27+4500=473=11*43.

The elementary discriminant-index identity says that this polynomial
discriminant equals the field discriminant times the square of
`[O_K:Z[B]]`. Since473 is squarefree, that index is1 and the field
discriminant is473. Finally

    N(B)=1,
    N(4B-1)=-64 F(1/4)=3.

Thus this entire cubic short-cap branch has a specified integral
basis, one physical generator and zero old-I contact. These data
could support a subsequent actual-wall or common-vector obstruction;
none is supplied here. The existence of the polynomial and its
admissible coarse conjugate strips does not realize a Coxeter cap,
let alone a compact polytope. This proposition neither changes the
seven arithmetic frontier groups nor closes the B16+I group.

**Verification status:** fully written manual argument, root is the
author; full outside-author review and root acceptance are pending.
There is no runtime, exact-enumeration or Lean credit from this note.
