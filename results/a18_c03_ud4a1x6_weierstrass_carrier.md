# Positive theorem: an equation-level `A18 C03 H11/28` K3 carrier

There is an explicit eight-moduli family of Jacobian elliptic K3 surfaces
whose very-general member has

\[
 \operatorname{NS}(X)=U+D_4(-1)+A_1(-1)^6
\]

and fibers `I0*+6I2+6I1`.  An exact determinant-minus-one marking identifies
its fiber-divisor basis with the scaled-dual Picard lattice and transports all
28 oriented walls and reflection matrices of the known finite-volume
noncompact `H11/28` Coxeter polytope.

## The eight-parameter equation

Take parameters `u2,u3,u4,s1,s2,s3,s5,s6`, with `s6 != 0`, and put

```text
u = 1 + 2t + u2*t^2 + u3*t^3 + u4*t^4 + t^5,
v2 = (-s2*s6^3 - 2*s5 + 3*s6^3*u2 + 3*s6*u4)/s6^3,
s4 = -(2*s1*s6^4 - s2*s6^3 - 2*s5 - 6*s6^4
       + 3*s6^3*u2 - 3*s6^2*u3 + 3*s6*u4)/s6,
S = 1 + s1*t + s2*t^2 + s3*t^3 + s4*t^4 + s5*t^5 + s6*t^6,
v = 3 + (12-3*s1)*t + v2*t^2
    + 3*(-s5+2*s6*u4)/s6^2*t^3 + 3/s6*t^4,
w = 1 + (6-2*s1)*t + (s2-3*u2+v2)*t^2 + t^3/s6^2,
A = -3*u^2 + S*v,
B = 2*u^3 - S*u*v + S^2*w.
```

Then

\[
 y^2=x^3+A(t)x+B(t).
\]

Here `deg(A)<=8`, `deg(B)<=12`, `t^2|A`, and `t^3|B`.  Direct expansion in
four universal indeterminates gives

\[
4A^3+27B^2=S^2\{-9u^2v^2+108u^3w
 +S(4v^3-54uvw)+27S^2w^2\}.
\]

The bracket is divisible by `t^6`; on the certified open it is `t^6 H(t)`
with `deg(H)=6`.  Thus zero is `I0*`, the six roots of `S` are `I2`, the six
roots of `H` are `I1`, and infinity is smooth.

## A nonempty exact witness

For

```text
(u2,u3,u4,s1,s2,s3,s5,s6)=(-1,1,1,0,2,0,0,1)
```

one obtains

```text
u = t^5+t^4+t^3-t^2+2t+1
S = t^6+11t^4+2t^2+1
v = 3t^4+6t^3-2t^2+12t+3
w = t^3+3t^2+6t+1
A = 22t^8+78t^7-22t^6+132t^5+11t^4+36t^3-2t^2
B = -37t^12+173t^11+222t^10+687t^9-125t^8+346t^7
    -149t^6+17t^5-53t^4-6t^3
H = 79555t^6+107370t^5+92717t^4+97974t^3
    +35999t^2+18900t+940.
```

The certificate records nonzero discriminants of `S` and `H`, and nonzero
resultants `Res(S,H)`, `Res(A,H)`, and `Res(u,S)`.  These prove that the stated
Kodaira configuration occurs on a nonempty Zariski open.

This construction avoids the misleading global-two-torsion ansatz.  At the
witness, a two-torsion point would give a polynomial root of degree at most
four of `x^3+A(t)x+B(t)`.  The exact coefficient ideal in the five possible
root coefficients has Groebner basis `[1]`: there is no such point.

## Why the Picard lattice is exact

The coefficient map has tangent rank eight at the witness.  The infinitesimal
`PGL2` and Weierstrass-homothety orbit has rank four, and the combined tangent
matrix has rank twelve; one displayed `12 x 12` minor is

```text
-37705049539909728.
```

Hence the moduli image really has dimension eight.  The ten reducible-fiber
roots plus fiber and section give the rank-12 trivial lattice.  Since the
moduli dimension equals `20-12`, the very-general Picard rank is exactly 12,
and Shioda--Tate gives Mordell--Weil rank zero.  All component groups are
2-primary; nontrivial torsion would contain two-torsion.  The unit ideal
witness and local constancy of the two-torsion monodromy type exclude it on
the very-general open.  Therefore the trivial lattice is saturated and is
the full Neron--Severi lattice.

## Literal attachment to the Coxeter chamber

In the fiber basis

```text
[F, F+O, C1, C2-central, C3, C4, E1, ..., E6]
```

the Gram matrix is exactly `U+D4(-1)+A1(-1)^6`.  The certificate gives an
integral matrix `Q`, with `det(Q)=-1`, satisfying

```text
Q^t Gram(N) Q = Gram(U+D4(-1)+A1(-1)^6).
```

All 28 transported roots are integral in this basis; their types are

```text
8*(-2,div1) + 12*(-2,div2) + 8*(-4,div2).
```

Their Hermite normal form is the identity, so they generate the full Picard
lattice.  The transported strict interior has square 608 and minimum wall
pairing 2.  The inherited literal Vinberg audit has 1,252 finite vertices,
64 ideal vertices, 8,851 elliptic edges, 186 connected affine components,
and zero edge or parabolic-neighborhood failures.  Thus the 28-wall chamber
is finite-volume noncompact.

## Scope

This is a new project-level equation realization and exact wall attachment of
a known abstract Coxeter polytope.  The imported chamber is not the complete nef
chamber as far as this certificate asserts.  The imported `(-2)` mirrors are not all asserted
to be irreducible curves, the square-608 interior is not asserted ample, and
no surface-automorphism or worldwide-priority claim is made.

Artifacts:

- [generator](../scripts/certify_a18_c03_ud4a1x6_weierstrass_carrier.py)
- [exact JSON certificate](../data/enumerations/a18_c03_ud4a1x6_weierstrass_carrier.json)
- [focused regression test](../tests/test_a18_c03_ud4a1x6_weierstrass_carrier.py)
