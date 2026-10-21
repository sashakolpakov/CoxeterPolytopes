# An exact Fekete probe on `[-2,1]`

**Status (2026-08-10): proved, with an exhaustive bounded-factor
subtheorem.**  This is a project-native auxiliary-polynomial argument.  It
does not import a field list, reflection-group classification, or a numerical
degree theorem from the literature.

The main output is the exact polynomial

\[
 P_6(x)=x(x-1)(x+1)(x+2)(x^2+x-1),
 \qquad \|P_6\|_{[-2,1]}={45\over64}<1.       \tag{1}
\]

It has the smallest possible degree among nonzero integer polynomials of
sup norm less than one on this interval, and it is unique in that degree up
to sign.  A finite search in a stated factor grammar then gives a degree-20
auxiliary polynomial which turns one distinguished conjugate in `[2,14]`
into the explicit field-degree bound `23`.

## 1. Minimal polynomial certificate

Put

\[
 y=x+\frac12,\qquad z=y^2.
\]

For `-2<=x<=1`, one has `0<=z<=9/4`, and

\[
 P_6=(z-\tfrac14)(z-\tfrac54)(z-\tfrac94).
\]

The endpoints give absolute values `45/64` and `0`.  The two non-root
critical points are

\[
 z={5\over4}\mathbin{\pm}{\sqrt3\over3},
\]

and both have absolute value `2sqrt(3)/9`.  Since

\[
 {4\over27}<\left({45\over64}\right)^2,
\]

equation (1) follows exactly.

> **Theorem 1 (minimality and uniqueness).**  No nonzero polynomial in
> `Z[x]` of degree at most five has sup norm less than one on `[-2,1]`.
> In degree six the only such polynomials are `P_6` and `-P_6`.

**Proof.**  If `P in Z[x]` has `||P||<1`, its values at the four integers
`-2,-1,0,1` are integers of absolute value less than one, hence vanish.
Therefore

\[
 A=x(x-1)(x+1)(x+2)
\]

divides `P` in `Z[x]`.  Let `r_+,r_-` be the roots of
`q=x^2+x-1`.  Then `A(r_+)=A(r_-)=-1` and
`|r_+-r_-|=sqrt(5)`.

For degree at most five write `P=A(ax+b)`.  The two strict inequalities

\[
 |ar_++b|<1,\qquad |ar_-+b|<1
\]

give `|a|sqrt(5)<2`, so the integer `a` is zero; then `|b|<1` gives
`b=0`.

For degree six write `P=A(ax^2+bx+c)`.  Modulo `q`, the quotient is

\[
 (b-a)x+(a+c).
\]

The same difference argument gives `b=a`, and then the common integer value
forces `c=-a`.  Thus `P=aP_6`.  Since `||P_6||=45/64`, strict norm less than
one permits exactly `a=0,+1,-1`; discarding the zero polynomial proves the
claim.  `square`

One immediate consequence is worth separating from the later degree bound.
If a totally real algebraic integer has **all** conjugates in `[-2,1]`, then
the integer norm of `P_6(alpha)` has absolute value less than one.  Hence it
vanishes and `alpha` is one of

\[
 -2,-1,0,1,{ -1+\sqrt5\over2},{-1-\sqrt5\over2}.       \tag{2}
\]

This conclusion is obtained from (1), not from a pre-existing classification
of small-span algebraic integers.

## 2. The bounded factor census

We exhaust the following explicitly bounded grammar:

\[
 x^a(x+1)^b(x+2)^c(x-1)^d(x^2+x-1)^e,
 \quad
 1\le a,b,c,d\le4,\quad0\le e\le4.             \tag{3}
\]

There are exactly `4^4*5=1,280` rows.  The unique winner for normalized sup
norm is

\[
 \begin{aligned}
 P_{20}(x)
   &=[x(x+1)]^4[(x-1)(x+2)]^3[x^2+x-1]^3,\\
 \deg P_{20}&=20.                                      \tag{4}
 \end{aligned}
\]

In the `z` coordinate,

\[
 P_{20}=(z-\tfrac14)^4(z-\tfrac94)^3(z-\tfrac54)^3.
\]

Its non-root critical points are `z=3/4` and `z=37/20`.  The complete list
of nonzero endpoint/critical magnitudes is

\[
 {91125\over1048576},\qquad
 {27\over1024},\qquad
 {884736\over9765625}.
\]

Consequently

\[
 \boxed{\|P_{20}\|_{[-2,1]}={884736\over9765625}},
 \qquad
 \|P_{20}\|^{1/20}=0.8868612567109771154\ldots .       \tag{5}
\]

The bounded-family optimality check does not promote floating-point
maximization.  For every one of the other `1,279` exponent tuples, the
checker isolates every root of its degree-at-most-five logarithmic-derivative
numerator in a rational interval of width at most `10^-8`, chooses a rational
midpoint, and verifies exactly

\[
 |P(r)|^{20}>\left({884736\over9765625}\right)^{\deg P}.
                                                               \tag{6}
\]

Thus one rational evaluation already proves that row's true normalized sup
norm is worse.  This proves optimality **only within (3)**.  No claim is made
that `P_20` is the global integer-Chebyshev polynomial of degree 20, or that
the grammar (3) contains all useful auxiliary factors.

## 3. One exceptional conjugate

> **Theorem 2 (exact degree-23 probe).**  Let `alpha` be a totally real
> algebraic integer, let `K=Q(alpha)`, and suppose one distinguished
> embedding satisfies
>
> \[
>             2\le alpha_+\le14,
> \]
>
> while every other embedding `sigma` satisfies
>
> \[
>             -2\le\sigma(alpha)\le1.
> \]
>
> Then
>
> \[
>                         \boxed{[K:Q]\le23}.             \tag{7}
> \]

**Proof.**  All roots of `P_20` lie in `[-2,1]`, so
`P_20(alpha)` is nonzero.  It is an algebraic integer.  On `[2,14]`, put
`w=x(x+1)>=6`; the function
`w^4(w-1)^3(w-2)^3` is increasing.  Hence

\[
 |P_{20}(alpha_+)|\le
 M:=210^4\,209^3\,208^3
 =159773975758768250880000.
\]

If `d=[K:Q]`, its nonzero integer norm gives

\[
 1\le |N_{K/Q}(P_{20}(alpha))|
    \le M\left({884736\over9765625}\right)^{d-1}.       \tag{8}
\]

Exact integer cross-multiplication gives

\[
 M\left({884736\over9765625}\right)^{22}>1,
 \qquad
 M\left({884736\over9765625}\right)^{23}<1.
\]

Thus `d-1<=22`, proving (7).  `square`

The generator hypothesis is essential.  The same computation bounds
`[Q(alpha):Q]`; it does not automatically bound an ambient field properly
containing `Q(alpha)`.

## 4. Translation to an AG Gram/dessin packet

This probe becomes an algebro-geometric field bound under a concrete marked
spectral condition.

Suppose an exact-return surface, stable-map, correspondence, or Gram--Belyi
packet returns its Gram field `k` through `t` marked **integral** spectral
coordinates `alpha_1,...,alpha_t`, and:

1. the `alpha_j` generate `k`;
2. for each `j`, its distinguished real value lies in `[2,14]` and all its
   other conjugates lie in `[-2,1]`.

Then Theorem 2 gives

\[
 [\mathbf Q(\alpha_j):\mathbf Q]\le23,
 \qquad
 [k:\mathbf Q]\le23^t.                                \tag{9}
\]

If one marked coordinate is primitive for `k`, the sharper bound is `23`.
This is particularly natural when an AG correspondence acts integrally on a
bounded cohomology lattice and the Gram coordinate is a marked eigenvalue.
For a Belyi cross-ratio, integrality and the conjugate envelope must be proved;
they do not follow from being Belyi.

Equation (9) makes the next graph problem genuinely finite without importing
a field catalogue.  An algebraic integer of degree at most `23` whose
conjugates have absolute value at most `14` has minimal-polynomial
coefficients bounded explicitly by

\[
 |c_{d-j}|\le {d\choose j}14^j.                         \tag{10}
\]

Thus a checker can enumerate the finite polynomial list, reconstruct the
finite marked Gram packets, and apply exact inertia, integrality, cyclic-
product, and Vinberg endpoint tests.  A uniform component cap `c` obtained
from that finite local Gram census would, in the **compact simple** lane and
under its separate narrow-facet coverage hypotheses, feed the elementary
surface-incidence inequality `n<=4c+1`.  The noncompact lane additionally
requires an affine-link/cusp argument and is not inferred from the compact
face average.  What remains open is the universal geometric input: proving
that general maximal arithmetic finite-volume chambers supply the integral
packet, the `[-2,1]` conjugate envelope, bounded `t`, and the required
compact/cusp-aware local coverage.

So the logical boundary is exact:

```text
bounded AG integral spectral grammar + [-2,1] envelope
        -> our finite polynomial/Gram census
        -> a dimension bound if the local component census closes;

failure of those hypotheses
        -> the surface/dessin decoration is the unbounded branch.
```

This replaces the broken use of an auxiliary interval bound wherever these
stated hypotheses hold.  It does not validate a printed equation, silently
identify `Q(alpha)` with a larger field, and does not assert a global Coxeter
theorem before the packet and cusp-coverage obligations are proved.

## 5. Reproducibility

- exact checker:
  [`scripts/certify_fekete_interval_minus2_1.py`](../scripts/certify_fekete_interval_minus2_1.py);
- canonical archive:
  [`data/enumerations/fekete_interval_minus2_1_certificate.json`](../data/enumerations/fekete_interval_minus2_1_certificate.json);
- focused tests:
  [`tests/test_fekete_interval_minus2_1_certificate.py`](../tests/test_fekete_interval_minus2_1_certificate.py).

Run

```text
python3 scripts/certify_fekete_interval_minus2_1.py --check
python3 -m pytest -q tests/test_fekete_interval_minus2_1_certificate.py
```
