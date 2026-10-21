# Exact low-label calibration of the Gram--Belyi packet

**Status (2026-08-10): five exact normalized-polynomial minima, one assigned
orthogonality packet, and one rigorous open interval.**  The standard
`m=5` Coxeter relation is the arithmetic outlier: its Gram--Belyi leaf has
degree at least `29`.  This is a theorem, not a numerical guess.  An explicit
factored Belyi polynomial proves finiteness of the upper endpoint, but the
exact `m=5` minimum remains open.

This report calibrates the invariant defined in
[`gram_belyi_coxeter_dessin_theorem.md`](gram_belyi_coxeter_dessin_theorem.md).
For a nonorthogonal pair it uses

\[
 u=\frac{(r_i,r_j)^2}{r_i^2r_j^2},\qquad
 \lambda=\frac{u}{1+u},
\]

and `h_Bel(lambda)`, the least degree of a normalized Belyi polynomial
sending `lambda` to zero or one.  For a finite Coxeter relation,
`u=cos^2(pi/m)`.

## Certified table

| pair type | `u` | `lambda` | certified packet degree |
|---|---:|---:|---:|
| orthogonal, `m=2` | `0` | not intrinsically defined | assigned `2`; the auxiliary `lambda=1/2` itself has exact height `2` |
| `m=3` | `1/4` | `1/5` | **exactly `5`** |
| `m=4` | `1/2` | `1/3` | **exactly `3`** |
| `m=5` | `(3+sqrt(5))/8` | `(7+2sqrt(5))/29` | **at least `29`** and at most `3099548558797153417819718240423033019647072`; exact minimum open |
| `m=6` | `3/4` | `3/7` | **exactly `7`** |
| parallel | `1` | `1/2` | **exactly `2`** |
| dotted test `u=2` | `2` | `2/3` | **exactly `3`** |

Thus the packet degree is not monotone in the Coxeter label.  In particular,
the `m=5` coordinate is not a degree-five leaf in this polynomial
normalization.

## 1. Exact rational rows

The primary lower-bound input is Theorem 3 of Rodriguez,
[*Bounding the Degree of Belyi Polynomials*](https://arxiv.org/abs/1104.2027)
(Journal of Number Theory 133 (2013), DOI
[`10.1016/j.jnt.2012.12.019`](https://doi.org/10.1016/j.jnt.2012.12.019)).
If a nonzero algebraic number has nonzero `p`-adic valuation, its normalized
Belyi-polynomial height is at least `p`.

For `0<a<b`, the elementary polynomial

\[
 B_{a,b}(x)=
 \frac{b^b}{a^a(b-a)^{b-a}}x^a(1-x)^{b-a}
\tag{1.1}
\]

has degree `b`, sends `0` and `1` to zero, sends `a/b` to one, and has no
other finite critical value.  Consequently:

- `lambda=1/5` has lower bound `5` from `v_5(lambda)=-1`, and (1.1) has
  degree `5`;
- `lambda=1/3` has lower bound `3` and a degree-`3` witness;
- `lambda=3/7` has lower bound `7` and a degree-`7` witness;
- `lambda=1/2` has lower bound `2` and witness `4x(1-x)`;
- the dotted test `lambda=2/3` has lower bound `3` and witness
  `(27/4)x^2(1-x)`.

These are proved minima.  They are not merely explicit upper bounds.

For `m=2`, however, the general Gram--Belyi construction deliberately does
not apply `lambda=u/(1+u)` at `u=0`.  It records an orthogonality bit and the
fixed auxiliary degree-two packet.  The auxiliary number `1/2` happens to
have exact height two, but it must not be confused with a coordinate recovered
from the zero Gram entry.

## 2. Why the `m=5` leaf has degree at least 29

Exact algebra gives

\[
 u=\cos^2(\pi/5)=\frac{3+\sqrt5}{8},\qquad
 \lambda=\frac{u}{1+u}=\frac{7+2\sqrt5}{29},
\]

and

\[
                       29\lambda^2-14\lambda+1=0.
\tag{2.1}
\]

At the prime `29`, the coefficient points of (2.1) are

\[
                         (0,0),\ (1,0),\ (2,1).
\]

The lower Newton polygon therefore has slopes `0` and `1`; the two root
valuations are `0` and `-1`.  Choose the embedding into
`overline(Q)_29` for which `v_29(lambda)=-1`.  Any hypothetical normalized
Belyi polynomial and its equation `B(lambda) in {0,1}` survive under this
field embedding.  Rodriguez's theorem then forces

\[
                              \boxed{h_{\rm Bel}(\lambda)\ge 29}.
\tag{2.2}
\]

Equivalently, one may use Galois invariance of `h_Bel`: conjugating both the
polynomial and its marked preimage preserves degree and the normalization.

## 3. A finite, completely explicit upper endpoint for `m=5`

The following witness is intentionally kept as a composition.  Expanding it
would create no mathematical information.

Put

\[
 C(x)=\frac{29x^2-14x}{15},\qquad
 Q(y)=-(y-1)^{15}(15y+1).
\tag{3.1}
\]

For either root `lambda` of (2.1), the exact marked and critical values are

\[
\begin{array}{c|cccc}
x&0&1&\lambda&7/29\\ \hline
C(x)&0&1&-1/15&-49/435.
\end{array}
\tag{3.2}
\]

Moreover

\[
 Q'(y)=-240y(y-1)^{14},\quad
 Q(0)=1,\quad Q(1)=Q(-1/15)=0.
\tag{3.3}
\]

At the remaining critical value,

\[
 t=Q(-49/435)=-\frac ba,
\]

where

\[
 a=3^{15}5^{14}29^{16}
  =21916280382795726534721079037762451171875,
\]

\[
 b=2^{32}11^{30}
  =74944612079615317772145115975457330692096.
\tag{3.4}
\]

The integer `a` is odd.  Define

\[
 F(z)=-\left(\frac ab\right)^b(z-1)^a
                         \left(z+\frac ba\right)^b.
\tag{3.5}
\]

Then `F(0)=1` and `F(1)=F(-b/a)=0`.  Its derivative is, up to a nonzero
constant,

\[
 (a+b)z(z-1)^{a-1}(z+b/a)^{b-1},
\]

so `F` is itself a normalized Belyi polynomial.  The critical values created
by `C` map through `Q` to `t`, those created by `Q` map to `0` or `1`, and
those created by `F` already map to `0` or `1`.  Hence

\[
                              B(x)=F(Q(C(x)))
\tag{3.6}
\]

is a normalized Belyi polynomial, with

\[
 B(0)=0,\qquad B(1)=1,\qquad B(\lambda)=1.
\]

Its exact degree is

\[
 2\cdot16\cdot(a+b)
 =3099548558797153417819718240423033019647072.
\tag{3.7}
\]

This proves the stated finite interval.  It does **not** show that the upper
endpoint is remotely optimal, and it does **not** produce a degree-`29`
witness.  Determining whether (2.2) is sharp is a bounded low-degree dessin
problem left open here.

## 4. Immediate global consequence

Let `M=binom(f,2)`, and let `N_m` count facet pairs carrying the finite label
`m`.  Since every pair contributes at least two, the calibrated standard
relations sharpen the universal estimate to

\[
 \beta_{\rm GB}(P)\ge
 2M+3N_3+N_4+27N_5+5N_6,
\tag{4.1}
\]

before adding any further excess forced by other dotted coordinates.  A
rational dotted entry `u=2` contributes one above the baseline.  If all pairs
belong to the calibrated exact rows and there is no `m=5` pair, the analogous
sum with the dotted count is an equality.  Formula (4.1) is useful precisely
because an `m=5` edge is arithmetically much more expensive than the Coxeter
diagram alone suggests.

## Reproducibility and claim boundary

The exact producer and frozen archive are:

- [`certify_gram_belyi_standard_packets.py`](../scripts/certify_gram_belyi_standard_packets.py);
- [`gram_belyi_standard_packets.json`](../data/enumerations/gram_belyi_standard_packets.json).

The producer checks every rational witness, the `29`-adic Newton polygon, all
values in the factored `m=5` composition, the reduction `t=-b/a`, the final
degree, and the exact status labels.  It performs no broad dessin enumeration,
uses no floating-point recognition, and makes no claim about the exact `m=5`
minimum.
