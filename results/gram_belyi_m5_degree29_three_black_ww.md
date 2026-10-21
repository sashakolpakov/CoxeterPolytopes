# Exact three-black/`WW` no-go at degree 29

**Status (2026-08-10): proved for the complete `1,372`-record stratum;
global degree-`29` sharpness remains open.**

Put

\[
 \lambda=\frac{7+2\sqrt5}{29},\qquad
 29\lambda^2-14\lambda+1=0.
\]

The preceding low-vertex reduction showed that, after interchanging the two
finite Belyi values, any normalized degree-`29` Belyi polynomial with
`B(lambda)=0` has at least three black vertices, and that both normalized
source points are white when there are exactly three.  The present exact
calculation closes that last case:

> **Three-black theorem.**  There is no normalized degree-`29` Belyi
> polynomial whose marked black fibre contains `lambda`, has exactly three
> distinct points, and whose points `0,1` lie in the white fibre.

Consequently, any still-hypothetical degree-`29` witness has at least four
points in the finite fibre containing `lambda`.  This is not a complete
degree-`29` no-go: `42,605,329` marked ramification records with at least four
black vertices remain.

## 1. Why `1,372` passports are only `784` systems

Write the three black multiplicities as `(a,b,c)`, where `a` marks the root
`lambda`, `b<=c`, and

\[
                         a+b+c=29.
\tag{1}
\]

There are exactly

\[
             \sum_{a=1}^{27}\left\lfloor\frac{29-a}{2}\right\rfloor=196
\tag{2}
\]

marked black-degree records.  Riemann--Hurwitz leaves two possible white
partitions:

```text
(3,1^26)       or       (2,2,1^25).
```

For each black-degree record the ordered endpoint degrees give seven marked
passport records: `(1,1),(3,1),(1,3)` in the first partition and
`(1,1),(2,1),(1,2),(2,2)` in the second.

Several of these are the same coefficient system.  Requiring `x=0` to be a
non-root critical point covers both endpoint degrees `2` and `3`; the same is
true at `x=1`.  If neither endpoint is critical, the two white partitions
again have identical equations and are distinguished only by whether the
residual critical quadratic is separable.  Thus one black-degree record gives
only four systems:

| coefficient system | marked passports represented |
|---|---:|
| only `x=0` required critical | 2 |
| only `x=1` required critical | 2 |
| both endpoints critical | 1 |
| both endpoints simple | 2 |

Hence there are `588` endpoint-critical systems representing `980` records,
and `196` simple/simple systems representing `392` records: `784` systems and
`1,372` marked records in total.

## 2. Uniform elimination when an endpoint is critical

Let

\[
 F(x)=(x-\lambda)^a(x-r)^b(x-s)^c.
\tag{3}
\]

After removing the forced root factors from `F'`, its residual quadratic is

\[
 D(x)=29x^2-Ax+B,
\tag{4}
\]

where

\[
 \begin{aligned}
 A&=(b+c)\lambda+(a+c)r+(a+b)s,\\
 B&=ars+b\lambda s+c\lambda r.
 \end{aligned}
\tag{5}
\]

If both endpoints are critical, necessarily

\[
                         D(x)=29x(x-1).
\tag{6}
\]

The equation `A=29` solves linearly for `s`; after substitution, `B=0` is a
quadratic in `r`, while `F(0)=F(1)` has exact degree `29-a`.  A
degree-preserving reduction with gcd one therefore proves their
characteristic-zero resultant is nonzero.

If one endpoint `e` is critical, write

\[
                         D(e)=C(r)s+K(r).
\tag{7}
\]

Both `C` and `K` are linear.  The possible exceptional branch `C=K=0` is
excluded by its nonzero linear resultant.  On the remaining branch set
`s=-K/C`; the other root of `D` is `z=A/29-e`.  Clearing denominators in

\[
                    F(1-e)=F(e),\qquad F(z)=F(e)
\tag{8}
\]

gives two univariate equations of exact degrees `29-a` and `58`.  Full
finite-field polynomial gcds, with both degrees preserved, are one in all
`392` one-endpoint systems.  Together with the `196` both-endpoint systems,
this excludes all `588` endpoint-critical systems.  The witness-prime
histogram is

| split prime | systems separated |
|---:|---:|
| 101 | 525 |
| 109 | 63 |

This is a Sylvester-resultant certificate over the algebraic closure, not a
search for rational points over the displayed finite fields.

## 3. Uniform elimination when both endpoints are simple

Translate the marked black root to zero and scale the root of multiplicity
`b` to one.  Every three-black configuration is then represented by

\[
                         G(y)=y^a(y-1)^b(y-t)^c.
\tag{9}
\]

If the original endpoint coordinates become `y_0,y_1`, their ratio is fixed:

\[
 q=\frac{y_1}{y_0}=-\frac{1-\lambda}{\lambda}
   =29\lambda-13,qquad q^2+12q+16=0.
\tag{10}
\]

The residual derivative quadratic is

\[
 D_t(y)=29y^2-\bigl(a+c+(a+b)t\bigr)y+at.
\tag{11}
\]

Reduce `G` modulo `D_t`:

\[
                         G(y)\equiv R_1(t)y+R_0(t).
\tag{12}
\]

The two non-root critical points have one value if and only if
`S(t):=R_1(t)=0`; this statement includes the double-critical-point case.
Homogeneous quadratic division gives `deg S<=28`.  Every displayed
certificate retains degree `28`, so no leading term has been lost.

For a root of `S`, put

\[
                         H_t(y)=G(y)-R_0(t).
\tag{13}
\]

Two white endpoints with ratio `q` exist only if

\[
                         H_t(y)=H_t(qy)=0
\tag{14}
\]

for some `y`.  At a split prime, the certificate factors `S` completely.  For
every irreducible factor `f`, it works in the exact field
`F_p[t]/(f)` and computes the full polynomial gcd

\[
                  \gcd_y\bigl(H_t(y),H_t(qy)\bigr)=1.
\tag{15}
\]

All factors of all `196` equations pass.  Their witness-prime histogram is

| split prime | systems separated |
|---:|---:|
| 101 | 173 |
| 109 | 22 |
| 131 | 1 |

### Why the residue-field calculation is a proof

This step does not assume that a characteristic-zero solution has integral
coordinates without justification.  At the selected prime, all denominators
are units because `p!=29`, and the leading coefficient of `S` is a unit
because its reduction has full degree `28`.  Scaling makes `S` monic over the
local ring, so every root `t` is integral there.  The polynomial `H_t` is
monic of degree `29` in `y`, hence a common root `y` would also be integral.
Any characteristic-zero solution would therefore reduce to a solution over
the algebraic closure of the residue field.  Complete factorization of `S`
and the extension-field gcds in (15) show that no such reduction exists.

## 4. Reproducibility and scope

The exact producer and frozen archive are:

- [`certify_gram_belyi_m5_degree29_three_black_ww.py`](../scripts/certify_gram_belyi_m5_degree29_three_black_ww.py);
- [`gram_belyi_m5_degree29_three_black_ww.json`](../data/enumerations/gram_belyi_m5_degree29_three_black_ww.json).

A replay performs `784` bounded coefficient-system checks and no plane-tree
census.  The archive explicitly records

```text
status = EXACT_THREE_BLACK_WW_NO_GO
global_degree_29_status = OPEN
degree_29_witness_found = false
complete_degree_29_no_go_claimed = false
remaining_three_black_WW_marked_records = 0
remaining_degree_29_marked_records_with_at_least_four_black_vertices = 42605329
```

The next boundary is therefore the four-black stratum, not another pass over
the now-closed three-black passports.
