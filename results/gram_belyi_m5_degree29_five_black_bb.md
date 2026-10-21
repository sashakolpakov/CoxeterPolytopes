# Exact no-go for the degree-29 five-black `BB` block

**Status (2026-08-10): all `11,011` five-black `BB` master systems,
representing `55,055` marked passports, are impossible; the five-black and
global degree-`29` problems remain open.**

Let

\[
 \lambda=\frac{7+2\sqrt5}{29},\qquad
 29\lambda^2-14\lambda+1=0,
\]

and complement the Belyi polynomial if necessary so that `B(lambda)=0`.
The exact five-black quotient has `190,047` marked passports and `33,079`
binary endpoint-status master systems.  Its `BB` block consists of `11,011`
systems / `55,055` marked records.  The calculation below eliminates the
whole block.

> **Five-black `BB` theorem.**  No degree-`29` Shabat polynomial with exactly
> five distinct black roots can put the marked point `lambda` and both
> normalized source endpoints `0,1` in the black fibre.

This is a complete theorem for `BB`, not a complete five-black no-go.

## 1. Canonical critical equations

For a degree record `(a,b,c,d,e)`, the first three multiplicities belong to
the marked root and the ordered endpoints; the remaining two are unlabelled,
so `d<=e`.  Write

\[
 G(x)=(x-\lambda)^a x^b(x-1)^c(x-u)^d(x-v)^e,
 \qquad a+b+c+d+e=29.
\tag{1}
\]

There are exactly `11,011` such records.  Remove the forced black-root powers
from `G'`, leaving the residual quartic `D`, and write

\[
             G(x)\equiv R_3(u,v)x^3+R_2(u,v)x^2
                         +R_1(u,v)x+R_0(u,v)\pmod D.
\tag{2}
\]

Every non-root critical point has the same value if and only if

\[
                         R_1=R_2=R_3=0.
\tag{3}
\]

This scheme-theoretic formulation includes repeated residual critical
points and imposes no distinct-root inequations.  It therefore gives a
necessary closed system for every `BB` passport.

## 2. Fixed projective degree bounds

Before fixing roots, the coefficient `R_j` is homogeneous of degree `29-j`
in all five black roots.  For the possible degree-`28` part of `R_1` in only
`u,v`, set the other three roots to zero.  Then the residual derivative is
divisible by `x^2`, as is the remainder of `G`; its `x` coefficient vanishes.
Consequently

\[
                 \deg_{u,v}R_1\leq27,qquad
                 \deg_{u,v}R_2\leq27,qquad
                 \deg_{u,v}R_3\leq26.
\tag{4}
\]

Homogenize the three equations to the fixed degrees `(27,27,26)` in
projective coordinates `[u:v:w]`.  Using fixed structural bounds is
important: the argument below does not assume that a leading coefficient is
a unit after reduction.

## 3. Exact finite-field projective certificates

For each record, choose a split prime for `Q(sqrt(5))`, always different
from `2` and `29`.  Division by the leading coefficient of `D` introduces
only powers of `29`, so all three equations are integral at the selected
place.

The reduction is checked on the affine chart `w=1` by an exact standard
basis over `F_p[u,v]`.  For every record,

\[
                         (R_1,R_2,R_3)=(1).
\tag{5}
\]

Thus there is no affine point even over the algebraic closure of `F_p`.

At infinity, let `L_i(u,v)` be the fixed-degree leading homogeneous form of
`R_i`.  On the chart `u=1`, the exact univariate calculation gives

\[
                   \gcd\bigl(L_1(1,v),L_2(1,v),L_3(1,v)\bigr)=1.
\tag{6}
\]

The remaining point `[u:v]=[0:1]` is excluded because at least one
`L_i(0,1)` is nonzero.  Hence the three reduced homogenized equations have no
common point anywhere in `P^2`.

The first successful prime histogram is:

| split prime | systems certified |
|---:|---:|
| 101 | 10,589 |
| 109 | 401 |
| 131 | 21 |

## 4. Why this proves characteristic-zero nonexistence

Suppose (3) had a characteristic-zero solution.  Its affine point
`[u:v:1]` lies in the projective closure defined using (4).  Projective space
is proper, so after extending the selected valuation the point has an
integral projective representative and a nonzero reduction.  Because the
homogenized equations have integral coefficients, that reduction is a common
projective zero over the residue-field algebraic closure.  Equations (5) and
(6), together with the check at `[0:1]`, show that no such point exists.

This projective argument rules out escape at infinity directly; it does not
deduce nonexistence from an affine finite-field search or from an unproved
monicity assumption.

## 5. Exact quotient and remaining boundary

The independent five-black quotient is:

| endpoint block | marked passports | binary-status masters |
|---|---:|---:|
| `BB` | 55,055 | 11,011 |
| `BW` | 50,700 | 8,450 |
| `WB` | 50,700 | 8,450 |
| `WW` | 33,592 | 5,168 |
| **total** | **190,047** | **33,079** |

After the theorem, the precise boundary is

```text
five_black_BB_status = IMPOSSIBLE
five_black_status = PARTIAL
five_black_master_systems_remaining = 22,068
five_black_marked_passports_remaining = 134,992
global_marked_passports_eliminated = 103,398
global_marked_passports_remaining = 42,507,280
global_degree_29_status = OPEN
degree_29_witness_found = false
complete_five_black_no_go_claimed = false
complete_degree_29_no_go_claimed = false
```

The smallest shared black-degree family is `WW`: `1,292` black records
represent `5,168` binary-status systems / `33,592` marked passports.  The
algebraically overdetermined one-black/critical-white family has `8,450`
systems / `59,150` records across both orientations.  These are the next
uniform targets; neither is claimed solved here.

## 6. Reproducibility

The exact quotient and `BB` certificates are produced by:

- [`certify_gram_belyi_m5_degree29_five_black_quotient.py`](../scripts/certify_gram_belyi_m5_degree29_five_black_quotient.py);
- [`gram_belyi_m5_degree29_five_black_quotient.json`](../data/enumerations/gram_belyi_m5_degree29_five_black_quotient.json);
- [`certify_gram_belyi_m5_degree29_five_black_bb.py`](../scripts/certify_gram_belyi_m5_degree29_five_black_bb.py);
- [`gram_belyi_m5_degree29_five_black_bb.json`](../data/enumerations/gram_belyi_m5_degree29_five_black_bb.json).

The `BB` archive records the exact equation digest, structural degrees,
affine unit-basis result, infinity gcd, and missing-point evaluation for every
one of the `11,011` systems.  A replay is a bounded coefficient calculation;
it performs no plane-tree enumeration and uses no timeout as evidence.
