# Exact `CC` no-go in five-black `WW` partition `(3^2,1^23)`

**Status (2026-08-10): all `1,292` `CC` systems in the exact white
partition `(3^2,1^23)` are impossible.  After the overlap with the preceding
high-marked-degree theorem, this eliminates `1,280` new marked passports.
The `SS/SC/CS` statuses, the full partition, and degree `29` remain open.**

> **Partition-`33` `CC` theorem.**  No degree-`29` Shabat polynomial with
> five black vertices and white partition `(3^2,1^23)` can put both
> normalized source endpoints at the two critical white vertices.

This is a complete theorem for the exact `CC` status, not for the complete
white partition.

## 1. Fixed quadratic factor

For a canonical black-degree record `(a,b,c,d,e)`, write

\[
 G(x)=(x-\lambda)^a(x-u)^b(x-v)^c(x-w)^d(x-z)^e,
 \qquad b\le c\le d\le e.
\tag{1}
\]

In the `CC` status, the two degree-three white vertices are precisely the
normalized endpoints `0,1`.  Therefore their quadratic factor is fixed:

\[
 K(x)=x(x-1),
 \qquad
 D(x)=29K(x)^2=29x^2(x-1)^2,
\tag{2}
\]

where `D` is the residual derivative quartic.  The two critical values must
also agree:

\[
                         G(0)=G(1).
\tag{3}
\]

The four non-leading coefficient equations in (2) have degrees
`(4,3,2,1)` in `(u,v,w,z)`.  Equation (3) is then a single value test.  No
critical-root variables and no degree-29 quotient expansion are needed.

## 2. Successive exact finite algebra

At a good split prime, let `I_base` be the ideal of the four equations in
(2).  Singular first computes its exact standard basis.  For every one of
the `1,292` canonical black-degree records,

\[
             \dim I_{\rm base}=0,
 \qquad
             \dim_{\mathbf F_p}\mathbf F_p[u,v,w,z]/I_{\rm base}=24.
\tag{4}
\]

Thus the nominal Bezout number `4*3*2*1=24` is attained exactly.  The value
equation (3) is adjoined only after entering this 24-dimensional finite
algebra.  In every row the resulting standard basis is `[1]`; equivalently,
no point of the exact base tower satisfies the common-value equation.

The first successful prime histogram is:

| split prime | systems certified |
|---:|---:|
| 101 | 1,277 |
| 109 | 15 |

The predetermined preflight tranche `17<=a<=20` contained `41` records and
closed before the full run was launched.  The identical fixed-degree engine
then replayed all `1,292` rows in about `22` seconds.  This is a successive
finite-algebra computation, not a monolithic general `SS/SC/CS/CC` basis.

## 3. Projective boundary

The structural homogenization degrees are `(4,3,2,1,29)`.  Retain the two
critical coordinates while homogenizing the fixed conditions as

\[
                         t=0,\qquad s=h.
\]

At `h=0`, they collide at `0`, and (2) becomes

\[
                         D_\infty(x)=29x^4.
\]

For collided black-root clusters,

\[
 D_\infty(x)=\prod_\alpha(x-\alpha)^{n_\alpha-1}E(x),
 \qquad
 E(\alpha)=A_\alpha\prod_{\beta\ne\alpha}(\alpha-\beta)\ne0
\]

in characteristic `p>29`.  The order-four zero at `0` forces all five
labelled black roots to collide with the fixed marked root there, leaving
only the zero projective vector.  Hence the projective special fibre has no
point at infinity.  Properness upgrades every modular unit certificate to
characteristic-zero nonexistence.

## 4. Exact boundary

There is one `CC` marked passport per black-degree record.  The preceding
high-`a` theorem already removed all four statuses for `12` records, so the
new contribution is

\[
                         1,292-12=1,280.
\]

The exact live boundary becomes

```text
partition_(3^2,1^23)_CC = IMPOSSIBLE
partition_(3^2,1^23)_black_records_remaining = 1,280
partition_(3^2,1^23)_passports_remaining = 3,840  # SS/SC/CS
five_black_status_free_systems_remaining = 9,742
five_black_WW_passports_remaining = 28,388
five_black_passports_remaining = 129,788
global_passports_eliminated = 108,602
global_passports_remaining = 42,502,076
partition_(3^2,1^23)_status = PARTIAL
five_black_WW_status = PARTIAL
global_degree_29_status = OPEN
```

The status-free necessary-system count does not decrease: each black-degree
system still represents the other endpoint statuses or other exact white
partitions.  The raw marked-passport boundary does decrease exactly as shown.

The next algebraic targets are the one-critical-endpoint statuses `SC/CS`.
There `K=(x-1)(x-r)` or `K=x(x-r)`, so the base tower gains only one critical
parameter; it should be attacked before the general `SS` system.

## 5. Reproducibility

The exact certificate is produced by:

- [`certify_gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.py`](../scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.py);
- [`gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.json`](../data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.json).

The archive stores the degree-24 base-algebra check, exact endpoint unit
basis, equation digest, split prime, and projective-boundary gate for every
row.  No timeout, numerical root test, or plane-tree enumeration is promoted
as evidence.
