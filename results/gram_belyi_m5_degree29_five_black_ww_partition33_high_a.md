# Exact high-marked-degree no-go in five-black `WW` partition `(3^2,1^23)`

**Status (2026-08-10): the `12` black-degree systems with marked black
degree `a>=21`, representing `48` marked passports, are impossible.  The
remaining `1,280` systems / `5,120` passports in this exact partition remain
open.**

This is a rigorously quantified strict sublocus theorem.  It is neither a
complete `(3^2,1^23)` theorem nor a complete five-black `WW` no-go.

## 1. Exact factorized system and coverage

For a canonical record `(a,b,c,d,e)`, write

\[
 G(x)=(x-\lambda)^a(x-u)^b(x-v)^c(x-w)^d(x-z)^e,
 \qquad b\le c\le d\le e,
\]

where

\[
 29\lambda^2-14\lambda+1=0.
\]

The white partition `(3^2,1^23)` has two critical white vertices `t,s`, each
of ramification degree three.  Thus the residual derivative quartic and
common-value equations are

\[
 D(x)=29(x-t)^2(x-s)^2,
\tag{1}
\]

\[
       G(s)=G(t),\qquad G(0)=G(t),\qquad G(1)=G(t).
\tag{2}
\]

Equating the four non-leading coefficients of (1) and adjoining (2) gives
seven equations in `(u,v,w,z,t,s)`, with fixed ordinary degrees

\[
                         (4,3,2,1,29,29,29).
\tag{3}
\]

The exact endpoint-degree choices in this partition are
`SS/SC/CS/CC=1/1/1/1`.  The same closed system covers all four: an endpoint
is critical when it coincides with `t` or `s`, and otherwise it is simple.
Consequently one black-degree system represents four marked passports.
Collision and diagonal `t=s` loci remain in the equations, so the system is
a necessary closed enlargement of the genuine passport stratum.

## 2. Uniform projective-infinity gate

Homogenize (1)--(2) in

\[
                  [u:v:w:z:t:s:h]\in\mathbf P^6.
\]

At `h=0`, the endpoint equations force `G_infinity(t)=0`, while the
critical-value equation then forces `G_infinity(s)=0`.  Hence both `t` and
`s` are black-root clusters.

For a cluster `alpha`, let `n_alpha` be its number of labelled roots and
`A_alpha` its total positive black degree.  As in the preceding partition
theorem,

\[
 D_\infty(x)=\prod_\alpha(x-\alpha)^{n_\alpha-1}E(x),
 \qquad
 E(\alpha)=A_\alpha\prod_{\beta\ne\alpha}(\alpha-\beta)\ne0
\tag{4}
\]

at every selected characteristic `p>29`.

If `t!=s`, equation (1) gives order two at each, so (4) requires two
disjoint clusters with `n_t=n_s=3`, impossible with only five labelled black
roots.  If `t=s`, the order is four, so all five roots collide there.  The
fixed marked root has become `0` at infinity, forcing every projective
coordinate to vanish.  This is again impossible.

Therefore the projective special fibre has no point at infinity for every
black-degree record, not merely for the certified high-`a` sublocus.

## 3. Exact bounded census

The direct six-variable standard-basis calculation is uneven on the full
`1,292` records, so no broad run was launched.  The predetermined sublocus

\[
                           a\ge21
\]

contains exactly `12` canonical records.  For every one, an exact Singular
standard basis over a good split finite field is `[1]`.  Together with the
uniform infinity theorem and projective properness, this excludes any
characteristic-zero solution.

The first successful prime histogram is:

| split prime | systems certified |
|---:|---:|
| 101 | 11 |
| 131 | 1 |

The calculation retains all collision components and performs no rational-
point sampling.  Interrupted preflights outside `a>=21` are not used as
evidence.

## 4. Exact residual

The theorem eliminates

\[
                        12\cdot4=48
\]

marked passports.  Since the same black-degree records occur in other exact
white partitions, it removes no entire broad status-free master system.  The
new boundary is

```text
partition_(3^2,1^23)_systems_eliminated = 12
partition_(3^2,1^23)_systems_remaining = 1,280
partition_(3^2,1^23)_passports_remaining = 5,120
five_black_status_free_systems_remaining = 9,742
five_black_WW_passports_remaining = 29,668
five_black_passports_remaining = 131,068
global_passports_eliminated = 107,322
global_passports_remaining = 42,503,356
partition_(3^2,1^23)_status = PARTIAL
five_black_WW_status = PARTIAL
global_degree_29_status = OPEN
```

The concrete next completion step is to replace the dense direct basis by a
regular-chain or tower decomposition of the symmetric quadratic factor

\[
 K=x^2-\sigma x+\pi,\qquad D=29K^2,
\]

and impose that `G mod K` is constant.  A monolithic four-variable basis of
that eliminated system was also unfavorable in preflight; it is not itself
a certificate or a no-go claim.

## 5. Reproducibility

The exact archive is produced by:

- [`certify_gram_belyi_m5_degree29_five_black_ww_partition33_high_a.py`](../scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition33_high_a.py);
- [`gram_belyi_m5_degree29_five_black_ww_partition33_high_a.json`](../data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_high_a.json).

The run is a fixed `12`-system coefficient census under the `a>=21`
condition.  It is not a plane-tree search, and no timeout or memory failure
is promoted as evidence.
