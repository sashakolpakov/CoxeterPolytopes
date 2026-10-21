# Exact no-go for the five-black `WW` partition `(5,1^24)`

**Status (2026-08-10): all `1,292` black-degree systems for the white
partition `(5,1^24)`, representing `3,876` marked passports, are
impossible.  The five-black `WW`, five-black, and global degree-`29` problems
remain open.**

Let

\[
 \lambda=\frac{7+2\sqrt5}{29},\qquad
 29\lambda^2-14\lambda+1=0,
\]

and complement the Shabat polynomial if necessary so that
`G(lambda)=0`.  The preceding exact quotient and `BB` theorem leave `5,168`
binary-status `WW` masters / `33,592` `WW` marked passports in the
five-black stratum.  This calculation eliminates the complete exact white
partition `(5,1^24)` inside that block.

> **Five-black `WW` partition `(5,1^24)` theorem.**  No degree-`29` Shabat
> polynomial with exactly five black roots, one white vertex of degree `5`,
> and `24` simple white vertices can put both normalized source endpoints in
> the white fibre while the marked point `lambda` has the required value.

This includes all three possible ordered endpoint markings `SS`, `SC`, and
`CS`.  It is not a complete five-black `WW` no-go.

## 1. Why all white-endpoint statuses collapse

Normalize the marked black root to `0` and one other black root to `1`, and
write

\[
 G(y)=y^a(y-1)^b(y-u)^c(y-v)^d(y-w)^e.
\tag{1}
\]

Let `D` be the residual derivative quartic after removing the forced
black-root powers from `G'`, and divide

\[
              G=DQ+R,
       \qquad \deg Q=25,\quad \deg R\leq3.
\tag{2}
\]

The Shabat common-critical-value condition says that `R` is constant.  If
`alpha` is a simple white vertex, then `D(alpha)` is nonzero and
`G(alpha)-R=0`, so `Q(alpha)=0`.  If it is critical of ramification index
`r`, then

\[
 \operatorname{ord}_\alpha(G-R)=r,
 \qquad
 \operatorname{ord}_\alpha D=r-1,
\]

and again `Q(alpha)=0`.  Therefore every fixed white endpoint satisfies the
same quotient equation whether it is simple or critical.  This applies not
only to both endpoints in `WW`, but also to the single white endpoint in
`BW` and `WB`.  The necessary-system quotient is

| endpoint block | binary-status systems | status-free systems |
|---|---:|---:|
| `BW` | 8,450 | 4,225 |
| `WB` | 8,450 | 4,225 |
| `WW` | 5,168 | 1,292 |
| **total after `BB`** | **22,068** | **9,742** |

In particular, the `WW` row collapses from

```text
5,168 binary-status systems = 1,292 black records * 4 statuses
```

to `1,292` black-degree systems at the necessary-equation level.  Direct
standard-basis preflights of the full common-value system were not favorable;
no timeout is used as a mathematical result.

## 2. The tractable exact partition

For a canonical black degree record `(a,b,c,d,e)`, with the final four
entries nondecreasing, use the original endpoint coordinate and write

\[
 G(x)=(x-\lambda)^a(x-u)^b(x-v)^c(x-w)^d(x-z)^e,
 \qquad a+b+c+d+e=29.
\tag{3}
\]

There are exactly `1,292` such records.  The white partition `(5,1^24)`
forces the residual derivative to have one root of multiplicity four:

\[
 G'=\prod_{r\in\{\lambda,u,v,w,z\}}(x-r)^{a_r-1}D(x),
 \qquad D(x)=29(x-t)^4.
\tag{4}
\]

Equating the four non-leading coefficients in (4), and requiring both
normalized endpoints to be white, gives the closed six-equation system

\[
       D=29(x-t)^4,
 \qquad G(0)=G(t),
 \qquad G(1)=G(t)
\tag{5}
\]

in `(u,v,w,z,t)`.  Allowing `t=0` or `t=1` includes the `CS` and `SC`
markings; `t` different from both endpoints is the `SS` marking.  Thus one
system covers three marked passports per black record.  Collision and
degeneracy loci are retained, making (5) a necessary closed system rather
than an inequation-sensitive search.

## 3. Exact affine certificates

For each record, the script selects a split prime for `Q(sqrt(5))`, always
larger than `29`.  All coefficients are integral there.  Over
`F_p[u,v,w,z,t]`, Singular computes an exact standard basis of (5).  Every
one of the `1,292` ideals is the unit ideal.

The first successful split-prime histogram is:

| split prime | systems certified |
|---:|---:|
| 101 | 1,263 |
| 109 | 28 |
| 131 | 1 |

The archive records the prime, reduced value of `lambda`, exact equation
digest, fixed projective degrees, and unit-basis result for every row.  This
is algebraic-closure nonexistence in the affine special fibre, not a search
through rational points.

## 4. Uniform projective-infinity theorem

The four coefficient equations in (4) have fixed degrees `(4,3,2,1)`, and
the two endpoint equations have degree `29`.  Homogenize them in

\[
                    [u:v:w:z:t:s]\in\mathbf P^5,
\]

so the marked root becomes `lambda*s`.  At `s=0`, both endpoint equations
have the same leading form `-G_infinity(t)`.  Hence `t` must be a black-root
cluster among `0,u,v,w,z`.

Group collided labelled roots by their distinct values `alpha`.  Let
`n_alpha` be the number of labels in a cluster and `A_alpha` the sum of their
positive black degrees.  The residual derivative factors as

\[
 D_\infty(x)=
 \prod_\alpha (x-\alpha)^{n_\alpha-1}E(x),
\tag{6}
\]

with

\[
 E(\alpha)=A_\alpha
    \prod_{\beta\ne\alpha}(\alpha-\beta)\ne0.
\tag{7}
\]

The last inequality holds in every selected characteristic because
`1<=A_alpha<=29<p`.  Thus the order of `D_infinity` at the black root `t` is
exactly `n_t-1`.  Equation (4) gives order four, so `n_t=5`: all five labelled
roots collide at `t`.  The fixed root at infinity is `0`, forcing

```text
u = v = w = z = t = 0,
```

which is not a projective point.  Hence the special fibre has no point at
infinity for any record.  Properness now upgrades each affine finite-field
unit certificate to characteristic-zero nonexistence.

An independent audit also recomputed the representative
`p=101`, `(a,b,c,d,e)=(1,1,1,1,25)` system: the affine ideal and all five
infinity charts are unit ideals.

## 5. Exact count boundary

The partition has three ordered endpoint markings per black record:

| marking | exact degree choice count |
|---|---:|
| `SS` | 1 |
| `SC` | 1 |
| `CS` | 1 |
| `CC` | 0 |

Therefore the theorem eliminates

\[
                         1,292\cdot3=3,876
\]

marked passports.  It does not remove any entire binary-status master,
because each such master also represents other exact white partitions.  The
status-free quotient therefore reduces the number of distinct unresolved
necessary coefficient systems to

\[
                 4,225+4,225+1,292=9,742.
\]

The old count `22,068` remains useful as the conservative binary-status
ledger; the two counts have different meanings.  The new exact boundary is

```text
five_black_WW_partition_(5,1^24) = IMPOSSIBLE
five_black_WW_marked_passports_remaining = 29,716
five_black_binary_status_masters_remaining = 22,068  # conservative ledger
five_black_status_free_necessary_systems_remaining = 9,742
five_black_marked_passports_remaining = 131,116
global_marked_passports_eliminated = 107,274
global_marked_passports_remaining = 42,503,404
five_black_WW_status = PARTIAL
five_black_status = PARTIAL
global_degree_29_status = OPEN
degree_29_witness_found = false
```

The remaining exact `WW` rows are:

| white partition | marked passports remaining |
|---|---:|
| `(4,2,1^23)` | 9,044 |
| `(3^2,1^23)` | 5,168 |
| `(3,2^2,1^22)` | 10,336 |
| `(2^4,1^21)` | 5,168 |
| **total** | **29,716** |

The next smallest exact residual partitions are `(3^2,1^23)` and
`(2^4,1^21)`, each with `5,168` marked passports.  Their systems must also
enforce equality among all distinct critical values; they are not claimed
solved here.

## 6. Reproducibility

The certificate is produced by:

- [`certify_gram_belyi_m5_degree29_five_black_ww_partition5.py`](../scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition5.py);
- [`gram_belyi_m5_degree29_five_black_ww_partition5.json`](../data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition5.json).

The run is a finite census of `1,292` coefficient systems.  It performs no
plane-tree enumeration and promotes neither a timeout nor a numerical
calculation as evidence.
