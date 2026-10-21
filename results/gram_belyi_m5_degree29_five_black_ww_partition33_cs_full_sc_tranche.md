# Exact `CS` theorem and bounded `SC` tranche in `(3^2,1^23)`

**Status (2026-08-10): the `CS` orientation is impossible for every black
degree record.  The predetermined `SC` tranche `17<=a_lambda<=20` is also
impossible.  The other `1,239` unresolved `SC` records, all `SS` records, the
full white partition, and degree `29` remain open.**

Endpoint status is ordered as `(x=0,x=1)`.  Reflection is not used:
`lambda -> 1-lambda` does not preserve the marked quadratic
`29 lambda^2-14 lambda+1`.

> **One-critical tower theorem.**  There is no degree-`29` Shabat polynomial
> in the five-black white partition `(3^2,1^23)` with endpoint status `CS`.
> In status `SC`, there is no such polynomial in the exact finite tranche
> `17<=a_lambda<=20`.

The complete `CS` statement combines the `1,280` rows certified here with
the `12` records `a_lambda>=21` already excluded for all four statuses by the
preceding high-degree theorem.

## 1. Separate critical factors

For black degrees `(a,b,c,d,e)`, let

\[
 G(x)=(x-\lambda)^a(x-u)^b(x-v)^c(x-w)^d(x-z)^e
\]

and let `D` be the residual derivative quartic.  In the two orientations,

\[
\begin{array}{lll}
 CS:&K=x(x-r),&D=29K^2,\quad G(r)=G(0),\quad G(1)=G(0),\\
 SC:&K=(x-1)(x-r),&D=29K^2,\quad G(r)=G(1),\quad G(0)=G(1).
\end{array}
\tag{1}
\]

The four coefficient equations in `D=29K^2` have structural degrees
`(4,3,2,1)` in `(r,u,v,w,z)`.  They define a one-dimensional derivative
base at every accepted split prime.  Adjoining the critical-value equation
gives an exact zero-dimensional algebra of dimension

\[
                         696=4\cdot3\cdot2\cdot1\cdot29.
\tag{2}
\]

Adjoining the remaining simple-endpoint value equation then gives the unit
ideal.  This is a successive finite-algebra calculation, not the stalled
monolithic six-variable system.

The frozen archive contains `1,321` exact rows:

| orientation | condition | certified rows | status |
|---|---|---:|---|
| `CS` | `a_lambda<=20` | 1,280 | complete after the preceding 12 high-`a` rows |
| `SC` | `17<=a_lambda<=20` | 41 | exact tranche only |

The first successful-prime histogram is `101: 1,288`, `109: 33`.

## 2. Uniform projective boundary

Homogenize the two fixed endpoints and `lambda` with `h`, retaining the
critical coordinate `r`.  At `h=0`, the fixed marked black root and the
fixed endpoints coalesce at `0`.  The critical-value equations force the two
roots of `K` to be black-root clusters.

If `r!=0`, the two order-two roots of `D_infinity` force disjoint clusters
of three labelled black roots each, impossible with only five labels.  If
`r=0`, the order-four root forces all five labels to collide at `0`, giving
the forbidden zero projective vector.  Thus the structural projective fibre
has no point at infinity for every selected prime `p>29`.  Properness lifts
the modular unit certificates to characteristic-zero nonexistence.

## 3. Exact ledger boundary

The preceding complete `CC` theorem left `3,840` marked passports: `1,280`
in each of `SS/SC/CS`.  The present certificate removes all `1,280` `CS`
rows and `41` `SC` rows.  Therefore

```text
partition_(3^2,1^23)_CS = IMPOSSIBLE
partition_(3^2,1^23)_SC_17<=a<=20 = IMPOSSIBLE
partition_(3^2,1^23)_SS_remaining = 1,280
partition_(3^2,1^23)_SC_remaining = 1,239
partition_(3^2,1^23)_passports_remaining = 2,519
five_black_status_free_systems_remaining = 9,742
five_black_WW_passports_remaining = 27,067
five_black_passports_remaining = 128,467
global_passports_eliminated = 109,923
global_passports_remaining = 42,500,755
global_degree_29_status = OPEN
```

The status-free necessary-system count does not fall: every surviving black
system still represents `SS`, unresolved `SC`, or another exact white
partition.

A broader live pass was deliberately interrupted after the strategic pivot.
Rows beyond the frozen `SC` tranche are not stored, counted, or promoted as
evidence.  No timeout or interruption is a mathematical claim.

## 4. Reproducibility

- producer:
  [`certify_gram_belyi_m5_degree29_five_black_ww_partition33_one_critical_towers.py`](../scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition33_one_critical_towers.py);
- archive:
  [`gram_belyi_m5_degree29_five_black_ww_partition33_one_critical_towers.json`](../data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_one_critical_towers.json).

Run with at most four orchestration workers:

```bash
python3 scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition33_one_critical_towers.py --check --workers 4
```

The heavy algebra is executed in isolated, sequential-per-row Singular
children.  The archive stores every black degree, orientation, split prime,
equation digest, exact intermediate dimensions, final unit basis, and the
uniform infinity gate.
