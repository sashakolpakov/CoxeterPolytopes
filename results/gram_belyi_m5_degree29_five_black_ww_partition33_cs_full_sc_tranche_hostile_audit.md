# Hostile audit: five-black `(3^2,1^23)` `CS` full and `SC` tranche

**Verdict (2026-08-10): PASS.** No mathematical, implementation, scope,
resource, ledger, or claim-boundary blocker was found. The frozen artifacts
certify exactly the `1,280` previously unresolved `CS` records and exactly the
`41` `SC` records with `17<=a_lambda<=20`. Together with the preceding
12-record all-status theorem, `CS` is complete. The other `1,239` `SC`
records, all `SS` records, the full partition, and degree `29` remain open.

## Frozen artifacts

| artifact | SHA-256 |
|---|---|
| `scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition33_one_critical_towers.py` | `374ad18c61b3c183ffeacdf4a5bce59a2fe8416e976365e433f856ee434dd361` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_one_critical_towers.json` | `53063a01bd7b83dfc1df35d2f9384c82eeda02d4704b812ec2415d16fd44a797` |
| `results/gram_belyi_m5_degree29_five_black_ww_partition33_cs_full_sc_tranche.md` | `bf729dc5c97d09973d7c5796383e44e3e8bed723a37d99dfcb84ec908a8d4030` |
| `tests/test_gram_belyi_m5_degree29_five_black_ww_partition33_one_critical_towers.py` | `528ffbf4d9e27964bff6b8ae3b2025e8eb5589e14234405941230a0cc8c7bc76` |

Direct dependencies were pinned as follows:

| dependency | SHA-256 |
|---|---|
| `scripts/certify_gram_belyi_m5_degree29_reduction.py` | `5bd3eae89b0b45f4f8b279fd1298ba0c0598d2f154cde9b26dd1f5d778132837` |
| `scripts/certify_gram_belyi_m5_degree29_five_black_quotient.py` | `deaaed9f90e550ecd8a2256dc65a2bdc4ca57aea7e18354b7dabfff0d344936f` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_quotient.json` | `89e2d1fafc10ac5ad91edd6e20ed5e96f6c238b404b8bb303a054f330cbc0e87` |
| `scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition5.py` | `a8c3f46b126424a80f02f0d4aab11f15d2fe9425cbbeed8236a861f0a966ab10` |
| `scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition33_high_a.py` | `f0cd32af85947d626e2ad74e193dc56d6a414ad72937abe62342f1330dc2ec36` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_high_a.json` | `011c3a48880f43fb282ce0bd4416ec875c3f6f7dc4b8ffb96909dc1b46d60ba7` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.json` | `6bc38604ce8c64519be1b71debdf48fddd9083616aee9409200334958e9a3630` |

The canonical compact `.rows` serialization, including the final newline,
has SHA-256
`8a51575e6b1c7dd08d9cecec845fc4241625ec2a1462885fb82b8a776337174e`.

## 1. Orientations and necessary equations

Endpoint status is ordered as `(x=0,x=1)`. For

\[
G(x)=(x-\lambda)^a(x-u)^b(x-v)^c(x-w)^d(x-z)^e,
\]

let `D` be the residual derivative quartic, of leading coefficient `29`.
In `CS`, zero is critical and one is simple, while in `SC` zero is simple
and one is critical. Thus the necessary systems are separately

\[
\begin{array}{lll}
CS:&K=x(x-r),&D=29K^2,\quad G(r)=G(0),\quad G(1)=G(0),\\
SC:&K=(x-1)(x-r),&D=29K^2,\quad G(r)=G(1),\quad G(0)=G(1).
\end{array}
\]

Independent expansion gives, from leading to constant coefficient,

```text
CS: 29, -58r, 29r^2, 0, 0
SC: 29, -58(1+r), 29(1+4r+r^2), -58r(1+r), 29r^2.
```

These are exactly the producer equations. The common critical value and
remaining simple endpoint value are also adjoined in the correct order.
Root collisions and forbidden endpoint coincidences are retained, enlarging
the necessary system; excluding it is therefore safe.

Reflection cannot identify the orientations. If lambda satisfies
`29lambda^2-14lambda+1=0`, substituting `1-lambda` leaves
`15(1-2lambda)`, not zero. The producer correctly sets
`reflection_symmetry_used=false`.

All `1,321` stored equation hashes were independently reconstructed and
matched; all are distinct. Every stored residue satisfies

\[
s^2=5\pmod p,
\qquad 29\lambda^2-14\lambda+1=0\pmod p,
\]

and the only selected primes, `101` and `109`, are split and exceed `29`.

## 2. Exact scope and row order

An enumeration independent of the producer formed every positive tuple

\[
(a,b,c,d,e),\qquad a+b+c+d+e=29,\qquad b\le c\le d\le e.
\]

It gives `1,292` records: exactly `12` have `a>=21`, exactly `1,280` have
`a<=20`, and exactly `41` have `17<=a<=20`. The archive is precisely

```text
rows 1 ... 1280: CS, all canonical a<=20 records in order
rows 1281 ... 1321: SC, exactly the 41 canonical 17<=a<=20 records in order.
```

The 12 high-`a` records agree exactly with the prior all-status archive. The
producer constructs its normal task list as those `1,280` `CS` tasks followed
by those `41` `SC` tasks. `ThreadPoolExecutor.map` preserves that input order.
No default or `--check` path schedules any other `SC` record. The transient
`--tranche-only` path instead runs both orientations on the same 41 records
and cannot be combined with `--check`; it does not broaden the theorem.

## 3. Successive dimension-696 algebra

The four derivative equations have structural degrees `(4,3,2,1)` in
`(r,u,v,w,z)` and define a one-dimensional base. The first value equation has
exact degree `29` because `G(r)` has leading term `r^29`. For every archived
row Singular recomputed

```text
DERIVATIVE_BASE
1
CRITICAL_VALUE_BASE
0
696
UNIT
```

Thus adjoining the critical-value equation gives a zero-dimensional affine
algebra of exact vector-space dimension

\[
696=4\cdot3\cdot2\cdot1\cdot29.
\]

Adjoining the remaining simple-endpoint value equation makes the standard
basis `[1]`. Computing `std(B+ideal(f4))` and then
`std(C+ideal(f5))` is ideal-theoretically identical to adjoining the two
equations successively; using standard bases as generators loses no
components or nilpotents. The final unit ideal excludes every affine point
over the algebraic closure, not only rational points.

## 4. Uniform projective gate and characteristic zero

Homogenize in `[r:u:v:w:z:h]`, replacing the marked root by `lambda*h` and
the fixed endpoint one by `h`. The structural degrees are
`(4,3,2,1,29,29)`. At `h=0`, both orientations have

\[
D_\infty(x)=29x^2(x-r)^2.
\]

The critical-value equation becomes `G_infinity(r)=0`; hence `r` is a
black-root cluster, while zero is already a cluster because the marked root
specializes there. For a collision cluster `alpha` containing `n_alpha`
of the five labelled black roots and total positive degree `A_alpha`,

\[
D_\infty(x)=\prod_\alpha(x-\alpha)^{n_\alpha-1}E(x),
\qquad
E(\alpha)=A_\alpha\prod_{\beta\ne\alpha}(\alpha-\beta)\ne0
\]

because `1<=A_alpha<=29<p`.

- If `r!=0`, the order-two zeros at zero and `r` force two disjoint
  three-label clusters, requiring six labels.
- If `r=0`, the order-four zero forces all five labels to collide at zero;
  then `r=u=v=w=z=h=0`, the forbidden zero projective vector.

These cases exhaust infinity, including the diagonal and every collision.
As a computational cross-check, the independently homogenized first five
equations were unit on all five ordinary infinity charts
`r,u,v,w,z=1` for both orientations and the representative records
`(1,1,1,1,25)` and `(20,2,2,2,3)` at `p=101`.

At the selected primes, `29` is a unit and lambda is integral in the chosen
localization of `Q(sqrt(5))`. The homogenized equations define a closed
projective scheme over that valuation ring. Properness would specialize any
characteristic-zero point, after finite extension, to either an affine
special-fibre point or a point at infinity; both have been excluded. Hence
the modular certificates rigorously imply characteristic-zero nonexistence.

## 5. Ledger and claim boundary

The preceding `CC` theorem left `3,840` marked passports, `1,280` in each of
`SS`, `SC`, and `CS`. The present rows are new and disjoint from the 12
previous high-`a` records, so the exact update is

```text
newly eliminated:                    1,280 CS + 41 SC = 1,321
partition passports:                  3,840 - 1,321 = 2,519
remaining partition statuses:                 SS 1,280; SC 1,239
five-black WW passports:             28,388 - 1,321 = 27,067
five-black passports:               129,788 - 1,321 = 128,467
global passports eliminated:        108,602 + 1,321 = 109,923
global passports remaining:      42,502,076 - 1,321 = 42,500,755
```

The `9,742` status-free necessary-system and `22,068` conservative binary-
status ledgers correctly remain unchanged. Every black record still supports
`SS`, and `1,239` also support unresolved `SC`.

The archive claims complete `CS` only after adjoining the prior 12 high-`a`
records, and complete `SC` only for the named tranche. It leaves
`complete_SC_no_go_claimed`, complete-partition, five-black-WW, and global
no-go flags false; five-black `WW` remains partial and degree `29` remains
open. An interrupted broader exploratory pass is not stored or counted.

## 6. Full replay, resources, and hygiene

One fresh in-process `runpy` replay executed the complete unchanged command
equivalent to

```text
--check --workers 4
```

It regenerated all `1,321` rows, ended after exactly `41` `SC` rows,
reproduced histogram

```text
CS: 101:1248, 109:32
SC: 101:40,   109:1
all: 101:1288, 109:33
```

and matched archive SHA-256
`53063a01bd7b83dfc1df35d2f9384c82eeda02d4704b812ec2415d16fd44a797`.
It exited successfully in `419.509` seconds. Darwin `resource.getrusage`
reported

```text
RUSAGE_SELF max RSS                    40,828,928 bytes = 38.938 MiB
RUSAGE_CHILDREN max RSS per maximum    13,123,584 bytes = 12.516 MiB
self + four child maxima               93,323,264 bytes = 89.000 MiB
```

At most four Singular children run concurrently. Multiplying the child
maximum by four and adding the parent maximum is a deliberately conservative
simultaneous bound, far below `0.6 GiB` (`614.4 MiB`).

The focused suite passed:

```text
6 passed in 0.20s
```

Independent checks found canonical duplicate-free JSON, valid Python syntax,
valid UTF-8, final LF, and no NUL, carriage return, or other forbidden C0
control byte in the four frozen artifacts. Their post-replay hashes remain
exactly the pinned values above.

The archive stores equation digests and standard-basis summaries rather than
full Bezout combinations. That is not a blocker: the unchanged full replay
reconstructed every equation, recomputed every successive basis, exercised
the ordered four-worker path, and matched the canonical archive bytes.

```text
hostile_audit_verdict = PASS
mathematical_blockers = 0
implementation_blockers = 0
scope_blockers = 0
resource_blockers = 0
partition_(3^2,1^23)_CS = IMPOSSIBLE
partition_(3^2,1^23)_SC_17_le_a_le_20 = IMPOSSIBLE
partition_(3^2,1^23)_SC_other = OPEN
partition_(3^2,1^23)_SS = OPEN
five_black_WW_status = PARTIAL
global_degree_29_status = OPEN
```
