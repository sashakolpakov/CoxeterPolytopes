# Hostile audit: five-black `WW` partition `(3^2,1^23)`, `a_lambda>=21`

**Verdict (2026-08-10): PASS.** No mathematical, implementation, coverage,
resource, or claim-boundary blocker was found. The frozen theorem excludes
exactly the `12` black-degree records with marked degree `a_lambda>=21`,
representing `48` marked passports. The other `1,280` records in this exact
partition remain open.

## Frozen artifacts

| artifact | SHA-256 |
|---|---|
| `scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition33_high_a.py` | `f0cd32af85947d626e2ad74e193dc56d6a414ad72937abe62342f1330dc2ec36` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_high_a.json` | `011c3a48880f43fb282ce0bd4416ec875c3f6f7dc4b8ffb96909dc1b46d60ba7` |
| `results/gram_belyi_m5_degree29_five_black_ww_partition33_high_a.md` | `c63ab0c7e9911c1dcad83367e0dae4f76c6f3d7833478b0de90eb9a41aaeeb1a` |
| `tests/test_gram_belyi_m5_degree29_five_black_ww_partition33_high_a.py` | `9a71568de465dc1d9ee38e4246e0adb7dd81e5d348acf4cd720fe364a8946524` |

The direct dependencies were also pinned:

| dependency | SHA-256 |
|---|---|
| `scripts/certify_gram_belyi_m5_degree29_reduction.py` | `5bd3eae89b0b45f4f8b279fd1298ba0c0598d2f154cde9b26dd1f5d778132837` |
| `scripts/certify_gram_belyi_m5_degree29_five_black_quotient.py` | `deaaed9f90e550ecd8a2256dc65a2bdc4ca57aea7e18354b7dabfff0d344936f` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_quotient.json` | `89e2d1fafc10ac5ad91edd6e20ed5e96f6c238b404b8bb303a054f330cbc0e87` |
| `scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition5.py` | `a8c3f46b126424a80f02f0d4aab11f15d2fe9425cbbeed8236a861f0a966ab10` |

The canonical `.rows` array, including the final newline from `jq -cS`, has
SHA-256
`b706e6b4f3025642081fe589feeab125658cc84f7a2c43e670938ad8d3e442c4`.

## 1. Independent 12-record coverage

A direct enumeration, independent of the producer's record generator, used

\[
 (a,b,c,d,e),\qquad a+b+c+d+e=29,qquad b\le c\le d\le e,
 \qquad a\ge21.
\]

It gives exactly

```text
a=21: (1,1,1,5), (1,1,2,4), (1,1,3,3), (1,2,2,3), (2,2,2,2)
a=22: (1,1,1,4), (1,1,2,3), (1,2,2,2)
a=23: (1,1,1,3), (1,1,2,2)
a=24: (1,1,1,2)
a=25: (1,1,1,1)
```

These `5+3+2+1+1=12` records are unique and occur in exactly this archive
order. No record with `a>=21` is missing and no record with `a<21` is
included.

For the white partition `(3^2,1^23)`, the simple degree occurs `23` times
and the critical degree `3` occurs twice. Quotienting equal-valency cycle
identities gives

```text
SS = 1, SC = 1, CS = 1, CC = 1.
```

The `CC` choice exists because there are two critical vertices. Hence each
black record represents four ordered endpoint markings and the quantified
sublocus contains exactly `12*4=48` marked passports.

## 2. Exact equations and residual arithmetic

For

\[
 G(x)=(x-\lambda)^a(x-u)^b(x-v)^c(x-w)^d(x-z)^e,
\]

remove the forced black-root powers from `G'`, leaving the residual quartic
`D` with leading coefficient `29`. The exact partition has two critical
white vertices `t,s`, each of ramification degree three, so necessarily

\[
                         D(x)=29(x-t)^2(x-s)^2.
\]

The target coefficients, from constant through leading, are

\[
 29t^2s^2, -58ts(t+s),
 29(t^2+4ts+s^2), -58(t+s), 29.
\]

The two critical points and both normalized endpoints must have the same
white value, giving

\[
 G(s)=G(t),\qquad G(0)=G(t),\qquad G(1)=G(t).
\]

Thus the producer sends seven exact equations in `(u,v,w,z,t,s)` to
Singular. Independent construction of all `12` archived rows reproduced the
exact degree vector

\[
                         (4,3,2,1,29,29,29)
\]

and leading coefficient `29` for `D` in every selected characteristic. The
four derivative-equation term counts were uniformly `(6,12,13,7)`. Across
the rows, the three value-equation term counts ranged respectively over

```text
G(s)-G(t): 830 ... 3,562
G(0)-G(t): 415 ... 1,781
G(1)-G(t): 430 ... 1,862
```

The system retains root collisions and the diagonal `t=s`. These are closed
degeneracy enlargements, so excluding the system remains a valid no-go for
the genuine stratum where the roots are distinct.

## 3. Affine unit certificates

Every selected prime is a split prime for `Q(sqrt(5))` and exceeds `29`.
Every stored residue satisfies

\[
 r^2=5\pmod p,
 \qquad
 29\lambda^2-14\lambda+1=0\pmod p.
\]

For each row, Singular computes `std(I)` in
`F_p[u,v,w,z,t,s]`. The producer accepts a prime only when the exact basis is
the singleton `[1]`; by the weak Nullstellensatz this excludes affine points
over the algebraic closure. A nonunit early reduction is only inconclusive
and is never promoted to a certificate.

All `12` frozen rows store the canonical record, valid residue embedding,
`affine_standard_basis=[1]`, `singular_stdout="UNIT"`, and a distinct valid
equation SHA-256. The first-successful-prime histogram is

```text
101: 11
131:  1
```

The sole fallback row is `(22,1,1,2,3)`. During replay the producer retained
the earlier nonunit reductions as inconclusive and continued to the exact
unit certificate at `p=131`.

## 4. Uniform projective-infinity exclusion

Homogenize the seven equations in

\[
                 [u:v:w:z:t:s:h]\in\mathbf P^6,
\]

where the marked root is `lambda*h`. At `h=0`, both endpoint equations give
`G_infinity(t)=0`, while `G(s)=G(t)` then gives
`G_infinity(s)=0`. Thus `t` and `s` are both black-root cluster values among
`0,u,v,w,z`.

For each distinct cluster `alpha`, let `n_alpha` be its number of labelled
black roots and `A_alpha` the sum of their positive degrees. Comparing the
ordinary derivative factorization with the labelled forced factors gives

\[
 D_\infty(x)=
 \prod_\alpha(x-\alpha)^{n_\alpha-1}E(x),
 \qquad
 E(\alpha)=A_\alpha
 \prod_{\beta\ne\alpha}(\alpha-\beta).
\]

Since `1<=A_alpha<=29<p`, `E(alpha)` is nonzero.

- If `t!=s`, the target quartic has order two at each root. Therefore
  `n_t=n_s=3`. These are disjoint clusters and would require six labelled
  black roots, but only five exist.
- If `t=s`, the target is `29(x-t)^4`, so `n_t=5`. All five labelled roots
  collide at `t`; the fixed marked root is `0` at infinity, forcing
  `u=v=w=z=t=s=0`, which is not a projective point.

These two cases exhaust infinity, including every collision and diagonal
degeneracy. The proof is uniform over all selected primes and all `1,292`
black records, not only the quantified high-`a` sublocus.

As a computational cross-check, an independently generated
`p=101`, `(25,1,1,1,1)` system returned `UNIT` for the affine ideal and all
six infinity charts `u=1`, `v=1`, `w=1`, `z=1`, `t=1`, `s=1`.

At a chosen prime of `Q(lambda)`, all coefficients are locally integral.
The homogenized equations define a closed projective scheme over the local
valuation ring. Properness would specialize any characteristic-zero point to
the special fibre, but the affine unit calculation and the infinity theorem
make that special fibre empty. Hence no solution can escape through
nonintegral coordinates.

## 5. Count arithmetic and claim boundary

The exact residual arithmetic is

\[
 1292-12=1280,
 \qquad
 5168-48=5120.
\]

This partition refinement does not remove a complete broad status-free
master, because the same black record persists in the other exact white
partitions. Therefore the status-free and binary-status ledgers correctly
remain `9,742` and `22,068`. The raw ledgers update by `48`:

\[
 29716-48=29668,
 \qquad
 131116-48=131068,
\]

\[
 107274+48=107322,
 \qquad
 42503404-48=42503356.
\]

The archive correctly leaves `complete_partition_no_go_claimed=false`, all
five-black and global complete-no-go flags false, the five-black `WW` status
partial, and the global degree-`29` status open.

## 6. Full replay and resource bound

The producer was run once through its complete `--check` path in a fresh
in-process `runpy` wrapper. It regenerated all `12` rows, tried fallback
primes in order, serialized the complete archive, and matched SHA-256
`011c3a48880f43fb282ce0bd4416ec875c3f6f7dc4b8ffb96909dc1b46d60ba7`.
The successful summary was

```text
status = EXACT_FIVE_BLACK_WW_PARTITION_33_HIGH_A_NO_GO
systems eliminated = 12
systems remaining = 1,280
marked passports eliminated = 48
marked passports remaining = 5,120
prime histogram = 101:11, 131:1
resource replay status = SUCCESS
```

The replay completed in `215.995` seconds. Darwin `resource.getrusage`, for
which `ru_maxrss` is measured in bytes, reported

```text
RUSAGE_SELF max RSS       =  36,044,800 bytes =  34.375 MiB
RUSAGE_CHILDREN max RSS   = 287,981,568 bytes = 274.641 MiB
conservative sum          = 324,026,368 bytes = 309.016 MiB
```

The producer launches only one Singular child at a time. The parent maximum
plus the largest-child maximum is therefore a conservative simultaneous RSS
bound, comfortably below `0.6 GiB` (`614.4 MiB`).

## 7. Tests, hashes, bytes, and environment

The focused suite completed with

```text
6 passed in 1.38s
```

Independent scans of the producer, archive, theorem report, and test found
valid UTF-8, final LF, valid Python or JSON syntax where applicable, and no
NUL, backspace, vertical-tab, form-feed, carriage-return, or other forbidden
C0 control bytes. Final hashes after replay and tests still equal the frozen
values above.

The audited environment was Python `3.14.6`, python-flint `0.9.0`, Singular
`4.4.1` (`44105`, arm64-Darwin), and FLINT `3.6.0`.

## 8. Auditability note

The archive stores equation hashes and unit-basis summaries rather than full
Singular Bezout combinations. Archive-only inspection therefore does not
re-prove each unit assertion. This is not a blocker: the unchanged full
replay reconstructed every equation and standard basis, exercised the prime
fallback logic, and matched the canonical archive bytes.

```text
hostile_audit_verdict = PASS
mathematical_blockers = 0
implementation_blockers = 0
resource_blockers = 0
partition_(3^2,1^23)_a_lambda_ge_21 = IMPOSSIBLE
partition_(3^2,1^23)_status = PARTIAL
five_black_WW_status = PARTIAL
global_degree_29_status = OPEN
```
