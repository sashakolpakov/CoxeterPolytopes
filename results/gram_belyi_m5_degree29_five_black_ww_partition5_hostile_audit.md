# Hostile audit: degree-29 five-black `WW` partition `(5,1^24)`

**Verdict (2026-08-10): PASS.** No mathematical, implementation, coverage,
resource, or claim-boundary blocker was found. The frozen certificate excludes
all `1,292` black-degree refinements of the exact white partition `(5,1^24)`,
representing `3,876` marked passports. It does not claim a complete
five-black `WW`, five-black, or degree-`29` no-go.

## Frozen artifacts

| artifact | SHA-256 |
|---|---|
| `scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition5.py` | `a8c3f46b126424a80f02f0d4aab11f15d2fe9425cbbeed8236a861f0a966ab10` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition5.json` | `875429b35c8601375b6d2063686d606b69bb819fb120f02f1094c3e59fcdbb91` |
| `results/gram_belyi_m5_degree29_five_black_ww_partition5.md` | `bc1ffdae976413d99b8d958bf616325d41f892b989a5659b025ea5c51ddec95d` |
| `tests/test_gram_belyi_m5_degree29_five_black_ww_partition5.py` | `ab5fc49f09bfb30247dbb4539fff36a52725e3199917c3441e78002f89015fd0` |

The producer's direct dependencies were also pinned:

| dependency | SHA-256 |
|---|---|
| `scripts/certify_gram_belyi_m5_degree29_reduction.py` | `5bd3eae89b0b45f4f8b279fd1298ba0c0598d2f154cde9b26dd1f5d778132837` |
| `scripts/certify_gram_belyi_m5_degree29_five_black_quotient.py` | `deaaed9f90e550ecd8a2256dc65a2bdc4ca57aea7e18354b7dabfff0d344936f` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_quotient.json` | `89e2d1fafc10ac5ad91edd6e20ed5e96f6c238b404b8bb303a054f330cbc0e87` |
| `scripts/certify_gram_belyi_m5_degree29_four_black_ww_simple.py` | `6528cf01fcd16f90ef65389226dfd1e094794d4f48490e97b6eaa86b8beeb8b1` |

The canonical `.rows` array, serialized by `jq -cS` with its final newline,
has SHA-256
`b25f447bbce3f63c674f24f4639c9db8350a68a63d6f28ec8a938634eadbec36`.

## 1. Independent coverage and status-free quotient

A direct enumeration, without calling the producer's record generator, used

\[
 (a,b,c,d,e),\qquad a+b+c+d+e=29,\qquad b\le c\le d\le e,
\]

where `a` is the degree at the marked black root. It produced `1,292` unique
records in exactly the archive order. Its compact record-list digest is
`fc18d3c8fc308091820ef68d9ff6642cff96d4660eed39ba684ba6afffeccea0`,
agreeing with the independently frozen five-black quotient.

The quotient-polynomial argument applies to every selected white endpoint,
not only to `WW`. If

\[
                 G=DQ+R,\qquad \deg Q=25,
\]

and the residual critical points have common value `R`, then a simple white
point `alpha` has `D(alpha)!=0` and hence `Q(alpha)=0`. At a critical white
point of ramification index `r`,

\[
 \operatorname{ord}_\alpha(G-R)=r,
 \qquad \operatorname{ord}_\alpha D=r-1,
\]

so `Q` again has a root at `alpha`. Thus simple/critical endpoint statuses
give the same necessary quotient equation. The exact remaining-system
distinction is consequently

| endpoint block | conservative binary-status systems | status-free necessary systems |
|---|---:|---:|
| `BW` | 8,450 | 4,225 |
| `WB` | 8,450 | 4,225 |
| `WW` | 5,168 | 1,292 |
| **total** | **22,068** | **9,742** |

The archive and report correctly retain both numbers with different labels.
The partition-specific theorem removes no entire broad status-free system,
since the same black-degree system still represents the other exact white
partitions.

## 2. Exact partition equations

For one black-degree record write

\[
 G(x)=(x-\lambda)^a(x-u)^b(x-v)^c(x-w)^d(x-z)^e.
\]

Removing the forced black-root powers gives

\[
 G'=\prod_i(x-r_i)^{a_i-1}D(x),
 \qquad
 D(x)=\sum_i a_i\prod_{j\ne i}(x-r_j),
\]

where `D` has degree four and leading coefficient `29`. For the white
partition `(5,1^24)`, the only residual critical point has derivative
multiplicity four, so necessarily

\[
                         D(x)=29(x-t)^4.
\]

Equating its four non-leading coefficients and imposing

\[
                         G(0)=G(t)=G(1)
\]

gives six closed equations in `(u,v,w,z,t)`. These equations are necessary
for every genuine configuration. Omitting root-distinctness, fibre
inequations, and `t!=0,1` only enlarges the system and is safe for a no-go.
Allowing `t=0` or `t=1` covers `CS` or `SC`; `SS` has `t` at neither
endpoint. `CC` is impossible because this partition has only one critical
white vertex.

An independent expansion reproduced the exact total degrees

```text
derivative coefficient equations: 4, 3, 2, 1
endpoint equations:              29, 29
```

For `p=101` and weights `(1,1,1,1,25)`, the six term counts were
`(6,11,11,6,415,623)`, independently matching the intended equations.

## 3. Affine finite-field criterion

Every selected prime is prime, exceeds `29`, and splits `Q(sqrt(5))`.
The stored values satisfy both

\[
 s^2=5\pmod p,
 \qquad
 29\lambda^2-14\lambda+1=0\pmod p.
\]

For each row, Singular computes `std(I)` over
`F_p[u,v,w,z,t]`. The exact test

```text
size(J) == 1 and J[1] == 1
```

is precisely the unit-ideal criterion. By the weak Nullstellensatz it rules
out affine points over the algebraic closure, not only `F_p`-rational points.
A nonunit result at an early prime is treated as inconclusive and the next
split prime is tried; it is never used as negative evidence.

The archive contains `1,292` rows, all with the exact canonical record,
valid residue embedding, `affine_standard_basis=[1]`,
`singular_stdout="UNIT"`, and the uniform infinity flag. All `1,292`
equation digests are distinct valid SHA-256 strings. The first-successful
prime histogram is

```text
101: 1,263
109:    28
131:     1
```

## 4. Complete projective-infinity argument

Homogenize in

\[
                   [u:v:w:z:t:s]\in\mathbf P^5,
\]

with the marked root represented by `lambda*s`. The coefficient equations
have fixed degrees `(4,3,2,1)` and the endpoint equations degree `29`. At
`s=0`, both endpoint equations have the same leading form
`-G_infinity(t)`. Therefore `t` is a black-root cluster among
`0,u,v,w,z`.

Group collided labelled roots by distinct values `alpha`. Let `n_alpha` be
the number of root labels in that cluster and `A_alpha` the sum of their
positive black degrees. Combining the two expressions for the derivative
gives

\[
 D_\infty(x)=
 \prod_\alpha (x-\alpha)^{n_\alpha-1}E(x),
\]

where

\[
 E(\alpha)=A_\alpha
 \prod_{\beta\ne\alpha}(\alpha-\beta).
\]

For every selected characteristic, `1<=A_alpha<=29<p`, so this value is
nonzero. Hence the order of `D_infinity` at the black root `t` is exactly
`n_t-1`. The equation `D_infinity=29(x-t)^4` forces `n_t=5`. All five
labelled roots therefore collide at `t`; because the fixed marked root is
`0` at infinity, this gives

```text
u = v = w = z = t = 0,
```

which is not a projective point. This exhausts infinity uniformly for every
record; no genericity or root-distinctness assumption is being made on the
special fibre.

As a computational cross-check, an independently generated representative
at `p=101`, weights `(1,1,1,1,25)`, returned `UNIT` for the affine ideal and
for each of the five infinity charts `u=1`, `v=1`, `w=1`, `z=1`, `t=1`.
Its two endpoint leading forms were identical, as required.

At a selected prime of `Q(lambda)`, `lambda` and every equation coefficient
are locally integral. The homogenized equations define a closed projective
scheme over the local valuation ring. If its characteristic-zero generic
fibre had a point, properness would specialize a point to the geometric
special fibre. The affine unit calculation and the infinity theorem show
that this special fibre is empty. Thus no characteristic-zero solution can
escape through nonintegral coordinates.

## 5. Independent count arithmetic and claim boundary

The five exact white partitions have ordered `SS/SC/CS/CC` multipliers

| white partition | `SS` | `SC` | `CS` | `CC` | total |
|---|---:|---:|---:|---:|---:|
| `(5,1^24)` | 1 | 1 | 1 | 0 | 3 |
| `(4,2,1^23)` | 1 | 2 | 2 | 2 | 7 |
| `(3^2,1^23)` | 1 | 1 | 1 | 1 | 4 |
| `(3,2^2,1^22)` | 1 | 2 | 2 | 3 | 8 |
| `(2^4,1^21)` | 1 | 1 | 1 | 1 | 4 |
| **total** | **5** | **7** | **7** | **7** | **26** |

Therefore

\[
 1292\cdot3=3876,
 \qquad
 33592-3876=1292\cdot23=29716.
\]

The remaining exact-partition rows independently reproduce

```text
(4,2,1^23):    9,044
(3^2,1^23):    5,168
(3,2^2,1^22): 10,336
(2^4,1^21):    5,168
total:         29,716
```

The wider ledgers also agree:

\[
 134992-3876=131116,
 \qquad
 103398+3876=107274,
 \qquad
 42507280-3876=42503404.
\]

It is correct to mark only the exact `(5,1^24)` `WW` partition impossible.
The archive keeps `five_black_ww_status="PARTIAL"`, all complete-no-go flags
false, the global degree-`29` status open, and `degree_29_witness_found=false`.

## 6. Full replay and resource bound

The unchanged producer was first run through its `--check` path with progress
every `100` records. It regenerated every row, matched the frozen archive
byte-for-byte, and printed

```text
status = EXACT_FIVE_BLACK_WW_PARTITION_5_NO_GO
systems = 1,292
marked passports = 3,876
prime histogram = 101:1263, 109:28, 131:1
archive SHA-256 = 875429b35c8601375b6d2063686d606b69bb819fb120f02f1094c3e59fcdbb91
SCRIPT_MAIN_STATUS = 0
```

The terminal replay took `189.96` seconds wall time (`178.17` user,
`8.59` system). After the successful child exited, macOS `/usr/bin/time -l`
was denied its auxiliary `sysctl kern.clockrate` query by the sandbox and
therefore omitted its RSS line. This post-run timer failure is not a producer
failure.

A second full `--check` replay was consequently run in a fresh in-process
`runpy` wrapper solely to measure `resource.getrusage`. It again matched the
same archive hash and completed in `195.998` seconds. On Darwin, where
`ru_maxrss` is measured in bytes, it reported

```text
RUSAGE_SELF max RSS       = 39,632,896 bytes = 37.797 MiB
RUSAGE_CHILDREN max RSS   = 16,334,848 bytes = 15.578 MiB
conservative sum          = 55,967,744 bytes = 53.375 MiB
```

The producer launches only one Singular child at a time, so the sum of the
parent maximum and largest-child maximum is a conservative simultaneous RSS
bound. It is far below `0.6 GiB = 644,245,094.4 bytes`.

## 7. Focused tests, bytes, and environment

The focused suite completed with

```text
7 passed in 0.96s
```

Independent scans of the producer, archive, theorem report, and focused test
confirmed valid UTF-8, final LF, valid Python or JSON syntax where applicable,
and no NUL, backspace, vertical-tab, form-feed, carriage-return, or other
forbidden C0 control bytes. Final hashes after both replays and the tests
still equal the frozen values above.

The audited environment was Python `3.14.6`, python-flint `0.9.0`, Singular
`4.4.1` (`44105`, arm64-Darwin), and FLINT `3.6.0`.

## 8. Auditability note

The archive stores equation digests and unit-basis summaries rather than
full Singular Bezout combinations. Archive-only inspection therefore does
not re-prove every unit assertion. This is not a blocker here: two unchanged
full replays reconstructed every equation, tried primes in the prescribed
order, recomputed every exact standard basis, and matched the canonical
archive bytes.

```text
hostile_audit_verdict = PASS
mathematical_blockers = 0
implementation_blockers = 0
resource_blockers = 0
five_black_WW_partition_(5,1^24) = IMPOSSIBLE
five_black_WW_status = PARTIAL
five_black_status = PARTIAL
global_degree_29_status = OPEN
```
