# Hostile audit: five-black `WW` partition `(3^2,1^23)`, full `CC` tower

**Verdict (2026-08-10): PASS.** No mathematical, implementation, coverage,
resource, or claim-boundary blocker was found. The frozen theorem excludes the
`CC` refinement for all `1,292` canonical black-degree records. Accounting for
the `12` records already removed in every endpoint status by the high-marked-
degree theorem, this eliminates exactly `1,280` new marked passports. The
other endpoint statuses and the full degree-`29` problem remain open.

## Frozen artifacts

| artifact | SHA-256 |
|---|---|
| `scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.py` | `6b0ad80d811970156a38389aa1f1302e4f91419836bdeb3500633bf118eb53a7` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.json` | `6bc38604ce8c64519be1b71debdf48fddd9083616aee9409200334958e9a3630` |
| `results/gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.md` | `e31b7af6b356dff1057eb05461f74adcf467ee1a010383971eda4ff0265fe567` |
| `tests/test_gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.py` | `6c3df50d1c4bb00302f6c4aeb4a4fababac542f558c87ed9f63ef7b5a575b9de` |

Direct dependencies were pinned as follows:

| dependency | SHA-256 |
|---|---|
| `scripts/certify_gram_belyi_m5_degree29_reduction.py` | `5bd3eae89b0b45f4f8b279fd1298ba0c0598d2f154cde9b26dd1f5d778132837` |
| `scripts/certify_gram_belyi_m5_degree29_five_black_quotient.py` | `deaaed9f90e550ecd8a2256dc65a2bdc4ca57aea7e18354b7dabfff0d344936f` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_quotient.json` | `89e2d1fafc10ac5ad91edd6e20ed5e96f6c238b404b8bb303a054f330cbc0e87` |
| `scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition5.py` | `a8c3f46b126424a80f02f0d4aab11f15d2fe9425cbbeed8236a861f0a966ab10` |
| `scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition33_high_a.py` | `f0cd32af85947d626e2ad74e193dc56d6a414ad72937abe62342f1330dc2ec36` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_high_a.json` | `011c3a48880f43fb282ce0bd4416ec875c3f6f7dc4b8ffb96909dc1b46d60ba7` |

The canonical compact `.rows` serialization, including its final newline,
has SHA-256
`cbc1f8d51c689748c99915d855f43004a2099fc6230101e24840553d662586ab`.

## 1. Independent coverage

An enumeration independent of the producer formed all positive tuples

\[
 (a,b,c,d,e),\qquad a+b+c+d+e=29,\qquad b\le c\le d\le e.
\]

It gives `1,292` distinct records, exactly equal to the archived rows in
canonical order. Exactly `41` satisfy `17<=a<=20`. Exactly `12` satisfy
`a>=21`, and that list agrees row-for-row with the preceding high-`a`
archive.

The white partition `(3^2,1^23)` has exactly two critical white vertices.
In status `CC`, the distinct normalized endpoints are those two vertices;
after quotienting equal-valency identities there is one `CC` refinement per
black record. Thus the theorem covers all `1,292` `CC` marked passports,
with no missing orientation or extra multiplier.

## 2. Equations and exact base algebra

For

\[
G(x)=(x-\lambda)^a(x-u)^b(x-v)^c(x-w)^d(x-z)^e,
\]

removing the forced black-root powers from `G'` gives the residual quartic

\[
D(x)=\sum_i a_i\prod_{j\ne i}(x-r_j),\qquad \operatorname{lc}(D)=29.
\]

Both critical roots are fixed by `CC`, hence necessarily

\[
K=x(x-1),\qquad D=29K^2=29x^4-58x^3+29x^2,
\qquad G(0)=G(1).
\]

Consequently the four stored base equations are exactly
`D[0]=0`, `D[1]=0`, `D[2]-29=0`, and `D[3]+58=0`. Their exact total degrees
in `(u,v,w,z)` are `(4,3,2,1)`; the common-value equation has structural
homogenization degree `29`. Independent reconstruction reproduced every one
of the `1,292` equation SHA-256 values, with `1,292` distinct digests, and
found the degree vector `(4,3,2,1)` in every row.

For every row Singular reports

```text
BASE_DIM
0
24
UNIT
```

so the four-equation base quotient is zero-dimensional of exact vector-space
dimension `24=4*3*2*1`. Adjoining `G(0)-G(1)` to that exact finite algebra
gives the unit ideal. This tests the whole finite scheme, including
nonreduced components; it is not a numerical or point-sampling argument.
Allowing forbidden black-root collisions only enlarges the necessary system,
so the absence of saturations cannot invalidate the no-go.

All archived residues satisfy `sqrt5^2=5` and
`29*lambda^2-14*lambda+1=0` modulo their stored prime. The selected primes
are `101` and `109`, both split, larger than `29`, and locally good.

## 3. Projective gate and characteristic zero

Homogenize in `[u:v:w:z:h]`, replacing the marked root by `lambda*h` and
the endpoints by `0,h`. The base equations have structural degrees
`(4,3,2,1)` and the value equation degree `29`. At `h=0`, the base equations
give

\[
D_\infty(x)=29x^4.
\]

Group the five labelled black roots into collision clusters `alpha`, with
cluster size `n_alpha` and positive total black degree `A_alpha`. Directly
from the residual-derivative formula,

\[
D_\infty(x)=\prod_\alpha(x-\alpha)^{n_\alpha-1}E(x),\qquad
E(\alpha)=A_\alpha\prod_{\beta\ne\alpha}(\alpha-\beta).
\]

Because `1<=A_alpha<=29<p`, the displayed value is nonzero. The fixed
marked root is at `0` on the boundary, and the order-four zero of
`D_infinity` there therefore forces `n_0=5`. All five labels collide at
zero, so `u=v=w=z=0`; together with `h=0` this is the forbidden zero
projective vector. Hence there is no point at infinity. The degree-29 value
equation may have zero top form on this boundary, but that is harmless: the
four base equations already make the boundary empty.

At `p>29`, `29` is a unit and `lambda` is integral in the chosen localization
of `Q(sqrt(5))`; all homogenized coefficients are locally integral. The
closed projective model is proper. Any characteristic-zero point would,
after finite valuation extension, specialize either to an affine special-
fibre point (excluded by `[1]`) or to infinity (excluded above). Thus the
finite-field unit certificates rigorously imply characteristic-zero
nonexistence; nonintegral affine coordinates cannot escape the argument.

## 4. Ledger and claim boundary

The preceding high-`a` theorem removed all four statuses for `12` records,
so the new contribution of the full `CC` theorem is

\[
1292-12=1280.
\]

Starting from the prior frozen boundary, all updates recompute exactly:

```text
partition marked passports:       5,120 - 1,280 = 3,840
partition black records remaining:                 1,280
five-black WW passports:         29,668 - 1,280 = 28,388
five-black passports:           131,068 - 1,280 = 129,788
global passports eliminated:    107,322 + 1,280 = 108,602
global passports remaining:  42,503,356 - 1,280 = 42,502,076
```

The status-free necessary-system ledger `9,742` and conservative binary-
status ledger `22,068` correctly do not decrease, because `SS/SC/CS` remain.
The archive claims a complete `CC` no-go only. It leaves the complete
partition and five-black/global no-go flags false, five-black `WW` partial,
and degree `29` open.

## 5. Full replay, resources, and hygiene

A fresh in-process `runpy` invocation exercised the complete unchanged
`--check` path. It regenerated all `1,292` rows, reproduced the prime
histogram `101:1277, 109:15`, matched archive SHA-256
`6bc38604ce8c64519be1b71debdf48fddd9083616aee9409200334958e9a3630`,
and exited successfully in `22.636` seconds. Darwin `resource.getrusage`
reported:

```text
RUSAGE_SELF max RSS       = 38,846,464 bytes = 37.047 MiB
RUSAGE_CHILDREN max RSS   =  5,488,640 bytes =  5.234 MiB
conservative sum          = 44,335,104 bytes = 42.281 MiB
```

Only one Singular child runs at a time, so the sum is a conservative
simultaneous bound and is far below `0.6 GiB` (`614.4 MiB`).

An immediately preceding full replay also reached the certifier's successful
archive summary, but its wrapper then raised `NameError` while formatting
the resource line because the wrapper used lowercase `true`. This occurred
after `runpy` returned and did not affect the certifier or artifacts; the
corrected exit-zero replay above is the terminal correctness and resource
measurement.

The focused suite passed:

```text
6 passed in 0.20s
```

Independent checks found canonical JSON with no duplicate keys, valid Python
syntax, valid UTF-8, final LF, and no forbidden C0 control bytes in the four
frozen artifacts. Their post-replay hashes remain exactly those pinned
above. The environment was Python `3.14.6`, python-flint `0.9.0`, Singular
`4.4.1` (`44105`, arm64-Darwin), and FLINT `3.6.0`.

The archive stores rowwise equation hashes and standard-basis summaries, not
full Bezout combinations. That is not a blocker because the unchanged full
replay reconstructed every exact equation and recomputed every basis.

```text
hostile_audit_verdict = PASS
mathematical_blockers = 0
implementation_blockers = 0
resource_blockers = 0
partition_(3^2,1^23)_CC = IMPOSSIBLE
partition_(3^2,1^23)_status = PARTIAL
five_black_WW_status = PARTIAL
global_degree_29_status = OPEN
```
