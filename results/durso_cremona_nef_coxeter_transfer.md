# D'Urso's Cremona--nef domains: three Coxeter rows and a repaired `H13/15` lineage

**Status (2026-08-09): exact project replay of a published 2025
surface-to-polytope theorem; one project-new algebro-geometric interpretation,
but no new abstract polytope and no strict inverse `L4`.**

Luíze D'Urso studies the blow-up `S_n` of `P2` at `n` very general points and
the Cremona action on its nef cone.  Her
[`v2` paper](https://arxiv.org/abs/2412.15460), Theorem 6.2 and Proposition
6.3, gives a rational-polyhedral fundamental cone for the `K_S`-negative
nef-effective locus and proves that its hyperbolic projectivization is a
Coxeter polytope exactly when

```text
n = 10, 11, or 13.
```

The present certificate reconstructs the physical divisor vectors, checks
every facet, enumerates every extreme ray exactly, identifies the three
weighted diagrams in the pinned CoxIter corpus, and separates the new
`H13/15` interpretation from the different `H13/15` chamber already carried
by this project.

The resulting exact positive ledger is:

| project row | dimension/facets | finite | ideal | hyperideal | abstract lineage |
|---|---:|---:|---:|---:|---|
| `DURSO-CREM-P10` | `H10/12` | 17 | 2 | **0** | Vinberg's odd-unimodular `P10` |
| `DURSO-CREM-P11` | `H11/13` | 26 | 2 | **0** | Vinberg's odd-unimodular `P11` |
| `DURSO-CREM-TUM13-01` | `H13/15` | 44 | 2 | **0** | Tumarkin's `13-tum04_13_01` pyramid |

Thus all three are genuine finite-volume noncompact hyperbolic Coxeter
polytopes.  The first two already had complete rational-surface nef carriers
in this project.  The third is a different published abstract chamber and a
new project-level algebro-geometric interpretation: it is a Cremona
fundamental sector of the `K`-negative nef-effective cone of a very-general
thirteen-point plane blow-up.

## 1. The divisor walls and the exact angle criterion

Use the total-transform basis

```text
NS(S_n) = <H,E1,...,En> = I_(1,n),
H^2=1,  Ei^2=-1.
```

The `n+2` physical wall vectors are

```text
Ri = Ei-E(i+1),                  1 <= i < n,
Rn = En,
L  = H-E1-E2-E3,
A  = -K_S = 3H-E1-...-En.
```

They span `NS(S_n)` at index one.  All nonzero pairings are one: the `R_i`
form a chain, `L.R3=1`, and `Rn.A=1`.  The norms are

```text
Ri^2=L^2=-2 for i<n,
Rn^2=-1,
A^2=9-n.
```

All angles except the last are therefore `pi/3` or `pi/4`.  At the terminal
pair,

```text
cos(theta_n)^2 = (Rn.A)^2/(Rn^2 A^2) = 1/(n-9).
```

This gives a complete elementary iff:

| `n` | terminal angle | decision |
|---:|---|---|
| 10 | `0` | parallel Coxeter walls |
| 11 | `pi/4` | Coxeter |
| 12 | strictly between `pi/4` and `pi/3` | not `pi/m` |
| 13 | `pi/3` | Coxeter |
| `>=14` | strictly between `pi/3` and `pi/2` | not `pi/m` |

This proves the all-`n>=10` Coxeter criterion without a numerical inverse
cosine test.

## 2. Finite volume is checked independently

The source proves that the fundamental polytope has finite volume for every
`n>=10`; it is not inferred from finite-sidedness.  The project replay gives
an independent exact double-description check in the three Coxeter rows and
in the adjacent non-Coxeter controls `n=12,14`.

For each cone, the wall matrix has rank `n+1` and inertia `(1,n,1)`.  Every
inequality supports an `n`-dimensional facet.  The complete primitive ray
censuses are

| `n` | finite | ideal | hyperideal | Coxeter? |
|---:|---:|---:|---:|---|
| 10 | 17 | 2 | 0 | yes |
| 11 | 26 | 2 | 0 | yes |
| 12 | 35 | 2 | 0 | no |
| 13 | 44 | 2 | 0 | yes |
| 14 | 53 | 2 | 0 | no |

Because every extreme ray is finite or ideal, the closed Klein section lies
inside the closed hyperbolic ball and the Coxeter rows have finite volume.
They are noncompact because each has two primitive isotropic rays.

## 3. Exact abstract lineage

The normalized relation graphs were compared with the pinned CoxIter snapshot
at commit `b800be48d8240fbabd2d1c9542e61054584ad7fe`.

- `H10/12` is literally the `group units qf` record for Vinberg `P10`.
- `H11/13` is literally the corresponding Vinberg `P11` record.
- `H13/15` is literally CoxIter's `13-tum04_13_01`, from Tumarkin's
  classification of finite-volume `n+2`-facet polytopes
  ([primary paper](https://arxiv.org/abs/math/0301133)).

The last row is **not** the project's previously realized Vinberg `P13`,
which is `13-tum04_13_02`.  The weighted graphs are nonisomorphic, and their
ray splits already distinguish them:

```text
D'Urso / Tumarkin 01: 44 finite + 2 ideal + 0 hyperideal,
Vinberg P13 / Tumarkin 02: 43 finite + 3 ideal + 0 hyperideal.
```

So this is a real new ledger item and interpretation, not a duplicate hidden
by the shared `H13/15` header.  The abstract chamber itself has been known
since Tumarkin's 2004 publication.

## 4. Exact algebraic meaning—and its limit

D'Urso's cone is not asserted to be the entire nef cone of `S_n`, and the
Cremona action is not the automorphism group of a fixed very-general blow-up.
The exact statement is

```text
a rational-polyhedral fundamental sector for the Cremona action on
Nef^e(S_n) intersect {D.K_S <= 0}.
```

There is a second important distinction.  The reflections in every physical
divisor wall preserve `NS(S_n)` for `n=10,11`.  In `H13`, the primitive
canonical wall has square `-4` and divisibility one, so

```text
x |-> x + (x.A)/2 A
```

does **not** preserve `I_(1,13)`.  The `H13` chamber is an exact Coxeter
sector cut out by divisor classes, but it is not an integral
divisor-reflection group on the raw Neron--Severi lattice.  This prevents an
incorrect strengthening of the new interpretation.

Conversely, a marked copy of one of these three diagrams together with its
odd-unimodular Cremona presentation recovers the blow-up carrier and the same
fundamental sector.  That is a useful three-row `P -> surface -> same P`
recognition.  It is not the requested `L4`: the AG theorem supplies the cone,
but the elementary Coxeter-angle calculation itself selects the three rows.
No independent AG classification is being pushed back to a new theorem about
an a-priori class of polytopes.

## 5. Reproduction

The exact certificate is
[`certify_durso_cremona_nef_coxeter_transfer.py`](../scripts/certify_durso_cremona_nef_coxeter_transfer.py),
with frozen archive
[`durso_cremona_nef_coxeter_transfer.json`](../data/enumerations/durso_cremona_nef_coxeter_transfer.json)
and focused tests
[`test_durso_cremona_nef_coxeter_transfer.py`](../tests/test_durso_cremona_nef_coxeter_transfer.py).

The optional hostile source replay verifies both the current arXiv `v2` PDF
SHA-256 and all four CoxIter graph files before rebuilding the same
deterministic payload.
