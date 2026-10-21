# Cubic-seven `Lambda22` prefix through H=100

## Outcome

Let

\[
K=\mathbb Q(\theta),\qquad \theta^3+\theta^2-2\theta-1=0,
\]

with the largest real root as the distinguished embedding, and consider

\[
L=\langle-\theta\rangle\oplus\Lambda_{22}.
\]

The exact simultaneous-height Vinberg prefix contains **no reflecting wall
at height zero and no positive-height wall through `H=100`**.  The control
enumerator exhausts 30 locally admissible controls at 24 exact heights.  An
exact rank-three Gram scan followed by the required discriminant-module gate
rejects every one.

This is only a bounded zero-wall prefix.  With no walls there is no closed
polyhedral chamber on which an extreme-ray finite-volume or compactness test
can be run.  In particular, this result makes no reflectivity or global
nonreflectivity claim.

## Authoritative integral spatial lattice

The vendored 22-by-22 integral Gram matrix is the exact lower-triangular
matrix on the [Nebe--Sloane `LAMBDA22` catalogue
page](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/LAMBDA22.html).
No floating Euclidean basis is used.

| item | exact value |
|---|---|
| source HTML SHA-256 | `b4e4f08d2c4a4926599b6551ab852b578877f7c3cbebbfed57ea8b8924e3aab4` |
| vendored file SHA-256 | `4f3dbdc96e7aead7fae09a834abdbe2d5319620ece24656720f13184dd85b2f1` |
| normalized full Gram stream SHA-256 | `53e1adff6866521330edfeb17d79e9813398c9f3dc73503713761e4153f46495` |
| determinant | 12 |
| Smith diagonal | `1^20,2,6` |
| minimum | 4 |
| kissing number | 49,896 |

Independent PARI `qfminim` on the exact 22-dimensional Gram returns minimum
4, signed count 49,896, and 24,948 representatives modulo sign.  The raw
representative-matrix stream has SHA-256
`9b752381c953619264cf7eb1aba9b73480e6de94553a68db65eec8bfda30f428`.
This is a complete shortest-vector enumeration; no shell-orbit transitivity
is assumed.

## Complete local mirror norms

For a primitive vector `r=(x,y)`, put

\[
q(r)=-\theta x^2+y^2,
\qquad I(r)=(\theta x,(y,\Lambda_{22})).
\]

Reflection integrality gives `(q(r)) | 2I(r)`.  Away from the inert primes
`(2)` and `(3)`, local unimodularity and primitivity force zero valuation.
The Smith diagonal has elementary-divisor exponent one at both discriminant
primes.

At `(3)`, the factor 2 is a unit, so `v_3(q)<=1`.  At `(2)`, an odd `x`
makes `q` odd because `Lambda22` is even.  If `x` is even, primitive `y` has
`v_2((y,Lambda22))<=1`, while the factor 2 in `2I` permits
`v_2(q)<=2`.  Thus the exponent-two classes must not be discarded:

\[
v_2(q)\in\{0,1,2\},\qquad v_3(q)\in\{0,1\}.
\]

The field has class number one and `U^+=U^2`, so the complete canonical
totally-positive norm alphabet is

\[
\boxed{1,2,4,3,6,12}.
\]

For `m=1,2,3,6`, define the exact pairing-radical preimage

\[
N_m=\{y\in\Lambda_{22}:(y,\Lambda_{22})\subset m\mathbb Z\}.
\]

Smith decomposition followed by canonical column HNF gives the following
indices.  Independent PARI enumeration gives both the ordinary minimum and
the minimum of the 66-dimensional `O_K` trace lattice.

| `m` | index `[Lambda22:N_m]` | ordinary min | signed shortest count | trace min |
|---:|---:|---:|---:|---:|
| 1 | 1 | 4 | 49,896 | 12 |
| 2 | `2^20` | 12 | 8,448 | 36 |
| 3 | `3^21` | 24 | 1,782 | 72 |
| 6 | `2^20 3^21` | 96 | 1,782 | 288 |

Consequently the spatial coefficient gates are:

| mirror norm | coefficient lattice | coefficient minimum |
|---:|---|---:|
| 1, 2 | `N_1=Lambda22` | 4 |
| 4 | `N_2` | 12 |
| 3, 6 | `N_3` | 24 |
| 12 | `N_6` | 96 |

The six height-zero target traces are respectively `3,6,12,9,18,36`, all
strictly below the corresponding trace minima `12,12,36,72,72,288`.
Therefore this host is rootless for every complete mirror-norm class.

## Exact H=100 control census

All archimedean control boxes, time-coordinate reflection divisibility, and
the two positive-definite conjugate targets are tested exactly.

| norm | 1 | 2 | 4 | 3 | 6 | 12 | total |
|---|---:|---:|---:|---:|---:|---:|---:|
| controls | 5 | 16 | 5 | 1 | 2 | 1 | 30 |

These controls occupy 24 distinct exact heights.  The final one in physical
height order has norm 6, control `6+9theta+3theta^2`, and height
`15+87/2 theta+39/2 theta^2`.  Every control with height at most 100 is
included.

## Complete necessary Gram scan

Write a prospective spatial vector as

\[
y=u+\theta v+\theta^2w
\]

and abbreviate

\[
(U,V,W,P,Q,R)=(u^2,v^2,w^2,u\!\cdot\!v,u\!\cdot\!w,v\!\cdot\!w).
\]

For target `y^2=A+B theta+C theta^2`, coefficient comparison gives

\[
A=U+2R-W,\quad B=2P+4R-W,\quad C=2Q+V-2R+3W.
\]

The inverse trace Gram has diagonal `1,2/7,2/7`.  Hence, with
`T=Tr(y^2)`, the finite bounds `U<=T` and `V,W<=2T/7` are complete.  The
enumerator retains every integral positive-semidefinite Gram, including
rank-one and rank-two cases, and checks the exact minimum of its integral
kernel quotient.

Across all 30 controls there are 89 ambient minimum-four Gram types, with
rank histogram

```text
rank 1 / rank 2 / rank 3 = 3 / 6 / 80.
```

Only the five norm-4 controls and the single norm-12 control have any ambient
type.  After imposing the pairing congruence, 11 types remain, with histogram
`3 / 6 / 2`.  The semidefinite cases are deliberately retained.

The first ambient necessary type occurs for norm 4 at height
`2theta+theta^2`:

```text
(U,V,W,P,Q,R;det) = (4,4,0,4,0,0;0), rank 1.
```

It fails `min(N_2)=12`.  At height `10+29theta+13theta^2`, five
congruence-compatible norm-4 types remain, including three degenerate types
and two positive-definite types.  Every one has a nonzero `U`, `V`, or `W`
below 12, so all fail the same exact gate.

For norm 12, 78 ambient types reduce to the single congruence-compatible
semidefinite type

```text
(U,V,W,P,Q,R;det) = (12,144,36,36,18,72;0), rank 2.
```

Its nonzero coefficient squares 12 and 36 are below `min(N_6)=96`.  Thus all
11 congruence-compatible types have an explicit diagonal minimum
contradiction, and the fully gated type count is zero.

Therefore every one of the 30 exact controls is unrealized, and the bounded
prefix has zero positive walls through `H=100`.

## Reproducibility and scope

- Vendored Gram: `data/lattices/lambda22_nebe_sloane_gram.txt`
- Backend: `src/coxeter_polytopes/search/cubic7_lambda22.py`
- Generator: `scripts/audit_cubic7_lambda22_h100.py`
- Tests: `tests/test_cubic7_lambda22_h100.py`
- Archive: `data/enumerations/cubic7_lambda22_h100.json`

The archive stores all 30 exact controls, all 89 ambient types, all 11
pairing-congruence types (with ranks and determinants), every radical gate,
the exact frontier, and raw PARI stream hashes.  Its scope is strictly
`H<=100` at the distinguished embedding.
