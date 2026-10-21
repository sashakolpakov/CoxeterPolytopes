# Silver C5 saturation: an exact rank-20 K3 carrier

## Result

The diagram-first algebro-geometric construction succeeds at the **carrier**
level.

Start with the order-five arithmetic Lanner cycle over

\[
F=\mathbf Q(\sqrt2),\qquad \mathcal O_F=\mathbf Z[\sqrt2],
\]

having one edge labelled `4` and four edges labelled `3`.  Its free trace
lattice has type

\[
E_8(-1)\oplus U(2).
\]

There is exactly one nonzero even index-two overlattice preserved by
multiplication by \(\sqrt2\).  It has type

\[
E_8(-1)\oplus U.
\]

Orthogonally adjoining either of two positive rank-five field blocks gives
two presentations of one rank-20 lattice:

\[
\begin{aligned}
C_5^{\mathrm{sat}}\oplus B_5
  &\rightsquigarrow U\oplus E_8(-1)\oplus D_{10}(-1),\\
C_5^{\mathrm{sat}}\oplus(F_4\oplus A_1)
  &\rightsquigarrow U\oplus E_8(-1)^2\oplus A_1(-1)^2.
\end{aligned}
\]

Both are even hyperbolic 2-elementary lattices with invariants

\[
(r,a,\delta)=(20,2,1).
\]

The second presentation is exactly the Picard lattice in Alexeev's
rank-20 example.  It embeds primitively in the K3 lattice with orthogonal
complement

\[
T=\langle2\rangle\oplus\langle2\rangle.
\]

This is a new exact bridge from the prescribed Coxeter missing face to a
specific algebraic surface lattice.  It is **not** yet a new Coxeter
polytope: saturation creates an additional \((-2)\)-root, so the complete
wall arrangement and Vinberg's finite-volume/compactness tests remain to be
run.

## 1. Trace convention

Use the totally negative codifferent scalar

\[
\alpha=-\frac{2+\sqrt2}{4}.
\]

For \(h=a+c\sqrt2\), put

\[
Q=\operatorname{Tr}(\alpha h)=-(a+c),\qquad
B=\operatorname{Tr}(\alpha\sqrt2h)=-(a+2c).
\]

Thus the field value is recovered from two divisor pairings by

\[
h=(B-2Q)+(Q-B)\sqrt2.
\]

In the interleaved basis

\[
e_0,\sqrt2e_0,e_1,\sqrt2e_1,\ldots,e_4,\sqrt2e_4,
\]

multiplication by \(\sqrt2\) is the integral block action

\[
R_0=\begin{pmatrix}0&2\\1&0\end{pmatrix}^{\oplus5},
\qquad R_0^2=2I.
\]

It is self-adjoint for the trace pairing.

## 2. The unique stable even saturation

Let \(G_0\) be the free C5 trace Gram matrix.  Exact Smith reduction and
inertia give

\[
\det G_0=-4,\qquad
\operatorname{inertia}(G_0)=(1,9,0),\qquad
A_{G_0}\simeq(\mathbf Z/2)^2.
\]

Every denominator-two discriminant class can be represented by \(b/2\)
with \(b\in\{0,1\}^{10}\).  Exhausting all \(2^{10}\) binary vectors leaves
the following four dual classes.  Supports are zero-indexed trace
coordinates.

| support of `b` | `(b/2)^2` | even overlattice | `sqrt(2)`-stable |
|---|---:|---:|---:|
| `[]` | `0` | trivial | yes |
| `[3,5,9]` | `-2` | yes | **yes** |
| `[1,2,4,7,8]` | `-1` | no | no |
| `[1,2,3,4,5,7,8,9]` | `0` | yes | no |

Hence the unique nontrivial stable even class is

\[
x=\frac{\sqrt2e_1+\sqrt2e_2+\sqrt2e_4}{2}.
\]

Replacing the last free trace basis vector by \(x\) gives an integral even
Gram matrix \(G_{\mathrm{sat}}\) with

\[
\det G_{\mathrm{sat}}=-1,
\qquad
\operatorname{inertia}(G_{\mathrm{sat}})=(1,9,0).
\]

The certificate contains a determinant-one integral change of basis whose
congruence is literally \(E_8(-1)\oplus U\); the conclusion therefore does
not rest only on classification by rank and determinant.  Conjugating
\(R_0\) through the index-two basis change remains integral and verifies

\[
R^2=2I,\qquad R^{t}G_{\mathrm{sat}}=G_{\mathrm{sat}}R.
\]

## 3. The extra root is the wall-closure warning

The saturation generator itself satisfies

\[
x^2=-2,
\qquad
x\cdot Rx=-1,
\qquad
Rx=e_1+e_2+e_4,
\qquad
(Rx)^2=-4.
\]

The trace-recovery formula gives

\[
h(x,x)=3-\sqrt2.
\]

This value is positive at both real embeddings, but it is not the original
field-root norm `2`.  Thus saturation does exactly what a divisor-lattice
construction must be expected to do: it adds integral divisor roots that
were invisible in the initial norm-two field diagram.  On a K3 surface one
of \(x\) and \(-x\) is effective by Riemann--Roch.  It cannot be discarded
from the nef/Weyl-wall analysis merely because it was absent from the seed
Coxeter labels.

This is the decisive next gate: enumerate the full \((-2)\)-root closure in
the chosen \(\sqrt2\)-eigenslice, orient the effective roots coherently, and
only then apply Vinberg's exact edge-endpoint and full parabolic-neighborhood
tests.

## 4. Two finite completions

The following identifications are certified by explicit unimodular
congruences.

| field block | field determinant | trace lattice | trace determinant |
|---|---:|---|---:|
| `B5` with labels `(4,3,3,3)` | `2` | `D10(-1)` | `4` |
| `F4 + A1` | `2` | `E8(-1) + A1(-1)^2` | `4` |

Each field block is positive definite at both real embeddings.  Adding it
orthogonally to the saturated C5 block gives a rank-ten field form with

\[
\det_F=-4-4\sqrt2,
\]

identity-embedding inertia \((9,1,0)\), and conjugate-embedding inertia
\((10,0,0)\).  Its saturated trace lattice has

\[
\det=-4,\qquad
\operatorname{inertia}=(1,19,0),qquad
A_N\simeq(\mathbf Z/2)^2,qquad
\delta=1.
\]

Multiplication by \(\sqrt2\) is integral and self-adjoint on both completed
lattices and still squares to `2`.

The two displayed trace presentations are isometric.  This also follows
independently from Nikulin's uniqueness theorem for indefinite even
2-elementary lattices with fixed signature and \((r,a,\delta)\).

## 5. Explicit primitive K3 embedding

Use the Alexeev-aligned presentation

\[
N=U\oplus E_8(-1)^2\oplus A_1(-1)^2
\]

and let \(a_1,a_2\) be its last two roots.  Put

\[
T=\mathbf Zt_1\oplus\mathbf Zt_2,
\qquad t_i^2=2,
\qquad t_1t_2=0.
\]

The two discriminant generators have opposite quadratic values:

\[
q_N(a_i/2)=-\frac12,
\qquad
q_T(t_i/2)=\frac12.
\]

Adjoin the two glue vectors

\[
g_i=\frac{a_i+t_i}{2}.
\]

The resulting index-four overlattice of \(N\oplus T\) is integral and even.
Its exact invariants are

\[
\det=-1,qquad
\operatorname{inertia}=(3,19,0).
\]

It is therefore the unique even unimodular lattice
\(II_{3,19}=U^3\oplus E_8(-1)^2\).  In the stored gluing basis the first
twenty vectors are literally a basis of \(N\), so the embedding is
primitive; the certified orthogonal-complement basis is

\[
t_i=2g_i-a_i,
\]

with Gram matrix `diag(2,2)`.

Consequently \(N\) is a genuine Picard-rank-20 K3 lattice.  Taking an
oriented orthogonal basis of \(T\) supplies the singular-K3 period.  The
self-adjoint \(\sqrt2\)-action is certified on \(N\); no claim is made here
that it is induced by a surface automorphism or an algebraic self-map.

## 6. Identification with Alexeev's surface

Alexeev states that an indefinite even 2-elementary lattice is determined by
its signature and \((r,a,\delta)\), and presents the rank-20 case as

\[
(20,2,1)\simeq U\oplus E_8(-1)^2\oplus A_1(-1)^2.
\]

He calls the associated complex K3 surface unique and computes its full
root-reflection chamber and elliptic pencils.  See
[Alexeev, arXiv:2209.09110v4](https://arxiv.org/abs/2209.09110), especially
the introduction and the subsection on `(20,2,1)`.  The v4 ancillary Coxeter
matrix used elsewhere in this repository has SHA-256

```text
4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457
```

Thus this inverse design did not merely find some abstract K3-admissible
lattice: it landed on the distinguished rank-20 reflective K3 lattice already
carrying Alexeev's 50-wall full root chamber.

What is new in the present calculation is the route into that lattice:
the one-`4` C5 missing-face block, its forced stable saturation, and an
integral self-adjoint \(\sqrt2\)-action are reconstructed directly from the
desired Coxeter arithmetic.

## 7. Exact scope

Proved:

- the complete discriminant-class enumeration of the C5 trace block;
- uniqueness of its nonzero even \(\sqrt2\)-stable index-two saturation;
- explicit integral isometries to `E8(-1)+U`, `D10(-1)`, and
  `E8(-1)+A1(-1)^2`;
- two exact rank-20 carrier presentations with invariants `(20,2,1)`;
- integrality, self-adjointness, and square `2` of the \(\sqrt2\)-action;
- an explicit primitive K3 embedding with complement `<2>+<2>`; and
- identification of the Picard lattice with Alexeev's rank-20 lattice.

Not proved:

- effectivity and common nef orientation of the proposed seed roots;
- closure under all \((-2)\)-roots introduced by saturation;
- acuteness or completeness of an eigenslice wall system;
- finite volume under Vinberg's full-neighborhood criterion; or
- compactness or a new hyperbolic Coxeter polytope.

The carrier is therefore a rigorous algebro-geometric search object, not a
finished example.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_silver_k3_carrier.py
PYTHONPATH=src:. pytest -q \
  tests/test_silver_k3_carrier.py \
  tests/test_silver_k3_carrier_certificate.py
```

Frozen certificate:

```text
data/candidates/silver_k3_carrier.json
SHA-256 f26310dc2c9475653f8d64fdabe7ba15abe312953621f5f757fc7fb344dec0b9
```
