# Supersingular K3 surfaces do not extend the golden eigen-nef factory to H10

## Result

Let `X` be a supersingular K3 surface over an algebraically closed field of
arbitrary positive characteristic, and put

\[
L=\operatorname{NS}(X).
\]

There is no self-adjoint integral endomorphism

\[
\tau\colon L\longrightarrow L,
\qquad \tau^2-\tau-1=0.
\]

Consequently, although a supersingular K3 surface has Picard rank 22, it cannot
carry the rank-11 golden module needed for a quadratic totally-real eigenslice
of signature `(1,10)`.  Together with the Picard-rank bound in characteristic
zero and for finite-height K3 surfaces, the integral golden full-NS eigen-nef
construction therefore has the same sharp dimensional ceiling in every
characteristic:

\[
\boxed{d\leq 9.}
\]

This is an obstruction to this particular algebro-geometric factory, not a
non-existence theorem for compact Coxeter polytopes in H10.  It does not exclude
a rational correspondence that fails to preserve the full integral NS lattice,
an action only on a finite-index sublattice, or a modified slice with additional
NS walls.

## Input from supersingular K3 theory

For a supersingular K3 surface of Artin invariant `sigma`, its Neron--Severi
lattice is an even hyperbolic lattice of rank 22 and discriminant

\[
-p^{2\sigma}.
\]

For odd `p` it is `p`-elementary.  For `p=2` it is 2-elementary of type I: its
discriminant quadratic form takes only integral values.  Conversely these
properties determine the lattice.  This is the Rudakov--Shafarevich
classification in the form recalled by Kondo--Shimada, Theorem 3.1 in
[arXiv:1212.0269](https://arxiv.org/abs/1212.0269).

Only these lattice properties are used below.

## The local golden trace form

Write

\[
\mathcal O=\mathbb Z[\phi],\qquad \phi^2-\phi-1=0.
\]

If `tau` existed, self-adjointness would make `L` a torsion-free
`O`-module.  Locally at 2 it is a finite torsion-free module over the discrete
valuation ring

\[
\mathcal O_2=\mathcal O\otimes\mathbb Z_2.
\]

It is therefore free.  Its `Z_2`-rank is 22, so its `O_2`-rank is 11.

The prime 2 is inert in `Q(sqrt(5))`, so

\[
\mathcal O_2/2\mathcal O_2\cong\mathbb F_4.
\]

Put `K_2=Q_2(phi)`, the unramified quadratic extension of `Q_2` with
ring of integers `O_2`.

Self-adjointness says that \(B(ax,y)=B(x,ay)\) for
\(a\in\mathcal O_2\).  Let \(\delta=\sqrt5\), a generator of the different.
Trace duality then gives a unique symmetric `O_2`-bilinear form `h` on the
free rank-11 module `M` such
that

\[
B(x,y)=\operatorname{Tr}_{K_2/\mathbb Q_2}
       (\delta^{-1}h(x,y)).
\]

Indeed, \(\delta^{-1}\mathcal O_2\) is the trace-dual of
\(\mathcal O_2\); the integrality of `B` therefore implies
\(h(M,M)\subset\mathcal O_2\).  Since `delta` is a 2-adic unit, its presence
does not affect the following parity argument.  The same trace-dual
description identifies the local discriminant module with the cokernel of

\[
h\colon M\longrightarrow M^*.
\]

More explicitly, if

\[
M^\#=\{z\in M\otimes_{\mathbb Z_2}\mathbb Q_2:
       h(z,M)\subset\mathcal O_2\},
\]

then perfect trace duality gives

\[
L^\vee\otimes\mathbb Z_2=M^\#,
\qquad
(L^\vee/L)_2=M^\#/M\cong\operatorname{coker}(h).
\]

### Evenness lemma

The evenness of `B` implies

\[
h(x,x)\in2\mathcal O_2
\qquad\text{for every }x\in M.
\]

Indeed, apply evenness both to `x` and to `phi x`.  If `u` is the residue of
`h(x,x)` and `c` is the nonzero residue of `delta^(-1)`, reduction modulo 2
gives

\[
\operatorname{Tr}_{\mathbb F_4/\mathbb F_2}(cu)=0,
\qquad
\operatorname{Tr}_{\mathbb F_4/\mathbb F_2}(c\phi^2u)=0.
\]

The two trace functionals are independent because `phi^2` is not in `F2`, so
`u=0`.

It follows that the reduction of `h` modulo 2 is alternating on the
11-dimensional vector space `M/2M` over `F4`.

### Local parity lemma

The preceding argument has a useful lattice-theoretic form independent of K3
surfaces.

> Let an even `Z_2`-lattice of rank `2n` carry a self-adjoint action of the
> unramified golden order `O_2`.  If it is 2-adically unimodular, then `n` is
> even.  The same conclusion holds if its discriminant is 2-elementary of
> type I.

In the unimodular case, `h mod 2` is nondegenerate alternating, so its
dimension `n` is even.  In the type-I 2-elementary case, the rank of `h mod 2`
is even by the evenness lemma, while the argument with the discriminant form
below shows that its nullity is even.  Their sum `n` is therefore even.

The supersingular-K3 application has `n=11`, so either local alternative is
impossible.  The next two sections spell out both contradictions in terms of
the geometric invariants.

## Odd characteristic

An alternating form in odd dimension is degenerate.  Hence `det(h)` is
divisible by 2 in `O_2`.  Because the local extension is unramified and its
different is a unit,

\[
v_2(\det_{\mathbb Z_2}B)
   =2v_{\mathcal O_2}(\det h).
\]

Thus the determinant of the underlying integral trace lattice is divisible
by

\[
N_{K_2/\mathbb Q_2}(2)=4.
\]

But a supersingular K3 lattice in odd characteristic has determinant
`-p^(2 sigma)`, which is odd.  This is a contradiction.

## Characteristic two

Assume now that `p=2`.  Let

\[
A=L^\vee/L\cong(\mathbb Z/2)^{2\sigma}.
\]

Self-adjointness extends `tau` from `L` to `L^vee` and hence to `A`.  Its
polynomial reduces to the irreducible polynomial `t^2+t+1` over `F2`.  Thus
`A` is an `F4`-vector space of dimension `sigma`.

The trace-dual identification above gives \(A=\operatorname{coker}(h)\).
Because `L` is
2-elementary, the elementary divisors of `h` over the discrete valuation ring
`O_2` are units or 2.  Since \(\dim_{\mathbb F_4}A=\sigma\), exactly `sigma`
of them are 2.  The nullity of `h mod 2` is therefore `sigma`, and

\[
\operatorname{rank}_{\mathbb F_4}(h\bmod2)=11-\sigma.
\]

The reduction is alternating, so its rank is even.  Hence

\[
\sigma\equiv1\pmod2. \tag{1}
\]

On the other hand, define the `F2`-valued discriminant pairing by

\[
b_A([x],[y])=2B(x,y)\pmod2,
\]

using \((\tfrac12\mathbb Z_2)/\mathbb Z_2\cong\mathbb F_2\).  Type I says
that the discriminant quadratic form is integral, so `B(x,x)` is integral for
\(x\in L^\vee\).  Hence the nondegenerate symmetric form `b_A` satisfies
`b_A(x,x)=0`; it is alternating.
Self-adjointness of `tau` gives

\[
b_A(\tau x,y)=b_A(x,\tau y).
\]

Finite-field trace duality defines `beta` uniquely by

\[
b_A(ax,y)=\operatorname{Tr}_{\mathbb F_4/\mathbb F_2}
                 (a\beta(x,y))
\qquad(a\in\mathbb F_4).
\]

The balance relation makes `beta` `F4`-bilinear; symmetry and nondegeneracy of
`b_A` make it symmetric and nondegenerate.

For every `x`, alternation applied to `x` and to `tau x` gives

\[
\operatorname{Tr}(\beta(x,x))=0,
\qquad
\operatorname{Tr}(\phi^2\beta(x,x))=0.
\]

Therefore `beta(x,x)=0`: `beta` itself is alternating.  A nondegenerate
alternating form over `F4` has even dimension, so

\[
\sigma\equiv0\pmod2. \tag{2}
\]

Equations (1) and (2) contradict each other.

## Why this matters for the search

The apparent positive-characteristic escape from the complex Picard-rank bound
is closed within the integral golden full-NS factory:

- complex K3: `rho <= 20`, so a quadratic eigenspace has dimension at most 10
  and its hyperbolic projectivization has `d <= 9`;
- finite-height K3 in positive characteristic: if the formal Brauer group has
  height `h`, the Artin--Mazur--Igusa inequality gives
  `rho <= 22-2h <= 20`, hence again `d <= 9`
  ([Ito, arXiv:1611.04721, Introduction](https://arxiv.org/abs/1611.04721));
  and
- supersingular K3: `rho=22`, but the required self-adjoint golden action is
  excluded by the theorem above.

Thus this precise K3 construction cannot produce a compact H10 example.  The
remaining geometric routes include a genuinely different algebraic object
(for example an OG10-type irreducible symplectic variety), a different field
construction, or a weaker correspondence/sublattice setup requiring a new
wall theorem.  Ordinary rank-22 Picard lattices still yield many
high-dimensional **noncompact** reflection polytopes; the obstruction is
specifically to the integral anisotropic golden full-NS eigenslice used here
to evade rational isotropy.

## Reproducibility

The finite-field core and a representative failed trace lattice are checked by
`scripts/certify_supersingular_golden_h10_obstruction.py`.  Run it from a clean
checkout as

```sh
PYTHONPATH=src python3 scripts/certify_supersingular_golden_h10_obstruction.py
```

The certificate is
only a guard against parity and convention errors; the proof above is the
mathematical completeness argument.
