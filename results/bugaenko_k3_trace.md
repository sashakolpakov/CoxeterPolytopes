# K3 Trace Realization of Bugaenko's Field Lattice

Status: exact derived result; the field-lattice normalization has been checked
against a pinned Vinberg regression fixture and independent secondary sources.

## Exact input

Bugaenko's compact \(\mathbb H^8\) lattice is

\[
L_8=\langle-2\phi\rangle\oplus E_8
\]

over \(\mathcal O_F=\mathbb Z[\phi]\), with \(F=\mathbb Q(\sqrt5)\).  The
exact ambient Gram matrix and eleven primitive roots are archived in
`data/known_examples/bugaenko_compact_h8.json`.  They agree with Bottinelli's
pinned Vinberg-algorithm regression fixture and Guglielmetti's lattice
normalization and Coxeter transcript.  The 1992 primary chapter is
access-gated in the current environment, a provenance limitation recorded
separately from the exact calculation.

## Trace construction

Take the totally negative codifferent element

\[
\alpha=-\frac{\phi}{\sqrt5}=-\frac25-\frac15\phi.
\]

On the underlying rank-18 \(\mathbb Z\)-module define

\[
B(x,y)=\operatorname{Tr}_{F/\mathbb Q}
\bigl(\alpha (x,y)_{L_8}\bigr).
\]

In the interleaved basis \(e_i,\phi e_i\), exact computation gives an even
integral lattice with

\[
\operatorname{inertia}(B)=(1,17),\qquad \det B=-4.
\]

The determinant-one basis change in
`scripts/bugaenko_trace_lattice.py` verifies

\[
(N,B)\simeq
\langle2\rangle\oplus\langle-2\rangle
\oplus E_8(-1)\oplus E_8(-1).
\]

## Primitive K3 embedding

Write the K3 lattice as

\[
\Lambda_{K3}=U_1\oplus U_2\oplus U_3
\oplus E_8(-1)\oplus E_8(-1),
\]

with bases \(e_i,f_i\) for \(U_i\).  Embed the first two summands of \(N\) by

\[
\langle2\rangle\ni 1\longmapsto e_1+f_1,
\qquad
\langle-2\rangle\ni 1\longmapsto e_2-f_2,
\]

and embed both \(E_8(-1)\) summands identically.  The script verifies the Gram
matrix and gives an explicit extension of these 18 columns to a determinant
\(-1\) basis of \(\Lambda_{K3}\), proving primitivity.  The orthogonal
complement is

\[
N^\perp\simeq\langle-2\rangle\oplus\langle2\rangle\oplus U,
\]

of signature \((2,2)\).

By the period theory of lattice-polarized K3 surfaces, this primitive
hyperbolic embedding yields a nonempty 2-dimensional K3 moduli stratum, and a
generic period has Néron--Severi lattice exactly \(N\).

## Consequence and continuation

The arithmetic field lattice passes the trace-lattice parity, signature,
discriminant, and primitive K3 embedding gates.  The subsequent exact
five-shell and Vinberg-completion argument now proves, existentially up to a
K3 Weyl-conjugate marking, that its eigen-nef section is Bugaenko's 11-wall
chamber.  It also enumerates every K3 root hyperplane meeting the closed
section.  See `results/bugaenko_k3_eigen_nef.md`.

The extra structure is a self-adjoint ring embedding
\(\mathbb Z[\phi]\hookrightarrow\operatorname{End}_{\mathbb Z}\operatorname{NS}(X)\),
not an action by automorphisms: multiplication by \(\phi\) does not preserve
the intersection form.  The result is existential as a marked K3 construction;
no explicit projective equation or Weyl word has been archived.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/bugaenko_trace_lattice.py
```

The positive and negative assertions are also covered by
`tests/test_bugaenko_trace.py`.
