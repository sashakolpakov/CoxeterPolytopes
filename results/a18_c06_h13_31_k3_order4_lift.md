# A genuine order-four K3 lift of the 31-wall hyperbolic 13-polytope

## Theorem

Let (N) be Alexeev's even 2-elementary K3 Picard lattice with invariants
((r,a,\delta)=(18,4,0)), and let (A) be the class-`C06` weighted-diagram
involution of its complete 48-wall mixed-root reflection chamber.  There is a
complex projective K3 surface (X) and a non-symplectic automorphism
(g\in\operatorname{Aut}(X)) of order four such that

1. `NS(X)` is isometric to (N);
2. (g^*|_{NS(X)}=A), so the action on the Picard lattice has order two and
   fixed rank 14; and
3. the section of Alexeev's chamber by (NS(X)^A\otimes\mathbb R) is the
   already certified finite-volume, noncompact Coxeter polytope with 31
   facets in (\mathbb H^{13}).

Thus `A18-C2c-H13-31` is no longer merely a diagram-symmetry section of a K3
Picard lattice: it has a genuine K3-automorphism realization.  This is a new
project-level algebraic realization of the existing chamber, not a claim that
the abstract polytope, its arithmetic lattice, or dimension 13 is new.

## Exact discriminant lift

The archived `C06` Picard action satisfies

\[
 A^2=1,\qquad AGA^t=G,\qquad \operatorname{rk}N^A=14.
\]

The discriminant group is ((\mathbb Z/2)^4).  In the frozen half-vector
basis in the certificate, the induced action is

\[
 M=\begin{pmatrix}
 0&0&1&0\\
 0&0&0&1\\
 1&0&0&0\\
 0&1&0&0
 \end{pmatrix}.
\]

Take

\[
 T=U(2)\oplus U(2),\qquad
 G_T=\begin{pmatrix}
 0&2&0&0\\2&0&0&0\\0&0&0&2\\0&0&2&0
 \end{pmatrix}
\]

and the order-four isometry

\[
 B=\begin{pmatrix}
 0&0&-1&0\\
 0&0&0&-1\\
 1&0&0&0\\
 0&1&0&0
 \end{pmatrix}.
\]

Exactly,

\[
 B^2=-1,\qquad B^4=1,\qquad BG_TB^t=G_T.
\]

The discriminant identification is represented over (mathbb F_2) by

\[
 P=\begin{pmatrix}
 0&0&1&1\\
 0&0&0&1\\
 1&1&0&0\\
 0&1&0&0
 \end{pmatrix}.
\]

The certificate checks all 16 discriminant classes and proves both

\[
 q_N(x)=-q_T(xP),\qquad MP=P(B\bmod 2).
\]

Here the first equality is written as equality modulo two; because
(\delta=0), the discriminant values are integral and negation has the same
value modulo two.

Glue (N\oplus T) along the graph of (P).  An exact Hermite-normal-form
replay gives a rank-22 lattice that is integral, even, unimodular, and has
signature ((3,19)).  It is therefore the K3 lattice.  Equivariance makes
(A\oplus B) integral on the glued lattice.  Its exact characteristic
polynomial is

\[
 (x-1)^{14}(x+1)^4(x^2+1)^2,
\]

so its order is exactly four.

## Period, exact Picard lattice, and the Kähler chamber

In row coordinates for (T\otimes\mathbb C), take

\[
 \omega=(1,\sqrt2,i,i\sqrt2).
\]

Then

\[
 \omega B=i\omega,\qquad \omega^2=0,\qquad
 \omega\overline\omega=8\sqrt2>0.
\]

If a rational vector ((a,b,c,d)\in T\otimes\mathbb Q) is orthogonal to
(\omega), its real and imaginary equations are

\[
 a\sqrt2+b=0,\qquad c\sqrt2+d=0.
\]

Hence it is zero.  The corresponding period has no extra rational algebraic
class, and the graph gluing is primitive, so the resulting K3 surface has
`NS(X)=N` exactly.

There is also an exact Kähler-class check.  The `C06` action permutes all 48
source walls.  It fixes a rational chamber class of square 385 whose pairing
with every source wall is positive, with minimum one.  Its double is the
primitive integral class

```text
(128,-165,-308,-284,-258,-223,-186,-126,-64,
 93,208,232,258,293,330,390,194,258)
```

of square 1540.  The complete mixed-root chamber lies inside a unique
`(-2)`-Weyl chamber, and the action preserves it.  After choosing the marking
with this Weyl chamber as the Kähler chamber, the displayed integral class is
ample.  The glued order-four action is therefore an effective Hodge isometry
preserving a Kähler class.  Strong global Torelli realizes it by a genuine K3
automorphism.  Since its action on (H^2) has order four, so does the
automorphism; since it acts on (H^{2,0}) by (i), it is non-symplectic.

## Hyperbolic conclusion uses Vinberg's criterion

The algebraic lift does not replace the volume proof.  The pinned Vinberg
certificate independently checks the literal edge-endpoint and affine
zero-neighborhood criteria for the 31-wall section:

- 28,766 elliptic edge diagrams;
- 3,433 finite vertices and 238 ideal vertices;
- every elliptic edge has exactly two finite-or-ideal endpoints;
- every maximal affine zero-neighborhood has the required cusp rank; and
- zero bad edges and zero bad parabolic neighborhoods.

Consequently the (\mathbb H^{13}) polytope has finite volume and is
noncompact.  The positive cusp count, not the finiteness of the wall list,
proves noncompactness.

## Scope and lineage

Purely non-symplectic order-four K3 automorphisms, including the case in
which the square acts trivially on the Picard lattice, are part of the
classical framework of Artebani--Sarti,
[*Symmetries of order four on K3 surfaces*](https://arxiv.org/abs/1102.4436).
Accordingly, no novelty is claimed for the automorphism order or its general
cohomological type.  The new project-level statement is the explicit `C06`
action on this Picard lattice, its equivariant discriminant glue, and its
preservation of the particular 48-wall chamber whose fixed section is the
31-wall hyperbolic polytope.

The result closes the previously explicit *geometric-lift* gap for the
single class `A18_4_0_C06`.  It does not turn the 48-wall mixed-root chamber
into the literal nef cone: 24 of its walls have norm `-4`.  The precise
geometric object is an automorphism-invariant reflection chamber in the
divisor lattice and its fixed hyperbolic section.  No quotient-surface nef
interpretation is asserted.

Order-four non-symplectic K3 automorphisms themselves are classical; see
Artebani--Sarti, [*Symmetries of order four on K3
surfaces*](https://arxiv.org/abs/1102.4436).  No novelty is claimed for that
general automorphism type.  The present contribution is the exact `C06`
Picard action, its discriminant-compatible lift, preservation of this
specific 48-wall chamber, and the resulting `H13/31` fixed section.

The abstract `H13/31` chamber and its finite-volume proof were already in the
project ledger.  Its maximal reflection-stable overclosure is the known
determinant-four lattice with the 15-wall Tumarkin/Scharlau--Kirschmer full
chamber.  The worldwide priority of the literal 31-wall finite-index subgroup
chamber remains open.  The new statement here is exactly the K3
automorphism realization.

## Reproduction

The permanent files are:

- [exact lift certificate](../data/enumerations/a18_c06_h13_31_k3_order4_lift.json);
- [certificate generator](../scripts/certify_a18_c06_h13_31_k3_order4_lift.py);
- [regression test](../tests/test_a18_c06_h13_31_k3_order4_lift.py);
- [original construction and subgroup certificate](../data/enumerations/alexeev_remaining_high_fixed_section_construction.json);
- [Vinberg finite-volume certificate](../data/enumerations/alexeev_remaining_high_fixed_section_vinberg.json); and
- [global abstract/arithmetic lineage certificate](../data/enumerations/a18_c06_h13_31_global_lineage.json).

Run from the repository root:

```bash
PYTHONPATH=src:. python3 scripts/certify_a18_c06_h13_31_k3_order4_lift.py
PYTHONPATH=src:. pytest -q tests/test_a18_c06_h13_31_k3_order4_lift.py
```

The certificate payload SHA-256 is
`c6dd715b4667021a7fa7a01189ed37b68ec20617f1ca58cf0cf141432278845e`.
The frozen JSON file SHA-256 is
`f6c271e87188983fd4acc0d9d528008d64f8afab988b4333a51db7f3af809f6c`.
