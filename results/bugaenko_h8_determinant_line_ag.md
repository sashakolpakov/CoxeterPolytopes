# Bugaenko compact `H8` as a canonical determinant-line AG arrangement

**Certified status (2026-08-08):** starting from the archived Bugaenko
polytope lattice and its exact compact Vinberg chamber, the coordinate-free
determinant-line compiler produces a unique even glue package of signature
`(2,8)`.  Its proper orthogonal Shimura stack contains the Bugaenko chamber
as a canonical real slice, its special mirrors recover all eleven facets, and
its stacky codimension-two strata recover every Coxeter label.  This is a new
project-level algebro-geometric realization of the known Bugaenko polytope,
not a new abstract Coxeter polytope and not an automorphic-principality claim.
The separate exact Whittaker calculation now proves that the pure norm-two
principal part and complete first compensator shell are obstructed.  The
general determinant safe-span theorem nevertheless proves abstractly that a
larger finite slice-disjoint compensation exists; with Bruinier's stated
Hilbert integral q-structure it can be chosen nonnegative integral and gives
an effective carrier.  Its least support and quantitative invariants remain
unknown.

The finite hidden-sheet gate gives the important exact split

\[
  320=11\ \text{facet mirrors}
      +309\ \text{nonfacet lower-boundary contacts}
      +0\ \text{interior cuts}.
\]

Thus the full numerical Heegner type does **not** pull back to a facet-only
divisor on the closed chamber.  The 309 extra sheets do **not** cut or refine
the open chamber.

## 1. Input and reproducibility boundary

The sole Coxeter-lattice input is
[`bugaenko_compact_h8.json`](../data/known_examples/bugaenko_compact_h8.json),
SHA-256
`c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`.
It supplies the rank-nine `Z[phi]` Gram matrix, eleven simple roots, all 41
finite vertices, no cusps, the `f`-vector

\[
 (41,164,316,374,294,156,54,11,1),
\]

and a verified exact Vinberg completion.  The compiler does not import the
hard-coded Bugaenko example module.

For the compact mirror test, the compiler independently repeats finite-root
closure at every vertex directly from this input.  It then compares the
answer with the existing
[`bugaenko_k3_eigen_nef.json`](../data/trace_lattices/bugaenko_k3_eigen_nef.json)
archive, SHA-256
`e77bb81bcc738f4f542676748dbcf09418cab830da579a600fd1e309f03de782`.
The two enumerations agree exactly, including the mirror-vector digest
`815c0c4cad8e4ab9da32ef491c104ac71a781d4f2b4300f988ecb484506da69e`.

The replay is
[`certify_bugaenko_h8_determinant_line_ag.py`](../scripts/certify_bugaenko_h8_determinant_line_ag.py),
the machine certificate is
[`bugaenko_h8_determinant_line_ag.json`](../data/enumerations/bugaenko_h8_determinant_line_ag.json),
and the focused tests are
[`test_bugaenko_h8_determinant_line_ag.py`](../tests/test_bugaenko_h8_determinant_line_ag.py).
The certificate SHA-256 is
`dce7d52e301b5e215a6f4c58f55ea9fca9114b77c2b84fc929c0c36ac21e6aad`.

## 2. Canonical arithmetic extension

Put `O=Z[phi]`, where `phi^2=phi+1`, and negate the archived form.  The
normalized hyperbolic lattice is

\[
 L=\langle2\phi\rangle\perp E_8(-1),\qquad
 \det L=2\phi.
\]

Its signatures at the identity and conjugate embeddings are respectively
`(1,8)` and `(0,9)`.  For `n=8`, its coordinate-free determinant line is

\[
 K_L=\det_O L,
 \qquad
 b_K(\wedge v_i,\wedge w_i)=2(-1)^8\det(b_L(v_i,w_j)).
\]

In the generator induced by the archived basis, `K_L=<4 phi>`.  This
description is basis-independent: changing the basis by `A` scales the
determinant-line generator by `det(A)` and its Gram entry by `det(A)^2`.
Consequently

\[
 \widehat L_0=L\perp K_L
   =\langle2\phi\rangle\perp E_8(-1)\perp\langle4\phi\rangle
\]

is even integral, has signatures `(2,8)` and `(0,10)`, and has determinant
`8 phi^2`.  Since `phi` is a unit, its determinant ideal is `(8)` and its
discriminant module has order 64.  Definiteness at the conjugate place proves
that the quadratic space is anisotropic over `Q(sqrt(5))`; hence its
arithmetic type-IV quotient is compact and its coarse quotient is
projective.

### Complete gluing enumeration

The `E8(-1)` summand is unimodular, so

\[
 A_{\widehat L_0}
  \cong O/(2\phi)\oplus O/(4\phi)
  \cong O/(2)\,u\oplus O/(4)\,v.
\]

As an abelian group this is `(Z/2)^2 + (Z/4)^2`.  Its finite quadratic form is

\[
 q(u)=\frac1{4\phi},\qquad q(v)=\frac1{8\phi},\qquad
 q(au+bv)=\frac{2a^2+b^2}{8\phi}\pmod O.
\]

The replay enumerates all `4*16=64` elements using exact pair arithmetic
modulo 2, 4, and 8.  The zero element is the only isotropic element.
Therefore the zero submodule is the only totally isotropic `O`-submodule and

\[
 \boxed{\widehat L_0\text{ is the unique even reflection-stable determinant-line
 gluing.}}
\]

This conclusion is exhaustive, not a bounded search.  Its enumeration digest
is `bbb3975b9a818fe3b4a5bbecb88050aade28009a330aa8450282e7fe64356047`.

## 3. Roots, Heegner type, and the canonical real slice

Every archived wall vector becomes a primitive square `-2` root of `L` and
extends as `(r_i,0)` in `Lhat_0`.  Exact matrix checks prove that all eleven
reflections are integral involutive isometries.  Every root has divisibility
ideal `O`, witnessed by an explicit unit pairing in the certificate.
Moreover every root reflection acts trivially on the discriminant module:
for `x in L^vee`,

\[
 s_r(x)-x=(x,r)r\in L,
\]

and it fixes the determinant line pointwise.

Choose `e in K_L` with `e^2=4 phi` at the identity embedding.  Then

\[
 j_e:\mathbb H(L)\longrightarrow\mathcal D^+(\widehat L_0),
 \qquad [x]\longmapsto\operatorname{span}^+(e,x)
\]

is the canonical real slice, and for every facet root

\[
 j_e^{-1}(\mathcal D_{(r_i,0)})=r_i^\perp.
\]

All eleven roots have one numerical special-divisor type: norm `-2`,
divisibility one, and discriminant coset zero.  After reversing the sign to
Bruinier's `(8,2)` convention, this is `(m,mu)=(1,0)`.  Numerical equality of
type does not prove that the eleven mirrors form one global arithmetic orbit
or one irreducible divisor component; no such claim is made.

## 4. Stacky recovery of the Coxeter labels

The 55 facet pairs split as follows.

| relation | pair/ridge count | distinguished generic inertia order |
|---|---:|---:|
| `I2(2)` | 44 | 4 |
| `I2(3)` | 8 | 6 |
| `I2(5)` | 2 | 10 |
| ultraparallel | 1, facets 9 and 10 | no real-slice crossing |

The 54 finite pairs equal the 54 ridges reconstructed from the complete
vertex list; this is not merely a count equality.  At a generic two-mirror
stratum the arithmetic stack's inertia contains the displayed distinguished
dihedral group.  It may contain further stabilizers, so equality with the
whole inertia is not asserted.  A neat coarse cover removes this torsion;
the labelled Coxeter data therefore live naturally on the stack.

## 5. The finite hidden-sheet gate

Because the chamber is a complete compact reflection chamber, a lattice root
mirror cannot cut its interior.  Any such mirror meeting the closed chamber
therefore contains a boundary face, and that face contains a vertex.  At each
of the 41 vertices the mirrors through it are exactly the roots in the finite
standard parabolic generated by the eight incident facets.  Exact reflection
closure of those 41 finite root systems is consequently a complete
enumeration, with no metric cutoff.

The contact dimensions of the 309 nonfacet mirrors are

| contact-face dimension | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|---:|
| nonfacet mirrors | 88 | 43 | 26 | 17 | 96 | 25 | 14 |

The remaining eleven mirrors contact seven-dimensional faces and are exactly
the facets.  Therefore:

- the complete `(m,mu)=(1,0)` Heegner shell has 320 restricted mirrors
  meeting the closed chamber;
- its 309 extras touch only lower-dimensional boundary faces;
- none cuts the open chamber; and
- the global numerical Heegner divisor is not a facet-only divisor on the
  closed real slice.

The last distinction matters.  The chamber remains a complete Coxeter
chamber, but selecting the entire numerical Heegner type remembers more of
the ambient reflection arrangement than just its facet boundary.

## 6. The associated algebraic object and an associated surface

Let `Gamma` be an arithmetic subgroup of `O^+(Lhat_0)` containing the eleven
extended reflections, for example a suitable effective version of the full
integral group.  The canonical carrier is the proper arithmetic
Deligne--Mumford stack

\[
 \mathscr X_\Gamma=[\mathcal D^+(\widehat L_0)/\Gamma],
 \qquad \dim_\mathbb C\mathscr X_\Gamma=8.
\]

Its coarse space is projective.  After passage to a neat finite-index cover,
six sufficiently general ample divisors have a smooth projective
complete-intersection surface.  Restricting the special divisors gives
algebraic curves on that surface.  Alternatively, a stacky complete
intersection retains the dihedral inertia labels.

The eightfold stack is canonical relative to the marked even lattice and a
choice of arithmetic group.  The extracted surface is not unique: choosing a
neat level, polarization powers, and six sections is auxiliary.  Thus this is
an unconditional `P -> AG arrangement -> surface` construction, but not yet
a canonical classical surface whose standard invariants classify `P`.

## 7. Exact transfer criterion exposed by this calibration

For this marked arithmetic polytope, the following statements are now exact.

1. **Existence of the projective AG carrier:** equivalent here to the verified
   one-positive-place signature and determinant-line construction; it passes.
2. **Integral extension of every facet reflection:** decided by the root and
   discriminant action checks; it passes for all eleven roots.
3. **Recovery of mirrors and angles:** the real-slice identity plus stacky
   inertia recovers all facets and all labels; it passes.
4. **Recovery of actual faces:** negative-definite divisor intersections give
   candidates, while the archived oriented halfspaces and vertex feasibility
   select the true faces; it passes because the full chamber is certified.
5. **Facet-only closed-slice divisor for a chosen collection of Heegner
   types:** equivalent to the finite hidden-sheet set being exactly the facet
   set.  It fails for the full numerical type `(1,0)` by the 309 boundary
   contacts.
6. **No refinement of the open chamber:** equivalent to zero interior-cutting
   sheets.  It passes.
7. **Automorphic principality:** requires the separate vector-valued
   Hilbert-modular/Whittaker obstruction pairing to vanish.  It is not
   implied by any preceding check.  The follow-up calculation proves that
   the pure standard norm-two column is obstructed by a nonzero weight-`(5,5)`
   Poincare series.  The complete globally minimal compensator shell is also
   obstructed: its divisor is disjoint from the slice, but the corresponding
   target-plus-shell `2x2` matrix is invertible.  A larger finite safe support
   exists, and exact coefficient-rank enumeration gives a terminating
   degree-two-field trace-height cutoff `max(R_det,B_CLW)`.  Its numerical
   least support, `R_det`, and obstruction matrix remain undecided; once that
   matrix is known, the quantitative compression bounds conic denominators,
   weight, and degree.  The facet-separating level is computable from the
   marked chamber.  See the [`exact Whittaker gate`](bugaenko_h8_whittaker_gate.md),
   [`first compensator theorem`](bugaenko_h8_whittaker_first_compensator.md),
   [`orbit-square theorem`](orbit_square_quantitative_safe_height.md), and
   [`weight/degree theorem`](borcherds_weight_degree_quantitative_bound.md).

What transfers back to the polytope is therefore its marked mirrors, complete
oriented face data, and Coxeter labels.  What does not transfer automatically
is a facet-only global divisor, an automorphic product, a canonical surface,
a nef cone, a divisor intersection matrix equal to the Lorentzian Gram, or a
new compactness/dimension theorem.  Those are precisely the remaining gates
for turning this uniform association into the stronger necessary-and-
sufficient surface classification envisioned in the inverse program.
