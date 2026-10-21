# Borcherds inverse leverage: exact support, compensation, and the compact gap

**Status (2026-08-09): theorem-level audit and no-go.**  The determinant-line
construction gives a genuine `P -> type-IV` carrier, but it does not presently
give a compact polytope theorem.  Exact reflective support is not automatic,
and uncontrolled effective compensation is too weak for algebraic geometry to
restrict the input polytope.  The precise completion target is a bounded-slope
Hilbert-reflective theorem, not merely existence of another compensated
divisor.

This note distinguishes four assertions which must not be conflated:

1. the facet mirrors occur as traces of special divisors;
2. a weakly holomorphic Whittaker relation with the **pure** facet principal
   part exists;
3. a holomorphic automorphic section exists after adding effective divisors
   disjoint from the real slice; and
4. its global divisor is reflective, or has uniformly bounded degree relative
   to the selected facets.

The first is unconditional for the marked determinant extension.  The second
can fail.  A complex compensated Whittaker/Green relation is automatic; the
third requires the rational/integral q-structure with the arithmetic scope
stated in Section 3 below.  Only the fourth is strong enough for known
reflective-form classification or a facet bound, and it is open for the
compact high-dimensional branch.

## 1. Exact reflective support is not automatic

There is a direct determinant-line counterexample in the rational
finite-volume, noncompact norm-two class.  In the sign convention
`signature=(1,n)`, put

\[
 S=U\oplus E_8(-1)\oplus D_4(-1),\qquad \operatorname{sig}(S)=(1,13).
 \tag{1.1}
\]

Wang's classification of hyperbolic `2`-reflective lattices lists the sign
opposite of `S`; hence the complete `(-2)`-reflection chamber of `S` has
finite volume.  Since `det(S)=-4`, the canonical determinant line is

\[
 K_S=\langle 2(-1)^{13}\det(S)\rangle=\langle8\rangle,
 \qquad \widehat S=S\oplus\langle8\rangle,
 \tag{1.2}
\]

and `widehat S` has signature `(2,13)`.  Wang's complete classification of
`2`-reflective lattices in signature `(n,2)`, `n>=8`, says that none exists
for `n=13`.  Therefore `widehat S` admits no nonconstant holomorphic modular
form whose zero divisor is supported on its `(-2)` mirrors.  In particular,
the determinant extension of this finite-volume Coxeter chamber has no exact
global `2`-reflective automorphic divisor carrying its facet mirrors.

This is stronger and more directly relevant than the often-cited statement
that eighteen hyperbolic `2`-reflective lattices have no automorphic
correction.  In that statement Wang defines an automorphic correction using
`U(m)\oplus S`, a rank-two extension.  It does **not**, by itself, rule out the
rank-one determinant extension `S\oplus K_S`.  Equations (1.1)--(1.2), combined
with the later forty-two-lattice classification, do rule it out for this
specific determinant extension.

There is also a useful positive control.  For

\[
 S_0=II_{1,9}=U\oplus E_8(-1),
 \qquad \widehat S_0=S_0\oplus\langle2\rangle,
 \tag{1.3}
\]

the two even lattices

\[
 U\oplus E_8(-1)\oplus\langle2\rangle
 \quad\text{and}\quad
 2U\oplus E_7(-1)
 \tag{1.4}
\]

have the same signature and discriminant form and are isometric by the
standard uniqueness theorem for indefinite even lattices in this rank.
Wang's list identifies `2U+E7(-1)` as **complete** `2`-reflective.  Its
complete reflective divisor pulls back on the determinant real slice to the
complete `(-2)` arrangement of `S_0`: every root with a nonzero auxiliary-line
component misses the whole slice.  The weight formula gives weight `165` for
the multiplicity-one complete form.  Thus the classical ten-facet `E10`
finite-volume noncompact `H9` simplex has a positive determinant-line
reflective realization.  This is a different realization of a known
polytope, not a new abstract polytope.  The complete divisor also contains
non-simple `S_0` root mirrors which can touch lower-dimensional boundary
faces, so this is an exact complete-arrangement/same-chamber return, not a
claim that one coarse divisor column has facet-only closed-slice support.

The pair (1.1), (1.3) proves the right qualitative statement:

\[
 \boxed{\text{determinant-line exact reflectivity is selective, not automatic}.}
 \tag{1.5}
\]

It does not settle the compact totally-real class.  The negative example is
rational and noncompact, while the positive example is likewise noncompact.
For Bugaenko's compact `H8` chamber the project certificate only proves that
the pure standard norm-two Whittaker column is obstructed.  It does not
exclude an arbitrary reflective modular form, because no converse theorem is
known in the compact anisotropic Hilbert setting.

## 2. Why rational reflective classifications do not bound compact `H4+`

A compact arithmetic hyperbolic chamber of dimension at least four cannot be
rational: by Meyer, an indefinite rational quadratic form of rank at least
five is isotropic.  Its field of definition therefore has degree `d>1`.
After the determinant extension, the relevant Shimura datum is

\[
 \operatorname{Res}_{F/\mathbb Q}SO(\widehat V),
 \qquad
 \operatorname{sig}(\widehat V)=((2,n),(0,n+2),\ldots,(0,n+2)).
 \tag{2.1}
\]

It has one noncompact type-IV factor and a projective quotient.  This is not
the classical integral `O(2,n)` situation classified by Ma and Wang.  The
trace quadratic lattice has signature

\[
 \left(2,\;n+(d-1)(n+2)\right),
 \tag{2.2}
\]

but the restriction-of-scalars arithmetic group is an arithmetic lattice in
the proper subgroup `Res SO`, not a finite-index subgroup of the full
orthogonal group of (2.2).  Consequently the classical `2`-reflective bounds
and forty-two-lattice list cannot be applied to the compact Hilbert-orthogonal
quotient.

No primary-source classification or finiteness theorem was found for
reflective Hilbert-orthogonal modular forms on the anisotropic varieties
(2.1).  This is an independent missing theorem, even if exact support could
be supplied.

## 3. What compensation really gives

At a direct determinant level, an auxiliary-nonzero discriminant component
has no vector orthogonal to the auxiliary line.  Its special divisor is
therefore disjoint from the complete canonical real slice.  The project's
safe-column argument proves that these Fourier columns span the Whittaker
obstruction dual.  Hence every finite auxiliary-zero target principal part
has a finite complex slice-disjoint correction.  With a real q-structure the
safe columns positively span; after denominator clearing an integral
q-structure gives an effective Borcherds divisor.

The arithmetic scope is:

| branch | compensated carrier status | unresolved issue |
|---|---|---|
| `F=Q`, interior | integral q-structure is supplied by McGraw, so effective interior compensation is unconditional | a noncompactification Borcherds product may have boundary poles; cusp-uniform exact-sheet bounds are also needed for a noncompact chamber |
| `F/Q` totally real, degree `>1` | the quotient is compact, so an integral weakly holomorphic input gives a holomorphic section with no cusp boundary | Bruinier states that McGraw's integral-basis argument adapts to the Hilbert--Weil space, but does not write the adaptation; the integral upgrade remains conditional on that input |
| fixed compact `P` | target-sheet separation, finite compensation, and a finite search cutoff are available per input | level, Hilbert detection, obstruction-matrix height, weight, and volume are not uniform over an independently defined class |

This is exact **on the marked real slice**, not exact in global divisor
support.  Known reflective-form classification theorems require the global
zero divisor to be supported on reflective mirrors.  They cannot be applied
after arbitrary safe Heegner compensation merely because the compensation
misses one real slice.

## 4. Ample-completion no-go

The failure of uncontrolled compensation to give AG leverage is structural,
not a missing clever BMY calculation.

### Proposition 4.1 (unbounded compensation is AG-tautological)

Let `X` be a projective integral variety, let `A` be an ample Cartier divisor,
and let `T` be any effective Cartier divisor.  For every sufficiently large
integer `k` there is an effective Cartier divisor `C_k` such that

\[
                  T+C_k\sim kA.
 \tag{4.1}
\]

#### Proof

For `k` sufficiently large, the Hilbert polynomial of
`O_X(kA-T)` is positive and higher cohomology vanishes by Serre vanishing.
Thus `H^0(X,O_X(kA-T))` is nonzero.  The divisor of a nonzero section is an
effective `C_k` satisfying (4.1).  `square`

On a compact determinant-line carrier, take `A=lambda` after passing to a
level on which the automorphic line is integral and take `T` to be the sum of
the separated facet components.  Proposition 4.1 says that the bare AG
assertion

```text
there is an effective divisor in a power of lambda containing every facet
```

holds for **every** finite target arrangement, independently of Coxeter
geometry.  Borcherds compensation adds arithmetic special support and real
slice disjointness, which are meaningful carrier data, but ordinary
projective classification cannot turn unbounded `C_k` and `k` into a bound on
`P`.

The same obstruction appears numerically.  On an exact separated level,

\[
 \#\operatorname{Facets}(P)q_{\min}\le w\lambda^n.
 \tag{4.2}
\]

An upper estimate for `w lambda^n` which is linear in the total target
multiplicity gives only `s<=constant*s`; it cannot bound the facet number
`s`.  A useful theorem must bound the normalized slope

\[
 U(P)=\frac{w\lambda^n}{q_{\min}}
 \tag{4.3}
\]

independently of `s`, or must force the global auxiliary divisor into a
separately classified bounded family.  Ordinary BMY and its natural
coefficient-one/orbifold mirror variants have already been computed in the
companion reports and do neither.

## 5. Exact completion theorem needed

The Borcherds route would become genuine
`POLYTOPE -> AG -> POLYTOPE` leverage on a non-singleton compact class only
after both of the following are proved.

1. **Bounded carrier theorem.**  For an independently defined compact
   arithmetic Coxeter class `C`, every `P in C` has either

   - a globally reflective determinant-line form with uniformly bounded
     mirror multiplicity, or
   - a positive safe compensation with a class-uniform, facet-count-independent
     bound on (4.3),

   with the Hilbert integral q-structure, level, detection height, matrix
   height, and Prasad/local volume all controlled.

2. **Hilbert-reflective AG theorem.**  For the anisotropic
   restriction-of-scalars varieties (2.1), globally reflective forms of the
   bounded multiplicity in item 1, or compensated divisors of the bounded
   slope in (4.3), satisfy an independent finiteness/dimension theorem with an
   effective list of terminal survivors.

Then (4.2) gives a facet bound directly, while a large-weight Whittaker
obstruction trapped inside the same bounded support could give a dimension
bound.  The terminal list must still be run through literal Vinberg
compactness/finite-volume tests.  Neither item is presently available.

The frozen inverse verdict is therefore

```text
C        = PASS  (marked compact arithmetic determinant-eligible class)
Carrier  = conditional PASS for effective Hilbert compensation;
           unconditional at the complex Whittaker/Green level
Return   = marked exact-real-slice PASS, not an intrinsic nef/Mori return
AG       = FAIL  (no applicable compact Hilbert-reflective classification)
Pushback = FAIL  (ample completion and linear target mass make the current
                   effective-divisor statement nonrestrictive)
```

This is a rigorous carrier theorem plus a rigorous missing-lemma/no-go result,
not a new compact dimension or facet bound.  The `E10` row is a positive
known-polytope calibration; no new compact survivor follows from this audit.

## 6. Primary sources and project certificates

1. J. H. Bruinier,
   [*Regularized theta lifts for orthogonal groups over totally real fields*](https://arxiv.org/abs/0908.3076),
   especially Proposition 4.5, Theorem 6.8, and the integral-basis sentence in
   the proof of Theorem 7.1.
2. H. Wang,
   [*The classification of 2-reflective modular forms*](https://arxiv.org/abs/1906.10459),
   especially Theorem 7.2 and its definition of `U(m)+S` automorphic
   correction.
3. H. Wang,
   [*2-reflective lattices of signature (n,2) with n>=8*](https://arxiv.org/abs/2301.11536),
   Theorem 1.2 and Corollary 4.7: exactly forty-two `2`-reflective lattices in
   the stated range, none in signature `(13,2)`, and the complete
   `2U+E7(-1)` positive row.
4. S. Ma,
   [*Finiteness of 2-reflective lattices of signature (2,n)*](https://arxiv.org/abs/1409.2969),
   for the classical rational `O(2,n)` finiteness theorem and its scope.
5. I. Metzler,
   [*A new converse theorem for Borcherds products*](https://arxiv.org/abs/2512.20701),
   whose converse assumes a hyperbolic split and therefore does not cover the
   compact anisotropic Hilbert branch.
6. Project details:
   [determinant transfer](polytope_to_automorphic_divisor_transfer_theorem.md),
   [Bugaenko Whittaker gate](bugaenko_h8_whittaker_gate.md),
   [compact BMY no-go](compact_shimura_ci_bmy_feedback_no_go.md), and
   [mirror-stratum BMY no-go](compact_mirror_multiple_strata_bmy_no_go.md).
