# Exact AGCat return and global finiteness of maximal arithmetic Coxeter chambers

**Status (2026-08-12): theorem, pending one final independent composition
audit.** This report gives the requested composite result. It constructs an
exact functor from arithmetic hyperbolic Coxeter chambers to a genuine
real-framed algebraic category, uses integral-point geometry on that carrier
to prove that only finitely many ground fields occur, and pushes this back to
finiteness of reflection-maximal chambers and reflection groups.

The proof is **AG-mediated**, not AG-only. Two classical geometric inputs are
stated explicitly: the universal dimension ceiling for finite-volume
hyperbolic reflection polyhedra and the bounded connected-wall-frame lemma.
No conclusion of the Agol--Belolipetsky--Storm--Whyte finiteness theorem, no
global reflective-lattice finiteness theorem, and no pre-existing global
ground-field finiteness theorem is used.

## 1. Source and target categories

Let \(\mathsf{Cox}^{\mathrm{mark}}_{\mathrm{arith}}\) be the groupoid of
finite-sided finite-volume arithmetic hyperbolic Coxeter chambers \(P\), with
their minimal Vinberg field \(k\), distinguished real place, ordered inward
facet-root rays, and chosen positive component. Morphisms are marked
field-compatible projective similitudes. Arithmetic reflection groups are in
the standard admissible Vinberg quadratic/lattice-form model; this is the
model canonically recovered from their reflection representation.  More
precisely, Vinberg's structural reflection lemma (recorded as ABSW Lemma 5.2,
citing Vinberg Lemma 7) says that an arithmetic hyperbolic lattice containing
real hyperplane reflections is defined by such a quadratic form.  This is a
structural reduction, not a field, dimension, or finiteness theorem.

Let \(\mathsf{AGCat}^{\mathbf R}_{\mathrm{arr}}\) be the groupoid of tuples,
defined over a displayed totally real base field \(k\),

\[
 (\mathbf P(V),Q,H_1,\ldots,H_f,\sigma_0,C),                         \tag{1.1}
\]

where \(Q\) is a smooth projective quadric, the \(H_i\) are ordered polar
hyperplanes, \(\sigma_0\) is a distinguished real place, and \(C\) is a
chosen component with oriented sides. This is an algebraic arrangement
category with finite real descent/framing data. The framing is necessary:
bare complex hyperplanes do not remember which half-spaces define \(P\).
Morphisms are field-isomorphism-compatible projective isomorphisms carrying
\(Q\), every ordered \(H_i\), the distinguished place \(\sigma_0\), and the
oriented component \(C\) to the corresponding target data. The return
functor on morphisms is restriction of that projective isomorphism to the
recovered hyperbolic half-space chamber.

## 2. Exact global functor

For inward facet normals \(r_1,\ldots,r_f\), define

\[
 F_{\mathrm{arr}}(P)=
 \bigl(\mathbf P(V),Q_q,\mathbf P(r_1^\perp),\ldots,
       \mathbf P(r_f^\perp),\sigma_0,C_P\bigr).                       \tag{2.1}
\]

> **Theorem 2.1 (exact AGCat carrier).** The assignment (2.1) is a faithful
> functor and admits an exact return functor on its image:
> \[
>                  R_{\mathrm{arr}}F_{\mathrm{arr}}(P)\cong P.        \tag{2.2}
> \]

Indeed, \(Q_q\) determines \(q\) up to common scale. Quadric polarity sends
\(H_i\) back to the line \(kr_i\); the real framing selects its inward ray.
Thus the carrier recovers all normalized pairings, finite angles, parallel
relations, ultraparallel distances, and the identical half-space
intersection. Passing to unmarked objects is the finite quotient by facet
permutations.

The regular coordinates on the diagonal-torus quotient of this arrangement
are the closed Cartan products

\[
 b_C=a_{i_1i_2}a_{i_2i_3}\cdots a_{i_si_1},\qquad
 a_{ij}=\frac{2q(r_i,r_j)}{q(r_i,r_i)}.                              \tag{2.3}
\]

They are invariant under independent root rescaling. If \(T_i=1-s_i\), then

\[
 b_C=\operatorname {tr}(T_{i_1}\cdots T_{i_s})\in\mathcal O_k,       \tag{2.4}
\]

also for nonfree projective \(\mathcal O_k\)-modules. Hence Vinberg
arithmeticity is literally an integral-point condition on the exact
algebraic arrangement carrier, together with one Lorentzian and all other
definite real-place conditions.

## 3. AG-arrangement pushback: finitely many fields

> **Theorem 3.1 (global finite field set).** In hyperbolic dimensions
> \(n\ge2\), the minimal Vinberg fields of finite-volume arithmetic
> hyperbolic Coxeter chambers form a finite set
> \(\mathscr K_{\mathrm{Cox}}\).

This theorem is obtained without reflection maximality and without a global
reflective-lattice or reflection-group finiteness theorem.

For dimensions \(n\ge6\), a minimal-visual-aperture facet and the
face-average/Hodge argument produce a full-field arrangement packet of
minimality \(14\). Its integral point lies in one of finitely many
Lanner/triangle/quadrilateral or transition charts. Each transition chart
contains a strict-capacity tail or fork determinant coordinate. Definite
nongeometric conjugates and relative Fekete theory make its integral locus
finite. This is the audited
[global minimality-14 field funnel](global_minimality14_integral_arrangement_field_funnel.md).

In dimensions four and five, exact codimension-one/two descent reduces the
field to dimension three or two, or to one of those finite transition
packets. Dimension two is closed by the exact real-DM orbifold-curve carrier
and the Long--Maclachlan--Reid genus-zero arithmetic theorem.

The formerly open compact-three-dimensional case is closed by the audited
[general-fork theorem](compact_h3_general_fork_arrangement_field_finiteness.md).
For a transverse three-wall packet, put

\[
 \beta=t_{12}+t_{13}+t_{23}-z-4,\quad
 t_{ij}=a_{ij}a_{ji},\quad z=a_{12}a_{23}a_{31}.                       \tag{3.1}
\]

This is an integral rescaling-invariant regular function on the arrangement
quotient. The selected packet has signature \((1,2)\), while every
nonselected packet is definite, so

\[
                 0<\beta_0<256,\qquad \beta^\tau<0,qquad
                 k=\mathbf Q(\beta).                                  \tag{3.2}
\]

For Coxeter labels \(m,n\), let
\(s_d=4\sin^2(\pi/d)\) and \(E=\mathbf Q(s_m,s_n)\). The conjugates obey

\[
 -\frac{s_m^\rho s_n^\rho}{4}\le\beta^\tau<0.                       \tag{3.3}
\]

The product of interval capacities is strictly below one. The corrected
norm estimate, including the exceptional selected place, gives a uniform
label cutoff \(m,n<10^6\); fixed-label Fekete finiteness then finishes
compact \(\mathbf H^3\). Noncompact chambers have \(k=\mathbf Q\). Combining
these branches proves Theorem 3.1.

## 4. From finite fields to finite native hosts

The downstream theorem is independently reconstructed in
[the bounded-frame native-host report](ag_field_funnel_bounded_frame_native_host_closure.md).

> **Theorem 4.1 (fixed field and rank).** Fix a totally real field \(k\) and
> rank \(m\ge3\). Up to quadratic similarity, only finitely many full native
> integral admissible arithmetic hosts can support a finite-covolume
> arithmetic reflection group in the Vinberg model; hence only finitely many
> reflection-maximal chambers occur for this \((k,m)\).

Here is the mechanism, included to make the noncircularity visible.

1. Arithmeticity gives a full \(W(P)\)-stable projective lattice \(L\) by
   taking a normal core and summing finitely many translates. Then
   \(W(P)\subseteq R(L)\subseteq O^+(L)\), all with finite index.
2. The geometric connected-frame lemma selects \(m\) wall reflections
   spanning the space, with connected nonorthogonality graph and
   \[
      c_{ij}=\frac{4q(\delta_i,\delta_j)^2}
                    {q(\delta_i,\delta_i)q(\delta_j,\delta_j)}<62^2.  \tag{4.1}
   \]
3. Length-two and longer simple cycle traces are algebraic integers. They
   are bounded by \(62^m\) at the selected place and by \(2^m\) elsewhere.
   Minkowski discreteness leaves finitely many cycle packets.
4. On a connected graph, the complete simple-cycle packet determines the
   Cartan matrix up to diagonal root rescaling. Thus only finitely many
   \(k\)-conjugacy classes of frame reflection representations occur.
5. Each representation is absolutely irreducible, so Burnside gives
   \(k[s_1,\ldots,s_m]=\operatorname {End}_k(V)\). For one stable lattice
   \(L_0\), the order
   \(\Lambda=\mathcal O_k[s_1,\ldots,s_m]\) is full in
   \(B=\operatorname {End}_{\mathcal O_k}(L_0)\). A nonzero conductor ideal
   \(\mathfrak c\) satisfies \(\mathfrak cB\subseteq\Lambda\).
6. For any \(\Lambda\)-stable full projective lattice \(L\), put \(M=BL\).
   Then
   \[
                         \mathfrak cM\subseteq L\subseteq M.          \tag{4.2}
   \]
   Matrix units show that a \(B\)-stable lattice is \(\mathfrak aL_0\) for
   a fractional ideal \(\mathfrak a\). The class group is finite, and for
   fixed \(M\), (4.2) has only finitely many submodules because
   \(M/\mathfrak cM\) is finite.

This proves the fixed-data theorem directly for projective modules. It is an
independent reconstruction and simplification of the arithmetic core behind
Nikulin's fixed-data argument; it does not invoke that theorem's finiteness
conclusion.

## 5. Dimension, maximality, and the global theorem

Vinberg's compact theorem and the Khovanskii--Prokhorov finite-volume theorem
give the purely geometric ceiling

\[
                              n\le995.                              \tag{5.1}
\]

This theorem applies to all finite-volume hyperbolic Coxeter polyhedra and
assumes neither arithmeticity nor any field or reflective-lattice
finiteness. Thus Theorem 3.1 and (5.1) leave only finitely many pairs

\[
              (k,m)\in\mathscr K_{\mathrm{Cox}}
                         \times\{3,\ldots,996\}.                       \tag{5.2}
\]

For a reflection-maximal input, the inclusion
\(W(P)\subseteq R(L)\) from Section 4 is equality. If maximality is only
among arithmetic reflection groups, equality still holds because \(R(L)\)
is a finite-index arithmetic reflection overgroup. Finally, the complete
locally finite mirror arrangement of \(R(L)\) has one chamber orbit: a
generic finite gallery is crossed by successive reflections.

> **Theorem 5.1 (AG-mediated global finiteness).** In hyperbolic dimensions
> \(n\ge2\), there are only finitely many unmarked normalized
> reflection-maximal finite-volume arithmetic hyperbolic Coxeter chambers,
> and hence only finitely many conjugacy classes of their side-reflection
> groups. The same conclusion holds for maximality among arithmetic
> reflection groups. Equivalently, this is the restriction of the exact
> functor on \(\mathsf{Cox}^{\mathrm{mark}}_{\mathrm{arith}}\) to its
> reflection-maximal full subgroupoid.

The complete logical chain is

\[
\boxed{
\begin{gathered}
P\longmapsto F_{\mathrm{arr}}(P)\in\mathsf{AGCat}^{\mathbf R}_{\mathrm{arr}},
\qquad R_{\mathrm{arr}}F_{\mathrm{arr}}(P)=P,\\
\text{integral arrangement coordinates + real-place positivity + Fekete/orbicurves}
\Longrightarrow |\mathscr K_{\mathrm{Cox}}|<\infty,\\
\text{geometric }n\le995+\text{ connected }62\text{-frame}
\Longrightarrow\text{ finitely many native hosts},\\
\text{reflection maximality}\Longrightarrow W(P)=R(L)
\Longrightarrow\text{ finitely many chambers and groups.}
\end{gathered}}                                                     \tag{5.3}
\]

## 6. Exact attribution and limitation

The result answers the functor-and-pushback question affirmatively, but the
attribution matters.

| ingredient | role |
|---|---|
| Vinberg structural reflection lemma | reduces every arithmetic reflection input to the admissible quadratic/lattice-form model; no finiteness conclusion |
| exact real-framed quadric arrangement | functorial AG encoding and exact return |
| arrangement cycle integrality, contracting determinant charts, Fekete | new global finite-field pushback |
| real-DM genus-zero orbifold curve + Long--Maclachlan--Reid | dimension-two and polygon branch |
| Vinberg; Khovanskii--Prokhorov | universal geometric dimension ceiling |
| Nikulin geometric appendix | connected bounded frame with \(A=62\) |
| cycle invariant theory + full-order conductor sandwich | reconstructed finite native hosts for fixed \((k,m)\) |
| reflection maximality | one full reflection group/chamber orbit per host |

The standard direct moduli implication

\[
 \text{reflection maximality}\Longrightarrow
 \text{uniformly bounded Hilbert/log degree of the facet divisor}    \tag{6.1}
\]

is still not known and is not claimed. BAB, KSBA, K-moduli, root stacks,
and automorphic compactifications do not presently prove (6.1). The theorem
instead uses AG where it has actual leverage--the integral arrangement
field funnel--and then uses transparent geometric and arithmetic bridges to
push the finite field output all the way back to the original chambers.

## 7. Frozen evidence to be included in the terminal audit

- [exact AGCat synthesis](exact_agcat_functor_and_ag_finiteness_pushback_synthesis.md);
- [dimension-six-plus field theorem](global_minimality14_integral_arrangement_field_funnel.md)
  and its [hostile audit](global_minimality14_integral_arrangement_field_funnel_hostile_audit.md);
- [low-dimensional descent](low_dimensional_minimality14_field_descent_to_h3.md);
- [compact \(\mathbf H^3\) theorem](compact_h3_general_fork_arrangement_field_finiteness.md)
  and its [independent hostile audit](compact_h3_general_fork_arrangement_field_finiteness_independent_hostile_audit.md);
- [finite-field-to-native-host closure](ag_field_funnel_bounded_frame_native_host_closure.md)
  and its [independent hostile audit](ag_field_funnel_bounded_frame_native_host_closure_independent_hostile_audit.md).
