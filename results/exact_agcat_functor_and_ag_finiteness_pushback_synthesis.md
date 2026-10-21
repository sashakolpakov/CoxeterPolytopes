# Exact AGCat functor and global finiteness pushback: terminal synthesis

**Status (2026-08-12): terminal PASS.**  The low-dimensional descent and the
assembled global composition have independent hostile PASS, and the full
relevant regression suite reports `31 passed`.
Every finite-volume arithmetic hyperbolic Coxeter chamber in dimension
`n>=2` in the
standard Vinberg quadratic/lattice model has an exact functorial carrier in
a genuine real-framed algebraic arrangement category.  Integral-point
arithmetic on that carrier gives a finite set of ground fields in all
dimensions.  A universal geometric dimension ceiling, a bounded connected
wall frame, a reconstructed fixed-data order argument, and reflection
maximality then give finitely many maximal chambers and reflection groups.

The theorem is **AG-mediated, not AG-only**.  The exact arrangement and its
integral coordinates are indispensable in the finite-field step.  The
dimension ceiling `n<=995` and frame constant `A=62` are classical
hyperbolic Coxeter geometry.  The fixed-field/fixed-rank host proof is
reconstructed from cycle invariant theory and orders in matrix algebras; it
does not cite the conclusion of a reflective-lattice finiteness theorem.

## 1. The answer in one diagram

\[
\boxed{
\begin{gathered}
 P\longmapsto
 F_{\rm arr}(P)=
 (\mathbf P(V),Q_q,H_1,\ldots,H_f,\sigma_0,C_P),
 \qquad R_{\rm arr}F_{\rm arr}(P)\cong P,\\[2mm]
 \text{integral Cartan cycles on the arrangement}
 +\text{ real-place definiteness}
 +\text{ Fekete/orbifold curves}
 \Longrightarrow |\mathscr K_{\rm Cox}|<\infty,\\[2mm]
 n\le995+\text{ connected }62\text{-frame}
 +\text{ finite full orders and projective modules}
 \Longrightarrow\text{ finitely many native hosts},\\[2mm]
 \text{reflection maximality}
 \Longrightarrow W(P)=R(L)
 \Longrightarrow\text{ finitely many chambers and groups.}
\end{gathered}}
\tag{1.1}
\]

Thus the requested functor exists, its return is exact, and algebraic-
arithmetic information extracted from its image is genuinely used to prove a
global theorem about the original maximal Coxeter chambers.

## 2. Source and target categories

Let \(\mathsf{Cox}^{\mathrm{mark}}_{\mathrm{arith}}\) be the groupoid
of finite-sided finite-volume arithmetic hyperbolic Coxeter chambers of
dimension `n>=2`, with minimal Vinberg field \(k\), distinguished real
embedding, ordered inward facet-root rays, and selected positive component.
Morphisms are marked, field-isomorphism-compatible projective similitudes.
The structural reflection lemma quoted as ABSW Lemma 5.2 from Vinberg Lemma 7
reduces arithmetic reflection inputs to this admissible quadratic/lattice
model; this is a model reduction, not a finiteness theorem.  The terminal
theorem restricts this source to its reflection-maximal full subgroupoid.
Maximality may mean maximal among all discrete hyperbolic reflection groups
or merely among arithmetic reflection groups; the two lead to the same
terminal argument.

Let \(\mathsf{AGCat}^{\mathbf R}_{\mathrm{arr}}\) be the groupoid of

\[
 (\mathbf P(V),Q,H_1,\ldots,H_f,\sigma_0,C),             \tag{2.1}
\]

defined over a displayed totally real field \(k\), where \(Q\) is a smooth
projective quadric, the \(H_i\) are ordered polar hyperplanes, \(\sigma_0\)
is a distinguished real place, and \(C\) is a chosen real component with its
oriented sides.  Morphisms are field-isomorphism-compatible projective
isomorphisms preserving \(Q\), each ordered \(H_i\), \(\sigma_0\), and
\(C\); the return functor on morphisms restricts them to the recovered
hyperbolic chambers.  This is an algebraic quadric--hyperplane arrangement
with finite real descent/framing data.  The framing is necessary: the bare
complex arrangement does not remember which half-spaces define the chamber.

For inward facet normals \(r_1,\ldots,r_f\), define

\[
 F_{\rm arr}(P)=
 \bigl(\mathbf P(V),Q_q,\mathbf P(r_1^\perp),\ldots,
       \mathbf P(r_f^\perp),\sigma_0,C_P\bigr).           \tag{2.2}
\]

> **Theorem 2.1 (exact global AGCat carrier).**  The assignment (2.2) is a
> faithful functor and has a return functor on its image satisfying
> \[
>                         R_{\rm arr}F_{\rm arr}(P)\cong P.             \tag{2.3}
> \]

A smooth projective quadric determines its quadratic form up to common
scale.  Quadric polarity returns each normal line \(kr_i\), and the real
framing selects its inward ray.  The normalized pairings then return every
finite angle, parallel relation, ultraparallel distance, and the identical
half-space intersection.  Passing to unmarked objects is the finite quotient
by facet permutations.

## 3. Arithmeticity is an integral-point condition on the carrier

Put

\[
 a_{ij}=\frac{2q(r_i,r_j)}{q(r_i,r_i)},\qquad
 b_C=a_{i_1i_2}a_{i_2i_3}\cdots a_{i_si_1}.              \tag{3.1}
\]

Independent rescaling \(r_i\mapsto\lambda_i r_i\) changes
\(a_{ij}\mapsto\lambda_j\lambda_i^{-1}a_{ij}\), while every closed product
\(b_C\) is invariant.  These are regular invariant coordinates on the
diagonal-torus quotient of the marked arrangement chart.

If \(s_i\) is a facet reflection and \(T_i=1-s_i\), projective-module trace
gives

\[
 b_C=\operatorname {tr}(T_{i_1}\cdots T_{i_s})\in\mathcal O_k.        \tag{3.2}
\]

This remains true for a nonfree finite projective
\(\mathcal O_k\)-module; no individual directed Cartan entry is declared
integral.  Vinberg arithmeticity is therefore literally an integral-point
condition in the exact arrangement coordinates, together with one
Lorentzian and all other definite real-place conditions.

## 4. The all-dimensional finite-field pushback

> **Theorem 4.1 (global finite Vinberg field set).**  The minimal Vinberg
> fields of finite-volume arithmetic hyperbolic Coxeter chambers form a
> finite set \(\mathscr K_{\mathrm{Cox}}\).  Reflection maximality is not
> needed.

The proof is assembled from four exact branches.

### 4.1 Dimensions at least six

A facet of minimal visual aperture has closed-star normalized pairings below
14.  The face-average/Dehn--Sommerville--Hodge argument and exact four-node
support reconstruction force either a rational, Lanner, triangle, or
quadrilateral branch, or one of finitely many transition arrangement charts.
On each transition chart an integral determinant coordinate has contracting
intervals at every nonselected embedding.  Relative Fekete theory makes the
integral locus finite.  This is the audited
[global minimality-14 field funnel](global_minimality14_integral_arrangement_field_funnel.md).

### 4.2 Dimensions four and five

Exact codimension-one/two descent sends the full field to a compact
three-dimensional or two-dimensional chamber unless it already lies in a
finite rational/Lanner/transition packet.  The descent retains literal
normal restrictions only in the orthogonal branches and keeps transverse
information in the four-wall arrangement charts otherwise.  See the
[low-dimensional descent theorem](low_dimensional_minimality14_field_descent_to_h3.md)
and its terminal independent audit listed in Section 9.

### 4.3 Compact dimension three

The previously missing compact \(\mathbf H^3\) branch is closed by a
three-wall general-fork chart.  With walls 1 and 2 ultraparallel and wall 3
meeting them at labels \(m,n\), define the invariant integral coordinate

\[
 \beta=t_{12}+t_{13}+t_{23}-z-4,\qquad
 t_{ij}=a_{ij}a_{ji},\quad z=a_{12}a_{23}a_{31}.          \tag{4.1}
\]

The selected packet has signature \((1,2)\), whereas every nonselected
packet is definite.  Hence

\[
                 0<\beta_0<256,\qquad \beta^\tau<0,\qquad
                 k=\mathbf Q(\beta).                    \tag{4.2}
\]

For \(s_d=4\sin^2(\pi/d)\) and
\(E=\mathbf Q(s_m,s_n)\), every nonselected extension satisfies

\[
 -\frac{s_m^\rho s_n^\rho}{4}\le\beta^\tau<0.          \tag{4.3}
\]

The interval-capacity product is below one.  The corrected norm estimate,
including the unique exceptional selected embedding, gives

\[
 1\le |N_{k/\mathbf Q}(\beta)|
 <\frac{1024}{s_{m,0}s_{n,0}}
   \left(\frac{|N_E(s_ms_n)|}{4^{[E:\mathbf Q]}}\right)^{[k:E]},     \tag{4.4}
\]

which forces \(m,n<10^6\).  Fixed-label Fekete finiteness then closes the
fork locus.  The all-perpendicular branch is an exact arithmetic Coxeter
polygon; the remaining incidence signs give the already finite tail and
small packets.  See the
[compact \(\mathbf H^3\) theorem](compact_h3_general_fork_arrangement_field_finiteness.md)
and its
[independent hostile audit](compact_h3_general_fork_arrangement_field_finiteness_independent_hostile_audit.md).

### 4.4 Dimension two and noncompact chambers

Doubling an arithmetic Coxeter polygon gives an exact real algebraic
genus-zero orbifold log curve.  Passing to a maximal arithmetic Fuchsian
over-lattice preserves the adjoint/Vinberg field and genus zero.
Long--Maclachlan--Reid finiteness then leaves finitely many fields.  This is
an AG-mediated arithmetic--spectral input, not a consequence of
\(\mathcal M_{0,r}\) alone.

For a noncompact chamber, a \(k\)-rational cusp vector remains isotropic at
every real embedding, contradicting definiteness at a nonselected place.
Thus its field is \(\mathbf Q\).  Together, Sections 4.1--4.4 prove Theorem
4.1.

## 5. Finite fields imply finitely many native hosts

The downstream proof is reconstructed in the
[bounded-frame native-host theorem](ag_field_funnel_bounded_frame_native_host_closure.md)
and independently audited.  Fix a totally real field \(k\) and
\(m=n+1\ge3\).

1. A normal-core/finite-translate sum turns the arithmetic chamber group
   \(W(P)\) into a group preserving one full projective lattice \(L\).  If
   \(R(L)\) is generated by every \(q\)-reflection preserving \(L\), then
   \[
                         W(P)\subseteq R(L)\subseteq O^+(q,L)          \tag{5.1}
   \]
   with finite index.
2. The imported geometric frame lemma selects \(m\) wall normals spanning
   \(V\), with connected nonorthogonality graph and
   \[
     \frac{4q(\delta_i,\delta_j)^2}
          {q(\delta_i,\delta_i)q(\delta_j,\delta_j)}<62^2.             \tag{5.2}
   \]
3. Projective traces make every simple directed cycle product integral.
   It is bounded by \(62^m\) at the selected place and by \(2^m\) at every
   other place.  Minkowski discreteness leaves finitely many complete cycle
   packets.
4. On a connected graph, the complete simple-cycle packet determines the
   Cartan matrix up to diagonal root rescaling.  Hence only finitely many
   \(k\)-conjugacy classes of frame reflection representations occur.
5. Each representation is absolutely irreducible.  Burnside gives
   \(k[s_1,\ldots,s_m]=\operatorname {End}_k(V)\), and the invariant form is
   unique up to \(k^\times\)-scale.
6. For one stable projective lattice \(L_0\), put
   \[
   \Lambda=\mathcal O_k[s_1,\ldots,s_m]
   \subset B=\operatorname {End}_{\mathcal O_k}(L_0).
   \]
   This is a full order.  A nonzero conductor ideal \(\mathfrak c\) has
   \(\mathfrak cB\subseteq\Lambda\).  For any \(\Lambda\)-stable full
   projective lattice \(L\), setting \(M=BL\) gives
   \[
                            \mathfrak cM\subseteq L\subseteq M.       \tag{5.3}
   \]
   Matrix units force \(M=\mathfrak aL_0\) for a fractional ideal
   \(\mathfrak a\).  The ideal class group is finite, and for fixed \(M\)
   the finite module \(M/\mathfrak cM\) has only finitely many submodules.

Thus only finitely many full native integral admissible arithmetic reflection
hosts occur for fixed \((k,m)\).  This handles nonfree modules, ideal classes,
units, and denominators directly.  It replaces rather than invokes Nikulin's
fixed-field/fixed-rank reflective-lattice finiteness conclusion.

## 6. Closing the dimension and maximality quantifiers

Vinberg's compact theorem and the Khovanskii--Prokhorov finite-volume
noncompact theorem give the universal geometric ceiling

\[
                              n\le995.                                \tag{6.1}
\]

This is a theorem about hyperbolic Coxeter polyhedra, independent of
arithmetic fields and reflective-lattice finiteness.  Theorem 4.1 therefore
leaves only finitely many pairs

\[
       (k,m)\in\mathscr K_{\mathrm{Cox}}\times\{3,\ldots,996\}.       \tag{6.2}
\]

For a reflection-maximal chamber, (5.1) is equality:

\[
                              W(P)=R(L).                              \tag{6.3}
\]

Maximality merely among arithmetic reflection groups is enough because
\(R(L)\) is an arithmetic finite-index reflection overgroup.  Conversely,
any discrete lattice overgroup of an arithmetic lattice is finite-index and
therefore arithmetic.

The complete locally finite mirror arrangement of \(R(L)\) has one chamber
orbit: a generic compact path between chamber interiors crosses finitely many
mirrors, and successive reflections give a finite gallery.  Each retained
native host therefore contributes at most one unmarked normalized chamber.

> **Theorem 6.1 (global AG-mediated finiteness).**  In hyperbolic dimensions
> `n>=2`, there are only finitely many unmarked normalized
> reflection-maximal finite-volume arithmetic
> hyperbolic Coxeter chambers, and hence finitely many conjugacy classes of
> their side-reflection groups.  The same conclusion holds for maximality
> among arithmetic reflection groups.  Equivalently, this is the restriction
> of the exact functor on
> \(\mathsf{Cox}^{\mathrm{mark}}_{\mathrm{arith}}\) to its
> reflection-maximal full subgroupoid.

An existential dimension-free facet bound follows immediately by taking the
maximum facet count over this finite set.  No facet bound is assumed in the
proof.

## 7. Exact attribution: what AG did and did not do

| ingredient | exact role |
|---|---|
| exact real-framed quadric arrangement | functorial algebraic carrier and identical return |
| integral Cartan-cycle functions and determinant charts | strict-capacity integral-point field funnel |
| real-DM genus-zero orbifold curve | dimension-two/polygon field branch |
| relative Fekete and Long--Maclachlan--Reid | arithmetic finiteness on those AG charts |
| ABSW Lemma 5.2 citing Vinberg Lemma 7 | structural reduction of arithmetic reflection inputs to the admissible quadratic/lattice model; no finiteness conclusion |
| Vinberg and Khovanskii--Prokhorov | universal geometric dimension ceiling |
| Nikulin geometric appendix | connected bounded frame with \(A=62\) |
| cycle invariant theory and conductor sandwich | reconstructed fixed-\((k,m)\) native-host finiteness |
| reflection maximality | equality with the full host and one chamber orbit |

The theorem is not attributed to BAB, KSBA, K-moduli, or a weak blow-up
encoding.  Nor is it attributed to the direct implication

\[
 \text{reflection maximality}\Longrightarrow
 \text{uniformly bounded Hilbert/log degree of the facet divisor},          \tag{7.1}
\]

which remains open and is called **AG-DEG**.  Likewise, the former
**MAX-ARR** proposal--a finite-type bounded arrangement envelope with
zero-dimensional exact integral Gram image--remains an interesting stronger
direct-moduli theorem.  It is no longer needed for global Coxeter finiteness,
because the field funnel plus geometric frame/order closure bypasses it.

The proof is therefore accurately described as

\[
 \text{AG arrangement/arithmetic geometry}
 +\text{ hyperbolic Coxeter geometry}
 +\text{ arithmetic of orders}.                           \tag{7.2}
\]

Calling it AG-only would be false; saying AG provides no global pushback is
also now false.

## 8. Relation to the surface and K3 inverse-carrier programme

The global arrangement carrier is not a claim that every chamber is the nef
cone of a surface.  The earlier three-branch direct-rational/weak-log/K3
Property-4 disjunction is still false, and the Borcherds \(\mathbf H^{21}\)
control still requires mixed characteristic-two curve/action data for a
surface interpretation.  Those surface-specific universality questions
remain open exactly as recorded in the
[canonical four-property audit](canonical_chamber_four_property_completion_audit.md)
and its [global AGCat scope addendum](canonical_chamber_four_property_completion_global_agcat_addendum.md).

What changes is the global finiteness conclusion: an exhaustive surface
compiler is no longer needed to prove that maximal arithmetic Coxeter
chambers are finite.  The exact quadric arrangement already supplies a
universal algebraic carrier, and its integral-point field pushback is enough
after the transparent geometric and order-theoretic bridges above.

The earlier scoped results remain useful calibrations:

- in \(\mathbf H^2\), the real-DM orbifold curve itself gives a direct
  quotient-curve proof;
- on the Lambert stable-pair curve, an explicit contracting integral
  coordinate makes the entire arithmetic locus finite; and
- K3/Mukai and rational-surface carriers give sharper bounds on their own
  intrinsic subclasses.

They are no longer the limit of the all-dimensional theorem.

## 9. Frozen evidence

| artifact | SHA-256 |
|---|---|
| [global composite theorem](agcat_global_maximal_arithmetic_coxeter_finiteness.md) | `441deb5a379fc59c6b4c1e1e86212f5f84db4e603a4b37259200db43a4508dfd` |
| [independent global composition audit](agcat_global_maximal_arithmetic_coxeter_finiteness_independent_hostile_audit.md) | `6fbf3f5954d66390e82916b3fd68720ba58ed106c681e7ff7196f8e30667a44e` |
| [dimension-two orbifold-curve theorem](arithmetic_coxeter_polygon_orbifold_curve_ag_finiteness.md) | `68ab8d85002106fde4493540ee6c6409247211cdda5364a3181000ae0135b659` |
| [dimension-two hostile audit](arithmetic_coxeter_polygon_orbifold_curve_ag_finiteness_hostile_audit.md) | `7be0ea8da1eb6dcc97ced1b68f40a47c9eb3ae382b47ab054e77735df4afcfc9` |
| [dimension-six-plus field theorem](global_minimality14_integral_arrangement_field_funnel.md) | `05acee48bd7a6c03e92febc32434cc38bcd3581c0e9044dafec3f07f9dbab678` |
| [dimension-six-plus hostile audit](global_minimality14_integral_arrangement_field_funnel_hostile_audit.md) | `e5e2754bf886d922ba50173e8afef67ff7eceffc62b706bb3c4dea4a309cac69` |
| [low-dimensional descent](low_dimensional_minimality14_field_descent_to_h3.md) | `8157c0779705f5f265e75126d10361441bf2e2b504e522d917cbda848b75f019` |
| [low-dimensional independent audit](low_dimensional_minimality14_field_descent_to_h3_independent_hostile_audit.md) | `74c32b41b593419951ecb6db0c43c3364134cc033e77038175156bffb385b5da` |
| [compact \(\mathbf H^3\) general-fork theorem](compact_h3_general_fork_arrangement_field_finiteness.md) | `4a971c64460a7d234114ecb12b793e613b080945850098da37b903d2c50876d5` |
| [compact \(\mathbf H^3\) independent audit](compact_h3_general_fork_arrangement_field_finiteness_independent_hostile_audit.md) | `7bded0c2dc46c1652e06dd15ab691685736db4b327c6fb8a218ca84066cef996` |
| [finite-field-to-native-host closure](ag_field_funnel_bounded_frame_native_host_closure.md) | `7ba392dd833358379ea92d1385909a52b2dc65417a9cb4ccb856170118df02a3` |
| [native-host independent audit](ag_field_funnel_bounded_frame_native_host_closure_independent_hostile_audit.md) | `0b05127f62a89b16060d9cbaacde93009a56ac056e08bb7e16e0b8f79322c3b1` |
| [independent composition certifier](../scripts/certify_agcat_global_maximal_arithmetic_coxeter_finiteness_independent_audit.py) | `8acee124164a974157ea19a32b5a8d765ad2252f15078fd9882909089bf189c7` |
| [independent composition tests](../tests/test_agcat_global_maximal_arithmetic_coxeter_finiteness_independent_audit.py) | `b534548ed8816e271f9a623666d756f7505e3ecb17d59730ae0100fc02564785` |

The independent composition verdict is **PASS after four repaired
scope/interface defects**, with no remaining mathematical defect.  The full
owner and independent arithmetic regression suite reports `31 passed`.
