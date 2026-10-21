# Polytope-to-automorphic-divisor transfer: theorem, obstruction, and missing converse

**Status (2026-08-09):** strengthened theorem/conjecture split.  The
determinant-line construction, real-slice mirror identity, finite hidden-sheet
test, stacky recovery of finite Coxeter labels, and Bruinier's finite
Whittaker obstruction are proved.  On a direct determinant level, the
auxiliary-nonzero coefficient columns now provably span the whole cusp
obstruction dual.  Hence every finite target principal part has a finite
slice-disjoint complex Whittaker compensation.  More strongly, over a real
q-structure the **nonnegative cone** of those safe columns is the whole real
obstruction dual: a diagonal restriction and the
Knopp--Kohnen--Pribitkin sign-change theorem rule out a separating cusp form.
Thus an integral q-structure gives nonnegative integral compensation after
scaling, and hence an effective interior Borcherds divisor.  This arithmetic
upgrade is unconditional for `F=Q` by McGraw and, for general totally real
`F`, conditional on the Hilbert--Weil integral-basis adaptation explicitly
stated but not proved in detail by Bruinier.  A deep discriminant-kernel level
also separates the target facet sheets on the closed chamber.

The positive cone is now quantitatively bounded at a fixed level.  If
\(\ell\) is the trace discriminant-module level and \(K\) the diagonal
weight, then over \(F=\mathbb Q\) the safe columns with
\(\ell m\ll_\varepsilon(8K\ell^2)^{2+\varepsilon}\) already positively
span.  For higher-degree totally real \(F\), the same bound must be enlarged
by the finite safe-detection cutoff \(R_{\rm det}\); that cutoff is found by
exact coefficient-rank enumeration but has no proved formula in \(K,\ell\).
The sheet-separating integer is computable per polytope and its new
discriminant level divides \(N_{\rm sep}^2\ell_0\).

This is not automatic *bounded strongly reflective* correction.  Trace
height is now controlled per level (polynomially over `Q`).  At a finite
support, Proposition 5.7 controls conic denominators, Borcherds weight, and
global divisor degree from the integral obstruction-column height and the
arithmetic volume.  Those latter quantities are not yet uniform over the
compact polytope class; for `F=Q`, poles at compactification cusps are also
not excluded.
On a neat compact exact-support cover, facets do inject into distinct divisor
components and hence satisfy the qualitative inequality
`facets(P)*q_min <= w*lambda^n`.  The missing input is a uniform upper bound
for its right side, not another sheet-count lemma.  Consequently this report
still asserts neither a new hyperbolic Coxeter polytope nor a new dimension or
facet bound.

## 1. Input and conventions

Let $F$ be a totally real field with ring of integers $\mathcal O_F$, and fix
a real embedding $\sigma _0:F\hookrightarrow\mathbb R$.  A **marked arithmetic
Coxeter realization** is an even integral $\mathcal O_F$-lattice $(L,b)$ of
rank $n+1$ such that

\[
 \operatorname{sig}(L_{\sigma _0})=(1,n),\qquad
 L_\sigma\text{ is negative definite for }\sigma\ne\sigma _0,
 \tag{1.1}
\]

together with primitive negative vectors $r_1,\ldots,r_s\in L$, oriented
halfspaces, and a compact Coxeter polytope

\[
 P=\{[x]\in\mathbb H(L_{\sigma _0}):(x,r_i)\geq0\text{ for all }i\}.
 \tag{1.2}
\]

Each reflection $s_{r_i}$ is required to preserve $L$.  Replacing an initial
integral realization by a positive rational multiple makes it even without
changing its hyperplanes or reflections.  All claims below are canonical
relative to this marked even realization.  They are not claims that an
unmarked real polytope canonically determines an integral lattice.

We use signature $(2,n)$ for the type-IV quadratic space.  To apply
Bruinier's convention, replace the resulting form by its negative; it then
has signature $(n,2)$ at $\sigma _0$ and is positive definite at every other
real place.

## 2. The canonical determinant-line extension

### Theorem 2.1 (determinant-line type-IV extension)

Let $(L,b)$ satisfy (1.1), and put $K_L=\det_{\mathcal O_F}L$.  Equip this
rank-one projective module with the integral bilinear form

\[
 b_K(v_1\wedge\cdots\wedge v_{n+1},
     w_1\wedge\cdots\wedge w_{n+1})
 =2(-1)^n\det\bigl(b(v_i,w_j)\bigr).
 \tag{2.1}
\]

Then

\[
 \widehat L_0=L\perp K_L
 \tag{2.2}
\]

is an even integral lattice of signature $(2,n)$ at $\sigma _0$ and is
negative definite at every other real embedding.  The Coxeter group $W(P)$
acts isometrically on $\widehat L_0$ by its given action on $L$ and the
trivial action on $K_L$.

The even $W(P)$-stable overlattices $\widehat L\supseteq\widehat L_0$ are a
finite, effectively enumerable set: they correspond to the
$W(P)$-stable isotropic $\mathcal O_F$-submodules of the finite discriminant
quadratic module of $\widehat L_0$.

#### Proof

In a local basis, let $\delta=\det(b|L)$.  At $\sigma _0$ its sign is
$(-1)^n$, so $(-1)^n\delta>0$.  At any other real place, $L$ is negative
definite of rank $n+1$, so $(-1)^n\delta<0$.  Thus (2.1) adds one positive
direction at $\sigma _0$ and one negative direction everywhere else.  The
factor $2$ makes the rank-one summand even, and determinant functoriality
makes (2.1) independent of a choice of basis.  The trivial extension of the
$W(P)$-action is isometric.  The final assertion is the usual
overlattice--isotropic-submodule correspondence, with $W(P)$-stability imposed
on the finite discriminant module.  $\square$

### Remark 2.2 (what has and has not become finite)

Theorem 2.1 gives one distinguished line and finitely many gluings of that
line.  It does **not** make all sign-correct one-line extensions finite.
Indeed, if $a\in F^\times$ is positive at $\sigma _0$ and negative at all
other real places, then $ap$ has the same signs for every positive rational
prime $p$.  Infinitely many such $p$ give distinct classes in
$F^\times/F^{\times2}$, as seen from their odd valuations at unramified prime
ideals.  A search over all one-lines therefore requires an independent
height, determinant, or level bound.  The determinant line avoids that
infinite search; it does not classify it.

## 3. The algebraic mirror carrier and its exact real slice

Fix one $W(P)$-stable overlattice $\widehat L$ from Theorem 2.1, and put
$\widehat V=\widehat L\otimes F$.  Let $\mathcal D^+(\widehat V)$ be one
component of the type-IV domain of oriented positive two-planes in
$\widehat V_{\sigma _0}$.  Choose an arithmetic group
$\Gamma\subset O^+(\widehat L)$ containing the extended wall reflections.
Write

\[
 \mathscr X_\Gamma=[\mathcal D^+(\widehat V)/\Gamma]
 \tag{3.1}
\]

for the effective arithmetic Deligne--Mumford stack, and $X_\Gamma$ for its
coarse quotient.  A neat normal subgroup gives a finite smooth cover of the
same arrangement.

If $[F:\mathbb Q]>1$, then $\widehat V$ is anisotropic over $F$: a nonzero
$F$-isotropic vector would remain isotropic at a definite real place.
Consequently $X_\Gamma$ is projective and the stack is proper.  This is the
compact case in Bruinier's setup.  For $F=\mathbb Q$, the open quotient is in
general only quasi-projective and one must use a Baily--Borel or toroidal
compactification; no compactness assertion is made in that case.

Choose $0\ne e\in (K_L)_F$ and its sign so that $e^2>0$ at $\sigma _0$.
There is a totally real embedding

\[
 j_e:\mathbb H(L_{\sigma _0})\longrightarrow\mathcal D^+(\widehat V),
 \qquad [x]\longmapsto\operatorname{span}^{+}(e,x).
 \tag{3.2}
\]

For a negative vector $v\in\widehat V_{\sigma _0}$ let

\[
 \mathcal D_v=\{Z\in\mathcal D^+(\widehat V):Z\perp v\}.
 \tag{3.3}
\]

### Theorem 3.1 (mirror identity and elliptic-candidate recovery)

For every $v\in\widehat V_{\sigma _0}$,

\[
 j_e^{-1}(\mathcal D_v)
 =\{[x]:(v,e)=0\text{ and }(v,x)=0\}.
 \tag{3.4}
\]

In particular, for each facet root $r_i\in L$,

\[
 j_e^{-1}(\mathcal D_{r_i})=r_i^\perp\cap\mathbb H(L_{\sigma _0}),
 \tag{3.5}
\]

so the marked facet mirror is the real trace of a special divisor mirror.
For a linearly independent set $I$ of roots, the common mirrors meet the
hyperbolic real slice precisely when $\langle r_i:i\in I\rangle$ is negative
definite.  Thus common Heegner intersections recover the elliptic
**candidate** spans of the Coxeter system.

They do not, by themselves, recover the faces of $P$.  Such an intersection
is an actual face exactly when its oriented real-slice intersection also
satisfies every remaining inequality $(x,r_j)\geq0$ from (1.2).  Accordingly,
the face poset is obtained from the divisor-intersection candidates only
after the complete oriented-halfspace test.  This distinction prevents an
elliptic subdiagram lying outside the chamber from being counted as a face.

#### Proof

The plane $j_e(x)$ is perpendicular to $v$ exactly when both of its spanning
vectors are perpendicular to $v$, which proves (3.4) and (3.5).  In a
Lorentzian space, the orthogonal complement of an independent root span
contains a positive vector exactly when that root span is negative definite.
The final statement is simply the additional requirement that this positive
ray lie in all the chosen halfspaces defining $P$.  $\square$

### Theorem 3.2 (finite hidden-sheet gate)

Let $\mathscr T$ be a finite set of norm and discriminant-coset types of
special divisors.  Define

\[
 \mathcal R_{\mathscr T}(P)=
 \left\{
  \begin{array}{c|c}
  v\in\widehat L^\vee\cap e^\perp &
  \operatorname{type}(v)\in\mathscr T,\quad
  v^\perp\cap P\ne\varnothing
  \end{array}
 \right\}/\{\pm1\}.
 \tag{3.6}
\]

Then $\mathcal R_{\mathscr T}(P)$ is finite and can be enumerated by a
positive-definite lattice bound.  The exact comparison separates three
cases: facet sheets, nonfacet sheets touching only a proper face of
$\partial P$, and sheets cutting $P^\circ$.  An interior-cutting sheet
invalidates chamber completeness.  A boundary-touching sheet does not cut the
chamber, but it must be retained in the divisor-stratification and it
invalidates any stronger claim that the global divisor restricts to the facet
sheets alone.  Thus the same enumeration decides both the chamber gate and,
when requested, the stronger facet-only arrangement gate without conflating
them.

#### Proof

Choose a continuous positive-definite majorant of the quadratic form along
the compact set $P$.  If a vector of one of the finitely many fixed norms has
an orthogonal hyperplane meeting $P$, its majorant norm is uniformly bounded.
The dual lattice contains only finitely many vectors of bounded
positive-definite norm.  Enumeration in the finitely many specified cosets,
followed by the exact tests in (3.6), therefore terminates.  $\square$

The compactness of $P$ is essential in this formulation.  For a finite-volume
noncompact chamber one would need separate certified cusp-neighborhood bounds.

## 4. The stack remembers Coxeter labels

A neat quotient is useful for intersection theory, but it erases the torsion
that carries the angle label.  The arithmetic stack retains it.

### Proposition 4.1 (dihedral inertia at a two-wall stratum)

Suppose the marked walls $r_i^\perp$ and $r_j^\perp$ meet at Coxeter angle
$\pi/m_{ij}$.  At a generic point of the corresponding codimension-two
mirror intersection in $\mathscr X_\Gamma$, the two distinguished reflection
elements generate an inertia subgroup

\[
 \langle s_{r_i},s_{r_j}\rangle\cong I_2(m_{ij}),
 \qquad |I_2(m_{ij})|=2m_{ij}.
 \tag{4.1}
\]

Thus the order of the subgroup generated by the two marked branch
reflections recovers $m_{ij}$, not merely the fact that the two divisors
intersect.  Extra arithmetic stabilizers may enlarge total inertia, so it is
the distinguished two-reflection subgroup, rather than the unmarked order of
the full stabilizer, that is the invariant.

This stacky statement concerns intersecting walls and their finite labels.
Parallel and ultraparallel mirrors have no codimension-two crossing; their
distinction, and the distance attached to a dotted edge, is recovered from
the marked Lorentzian Gram on the real slice rather than from local inertia.

#### Proof

Both reflections fix the mirror intersection.  On its two-dimensional normal
slice they are the two real reflections through lines meeting at angle
$\pi/m_{ij}$; their product is a rotation of order $m_{ij}$.  Their faithful
normal action is the dihedral group of order $2m_{ij}$.  $\square$

More generally, a spherical collection of marked facet reflections gives the
corresponding finite Coxeter subgroup in the inertia of its mirror stratum.
On a stacky complete-intersection surface, the same local group survives at
the induced intersection point.  If this dihedral group is the full generic
stabilizer, orbifold integration weights the zero-dimensional stratum by
$1/(2m_{ij})$; with extra stabilizer it uses $1/|I_z|$, while the distinguished
dihedral subgroup still records $m_{ij}$.  This is the entry point for
label-sensitive orbifold Chern calculations.

## 5. The exact finite Bruinier obstruction

Put $M=-\widehat L$ in Bruinier's sign convention.  Then $M$ has signature
$(n,2)$ at $\sigma _0$ and is positive definite at all other real places.
Let $A_M$ be its finite discriminant module and $\rho_M$ its Weil
representation.  Bruinier's harmonic Whittaker input has mixed weight

\[
 k=\left(\frac{2-n}{2},\frac{n+2}{2},\ldots,
                   \frac{n+2}{2}\right),
 \tag{5.1}
\]

and the obstruction cusp forms have parallel weight

\[
 \kappa=\left(\frac{n+2}{2},\ldots,\frac{n+2}{2}\right).
 \tag{5.2}
\]

Let $S$ be a finite set of allowed pairs $(m,\mu)$ indexing special divisors,
and let $c=(c(m,\mu))_{(m,\mu)\in S}$.  If

\[
 g=\sum_{\mu\in A_M}\sum_m b_g(m,\mu)
       q^m\mathfrak e_\mu
 \in S_{\kappa,\rho_M},
 \tag{5.3}
\]

define

\[
 \operatorname{ob}_S(c)(g)
   =\sum_{(m,\mu)\in S}c(m,\mu)b_g(m,\mu).
 \tag{5.4}
\]

### Theorem 5.1 (finite necessary-and-sufficient Whittaker gate)

The prescribed finite principal part $c$ is the principal part of a weakly
holomorphic Whittaker form if and only if

\[
 \operatorname{ob}_S(c)=0
 \quad\text{in }S_{\kappa,\rho_M}^{\vee}.
 \tag{5.5}
\]

Equivalently, after choosing a basis of the finite-dimensional cusp space,
(5.5) is a finite matrix-kernel calculation.  If the coefficients $c(m,\mu)$
are integral, Bruinier's lift gives a meromorphic orthogonal modular form
with finite-order multiplier and divisor

\[
 \operatorname{div}(\Psi_c)
   =\sum_{(m,\mu)\in S}c(m,\mu)Z(m,\mu).
 \tag{5.6}
\]

#### Source and proof

Bruinier proves the exact sequence

\[
 0\longrightarrow M_k^!
 \longrightarrow H_k\xrightarrow{\,\xi_k\,}
 S_{\kappa,\rho_M}\longrightarrow0
 \tag{5.7}
\]

and the coefficient-pairing criterion (5.5) in Proposition 4.5 of
*Regularized theta lifts for orthogonal groups over totally real fields*.
The divisor statement is Theorem 6.8 (and Theorem 1.3); finite order of the
multiplier in the low-dimensional cases is completed in Corollary 7.4.

### Scope warning 5.2

Theorem 5.1 is necessary and sufficient for a **weakly holomorphic
Whittaker/Borcherds carrier with the prescribed principal part**.  It is not
known to be necessary for an arbitrary meromorphic automorphic form with the
same special-divisor support in the compact totally-real case.  Such a
converse would require, for example, injectivity of the relevant
Kudla--Millson lift.

Metzler's current converse theorem proves this implication for rational
orthogonal lattices of signature $(m^+,2)$, $m^+>3$, that split a hyperbolic
plane.  The hyperbolic-split hypothesis is essential there.  The present
totally-real compact spaces are anisotropic, hence do not split a hyperbolic
plane over $F$, so that theorem does not close this gap.

Nor can one assert that hyperbolic reflectivity automatically supplies an
automorphic correction.  There are two distinct statements here.  Wang's
Theorem 7.2 gives eighteen hyperbolically $2$-reflective lattices with no
automorphic correction in his sense, but that sense uses the rank-two
extension $U(m)\oplus S$; it is not by itself a counterexample to the
rank-one determinant extension used in this report.

A direct determinant-line counterexample follows instead by combining that
hyperbolic list with Wang's later complete classification of $2$-reflective
lattices of signature $(n,2)$, $n\geq8$.  Take

\[
 S=U\oplus E_8(-1)\oplus D_4(-1),\qquad
 \widehat S=S\oplus\langle8\rangle.
 \tag{5.8}
\]

The complete $(-2)$-reflection chamber of $S$ has finite volume, whereas
$\widehat S$ has signature $(2,13)$ and no lattice of signature $(2,13)$ is
$2$-reflective.  Hence no holomorphic modular form on the determinant
extension has zero divisor supported on its $(-2)$ mirrors while carrying
the facet mirrors.  Thus exact determinant-line reflectivity is genuinely
nonautomatic in the rational finite-volume class.  The example is
noncompact, so it does not settle the compact totally-real branch.  The full
scope audit, including the positive $II_{1,9}$ control and the distinction
between exact and compensated support, is in
[`borcherds_compact_inverse_leverage_audit.md`](borcherds_compact_inverse_leverage_audit.md).

### Theorem 5.3 (automatic finite slice-disjoint compensation)

Use the **unglued direct** determinant extension

\[
                 \widehat L_0=L\perp K_L.
\tag{5.9}
\]

One may replace it by the finite-index direct level
`N*widehat L_0=N*L perpendicular N*K_L`, `N>=2`, without changing the
rational type-IV domain or the canonical real slice.  In Bruinier's sign
convention let

\[
 A_N=A_{N L}\oplus A_{N K_L},\qquad
 V_N=S_{\kappa,\rho_{A_N}}.
\tag{5.10}
\]

Call a Whittaker index `(m,(alpha,beta))` **safe** when
`beta!=0` in `A_(N K_L)`.  Then:

1. every special-divisor vector of a safe type has nonzero auxiliary-line
   component, so its divisor is disjoint from the entire canonical real
   slice, not merely from `P`;
2. the Fourier-coefficient functionals at all safe indices span `V_N^*`;
3. consequently, for every finite target principal part `c_T` supported in
   auxiliary-zero components, there is a finite safe correction `c_D` such
   that `c_T+c_D` annihilates every cusp form.  At most `dim V_N`
   compensator columns are needed over `C`.

The spanning proof is short and exact.  If a cusp form is annihilated by all
safe coefficient evaluations, every one of its components with `beta!=0`
has identically zero Fourier expansion.  Its values therefore lie in

\[
 U=\mathbb C[A_{N L}]\otimes\mathbb C\mathfrak e_0.
\tag{5.11}
\]

The finite Weil representation factors over the orthogonal direct sum, while

\[
 \rho_{A_{N K_L}}(S)\mathfrak e_0
   =c\sum_{\beta\in A_{N K_L}}\mathfrak e_\beta,
 \qquad c\ne0.
\tag{5.12}
\]

The auxiliary discriminant group is nontrivial because
`K_L/N*K_L` injects into `A_(N K_L)`.  Hence
`U intersect rho(S)U=0`.  Modularity under `S` forces the cusp form to be
zero.  The annihilator of the safe columns is thus zero, and finite
dimensionality gives assertion 2 and the finite correction in assertion 3.
The central `mu <-> -mu` symmetry and the finite unit action preserve
`beta!=0`, so symmetrizing does not spoil safety.

There are two coefficient levels which must not be conflated.

- The argument above unconditionally gives a complex weakly holomorphic
  Whittaker combination and its regularized Green relation.
- In the proof of his Theorem 7.1, Bruinier explicitly states that, by
  adapting McGraw's argument, `V_N` has a basis of cusp forms with integral
  Fourier coefficients.  Relative to that published q-integrality input,
  all obstruction columns are integral.  A nonzero integral minor of safe
  columns gives a rational correction; clearing denominators gives an
  integral principal part with a positive integral multiple of `c_T` still
  nonzero.  Bruinier's Theorem 6.8 then produces a meromorphic modular form
  with finite-order multiplier and the prescribed target-plus-safe divisor.

The second bullet relies on Bruinier's stated but briefly sketched Hilbert
q-integrality adaptation; without it, complex linear algebra alone does not
license the integral Borcherds-product conclusion.

Safe *formal* indices are enough for the span theorem.  Some corresponding
special divisors may be empty.  Empty divisors are harmless for slice
separation, but finite dimensionality alone does not prove that the
coefficient columns belonging only to **nonempty** safe divisors span
`V_N^*`.

### Theorem 5.4 (positive safe cone and effective carrier)

Assume that `V_N` has a real q-structure: it has a basis with real Fourier
coefficients in the standard discriminant basis.  Let `V_(N,R)` be that real
form and let `ell_s in V_(N,R)^*` run through every safe Fourier-coefficient
column.  Then

\[
 \operatorname{cone}_{\mathbb R_{\geq0}}\{\ell_s:s\text{ safe}\}
                  =V_{N,\mathbb R}^{*}.
\tag{5.13}
\]

One finite set of safe columns already has full positive cone.  Consequently
every target principal part with nonnegative real coefficients has a finite
nonnegative real safe correction; at most `dim V_N` correction columns are
needed.  With a rational q-structure and rational target the correction can
be chosen rational.  With integral data, clearing denominators gives a
positive integral multiple of the target plus nonnegative integral safe
coefficients.

Here is the sign argument, formulated without assuming that the initial
infinite cone is closed.  Suppose a nonzero real cusp form `g` had every safe
Fourier coefficient nonnegative.  Theorem 5.3 implies that some component
`g_(alpha,beta)` with `beta!=0` is nonzero.  Pass to the finite-index kernel
of the finite Weil representation and restrict that scalar component to the
diagonal Hilbert modular curve:

\[
 h(z)=g_{\alpha,\beta}(z,\ldots,z)
     =\sum_{m\gg0}b_g(m,(\alpha,\beta))
          e(\operatorname{Tr}_{F/\mathbb Q}(m)z).
\tag{5.14}
\]

This is a scalar cusp form of positive real weight
`sum_j kappa_j`, with unitary multiplier, on a finite-index subgroup of
`SL_2(Z)`.  It is cuspidal at every subgroup cusp because each rational
diagonal cusp is a Hilbert cusp.  For fixed trace there are only finitely many
totally positive indices in the relevant shifted lattice, so every Fourier
coefficient of `h` is a finite sum of nonnegative real numbers.  A nonzero
safe coefficient makes `h` nonzero.  This contradicts the
Knopp--Kohnen--Pribitkin theorem, which forces infinitely many sign changes
for a nonzero cusp form of positive real weight with real coefficients and a
unitary multiplier on a finitely generated Fuchsian group of the first kind.
Fractional exponents are absorbed by the cusp width and multiplier.

For complete finiteness, the open sets

\[
 \{g\in S(V_{N,\mathbb R}):\ell_s(g)<0\},\qquad s\text{ safe},
\tag{5.15}
\]

cover the unit sphere by the preceding argument.  Choose a finite subcover.
The polar of the finitely generated, hence closed, cone of its columns is
zero, so that cone is the whole dual.  Conic Caratheodory gives the
`dim V_N` support bound.  For rational columns and target, the resulting
finite feasible polyhedron has a rational point.

McGraw proves the required integral q-structure for the classical
`F=Q` Weil representation.  For general totally real `F`, Bruinier states in
the proof of his Theorem 7.1 that McGraw's argument adapts to give an integral
basis, but does not supply that adaptation.  Thus (5.13) is unconditional
whenever a real q-structure is supplied, while its integral Borcherds upgrade
is unconditional for `F=Q` and conditional on Bruinier's stated input in the
general Hilbert--Weil setting.

Empty safe divisor types do not obstruct effectivity: their geometric terms
are simply zero.  Under the integral hypothesis, Bruinier's lift therefore
has an effective divisor on the type-IV domain, with the positive target
present.  If `[F:Q]>1`, the quotient is anisotropic and projective, so there
is no cusp boundary carrying an unrecorded pole; after killing the finite
multiplier by a power, this is a holomorphic automorphic section with a
nonzero effective divisor.  For `F=Q`, only interior effectivity is asserted:
the Borcherds product may still have a pole on a compactification boundary.

### Theorem 5.5 (orbit-square quantitative safe-height bound)

Let `ell` be the level of the finite quadratic module of the underlying even
trace lattice and put

```text
K = [F:Q]*(n+2)/2,
r(m) = ell*Tr_F/Q(m).
```

For a safe component `x`, let `O_x=(Z/ell Z)^times*x`.  The orbit-square

```text
H_O(z) = sum_(y in O_x) g_y(ell*z,...,ell*z)^2
```

is a scalar cusp form of integral weight `2K` on
`Gamma_0(4*ell^2)` with Dirichlet nebentypus.  Indeed, for
`gamma=(a b; c d)` in that group,

\[
 \delta=\operatorname{diag}(\ell,1)\gamma
        \operatorname{diag}(\ell,1)^{-1}
 =\begin{pmatrix}a&\ell b\\c/\ell&d\end{pmatrix}
 \in\Gamma_0(\ell),
\tag{5.16}
\]

and Strömberg's arXiv Lemma 5.13 (journal Lemma 5.6) gives

\[
 \rho(\delta)e_x=e(\ell bdQ(x))\chi_A(\delta)e_{dx}
                 =\chi_A(\delta)e_{dx}.
\tag{5.17}
\]

The unit orbit removes the permutation.  The square of `chi_A` is a
Dirichlet character: in even signature `chi_A` already is one, while in odd
signature it is a theta multiplier or its conjugate times one.  Diagonal
restriction is cuspidal on the finite Weil kernel because every rational
diagonal cusp is a Hilbert cusp; after rescaling, intersect the conjugated
kernel with `Gamma_0(4*ell^2)`.  This remains finite index, so vanishing at
its rational cusps plus modularity on the larger group gives cuspidality at
every cusp there.  A real q-structure gives real coefficients, and on the
imaginary axis a sum of squares is nonzero whenever one restricted component
is nonzero.

Chen--Lau--Wang therefore supplies

\[
 B_{\rm CLW}=\left\lceil
 C_\varepsilon(8K\ell^2)^{2+\varepsilon}\right\rceil.
\tag{5.18}
\]

For `F=Q`, the safe columns with `r(m)<=B_CLW` already have positive cone
`V_(N,R)^*`.  For general totally real `F`, define

\[
 R_{\rm det}=\min\{R:
 \operatorname{span}\{\ell_{m,x}:x\text{ safe},\ r(m)\le R\}
 =V_{N,\mathbb R}^*\}.
\tag{5.19}
\]

This cutoff is finite by Theorem 5.3.  The safe columns with

\[
 r(m)\le\max\{R_{\rm det},B_{\rm CLW}\}
\tag{5.20}
\]

have full positive cone.  The term `R_det` is essential to the present proof:
a nonzero Hilbert cusp form can vanish on the diagonal and equal-trace
coefficients can cancel.  Exact q-expansion rank enumeration finds it for a
fixed cusp space, but no bound in `K,ell` is proved.  In both field cases,
conic Caratheodory retains the `dim V_N` support bound and rational data give
a rational correction.

Thus the old scalar-descent gap is closed polynomially over `Q` and by a
finite per-space certificate in general.  By itself this does not bound the
conditioning of the resulting finite obstruction matrix.  Proposition 5.7
below shows that this one matrix height controls the conic denominator,
automorphic weight, and global divisor degree.  The full height proof is in
[`orbit_square_quantitative_safe_height.md`](orbit_square_quantitative_safe_height.md).

### Proposition 5.6 (explicit per-polytope target-sheet separation)

For the finitely many facet roots, a direct level `N=N_sep(P)` can be
computed so that the target cosets `+/-r_i mod N*widehat L_0` contain no
unwanted same-norm mirror meeting the closed chamber `P`.

Indeed, compactness and the positive-definite-majorant argument of Theorem
3.2 make the set of vectors in `widehat L_0` of a target norm whose mirror
meets `P` finite.  In a fixed integral basis, let `Delta(P)` be the finite set
of all nonzero differences between an unwanted vector and `+/-r_i`, together
with the differences needed to separate distinct desired facets modulo sign.
Then

```text
N_sep(P) = 1 + max{|d_j| : d in Delta(P), j a basis coordinate},
```

with `N_sep(P)=1` when `Delta(P)` is empty, keeps every forbidden difference
out of `N_sep(P)*widehat L_0`.  Since `r_i` belongs to
`(N_sep(P)*widehat L_0)^vee`, its special divisor is defined at that level.
Any vector in its target coset that meets `P` belongs to the finite set just
separated, and therefore defines the chosen facet mirror.  The `+/-` pair
defines the same divisor and is treated together.  If `ell_0` is the trace
discriminant-module level before scaling, then

```text
ell(N_sep(P)*widehat L_0) divides N_sep(P)^2*ell_0,
```

because `(N*M)^vee=N^(-1)*M^vee`.  This closes the level algorithm for each
fixed marked `P`; it is not a uniform bound on `N_sep(P)` over a class.  The
full coordinate audit is in
[`explicit_sheet_separating_level_bound.md`](explicit_sheet_separating_level_bound.md).

### Proposition 5.7 (weight--degree compression at bounded matrix height)

Choose an integral q-basis of the obstruction cusp space of rank `d`.  Let
`b_T` be the integral target obstruction column, let `S` be a finite safe set
whose positive cone contains `-b_T`, and put

```text
R   = max_(s in S) ||b_s||_infinity,
R_T = ||b_T||_infinity,
C_T = total target principal-part multiplicity.
```

There is a correction using `t<=d` linearly independent safe columns and
nonnegative integers `q,y_1,...,y_t` such that

\[
 qb_T+\sum_i y_i b_{s_i}=0,
\tag{5.21}
\]

with

\[
 q\le t^{t/2}R^t,
 \qquad
 y_i\le t^{t/2}R^{t-1}R_T,
 \qquad
 \sum_i y_i\le t^{t/2+1}R^{t-1}R_T.
\tag{5.22}
\]

Indeed, conic Caratheodory gives independent columns.  Select a nonsingular
`t`-row minor, apply Cramer's rule, clear its determinant, and apply
Hadamard's inequality.  Thus denominators are bounded minors once the
integral matrix height `R` is controlled; they are not a further independent
unknown.

If `E(s)=-B(s)=deg_lambda Z(s)/lambda^n` is Bruinier's normalized
Eisenstein coefficient and `E_H` is its maximum on the finite target and
safe support, the resulting Borcherds form satisfies

\[
 w\le E_Ht^{t/2}R^{t-1}(RC_T+tR_T),
\qquad
 w\lambda^n\le
 E_H\lambda^n t^{t/2}R^{t-1}(RC_T+tR_T).
\tag{5.23}
\]

These follow from Bruinier's exact identities

\[
 w=-\sum_s c_sB(s),
 \qquad
 w\lambda^n=\sum_sc_s\deg_\lambda Z(s).
\tag{5.24}
\]

For `F=Q`, unitarity in the vector Eisenstein Kloosterman expansion gives,
for `kappa=(n+2)/2>2`, the explicit estimate

\[
E(m,\mu)\le
 \frac{2(2\pi)^\kappa\zeta(\kappa-1)}{\Gamma(\kappa)}m^{\kappa-1}.
\tag{5.25}
\]

This rational closure is principally a noncompact/low-dimensional
calibration, not the compact high-dimensional answer: Meyer makes every
indefinite rational form of rank at least five isotropic, so a compact
arithmetic hyperbolic polytope over `Q` has dimension at most three.  The
compact `H8+` route therefore lies in degree greater than one, exactly where
Hilbert q-integrality, `R_det`, and `R` remain.

Moreover `lambda^n=2 vol_HM(Gamma)`, and the
Gritsenko--Hulek--Sankaran determinant/local-density formula computes this
volume and its subgroup-index scaling.  Over a general totally real field,
Fourier inversion gives a finite per-carrier Eisenstein majorant, while
Prasad's formula makes the arithmetic volume explicit; a class-uniform
bound still requires uniform Hilbert q-integrality, `R_det`, obstruction
height `R`, and field/local-parahoric data.

Bounded support alone cannot replace `R`: the integral columns
`(1,0),(-Q,-1),(0,1)` positively span `R^2`, but representing `(0,-1)`
forces coefficient at least `Q`.  This is an abstract conditioning no-go,
not a claimed family of Weil cusp spaces.  The complete proof and exact
scope are in
[`borcherds_weight_degree_quantitative_bound.md`](borcherds_weight_degree_quantitative_bound.md).

Combining Theorems 5.3--5.5 and Propositions 5.6--5.7 gives the following universal
carrier statement for the marked totally-real input class of Section 1:

```text
compact arithmetic P
  -> a deep direct determinant-level discriminant-kernel cover
  -> a finite slice-disjoint compensated weak-Whittaker relation
  ->, with the applicable integral q-structure, an effective interior
     Borcherds divisor
     whose target trace on closed P is exactly a positive multiple of its
     facet mirrors and whose compensating divisors miss the whole real slice.
```

For a general totally real field the last arrow has the conditional scope in
Theorem 5.4; in the anisotropic degree-greater-than-one case it has no boundary
cusp caveat once that q-structure is supplied.  Over `Q` it is unconditional
in the interior but may acquire boundary poles on compactification.

The deep discriminant-kernel cover need not contain the wall reflections.
Thus Proposition 4.1's inertia recovery stays on the original reflection
stack, while the deep-level carrier retains labels through the marked Gram;
no claim is made that the dihedral inertia itself descends to that cover.

The detailed proof and hostile coefficient audit are frozen in
[`determinant_line_slice_disjoint_obstruction_span.md`](determinant_line_slice_disjoint_obstruction_span.md).

## 6. The algebraic surface shadow and the canonical dividend

Assume for this section that $n\geq2$, that $X=X_\Gamma$ has been replaced by
a smooth projective neat cover, and let $\lambda$ denote the automorphic line
class.
Choose a sufficiently ample divisor $H$ and a general complete-intersection
surface

\[
 \Sigma=H_1\cap\cdots\cap H_{n-2},\qquad H_i\sim H,
 \tag{6.1}
\]

transverse to the selected special divisors.  Put $C_i=Z_i|_\Sigma$ and
$C=\sum_iC_i$.  This is an unconditional algebraic surface produced from the
polytope-derived orthogonal variety.  It is an arrangement carrier: it does
not assert that the $C_i$ form a Lorentzian Neron--Severi root basis or that
their halfspaces are the nef cone of $\Sigma$.

Adjunction gives

\[
 K_\Sigma=(K_X+(n-2)H)|_\Sigma,
 \tag{6.2}
\]

and

\[
 c_2(\Sigma)=
 \left(c_2(X)+(n-2)HK_X+\binom{n-1}{2}H^2\right)|_\Sigma.
 \tag{6.3}
\]

When $(\Sigma,C)$ is an SNC pair,

\[
 \bar c_1^2=(K_\Sigma+C)^2,
 \tag{6.4}
\]

\[
 \bar c_2=e(\Sigma\setminus C)
 =c_2(\Sigma)+K_\Sigma C+\sum_iC_i^2
                  +\sum_{i<j}C_iC_j.
 \tag{6.5}
\]

If the induced arrangement is not SNC, one applies the same framework on a
log resolution and includes the exceptional-divisor correction terms;
(6.4)--(6.5) are not asserted unchanged after such a resolution.

Under the standard log-nef/log-general-type hypotheses, logarithmic or
orbifold Bogomolov--Miyaoka--Yau gives

\[
 \bar c_1^2\leq3\bar c_2.
 \tag{6.6}
\]

Langer's orbifold form also permits quotient boundary coefficients and local
orbifold Euler terms, so Proposition 4.1 can in principle insert the labels
$m_{ij}$ into (6.6).

There is also a sharper canonical-class consequence when the divisor is the
entire reflection branch.  Let $B$ be the full branch divisor on a compact
coarse orthogonal quotient.  In rational Picard class,

\[
 K_X\equiv n\lambda-\frac12B.
 \tag{6.7}
\]

If a strongly reflective modular form of weight $w$ vanishes simply on all
reflection mirrors and on no other divisor, then

\[
 B\equiv2w\lambda,
 \qquad K_X\equiv(n-w)\lambda.
 \tag{6.8}
\]

In particular $w>n$ is the strongly reflective uniruledness criterion of
Gritsenko--Hulek; $w=n$ gives a numerically trivial canonical class in this
model, and $w<n$ lies on the positive-canonical side.

Equations (6.4)--(6.8) are genuine algebro-geometric invariants of the
polytope-derived divisor arrangement.  They do **not** presently give a
dimension bound, and log BMY by itself does not yet give a facet bound: the
quadratic terms $C^2$ and the pairwise intersections occur on both sides of
(6.6).  However, the exact-support level of Proposition 5.6 gives a direct
sheet-to-degree comparison.  Pass further to a finite neat cover of the
compact anisotropic quotient and let $\lambda=\mathcal M_1$.  Its Petersson
curvature is the invariant Kähler form, so $\lambda$ is an ample integral
line bundle.  If the effective carrier divisor $D$ contains every facet
target and $q_{\min}$ is the minimum multiplicity of a component carrying a
marked facet, exact support injects the facets into distinct irreducible
components of $D$: a component shared by facets $i\ne j$ would force
$Z_i\cap P$ to contain $F_j$, contrary to Proposition 5.6.  Therefore

\[
 \#\operatorname{Facets}(P)\,q_{\min}
 \leq D\lambda^{n-1}
 =w\lambda^n
 =w\operatorname{vol}_{\Omega}(X).
 \tag{6.9}
\]

Every nonzero irreducible divisor has positive integral $\lambda$-degree.
The equalities are Bruinier's weight--degree formula, after taking a power
to kill the finite multiplier.  Thus the qualitative sheet/degree bridge is
closed on the exact compact carrier.  It is false on an unseparated coarse
quotient, where several walls may remain sheets of one Heegner component.

The positive-cone theorem supplies the effective input required by (6.9),
but it does not close the quantitative problem.  Under the arithmetic
q-structure in Theorem 5.4 it supplies an effective interior divisor; in the
compact anisotropic case there is no cusp boundary either.  Proposition 5.6
computes the sheet-separating level for each fixed `P`, and Theorem 5.5 bounds
the largest scaled trace index polynomially over `Q`; for higher-degree
totally real fields it gives the exact finite cutoff
`max(R_det,B_CLW)`, with no proved formula for `R_det`.  Neither result is
uniform over a polytope class.  Proposition 5.7 does control conic
denominators, automorphic weight, and global Heegner degree once the finite
integral obstruction-column height `R` and arithmetic volume are bounded.
Passing to a
sheet-separating and neat level can itself increase discriminant,
$\lambda^n$, degree, and cusp-form dimension.  Hence (6.9) is not yet a
numerical wall bound intrinsic to the input class, and the effective divisor
cannot yet be inserted into (6.6)--(6.8) as a uniformly bounded log boundary.
What remains is a class-uniform level/matrix-height theorem (including a
Hilbert detection bound), field/local-volume control, and hence--by
Proposition 5.7--a noncircular bound for weight and top automorphic degree.
The full proof and coarse-quotient warning are in
[the sheet-to-degree report](facet_sheet_automorphic_degree_injection.md).

## 7. Exact theorem/conjecture boundary

For a **fixed canonical level and fixed prescribed support**, the maximal
directly decidable success locus of the original finite calculation is the
following.  It remains a calibration domain, not the independently defined
class required for a transferred bound.  Call a marked polytope
**determinant--Heegner-complete and automorphically correctable on that
support** if:

1. one of the finitely many determinant-line overlattices of Theorem 2.1
   passes the hidden-sheet gate of Theorem 3.2;
2. the selected facet Heegner types exhaust the full reflection branch shell
   that is being claimed, with no omitted branch component;
3. the integral facet multiplicity vector satisfies the finite obstruction
   equation (5.5).

This definition is surface-free and membership is decidable from the marked
polytope, its arithmetic lattice and discriminant module, and vector-valued
Hilbert modular forms.  But because it selects inputs by success of a
specified automorphic gate, it supplies no universality leverage by itself.
For this fixed-support success locus the following route is proved:

\[
 \boxed{
 P\longmapsto
 (\mathscr X_\Gamma,\{Z_i\})
 \longmapsto
 (\Sigma,\{C_i\}),
 }
 \tag{7.1}
\]

with exact marked mirror recovery on the real slice, a finite hidden-sheet
test, dihedral recovery of finite Coxeter labels on the stack, a meromorphic modular
form with the prescribed divisor, and the canonical/log-Chern identities of
Section 6.

Theorems 5.3--5.5 now prove that, after the explicit per-polytope direct level,
a bounded-height finite complex compensation always exists and a real
q-structure makes it nonnegative.  Over `Q` its scaled trace height has the
polynomial bound (5.18); over a higher-degree totally real field its exact
finite bound additionally contains the rank-enumerated term `R_det`.  The
integral effective Borcherds
upgrade is unconditional over `Q` in the interior and, for general totally
real `F`, follows from Bruinier's stated Hilbert q-integrality input; the
degree-greater-than-one anisotropic quotient then has no cusp boundary.  Thus
failure of one **fixed** support is no longer evidence that no compensated
carrier exists.  In particular there is still no
implication

\[
 \text{failure of (5.5)}\Longrightarrow
 \text{nonexistence of }P,
 \tag{7.2}
\]

The automatic carrier explains why this implication cannot come from the
finite cusp obstruction alone: a solution may escape a prescribed support,
and the per-input height and level bounds need not be uniform over a class.
What is still missing is a **class-uniform bounded-height and bounded-degree**
effective divisor package that an independent AG inequality or classification
can use, together with a written Hilbert q-integrality proof where that
arithmetic upgrade is needed.

### 7.1 Exact compact calibration: Bugaenko's $H^8$ chamber

The canonical compiler has now been run, without a supplied surface, on the
marked arithmetic Bugaenko polytope.  With $O=\mathbf Z[\phi]$ it reconstructs

\[
 L=\langle2\phi\rangle\perp E_8(-1),\qquad
 K_L=\langle4\phi\rangle,\qquad
 \widehat L=\langle2\phi\rangle\perp E_8(-1)\perp\langle4\phi\rangle.
 \tag{7.3}
\]

The signatures are $(2,8)$ and $(0,10)$ at the two real embeddings.  Moreover

\[
 A_{\widehat L}\cong O/(2)\oplus O/(4),\qquad |A_{\widehat L}|=64,
 \tag{7.4}
\]

and exhaustive evaluation of
$q(au+bv)=(2a^2+b^2)/(8\phi)$ finds no nonzero isotropic element.  Hence the
direct sum is the unique even determinant-line gluing; this is an exhaustive
finite conclusion, not a denominator cutoff.

All eleven wall vectors extend primitively with square $-2$, divisibility one,
and numerical Heegner type $(m,\mu)=(1,0)$.  Stack crossings recover
$44I_2(2)+8I_2(3)+2I_2(5)$, while the remaining facet pair is ultraparallel
and stays in the marked Gram.  Independent finite-parabolic root closure at
all 41 finite vertices gives the exact closed-slice split

\[
 320=11\ \text{facets}+309\ \text{lower-face contacts}
     +0\ \text{interior cuts}.
 \tag{7.5}
\]

Thus the complete numerical Heegner shell is not the facet-only divisor on
the closed chamber.  Nevertheless it does not cut the open chamber.  This is
the first exact compact $P\to$ projective-AG-arrangement calibration of the
general theorem.

The Whittaker gate is now exact as well.  In Bruinier's normalization all
eleven facets collapse to the single standard full-group/discriminant-module
column

\[
 (m_0,\mu_0)=(\phi/\sqrt5,0),
 \qquad k=(-3,5),\quad\kappa=(5,5).
 \tag{7.6}
\]

For the corresponding Hilbert--Weil Poincare series, exact cancellation of
the ideal-norm-four row and certified Kloosterman--Bessel tail estimates give

\[
 \left|a_{P_{m_0,0}}(m_0,0)-1\right|
 <0.742149649308345<1.
 \tag{7.7}
\]

Hence $P_{m_0,0}\ne0$, and the pure norm-two Whittaker principal part is
obstructed.  This proves that pure-shell automatic correction already fails
for the first compact benchmark.  It does not exclude a finite principal
part with compensating Heegner types, a refinement at finite level, or an
arbitrary meromorphic automorphic form.  The statement that the eleven
facets are not separately selectable is scoped to the standard column
$Z(m_0,0)$, which also contains the 309 lower-face contacts on the closed
slice.

The first compensator support is now closed exactly.  Enumerating all 64
discriminant components shows that the complete globally minimal positive
index shell, modulo positive units and component sign, is one
three-component orbit represented by

\[
 (m_1,\mu_1)=
 \left(\frac{\phi}{4\sqrt5},(\phi,2\phi^2)\right),
 \qquad N(m_1)=\frac1{80}.
 \tag{7.8}
\]

An explicit vector proves that its Heegner divisor is nonempty.  Every
component has nonzero coordinate on the determinant-line factor modulo four,
so every defining vector has nonzero auxiliary part.  The complete divisor
is therefore disjoint from the canonical real slice:

\[
 0\text{ facets}+0\text{ lower-face contacts}+0\text{ interior cuts}.
 \tag{7.9}
\]

After the exact norm similarity scaling, certified coefficient bounds are
`<3/4` and `<1/500` on the two diagonal errors and `<1/6` on both cross
errors.  The two row-dominance margins are `1/12` and `1247/1500`, so the
target-plus-shell `2x2` Poincare matrix is invertible.  Consequently no
Whittaker principal part on this complete first support can have nonzero
target coefficient.  Theorem 5.3 now proves abstractly that some finite
auxiliary-nonzero formal support cancels the target obstruction, and the
explicit per-polytope level can separate the facet trace.  Conditional here
on Bruinier's stated
Hilbert integral q-structure, Theorem 5.4 further proves that some such
correction can be chosen nonnegative integral, hence effective.  What remains
open is the numerical least support and the actual value of its
degree-two-field detection cutoff `R_det`, the corresponding obstruction
matrix, weight and global degree, and a nonempty-only support statement if
degree bookkeeping requires one.  Once that finite matrix is known,
Proposition 5.7 bounds its conic denominators, weight, and degree rather than
leaving three further existence problems.  The sheet-separating level itself is now computable from this
fixed marked chamber.  Thus the `2x2` result is a sharp first-shell
obstruction, not a global no-compensation theorem.

This remains a calibration, not a new abstract polytope or a transferred
dimension/facet bound.  The replays and certificates are in the
[`Bugaenko determinant-line report`](bugaenko_h8_determinant_line_ag.md) and
the [`exact Whittaker report`](bugaenko_h8_whittaker_gate.md), plus the
[`first compensator report`](bugaenko_h8_whittaker_first_compensator.md).

## 8. The proved carrier bridge, two missing terminal controls, and the proved return

Theorems 5.3--5.5 and Proposition 5.6 close the qualitative carrier part of
the old Missing Lemma A: target-preserving finite slice-disjoint compensation
is automatic on a sufficiently deep direct determinant level, and a real
q-structure makes it nonnegative.  An integral q-structure upgrades this to
an effective interior Borcherds divisor; in the compact anisotropic case it
has no boundary-cusp caveat.  The route becomes theorem-level leverage on
polytopes only after the following **quantitative** controls are proved.  They
are stated separately so that no desired conclusion is hidden inside the
definition of a carrier.

### Missing Lemma A+: bounded-height/degree carrier or exhaustive alternative

Find a non-singleton class $\mathcal C_0$ of compact arithmetic Coxeter
polytopes defined solely by intrinsic polytope/arithmetic data.  For every
$P\in\mathcal C_0$, strengthen the automatic effective interior carrier by producing
a determinant-level package and support $S(P)$ for which:

1. the compensator indices and discriminant level have a uniform explicit
   bound in the intrinsic parameters of $\mathcal C_0$ (the qualitative
   theorem bounds support cardinality only after such columns are found);
2. the finite integral obstruction matrix has a uniform height bound;
   Proposition 5.7 then supplies explicit automorphic-weight and global-degree
   bounds from its minors, the Eisenstein coefficients, and the arithmetic
   volume.  For $F=\mathbb Q$, compactification-boundary orders are controlled
   as well;
3. the needed rational/integral q-structure is proved in the relevant field
   family, rather than merely quoted as an adaptation; if the AG argument
   counts compensators, the nonempty retained types and their degrees are
   controlled; and
4. the resulting divisor package is sufficient for an AG inequality or
   classification that transfers back to the same $P$;

or place $P$ in a jointly exhaustive non-Borcherds carrier branch.

The sign and fixed-level height conditions are no longer missing.  If `b_T` is the target
obstruction column and `b_s` are safe columns, Theorem 5.4 proves, over a real
q-structure,

\[
          -b_T\in\operatorname{cone}_{\mathbb Q_{\ge0}}
                    \{b_s:s\text{ safe}\},
\tag{8.1}
\]

after fixing which sign means a zero rather than a pole, with $\mathbb Q$
replaced by $\mathbb R$ for a merely real structure.  A rational/integral
q-structure gives the displayed rational cone and denominator clearing.  The
theorem gives a cardinality bound `<=dim V_N`.  Theorem 5.5 also bounds the
largest scaled trace index by `B_CLW` over `Q`, and by
`max(R_det,B_CLW)` in higher degree.  What remains in A+ is uniformity over
the independently defined class: `N_sep(P)` may grow with `P`, no formula
for the Hilbert detection term `R_det` is proved, and the height `R` of the
finite integral obstruction columns is not uniformly bounded.  Proposition
5.7 proves that once `R` is bounded, Cramer minors control denominator
clearing and Bruinier's Eisenstein/degree identity controls both `w` and
`w*lambda^n`.  Over `Q` the needed Eisenstein coefficient and
Hirzebruch--Mumford volume formulas are explicit.  Thus weight, denominator,
and degree are no longer three independent missing mechanisms; the residual
problem is uniform matrix/level/volume control.

Pure-shell automatic correction cannot replace A+: it fails for the
Bugaenko chamber by (7.7), in addition to Wang's rational no-correction
examples.  Bruinier's 2017 prescribed-divisor theorem likewise produces
*some* nonzero-weight Borcherds product on an infinite admissible support in
the rational orthogonal setting, but does not supply the target-preserving
weight and degree bounds needed here.  Theorem 5.5 says how high a safe
positive spanning set must be sought at a fixed level over `Q`, and gives an
exact rank-enumerated stopping rule in higher degree.  Proposition 5.7 turns
that finite set into a weight/degree bound conditional on its integral
column height.  Neither result yet supplies the class-uniform matrix bound.

### Missing Lemma B: uniform finite-support obstruction

For fixed $F$, discriminant representation $A$, and finite allowed shell $S$,
normalize the Hilbert Poincare cusp series $G_s$ seeded at $s\in S$.  Prove an
explicit estimate making the coefficient matrix

\[
 \bigl(b_{G_s}(t)\bigr)_{s,t\in S}
 \tag{8.2}
\]

diagonally dominant, hence invertible, once the parallel obstruction weight
$\kappa=(n+2)/2$ is sufficiently large.  Analytically, the required input is
an explicit uniform bound for the off-diagonal Kloosterman--Bessel tail.

This would make $\operatorname{ob}_S$ injective in the proved range and would
exclude every nonzero $S$-supported Whittaker principal part there.  Theorem
5.3 shows that a target-preserving solution can then escape any such fixed
support.  Only a uniform height/level bound from Missing Lemma A+ can trap the
automatic solution inside the injective range and transfer the obstruction
back to $\mathcal C_0$.  No threshold, and therefore no dimension bound, is
asserted here.

### Lemma C (proved qualitatively): exact-level sheet-to-degree injection

On the exact-support level of Proposition 5.6, followed by a neat compact
cover, the map

\[
 \{\text{facets of }P\}\longrightarrow
 \{\text{irreducible components of the effective carrier divisor}\}
 \tag{8.3}
\]

is injective, and every image component contributes a positive integer to
its $\lambda$-degree.  Hence (6.9) gives

\[
 \#\operatorname{Facets}(P)
 \leq w\lambda^n/q_{\min}.
 \tag{8.4}
\]

The remaining quantitative problem is no longer a separate sheet-count
lemma: it is the level/weight/top-degree upper bound already demanded by
Missing Lemma A+.  At a coarse unseparated quotient the map need not be
injective, so its fixed divisor degree cannot be substituted into (8.4).
Orbifold Chern refinements using the factors $1/(2m_{ij})$ remain open if a
label-sensitive, rather than incidence-only, inequality is desired.

### Logical completion paths

- **Dimension path:** Missing Lemma A+ turns the qualitative positive carrier
  into a uniformly bounded list of allowed supports, while Missing Lemma B
  makes every corresponding
  coefficient matrix injective above an explicit threshold.  Their
  incompatibility would give a genuine nonexistence theorem for the
  independently defined class $\mathcal C_0$.
- **Facet path:** a Missing-Lemma-A+ bound on $w\lambda^n$ now feeds directly
  into (8.4), giving a wall bound for one Coxeter chamber.  A quantitative
  orbifold log-Chern inequality is needed only for a sharper label-sensitive
  refinement.
- **Classification path:** after a genuine dimension/discriminant bound,
  enumerate the finite determinant-line glue packages, run the hidden-sheet
  and Whittaker gates, and finish each survivor with literal Coxeter and
  Vinberg compactness/finite-volume tests.  Automorphic divisors do not replace
  those terminal hyperbolic tests.

## 9. Primary references

1. W. L. Baily Jr. and A. Borel,
   [*Compactification of arithmetic quotients of bounded symmetric domains*](https://doi.org/10.2307/1970449),
   Ann. of Math. **84** (1966), 442--528.  This supplies the algebraic
   quotient and projective minimal compactification.
2. J. H. Bruinier,
   [*Regularized theta lifts for orthogonal groups over totally real fields*](https://arxiv.org/abs/0908.3076),
   J. Reine Angew. Math. **672** (2012), 177--222,
   [doi:10.1515/CRELLE.2011.163](https://doi.org/10.1515/CRELLE.2011.163).
   See Proposition 4.5, Theorems 1.3 and 6.8, Corollary 7.4, and Section 2
   for the projectivity/anisotropy criterion.
3. J. H. Bruinier,
   [*Borcherds products with prescribed divisor*](https://doi.org/10.1112/blms.12090),
   Bull. London Math. Soc. **49** (2017), 979--987
   ([arXiv](https://arxiv.org/abs/1607.08713)).  Its Theorem 1 gives a
   nonzero-weight product supported on an infinite admissible set; it does
   not force a selected target coefficient.
4. W. J. McGraw,
   [*The rationality of vector valued modular forms associated with the Weil
   representation*](https://doi.org/10.1007/s00208-003-0413-1), Math. Ann.
   **326** (2003), 105--122.  This proves the integral Fourier-basis theorem
   in the classical rational setting; Bruinier invokes an adaptation for the
   Hilbert--Weil space without reproducing its proof.
5. M. Knopp, W. Kohnen, and W. Pribitkin,
   [*On the signs of Fourier coefficients of cusp forms*](https://doi.org/10.1023/A:1026207515396),
   Ramanujan J. **7** (2003), 269--277.  Their theorem applies to positive
   real weight with a unitary multiplier on a finitely generated Fuchsian
   group of the first kind and supplies the contradiction in Theorem 5.4.
6. V. Gritsenko, K. Hulek and G. K. Sankaran,
   [*Hirzebruch--Mumford proportionality and locally symmetric varieties of
   orthogonal type*](https://arxiv.org/abs/math/0609774),
   Doc. Math. **13** (2008), 1--19.  This is a primary source for canonical
   bundle and branch-divisor calculations on orthogonal modular varieties.
7. V. Gritsenko and K. Hulek,
   [*Uniruledness of orthogonal modular varieties*](https://arxiv.org/abs/1202.3361),
   J. Algebraic Geom. **23** (2014), 711--725,
   [doi:10.1090/S1056-3911-2014-00632-9](https://doi.org/10.1090/S1056-3911-2014-00632-9).
   See Theorem 2.1 for the strongly reflective weight criterion.
8. A. Langer,
   [*Logarithmic orbifold Euler numbers of surfaces with applications*](https://arxiv.org/abs/math/0012180),
   Proc. London Math. Soc. **86** (2003), 358--396,
   [doi:10.1112/S0024611502013874](https://doi.org/10.1112/S0024611502013874).
   This supplies the orbifold logarithmic BMY framework.
9. H. Wang,
   [*The classification of 2-reflective modular forms*](https://arxiv.org/abs/1906.10459),
   J. Eur. Math. Soc. **26** (2024), 111--151,
   [doi:10.4171/JEMS/1358](https://doi.org/10.4171/JEMS/1358).
   See Theorem 1.1 for the high-rank classification and Theorem 7.2 for the
   eighteen hyperbolically reflective lattices without automorphic correction.
10. I. Metzler,
   [*A new converse theorem for Borcherds products*](https://arxiv.org/abs/2512.20701)
   (2025 preprint), especially the hyperbolic-split converse in Section 6 and
   the warning that the split hypothesis is essential.
11. M. Zhang and Y. Zhang,
   [*On the non-vanishing of Hilbert Poincare series*](https://arxiv.org/abs/2312.04169),
   J. Number Theory (2026).  Its Hilbert Kloosterman--Bessel coefficient
   formula underlies the explicit Bugaenko diagonal estimate.
12. F. Strömberg,
   [*Weil representations associated to finite quadratic modules*](https://arxiv.org/abs/1108.0202),
   Math. Z. **275** (2013), 509--527.  Its explicit finite-Weil matrix
   formula supplies the kernel-size bounds used in that calibration and,
   especially in arXiv Lemma 5.13 (journal Lemma 5.6), the monomial
   `Gamma_0(ell)` action used in
   Theorem 5.5.
13. G. Chen, Y.-K. Lau, and Y. Wang,
   [*On the first sign change of Fourier coefficients of cusp forms*](https://preprint.press.jhu.edu/ajm/sites/default/files/AJM-chen-lau-wang.pdf),
   Theorem 1.  This supplies the `(kN)^(2+epsilon)` scalar first-sign-change
   bound used after orbit-square scalarization.
14. M. Schwagenscheidt,
   [*Eisenstein series for the Weil representation*](https://arxiv.org/abs/1803.10550),
   J. Number Theory **193** (2018), 74--90.  Proposition 2.2 and unitarity
   give the coefficient majorant (5.25); Theorem 1.4 gives the exact local
   representation-number formula.
15. V. Gritsenko, K. Hulek, and G. K. Sankaran,
   [*The Hirzebruch--Mumford volume for the orthogonal group and
   applications*](https://arxiv.org/abs/math/0512595), Doc. Math. **12**
   (2007), 215--241.  Their main formula and Proposition 1.2 give the
   rational determinant/local-density volume and `lambda^n=2 vol_HM`.
16. G. Prasad,
   [*Volumes of S-arithmetic quotients of semi-simple groups*](https://doi.org/10.1007/BF02698841),
   Publ. Math. IHES **69** (1989), 91--117.  This supplies the general
   arithmetic covolume framework used in the totally-real branch.
17. H. Wang,
   [*$2$-reflective lattices of signature $(n,2)$ with $n\geq8$*](https://arxiv.org/abs/2301.11536).
   Theorem 1.2 gives the complete forty-two-lattice list and, in particular,
   no $2$-reflective lattice of signature $(13,2)$; its complete-reflective
   corollary contains the positive $2U\oplus E_7$ control.

## 10. Frozen conclusion

Every marked compact arithmetic Coxeter polytope of the stated totally-real
type has a canonical determinant-line path to a type-IV arithmetic stack.
Its facet mirrors occur as exact real traces of special divisors; a finite
calculation detects hidden arithmetic sheets; marked stack inertia recovers
the finite Coxeter edge labels; and a finite vector-valued Hilbert cusp-form pairing
decides whether a prescribed finite divisor is carried by Bruinier's
Whittaker lift.  General complete intersections then produce algebraic
surfaces with computable canonical and log-Chern data.

The finite-support existence and sign problems are now solved at their proper
arithmetic levels: auxiliary-nonzero components span the cusp obstruction
dual, their nonnegative cone is the whole real obstruction dual, and a deep
level isolates the target facet sheets.  Thus a real q-structure gives a
positive finite correction, while an integral q-structure gives an effective
interior Borcherds divisor after scaling.  The latter is unconditional over
`Q` by McGraw and conditional for general totally real `F` on Bruinier's
stated Hilbert adaptation; the anisotropic degree-greater-than-one quotient
then has no boundary cusp.

Exact level separation plus ampleness of the weight-one bundle also closes
qualitative Lemma C:
\[
 \#\operatorname{Facets}(P)\,q_{\min}
 \leq w\lambda^n.
\]
What is not known is exactly the quantitative input that would create
polytope-level leverage: a **class-uniform** bound on sheet-separating level,
the higher-degree detection term, the finite integral obstruction-column
height, and the arithmetic volume (or an exhaustive alternative carrier),
plus a compatible large-weight obstruction with explicit constants.
Proposition 5.7 then bounds conic denominators, weight, and the top
automorphic degree $w\lambda^n$ from those data.  At each
fixed level the scaled trace height is now bounded polynomially over `Q` and
by an exact finite rank-enumerated cutoff in higher degree.  Bugaenko proves
that the pure and first-shell bounded
versions fail, while Theorems 5.3--5.5 show that some higher positive solution
must escape those supports at the relevant q-structure level.  The corrected
remaining demands are Missing Lemmas A+ and B; the quantitative content of
old Lemma C and the weight/degree identities have folded into A+'s uniform
level/detection/matrix/volume bound.
Until one of the stated completion paths is closed, this transfer is a
rigorous universal complex Whittaker/Green carrier and, with the stated
arithmetic scope, an effective exact-slice automorphic carrier.  It is not a
new bound or a nonexistence theorem for compact arithmetic Coxeter polytopes.
