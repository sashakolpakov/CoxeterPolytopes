# Slice-disjoint Fourier coefficients span the Whittaker obstruction space

**Status (2026-08-09):** proved at the complex weak-Whittaker/Green-current
level, with an exact finite-level separation of the selected facet sheets.
At the real/integral q-structure level, a diagonal sign-change argument also
proves that the safe obstruction columns positively span, giving an effective
interior Borcherds divisor.  This arithmetic upgrade is unconditional for
\(F=\mathbf Q\) by McGraw.  For general totally real \(F\), it uses an integral
Fourier-structure assertion that Bruinier states, but does not prove in detail,
by saying that McGraw's argument adapts to the Hilbert--Weil setting.  We
therefore record that upgrade separately rather than folding it into the
unconditional theorem.

An orbit-square scalarization now makes the positive-cone theorem
quantitative at a fixed discriminant level.  Over \(F=\mathbf Q\), safe
columns of scaled index

\[
 \ell m\ll_\varepsilon(8K\ell^2)^{2+\varepsilon}
\]

already positively span, where \(K\) is the diagonal weight and \(\ell\) the
finite quadratic-module level.  Over a higher-degree totally real field the
same statement requires the additional finite safe-detection cutoff
\(R_{\rm det}\); it is obtained by exact coefficient-rank enumeration but is
not bounded here in terms of \(K\) and \(\ell\).  Thus this is a closed
polynomial height theorem over \(\mathbf Q\) and a per-space finite
certificate in general, not yet a uniform totally-real class bound.

This result closes a precise linear-algebraic part of Missing Lemma A in
[`polytope_to_automorphic_divisor_transfer_theorem.md`](polytope_to_automorphic_divisor_transfer_theorem.md).
It does **not** give a uniform class bound for the sheet-separating level,
the resulting automorphic weight, or top degree.  The separate exact-level sheet-to-degree theorem now
injects facets into divisor components once the arithmetic effective upgrade
is available.  Effectivity belongs only to that arithmetic upgrade, not to
the unconditional complex statement.

## 1. Setup

Let \(F\) be totally real and let

\[
 M=M_0\perp K
 \tag{1.1}
\]

be the sign-reversed, **unglued** determinant-line extension used in
Bruinier's convention.  Thus \(M\) has signature \((n,2)\) at one real place
and is positive definite at all other real places, while \(K\) is the
rank-one determinant summand.  Write trace duals as \(M'\), \(M_0'\), and
\(K'\).  Orthogonality gives a literal decomposition of finite quadratic
modules

\[
 A_M=M'/M=A_0\oplus B,
 \qquad A_0=M_0'/M_0,\quad B=K'/K.
 \tag{1.2}
\]

The finite bilinear pairings on \(A_0\) and \(B\) are nondegenerate.  Let
\(e\in K_F\) be the vector used in the canonical real slice; every plane in
that slice contains \(e\).

For \(\mu=(\alpha,\beta)\in A_0\oplus B\), call \(\mu\) **safe** if
\(\beta\ne0\).  Let

\[
 \mathcal S=S_{\kappa,\rho_M}
 \tag{1.3}
\]

be Bruinier's finite-dimensional obstruction cusp space.  If

\[
 g=\sum_{\mu\in A_M}\sum_{m\gg0}
 b_g(m,\mu)e(\operatorname{tr}(m\tau))\chi_\mu\in\mathcal S,
 \tag{1.4}
\]

denote the coefficient functional \(g\mapsto b_g(m,\mu)\) by
\(\ell_{m,\mu}\).

The direct-sum hypothesis is important.  An arbitrary isotropic gluing of
\(M_0\perp K\) need not retain a well-defined nonzero auxiliary coordinate.
There is no loss of existence here: the zero gluing, namely the canonical
unglued determinant extension, is always one of the allowed extensions.

## 2. The safe-component spanning theorem

### Theorem 2.1 (safe coefficients separate cusp forms)

Assume \(B\ne0\).  Then

\[
 \operatorname{span}_{\mathbf C}
 \{\ell_{m,(\alpha,\beta)}:\ m\gg0,\ \beta\ne0\}
 =\mathcal S^*.
 \tag{2.1}
\]

Consequently a finite subset of safe coefficient functionals spans
\(\mathcal S^*\).  It can be chosen with at most \(\dim\mathcal S\) members,
but this statement supplies no a priori bound for their indices \(m\).

#### Proof

It suffices to prove that the common annihilator of the displayed
functionals is zero.  Suppose that \(g\in\mathcal S\) has

\[
 b_g(m,(\alpha,\beta))=0
 \quad\text{for every admissible }m\gg0\text{ and every }\beta\ne0.
 \tag{2.2}
\]

Because \(g\) is cuspidal, these are all Fourier terms of its safe component
functions.  Hence

\[
 g_{\alpha,\beta}(\tau)=0\qquad(\beta\ne0).
 \tag{2.3}
\]

Bruinier's formula for the Weil action of
\(S=\left(\begin{smallmatrix}0&-1\\1&0\end{smallmatrix}\right)\) is, up to a
nonzero scalar independent of \(\mu,\nu\),

\[
 \rho_M(S)\chi_\mu
 =|A_M|^{-1/2}\sum_{\nu\in A_M}
 e\!\left(-\operatorname{tr}(\mu,\nu)\right)\chi_\nu.
 \tag{2.4}
\]

Fix any \(b\in B\setminus\{0\}\).  The \((a,b)\)-component of
\(g(S\tau)\) is zero by (2.3).  Applying modularity and using (1.2), it is a
nonzero automorphy scalar times

\[
 \sum_{\alpha\in A_0}
 e\!\left(-\operatorname{tr}(a,\alpha)\right)g_{\alpha,0}(\tau).
 \tag{2.5}
\]

This vanishes for every \(a\in A_0\).  The matrix in (2.5) is the finite
Fourier transform of the nondegenerate finite bilinear module \(A_0\), hence
is invertible.  Therefore every \(g_{\alpha,0}\) vanishes.  Together with
(2.3), this gives \(g=0\).  The annihilator criterion proves (2.1), and
finite dimensionality gives a finite spanning subset. \(\square\)

### Lemma 2.2 (the determinant line always has a safe component)

For the canonical determinant line with form

\[
 b_K=2(-1)^n\det(b_L),
 \tag{2.6}
\]

one has \(B=K'/K\ne0\).

#### Proof

Localize at a prime \(\mathfrak p\mid2\).  The projective lattice \(L\) and
its determinant line are free locally.  In local generators the determinant
line Gram entry is

\[
 a=2(-1)^n\det(G_L).
 \tag{2.7}
\]

The entries of \(G_L\) are integral, so \(a\) is a nonunit at
\(\mathfrak p\).  The ordinary \(\mathcal O_{F,\mathfrak p}\)-dual
\(a^{-1}\mathcal O_{F,\mathfrak p}\) therefore strictly contains the local
line.  The trace dual contains the ordinary integral dual, so \(K'/K\) is
nonzero. \(\square\)

### Orbit audit

The safe set is stable under \(\mu\mapsto-\mu\) and under multiplication by
totally positive units: a unit cannot send a nonzero element of \(B\) to
zero.  Hence the central \(\mu\leftrightarrow-\mu\) relation and Bruinier's
unit-orbit convention only remove redundant columns; they do not invalidate
the separation proof.  One may equivalently take complete sign/unit orbits
of safe indices.

## 3. Why safe Heegner divisors miss the real slice

### Proposition 3.1 (whole-divisor slice disjointness)

Work on the discriminant-kernel quotient, or on any quotient whose action on
\(A_M=A_0\oplus B\) preserves the subset \(\beta\ne0\).  If
\(\mu=(\alpha,\beta)\) is safe, then

\[
 j_e(P)\cap Z(m,\mu)=\varnothing.
 \tag{3.1}
\]

This includes the possibility that \(Z(m,\mu)=0\) globally.

#### Proof

Every vector \(v\) defining a sheet of \(Z(m,\mu)\) on the indicated cover
has finite congruence class \(\mu\).  Its \(K_F\)-component cannot be zero,
since zero has auxiliary class \(0\), whereas \(\beta\ne0\).  Because \(K_F\)
has rank one and is nondegenerate,

\[
 (v,e)\ne0.
 \tag{3.2}
\]

Every two-plane in the canonical real slice contains \(e\), while a plane
on the mirror \(\mathcal D_v\) must be perpendicular to \(v\).  Thus no such
plane lies on \(\mathcal D_v\).  The safe subset is invariant under the
quotient group by hypothesis, so the same argument applies to every sheet in
the quotient pullback. \(\square\)

Passing blindly to a quotient that mixes the zero and nonzero auxiliary
classes would destroy this conclusion.  The discriminant-kernel cover is the
canonical way to avoid that problem.  It is finite and remains projective in
the anisotropic totally-real case.

## 4. A deep-level theorem isolates exactly the chosen facets

The safe-span theorem does not by itself remove other zero-auxiliary sheets
from a target Heegner class.  A finite congruence cover does.

### Theorem 4.1 (finite-level exact closed-slice support)

Let \(r_1,\ldots,r_s\in M_0\subset M\) be the marked facet vectors and let
\(P\) be compact.  There is a positive integer \(N\) such that, for

\[
 M_N=NM=NM_0\perp NK,
 \qquad \mu_i=r_i+M_N\in M_N'/M_N,
 \tag{4.1}
\]

the following hold.

1. Every \(r_i\) belongs to \(M_N'\).
2. The classes \(\mu_i\) are pairwise distinct modulo sign, after duplicate
   facet mirrors have been removed.
3. On the discriminant-kernel cover for \(M_N\), the support of
   \(Z(Q(r_i),\mu_i)\) pulled back to the **closed** chamber \(P\) is exactly
   the facet mirror \(r_i^\perp\cap P\).  In particular there are no extra
   lower-face contacts and no interior-cutting sheets in that Heegner class.
4. The safe-component theorem and Proposition 3.1 apply to
   \(A_{M_N}=A_{NM_0}\oplus A_{NK}\), and \(A_{NK}\ne0\).

The claim is about support.  It does not silently choose a normalization for
the cycle multiplicity of \(Z(m,\mu)=Z(m,-\mu)\).

#### Proof

Let \(T=\{Q(r_i):1\le i\le s\}\).  The set

\[
 \mathcal U_T=
 \{v\in M:Q(v)\in T,\ \mathcal D_v\cap j_e(P)\ne\varnothing\}/\{\pm1\}
 \tag{4.2}
\]

is finite.  Indeed, compactness of \(P\) gives a positive-definite majorant
bound at the distinguished real place for vectors of any fixed norm whose
mirrors meet \(P\).  At every other real place the quadratic space is
positive definite and the conjugates of \(Q(v)\) are fixed by \(T\).  The
Minkowski embedding of \(M\) is discrete, so the resulting bounded set
contains only finitely many lattice vectors.

Choose a \(\mathbf Z\)-basis of the underlying free abelian group of \(M\).
There are only finitely many nonzero differences

\[
 v-r_i,\quad v+r_i
 \quad(v\in\mathcal U_T,\ v^\perp\ne r_i^\perp,
       \ Q(v)=Q(r_i)),
 \qquad r_i-r_j,\quad r_i+r_j
 \tag{4.3}
\]

that must be kept out of \(NM\).  Since
\(\bigcap_{N\ge1}NM=\{0\}\), one can choose one sufficiently large \(N\)
for which none of the differences in (4.3) belongs to \(NM\).  This proves
the sign-distinctness and the desired congruence separation.

This choice is effective.  In the fixed \(\mathbf Z\)-basis, let \(C\) be
the largest absolute coordinate among the nonzero differences in (4.3).
Then

\[
                 N_{\rm sep}=C+1
 \tag{4.3a}
\]

works, with \(N_{\rm sep}=1\) if the difference set is empty.  A nonzero
vector in \(N_{\rm sep}M\) has a nonzero coordinate of absolute value at
least \(N_{\rm sep}\).  If \(\ell_0\) is the level of the
base trace discriminant module and \(\ell_N\) that of \(NM\), then

\[
                 \ell_N\mid N^2\ell_0.
 \tag{4.3b}
\]

Indeed \((NM)'=N^{-1}M'\) and
\(N^2\ell_0Q(y/N)=\ell_0Q(y)\in\mathbf Z\).
The basis-coordinate construction and its exact scope are isolated in
[`explicit_sheet_separating_level_bound.md`](explicit_sheet_separating_level_bound.md).

Since \(r_i\in M\) and \(M\) is integral,
\((r_i,NM)\subset\mathbf Z\), so \(r_i\in M_N'\).  Any sheet of
\(Z(Q(r_i),\mu_i)\) meeting \(P\) is defined by a vector

\[
 v\in r_i+NM,\qquad Q(v)=Q(r_i).
 \tag{4.4}
\]

It belongs to \(\mathcal U_T\), and the choice of \(N\) forces its mirror to
be the mirror of \(r_i\).  Conversely \(r_i\) itself supplies that mirror.
This proves the exact support claim, including boundary contacts because
\(P\) was kept closed.

Finally, (4.1) remains an unglued orthogonal sum.  The proof of Lemma 2.2,
or simply localization at a prime dividing \(2N\), gives
\(A_{NK}\ne0\).  Theorems 2.1 and Proposition 3.1 therefore apply at this
level. \(\square\)

Wall reflections need not lie in the discriminant-kernel group of this deep
cover.  The theorem makes no such claim.  The marked Gram matrix still
records the Coxeter labels externally; this finite-level separation is not a
replacement for the stack-inertia construction on a quotient containing the
reflections.

## 5. The two theorem levels

### Level I: unconditional complex weak-Whittaker compensation

Let \(p_{\mathrm{fac}}\) be any finite complex principal part supported on
the separated target classes of Theorem 4.1.  Its obstruction is an element

\[
 \operatorname{ob}(p_{\mathrm{fac}})\in\mathcal S^*.
 \tag{5.1}
\]

By Theorem 2.1 there are finitely many safe indices \(t_1,\ldots,t_a\) and
complex numbers \(c_1,\ldots,c_a\) such that

\[
 \operatorname{ob}\!\left(
 p_{\mathrm{fac}}+\sum_{j=1}^a c_j f_{t_j}
 \right)=0.
 \tag{5.2}
\]

Bruinier's Proposition 4.5 then gives a weakly holomorphic Whittaker form.
Its target Heegner support pulls back to exactly the chosen facets, and all
compensating Heegner divisors miss the real slice.  For arbitrary complex
coefficients, the unconditional geometric output is the regularized theta
lift/Green-current relation with a complex special-divisor combination.

This is a genuine universal finite-support **existence** theorem.  The
quantitative refinement and its exact field-degree scope are given in
Section 6 below.

### Level II: integral meromorphic Borcherds upgrade

If \(\mathcal S\) has a basis whose Fourier coefficients in the standard
\(\chi_\mu\) basis are rational, the finite full-rank safe coefficient
matrix in Theorem 2.1 can be chosen rational.  For an integral target column,
(5.2) then has a rational solution.  Clearing denominators gives an integral
weakly holomorphic Whittaker principal part with nonzero (possibly scaled)
target multiplicity.  Bruinier's Theorem 6.8 produces a meromorphic modular
form with that divisor and finite-order multiplier.

Bruinier writes in the proof of his Theorem 7.1:

> Adapting the argument of McGraw, it can be proved that
> \(S_{\kappa,\rho_L}\) has a basis of cusp forms with coefficients in
> \(\mathbf Z\).

If that published assertion is accepted as the Hilbert--Weil integral
q-structure theorem, the preceding paragraph upgrades Level I to an actual
integral meromorphic carrier on the deep cover.  However, McGraw's proved
theorem is for the classical \(F=\mathbf Q\) Weil representation, and
Bruinier does not include the adaptation in this paper.  Our archive contains
no independent proof of the Hilbert adaptation.  Accordingly:

- for \(F=\mathbf Q\), McGraw's theorem supplies the arithmetic upgrade;
- for general totally real \(F\), the meromorphic upgrade is conditional on
  the quoted Hilbert integral-structure assertion until that adaptation is
  written out or replaced by a precise q-expansion theorem.

Complex spanning alone is not enough: a complex linear kernel can fail to
contain any nonzero rational vector.  Nor can multiplying arbitrary complex
coefficients make them integral.

### Theorem 5.1 (the safe columns positively span over a real q-structure)

Assume that \(\mathcal S\) has a basis \(g_1,\ldots,g_d\) whose Fourier
coefficients are real.  Identify the real form

\[
 \mathcal S_{\mathbf R}=\operatorname{span}_{\mathbf R}\{g_1,\ldots,g_d\}
 \tag{5.3}
\]

and its dual with \(\mathbf R^d\).  Let \(C_{\rm safe}\) be the convex cone
generated by all safe coefficient columns.  Then

\[
                       C_{\rm safe}=\mathcal S_{\mathbf R}^*.
 \tag{5.4}
\]

More strongly, one finite set of safe columns already has positive cone equal
to the whole dual.  Hence for every integral target principal part there is a
finite nonnegative real safe correction.  If the chosen q-structure is
integral, the correction can be chosen rational and, after scaling the target,
nonnegative integral.  It may be chosen with at most \(d\) nonzero safe
coefficients.

#### Proof: diagonal restriction and sign change

Suppose that a nonzero \(g\in\mathcal S_{\mathbf R}\) had

\[
 b_g(m,(\alpha,\beta))\geq0
 \qquad\text{at every safe index.}
 \tag{5.5}
\]

Theorem 2.1 implies that some safe component \(g_\mu\),
\(\mu=(\alpha,\beta)\), \(\beta\ne0\), is nonzero.  Restrict this component
to the diagonal:

\[
 h_\mu(z)=g_\mu(z,\ldots,z)
 =\sum_{m\gg0}b_g(m,\mu)e(\operatorname{Tr}_{F/\mathbf Q}(m)z).
 \tag{5.6}
\]

Restrict the Hilbert metaplectic group to the diagonal copy of
\(\operatorname{Mp}_2(\mathbf Z)\), and then pass to the finite-index kernel
of the finite Weil representation.  On this subgroup an individual component
transforms as a scalar cusp form of positive real weight

\[
                 \kappa_\Delta=\sum_{j=1}^{[F:\mathbf Q]}\kappa_j>0
 \tag{5.7}
\]

with a unitary multiplier.  The subgroup is finitely generated and of
Fuchsian type of the first kind.  The restriction is cuspidal at every one of
its cusps: a rational cusp of the diagonal subgroup is a Hilbert cusp, and
after applying a diagonal modular transformation each component is a linear
combination of components of the original Hilbert cusp form.

The coefficient of a fixed diagonal exponent is a **finite** sum.  Indeed,
the totally positive \(m\) of fixed trace lie in a bounded subset of the
Minkowski space of the shifted lattice
\(\partial_F^{-1}+Q(\mu)\), hence form a finite set.  Thus (5.5) makes every
Fourier coefficient of \(h_\mu\) nonnegative.  Since \(g_\mu\ne0\), at least
one summand is positive, so no cancellation can make \(h_\mu\) zero.  The
possible fractional shift \(\operatorname{Tr}Q(\mu)\bmod\mathbf Z\) is simply
the cusp parameter in the expansion

\[
 \sum_{r+\theta>0}a(r)e^{2\pi i(r+\theta)z/\lambda};
 \tag{5.8}
\]

it is included in the Knopp--Kohnen--Pribitkin theorem.

That theorem says that a nonzero cusp form of positive real weight with
unitary multiplier on such a group, if its Fourier coefficients are real,
has infinitely many positive and infinitely many negative coefficients.
This contradicts (5.5).  Therefore every nonzero real cusp form has at least
one negative safe coefficient.  Applying the same conclusion to \(-g\) also
gives at least one positive safe coefficient.

For completeness, this pointwise assertion gives one **finite** positive
spanning set.  On the unit sphere in \(\mathcal S_{\mathbf R}\), the open sets

\[
 U_s=\{g:b_g(s)<0\},\qquad s\text{ safe},
 \tag{5.9}
\]

cover the sphere.  Compactness supplies a finite subcover indexed by
\(S_0\).  If the cone generated by the columns in \(S_0\) were not the whole
dual, finite-dimensional separation would give a nonzero \(g\) nonnegative
on every one of those columns, contradicting the subcover.  Hence that finite
cone is the whole dual.

Conic Caratheodory now expresses any target functional using at most \(d\)
columns of \(S_0\).  If the basis and target coefficients are integral, the
resulting nonempty rational polyhedron

\[
 \{x\geq0:B_{S_0}x=-b_{\rm target}\}
 \tag{5.10}
\]

has a rational point (equivalently, take a rational basic feasible solution).
Clearing denominators proves the integral assertion. \(\square\)

### Corollary 5.2 (effective exact-slice carrier, with arithmetic scope)

Under the integral q-structure hypothesis, take the target facet
multiplicities to be positive.  Theorem 5.1 gives a positive integral multiple
of that target plus nonnegative integral safe compensators in the Whittaker
kernel.  Bruinier's Theorem 6.8 therefore gives a meromorphic modular form
whose interior divisor is effective, whose pullback to the closed real
chamber is exactly the selected positive facet multiple, and whose other
components miss the whole real slice.

The modular form has no pole on the type-IV domain.  When \([F:\mathbf Q]>1\)
the orthogonal quotient is anisotropic and projective, so there is no cusp
boundary at which an unrecorded pole could occur.  For \(F=\mathbf Q\), the
effective **interior** divisor does not exclude poles on a Baily--Borel or
toroidal boundary.

McGraw makes Theorem 5.1 and Corollary 5.2 unconditional for \(F=\mathbf Q\).
For a general totally real field they are conditional on Bruinier's quoted
Hilbert integral-basis assertion.  The orbit-square theorem below bounds the
largest **scaled trace index** at a fixed level with an essential
field-degree split.  It does not bound the conic denominators, automorphic
weight, or global Heegner degree.

## 6. Quantitative orbit-square theorem

Let \(A\) be the finite quadratic module of the underlying even trace
\(\mathbf Z\)-lattice and let \(\ell\) be its level.  Put

\[
 K=\sum_j\kappa_j=\frac{[F:\mathbf Q](n+2)}2,
 \qquad r(m)=\ell\operatorname{Tr}_{F/\mathbf Q}(m)\in\mathbf Z_{>0}.
 \tag{6.1}
\]

For a safe component \(x\in A\), let
\(\mathcal O_x=(\mathbf Z/\ell\mathbf Z)^\times x\) and set

\[
 H_{\mathcal O_x}(z)=
 \sum_{y\in\mathcal O_x}g_y(\ell z,\ldots,\ell z)^2.
 \tag{6.2}
\]

### Theorem 6.1 (scalarization and bounded safe cone)

The function (6.2) is a scalar cusp form of integral weight \(2K\) on
\(\Gamma_0(4\ell^2)\) with Dirichlet nebentypus.  It has real coefficients
for \(g\in\mathcal S_{\mathbf R}\), its orbit remains safe, and it is nonzero
whenever one of the restricted components in its orbit is nonzero.

For every \(\varepsilon>0\), let

\[
 B_{\rm CLW}=\left\lceil
 C_\varepsilon(8K\ell^2)^{2+\varepsilon}\right\rceil
 \tag{6.3}
\]

be the Chen--Lau--Wang first-sign-change bound.  If \(F=\mathbf Q\), the
finite safe columns with \(r(m)\le B_{\rm CLW}\) already have positive cone
\(\mathcal S_{\mathbf R}^*\).

For general totally real \(F\), define

\[
 R_{\rm det}=\min\left\{R:
 \operatorname{span}_{\mathbf R}
 \{\ell_{m,x}:x\text{ safe},\ r(m)\le R\}
 =\mathcal S_{\mathbf R}^*\right\}.
 \tag{6.4}
\]

This integer is finite by Theorem 2.1.  The columns with

\[
 r(m)\le R_*:=\max\{R_{\rm det},B_{\rm CLW}\}
 \tag{6.5}
\]

have positive cone \(\mathcal S_{\mathbf R}^*\).  In either case a correction
uses at most \(\dim\mathcal S\) columns; rational data give a rational
solution and denominator clearing gives a nonnegative integral correction
after scaling the target.

#### Proof and hostile audit

For \(\gamma\in\Gamma_0(4\ell^2)\), conjugation by
\(\operatorname{diag}(\ell,1)\) gives

\[
 \delta=\begin{pmatrix}a&\ell b\\c/\ell&d\end{pmatrix}
 \in\Gamma_0(\ell).
 \tag{6.6}
\]

Strömberg's arXiv Lemma 5.13 (journal Lemma 5.6) gives

\[
 \rho_A(\delta)e_x
 =e(\ell bdQ(x))\chi_A(\delta)e_{dx}
 =\chi_A(\delta)e_{dx}.
 \tag{6.7}
\]

Thus unit-orbit summation removes the permutation.  In even signature
\(\chi_A\) is a Dirichlet character; in odd signature it is a theta
multiplier or its conjugate times one, and its square is again Dirichlet.
Diagonal restriction is cuspidal on the finite Weil kernel because every
rational diagonal cusp is a Hilbert cusp.  Rescaling conjugates that
finite-index group; intersecting it with \(\Gamma_0(4\ell^2)\) still has
finite index.  Since (6.2) is modular on \(\Gamma_0(4\ell^2)\), its vanishing
at all rational cusps for the intersection proves cuspidality at every cusp
of the larger group.  Holomorphy is
preserved by restriction, rescaling, squaring, and finite summation.  On the
imaginary axis the component values are real, so a sum of squares cannot
vanish identically unless all restricted components vanish.

Over \(F=\mathbf Q\), diagonal restriction is the identity.  If a nonzero
real \(g\) were nonnegative on every safe column through
\(B_{\rm CLW}\), Theorem 2.1 would give a nonzero safe orbit-square.  Every
coefficient of its square-sum through \(B_{\rm CLW}\) is a sum of products
of bounded nonnegative safe coefficients.  Chen--Lau--Wang forces an
opposite-sign pair in that range, a contradiction.

For general \(F\), diagonal restriction may have a kernel and equal-trace
coefficients may cancel; this is why (6.4) is necessary.  If a nonzero
\(g\) were nonnegative on all columns through \(R_*\), linear spanning at
\(R_{\rm det}\) would make one bounded coefficient positive.  Its diagonal
trace coefficient is a finite sum of bounded nonnegative coefficients and is
therefore positive.  The same square-sum contradiction applies.  The polar
of the finite bounded cone is zero, hence the cone is the full dual.  Conic
Caratheodory and rational polyhedral linear algebra finish the proof.

The full derivation is recorded in
[`orbit_square_quantitative_safe_height.md`](orbit_square_quantitative_safe_height.md).

## 7. Nonempty divisors, signs, and what remains open

Bruinier defines \(Z(m,\mu)=0\) when the relevant global norm orbit is empty.
Theorem 2.1 ranges over all admissible Fourier indices, so some of its finite
spanning columns may have zero global divisor.  Such a column is vacuously
slice-disjoint and is harmless at Level I, but the proof does **not** show
that coefficient functionals attached only to nonempty safe divisors span
\(\mathcal S^*\).

Every safe discriminant coset does contain vectors of totally positive norm
after a suitable lattice translation: choose a lattice vector whose norm is
positive at the distinguished place (the other places are definite) and add
large multiples to a fixed coset representative.  Thus nonempty safe types
exist abundantly.  This fact does not imply that an arbitrary cusp form has a
nonzero coefficient at one of those represented indices.  Proving the
nonempty-only spanning statement would require an additional representation
or theta-lift injectivity theorem; the finite Fourier-transform argument
alone cannot supply it.

No sign or effectivity statement follows from the unconditional complex
relation (5.2): its coefficients can have either sign or be non-real.  The
real/integral q-structure hypothesis changes this conclusion through Theorem
5.1 and gives an effective interior divisor.  Theorem 6.1 controls trace
height polynomially over \(\mathbf Q\), and controls it in general after the
per-space detection term \(R_{\rm det}\).  It does not control the height of
the resulting finite integral obstruction matrix.  Once that height and the
arithmetic volume are bounded, Cramer--Hadamard minors and Bruinier's
Eisenstein/degree identity bound conic denominators, automorphic weight, and
top degree; see
[`borcherds_weight_degree_quantitative_bound.md`](borcherds_weight_degree_quantitative_bound.md).
Neither the matrix height, \(R_{\rm det}\), nor the sheet-separating level is
uniformly bounded in a totally-real polytope class.  Thus the effective divisor cannot yet be
inserted into a quantitative log-BMY or canonical-divisor argument.
Qualitative
Missing Lemma C is now closed separately: on the exact-support neat compact
cover,
\[
 \#\operatorname{Facets}(P)\,q_{\min}\leq w\lambda^n.
\]
The still-missing input is a uniform upper bound for the right side; see
[the sheet-to-degree report](facet_sheet_automorphic_degree_injection.md).

What is newly closed is exact and useful:

\[
 \boxed{
 \begin{gathered}
 \text{compact marked }P
 \ \longrightarrow\ \text{explicit per-}P\text{ determinant-line cover with exact facet
 support}\ \\
 \ \longrightarrow\ \text{finite slice-disjoint complex obstruction
 compensation.}
 \end{gathered}}
 \tag{7.1}
\]

With the real/integral q-structure, Theorem 5.1 upgrades the last line of
(7.1) to a nonnegative integral correction and Corollary 5.2 supplies the
effective meromorphic carrier.  Thus the remaining arithmetic step in the
general totally-real case is precisely the integral Hilbert--Weil q-structure
issue above, not a further geometric hidden-sheet problem.

## 8. Primary sources

1. J. H. Bruinier,
   [*Regularized theta lifts for orthogonal groups over totally real fields*](https://arxiv.org/abs/0908.3076),
   J. Reine Angew. Math. **672** (2012), 177--222.  Equations (3.7)--(3.10)
   give the finite Fourier-transform formula for the Weil \(S\)-action;
   Proposition 4.5 is the coefficient obstruction pairing; Theorem 6.8 is
   the integral meromorphic lift; and the proof of Theorem 7.1 contains the
   quoted Hilbert integral-basis assertion.
2. W. J. McGraw,
   [*The rationality of vector valued modular forms associated with the Weil
   representation*](https://doi.org/10.1007/s00208-003-0413-1),
   Math. Ann. **326** (2003), 105--122.  This proves the integral Fourier-basis
   theorem for the classical Weil representation over \(\mathbf Q\); it is
   not itself the omitted totally-real Hilbert adaptation.
3. M. Knopp, W. Kohnen, and W. Pribitkin,
   [*On the signs of Fourier coefficients of cusp forms*](https://doi.org/10.1023/A:1026207515396),
   Ramanujan J. **7** (2003), 269--277.  Their theorem applies to positive
   real weight and a unitary multiplier on a finitely generated Fuchsian group
   of the first kind, and permits the shifted Fourier expansion occurring in
   (5.8).
4. F. Strömberg,
   [*Weil representations associated to finite quadratic modules*](https://arxiv.org/abs/1108.0202),
   especially arXiv Lemma 5.13 (journal Lemma 5.6) and the preceding identification of the
   \(\Gamma_0(\ell)\) multiplier.
5. G. Chen, Y.-K. Lau, and Y. Wang,
   [*On the first sign change of Fourier coefficients of cusp forms*](https://preprint.press.jhu.edu/ajm/sites/default/files/AJM-chen-lau-wang.pdf),
   Theorem 1.  This gives a \((kN)^{2+\varepsilon}\) first-sign-change bound
   for integral or half-integral weight forms on \(\Gamma_0(N)\) with scalar
   nebentypus.  The orbit-square construction in Section 6 supplies the
   required scalar form.

## 9. Novelty statement

The invertibility of the full finite Fourier transform, the ordinary
q-expansion separation principle, and the Knopp--Kohnen--Pribitkin sign-change
theorem are standard.  The project contribution is their combination with
the determinant-line auxiliary discriminant factor: the proof that this
factor is automatically nontrivial, the compact deep-level congruence
argument that simultaneously isolates every selected facet on the closed
real chamber, and the diagonal/orbit-square arguments converting sign change
into positive-cone generation by slice-disjoint columns with the quantitative
field-degree split of Theorem 6.1.  No claim is made
that the underlying modular-form ingredients are new, or that a new polytope
or dimension bound follows from this lemma alone.
