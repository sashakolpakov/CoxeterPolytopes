# Orbit-square scalarization and a quantitative safe-height theorem

**Status (2026-08-09): proved, with an essential field-degree split.**
For the classical trace representation over \(F=\mathbf Q\), the
Chen--Lau--Wang first-sign-change theorem gives a polynomial bound for the
largest scaled Fourier index needed in the positive safe cone.  For a
general totally real field, the same scalar argument gives the same
polynomial sign bound after one also supplies a finite safe-detection cutoff;
such a cutoff always exists by the safe-span theorem and is found by exact
coefficient-rank enumeration, but no bound for it in terms of weight and
discriminant level alone is proved here.  This distinction is forced by the
possible kernel of diagonal restriction for Hilbert modular forms.

The theorem bounds compensator **trace height at an already chosen
discriminant level**.  It does not by itself bound the sheet-separating level
uniformly in a class of polytopes, the denominators of a rational conic
solution, the resulting Borcherds weight, or its global divisor degree.

## 1. Setup and exponent normalization

Let

\[
 V=S_{\kappa,\rho_A}
\]

be the Hilbert--Weil cusp obstruction space for the unglued direct
determinant extension over a totally real field \(F\) of degree \(d\).  The
parallel weight is \(\kappa_j=(n+2)/2\), and the diagonal weight is

\[
 K=\sum_{j=1}^d\kappa_j=\frac{d(n+2)}2.
 \tag{1.1}
\]

Here \(A=M'/M\) is also the discriminant module of the underlying even trace
\(\mathbf Z\)-lattice, with quadratic form
\(Q_{\mathbf Z}=\operatorname{Tr}_{F/\mathbf Q}Q\).  Let \(\ell\) be its
level, the least positive integer for which
\(\ell Q_{\mathbf Z}(x)\in\mathbf Z\) for every \(x\in A\).  Bruinier's
formulas for diagonal \(S\) and \(T=T_1\) are exactly the classical Weil
formulas for this finite quadratic module.

For

\[
 g=\sum_{x\in A}\sum_{m\gg0}b_g(m,x)
       e(\operatorname{Tr}(m\tau))\,\mathfrak e_x,
 \tag{1.2}
\]

put

\[
 r(m)=\ell\operatorname{Tr}_{F/\mathbf Q}(m)\in\mathbf Z_{>0}.
 \tag{1.3}
\]

For a fixed trace, only finitely many totally positive \(m\) occur: they lie
in a compact simplex in the Minkowski embedding of the relevant shifted
lattice.  Thus the diagonal component

\[
 g_x^\Delta(z)=g_x(z,\ldots,z)
\]

has well-defined real coefficients whenever \(g\) belongs to a real
q-structure.

Write \(A=A_0\oplus B\), where \(B\) is the auxiliary determinant-line
factor, and call \(x=(\alpha,\beta)\) safe when \(\beta\ne0\).

## 2. Why a single component does not descend to \(\Gamma_0\)

An individual nonisotropic component satisfies

\[
 g_x^\Delta(z+1)=e(Q_{\mathbf Z}(x))g_x^\Delta(z).
 \tag{2.1}
\]

Since every \(\Gamma_0(N)\) contains \(T\), while a Dirichlet nebentypus is
trivial on \(T\), (2.1) rules out direct scalar descent whenever
\(Q_{\mathbf Z}(x)\notin\mathbf Z\).  Passing to the Weil kernel only gives a
principal- or \(\Gamma_1\)-type subgroup.  This is the exact failure in the
naive scalarization.

## 3. The orbit-square scalarization

Let

\[
 \mathcal O_x=\{u x:u\in(\mathbf Z/\ell\mathbf Z)^\times\}
 \tag{3.1}
\]

and define

\[
 H_{\mathcal O_x}(z)=
   \sum_{y\in\mathcal O_x}\bigl(g_y^\Delta(\ell z)\bigr)^2.
 \tag{3.2}
\]

### Theorem 3.1 (controlled scalar cusp form)

For every \(g\in V\), the function (3.2) is a scalar holomorphic cusp form
of integral weight \(2K\) on

\[
 \Gamma_0(N_*),\qquad N_*=4\ell^2,
 \tag{3.3}
\]

with a Dirichlet nebentypus.  If the component q-expansions of \(g\) are
real, then the Fourier coefficients of \(H_{\mathcal O_x}\) are real.  The
unit orbit of a safe \(x\) consists entirely of safe elements.  If at least
one diagonal component \(g_y^\Delta\), \(y\in\mathcal O_x\), is nonzero and
has real coefficients, then \(H_{\mathcal O_x}\ne0\).

#### Proof

Let \(A_\ell=\operatorname{diag}(\ell,1)\).  For
\(\gamma=\left(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right)
\in\Gamma_0(4\ell^2)\),

\[
 \delta=A_\ell\gamma A_\ell^{-1}
 =\begin{pmatrix}a&\ell b\\c/\ell&d\end{pmatrix}
 \in\Gamma_0(\ell).
 \tag{3.4}
\]

Strömberg's arXiv Lemma 5.13 (journal Lemma 5.6) gives, for the Weil representation of an arbitrary
finite quadratic module of level \(\ell\),

\[
 \rho_A(\delta)\mathfrak e_y
 =e(\ell b d Q_{\mathbf Z}(y))\chi_A(\delta)
       \mathfrak e_{dy}
 =\chi_A(\delta)\mathfrak e_{dy}.
 \tag{3.5}
\]

The last equality is the definition of the level.  Scalar slashing commutes
with the rescaling because

\[
 \delta(\ell z)=\ell\gamma z,qquad
 (c/\ell)(\ell z)+d=cz+d.
 \tag{3.6}
\]

Multiplication by \(d\), a unit modulo \(\ell\), permutes
\(\mathcal O_x\).  Squaring and summing therefore makes (3.2) transform by
the common multiplier \(\chi_A(\delta)^2\).

Strömberg proves that on \(\Gamma_0(\ell)\), \(\chi_A\) is a Dirichlet
character in even signature, and is the theta multiplier (or its conjugate)
times a Dirichlet character in odd signature; in the odd case \(4\mid\ell\).
His formulas (5.5)--(5.6) give
\(v_\theta(\delta)^2=\left(\frac{-1}{d}\right)\).  Hence
\(\chi_A(\delta)^2\) depends only on \(d\) and is a Dirichlet character after
restriction to (3.3).  This also proves that the weight is the integral
weight \(2K\).

Diagonal restriction of a Hilbert cusp form is cuspidal on the finite-index
kernel of the diagonal Weil representation: every rational diagonal cusp is
a Hilbert cusp.  Rescaling conjugates that finite-index group, and (3.2) is
cuspidal on its intersection with \(\Gamma_0(4\ell^2)\), which is finite
index in \(\Gamma_0(4\ell^2)\).  Since (3.2) is modular on (3.3), vanishing at
every rational cusp for that finite-index intersection implies vanishing at
every cusp of (3.3), with only the cusp width changed.  Holomorphy is
immediate from restriction, rescaling, multiplication, and finite summation.

Finally, unit multiplication is an automorphism of \(A=A_0\oplus B\), so it
cannot send \(\beta\ne0\) to zero.  On the imaginary axis all component
values in a real q-structure are real.  A sum of their squares can vanish
identically only if every component in the orbit vanishes.  This proves the
last assertion. \(\square\)

## 4. Quantitative positive-cone bounds

Chen--Lau--Wang's Theorem 1 says that for every \(\varepsilon>0\), a nonzero
scalar cusp form of real Fourier coefficients, weight \(k\), level \(N\),
and Dirichlet nebentypus has two coefficients of opposite sign at indices

\[
 j\le C_\varepsilon(kN)^{2+\varepsilon}.
 \tag{4.1}
\]

Put

\[
 B_{\rm CLW}(K,\ell;\varepsilon)=
 \left\lceil C_\varepsilon(8K\ell^2)^{2+\varepsilon}\right\rceil,
 \tag{4.2}
\]

using \(k=2K\) and \(N=4\ell^2\).  The published theorem uses
\(\ll_\varepsilon\), so (4.2) is a controlled polynomial bound rather than a
printed numerical constant.

### Theorem 4.1 (closed polynomial bound over \(\mathbf Q\))

Assume \(F=\mathbf Q\) and a real q-structure on \(V\).  The positive cone
of the finite set of safe coefficient columns satisfying

\[
                  \ell m\le B_{\rm CLW}
 \tag{4.3}
\]

is all of \(V_{\mathbf R}^*\).  Consequently every real target obstruction
has a nonnegative safe correction using only indices (4.3), with at most
\(\dim V\) nonzero coefficients.  For a rational q-structure and rational
target the correction is rational; clearing denominators gives an integral
correction after multiplying the target.

#### Proof

Suppose a nonzero \(g\in V_{\mathbf R}\) were nonnegative on every column in
(4.3).  Safe linear spanning supplies a nonzero safe component \(g_x\).
For \(F=\mathbf Q\), diagonal restriction is the identity, so Theorem 3.1
gives a nonzero \(H_{\mathcal O_x}\).  Every coefficient of this square-sum
through index \(B_{\rm CLW}\) is a sum of products of safe coefficients whose
scaled indices are individually at most \(B_{\rm CLW}\); it is therefore
nonnegative.  This contradicts (4.1).  Hence the polar of the finite cone is
zero, and the finite-dimensional bipolar theorem gives the whole dual.
Conic Caratheodory and rational polyhedral linear algebra give the remaining
claims. \(\square\)

For a general totally real field, diagonal restriction need not be
injective.  Define the finite safe-detection cutoff

\[
 R_{\rm det}=min\left\{R:
 \operatorname{span}_{\mathbf R}
 \{\ell_{m,x}:x\text{ safe},\ r(m)\le R\}=V_{\mathbf R}^*
 \right\}.
 \tag{4.4}
\]

It is finite by the safe-span theorem.  Given exact q-expansions, it is found
by increasing \(R\) until the coefficient matrix has full rank.  This is a
per-space finite certificate, not a bound in \(K\) and \(\ell\).

### Theorem 4.2 (totally-real bound with the necessary detection term)

For arbitrary totally real \(F\), put

\[
 R_*=max\{R_{\rm det},B_{\rm CLW}(K,\ell;\varepsilon)\}.
 \tag{4.5}
\]

Then the safe coefficient columns with
\(\ell\operatorname{Tr}(m)\le R_*\) have positive cone
\(V_{\mathbf R}^*\), with the same support-cardinality and rational/integral
consequences as Theorem 4.1.

#### Proof

If a nonzero \(g\) were nonnegative on all columns through \(R_*\), the
spanning property in (4.4) would make at least one such coefficient strictly
positive.  All indices in the same safe component with the same trace are
also in the bounded set and are nonnegative, so the corresponding diagonal
trace coefficient is positive.  Thus some safe diagonal component, and
hence its orbit-square, is nonzero.  Every orbit-square coefficient through
\(B_{\rm CLW}\) is nonnegative, contradicting (4.1).  The cone argument is
then identical to Theorem 4.1. \(\square\)

The term \(R_{\rm det}\) cannot simply be deleted from this proof: nonzero
Hilbert cusp forms can vanish on the diagonal, and coefficients with equal
trace can cancel.  A uniform Hilbert Sturm/detection theorem for the present
finite Weil kernels would be precisely what replaces (4.4) by a formula in
field, weight, and level data.

## 5. Interaction with the sheet-separating level

The deep level itself is effective for each fixed polytope.  In a fixed
\(\mathbf Z\)-basis, enumerate the finite nonzero difference set used to
separate target facets from all same-norm sheets meeting the closed chamber,
and let \(C\) be the maximum absolute coordinate of those differences.  Then

\[
                   N_{\rm sep}=C+1
 \tag{5.1}
\]

keeps every nonzero difference out of \(N_{\rm sep}M\), with
\(N_{\rm sep}=1\) if the difference set is empty.  If \(\ell_0\) is
the level of the base trace discriminant module, the inherited form on the
sublattice \(N_{\rm sep}M\) has discriminant level \(\ell\) satisfying

\[
                   \ell\mid N_{\rm sep}^2\ell_0.
 \tag{5.2}
\]

Indeed \((N M)'=N^{-1}M'\), so
\(N^2\ell_0 Q(y/N)=\ell_0Q(y)\in\mathbf Z\).
The exact difference set and basis-coordinate proof are recorded in
[`explicit_sheet_separating_level_bound.md`](explicit_sheet_separating_level_bound.md).

Equations (4.2), (5.1), and (5.2) make the \(F=\mathbf Q\) compensator
trace-height bound fully computable for each fixed marked compact polytope.
They are not a uniform class bound: the finite mirror list and
\(N_{\rm sep}\) can grow with \(P\).  Over \([F:\mathbf Q]>1\),
\(R_{\rm det}\) is an additional uncontrolled term.  In every field, conic
denominators, the Borcherds weight, and global degree are not controlled by
the orbit-square cutoff alone.  The separate Cramer--Hadamard compression
shows that a bound for the finite integral obstruction-column height and
arithmetic volume controls those quantities; those inputs, together with
nonempty Heegner degrees when needed, still require class-uniform bounds before
the carrier can yield a new dimension or facet theorem.  See
[`borcherds_weight_degree_quantitative_bound.md`](borcherds_weight_degree_quantitative_bound.md).

## 6. Primary sources

1. J. H. Bruinier,
   [*Regularized theta lifts for orthogonal groups over totally real fields*](https://arxiv.org/abs/0908.3076),
   especially the trace-lattice and Weil formulas in Section 3 and the cusp
   space used in Proposition 4.5.
2. F. Strömberg,
   [*Weil representations associated to finite quadratic modules*](https://arxiv.org/abs/1108.0202),
   arXiv Lemma 5.13 (journal Lemma 5.6) for the
   \(\Gamma_0(\ell)\) monomial action and the discussion
   immediately before it identifying \(\chi_A\) as a Dirichlet character or
   theta multiplier times one.
3. G. Chen, Y.-K. Lau, and Y. Wang,
   [*On the first sign change of Fourier coefficients of cusp forms*](https://preprint.press.jhu.edu/ajm/sites/default/files/AJM-chen-lau-wang.pdf),
   Theorem 1.
