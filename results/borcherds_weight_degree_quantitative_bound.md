# Quantitative Borcherds weight and divisor degree after bounded safe support

**Status (2026-08-09): proved reduction and a closed fixed-carrier
rational-field formula once the finite obstruction matrix is enumerated;
the totally-real class-uniformity gap is isolated exactly.**  Bruinier's weight and
degree identities show that automorphic weight and top divisor degree are not
independent parameters once the integral principal part is known.  A
Cramer--Hadamard argument bounds the coefficients and denominator of a
positive safe correction from the height of one finite obstruction matrix.
For $F=\mathbf Q$, a direct Kloosterman estimate gives an explicit universal
bound for the required Eisenstein coefficients, and the
Hirzebruch--Mumford formula computes the remaining automorphic volume from
the lattice and subgroup index.

This does **not** yet give a new facet or dimension bound.  Bounded support
cardinality and index do not by themselves bound the conditioning of the
obstruction matrix.  In the compact totally-real branch one still needs a
uniform effective integral q-lattice, a bound for the safe-detection cutoff,
and a bound for the associated finite obstruction-column height.  Moreover,
an upper bound which grows linearly with the number of target facets is
circular when inserted into the facet inequality.

## 1. Setup and Bruinier's exact identities

Let $M$ be the even lattice in the determinant-line carrier, let
$X=\Gamma\backslash\mathcal D_M$, and put

\[
 \kappa=(n+2)/2,
 \qquad \lambda=\mathcal M_1,
 \qquad V_X=\int_X\Omega^n=\lambda^n.
\tag{1.1}
\]

For a totally positive Whittaker index $s=(m,\mu)$, write

\[
 E(s):=-B(m,\mu).
\tag{1.2}
\]

Bruinier's Eisenstein series is normalized as

\[
 E_M(\tau,\kappa)=\mathfrak e_0+
 \sum_{\mu}\sum_{m>0}B(m,\mu)q^m\mathfrak e_\mu.
\tag{1.3}
\]

His Remark 6.4 (Remark `rem:coeff` in the source) gives

\[
 B(m,\mu)=-\frac{\deg_\lambda Z(m,\mu)}{V_X},
 \qquad
 E(s)=\frac{\deg_\lambda Z(s)}{V_X}\ge 0.
\tag{1.4}
\]

The value is zero when the special divisor is empty.  If

\[
 f=\sum_s c_s f_s,
 \qquad c_s\in\mathbf Z_{\ge0},
\tag{1.5}
\]

is weakly holomorphic, Theorem 6.8 produces a Borcherds form with

\[
 \operatorname{div}(\Psi_f)=\sum_sc_sZ(s),
 \qquad
 w(\Psi_f)=-B(f)=\sum_sc_sE(s).
\tag{1.6}
\]

Consequently, on the compact neat carrier,

\[
 \boxed{
 w=\sum_sc_s\frac{\deg_\lambda Z(s)}{V_X},
 \qquad
 w\lambda^n=\sum_sc_s\deg_\lambda Z(s).
 }
\tag{1.7}
\]

Thus there is no additional unknown ``Borcherds weight formula'' after the
principal-part coefficients and special-cycle degrees are controlled.  Over
$\mathbf Q$, (1.7) is an identity for the interior cycles; a compactified
algebraic statement must additionally include the boundary order of the
product.

## 2. Effective denominator clearing is a minor bound

Assume that the obstruction space has an integral q-lattice
$V_{\mathbf Z}$ of rank $d$, and choose an integral basis
$g_1,\ldots,g_d$.  The obstruction column of an index $s$ is

\[
 b_s=(b_{g_1}(s),\ldots,b_{g_d}(s))^t\in\mathbf Z^d.
\tag{2.1}
\]

Let the nonnegative integral target principal part be

\[
 p_T=\sum_{u\in T}a_u f_u,
 \qquad
 C_T=\sum_{u\in T}a_u,
 \qquad
 b_T=\sum_{u\in T}a_ub_u.
\tag{2.2}
\]

Let $S$ be a finite safe set satisfying

\[
 -b_T\in\operatorname{cone}_{\mathbf R_{\ge0}}
 \bigl\{b_s:s\in S\bigr\}.
\tag{2.3}
\]

Put

\[
 R=\max_{s\in S}\lVert b_s\rVert_\infty,
 \qquad R_T=\lVert b_T\rVert_\infty.
\tag{2.4}
\]

### Theorem 2.1 (Cramer--Hadamard correction bound)

If $b_T\ne0$, there are linearly independent
$s_1,\ldots,s_t\in S$, with $t\le d$, and integers

\[
 q\in\mathbf Z_{>0},\qquad y_1,\ldots,y_t\in\mathbf Z_{\ge0}
\]

such that

\[
 qb_T+\sum_{i=1}^t y_i b_{s_i}=0,
\tag{2.5}
\]

and

\[
 \boxed{
 q\le t^{t/2}R^t,
 \qquad
 y_i\le t^{t/2}R^{t-1}R_T,
 \qquad
 \sum_i y_i\le t^{t/2+1}R^{t-1}R_T.
 }
\tag{2.6}
\]

If $b_T=0$, take $q=1,t=0$.

#### Proof

A minimal conic representation of $-b_T$ uses $t\le d$ linearly
independent safe columns.  Let $A\in\mathbf Z^{d\times t}$ be their
matrix and write $Ax=-b_T$, $x\ge0$.  Choose $t$ rows $J$ for
which $A_J$ is nonsingular.  Cramer's rule gives

\[
 x_i=\frac{\det A_J(i\leftarrow-b_{T,J})}{\det A_J}.
\]

Set $q=|\det A_J|$ and $y_i=qx_i$.  These are nonnegative integers.
Every column of $A_J$ has Euclidean norm at most $\sqrt tR$, while
$b_{T,J}$ has norm at most $\sqrt tR_T$.  Hadamard's determinant
inequality gives (2.6).  The full equation (2.5) follows because $A$ has
column rank $t$ and the selected-row solution is the original solution.
\(\square\)

The same statement holds for a rational q-structure once a common
denominator is supplied: multiply all obstruction columns by it before
using (2.6).  McGraw proves the required integral Fourier structure for the
classical Weil representation.  For general totally real $F$, Bruinier
states the analogous Hilbert integral-basis result as an adaptation of
McGraw but does not give its proof or an effective denominator.

### Corollary 2.2 (weight and degree from a finite matrix)

Assume $b_T\ne0$; the unobstructed case has $q=1,t=0$ and is bounded
directly by the target term.

Let

\[
 E_H=\max\{E(s):s\in T\cup S\}.
\tag{2.7}
\]

The effective principal part

\[
 qp_T+\sum_{i=1}^ty_i f_{s_i}
\tag{2.8}
\]

has weight bounded by

\[
 \boxed{
 w\le E_H\,t^{t/2}R^{t-1}(RC_T+tR_T)
 }
\tag{2.9}
\]

and, on the compact carrier,

\[
 \boxed{
 w\lambda^n\le
 E_HV_X\,t^{t/2}R^{t-1}(RC_T+tR_T).
 }
\tag{2.10}
\]

Equivalently, replace $E_HV_X$ by the largest special-cycle degree in
the finite support.  Equations (1.6), (2.6), and additivity prove the
corollary.

Suppose, in addition, that the exact-support normalization identifies each
marked facet with a distinct target divisor component and that its coefficient
is at least $a_{\min}=\min_{u\in T}a_u$.  Then the least marked-facet
multiplicity after clearing is at least $q a_{\min}$.  Since $q\ge1$, the
sheet-to-degree inequality yields

\[
 \#\operatorname{Facets}(P)
 \le \frac{E_HV_X}{a_{\min}}
 \left(C_T+t^{t/2+1}R^{t-1}R_T\right).
\tag{2.11}
\]

This is useful only when the quantities on the right are controlled
independently of the facet count.  In particular, the crude estimate
$C_T=\#\operatorname{Facets}(P)$ merely returns a circular inequality.

## 3. Why bounded support alone is not a coefficient bound

There is a sharp elementary obstruction to omitting $R$.  For any
$Q\ge1$, take the three integral safe columns

\[
 v_1=(1,0),\qquad v_2=(-Q,-1),\qquad v_3=(0,1)
\tag{3.1}
\]

and target $b_T=(0,1)$.  Their nonnegative cone is all of
$\mathbf R^2$, but every solution

\[
 x_1v_1+x_2v_2+x_3v_3=-b_T,qquad x_i\ge0,
\]

satisfies

\[
 x_1=Qx_2,\qquad x_2-x_3=1,
\]

so $x_1\ge Q$.  Dimension, support cardinality, and target size stay
fixed while the least correction mass tends to infinity.

This is a no-go for an inference from finite-dimensional cone theory alone,
not a construction of Weil cusp spaces with this behavior.  Arithmetic may
bound $R$; such a bound has to be proved or computed.  Once it is known,
Theorem 2.1 shows that coefficient denominators and minors are **not** an
additional uncontrolled gap.

## 4. An explicit Eisenstein-coefficient bound over $\mathbf Q$

Assume $F=\mathbf Q$, $n>2$, and
$\kappa=(n+2)/2>2$.  Schwagenscheidt's Proposition 2.2, reproducing
Bruinier's vector Eisenstein expansion, gives

\[
 c_0(\mu,m)=\frac{(2\pi)^\kappa m^{\kappa-1}}{\Gamma(\kappa)}
 \sum_{c\in\mathbf Z\setminus\{0\}}|c|^{1-\kappa}
 H_c^*(0,0,\mu,m).
\tag{4.1}
\]

The Weil representation is unitary.  Since $H_c^*$ is $1/|c|$ times
a sum over at most $\varphi(|c|)$ unit residues, every summand has

\[
 |H_c^*|\le \varphi(|c|)/|c|\le1.
\tag{4.2}
\]

Therefore

\[
 |c_0(\mu,m)|\le
 \frac{2(2\pi)^\kappa\zeta(\kappa-1)}{\Gamma(\kappa)}m^{\kappa-1}.
\tag{4.3}
\]

Bruinier's Eisenstein series (1.3) has constant coefficient
$\mathfrak e_0$, whereas the convention in (4.1) may have twice that
constant coefficient.  Thus the deliberately larger bound

\[
 \boxed{
 E(m,\mu)\le C_\kappa m^{\kappa-1},\qquad
 C_\kappa=\frac{2(2\pi)^\kappa\zeta(\kappa-1)}{\Gamma(\kappa)}
 }
\tag{4.4}
\]

is valid in either normalization.  If every target and safe index is at
most $H$, substitute $E_H\le C_\kappa H^{\kappa-1}$ in
(2.9)--(2.11).

Combining (4.4) with the proved orbit-square cutoff gives a fully finite
fixed-carrier rational-field expression once the integral obstruction matrix,
and hence $R$, is enumerated.  If $\ell$ is the trace-discriminant level and

\[
 B_{\rm CLW}=
 \left\lceil C_\varepsilon(8\kappa\ell^2)^{2+\varepsilon}\right\rceil,
\]

then the safe indices may be taken with
$\ell m\le B_{\rm CLW}$.  Take

\[
 H=\max\left\{\max_{u\in T}m_u, B_{\rm CLW}/\ell\right\}.
\tag{4.5}
\]

The published first-sign-change theorem writes its constant with
$\ll_\varepsilon$, so this is polynomial control rather than a printed
numerical constant.

### Compactness scope of the rational theorem

The closed bound above is not the high-dimensional compact branch of the
project.  By Meyer's theorem, an indefinite rational quadratic form in at
least five variables is isotropic.  An arithmetic hyperbolic lattice of
quadratic-form type is cocompact exactly when the defining rational form is
anisotropic.  Hence a compact arithmetic Coxeter polytope over
$F=\mathbf Q$ has hyperbolic dimension at most three.  In particular, the
rational theorem is principally a noncompact or low-dimensional
calibration, and its determinant type-IV quotient is noncompact in the
high-dimensional range as well.

The actual compact $\mathbb H^8$-and-higher route necessarily uses a totally
real field of degree greater than one.  There the quotient is anisotropic,
but the unresolved inputs are precisely the Hilbert integral q-lattice, the
detection cutoff $R_{\rm det}$, and the finite obstruction-column height
$R$ described in Section 6.

## 5. The top automorphic volume is arithmetic data

For a rational lattice $M$ of rank $r=n+2$, Gritsenko--Hulek--Sankaran
prove

\[
 \operatorname{vol}_{HM}(O(M))=
 \frac{2}{g_{sp}^+(M)}|\det M|^{(r+1)/2}
 \prod_{j=1}^{r}\pi^{-j/2}\Gamma(j/2)
 \prod_p\alpha_p(M)^{-1}.
\tag{5.1}
\]

Here $g_{sp}^+(M)$ is the number of proper spinor genera and the local
densities $\alpha_p(M)$ are explicit from the $p$-adic Jordan
decompositions; their paper records the odd-prime and dyadic formulas.
Their Hirzebruch--Mumford proportionality calculation also gives

\[
 \boxed{\lambda^n=2\operatorname{vol}_{HM}(\Gamma).}
\tag{5.2}
\]

The volume scales by the projective subgroup index.  Thus (5.1), the
orientation factor, and the exact index of the discriminant-kernel/neat
subgroup give $V_X$ in (2.10) without a new geometric unknown.  At the
direct level $M_N=NM$,

\[
 |\det M_N|=N^{2r}|\det M|,
\tag{5.3}
\]

and the local Jordan data and subgroup index are finite exact calculations.
The explicit sheet-separation theorem supplies $N=N_{\rm sep}(P)$ for an
individual polytope.  It does not give a class-uniform bound for that number.

For a general totally real field, Prasad's volume formula likewise expresses
the covolume of a principal arithmetic subgroup through the field
discriminant, global exponents, and local parahoric factors, with subgroup
index scaling.  That makes the volume arithmetic rather than mysterious,
but a uniform bound requires the relevant field, lattice, parahoric, and
index data to be bounded.  The rational formula (5.1) is not silently
asserted over a general $F$.

## 6. A general totally-real coefficient majorant

The exact identities (1.4)--(1.7) hold over every totally real field in
Bruinier's setup.  There is also a simple analytic finite-height bound which
does not require writing the local Whittaker factors explicitly.  Normalize
Haar measure on $\mathcal O_F\backslash\mathbf R^d$ to have mass one and
fix $v\in\mathbf R_{>0}^d$.  Fourier inversion gives

\[
 |B(m,\mu)|\le
 \exp(2\pi\operatorname{Tr}(mv))
 \sup_{u\in\mathcal O_F\backslash\mathbf R^d}
 |E_{M,\mu}(u+iv,\kappa)|.
\tag{6.1}
\]

For $n>2$, the defining adelic Eisenstein series is absolutely convergent
at $s_0=n/2>1$, so the supremum in (6.1) is finite and has a convergent
absolute-series majorant.  Hence trace height
$\operatorname{Tr}(m)\le H$, with $v=(1,\ldots,1)$, gives

\[
 E_H\le M_E(F,M,\kappa)e^{2\pi H},
\tag{6.2}
\]

where $M_E$ is the displayed compact-torus Eisenstein supremum (or its
absolute-series majorant).  This is a rigorous per-carrier bound, but not a
closed uniform formula in $\ell$ and $|A_M|$ alone.  An explicit Hilbert
local-density/Whittaker estimate, or a uniform bound for $M_E$ in a stated
field family, would strengthen (6.2).

The orbit-square theorem gives safe trace height

\[
 R_*=\max\{R_{\rm det},B_{\rm CLW}\}
\]

over general $F$.  The additional $R_{\rm det}$ is necessary because
diagonal restriction can have a kernel.  Even after $R_*$ is bounded,
Theorem 2.1 still requires the finite obstruction height $R$.  These are
the two exact quantitative inputs still absent from a uniform compact
totally-real theorem.

## 7. What is closed and what remains

The following implication is now rigorous:

```text
bounded exact-support level
  + bounded safe/target index set
  + effective integral (or denominator-bounded rational) cusp q-lattice
  + bounded finite obstruction-column height R
  + arithmetic volume/local-density bound
    => explicit Borcherds weight and top-degree bound.
```

In particular, the automorphic weight, coefficient denominator, and top
degree should no longer be listed as three unrelated missing lemmas.
Denominators are controlled by obstruction minors, weight is the Eisenstein
linear functional (1.6), and top degree is (1.7).  The remaining hard input
is quantitative control of the finite modular matrices and arithmetic
level/volume across an independently defined polytope class.

This report does **not** prove any of the following:

1. a uniform bound for $R_{\rm det}$ over totally real fields;
2. Bruinier's stated Hilbert integral-basis adaptation with effective
   denominators;
3. a uniform bound for the obstruction-column height $R$;
4. a compactified boundary-order bound over $\mathbf Q$;
5. a noncircular facet bound or a new dimension bound.

## 8. Primary sources

1. J. H. Bruinier,
   [*Regularized theta lifts for orthogonal groups over totally real
   fields*](https://arxiv.org/abs/0908.3076), J. Reine Angew. Math.
   **672** (2012), 177--222.  See the Eisenstein expansion in Section 6,
   Theorem 6.8, and Remark 6.4 for (1.3)--(1.7).
2. M. Schwagenscheidt,
   [*Eisenstein series for the Weil representation*](https://arxiv.org/abs/1803.10550),
   J. Number Theory **193** (2018), 74--90.  Proposition 2.2 gives
   (4.1); Theorem 1.4 gives the exact local representation-number formula.
3. V. Gritsenko, K. Hulek, and G. K. Sankaran,
   [*The Hirzebruch--Mumford volume for the orthogonal group and
   applications*](https://arxiv.org/abs/math/0512595), Doc. Math.
   **12** (2007), 215--241.  The main volume formula gives (5.1), and
   Proposition 1.2 gives (5.2).
4. G. Prasad,
   [*Volumes of S-arithmetic quotients of semi-simple groups*](https://doi.org/10.1007/BF02698841),
   Publ. Math. IHES **69** (1989), 91--117.
5. W. J. McGraw,
   [*The rationality of vector valued modular forms associated with the
   Weil representation*](https://doi.org/10.1007/s00208-003-0413-1),
   Math. Ann. **326** (2003), 105--122.
