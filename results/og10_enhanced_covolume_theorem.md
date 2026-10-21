# Finite covolume from birational symmetries of the golden OG10 movable cone

> **Subsequent decision.**  The Bugaenko-facet reversal and the intrinsic H9
> Pell-fold involution both preserve the complete fixed chamber, and their
> product has trace `8+8*phi>10`.  The fixed-chamber stabilizer is therefore
> infinite, forcing both the H9 and H10 chambers to have infinite volume.
> The equivalence below then proves `Bir(X,tau)` infinite.  See
> [`og10_golden_h9_infinite_stabilizer.md`](og10_golden_h9_infinite_stabilizer.md).

## The theorem

Let `(X,eta,tau)` be the projective marked OG10 carrier constructed in
`og10_golden_h10_carrier.md`, re-marked as in
`og10_golden_complete_wall_chambers.md` so that its open eigen-movable cone
is a chamber `C` of the complete norm-two field-reflection arrangement in
the compatible hyperbolic space `H(V_phi^+)`.  Let `P=closure(C)` and define

\[
\operatorname{Bir}(X,\tau)=
\{f\in\operatorname{Bir}(X):f^*\tau=\tau f^*
  \text{ on }NS(X)\}.
\]

The commutation condition is numerical and uses the marked divisor-lattice
endomorphism.  It does not assert that `tau` is induced by a birational map
or extends to all of `H^2(X,Z)`.

Then

\[
\boxed{
 \operatorname{vol}_{hyp}(P)<\infty
 \iff P\text{ is compact}
 \iff |\operatorname{Bir}(X,\tau)|<\infty.}
\]

When these conditions hold, `P` has finitely many facets and is a compact
hyperbolic Coxeter polytope in dimension `10`.

For the present carrier this decision theorem is now evaluated negatively:
the complete chamber has infinite volume and is noncompact, its norm-two
reflection group has infinite arithmetic index, and the enhanced birational
centralizer is infinite.  Finite-sidedness remains undecided, and no general
H9 or H10 non-existence claim follows.

## 1. Arithmetic centralizer

Put

\[
F=\mathbf Q(\sqrt5),\qquad\mathcal O=\mathbf Z[\phi],
\]

and let `(M,H)` be the rank-11 field lattice

\[
H=\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}
  \oplus E_8\oplus\langle2\rangle.
\]

Its trace lattice is `(N,B)`, where

\[
B(x,y)=\operatorname{Tr}_{F/\mathbf Q}
\left(-\frac{\phi}{\sqrt5}H(x,y)\right),
\]

and `tau` is multiplication by `phi`.  Define the positive-component
centralizer

\[
\Gamma_\tau=\operatorname{Cent}_{O^+(N,B)}(\tau).
\]

The pair `(B,tau)` reconstructs `H`.  If `H(x,y)=a+b phi`, then

\[
B(x,y)=-(a+b),\qquad B(\tau x,y)=-(a+2b),
\]

and hence

\[
H(x,y)=
\bigl(B(\tau x,y)-2B(x,y)\bigr)
+\bigl(B(x,y)-B(\tau x,y)\bigr)\phi.
\]

Thus commuting with `tau` is equivalent to `O`-linearity, and a `B`-isometry
commuting with `tau` is an `H`-isometry.  Consequently

\[
\Gamma_\tau=O^+(M,H).
\]

Its action on `V_phi` is faithful: an `F`-linear transformation is determined
by either real embedding.

The arithmetic lattice `Gamma_tau` is uniform.  Indeed, `H` is anisotropic
over `F`.  If a nonzero `F`-vector were isotropic, its conjugate would be a
nonzero isotropic vector for the positive-definite conjugate form, which is
impossible.  Borel--Harish-Chandra and the anisotropy compactness criterion
therefore give

\[
\Gamma_\tau\backslash\mathbb H(V_\phi^+)\quad\text{compact}.
\]

This is the source of the important no-cusp consequence: any finite-covolume
norm-two reflection subgroup is automatically cocompact.

## 2. Reflection subgroup and chamber stabilizer

Let

\[
R_2=\{r\in M:H(r,r)=2\},\qquad
W_2=\langle s_r:r\in R_2\rangle.
\]

Every element of `Gamma_tau` permutes `R_2`, so

\[
W_2\triangleleft\Gamma_\tau.
\]

The complete-wall theorem identifies the mirrors of `W_2` with the complete
restriction of the square-`-2` stably prime-exceptional arrangement, hence with the
marked eigen-movable walls.  The arrangement is locally finite, and standard
reflection theory makes `W_2` simply transitive on its chambers.

For the selected chamber, put

\[
\Gamma_C=\{\gamma\in\Gamma_\tau:\gamma C=C\}.
\]

Every `Gamma_tau`-translate of `C` has a unique `W_2`-translate back to `C`.
Therefore

\[
\Gamma_\tau=W_2\rtimes\Gamma_C,
\qquad
\Gamma_\tau/W_2\cong\Gamma_C.
\]

Since `C` is a fundamental chamber for `W_2`, covolume multiplicity gives

\[
\operatorname{vol}(C)
=[\Gamma_\tau:W_2]\,
 \operatorname{vol}(\Gamma_\tau\backslash\mathbb H(V_\phi^+)),
\]

with infinite left side when the index is infinite.  Hence

\[
\operatorname{vol}(C)<\infty
\iff[\Gamma_\tau:W_2]<\infty
\iff|\Gamma_C|<\infty.
\]

Because `Gamma_tau` is uniform, finite index makes `W_2` uniform.  Local
finiteness of the mirror arrangement then forces the chamber closure `P` to
be compact.  Conversely compactness gives finite volume.  Thus

\[
\operatorname{vol}(C)<\infty
\iff P\text{ is compact}
\iff|\Gamma_C|<\infty.
\]

A compact chamber meets only finitely many locally finite mirrors, so it is
automatically finite-sided.

## 3. Passage to OG10 monodromy

The main scope issue is that `tau` lives only on `N`.  Define

\[
\Gamma_\tau^{mon}
=\{g|_N:g\in\operatorname{Mon}^2_{Hdg}(X),
      \ g(N)=N,\ g|_N\in\Gamma_\tau\}.
\]

This is a finite-index subgroup of `Gamma_tau`.  To see this, take the
discriminant kernel

\[
\Gamma_\tau^\#=
\ker(\Gamma_\tau\longrightarrow O(A_N)).
\]

It has finite index.  Since `N` is primitively embedded in the OG10
Beauville lattice, every element of this kernel extends by the identity on
the positive period complement
\(T=\langle2\rangle\oplus\langle66\rangle\) through the discriminant
gluing.
For the chosen positive component the extension preserves orientation and
fixes the period.  Onorati's theorem
`Mon^2(OG10)=O^+(Lambda_OG10)` makes it a Hodge monodromy operator.  Hence

\[
\Gamma_\tau^\#\subset\Gamma_\tau^{mon}
\subset\Gamma_\tau.
\]

The field reflection subgroup lies in this geometric finite-index subgroup.
For `H(r,r)=2`, put

\[
D=r,\qquad E=\tau^{-1}D.
\]

Then `D` and `E` are orthogonal square-`-2` stably prime-exceptional
divisorial wall classes, and

\[
s_r=(R_D R_E)|_N.
\]

The two Markman reflections are Hodge monodromy operators.  Therefore

\[
W_2\subset\Gamma_\tau^\#\subset\Gamma_\tau^{mon}.
\]

The first inclusion can also be seen directly on discriminant groups.  A
reflection in a square-`-2` vector `D` acts trivially on `A_N`, because for
`x in N^vee` one has `R_D(x)-x=(x,D)D in N`.  Hence the product `R_D R_E`,
and therefore every `s_r`, lies in the stable discriminant kernel.

Let

\[
(\Gamma_\tau^{mon})_C=Gamma_\tau^{mon}\cap\Gamma_C.
\]

It has finite index in `Gamma_C`, so it is finite exactly when `Gamma_C` is
finite.

## 4. The stabilizer is the enhanced birational group

The chosen marking satisfies

\[
C=FE_X\cap\mathbb H(V_\phi^+),
\]

where `FE_X=Mov^0(X)` is the fundamental exceptional chamber.  Take a
regular point `x in C`.  It lies on no ambient square-`-2` hyperplane: any
such hyperplane through `x` would meet `V_phi^+` and hence occur in the
complete restricted norm-two arrangement.

Nor can a nonzero ambient square-`-2` hyperplane contain all of `V_phi`.
Its rational normal would lie in the conjugate eigenspace, but the
irreducible polynomial `t^2-t-1` gives no nonzero rational eigenvector of
`tau`.

If `gamma in (Gamma_tau^mon)_C`, an ambient Hodge-monodromy extension sends
`FE_X` to another exceptional chamber containing `gamma x`.  But `gamma x`
also lies in `C subset FE_X` and is regular.  A regular point belongs to a
unique exceptional chamber, so the extension preserves `FE_X`.  Markman's
birational Torelli theorem then realizes it by a birational selfmap of `X`.
Its action on `N` commutes with `tau`.

Conversely, a member of `Bir(X,tau)` acts by Hodge monodromy, preserves the
movable cone, preserves `V_phi^+`, and hence stabilizes `C`.  Therefore

\[
\operatorname{im}\bigl(\operatorname{Bir}(X,\tau)\to O(N)\bigr)
=(\Gamma_\tau^{mon})_C.
\]

The kernel is finite.  An element acting trivially on `N` fixes an ample
class.  A birational selfmap fixing an ample class is regular, and therefore
belongs to the finite automorphism group of a polarized OG10 manifold.  It
follows that

\[
\operatorname{Bir}(X,\tau)
\quad\text{and}\quad\Gamma_C
\quad\text{are commensurable}.
\]

Combining this with the arithmetic calculation proves the boxed theorem.

There is also a compact quotient statement.  Normality and chamber
transitivity give

\[
\Gamma_\tau^{mon}
=W_2\rtimes(\Gamma_\tau^{mon})_C,
\]

and hence the orbifold identification

\[
\operatorname{im}(\operatorname{Bir}(X,\tau))\backslash P
\cong
\Gamma_\tau^{mon}\backslash\mathbb H(V_\phi^+).
\]

The right side is compact because `Gamma_tau^mon` has finite index in the
uniform arithmetic group `Gamma_tau`.  The left side is understood as the
closed Coxeter fundamental domain with its mirror-orbifold boundary.  The
quotient of the open chamber `C` is only the full-measure regular stratum; it
has the same finite measure, and its metric/orbifold completion is compact.
Thus the invariant measure always has finite mass after quotienting by
enhanced birational symmetries; the unquotiented chamber has finite mass
exactly when that symmetry group is finite.

Quantitatively,

\[
\operatorname{vol}(C)=
\left|\operatorname{im}
 (\operatorname{Bir}(X,\tau)\to O(N))\right|
\operatorname{vol}
 (\Gamma_\tau^{mon}\backslash\mathbb H(V_\phi^+)),
\]

where the right side is `+infinity` when the enhanced birational image is
infinite.

## 5. The invariant measure comes from divisor intersection algebra

The hyperbolic measure has a direct algebro-geometric formulation.  Since
`dim_C X=10`, the Fujiki relation is

\[
\operatorname{vol}_X(D)=c_X B(D,D)^5
\]

on big and nef classes.  On the positive movable cone the same formula
follows chamber-by-chamber on birational nef models, using birational
invariance of volume, and then on the walls by continuity.  Put

\[
\mathcal F(D)=-\frac1{10}\log\operatorname{vol}_X(D)
=-\frac12\log B(D,D)+\text{constant}.
\]

On the unit slice `B(D,D)=1`, tangent vectors satisfy `B(D,u)=0`, and direct
differentiation gives

\[
\operatorname{Hess}(\mathcal F)_D(u,v)=-B(u,v).
\]

This is exactly the hyperbolic metric on `H(V_phi^+)`.  On the
11-dimensional positive cone its Hessian Riemannian density is, up to a
constant,

\[
|\det(B|_{V_\phi})|^{1/2}B(D,D)^{-11/2}\,dD,
\]

and its radial/projective slice is hyperbolic volume.  The corresponding real
Monge--Ampere density is proportional to

\[
|\det(B|_{V_\phi})|B(D,D)^{-11}\,dD.
\]

Thus the invariant measure is encoded by the classical top-intersection
polynomial, or equivalently by asymptotic divisor-section growth.  What is
not encoded by the deformation-invariant Fujiki constant, Hodge numbers, or
Chern numbers is its **total mass on this chamber**: the constant cancels
from the Hessian, and finiteness is the global reflection-index/enhanced-
birational-group problem.

This answers the earlier intrinsic-invariant question sharply: the algebra
of the variety supplies the metric and measure, while the enhanced
birational symmetry group decides whether the selected Coxeter chamber has
finite total measure.

## 6. The enhancement is algebraic but not an ambient golden action

The divisor-lattice operator can be encoded by the algebra of the variety
itself.  Choose an ample divisor `A` and a rational divisor basis
`D_1,...,D_22` of `NS(X)_Q`.  Hodge--Riemann and the Fujiki relation make

\[
G_{ij}=\int_X D_iD_jA^8
\]

nondegenerate.  Let `Z_i in CH^9(X)_Q` be the dual linear combinations of
the complete-intersection cycles `D_j A^8`, so that

\[
\int_X D_k Z_i=\delta_{ki}.
\]

Then the divisor-supported correspondence

\[
\mathcal Z_\tau=
\sum_{i=1}^{22} Z_i\times\tau(D_i)
\in CH^{10}(X\times X)_{\mathbf Q}
\]

acts on `H^2(X,Q)` as `tau` on `NS(X)_Q` and as zero on the transcendental
period plane.  Thus `Bir(X,tau)` is the cohomological centralizer of an
explicit rational algebraic correspondence built from divisors and their
complete intersections.

Here the action convention is

\[
(\mathcal Z_\tau)_*(x)=
(p_2)_*\bigl(p_1^*x\cdot\mathcal Z_\tau\bigr).
\]

It kills the transcendental plane because, at maximal Picard rank,
\(T_{\mathbf C}=H^{2,0}\oplus H^{0,2}\), which pairs trivially with the
divisor complete intersections `D_j A^8`.

The cycle representative depends on `A` and the chosen basis, and no Chow-
level invariance is asserted.  Its cohomological operator is the canonical
projection to `NS_Q` followed by `tau`.

There is nevertheless an unavoidable ambient obstruction.  There is not
merely a lack of an ambient golden action: a self-adjoint one is impossible.

The OG10 Beauville lattice modulo `3` has a one-dimensional
radical.  Any self-adjoint endomorphism preserves that radical, but

\[
t^2-t-1
\]

is irreducible over `F_3`.  It cannot act on a nonzero one-dimensional
`F_3`-space while satisfying the golden polynomial.  Hence `tau` has no
integral self-adjoint extension to the full Beauville lattice.

The theorem avoids such an extension.  Only the finite-index discriminant
kernel and the paired reflection products are extended; these are genuine
Hodge monodromy operators.  An arbitrary non-self-adjoint integral extension
is neither required nor given a Beauville-form meaning.

Notice also that `Z_tau` does not evade the obstruction: it acts as zero on
the period plane, so the golden polynomial does not hold on all of `H^2`.

## 7. Primal--dual search consequence

There are now two exact routes to the volume decision:

1. **Primal:** complete the norm-two Vinberg chamber or exhibit a global
   hyperideal escape.
2. **Dual:** compute the centralizer of `tau` inside `Bir(X)`.

Proving `Bir(X,tau)` finite immediately makes the complete eigen-movable cone
a finite-sided compact Coxeter `10`-polytope.  Exhibiting one infinite-order
enhanced birational selfmap proves that the chamber has infinite volume and
is noncompact.

The current 171-wall and 321-wall calculations are bounded primal prefixes;
their surviving hyperideal rays do not decide the completed chamber.  The
dual group has not yet been computed.  The theorem is precisely the desired
primal--dual bridge between divisor geometry and hyperbolic covolume.

## Sources

The arithmetic inputs are Borel--Harish-Chandra's
[arithmetic-lattice theorem](https://doi.org/10.2307/1970210) and the
anisotropy compactness criterion.  Chamber normality and covolume use
Vinberg's [reflection-group theory](https://doi.org/10.1070/IM1971v005n05ABEH001203).
The geometric inputs are Markman's
[Corollary 5.7, Theorems 6.17 and 6.18(4)--(5), and Lemma
6.22](https://arxiv.org/abs/1101.4606), Mongardi--Onorati's
[Proposition 3.1](https://arxiv.org/abs/2010.12511), and Onorati's
[Theorem 5.4](https://arxiv.org/abs/2002.04129), which gives
`Mon^2(OG10)=O^+(Lambda)`.  Rapagnetta's
[OG10 Beauville--Fujiki calculation](https://arxiv.org/abs/math/0606409)
supplies the divisor-intersection volume polynomial.
