# Finite covolume from automorphisms of a golden K3 pair

**Status.**  This is a scoped theorem for the even golden K3 carrier of
`results/k3_eigen_nef_factory_theorem.md`.  It is a direct synthesis of
arithmetic reduction theory, the universal shell theorem, K3 Weyl theory,
discriminant gluing, and global Torelli.  It is not asserted for an arbitrary
self-adjoint correspondence or arbitrary linear section of a nef cone.

The conclusion is the enhanced-pair analogue of the usual Sterk theorem:

\[
\boxed{
 \operatorname{vol}_{\mathrm{hyp}} P_{X,\tau}<\infty
 \iff P_{X,\tau}\text{ is compact}
 \iff |\operatorname{Aut}(X,\tau)|<\infty .}
\]

Here `Aut(X,tau)` means automorphisms whose action on the Neron--Severi
lattice commutes with the specified golden endomorphism.  No finite-generation
claim for an ``eigen-Cox ring'' is made.

## 1. Hypotheses and notation

Put

\[
 F=\mathbf Q(\sqrt5),\qquad
 \mathcal O=\mathbf Z[\phi],\qquad
 \phi^2=\phi+1,\qquad
 \alpha=-\frac{\phi}{\sqrt5}.
\]

Let `(M,h,j,p)` be an even golden K3 carrier datum of hyperbolic dimension
`d` in the sense of the eigen-nef factory theorem.  Thus:

1. `M` is a free `O`-module of rank `d+1` and `h` is an even symmetric
   `O`-valued form;
2. at the distinguished embedding of `F`, `h` has signature `(d,1)`, and at
   the conjugate embedding it is positive definite;
3. the trace form
   \[
   q(x,y)=\operatorname{Tr}_{F/\mathbf Q}(\alpha h(x,y))
   \]
   makes `N=(M,q)` an even lattice with a primitive embedding in the K3
   lattice; and
4. `p` is a regular time-like vector selecting a component `C` of the
   complete norm-two mirror arrangement.

Multiplication by `phi` is a self-adjoint integral endomorphism

\[
 \tau^2-\tau-I=0
\]

of `N`.  Let `E` be its distinguished real eigenspace.  The restriction of
`q` to `E` has signature `(1,d)`, and its positive unit hyperboloid is denoted
by `H(E)`.

Let

\[
 R_2=\{r\in M:h(r,r)=2\},\qquad
 W_F=\langle s_r:r\in R_2\rangle,
 \quad s_r(x)=x-h(x,r)r.
\]

### Complete-arrangement lemma

The universal shell theorem identifies the mirrors of `W_F` on `E` with the
restriction of the **complete** K3 `(-2)`-root arrangement.  Here is the
part of the shell argument needed below.

If `delta in N` has `q(delta,delta)=-2` and its hyperplane cuts the positive
cone of `E`, then

\[
 t=h(\delta,\delta)=a+b\phi
\]

is positive at both real embeddings.  The trace equation is

\[
 \operatorname{Tr}_{F/\mathbf Q}(\alpha t)=-(a+b)=-2.
\]

Total positivity leaves the five scalar possibilities

\[
 5-3\phi,\quad4-2\phi,\quad3-\phi,\quad2,\quad1+\phi.
\]

Evenness of `h` leaves only `2` and `4-2phi=2(phi-1)^2`.  Multiplication by
the unit `phi-1` identifies their hyperplanes and turns the latter into a
norm-two field root.  Conversely, every `r` with `h(r,r)=2` has trace square
`-2`.  Thus, in both directions,

\[
\boxed{
 \{E\cap\delta^\perp:\delta^2=-2,\ \delta^\perp\cap\mathbb H(E)\ne\varnothing\}
 =\{r^\perp\cap\mathbb H(E):h(r,r)=2\}.}                            \tag{1.1}
\]

On the trace lattice, `s_r` is the product of the two orthogonal K3 root
reflections in `r` and `(phi-1)r`.  In particular, `W_F` is an integral
subgroup of the K3 Weyl group and every mirror crossing is geometrically
legitimate.  Equality (1.1), rather than merely the lifting of a selected
finite wall list, is what makes the automorphism argument below possible.

After applying the Weyl framing of the factory theorem, take a projective
marked K3 surface `X` with `NS(X)=N`.  If `w` is the chosen K3 Weyl isometry,
all enhanced data are conjugated simultaneously:

\[
 \tau'=w\tau w^{-1},\quad E'=wE,\quad
 \Gamma'=w\Gamma_\tau w^{-1},\quad W_F'=wW_Fw^{-1},\quad C'=wC.
\]

This is a change of the enhanced marking; `w` is not asserted to be a surface
automorphism.  Retain the unprimed notation for the conjugated data.  The
factory theorem, using the complete arrangement equality (1.1), then gives
the equality of interiors and, separately, of closures:

\[
 C=\operatorname{Amp}(X)\cap\mathbb H(E),\qquad
 \overline C=P_{X,\tau}:=
 \operatorname{Nef}(X)\cap\mathbb H(E).
\]

Finally define the numerical automorphism group of the enhanced pair by

\[
 \operatorname{Aut}(X,\tau)=
 \{f\in\operatorname{Aut}(X):f^*\tau=\tau f^*
   \text{ on }\operatorname{NS}(X)\}.
\]

The endomorphism `tau` can be represented rationally by a divisor-supported
cycle on \(X\times X\), but the displayed group depends only on its numerical
action on `NS(X)`.  It is not claimed to fix a chosen Chow-cycle
representative.

## 2. The arithmetic centralizer is uniform

Let

\[
 \Gamma_\tau=
 \operatorname{Cent}_{O^+(N,q)}(\tau)
\]

acting on `H(E)`.

The pair `(q,tau)` reconstructs the field form `h`.  Explicitly, if
`h(x,y)=a+b phi`, then

\[
 q(x,y)=-(a+b),\qquad q(\tau x,y)=-(a+2b),
\]

so

\[
 h(x,y)=
 \bigl(q(\tau x,y)-2q(x,y)\bigr)
 +\bigl(q(x,y)-q(\tau x,y)\bigr)\phi.
\]

It follows that commuting with `tau` is the same as being `O`-linear, and a
`q`-isometry commuting with `tau` is precisely an `h`-isometry.  An
`O`-linear map is determined by its matrix at either real embedding, so the
action on `E` is faithful.  Consequently `Gamma_tau` is exactly the
positive-component arithmetic group

\[
 O^+(M,h).
\]

This arithmetic lattice is cocompact.  Indeed, `h` is anisotropic over `F`:
if a nonzero `v in M tensor F` satisfied `h(v,v)=0`, applying the conjugate
real embedding would produce a nonzero isotropic vector for the
positive-definite conjugate form.  This is impossible.  Arithmetic reduction
theory and the anisotropy compactness criterion therefore give

\[
 \Gamma_\tau\backslash\mathbb H(E)\quad\text{compact}.
\]

This is the decisive difference from the full rational Neron--Severi
hyperbolic space.  It is also why a finite-volume golden reflection chamber
cannot have cusps: the ambient arithmetic lattice is uniform.

## 3. Reflection subgroup and chamber stabilizer

Every element of `Gamma_tau` preserves `h` and hence permutes `R_2`.
Therefore

\[
 W_F\triangleleft\Gamma_\tau.
\]

The norm-two arrangement is locally finite.  Its chambers are the connected
components of the complement of its mirrors in `H(E)`.  Standard reflection
group theory gives a simply transitive action of `W_F` on these chambers.
For the selected chamber `C`, put

\[
 \Gamma_C=\{\gamma\in\Gamma_\tau:\gamma C=C\}.
\]

For every `gamma in Gamma_tau`, the chamber `gamma C` has a unique translate
back to `C` by `W_F`.  Since `W_F` is normal and has trivial chamber
stabilizer,

\[
 \boxed{\Gamma_\tau=W_F\rtimes\Gamma_C,
 \qquad \Gamma_\tau/W_F\cong\Gamma_C.}                 \tag{3.1}
\]

Since `C` is a fundamental chamber for `W_F`, covering-volume multiplicity
gives

\[
 \operatorname{vol}(C)
 =[\Gamma_\tau:W_F]\,
 \operatorname{vol}(\Gamma_\tau\backslash\mathbb H(E)),             \tag{3.2}
\]

with infinite left side when the index is infinite.  Combining (3.1),
(3.2), and cocompactness yields

\[
\begin{aligned}
 \operatorname{vol}(C)<\infty
 &\iff [\Gamma_\tau:W_F]<\infty\\
 &\iff |\Gamma_C|<\infty\\
 &\iff \overline C\text{ is compact}.
\end{aligned}                                                       \tag{3.3}
\]

For the last equivalence, finite index in a uniform lattice makes `W_F`
uniform; a reflection chamber for a uniform group has no ideal escape.  The
reverse implication is immediate from finite hyperbolic volume of a compact
polytope.  Local finiteness then also shows that a compact chamber has only
finitely many facets.

Thus finite sidedness, finite volume, and compactness do not need to be
proved by a global wall list if one can instead prove that `Gamma_C` is
finite.

## 4. The chamber stabilizer is an automorphism group up to finite index

Two points require care: preserving the restricted chamber must imply
preserving the **full** K3 ample chamber, and an isometry of `NS(X)` must pass
the Hodge/discriminant extension gate.

### 4.1 No hidden ambient-chamber ambiguity

Let `D=Amp(X)` be the open ambient K3 root chamber and take `x in C`.
Regularity and the shell theorem imply that `x` lies on no ambient K3 root
hyperplane.  In particular, no K3 root is orthogonal to all of `E`: such a
root hyperplane would contain `x`, and would therefore meet the positive cone
of `E`, contrary to regularity of the complete restricted arrangement.

If `gamma in Gamma_C`, then

\[
 \gamma x\in C\subset D.
\]

On the other hand, `gamma D` is an ambient K3 root chamber containing
`gamma x`.  A regular point belongs to a unique ambient root chamber, so

\[
 \boxed{\gamma D=D.}                                                \tag{4.1}
\]

Thus every restricted-chamber stabilizer automatically preserves the full
ample cone.  No extra assumption about roots whose mirrors miss `E` is
needed.

### 4.2 Discriminant kernel and Torelli

Let

\[
 \Gamma_\tau^\#=
 \ker\bigl(\Gamma_\tau\longrightarrow O(A_N)\bigr),
 \qquad A_N=N^\vee/N.
\]

The discriminant group is finite, so `Gamma_tau^#` has finite index in
`Gamma_tau`.  Because `N` is primitively embedded in the K3 lattice, every
element of `Gamma_tau^#` extends by the identity on \(T=N^\perp\) to an integral
isometry of the full K3 lattice: triviality on `A_N` preserves the
discriminant gluing.

For

\[
 \gamma\in\Gamma_\tau^\#\cap\Gamma_C,
\]

the extension is the identity on the period plane, preserves its orientation,
and preserves the ample cone by (4.1).  Global Torelli therefore realizes it
as the action of an automorphism of `X`.  It commutes with `tau` by
construction.  Hence

\[
 \Gamma_\tau^\#\cap\Gamma_C
 \subset \operatorname{im}\bigl(
 \operatorname{Aut}(X,\tau)\to O(N)\bigr)
 \subset \Gamma_C.                                                  \tag{4.2}
\]

The left inclusion has finite index in `Gamma_C`.  The kernel of the map
from `Aut(X,tau)` to `O(N)` is finite: it is contained in the automorphism
group of a fixed polarization, or equivalently in the finite group of Hodge
isometries acting trivially on `N`.  Thus (4.2) proves that

\[
 \boxed{\Gamma_C\text{ and }\operatorname{Aut}(X,\tau)
 \text{ are commensurable}.}                                       \tag{4.3}
\]

Combining (3.3) and (4.3) proves the theorem:

\[
\boxed{
 \operatorname{vol}_{\mathrm{hyp}}P_{X,\tau}<\infty
 \iff P_{X,\tau}\text{ is compact}
 \iff |\operatorname{Aut}(X,\tau)|<\infty.}                        \tag{4.4}
\]

Equivalently, the golden norm-two reflection subgroup has finite index in
the arithmetic centralizer exactly when the automorphism group of the
enhanced K3 pair is finite.

## 5. Algebro-geometric meaning and the Cox-ring boundary

The chamber remains encoded by divisor algebra.  If

\[
 \operatorname{VisEff}_\tau(X)=
 \overline{\operatorname{cone}}
 \{(D,\mathord\cdot)|_E:H^0(X,\mathcal O_X(D))\ne0\}\subset E^*,
\]

then

\[
 \operatorname{Nef}(X)\cap E
 =\overline{\mathcal C_E^+}\cap
   \operatorname{VisEff}_\tau(X)^\vee.
\]

The hyperbolic metric is the Hessian metric of divisor-section growth:
on a big and nef class, `vol_X(D)=D^2`, and on `D^2=1`,

\[
 \operatorname{Hess}\!\left(-\tfrac12\log\operatorname{vol}_X\right)
 (u,v)=-u\mathbin\cdot v.
\]

The theorem therefore says intrinsically:

> The space of unit-volume golden eigen-nef polarizations has finite Hessian
> volume exactly when the numerical automorphism group of the enhanced pair
> is finite; in that event the space is compact.

This does **not** make finite generation of the ordinary Cox ring equivalent
to (4.4).  For example, a Picard-rank-20 singular K3 has infinite full
automorphism group and a non-finitely-generated Cox ring, while its enhanced
centralizer `Aut(X,tau)` could still be finite.  Nor is a canonical algebra
graded by the irrational eigenspace being asserted.  The ordinary Cox ring
encodes the visible inequalities through its degree support, while the
centralizer theorem supplies the finite-covolume criterion.

## 6. Search consequence

For the discriminant-11 carrier, further wall traversal is no longer the only
route to a decision.  There are now exact dual certificates:

1. prove `Aut(X,tau)` finite; then the complete eigen-nef chamber is
   automatically a finite-sided compact hyperbolic Coxeter 9-polytope;
2. exhibit an infinite-order element of `Aut(X,tau)`; then the chamber has
   infinite volume and is noncompact.

The observed D68/D27 return supplies one involution in `Aut(X,tau)`.  It
does not decide whether the whole centralizer is finite.  The correct
algebro-geometric next task is therefore to compute the centralizer of `tau`
inside Shimada's presented automorphism group of the discriminant-11 K3,
in parallel with the root-wall completion.  This is precisely the
primal--dual strategy: the divisor walls approach from one side and enhanced
automorphisms from the other.

For Bugaenko's dimension-eight carrier, compactness of the archived chamber
and (4.4) imply that the corresponding enhanced automorphism group is finite.

## 7. Scope and dependencies

The proof uses all of the following and should not be detached from them:

- the totally real quadratic structure with exactly one indefinite real
  place and a definite conjugate place;
- the even golden shell theorem, which identifies **every** restricted K3
  root mirror with an integral norm-two field reflection;
- a regular eigenspace meeting an ambient ample chamber;
- a primitive embedding of the trace lattice in the K3 lattice; and
- global Torelli.

For a general self-adjoint endomorphism, the centralizer need not be an
arithmetic lattice of finite covolume, need not be uniform, and its visible
walls need not be reflection mirrors.  Formula (4.4) is not claimed in that
generality.

The arithmetic ingredients are the Borel--Harish-Chandra lattice theorem and
the anisotropy compactness criterion.  The chamber decomposition and normal
reflection subgroup use Vinberg's reflection-group lemmas.  The geometric
ingredients are the K3 Weyl-chamber theorem and global Torelli.  The golden
specialization, shell equality, trace reconstruction, and Weyl framing are
proved in `results/k3_eigen_nef_factory_theorem.md`; the full-nef comparison
is in `results/k3_nef_covolume_dictionary.md`.

## Primary sources

1. A. Borel and Harish-Chandra, *Arithmetic subgroups of algebraic groups*,
   Ann. of Math. 75 (1962), 485--535,
   [doi:10.2307/1970210](https://doi.org/10.2307/1970210).
2. E. B. Vinberg, *Discrete linear groups generated by reflections*, Math.
   USSR-Izv. 5 (1971), 1083--1119,
   [doi:10.1070/IM1971v005n05ABEH001203](https://doi.org/10.1070/IM1971v005n05ABEH001203).
3. I. I. Piatetski-Shapiro and I. R. Shafarevich, *A Torelli theorem for
   algebraic surfaces of type K3*, Math. USSR-Izv. 5 (1971), 547--588,
   [doi:10.1070/IM1971v005n03ABEH001075](https://doi.org/10.1070/IM1971v005n03ABEH001075).
