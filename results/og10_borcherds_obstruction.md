# The Borcherds route for the golden OG10 carrier

## Result and scope

There are two different automorphic constructions one can attach to the
golden OG10 carrier.  They have sharply different outcomes.

1.  The ordinary rational type-IV extension of the Neron--Severi lattice
    admits **no reflective Borcherds product**.  This is an immediate exact
    application of Wang's classification in signature `(22,2)`.
2.  The real-quadratic form can be enlarged by one sign-changing
    `Z[phi]`-line so that Bruinier's totally-real Whittaker lift applies.
    For this enlargement the single-special-divisor principal part has a
    concrete cusp-form obstruction in parallel weight `(6,6)`, and every
    fixed finite compensator ansatz is an exact finite-dimensional kernel
    calculation.  The pure norm-two Poincare series is nonzero, so the
    single-term principal part is obstructed.  Moreover, complete
    discriminant enumeration finds two off-slice unit orbits at the global
    smallest index norm `9/880`, and the target plus both orbits has an
    invertible `3x3` coefficient matrix.  Thus the first complete compensator
    shell is also obstructed; larger shells remain open.

Neither statement decides whether the norm-two reflection group in
`H^10` has finite index.  These two non-implications are essential:

> **No ordinary reflective product does not imply nonreflectivity.**

> **A totally-real Whittaker lift with the desired divisor would not
> imply finite covolume.**

Thus the ordinary Borcherds route is closed, while the totally-real route is
a rigorous automorphic-divisor test rather than a replacement for Vinberg's
algorithm or the enhanced-birational-group criterion.

## 1. The ordinary type-IV extension

The OG10 carrier has trace lattice

\[
N=U\oplus T_{11}(-1)\oplus2E_8(-1)\oplus2A_1(-1),
\qquad
T_{11}=\begin{pmatrix}2&1\\1&6\end{pmatrix}.
\]

It has signature `(1,21)` and determinant `-44`.  Use Wang's convention,
in which an orthogonal modular lattice has signature `(l,2)`, and put

\[
M_0=U\oplus N(-1)
   \cong2U\oplus T_{11}\oplus2E_8\oplus2A_1.
\]

Then

\[
\operatorname{sig}(M_0)=(22,2),\qquad
\det(M_0)=44.
\]

The discriminant calculation is also explicit.  The vector

\[
g=T_{11}^{-1}e_2=\frac{-e_1+2e_2}{11}
\]

generates `A_T11` and has quadratic value `2/11` modulo `2Z`.  A half-root
generates each `A_A1` and has value `1/2`.  Therefore

\[
A_{M_0}\cong\mathbf Z/11\oplus(\mathbf Z/2)^2,
\qquad
q_{M_0}\cong
\left\langle\frac2{11}\right\rangle
\oplus\left\langle\frac12\right\rangle^2.
\]

Its Smith invariant factors are `(2,22)`.  The level is `44`: multiplying
`2/11` by `11` makes it even, while multiplying `1/2` by `4` is necessary
and sufficient to make it even.

For a Borcherds product on a lattice of signature `(l,2)`, the vector-valued
input has weight

\[
1-\frac l2=-10,
\]

and the dual obstruction space has weight

\[
1+\frac l2=12.
\]

The term `q^-1 e_0` is the principal part attached to norm-two vectors.
In particular, this is a level-44, 44-component vector-valued modular-form
problem, not a scalar level-11 problem.

### Classification obstruction

Wang's Theorem 1.2(2), restated and proved in detail as Theorem 3.5,
states that if an even lattice

\[
M=U\oplus K,
\qquad \operatorname{sig}(M)=(22,2),
\]

has a reflective Borcherds product, then

\[
M\cong2U\oplus D_{20}.
\]

The latter lattice has determinant `4`, whereas `M_0` has determinant
`44`.  Hence

\[
\boxed{M_0\text{ has no reflective Borcherds product}.}
\]

Wang's earlier complete classification of 2-reflective modular forms gives
the narrower conclusion independently: in signature `(2,n)` with `n>=14`,
the only 2-reflective lattices are `II_(2,18)`,
`2U + 2E8(-1) + A1(-1)`, and `II_(2,26)`.  Our signature `(2,22)` lattice
is not on the list, so it has no 2-reflective modular form.

The conclusion is robust under passing to a finite-index arithmetic
subgroup.  A reflective form can be symmetrized over cosets, and because
`M_0` splits `2U`, the resulting form for the stable orthogonal group is a
Borcherds product by the converse result used in Wang's paper.  For the
present project, however, the directly needed and completely sufficient
statement is already the nonexistence of a reflective Borcherds product.

## 2. Why this does not prove infinite index

Let

\[
\Gamma_\tau=O^+(H,\mathbf Z[\phi]),\qquad
W_2=\langle s_r:H(r,r)=2\rangle.
\]

The group `Gamma_tau` acts on real hyperbolic `10`-space.  It is not a
finite-index arithmetic subgroup of `O(M_0)`, whose Hermitian type-IV
domain has complex dimension `22`.  A modular form for `O(M_0)` therefore
belongs to a different ambient automorphic problem.

More fundamentally, there is no converse from hyperbolic reflectivity to
an automorphic correction.  Borcherds states this explicitly in the
introduction to *Reflection groups of Lorentzian lattices* and asks for a
converse in Problem 13.2.  Wang makes the failure concrete: his Theorem 7.2
lists 18 hyperbolically 2-reflective lattices of rank greater than five with
no automorphic correction; for example

\[
S=U\oplus E_8\oplus E_7
\]

is hyperbolically 2-reflective although `U + S` is not 2-reflective modular.
Consequently

\[
\text{no reflective product on }M_0
\quad\not\Longrightarrow\quad
[\Gamma_\tau:W_2]=\infty.
\]

This is a theorem-level obstruction to a method, not a nonreflectivity
theorem for the OG10 chamber.

## 3. Why `H` itself is not a Borcherds domain

At the two real embeddings of `F=Q(sqrt(5))`, the field form

\[
H=\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}
  \oplus E_8\oplus\langle2\rangle
\]

has signatures

\[
((10,1),(11,0)).
\]

Thus

\[
\operatorname{Res}_{F/\mathbf Q}O(H)(\mathbf R)
\cong O(10,1)\times O(11),
\]

and its symmetric space is real hyperbolic `10`-space.  For dimension
greater than two this is not Hermitian.  Classical Borcherds products and
Bruinier's totally-real extension require signature `(n,2)` at one real
place and a definite form at every other real place.  They therefore do not
apply directly to `H`.

## 4. A minimal totally-real type-IV enlargement

There is nevertheless a natural one-line repair.  Put

\[
\widehat H=H\oplus\langle2(1-\phi)\rangle
\]

over `O=Z[phi]`.  Since `phi'=1-phi`, the added coefficient has signs

\[
2(1-\phi)<0,\qquad
2(1-\phi')=2\phi>0.
\]

Therefore

\[
\operatorname{sig}(\widehat H)=((10,2),(12,0)),
\]

exactly Bruinier's totally-real orthogonal-Shimura signature with `n=10`.

Let `sqrt(5)=2phi-1`, let `partial=(sqrt(5))` be the different, and define
the quadratic form

\[
\widehat Q(x)=\frac{\phi}{2\sqrt5}\widehat H(x,x).
\]

The scalar `phi/sqrt(5)` is totally positive.  Since `widehat H` is even,
`widehat Q(O^12)` lies in `partial^-1`, which is precisely Bruinier's
integrality convention.

On the `H` summand the associated integral trace bilinear form is `-B`, so
its trace lattice is `N(-1)`.  If `e` generates the new line, then

\[
\widehat Q(e)
=\frac{\phi(1-\phi)}{\sqrt5}
=-\frac1{\sqrt5}.
\]

In the integral basis `(1,phi)` of that line, its trace bilinear Gram matrix
is

\[
\left(\operatorname{Tr}\left(-\frac{2e_i e_j}{\sqrt5}\right)\right)
=\begin{pmatrix}0&-2\\-2&-2\end{pmatrix}
\cong\langle2\rangle\oplus\langle-2\rangle.
\]

Thus the rank-24 trace lattice of the totally-real extension is

\[
M_F\cong N(-1)\oplus\langle2\rangle\oplus\langle-2\rangle,
\]

with

\[
\operatorname{sig}(M_F)=(22,2),\qquad \det(M_F)=176,
\]

and

\[
\begin{aligned}
A_{M_F}&\cong\mathbf Z/11\oplus(\mathbf Z/2)^4,\\
q_{M_F}&\cong
\left\langle\frac2{11}\right\rangle
\oplus\left\langle\frac12\right\rangle^3
\oplus\left\langle-\frac12\right\rangle.
\end{aligned}
\]

The Smith invariant factors are `(2,2,2,22)`, the group algebra of the
Weil representation has dimension `176`, and the underlying finite
quadratic module again has level `44`.

The ordinary rational group `O(M_F)` is not the totally-real group in
Bruinier's construction.  Wang's theorem also rules out a reflective
Borcherds product for the full rational orthogonal group of `M_F`, but it
does not rule out a Whittaker lift for the proper subgroup

\[
\operatorname{Res}_{F/\mathbf Q}\operatorname{GSpin}(\widehat H).
\]

## 5. The exact Whittaker obstruction test

For `n=10` and `[F:Q]=2`, Bruinier's input weight and dual cusp weight are

\[
k=\left(\frac{2-n}{2},\frac{n+2}{2}\right)=(-4,6),
\qquad
\kappa=(6,6).
\]

Let `L_F=(O^12,widehat Q)` and let `rho_LF` be Bruinier's Hilbert--Weil
representation of this `O`-lattice on `C[L_F'/L_F]=C[A_MF]`.  (Its
underlying vector space has dimension `176`, but its
`SL_2(O)`-action also uses the `O`-module structure; it is not merely the
classical Weil representation of the rank-24 trace lattice.)  The
obstruction space is the finite-dimensional vector-valued Hilbert cusp-form
space

\[
S_{(6,6),\rho_{L_F}}.
\]

For a norm-two field root `r`,

\[
m_0=\widehat Q(r)=\frac{\phi}{\sqrt5},
\qquad \mu=0.
\]

This is totally positive.  Let `f_(m0,0)` denote Bruinier's corresponding
harmonic Whittaker basis element.  If

\[
g=\sum_{m,\mu}b_g(m,\mu)q^m\mathfrak e_\mu
\in S_{(6,6),\rho_{L_F}},
\]

then Bruinier's Proposition 4.5 gives the exact criterion

\[
\boxed{
f_{m_0,0}\text{ is weakly holomorphic}
\iff
b_g(m_0,0)=0\quad\text{for every }g
\in S_{(6,6),\rho_{L_F}}.}
\]

This coefficient functional has now been evaluated without constructing a
full cusp-form basis.  The normalized diagonal coefficient of Bruinier's
Poincare series

\[
P_{m_0,0}=\xi_{(-4,6)}(f_{m_0,0})
\]

has identity term `1`; exact rational Kloosterman--Bessel bounds give

\[
\left|b_{P_{m_0,0}}(m_0,0)-1\right|<0.931739<1.
\]

Therefore `P_(m0,0)` is nonzero and there exists
`g in S_(6,6),rho_LF` with `b_g(m0,0) != 0`.  In particular the
single-term Whittaker principal part is obstructed.  See the
[exact nonvanishing theorem](og10_whittaker_nonvanishing.md), its
[`Fraction` certificate](../scripts/certify_og10_whittaker_nonvanishing.py),
and [test](../tests/test_og10_whittaker_nonvanishing.py).

Had the functional been zero, the coefficient of the Whittaker principal part would be the
integer `1`, so Bruinier's Theorem 6.8 (stated as Theorem 1.3 in the
introduction) produces a meromorphic modular form with divisor
`Z(m0,0)`.  The proved nonzero functional means that the **single-term
principal part** does not lift.  This alone does not exclude a different principal part with the same
geometric divisor: one may have to add terms whose special divisors are
empty or cancel.

More generally, a prescribed finite principal part

\[
f=\sum_s c_s f_{m_s,\mu_s}
\]

passes exactly when

\[
\sum_s c_s b_g(m_s,\mu_s)=0
\]

for every cusp form `g`.  Hence any finite compensator ansatz is an exact
matrix-kernel calculation, not a heuristic search.  Without a bound on the
allowed compensators, no finite calculation is a global nonexistence theorem
for all meromorphic forms with the same divisor.  The first complete finite
calculation is now available: among all classes whose auxiliary projection
is nonzero, the smallest totally positive index norm is `9/880`, there are
two divisor-relevant unit orbits at that norm, and the resulting target plus
two-compensator `3x3` matrix is strictly diagonally dominant after an exact
similarity scaling.  Hence it has no kernel.  See the
[first-compensator theorem](og10_whittaker_first_compensator.md),
[certificate](../scripts/certify_og10_whittaker_first_compensator.py), and
[test](../tests/test_og10_whittaker_first_compensator.py).

## 6. Why the special divisor restricts to the correct mirrors

Work on the universal cover of the lattice component represented by
`h=1`.  At the identity embedding, let `e` be the negative generator of the
added line.  A negative line `R x` in `H` determines the oriented negative
2-plane

\[
\iota([x])=\operatorname{span}_{\mathbf R}(e,x)
\]

in the type-IV domain of `widehat H`.  This embeds hyperbolic `10`-space as
a totally real, totally geodesic submanifold.

For `r in H`, the special divisor `r^perp` intersects this real slice in

\[
\{[x]:(x,r)=0\},
\]

the required hyperbolic reflection mirror.  Conversely, write a vector of
`widehat H` as `lambda=r+a e`.  If `a` is nonzero, then

\[
(e,\lambda)=a(e,e)\ne0,
\]

so `lambda^perp` misses the real slice.  Consequently, on this distinguished
component,

\[
\iota(\mathbb H^{10})\cap \operatorname{supp} Z(m_0,0)
=\bigcup_{H(r,r)=2}\iota(r^\perp).
\]

The restriction is essential: globally `Z(m0,0)` also has components from
vectors with nonzero `e`-coordinate.  For example, if `u` generates the
`<2>` summand of `H`, then `lambda=phi*u+phi*e` satisfies
`widehat H(lambda,lambda)=2`, but its divisor misses the real slice by the
preceding orthogonality test.

The single Whittaker coefficient therefore targets exactly the complete
norm-two mirror arrangement on the OG10 movable slice.

## 7. Why passing the test would still not prove finite covolume

Because the conjugate field form is definite, the orthogonal Shimura
variety for `widehat H` is compact.  Bruinier emphasizes that for degree
greater than one there are no cusps and no Borcherds product expansions at
cusps.  In particular the lift supplies no Weyl vector or denominator
expansion from which a fundamental chamber could be read.

The restriction of the resulting automorphic section would encode the
mirror divisor, but the same divisor can coexist with either a finite or an
infinite chamber stabilizer.  It therefore gives neither direction of the
desired index decision:

\[
\text{Whittaker lift exists}
\quad\not\Longrightarrow\quad
[\Gamma_\tau:W_2]<\infty,
\]

and failure of its obstruction test only rules out that automorphic
realization, not finite index.

For the present OG10 carrier the valid decision theorem remains

\[
[\Gamma_\tau:W_2]<\infty
\iff |\operatorname{Bir}(X,\tau)|<\infty
\iff \overline C\text{ is compact},
\]

as proved in `og10_enhanced_covolume_theorem.md`.

## Primary sources

- R. E. Borcherds,
  [*Reflection groups of Lorentzian lattices*](https://arxiv.org/abs/math/9909123),
  especially the introduction, Lemma 12.1, and Problems 13.1--13.2.
- H. Wang,
  [*The classification of 2-reflective modular forms*](https://arxiv.org/abs/1906.10459),
  Theorem 1.1 and Theorem 7.2.
- H. Wang,
  [*On the classification of reflective modular forms*](https://arxiv.org/abs/2301.12606),
  Theorem 1.2(2) and the detailed Theorem 3.5.
- J. H. Bruinier,
  [*Regularized theta lifts for orthogonal groups over totally real fields*](https://arxiv.org/abs/0908.3076),
  Proposition 4.5 and Theorem 6.8 (Introduction, Theorem 1.3).
- M. Zhang and Y. Zhang,
  [*On the non-vanishing of Hilbert Poincare series*](https://doi.org/10.1016/j.jnt.2025.11.003),
  Sections 2--3 and Theorem 3.2.

## Exact verification

The convention-sensitive lattice calculations are independently replayed by

```bash
PYTHONPATH=.:src python3 scripts/certify_og10_borcherds_obstruction.py
PYTHONPATH=.:src python3 scripts/certify_og10_whittaker_nonvanishing.py
PYTHONPATH=.:src pytest -q tests/test_og10_borcherds_obstruction.py tests/test_og10_whittaker_nonvanishing.py
```

The external classification and theta-lift theorems are cited inputs, not
computer-certified claims.
