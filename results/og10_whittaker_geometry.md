# A geometric decision test for the pure OG10 special divisor

## Outcome

Let

\[
F=\mathbf Q(\sqrt5),\qquad \mathcal O=\mathbf Z[\phi],\qquad
L=(\mathcal O^{12},\widehat Q),
\]

where

\[
\widehat H=H\oplus\langle 2(1-\phi)\rangle,
\qquad
\widehat Q(x)=\frac{\phi}{2\sqrt5}\widehat H(x,x).
\]

At the two real places this has signatures `((10,2),(12,0))`.  Let
`X=X_K` be a neat compact orthogonal Shimura variety for `L`, let

\[
\ell=c_1(\mathcal M_1)\in \operatorname{Pic}(X)_\mathbf Q
\]

be the ample automorphic class, and put

\[
m_0=\frac{\phi}{\sqrt5},\qquad D=Z(m_0,0).
\]

There are two different questions, and they must not be identified.

1. **Pure Whittaker question.**  Is Bruinier's harmonic Whittaker basis
   element `f_0=f_(m_0,0)` weakly holomorphic?
2. **Scalar-section question.**  Is there a scalar automorphic section, of
   rational weight and finite-order multiplier, whose divisor is exactly
   `D`?

Write

\[
g_0=\xi_k(f_0)\in S_{(6,6),\rho_L}.
\]

Bruinier's formulas give the exact distinction

\[
\boxed{f_0\text{ is weakly holomorphic}\iff g_0=0,}
\]

whereas the Kudla--Millson lift `Lambda` gives

\[
\boxed{
\text{a scalar section with divisor }D\text{ exists}
\iff \Lambda(g_0)=0.
}
\]

The first question is now decided by an exact vector-valued Poincare-series
estimate:

\[
\boxed{g_0\ne0.}
\]

Hence the pure one-term Whittaker input is obstructed.  The second question
remains open: nonzero `g_0` does not exclude `g_0` from the kernel of the
Kudla--Millson lift.

Consequently the converse from a divisor relation to a Whittaker lift needs
injectivity of the Kudla--Millson lift.  It is not a formal consequence of
Bruinier's Chow-valued modularity theorem.

There is, however, an exact one-number geometric decision test for the
second question.  Set

\[
a=\frac{D\ell^9}{\ell^{10}}=-B(m_0,0)>0
\]

and

\[
\Delta_D
=D^2\ell^8-\frac{(D\ell^9)^2}{\ell^{10}}.
\]

Then

\[
\boxed{
\Delta_D\le 0,
\qquad
\Delta_D=0
\iff [D]=a\ell\text{ in }\operatorname{Pic}(X)_\mathbf Q
\iff \text{the scalar section exists}.}
\]

Thus `Delta_D<0` proves nonexistence of a pure scalar automorphic section,
without constructing a basis of the 176-component cusp space.  The pure
Whittaker obstruction proves `g_0 != 0`, but the value of `Delta_D` has **not**
yet been computed for this lattice.

## 1. Whittaker vanishing versus Kudla--Millson vanishing

Bruinier defines

\[
\xi_k:\mathcal H_{k,\bar\rho_L}\longrightarrow
S_{(6,6),\rho_L}.
\]

His equations (4.17)--(4.18) show that `g_0=xi_k(f_0)` is, up to the
displayed nonzero normalizing constant, the vector-valued Hilbert Poincare
series with seed

\[
e(\operatorname{tr}(m_0\tau))\,\mathfrak e_0.
\]

Proposition 4.3 and Definition 4.4 give

\[
f_0\in M^!_{k,\bar\rho_L}\iff g_0=0.
\]

Equivalently, by Proposition 4.5,

\[
b_g(m_0,0)=0
\quad\text{for every }g\in S_{(6,6),\rho_L}.
\]

On the orthogonal side, Corollary 5.16 and Theorem 6.4 give, in cohomology,

\[
[D]=[\Lambda(g_0)]-B(m_0,0)\ell.
\]

Bruinier's Remark 6.5 says

\[
B(m_0,0)=-\frac{D\ell^9}{\ell^{10}}=-a,
\]

so

\[
\boxed{[D]-a\ell=[\Lambda(g_0)].}
\]

The form `Lambda(g_0)` is harmonic.  Hence its cohomology class vanishes if
and only if the form itself vanishes.  This proves

\[
[D]=a\ell\iff \Lambda(g_0)=0.
\]

It does **not** prove `g_0=0`: a nonzero `g_0` in the kernel of the
Kudla--Millson lift would produce the totally-real analogue of an exceptional
automorphic product.

The available blanket injectivity theorems do not cover the present case.
The field quadratic space is anisotropic because it is positive definite at
the conjugate real place, so it cannot split a hyperbolic plane over `F`.
Moreover `L` is not unimodular: its trace discriminant has order `176`.
Thus neither the even-unimodular theorem of Bruinier--Funke nor the
split-`U+U(N)` theorem applies.  A Rallis inner-product proof here would still
have to establish nonvanishing of the bad local zeta integrals.

## 2. The scalar-section/Picard equivalence

Work first at neat level, so `X` is smooth and projective.  The same statement
can be read on the Shimura stack at non-neat level, or obtained after a normal
neat cover and taking the norm over the finite deck group.

Suppose `Psi` is a scalar meromorphic automorphic form of rational weight `w`
and finite-order multiplier with

\[
\operatorname{div}(\Psi)=D.
\]

After taking a positive power that kills the multiplier and clears the
denominator of `w`, `Psi` becomes a rational section of an integral power of
`M_1`.  Therefore

\[
[D]=w\ell\quad\text{in }\operatorname{Pic}(X)_\mathbf Q.
\]

Intersecting with `ell^9` forces `w=a`.

Conversely, if `[D]=a ell` in `Pic(X)_Q`, then for some positive integer `q`
the canonical section of `O_X(qD)` transports to a holomorphic modular form
of weight `qa`, up to a torsion line bundle.  A further power kills that
torsion.  On the simply connected type-IV domain, taking the corresponding
root gives a weight-`a` automorphic form with finite-order multiplier and
divisor `D`.  Equivalently, one may state the conclusion without roots:
some positive multiple `qD` is the divisor of an ordinary integral-weight
modular form.

This is a **positive-weight section**, not a principal divisor.  Indeed, `D`
is nonzero and effective.  On a connected projective variety a meromorphic
function with effective divisor has no poles, hence is a regular function,
hence is constant.  The same argument after taking a power excludes weight
zero with a finite-order character.  Since `D ell^9>0`, the forced weight
`a` is strictly positive; a section with divisor `D` is automatically
holomorphic.

## 3. Hodge-index decision theorem

Put

\[
\alpha=[D]-a\ell.
\]

By construction, `alpha ell^9=0`.  The Hodge--Riemann bilinear relations for
primitive real `(1,1)` classes on the tenfold `X` give

\[
\alpha^2\ell^8\le 0,
\]

with equality exactly when `alpha=0` in `H^2(X,R)`.  The left side is
precisely `Delta_D`.

For these compact orthogonal Shimura varieties of dimension at least three,
Yuan--Zhang--Zhang use the vanishing of `H^1(X,O_X)` (Kumaresan and
Vogan--Zuckerman) to show that the divisor cycle-class map is injective up
to torsion.  Consequently

\[
\alpha=0\text{ in }H^2(X,\mathbf R)
\iff
\alpha=0\text{ in }\operatorname{Pic}(X)_\mathbf Q.
\]

This proves the boxed decision theorem in the Outcome.  It also identifies

\[
-\Delta_D
\]

with a positive constant times the squared Hodge norm of the primitive
harmonic form `Lambda(g_0)`.

No general Picard-rank-one theorem is available here.  Compactness removes
the orthogonal boundary; it does not remove cuspidal `(1,1)` cohomology.
That primitive cohomology is exactly where the Kudla--Millson image can live.

## 4. A finite genus-two coefficient calculation for `Delta_D`

The Hodge-index test can be evaluated by special-cycle intersection theory.
Let

\[
\mathcal Z^{(1)}(\tau)=-\ell+
\sum_{m,\mu}Z(m,\mu)q^m\mathfrak e_\mu.
\]

Yuan--Zhang--Zhang's product formula (Theorem 1.1) says, in the Chow
ring,

\[
\mathcal Z^{(1)}(\tau_1)\mathcal Z^{(1)}(\tau_2)
=
\mathcal Z^{(2)}
\!\left(\begin{smallmatrix}\tau_1&0\\0&\tau_2\end{smallmatrix}\right).
\]

Pairing with `ell^8` and applying the geometric Siegel--Weil formula expresses
the resulting genus-two volume series as the coherent genus-two
Hilbert--Siegel Eisenstein series of parallel weight `6`, for the finite
Schwartz function `1_(Lhat^2)`.  Thus:

* `D ell^9 / ell^10` is the normalized `(m_0,0)` coefficient of the
  genus-one Eisenstein series; and
* `D^2 ell^8 / ell^10` is the normalized coefficient of
  `q_1^(m_0) q_2^(m_0)` in the block-diagonal restriction of the genus-two
  Eisenstein series.

More explicitly, normalize the genus-one Eisenstein series to have constant
coefficient `1`, and call its `(m_0,0)` coefficient `B_1(m_0,0)`.  Normalize
the genus-two series compatibly, and call the sum of its coefficients with
both diagonal entries `m_0` by `C_2(m_0,m_0)`.  With Bruinier's sign
convention and the ample class `ell`,

\[
B_1(m_0,0)=-a,
\qquad
C_2(m_0,m_0)=\frac{D^2\ell^8}{\ell^{10}},
\]

and therefore

\[
\boxed{
\frac{\Delta_D}{\ell^{10}}
=C_2(m_0,m_0)-B_1(m_0,0)^2.}
\]

Only nine genus-two indices can enter the latter coefficient.  If `x,y in L`
both have `Qhat(x)=Qhat(y)=m_0`, put `h=Hhat(x,y) in O`.  Their moment matrix is

\[
T_h=m_0
\begin{pmatrix}
1&h/2\\
h/2&1
\end{pmatrix}.
\]

Total positive semidefiniteness is equivalent to

\[
|\sigma_1(h)|\le2,
\qquad
|\sigma_2(h)|\le2.
\]

The algebraic integers in `Z[phi]` satisfying these two bounds are exactly

\[
h\in\{0,\ \pm1,\ \pm\phi,\ \pm(\phi-1),\ \pm2\}.
\]

The first seven give positive-definite rank-two matrices; `h=+/-2` are the
rank-one excess-intersection terms.  Yuan--Zhang--Zhang's definition inserts
`c_1(L^vee)=-ell` for precisely these excess terms, so they must not be
discarded.

An exact coefficient computation therefore requires the following data.

1. The normalized genus-one Whittaker coefficient at `m_0`.
2. The nine normalized genus-two Whittaker coefficients at the matrices
   `T_h` above.
3. The local Jordan decompositions of `L` and compatible additive-character
   conventions.  Since

   \[
   \det(\widehat H)=4(\phi+3),\qquad
   N_{F/\mathbf Q}(\det\widehat H)=176,
   \]

   the nonspherical lattice data occur at the prime over `2` and at the
   chosen prime `(phi+3)` over `11`.  The prime over `5` must also be tracked
   in the Fourier-index/additive-character normalization because
   `m_0 in partial^{-1}` and `partial=(sqrt(5))`.
4. The standard archimedean holomorphic Whittaker factors of weight `(6,6)`
   and one common Tamagawa/volume normalization.  Away from the displayed
   places the factors are the unramified local Siegel series.

In Langlands' normalization the two Siegel--Weil points are

\[
s_1=5\quad(r=1),
\qquad
s_2=\frac92\quad(r=2).
\]

Local Siegel-series/representation-density formulas make this a finite exact
calculation, but the dyadic factor and the normalization have not been
implemented here.  Until that calculation is completed, the sign of
`Delta_D` remains open.

## 5. Why the totally real hyperbolic slice does not decide the Picard class

Let `e` generate the added negative line at the identity embedding.  On the
type-IV domain, the antiholomorphic isometry

\[
z\longmapsto s_e(\bar z)
\]

has the embedded hyperbolic ten-space as a component of its fixed locus:
projectively it consists of lines represented by `e+i t x`, with `x` a
negative vector of `H` and `t` chosen to make the vector isotropic.

The fixed locus of an antiholomorphic isometry of a Kahler manifold is
Lagrangian.  Hence

\[
\iota^*\Omega=0,
\qquad
\iota^*c_1(\mathcal M_1)=0
\quad\text{in real cohomology}.
\]

Moreover, `D` meets this locus nontransversely: its intersection is the
real codimension-one norm-two mirror arrangement, whereas a transverse
intersection with a complex divisor would have real codimension two.  Thus
ordinary restriction of the Picard class to the hyperbolic locus cannot see
the primitive class `[D]-a ell`.  Restricting a hypothetical section records
the mirror zero set and a flat phase/sign character, but it does not decide
`Delta_D`, Whittaker vanishing, or a reflection index.

## 6. The pure Whittaker input is now proved obstructed

The required vector-valued identity-term estimate has been carried out.  After
normalizing the diagonal coefficient of the Hilbert Poincare series

\[
g_0=P_{m_0,0}=\xi_k(f_0)
\]

to have identity contribution `1`, the exact rational certificate in
`scripts/certify_og10_whittaker_nonvanishing.py` proves

\[
\left|b_{g_0}(m_0,0)-1\right|<0.931739<1.
\]

It includes the full Hilbert unit orbit, every lower-left row, the unitarity
bound for the generalized Weil Kloosterman sums, and the sharper unit-row
value

\[
|\rho_L(S)_{0,0}|=|L^\vee/L|^{-1/2}=176^{-1/2}.
\]

Consequently

\[
\boxed{g_0\ne0.}
\]

Thus `f_0` is not weakly holomorphic and the principal part consisting only
of the coefficient `1` at `(m_0,0)` is obstructed.  Equivalently, the
functional `g -> b_g(m_0,0)` on the cusp space is nonzero.  The derivation and
exact tail certificate are recorded in
`results/og10_whittaker_nonvanishing.md`.

The complete globally minimal off-slice compensator shell has since been
excluded as well: its two unit orbits, together with the target, give an
invertible `3x3` coefficient matrix; see
`results/og10_whittaker_first_compensator.md`.  This still does not exclude
another Whittaker input using larger-index or multiple compensator shells
whose special divisors vanish or cancel.  More
importantly, it does not prove that no scalar automorphic section has divisor
`D`: that stronger conclusion requires

\[
\Lambda(g_0)\ne0,
\]

equivalently `Delta_D<0`, or an injectivity theorem applicable to this exact
anisotropic nonunimodular lattice.

## 7. What Rallis nonvanishing does, and does not, settle here

There is a useful stable-range theorem, but it stops one quantifier short of
the assertion needed for `g_0`.

Rallis's Corollary 1 to Proposition 1.1 applies to an anisotropic quadratic
space over a totally real field when its dimension is greater than `6`.  It
characterizes the cuspidal `SL_2` representations occurring in the **full**
global theta lift by local occurrence in the Weil representation.  In the
proof (pages 155--156), at every exceptional finite place Rallis chooses
Schwartz functions `phi_v, phi'_v` suitably so that the local matrix-coefficient
integral is nonzero.  Since a quadratic space of dimension greater than `6`
has no anisotropic finite localization, Remark 1.4 reduces the
representation-level occurrence conditions to the real places.

For our `12`-dimensional space the archimedean occurrence condition is the
expected one: the parallel weight-`6` discrete-series type occurs in the Weil
representation, and the Kudla--Millson archimedean factor is nonzero because
the indefinite signature has `p=10>1`.  Consequently Rallis theory says that
each relevant automorphic representation has a nonzero orthogonal theta lift
for **some** choice of finite Schwartz data.

The Bruinier--Kudla--Millson lift in this problem does not allow that choice.
Its finite Schwartz data are the fixed lattice/discriminant-module type

\[
\mathcal S_L=\bigotimes_{v<\infty}\mathcal S_{L,v},
\qquad
\dim \mathcal S_L=|L^\vee/L|=176,
\]

spanned by the characteristic functions of the discriminant cosets.  For an
irreducible constituent `pi` of `g_0`, its norm has the schematic form

\[
\|\Lambda(g_\pi)\|^2
=C_\infty\,L^S(5,\pi,\mathrm{std})
\prod_{v\in S_f} Z_v(g_{\pi,v},\Phi_{L,v}).
\]

Here the partial standard `L`-value is in its absolute-convergence region and
is nonzero, and `C_infty` is nonzero.  The remaining question is whether the
local doubling functional is nonzero on the **specified** lattice type at
every bad finite place.  Local occurrence

\[
\operatorname{Hom}_{SL_2(F_v)}(\omega_{V_v},\pi_v^\vee)\ne0
\]

only says that this functional is nonzero somewhere in the full Schwartz
space.  It does not say that its restriction to
\(\mathcal S_{L,v}\) is nonzero.
Allowing all `176` discriminant components does not remove this logical gap:
one must show that the projection of the whole local lattice type, rather
than merely one of its basis vectors, survives the local theta quotient.

Bruinier--Funke make exactly this limitation explicit.  For a nonunimodular
lattice they obtain the nonzero partial standard `L`-value times bad local
factors and state that the lift vanishes precisely when a bad factor
vanishes; they remove the finite bad factors only in the even-unimodular
specialization and note that finite-level vector-valued Hecke theory would be
needed for the general case.  Thus neither their theorem nor Rallis's older
theorem supplies fixed-vector injectivity for the present anisotropic lattice.

### Corollary 7.1 (variable-level primitive divisor theorem)

Let `Pi(g_0)` be the set of irreducible cuspidal automorphic representations
of `SL_2(A_F)` whose projection in the automorphic representation generated by
`g_0` is nonzero.  For every `pi in Pi(g_0)` there exist

1. a nonzero cusp vector `g_pi in pi` with archimedean weight `(6,6)`;
2. a factorizable finite Schwartz function
   `phi_f in S(V(A_(F,f)))`; and
3. a neat compact open
   `K' subset K cap Stab(phi_f)`

such that the Kudla--Millson lift with the standard archimedean
Kudla--Millson form satisfies

\[
\boxed{\Lambda_{K',\phi_f}(g_\pi)\ne0
\quad\text{in }H^{1,1}(X_{K'},\mathbf C).}
\]

The finite Schwartz function may be chosen `Q`-valued.  Consequently the
cohomology-valued special-divisor generating series on some finite algebraic
cover `X_(K') -> X_K` has nonzero cuspidal, hence primitive, projection.  In
particular, some rational linear combination of algebraic special divisors on
`X_(K')` has nonzero primitive `(1,1)` class.

Indeed, the Fourier coefficients of the cohomology-valued
Kudla--Millson theta series are the Poincare-dual classes of these special
divisors.  If every primitive coefficient vanished, its Petersson pairing
with every cusp vector, including `g_pi`, would vanish.

Here are the exact hypotheses used.  The field quadratic space `V` has
dimension `12>6` and is anisotropic over `F` because it is positive definite
at the second real place.  At the positive-definite place the weight-`6`
discrete series is the lowest type in the Weil representation; at the
signature `(10,2)` place Rallis's proof includes all discrete series because
both signs have dimension greater than one.  Thus every `pi in Pi(g_0)` meets
the archimedean occurrence conditions in Rallis's Corollary 1 to Proposition
1.1.  His proof on pages 155--156 chooses the finite Schwartz functions so
that every exceptional local integral is nonzero.  Bruinier--Funke's
archimedean calculation shows that the fixed Kudla--Millson factor is nonzero
for `p=10>1`.  The resulting theta form is closed and harmonic by Bruinier's
definition (6.10).  Finally, a finite Schwartz function is a finite complex
linear combination of characteristic functions of compact opens.  Since the
lift is linear, one characteristic summand still has nonzero lift; it is
`Q`-valued and fixed by an open compact subgroup, which can be intersected
with `K` and shrunk to be neat.

Thus a possible vanishing at the original `176`-component level would be a
genuine bad-local-type phenomenon, not representation-level disappearance of
the theta lift.

### Finite covers, cycle combinations, and the hyperbolic slice

The variable-level theorem has an entirely algebraic interpretation.  If
`phi_f` is fixed by `K'`, then every Fourier coefficient `Z(m,phi_f)` is a
finite `Q`-linear combination of orthogonal Shimura divisors on `X_(K')`.
After making both groups neat, the map

\[
p:X_{K'}\longrightarrow X_K
\]

is a finite algebraic cover.  Pullback refines a special divisor into its
level-`K'` components, while pushforward is, up to the standard stabilizer
multiplicities, the finite trace of the Schwartz function over the deck
cosets.

There is also an exact test for retaining the OG10 mirrors.  Write
`V=H direct_sum F e` and let `iota` be the distinguished hyperbolic slice in
the identity adelic component.  For the index `m_0`, restriction has the form

\[
\iota^{-1}Z(m_0,\phi_f)
=\sum_{\substack{[r]\subset H(F)\\H(r,r)=2}}
w_{\phi_f}([r])\,\mathcal M_r,
\tag{7.1}
\]

where `[r]` is an unoriented rational root line, `M_r` is its hyperbolic
mirror, and `w_(phi_f)([r])` is the corresponding orbit sum of the values at
`r` and `-r` (with the same stabilizer convention used in the definition of
the special divisor).  A vector `lambda=r+a e` with `a != 0` contributes no
mirror because `(e,lambda) != 0`; its type-IV divisor misses the slice.

Let `chi_L` be the characteristic function defining the original divisor
`D`.  Equation (7.1) gives a necessary and sufficient condition for the same
complete mirror arrangement, up to one common multiplicity `c`:

\[
w_{\phi_f}([r])=
\begin{cases}
c\,w_{\chi_L}([r]),&r\in H(\mathcal O),\\
0,&r\notin H(\mathcal O),
\end{cases}
\qquad H(r,r)=2.
\tag{7.2}
\]

A convenient stronger sufficient condition is

\[
\phi_f-c\chi_L=0
\quad\text{on }\{r\in H(F):H(r,r)=2\},
\tag{7.3}
\]

after symmetrizing under `r -> -r`.  Condition (7.3) permits arbitrary extra
Schwartz support off `H`; geometrically this adds algebraic divisors that are
disjoint from the hyperbolic slice.

Rallis's theorem does **not** assert that its finite Schwartz function lies
in the affine subspace (7.2), nor that its nonzero primitive coefficient has
index `m_0`.  Averaging over deck transformations can restore the complete
mirror weights, but it can also kill a primitive class lying in a nontrivial
deck representation.  Equivalently, a nonzero class upstairs may have zero
pushforward to `X_K`.

The theorem therefore proves nontrivial primitive special-divisor cohomology
somewhere in the tower, but it does not yet sharpen the fixed OG10 divisor
decision.  It does suggest a finite local refinement of the divisor-lattice
factory: at the bad places search inside the affine constraint (7.2) for
Schwartz data on which every local doubling functional is nonzero.  Success
would preserve the full norm-two mirror restriction while allowing extra
type-IV components disjoint from the real slice.  Even then, finite volume of
the resulting real Coxeter chamber would remain a separate Vinberg/nef-cone
condition.

The residual finite problem is small and concrete.  Decompose `g_0` into
irreducible Hecke constituents and, for each constituent, compute the local
doubling operator on \(\mathcal S_{L,v}\) at the primes over `2` and `11`; the
place over `5` must also be checked under the chosen different/additive-character
normalization.  If one constituent has nonzero local factors at all these
places, positivity of the Rallis norm gives

\[
\Lambda(g_0)\ne0,
\qquad
\Delta_D<0.
\]

This is substantially narrower than an unrestricted search for an
injectivity theorem, but it is still an unperformed local calculation.

## Primary sources

* J. H. Bruinier,
  [*Regularized theta lifts for orthogonal groups over totally real fields*](https://arxiv.org/abs/0908.3076),
  especially Proposition 4.3, Proposition 4.5, Corollary 5.16,
  Theorem 6.4, Remark 6.5, Theorem 6.8, and Theorem 7.1.
* X. Yuan, S.-W. Zhang, and W. Zhang,
  [*The Gross--Kohnen--Zagier theorem over totally real fields*](https://doi.org/10.1112/S0010437X08003734),
  especially Theorems 1.1--1.3, Proposition 2.6, and Section 4.1.
* S. Kudla,
  [*Special cycles and derivatives of Eisenstein series*](https://arxiv.org/abs/math/0308295),
  Theorem 4.1 and the compact volume formula (4.4).
* J. H. Bruinier and J. Funke,
  [*On the injectivity of the Kudla--Millson lift and surjectivity of the Borcherds lift*](https://arxiv.org/abs/math/0606178),
  especially Corollary 1.7, Theorem 1.8, and Corollary 4.12.
* I. Metzler and R. Zuffetti,
  [*Injectivity of the genus 1 Kudla--Millson lift on locally symmetric spaces*](https://arxiv.org/abs/2312.00572),
  Theorem 1.1 (whose splitting hypotheses do not hold here).
* S. Rallis,
  [*Injectivity properties of liftings associated to Weil representations*](http://www.numdam.org/item/CM_1984__52_2_139_0/),
  Compositio Math. 52 (1984), especially Proposition 1.1, its first
  corollary, Remark 1.4, and the proof on pages 155--156.
* M. Zhang and Y. Zhang,
  [*On the non-vanishing of Hilbert Poincare series*](https://arxiv.org/abs/2312.04169),
  especially the Fourier-coefficient formula in Section 2 and the Bessel
  bounds used in Theorem 3.2; the vector-valued adaptation is audited in
  `results/og10_whittaker_nonvanishing.md`.

## Scope

The Hodge-index/Eisenstein calculation decides the existence of a scalar
automorphic section with the single divisor `Z(m_0,0)`.  It does not by itself
decide whether that section is a Bruinier Whittaker lift, and neither outcome
implies anything about the index or covolume of the real-hyperbolic reflection
group.
