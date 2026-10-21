# Compact determinant-line Shimura surfaces: canonical sections and exact BMY cancellation

**Status (2026-08-09): proved compact-branch no-go and a closed surface-normalization
lemma.**  The automatic determinant-line construction does admit a complete-intersection
surface in a multiple of the automorphic line bundle bounded solely by the complex
dimension.  Thus the freedom to replace the surface polarization by arbitrarily high
multiples is not an essential obstruction on the compact type-IV carrier.

Nevertheless, ordinary BMY on that canonically bounded surface is an explicit strict
positive identity which contains neither the Coxeter facets nor their labels.  Restricting
an effective automorphic carrier divisor to the surface gives exactly the existing ambient
sheet-to-degree inequality; the complete-intersection multiple cancels.  Hence the surface
section itself supplies no new dimension or facet bound.  The remaining feedback gate is
precisely a class-uniform, noncircular upper bound for the P-derived automorphic top degree.

This note is scoped to the compact anisotropic determinant-line branch.  It makes no
finite-volume noncompact inference.

## 1. Setup

Let $P\subset\mathbb H^n$, $n\ge2$, be a marked compact arithmetic Coxeter
polytope to which the determinant-line compiler associates a type-IV arithmetic
quotient.  Pass to the P-computable exact-support level and then to a neat normal
finite-index subgroup.  Since the quadratic datum is anisotropic, the resulting
quotient

\[
 X=\Gamma'\backslash\mathcal D
\]

is a smooth compact projective $n$-fold.  Let $\lambda$ be its weight-one
automorphic line bundle.  It is ample, and numerically

\[
 K_X=n\lambda. \tag{1.1}
\]

Finite automorphy characters do not affect any formula below.  More concretely,
the very ample divisor used below is $2K_X+m\lambda$, so no linear equivalence
$K_X\simeq n\lambda$ is needed.

## 2. A dimension-uniform complete-intersection surface

### Theorem 2.1 -- canonical bounded polarization

Put

\[
 m_n=2+\binom{3n+1}{n},\qquad
 k_n=2n+m_n=2n+2+\binom{3n+1}{n}. \tag{2.1}
\]

Then $2K_X+m_n\lambda$ is very ample and has numerical class
$k_n\lambda$.  Consequently $n-2$ general members cut out a smooth
projective surface

\[
 \Sigma=H_1\cap\cdots\cap H_{n-2},\qquad [H_i]=k_n\lambda, \tag{2.2}
\]

transverse to every member of any fixed finite special-divisor arrangement.

#### Proof

Demailly's effective very-ampleness theorem says that $2K_X+mL$ is very
ample on every smooth projective $n$-fold when $L$ is ample and
$m\ge2+\binom{3n+1}{n}$.  Apply it to $L=\lambda$.  Equation (1.1)
gives the numerical class in (2.1), and Bertini plus stratified transversality
gives (2.2).  $\square$

This is enormous rather than sharp--for $n=8$, $k_8=1,081,593$--but it
closes the logical normalization issue: a P-independent multiple depending only
on $n$ exists.  One cannot blame the failure of BMY below merely on an
unbounded arbitrary choice $H\mapsto aH$.

## 3. Exact Chern calculation

### Theorem 3.1 -- proportionality formula

For any integer $k\ge1$, let $\Sigma_k$ be a smooth complete intersection
of $n-2$ members numerically equivalent to $k\lambda$.  Put

\[
 V_X=\int_X\lambda^n>0,\qquad r=n-2.
\]

Then

\[
 K_{\Sigma_k}^2=(n+rk)^2k^rV_X, \tag{3.1}
\]

\[
 c_2(\Sigma_k)=
 \left(
 \frac{n^2-n+2}{2}+rnk+\frac{r(r+1)}2k^2
 \right)k^rV_X, \tag{3.2}
\]

and therefore

\[
\boxed{
3c_2(\Sigma_k)-K_{\Sigma_k}^2
=\left(
 \frac{n^2-3n+6}{2}
 +n(n-2)k
 +\frac{(n-2)(n+1)}2k^2
 \right)k^{n-2}V_X.}
\tag{3.3}
\]

The right side is strictly positive for every $n\ge2$ and $k\ge1$.

#### Proof

The compact dual of the type-IV domain is the smooth quadric
$Q^n\subset\mathbf P^{n+1}$.  If $h$ is its hyperplane class, the normal
sequence gives

\[
 c(TQ^n)=\frac{(1+h)^{n+2}}{1+2h},
\]

so $c_1(TQ^n)=nh$ and
$c_2(TQ^n)=\frac{n^2-n+2}{2}h^2$.  Hirzebruch proportionality, with
$h\mapsto-\lambda$, yields

\[
 c_1(TX)=-n\lambda,\qquad
 c_2(TX)=\frac{n^2-n+2}{2}\lambda^2. \tag{3.4}
\]

Now divide $c(TX)|_{\Sigma_k}$ by
$(1+k\lambda)^r$, or equivalently use adjunction.  This gives (3.1) and
(3.2); subtraction gives (3.3).  Its constant coefficient is always positive
(the quadratic $n^2-3n+6$ has discriminant $-15$), and its linear and
quadratic coefficients are nonnegative for $n\ge2$.  $\square$

### Corollary 3.2 -- ordinary BMY is structurally silent

For the universal surface $\Sigma_{k_n}$, ordinary BMY is already the strict
identity (3.3).  It cannot bound $n$, the number of facets, any Coxeter label,
or the number of special-divisor components.  All input dependence occurs through
the positive scalar $V_X$, so even an exact value of the arithmetic volume does
not change the sign.

The asymptotic Chern slope is

\[
 \lim_{k\to\infty}\frac{K_{\Sigma_k}^2}{c_2(\Sigma_k)}
 =\frac{2(n-2)}{n-1}<2, \tag{3.5}
\]

uniformly far from the BMY boundary $3$.

If a **single smooth** automorphic divisor restricts to
$C\equiv w\lambda|_{\Sigma_k}$, then

\[
3\bar c_2-(K_{\Sigma_k}+C)^2
=(3c_2-K_{\Sigma_k}^2)+K_{\Sigma_k}C+2C^2>0. \tag{3.6}
\]

Thus adding a smooth positive boundary makes the inequality still less
restrictive.  For a reduced SNC union $C=\sum C_i$, the correction is

\[
 K_{\Sigma_k}C+2\sum_iC_i^2+\sum_{i<j}C_iC_j. \tag{3.7}
\]

For the actual separated mirror components, even that apparent possibility
closes negatively.

### Proposition 3.3 -- the SNC mirror-boundary branch is tautological

Let $Z_i$ be distinct mirror components on the neat exact-support cover and
let $C_i=Z_i|_{\Sigma_k}$.  **Assume the induced reduced curve arrangement is
SNC**, equivalently for this calculation that no codimension-two stratum met by
the surface lies on three or more selected mirror components.  Put

\[
 d_i=\lambda C_i>0.
\]

Then

\[
 C_i^2=-d_i. \tag{3.8}
\]

Indeed, Bruinier defines the positive weight-one bundle $\mathcal M_1$, our
$\lambda$, as the quotient of the **tautological** bundle $\mathcal O(-1)$
on the type-IV domain.  The equation $z\mapsto(z,r)$ of a mirror is
fiberwise linear on that tautological line, hence is a section of its dual.
Thus the mirror normal is $\lambda^{-1}$.  On the neat cover the stabilizer
cannot send $r$ to $-r$, since that would give the forbidden torsion
eigenvalue $-1$ on the normal line.  The bundle therefore descends with
factor one, not the coarse reflection-branch factor two.  As a consistency
check, adjunction gives $K_{Z_i}=(n-1)\lambda|_{Z_i}$, the canonical formula
for the dimension-$n-1$ type-IV mirror quotient.  Restriction to $C_i$ gives
(3.8).

Writing $a=n+(n-2)k$, equations (3.3), (3.7), and (3.8) give the exact
identity

\[
\boxed{
 3\bar c_2-(K_{\Sigma_k}+C)^2
 =\Delta_0+(a-2)\sum_i d_i+\sum_{i<j}C_iC_j,}
 \tag{3.9}
\]

where $\Delta_0=3c_2(\Sigma_k)-K_{\Sigma_k}^2>0$.  Since
$a-2\ge0$ for every $n\ge2,k\ge1$, and all pair intersections are
nonnegative, reduced log BMY is strictly satisfied for every such **SNC**
separated mirror arrangement, regardless of the number of facets or labels.

Moreover $K_{\Sigma_k}+C$ is ample: it has positive intersection with a
boundary component because

\[
 (K_{\Sigma_k}+C)C_i=(a-1)d_i+\sum_{j\ne i}C_iC_j>0,
\]

and with every other curve because $K_{\Sigma_k}$ is ample and $C$ is
effective.  Thus (3.9) is not evading a log-nef hypothesis.

The SNC hypothesis is substantive.  A general surface avoids ambient strata
of codimension at least three, but it cannot avoid a codimension-two stratum
contained in three or more mirrors.  Such strata occur when several rational
quadratic hyperplanes have the same codimension-two intersection (for example
inside a rank-two finite root subsystem), and deck translates can create them
away from the marked real chamber.  Exact real-slice support does not exclude
this.  In the general arrangement one must blow up the multiple points and
include the exceptional-divisor/local-orbifold corrections.  Those corrections
are **not** proved nonnegative here and remain a legitimate label-sensitive
feedback route.

The first local calculation is nevertheless exact.  At an ordinary point with
$b$ transverse mirror branches, one blowup with reduced total boundary gives

\[
 K_{\widetilde\Sigma}+\widetilde C
 =\pi^*(K_\Sigma+C)+(2-b)E. \tag{3.10}
\]

Thus a reduced boundary with $b\ge3$ is not log canonical, so ordinary
reduced log BMY cannot simply be applied across that point.  The finite next
calculation is to enumerate every $b$ from the P-derived special-divisor
arrangement and test the actual stack/orbifold coefficients (for reflection
branches, typically one-half) for log canonicity and log-nefness after
resolution.  More generally, for transverse branches of coefficients
$\alpha_1,\ldots,\alpha_b$, the exceptional discrepancy is
$1-\sum_i\alpha_i$, so the local pair is log canonical exactly when
$\sum_i\alpha_i\le2$.  This is exactly where a nontrivial label term could survive;
equation (3.9) proves that it cannot come from the SNC locus.

There is an additional label warning.  On the original arithmetic DM stack,
the reflection inertia and all dihedral labels are already part of the orbifold
tangent bundle.  Hirzebruch proportionality computes its orbifold Chern classes
by the same formula (3.4).  Hence the **natural quotient-stack orbifold BMY
number** again equals (3.3): the branch components and local inertia orders are
absorbed into the stack Chern number, rather than appearing as a separately
signed combinatorial penalty.  A label-sensitive inequality would need a new
log modification or boundary weighting, not merely the canonical quotient stack.

## 4. Surface degree is exactly the ambient degree inequality

Let an effective automorphic divisor $D\equiv w\lambda$ on the exact-support
neat cover contain a distinct component $Y_i$ over every marked facet, and let
$q_{\min}>0$ be the least such multiplicity.  Exact support makes the $Y_i$
pairwise distinct.  On $\Sigma_k$, set $C_i=Y_i|_{\Sigma_k}$.  Ampleness and
integrality give

\[
 C_i\lambda|_{\Sigma_k}
 =k^{n-2}Y_i\lambda^{n-1}\ge k^{n-2}.
\]

Consequently

\[
 \#\operatorname{Facets}(P)\,q_{\min}k^{n-2}
 \le D\lambda(k\lambda)^{n-2}
 =wk^{n-2}\lambda^n.
\]

Cancelling $k^{n-2}$ gives exactly

\[
\boxed{
 \#\operatorname{Facets}(P)
 \le \frac{w\lambda^n}{q_{\min}}.}
\tag{4.1}
\]

This is the ambient sheet-to-degree theorem, neither stronger nor weaker.
The canonical surface of Theorem 2.1 makes the association uniform but creates
no new numerical leverage.

## 5. The exact P-only condition that would produce feedback

For a surface-free compact arithmetic class $\mathcal C$, define the
**normalized exact-support carrier degree** of a selected automatic package by

\[
 U_{\mathrm{det}}(P)=\frac{w_P\lambda_P^n}{q_{\min}(P)}. \tag{5.1}
\]

All terms are arithmetic outputs of the determinant-line compiler once the
integral Hilbert q-structure and positive correction have been made effective.
The precise sufficient feedback condition is

\[
 \boxed{\sup_{P\in\mathcal C}U_{\mathrm{det}}(P)\le U_{\mathcal C}<\infty,}
 \tag{5.2}
\]

with the bound proved independently of the facet count.  Under (5.2), (4.1)
immediately returns the genuine theorem

\[
 \#\operatorname{Facets}(P)\le\lfloor U_{\mathcal C}\rfloor
 \quad\text{for every }P\in\mathcal C. \tag{5.3}
\]

The existing finite Cramer--Hadamard and Eisenstein/degree formulas reduce
(5.2) to the following **finite, P-first proof obligations**:

1. bound the exact sheet-separating level and the index of a normal neat cover
   from the raw class invariants, without using the facet count;
2. prove an effective integral Hilbert q-lattice (or a uniformly bounded
   denominator) for the obstruction space;
3. bound the higher-degree diagonal-detection cutoff and the integral
   obstruction-column height on the finite safe support;
4. bound the Prasad/local-density volume and cover index; and
5. keep the target mass in the weight estimate sublinear in, or independent of,
   the number of facets.

Each item terminates for a fixed $P$.  Their **class-uniform** bounds are not
proved.  This is the exact obstruction, rather than a missing choice of surface
or an unperformed BMY calculation.

For a label-sensitive refinement, one needs one additional theorem: on a common
normal exact-support cover retaining the finite deck-group inertia, enumerate
every codimension-two multiple-mirror stratum, resolve it, and prove that the
resulting exceptional/local-orbifold term has a uniformly negative contribution
per selected Coxeter stratum.  The natural quotient-stack boundary cannot do
this, because its orbifold Chern number is already (3.3); a deliberately chosen
log boundary or coarse resolution may still do so.

## 6. Position on the strict inverse-carrier spine

The strongest honest ledger is:

| flagship clause | compact determinant-line status |
|---|---|
| arbitrary $P$ in a surface-free compact arithmetic class | passes for the marked determinant-line-eligible class |
| finite exhaustive compiler | passes for determinant-line gluings, finite hidden sheets, exact level, and fixed-support obstruction; it is not an exhaustive list of every possible algebraic carrier species |
| intrinsic iff returning the same $P$ | passes for the marked real-slice divisor arrangement; it does not identify the surface nef/Mori cone with $P$ |
| essential surjectivity | passes only for the automatic arrangement shadow, not for all divisor-lattice/nef carriers |
| AG-native theorem back on all $P$ plus a positive survivor | fails: Theorems 2.1 and 3.1 prove that ordinary/natural-orbifold BMY supplies no bound; (5.2) is the remaining quantitative gate |

Thus this result is not flagship completion and should not be advertised as a
new compact polytope theorem.  Its value is that it closes the canonical-surface
question negatively and leaves a finite quantitative obligation.

## 7. Bugaenko calibration and the coherent-K3 alternative

For Bugaenko's compact $H^8/11$ chamber, (3.3) specializes to

\[
 (27k^2+48k+23)k^6V_X>0. \tag{7.1}
\]

At Demailly's universal $k_8=1,081,593$, the coefficient in parentheses is
$31,585,824,193,010$.  Hence Bugaenko remains a positive exact arrangement
carrier, but it does not make BMY sharp and no new abstract polytope follows.

The other proposed compact shortcut--declaring every NS eigencarrier to be a
full-$H^2$ Hodge-coherent K3 carrier--is false.  At $H^9$, the trace Picard
lattice has rank twenty and its K3 complement has rational rank two, hence
rank one over a real quadratic field.  One real eigenspace cannot contain the
two-dimensional period plane.  Existing NS-only rank-twenty packages therefore
show that primitive K3 embedding does not force Hodge coherence.  A valid
P-only class must supply a genuine complement-transfer criterion; simply adding
"coherent carrier exists" to the definition would select by carrier success and
would not prove the flagship essential-surjectivity clause.

## 8. Replay and sources

Run

```sh
python3 -m py_compile scripts/certify_compact_shimura_ci_bmy_no_go.py
PYTHONPATH=src:. python3 scripts/certify_compact_shimura_ci_bmy_no_go.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_compact_shimura_ci_bmy_no_go.py
```

The executable recomputes (3.1)--(3.3), proves coefficient positivity, and
checks the $H^8$ specialization.

Primary inputs:

1. J.-P. Demailly, *Effective bounds for very ample line bundles*, Invent.
   Math. **124** (1996), 243--261, Theorem 0.2(a).
2. F. Hirzebruch, *Automorphe Formen und der Satz von Riemann--Roch*, in
   *Symposium Internacional de Topologia Algebraica* (1958), for proportionality;
   equivalently use its standard compact locally symmetric formulation.
3. J. H. Bruinier, *Regularized theta lifts for orthogonal groups over totally
   real fields*, J. Reine Angew. Math. **672** (2012), for $\lambda$, its
   curvature, and the automorphic weight--degree identity.
4. The upstream exact-support and degree inputs are
   [`polytope_to_automorphic_divisor_transfer_theorem.md`](polytope_to_automorphic_divisor_transfer_theorem.md),
   [`explicit_sheet_separating_level_bound.md`](explicit_sheet_separating_level_bound.md),
   [`facet_sheet_automorphic_degree_injection.md`](facet_sheet_automorphic_degree_injection.md),
   and [`borcherds_weight_degree_quantitative_bound.md`](borcherds_weight_degree_quantitative_bound.md).
