# Codimension-two mirror strata: exact local BMY no-go

**Status (2026-08-09): proved local no-go with exact replay.**

Resolving ordinary multiple-mirror points does not repair the compact
determinant-line feedback problem.  For a reduced boundary, blowing up an
ordinary `b`-fold point increases the BMY defect by `(b-2)^2`.  For arbitrary
log-canonical fractional weights, Langer's local orbifold bound implies that
both the coefficient-one and crepant exceptional-divisor corrections are
nonnegative.  In the reflection coefficient `1/2` specialization, the exact
dihedral point terms for `I_2(m)` with `m=2,3,4` are respectively

```text
1/4, 3/16, 0.
```

For `m>=5` the half-weight pair is not log canonical, so orbifold BMY does not
apply.  The natural dihedral quotient stack remains log smooth, but its local
mass `1/(2m)` is already absorbed in the orbifold tangent bundle and
Hirzebruch proportionality.  Thus no codimension-two stratum supplies a
negative-per-label or negative-per-facet term.

This is a compact anisotropic no-go.  It gives no result for finite-volume
noncompact polytopes and no new abstract Coxeter-polytope bound.

## 1. Ordinary multiple points and weighted blowup

Let `Sigma` be a smooth projective surface and let

\[
             \Delta=\sum_{i=1}^b\alpha_i C_i,
             \qquad 0\le\alpha_i\le1,
\]

where the smooth curves `C_i` have distinct tangent directions at an ordinary
multiple point `p`.  Put

\[
                     A=\sum_i\alpha_i.
\]

Let `pi:Y->Sigma` be the blowup, let `C_i'` be the strict transforms, and give
the exceptional curve `E` coefficient `beta`.  Then

\[
 K_Y+\sum_i\alpha_iC_i'+\beta E
 =\pi^*(K_\Sigma+\Delta)+(1-A+\beta)E.                 \tag{1.1}
\]

Consequently

\[
 \delta\bar c_1^2=-(1-A+\beta)^2.                     \tag{1.2}
\]

Write `e_p=e_orb(p;Sigma,Delta)` for Langer's local orbifold Euler
number.  The exceptional open stratum is `P1` minus `b` points and each new
crossing has product weight.  Direct stratification therefore gives

\[
 \delta e_{\rm orb}=(1-\beta)(2-A)-e_p.               \tag{1.3}
\]

For the BMY defect `B=3e_orb-(K+Delta)^2`, equations (1.2)--(1.3) give the
exact local resolution correction

\[
 \boxed{\delta B
 =3\big((1-\beta)(2-A)-e_p\big)+(1-A+\beta)^2.}        \tag{1.4}
\]

This formula is independent of the global automorphic volume, complete-
intersection degree, or Coxeter combinatorics.

## 2. General log-canonical sign theorem

Langer proves for a log-canonical ordinary point that `A<=2` and

\[
                     e_p\le \left(1-\frac A2\right)^2
                     =\frac{(2-A)^2}{4}.              \tag{2.1}
\]

Two natural resolution conventions now have a forced sign.

### Theorem 2.1 -- coefficient-one exceptional divisor

Put `beta=1`.  Then

\[
 \delta B=(2-A)^2-3e_p
 \ge\frac{(2-A)^2}{4}\ge0.                            \tag{2.2}
\]

### Theorem 2.2 -- crepant exceptional coefficient

When `1<=A<=2`, put `beta=A-1`.  Equation (1.1) is crepant and

\[
 \delta B=3\big((2-A)^2-e_p\big)
 \ge\frac94(2-A)^2\ge0.                               \tag{2.3}
\]

Thus every log-canonical weighted ordinary point makes the resolved BMY
inequality unchanged or weaker.  There is no choice between these standard
exceptional coefficients that yields a uniformly negative correction.

For a **reduced** boundary with `b>=3`, the unresolved pair is not log
canonical.  Nevertheless its reduced log resolution has an even simpler
exact identity.  The complement is unchanged, hence `delta c2_log=0`, while

\[
 K_Y+C'_{\rm red}+E=\pi^*(K_\Sigma+C_{\rm red})+(2-b)E.
\]

Therefore

\[
 \boxed{\delta\big(3\bar c_2-\bar c_1^2\big)=(b-2)^2>0.} \tag{2.4}
\]

The multiple point contributes with the wrong sign for a non-existence
argument.

## 3. Exact half-weight dihedral values

On a neat exact-support cover, a rank-two reflection subsystem gives an
ordinary multiple-mirror point on a general complete-intersection surface.
If the full finite subsystem is present, type `I_2(m)` has `b=m` mirror
branches and local reflection group order `2m`.

Give every mirror coefficient `alpha=1/2`.  The pair is log canonical exactly
when `m/2<=2`, hence only for `m=2,3,4`.  Langer's exact three-line theorem and
ordinary-line bound give

| `m` | `A=m/2` | `e_p` | BMY point term `Phi_m` |
|---:|---:|---:|---:|
| 2 | 1 | `1/4` | `1/4` |
| 3 | `3/2` | `1/16` | `3/16` |
| 4 | 2 | 0 | 0 |

Here the point term after separating the global curve contributions is

\[
 \Phi_m
 =3\left(\frac m2+e_p-1\right)-\frac{m(m-1)}4.       \tag{3.1}
\]

All three values are nonnegative.  The corresponding coefficient-one
exceptional corrections from (2.2) are `1/4,1/16,0`; the crepant corrections
from (2.3) are `9/4,9/16,0`.

For completeness, let `Sigma_k` be the canonical complete-intersection
surface in a compact type-IV `n`-fold and put

\[
                 a=n+(n-2)k,\qquad d_i=\lambda C_i,
                 \qquad C_i^2=-d_i.
\]

For an ordinary half-weight arrangement with only log-canonical multiple
points, the exact global defect is

\[
 \boxed{
 3e_{\rm orb}(\Sigma_k,\tfrac12C)
 -(K_{\Sigma_k}+\tfrac12C)^2
 =\Delta_0+\frac12\left(a-\frac52\right)\sum_i d_i
 +\sum_p\Phi_{b_p}.}                                  \tag{3.2}
\]

For `n>=3`, `a>=4`, so every displayed term is nonnegative and `Delta_0>0`.
Orbifold BMY is then a strict positive identity for every number of mirror
components.  If an `m>=5` point occurs, the half-weight pair is non-lc and
Langer's BMY theorem is unavailable; failure of a hypothesis is not a
polytope obstruction.

## 4. Why natural dihedral inertia also cannot help

The half-weight arrangement on a neat cover is a deliberately chosen log
pair.  The natural arithmetic Deligne--Mumford stack is different.  At a
rank-two stratum its local group is `I_2(m)` of order `2m`.  Finite-cover
multiplicativity gives local orbifold mass

\[
                              \frac1{2m}.              \tag{4.1}
\]

But that inertia is part of the orbifold tangent bundle itself.  Hirzebruch
proportionality therefore computes the complete stack Chern numbers before
any separate mirror penalty is added.  Extracting `1/(2m)` a second time
would double count the same local stabilizer.  The natural-stack BMY defect
remains the positive proportionality identity of the upstream compact
surface theorem.

Thus the three available local constructions exhaust the obvious route:

```text
reduced resolution       -> positive +(b-2)^2 correction;
lc weighted resolution   -> nonnegative correction by Langer;
natural dihedral stack   -> label absorbed by proportionality.
```

None yields a negative amount per Coxeter ridge, stratum, or facet.

## 5. Feedback verdict and remaining completion gate

This closes the codimension-two loophole left in
[`compact_shimura_ci_bmy_feedback_no_go.md`](compact_shimura_ci_bmy_feedback_no_go.md):
multiple mirrors cannot turn determinant-line BMY into primal leverage.

The exact sheet inequality remains

\[
 \#\operatorname{Facets}(P)
 \le U_{\rm det}(P)
 =\frac{w_P\lambda_P^n}{q_{\min}(P)}.                 \tag{5.1}
\]

A genuine compact theorem still requires a class-uniform upper bound for
`U_det` proved without using the facet count, or a different AG invariant not
absorbed by locally symmetric proportionality.  The local multiplicity and
dihedral data provide neither.  Accordingly this report is a **sharp local
no-go**, not strict inverse-carrier completion: it produces no new theorem on
every compact arithmetic polytope and no positive endpoint beyond the already
known Bugaenko calibration.

## 6. Replay and source

The exact rational replay is
[`certify_compact_mirror_multiple_strata_bmy_no_go.py`](../scripts/certify_compact_mirror_multiple_strata_bmy_no_go.py),
with frozen output
[`compact_mirror_multiple_strata_bmy_no_go.json`](../data/enumerations/compact_mirror_multiple_strata_bmy_no_go.json)
and tests
[`test_compact_mirror_multiple_strata_bmy_no_go.py`](../tests/test_compact_mirror_multiple_strata_bmy_no_go.py).

The local orbifold input is Adrian Langer,
[*Logarithmic orbifold Euler numbers of surfaces with applications*](https://arxiv.org/abs/math/0012180),
especially Theorem 0.1, property (0.2.3), and Theorems 8.3 and 8.7.  The
global type-IV Chern and mirror-normal identities are proved in the upstream
compact complete-intersection report.

