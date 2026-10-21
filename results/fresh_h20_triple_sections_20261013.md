# Three-wall H20 sections: a continuous obstruction and 64 exact angle repairs

No new Coxeter polytope was found. This round proves a stronger, precisely restricted section obstruction:

> In the frozen Borcherds H21 chamber, every interior H20 section through the fixed center `h` whose normal lies in the real span of at most three original wall normals is non-Coxeter.

The result covers **all 1,521,520 unordered wall triples and all real coefficients satisfying `(e,h)=0`**, not a bounded coefficient grid. Two-wall boundary cases use the [previous all-parameter theorem](fresh_h20_global_20261013.md). Sections with three-wall normals away from `h`, or normals requiring four or more walls, are not excluded by this theorem.

A separate noncentral construction repaired a chosen actual dihedral angle exactly. It produced **64 algebraic parameter choices in 17 symmetric elliptic-triple pencils**, including 51 quartic parameters. All 64 finite-volume sections have another exact forbidden actual ridge. This second result covers the declared pencils and target labels 2–6 only.

Both results passed independent Python standard-library/Fraction verifiers, with no `numpy`, `sympy`, numerical geometry, or search-module imports. The checks include complete orbit/root denominators, strict timelike geometry, and rejection of forged polynomial and algebraic ridge points.

## Rational sections are not excluded merely by arithmetic dimension bounds

For a rational section normal, the projected wall normals and the quadratic form are rational. A finite-volume Coxeter group produced this way would be **quasi-arithmetic over Q**. Arithmeticity additionally requires integrality of Vinberg's cyclic products; rationality alone does not supply it. Bogachev–Yorov explicitly exhibit noncompact properly quasi-arithmetic examples over Q in Tables 5–6, and distinguish the integrality condition in Theorem 2.2. Their Theorem 1.5 even gives a compact properly quasi-arithmetic example over an isotropic rational form in H3. [*Quasi-arithmetic hyperbolic Coxeter prisms*](https://nvbogachev.github.io/Quasi_arithmetic_Coxeter_prisms-2024.pdf).

Thus the known exclusion of noncompact **arithmetic** H20 reflection groups does not invalidate a rational, potentially nonintegral section search. No candidate below was rejected merely for being rational. The angle-targeted search also explicitly includes nonrational coefficients.

## Input, symmetry denominator, and finite volume

The exact input remains [the 210 Borcherds walls](../data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json), with quadratic-form sign reversed so spacelike vectors have positive norm. Their primitive sum is the strict center

\[
 h^2=-14,\qquad (h,r_i)=-1\quad\text{or}\quad-8
\]

for the 42 short or 168 long walls respectively. Set

\[
 u_i=\frac{8r_i}{-(h,r_i)},\qquad (h,u_i)=-8.
\]

The exact ambient isometries saved in the previous round induce **49 orbits on all 1,521,520 unordered triples**. This count is independently reconstructed by orbit traversal; it is not inferred from an asserted symmetry-group order.

For each orbit, choose which of its three walls has the negative coefficient. The balanced calibration `e=u_i+u_j-2u_k` gives 147 sections, all with exact actual forbidden ridges. The continuous family is

\[
 e(t)=u_i+t u_j-(1+t)u_k,\qquad t>0.
\]

Every nonzero central three-wall normal is of this form after rescaling, possibly reversing its sign and reordering its support. Interchanging `i,j` replaces `t` by `1/t`, up to overall scale. Consequently **294 ordered half-pencils on `0<t≤1` cover every central three-wall section**. No ambient symmetry interchanging the two positive walls is assumed for this reciprocal reparameterization.

Each section contains `h` strictly; its normal is nonzero and lies in the positive-definite space `h⊥`. It is therefore a full-dimensional H20 section. Finite volume follows from the finite/ideal Klein-vertex section argument proved in the [previous report](fresh_h20_global_20261013.md). The noncentral cases below instead supply explicit strict timelike section centers and use the same volume argument.

## Continuous parameter certificates

The 294 half-pencils are excluded by **441 actual-ridge constraints**: 430 necessary Coxeter-angle polynomials and 11 everywhere-obtuse ridge certificates.

For `q=e²` and `v_a=(r_a,e)`, the matrix

\[
 H_{ab}(t)=q(t)(r_a,r_b)-v_a(t)v_b(t)
\]

is `q` times the Gram matrix of the projected normals. Every saved ridge pair has positive-definite projected two-by-two Gram matrix over the entire parameter interval. The archive gives a polynomial point `p(t)` with

\[
 (p,e)=0,\qquad p^2<0,\qquad(p,h)<0.
\]

Exactly the indicated inequalities vanish. Every other inequality is strictly negative, except explicitly listed redundant zero inequalities whose projected normals have exact nonnegative coefficients in the two ridge normals. Thus the two indicated normals define true local facets and an actual ridge; no raw pairwise Gram defect is used as a geometric rejection.

Strict polynomial signs on `0<t<1` are certified by nonnegative coefficients after

\[
 t=\frac{s}{1+s},\qquad s>0,
\]

clearing the positive denominator. The endpoint `t=1` is checked separately. Nonnegative cone coefficients use the corresponding nonstrict check.

For an acute candidate ridge, write its squared cosine as `N(t)/D(t)`. The certificates prove `0≤N/D<3/4`. Any Coxeter value must therefore satisfy

\[
 N(4N-D)(2N-D)(4N-3D)(16N^2-12ND+D^2)=0.
\]

The extra label-6 endpoint and the conjugate label-5 root enlarge this necessary zero set harmlessly. The polynomial gcd of all saved constraints has **no root in `0<t≤1`** for every half-pencil. Gcds need not be constant: the independent verifier uses an integer pseudo-remainder Sturm implementation to prove absence of roots in the interval. An everywhere-obtuse actual ridge instead contributes the impossible condition `1=0`.

Two initially remaining algebraic possibilities, one quartic and one octic for the support `{1,20,190}`, were used to select additional actual-ridge constraints. The final certificates close those whole half-pencils; there are no unresolved parameter intervals in the central theorem.

## A useful geometric construction: interpolate ambient face points

Projecting `h` orthogonally onto the candidate ridge plane can leave the chamber even when the ridge really exists. Eight saved certificates resolve that issue by an exact geometric construction.

For a ridge pair containing the positive support wall `i`, take strict relative-interior timelike points `z_1,z_2` of the actual ambient triple faces `(i,b,j)` and `(i,b,k)`. Their section-normal pairings have opposite signs throughout the parameter interval. Then

\[
 p(t)=(e(t),z_1)z_2-(e(t),z_2)z_1
\]

is a positive combination of future timelike vectors, lies in `e(t)⊥`, and vanishes on precisely the two desired ambient walls. The remaining inequalities are strict because the two original triple faces have different third walls. The archive saves the face points and verifies all resulting polynomial identities directly.

This supplies genuine ridges where the simpler orthogonal-projection witness failed; it does not enlarge the search cutoff.

## Noncentral construction: choose algebraic coefficients to repair an angle

Among the 49 triple orbits, 15 are positive definite. They contain 17 choices of a positive pair `i,j` with equal norms and equal pairings with the third wall `k`. We tested the symmetric pencils

\[
 e(t)=r_i+r_j-t r_k,\qquad t>0.
\]

These sections generally do not pass through `h`. The actual ambient triple face persists as a section ridge: the projected normals satisfy `r'_k=(r'_i+r'_j)/t`, so the third inequality is redundant locally. We solved exactly for its angle to be `π/m`, for `m=2,3,4,5,6`.

The independent Sturm census finds exactly:

| Repaired local label | Parameter choices |
|---|---:|
| 2 | 0 |
| 3 | 13 |
| 4 | 17 |
| 5 | 17 |
| 6 | 17 |
| **Total** | **64** |

There are 9 rational, 4 quadratic, and 51 quartic parameters. Roots of the conjugate label-5 angle equation are discarded unless they give the intended squared cosine. Completeness is independently checked for each declared pencil and label.

For example, on a mutually orthogonal short triple,

\[
 \cos\theta=\frac1{1+t^2}.
\]

Thus labels 3, 4, and 5 require respectively

\[
 t=1,\qquad t=\sqrt{\sqrt2-1},\qquad
 t=\sqrt{\sqrt5-2}.
\]

These are deliberate exact angle repairs, not approximate numerical tuning. For the first archived rational case, `e=r_2+r_3-r_1` repairs the persistent local ridge to label 3, but another actual ridge, at ambient walls 2 and 4, has positively paired outward normals with scaled Gram

\[
 \begin{pmatrix}8&2\\2&11\end{pmatrix}.
\]

Its exact timelike witness has square `−3828384`, so its obtuse dihedral angle conclusively rejects this candidate.

All 64 candidates have strict timelike interior certificates: 57 use the projection of `h`, and 7 use a small inward perturbation of the persistent triple-face point. Every candidate has another exact forbidden actual ridge. No candidate survives as a finite-volume Coxeter example. No claim is made about higher target labels, unequal positive coefficient ratios, or other noncentral triples.

## Artifacts and independent replay

- [49 triple orbits](../data/enumerations/fresh_h20_triple_orbits_20261013.json).
- [147 balanced calibration sections](../data/enumerations/fresh_h20_triple_central_20261013.json); [producer](../scripts/fresh_h20_triple_central_20261013.py).
- [294 complete continuous half-pencil certificates](../data/enumerations/fresh_h20_triple_pencils_20261013.json); [producer](../scripts/fresh_h20_triple_pencils_20261013.py).
- [64 noncentral algebraic angle repairs](../data/enumerations/fresh_h20_triple_angle_targets_20261013.json); [producer](../scripts/fresh_h20_triple_angle_targets_20261013.py).
- [Independent central verifier](../scripts/fresh_h20_triple_verify_20261013.py); [PASS record](../data/enumerations/fresh_h20_triple_verify_20261013.json).
- [Independent algebraic verifier](../scripts/fresh_h20_triple_angle_verify_20261013.py); [PASS record](../data/enumerations/fresh_h20_triple_angle_verify_20261013.json).

```sh
python scripts/fresh_h20_triple_verify_20261013.py
python scripts/fresh_h20_triple_angle_verify_20261013.py
```

The verifiers use only the standard library and exact rational polynomial arithmetic. They reconstruct the complete orbit and positive-root denominators, verify the strict centers and actual ridge geometry, and reject deliberately changed ridge-point coordinates. Seven Sturm calibration checks cover multiple roots, sign changes, irrational roots, and the excluded zero-parameter endpoint. All producer and verifier scripts compile.

A separate proof-and-code review confirmed the half-pencil coverage, reciprocal reparameterization without a symmetry assumption, actual-ridge logic, endpoint handling, and finite-volume section argument. It also checked a subtle cone-identity step: the residual vector lies in the positive-definite orthogonal complement of the certified timelike point, so its verified zero norm forces the full vector identity. The review checked the archived counts and reran the seven Sturm regressions; it did not repeat the large producer search.

The next available section directions are genuinely noncentral three-wall normals outside the tested symmetric coefficient pencils, or central normals requiring at least four walls. The present results establish neither a general H20 obstruction nor a new example.
