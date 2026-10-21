# Noncentral two-parameter H20 sections

No new Coxeter polytope was found. The new exact result closes a continuous family left open by the central-section theorem:

> Every interior H20 section of the frozen Borcherds H21 chamber whose normal lies in an **elliptic span of at most three original wall normals** is non-Coxeter.

There is no requirement that the section pass through the old center. The proof covers all real coefficient ratios, including nonrational ones. It covers **266,000 unordered elliptic triples in 15 symmetry orbits**, with 45 choices of the negative coefficient. Affine and Lorentzian three-wall spans, and larger supports, are outside this theorem. Two-wall boundary cases use the [earlier separate obstruction](fresh_h20_global_20261013.md).

The [independent standard-library verifier](../scripts/fresh_h20_bivariate_verify_20261013.py) passes. It reconstructs all 49 triple orbits, selects the 15 elliptic ones, and checks **58 actual-ridge certificates**, including an algebraic construction that repairs two chosen angles simultaneously. It rejects both a forged face point and a forged coefficient parameter.

## A genuine two-angle repair

Use the original one-based wall indices and set

\[
 e=r_1+s r_{48}-t r_{20},\qquad t=1-2s,
\]

where `s` is the unique root in `(0,1/2)` of

\[
 23s^4-52s^3+44s^2-16s+2=0.
\]

Numerically `s≈0.26275747346093226` and `t≈0.4744850530781355`. These approximations merely identify the archived rational isolating interval. The quartic is irreducible modulo 3, so the coefficient has degree four.

This normal comes from **two independent exact angle equations**. The actual ridge at walls 48/22 becomes right-angled when `2s+t=1`. On that line, the actual ridge inherited from the ambient face `{1,48,20}` has squared cosine

\[
 \frac{6s^2-8s+3}{2(5s^2-6s+2)}.
\]

Setting this to `cos²(π/8)` gives the quartic above; the selected real root gives the intended angle, not its conjugate. Thus the section has actual repaired angles **π/2 and π/8**, unequal positive coefficients, and does not pass through the original center. An exact strict timelike interior point proves that the section is full-dimensional. Finite volume is inherited from the finite-volume chamber by the Klein-section argument used in the [previous report](fresh_h20_global_20261013.md).

The resulting section is nevertheless non-Coxeter: the actual ridge at walls 1/19 has

\[
 0<\cos^2\theta
 =\frac{(2s-1)^2}{2(5s^2-6s+2)}<\frac14.
\]

Its acute angle lies strictly between π/3 and π/2. The archive supplies an exact timelike ridge point, strictly inside every other ambient wall inequality. This is a geometric rejection, not merely a forbidden raw Gram entry.

## Continuous coverage and actual ridges

An interior section cannot have all nonzero coefficients of its normal with the same sign: all ambient wall inequalities are strict at an interior point. After reversing and rescaling the normal, each three-wall case is therefore

\[
 e=r_i+s r_j-t r_k,\qquad s,t>0.
\]

The two parameters are independent. Choosing `i<j` does not restrict their ratio. Every elliptic triple has a strict timelike ambient face point. Perturbing that point in its positive-definite normal span gives a strict section interior: prescribe its three wall pairings proportional to `(-t,-t,-(1+s))`, which pairs to zero with `e`, and take a sufficiently small perturbation. This also explains why every parameter in the declared family gives a genuine full-dimensional finite-volume section.

For an auxiliary wall `l`, let `a` be one of the two positive-coefficient walls and let `b` be the other. The archive supplies strict ambient face points `z₁` on `{a,b,l}` and `z₂` on `{a,k,l}`. Their pairings satisfy

\[
 (e,z_1)>0,\qquad(e,z_2)<0.
\]

Consequently

\[
 p=(e,z_1)z_2-(e,z_2)z_1
\]

is a positive combination of future timelike vectors, lies in the section, and vanishes **only** on walls `a,l`. The section normal cannot belong to `span(r_a,r_l)`, since its pairing with `z₁` is nonzero. Hence the projected pair is independent and positive definite in `p⊥`; it defines two true facets meeting at an actual ridge. The verifier checks the point coordinates, square, time orientation, section identity and all 210 ambient inequalities exactly.

Write `v_a=(e,r_a)` and `H_ab=e²(r_a,r_b)-v_a v_b`, so `H` is a positive multiple of the projected Gram matrix. The 45 families divide into three elementary proof types:

| Families | Exact obstruction for every positive parameter pair |
|---:|---|
| 19 | One persistent actual ridge is always obtuse. |
| 14 | One persistent actual ridge always has squared cosine strictly between 0 and 1/4. |
| 11 | One of two persistent actual ridges is obtuse, depending on the signs of `v_i,v_j`. |
| 1 | Two sign regions have an obtuse ridge; their shared boundary has a forbidden acute ridge. |

The first 33 cases use sparse bivariate polynomials with nonnegative coefficients to prove the strict signs on the entire positive quadrant.

For each of the next 11, there is an auxiliary wall attached only to the other positive wall, and orthogonal to the chosen wall and the negative wall. Its section-normal pairing is strictly negative. Thus its ridge with wall `a` is obtuse whenever `v_a>0`. At least one of `v_i,v_j` is positive, because

\[
 v_i+s v_j
 =r_i^2+2s(r_i,r_j)+s^2r_j^2
   -t\big((r_i,r_k)+s(r_j,r_k)\big)>0.
\]

The first three terms form a strictly positive quadratic form, and the remaining terms are nonnegative by ambient acuteness. The verifier checks a completed-square identity and both actual-ridge certificates.

The exceptional ordered triple is `(i,j,k)=(1,48,20)`. Put `B=2s+t-1`. The projected pairings for actual ridges 48/22 and 48/19 are respectively `2B` and `-4tB`; one is positive whenever `B≠0`. On `B=0`, positivity of `s,t` gives exactly `0<s<1/2`. The actual 1/19 ridge then has the squared cosine displayed above. The upper bound follows from

\[
 2(5s^2-6s+2)-4(2s-1)^2=2s(2-3s)>0.
\]

This includes the entire dividing line, not just the quartic two-angle candidate.

## Artifacts and replay

- [Exact universal certificates and quartic candidate](../data/enumerations/fresh_h20_bivariate_elliptic_20261013.json).
- [Producer](../scripts/fresh_h20_bivariate_elliptic_20261013.py).
- [Independent verifier](../scripts/fresh_h20_bivariate_verify_20261013.py) and [PASS record](../data/enumerations/fresh_h20_bivariate_verify_20261013.json).

```sh
python scripts/fresh_h20_bivariate_verify_20261013.py
```

The producer uses symbolic arithmetic to discover and archive the certificates. The verifier imports no producer, SymPy, NumPy or project geometry code; it uses independent rational sparse-polynomial arithmetic and the existing standard-library algebraic-field checker. The angle repair is verified over an isolated real quartic root. No numerical finite-volume or Coxeter claim is made.

This theorem leaves **34 nonelliptic triple orbits** for a separate noncentral search. The absence of an ambient triple face in those cases removes the main geometric mechanism used here. The [separate nonelliptic scout](fresh_h20_bivariate_nonelliptic_20261013.md) constructed two exact quadratic sections repairing π/2 and π/3; both retain actual obtuse ridges. That bounded calculation does not extend this universal theorem to nonelliptic spans.
