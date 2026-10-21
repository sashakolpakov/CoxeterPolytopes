# Independent review of the one-omitted local classification

No issue was found in the completeness argument or the local arithmetic of the [rank-at-least-four classifier](../scripts/fresh_h20_omitted_factor_classification_20261014.py). A fresh primary `--verify` replay passes. In addition, a [separate Fraction/Q(√6) checker](../scripts/fresh_h20_omitted_factor_classification_20261014_local_review.py), importing no producer arithmetic, verifies all 175 accepted records, their 119 spatial-plane classes, and their exact shift minima. Its [saved certificate](fresh_h20_omitted_factor_classification_20261014_local_review.json) includes the 48 irrational records. This separate checker does not repeat the exhaustive support enumeration; that part receives a source/proof review and primary replay.

**Local completeness.** With exactly one local facet omitted from a rank-`n` affine simplex factor, the retained bounded section has dimension `n−1` and `n` actual facets. Its Coxeter Gram is therefore irreducible affine. Deleting the proposed omitted root leaves a positive-definite crystallographic forest `B`. For a unit spatial normal with contact vector `a`, the projected Gram is `H=B−aaᵀ`; the necessary rank and positive-kernel conditions are exactly `aᵀB⁻¹a=1` and a consistently signed, strictly positive vector `B⁻¹a`.

For `n≥4`, the target affine diagram has no triangle, and its finite labels are only 3 and 4. Three contacts of one sign would create a triangle, so each sign class has size at most two. Two of each sign would require every cross pair to be an old edge to avoid an obtuse projected pairing; this would put a four-cycle in the old forest. Hence support size is at most three. If the supported set had no old edge to its nonempty complement, the projected diagram would remain disconnected. A boundary edge therefore exists. Its squared cosine strictly increases, excluding an old label 4 and forcing an old label 3 to become label 4. This fixes the supported contact square to half its old root norm. Overall reversal permits the positive square root.

Every other supported contact then solves the displayed quadratic for one of the complete target cosine values 0, 1/4, and 1/2. Both quadratic branches are enumerated, including repeated roots. The subsequent full projected-pair checks, unit equation, and positive-kernel test remove extraneous squared-equation solutions. Since `B` is positive definite and the unit equation holds, `H` is positive semidefinite of corank one; a strictly positive kernel also precludes separate nonsingular components. There is no coefficient or angle-grid assumption in this reduction. The code uses the full inverse `B⁻¹` restricted to the support, not the inverse of the supported principal submatrix.

**Complete shift interval.** Write `D=Σ m_i r_i`, `e₀=Σ k_i r_i`, and normalize the cusp horosphere by `(D,x)=−1`. On the retained simplex, `y_i=−(r_i,x)≥0` and `Σ k_i y_i=−t`. The omitted inequality is `Σ m_i y_i≤1`. Its maximum occurs at a simplex vertex, so redundancy is equivalent to

\[
-L\le t<0,\qquad L=\min_i\frac{k_i}{m_i}.
\]

The endpoint is not discarded: the exact identity

\[
e_0-LD=\sum_i(k_i-Lm_i)r_i-Lm_o r_o
\]

has nonnegative retained coefficients. If all vanish, the section is the old omitted facet and misses the source interior. Otherwise the omitted inequality is redundant at the endpoint and the endpoint must be checked globally. These facts justify the interval passed to the separate global certificate.

The independent record replay rebuilds every old retained Gram and its positive rational LDL pivots; verifies `Bk=a`, `aᵀk=1`, the ambient unit normal, null direction and positive cusp scale; checks every projected pair, the positive kernel, the exact minimum `L`, and the endpoint identity; and reconstructs the 119 plane equivalence classes from proportional full-factor covectors. Its arithmetic uses the independent two-term basis `1, √6`, with exact sign comparisons.

```sh
python scripts/fresh_h20_omitted_factor_classification_20261014.py --verify
python scripts/fresh_h20_omitted_factor_classification_20261014_local_review.py --verify
```

This review covers the rank-at-least-four local classification and shift intervals. The global actual-ridge interpolation, endpoint exceptions, and separate rank-three classification have distinct reviews and are not claimed as independently replayed here.
