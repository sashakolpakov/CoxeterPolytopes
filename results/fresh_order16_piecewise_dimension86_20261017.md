# Piecewise ideal payments: dimension at most86

**Complete internally checked computer-assisted proof.** This argument
gives dimension at most **86** for every finite-sided,
finite-volume real hyperbolic Coxeter polytope, with no arithmetic,
compactness, simple-cusp, cusp-rank or facet-count hypothesis. The
compact-facet ambient ceiling79, compact-arithmetic ceiling21 and
arbitrary compact ceiling29 are unchanged.

The [source](../scripts/fresh_order16_piecewise_dimension86_20261017.py)
and [certificate](fresh_order16_piecewise_dimension86_20261017.json)
combine the [new local payment theorem](fresh_order16_finite_cross_payments_20261017.md)
with the completed [dimension91 proof](fresh_ideal7_cross_dimension91_20261017.md).
The [separate exact review](fresh_order16_piecewise_dimension86_20261017_review.md)
records independent scalar reconstruction and its written composition scope.
Verification status is stated below; no construction or sharpness is claimed.

## Actual face payments and costs

The local theorem retains the preceding unit finite payments and adds
finite half-payments between whole A1/D15, A2/D13 and A1/A15 components,
at specified fork marks or the unique middle of A15. The preceding
A1/D16 and A2/D14 half-payments remain. Every triangle receives at
least2 and every quadrilateral at least1. The total finite-vertex cost
remains at most9n for n<=108.

These payments reduce the order16 ideal costs to112,107,103,111 for
affine A,B,C,D. Complete induced embeddings cover all260 potentially
relevant preferred components from the inherited classification. The
cycle replacement and nine exceptional D stars are proved universally,
including arbitrary additional supporting whites; they do not depend
on seeing the same component types in one known small chamber.

For an affine factor of rank r, the pair count divided by r+1 is at most

    b(r) = r/2       if 1<=r<=14,
           7        if r=15,
           127/17   if r=16,
           217/29   if r>=17.

The large-rank envelope follows from the explicit all-rank affine
spine/fork formulas. It is not inferred from a finite sample.
At a cusp let d=n-1, V=product_i(r_i+1), and
h(t)=product_i(1+t+...+t^(r_i)). The inherited
[factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md)
gives total ideal cost at most V sum_i b(r_i), including at nonsimple cusps.

For G=sum_j h_j min(j,d-j), the normalized curvature from the
[shape inequality](fresh_cusp_link_shape_curvature_20261017.md) is
(5G-dV)/(2V). Its finite coefficient is

    beta_n=n(n-6)/8              for even n,
    beta_n=(n-1)(n-5)/8          for odd n.

The shape inequality bounds D5, the sum of5 minus the number of sides
over all two-faces, below by these vertex coefficients plus a strictly
positive global slack. Actual corner payments bound D5 above by the
finite and ideal costs. Comparing coefficients therefore suffices.

## Complete concentration for the piecewise budget

The inherited [all-rank convex-order theorem](fresh_distance14_unrestricted_dimension156_20261017.md)
states that for1<=a<=b, replacing(a,b) by(a-1,b+1) increases the spread
of the centered uniform sum and can only decrease normalized curvature.
It applies with the remaining factors fixed. A zero rank is omitted.

Use only these transfers:

1. For two ranks below14, concentrate them until one is0 or the other14.
   The budget r/2 is linear throughout, so total cost is preserved.
2. Replace two ranks15 by14 and16. Their budget increases by
   7+127/17-14=8/17, while curvature decreases. Thus curvature minus
   cost can only decrease; cost preservation is not required.
3. For two ranks greater than17, concentrate until one is17. Their
   budget is the constant217/29, so total cost is preserved.

Each nontrivial transfer strictly increases the integer sum of squared
ranks at fixed total d, so iteration terminates. The resulting shapes
have any number of14,16,17 parts, at most one smaller part in1,...,13,
at most one15, and at most one larger part of size at least18. These
shapes form a covering reduction of every positive integer partition;
it is harmless if some shapes are not themselves final outcomes of a
particular transfer sequence.

At every step curvature-minus-cost does not increase. Positive margin
for every reduced shape therefore proves positive margin for every
original cusp-rank partition. This argument explicitly handles the
nonlinear ranks15 and16; applying the previous single-cap reduction
without this adjustment would not be justified.

## Exact comparison and conclusion

The source checks all3355 reduced profiles in link dimensions86 through90:

| Ambient dimension | Profiles | Minimum curvature minus cost | Finite slack beta_n-9n |
|---:|---:|---:|---:|
|87|621|146417/911250|197/2|
|88|646|203099/607500|110|
|89|671|136309/253125|123|
|90|696|34528/50625|135|
|91|721|606679/708750|297/2|

The smallest ideal margin is146417/911250, at ranks(2,14,14,14,14,14,14),
ambient87. All ideal margins and finite slacks are strictly positive.
The resulting D5 lower bound exceeds its payment upper bound, excluding
dimensions87 through91. The completed91 proof excludes all larger
dimensions. Hence the ceiling is86.

At ambient86, the profile(1,14,14,14,14,14,14) has margin
-112409/1518750. This comparison therefore cannot yield a ceiling85.
The formal link profile is not a realized finite-volume polytope.

## Verification scope

The primary calculation uses exact integer convolution and a separate
inclusion-exclusion moment identity, compares all3355 profiles, checks
the finite coefficients and records the boundary obstruction. The
separate checker uses a different enumeration, sliding-window
convolution and a cumulative-tail identity without project arithmetic
imports. Primary generation and a separate fresh saved replay both
terminate exit0/PASS for all3355 profiles. The independent comparison
agrees on every polynomial, G, curvature, cost and margin, all five
dimension minima and finite slacks, and the boundary. Its linked report
records the final saved bindings and replay status. Separate full
mathematical/source reviews confirm the piecewise concentration,
actual-payment composition and completed91 tail, with no identified gap.

The local theorem has a separate full mathematical/source review,
including its complete embedding domain, universal support arguments,
finite budget and all-rank ideal envelope. The independent scalar
checker does not duplicate that embedding census or the inherited
complete geometric classification. External refereeing, publication
priority and proof-assistant formalization are not established here.
