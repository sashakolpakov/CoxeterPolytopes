# Finite cross-component payments and ideal cost seven and a half: dimension91

**Complete internally checked computer-assisted proof.** This argument
gives dimension at most **91** for every finite-sided,
finite-volume real hyperbolic Coxeter polytope. There is no arithmetic,
compactness, simple-cusp, cusp-factor-rank or facet-count hypothesis.
Verification and its limits are recorded below. The compact-facet ambient
ceiling79, compact-arithmetic ceiling21 and arbitrary compact ceiling29
are unchanged.

The [source](../scripts/fresh_ideal7_cross_dimension91_20261017.py) and
[certificate](fresh_ideal7_cross_dimension91_20261017.json) combine the
new [local payment theorem](fresh_ideal7_finite_cross_payments_20261017.md)
with the completed [dimension97 proof](fresh_distance8_ideal_dimension97_20261017.md).
The [independent review](fresh_ideal7_cross_dimension91_20261017_review.md)
reconstructs the scalar arithmetic separately. No construction or
sharpness assertion is involved.

## Payment of actual faces

Retain finite weight1 for connected marked pairs at distance at most9
or in a whole irreducible component of order at most17. Add weight1/2
between A1 and D16, or between A2 and D14, when the D mark is its fork
center or a short fork leaf. All component types refer to whole
components of the actual finite vertex diagram.

At ideal vertices use the intrinsic factor rules of the local theorem:
whole factors of order at most16; selected order17 pairs; and for larger
orders distance at most7 together with distances8/9 near a diagram feature.
The theorem proves that every triangle receives at least2 and every
quadrilateral at least1. Its complete preferred-component embedding
census handles the order17 exception. A Lorentzian radical relation
seals the exceptional19-root stars against outside support, so the two
finite half-payments remain valid in every actual ambient face star.
This step is not an inference from the known H17 example alone.

Thus the sum of all assigned corner weights is at least

    D5 = sum over two-faces F of (5 - number of sides of F).

Larger polygons cause no difficulty because weights are nonnegative.

## Finite and ideal cost bounds

The local theorem proves finite-vertex cost at most9n for every n<=108.
If a,b,c,d count A1,A2,D14,D16 components, its explicit bound is

    9n - 9a - 17b - 35c - 24d + (3/2)ad + 3bc <= 9n,
    a + 2b + 14c + 16d <= n.

The proof allows arbitrary other finite components. The extra finite
payments therefore cost no additional total budget in the dimensions
needed here.

Every irreducible ideal factor of order m costs at most

    min(15m/2, m(m-1)/2).

The proof covers every rank by path, cycle and fork counts; the saved
finite BFS range is a corroboration, not the justification of the tail.
For a cusp with affine ranks r_i, let d=n-1 and V=product_i(r_i+1).
The [actual factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md)
gives multiplicity V/(r_i+1) for each pair in factor i, including at
nonsimple cusps. Its total normalized cost is consequently at most

    sum_i min(15/2, r_i/2).

For h(t)=product_i(1+t+...+t^(r_i)) and G=sum_j h_j min(j,d-j), the
inherited [shape curvature inequality](fresh_cusp_link_shape_curvature_20261017.md)
has normalized ideal coefficient (5G-dV)/(2V). Its finite coefficient is

    beta_n = n(n-6)/8              for even n,
    beta_n = (n-1)(n-5)/8          for odd n.

It gives a lower bound for D5 with a strictly positive global slack.
Therefore coefficients at least the corresponding payment budgets would
already contradict the payment upper bound. The comparisons below are
strict at every vertex type.

## All-rank reduction and exact comparison

The inherited [convex-order transfer](fresh_distance14_unrestricted_dimension156_20261017.md)
applies to an arbitrary integer cap, here15. It preserves the cost
sum_i min(15/2,r_i/2) and can only decrease the curvature when a partition
is concentrated to copies of15, at most one smaller positive part and
at most one larger part. Thus those shapes suffice for all rank partitions.

The source enumerates all477 reduced shapes for link dimensions91 through96:

| Ambient dimension | Reduced shapes | Minimum curvature minus cost | Finite slack beta_n-9n |
|---:|---:|---:|---:|
|92|77|10121/4194304|161|
|93|78|1497773/6291456|175|
|94|79|3493651/8388608|188|
|95|80|1307621/2097152|405/2|
|96|81|3244935/4194304|216|
|97|82|13980081/14680064|231|

The overall minimum is10121/4194304 at ranks(1,15,15,15,15,15,15),
ambient dimension92. Every listed ideal margin and finite slack is
positive, contradicting the two bounds for D5. These dimensions cannot
occur. The completed97 theorem excludes all larger dimensions, so the
new ceiling is91.

At ambient91, the six-rank15 profile has margin-428955/2097152. This
specific scalar comparison therefore cannot prove a ceiling90. The
profile is a formal possible link shape, not a realized finite-volume
polytope or an existence result.

## Verification scope

The producer uses exact integer convolution and checks the curvature a
second way by an inclusion-exclusion moment formula. It binds the local
payment theorem, its mixed replacement input, the completed97 proof and
review, and the factor/shape/convex-order arguments. Its full generation
and fresh saved replay both terminate exit0 with
`PASS_FINITE_CROSS_IDEAL7_DIMENSION91` for all477 profiles.

The separate checker uses no project arithmetic: it enumerates the
reduced shapes differently, reconstructs h by sliding-window convolution,
and compares a cumulative-tail inclusion-exclusion formula for G and an
absolute-deviation curvature identity. Its independent exact comparison
passes every polynomial, G, curvature, cost and margin for all477 profiles,
all six minima, the finite coefficients and the boundary. Separate internal
mathematical/source reviews confirm the local payments and budgets,
actual factor multiplicities, cap15 concentration and completed97 tail,
with no identified gap. The linked review records its saved certificate
and final replay status.
These checks do not rerun the inherited complete core classifications.
External refereeing and proof-assistant formalization are not asserted.
