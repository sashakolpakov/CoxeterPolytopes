# Complete finite9/ideal6 payments and piecewise cusp budgets

This package combines the local arguments into a single payment rule for
arbitrary finite-sided finite-volume hyperbolic Coxeter polytopes. Every
triangle receives at least2 and every quadrilateral at least1. For n<=108
every finite vertex costs at most9n. At an ideal vertex with simplex-factor
ranks r_i and V link vertices, the cost is at most V*sum_i b(r_i), where

    b(r)=r/2 for1<=r<=13;
    b(14)=20/3; b(15)=111/16; b(16)=127/17;
    b(r)=95/13 forr>=17.

The [binder](../scripts/fresh_ideal6_combined_payments_20261017.py) and
[certificate](fresh_ideal6_combined_payments_20261017.json) retain the complete
source dependencies. A global dimension comparison is separate.

## The rule

At finite vertices retain unit payments for same-component marked pairs
of distance at most9 or whole-component order at most17, together with the
[broad finite half-coupons](fresh_broad_finite_half_coupons_20261017.md).
Those coupons apply to whole components, with one of rank at most4 and the
other a specified fork-marked D12,...,D16 or middle-marked A13,A15.

At ideal vertices pay every pair in a varying factor of order at most14.
For larger factors the basic unit rule is distance at most6, or distance7--9
with an endpoint within6 of a branch or nonsimple bond endpoint. At orders
15,16,17 add all induced marked images of the complete inherited preferred
components. Remove the nine D pairs at each of those orders: the two fork
centers and the eight pairs spanning the two fork triples. Their finite
replacements are proved separately, including all supporting whites.

The exact resulting A/B/C/D costs are

| Factor order | A | B | C | D |
|---:|---:|---:|---:|---:|
|15|90|100|95|96|
|16|96|107|103|111|
|17|102|111|108|127|

The [small-factor census](fresh_ideal6_small_factor_frontier_20261017.md)
enumerates every relevant induced labelled map of437 preferred contexts.
The new [direct affine-D14 transfer](fresh_affine_d14_long_pair_transfer_20261017.md)
handles all nine order15 pairs universally, so it handles all63 inherited
contexts affected by their deletion without assuming that their original
small cores are whole actual stars. The frozen
[order16](fresh_order16_finite_cross_payments_20261017.md) and
[order17](fresh_ideal7_finite_cross_payments_20261017.md) arguments provide
the other two nine-pair transfers. Every finite coupon they use is included
in the broader rule.

## Every face is paid

The [triangle theorem](fresh_ideal6_triangle_payments_20261017.md) applies
directly: its short paths and feature witnesses remain paid, and its sole
affine-E8 exception uses a whole factor of order at most9. None of the nine
removed D pairs at orders15--17 belongs to these short/feature payments.

For a quadrilateral, suppose first that some actual corner has a path of
length at most9. A finite corner pays immediately. An ideal factor of order
at most14 pays in full. An unpaid pair in a larger classical factor has
distance7,8 or9 and lies outside the relevant feature neighborhood.
The [mixed6 replacement](fresh_ideal_distance6_mixed_replacement_20261017.md)
handles distance7 cycles and tree distances7--9. For cycles of distance8
or9, the [preceding mixed7 theorem](fresh_ideal_distance7_mixed_replacement_20261017.md)
gives another finite-distance9 or ideal-distance7 payment when the factor
order is at least17. If that ideal-distance7 corner is still unpaid, the
new mixed6 theorem applies there. Since the whole face was assumed unpaid,
each such transfer directly contradicts that assumption; no iteration or
termination premise is needed. The antipodal order16 cycle case is handled
by its [finite A15-center coupons](fresh_cycle16_half_coupon_transfer_20261017.md).
No cycle of order15 has distance8 or9. These cases cover every short path.

Otherwise there is no short9 path at any corner, precisely the hypothesis
of the inherited preferred-component theorem. Its complete small-core
classification gives a preferred component of order at most12; the complete
high-core classification gives one of order at most17. A finite preferred
corner pays by its whole-component rule. At an ideal corner, order at most14
pays in full, and every possible image at orders15--17 was included above.
Only the nine removed D pairs can remain; their universal finite transfer
pays the face. Thus every quadrilateral receives at least1. Polygons with
at least5 sides need no positive lower payment for the D5 comparison.

## Vertex budgets and scope

The broad finite-coupon theorem proves the9n bound for every n<=108, using
only whole-component ranks, slack and at most three eligible large marks.
At factor order m<=14, binom(m,2)/m=(m-1)/2. The table gives the exceptional
normalized maxima20/3,111/16,127/17. For order at least18 the preferred
component bound allows no extra pair. The mixed6 all-rank spine formulas
and their12 linear branches give a common normalized bound95/13; cycles
cost6 per node. These are all-rank arguments, not extrapolation from a
finite graph sample.

The [actual cusp-factor theorem](fresh_cusp_factor_payment_framework_20261017.md)
counts each pair in a rank-r factor with multiplicity V/(r+1). This gives
the stated ideal cost V*sum_i b(r_i), without a simple-cusp, arithmetic or
minimum-factor-rank assumption.

The direct-D14 proof and this composition have received separate internal
mathematical reads. These check the complete large-D extension profiles,
arbitrary supporting spectators, whole-component sealing, the short-path
case split, and retention of every triangle and preferred-component payment.
The binder checks the exact nine removed D14 pairs against the local
certificate, all72 equality payments, all12 small-factor costs and all12
large-factor linear branches. Generation and saved replay are recorded in
the dimension80 verification package and its archive.
No external refereeing or proof-assistant formalization is claimed.
