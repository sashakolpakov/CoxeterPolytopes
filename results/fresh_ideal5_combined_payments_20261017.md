# Complete finite9 and ideal5 face payments

For every finite-sided finite-volume hyperbolic Coxeter polytope of dimension
n≤80, the rules below give each triangular two-face total corner payment at
least2 and each quadrilateral at least1. Finite-vertex payment is at most9n.
At an ideal vertex the exact factor costs have the piecewise envelope stated
below. This is a complete local-payment and counting composition, including
nonsimple cusps and arbitrary supporting white roots. The ambient dimension
comparison is separate.

The [checker](../scripts/fresh_ideal5_combined_payments_20261017.py) and
[certificate](fresh_ideal5_combined_payments_20261017.json) bind the frozen
proofs and check their interfaces. The complete core classification used
here is the finite9 classification. No claim of complete finite8 overlap
coverage is needed or made.

## The rules

At a finite vertex, pay1 when the marks lie in the same whole irreducible
component and either their graph distance is at most9 or that component has
order at most17. In addition, pay one half between distinct whole components
when one has rank at most5 and the other is one of:

* D11 through D16, marked at its fork center or either short fork leaf;
* A13 or A15, marked at its unique middle;
* A17, marked five steps from one endpoint and eleven from the other.

Each qualifying pair receives one half, with no duplicate coupon for multiple
descriptions. The unit and coupon conditions use the full actual finite
vertex diagram, not a selected subdiagram.

At an ideal vertex, use the full varying affine factor containing the two
marks, as supplied by the [actual factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md).
Pay every pair when its order is at most13. At orders14–17 use the exact
automorphism-invariant pair sets in the
[571-pattern census](fresh_ideal5_small_factor_payments_20261017.md).
At larger orders pay distance at most5, and distances6–9 if an endpoint is
within5 of a branch or a nonsimple-bond endpoint. All these ideal payments
are zero or one.

## Triangle coverage is retained

The [complete finite8/ideal5 triangle theorem](fresh_ideal5_finite8_triangle_payments_20261017.md)
already pays every triangle twice by finite distance at most8, ideal factors
of order at most9, or ideal distance at most5 and feature-near distances at
most8. The present finite rule contains that finite rule. The source checks
the inclusion of every elementary triangle pair in every order14–17 target.
Conceptually, cycles have no features, so their removed distances greater
than5 were never elementary triangle payments. The nine removed D pairs
have minimum distance9 at order14 and larger minimum distances in higher
orders, whereas the triangle theorem's feature-near paths have length at
most8. Factors of order at most13 are wholly retained. Thus no replacement
of one triangle corner by only one total quadrilateral payment is used.

## Quadrilateral coverage from the complete finite9 domain

Start with an actual quadrilateral. If it has a finite marked corner path of
length at most9, its unit payment remains. Suppose an ideal corner has such
a path. It pays immediately at distance at most5 or within the radius5
feature neighborhood, unless it is one of the explicitly removed pairs.

For a classical tree factor outside that neighborhood, the complete
[finite8/ideal5 local theorem](fresh_finite8_ideal5_local_replacement_20261017.md)
replaces marked distances6–9 by retained units or sealed finite coupons.
The same theorem treats cycles of order at least15, including the sealed
A1/A15 and A1/A17 equalities. The separate
[cycle14 theorem](fresh_cycle14_ideal5_transfer_20261017.md) treats its
remaining distances6 and7 using the strict spectator spectrum and the
sealed A13-middle coupons. Smaller factors are wholly paid.

All D pairs deliberately deleted from the basic rule are handled by the
specific transfers below. Therefore, when no corner has yet paid and no
transfer has already supplied a total unit, there can be no adjacent marked
path of length at most9. This is exactly the original hypothesis of the
complete [dimension108 core classification](fresh_selective_distance9_dimension108_20261017.md),
not a silently weakened no-short8 hypothesis.

That classification includes its complete all-label low-core supplement
and949 high-core overlap contexts. The low-core selected component bounds
are at most12. At a finite selected corner, the retained small17 rule pays.
At an ideal selected corner, the whole varying factor has the proved bound
at most17. Orders at most13 pay. For orders14–17 the actual marked principal
component gives an induced labelled embedding into the full factor; all
such maps are in the571-pattern census, with its extension-bound filter.
Thus its pair is in the union before the certified deletions. No arbitrary
provenance, favorable map, or minimum across incompatible realizations is
selected.

## Every deleted pair has a retained alternative

The cycle transfers above directly use the present short/feature rule and
the listed whole finite coupons. They do not leave an unexamined ideal
distance6 or7 obligation behind.

For order14 affine D13, the
[universal long-pair transfer](fresh_affine_d13_long_pair_transfer_20261017.md)
replaces all nine cross-fork pairs. Its immediate ideal alternatives are
paths of length at most3 or the two literal affine components of order12;
these all pay. Its other conclusions are sealed whole finite components or
small components of ranks1,2,4,5 paired with D11 or D13 fork marks, all in
the finite rule.

For order15 affine D14, the
[preceding universal transfer](fresh_affine_d14_long_pair_transfer_20261017.md)
has only two immediate whole-small ideal alternatives, both of order13.
They arise from the center deletion D12 plus the new root, with either an
endpoint4 or a second-node3 attachment. They are wholly paid. A possible
order14 affine component after another outer attachment is **not** an
immediate paid alternative in that proof: it remains among the bounded
nonextreme profiles and is eliminated by the strict determinant gap. The
only short ideal alternatives used are the adjacent-product or common
spectator paths of length1 or2, not an arbitrary distance6 pair. Its sealed
small-rank≤3 against D12/D14 coupons are included.

For orders16 and17, the nine affected preferred-core contexts at each order
are matched exactly to the saved
[order16](fresh_order16_finite_cross_payments_20261017.md) and
[order17](fresh_ideal7_finite_cross_payments_20261017.md) replacements.
Those are literal sealed-core arguments, not a new universal center-profile
extension. A nonorthogonal supporting white gives a retained length2 path;
otherwise the saved finite unit/coupon corners sum to at least1. The source
checks that every nonzero terminal payment in all18 cases is finite. There
is consequently no hidden order14 ideal payment inherited from their older
global rules. Their core matching is applied only after the original
no-short9 reduction has been restored.

These implications complete quadrilateral coverage. They are direct
contradictions to a face whose total new payment is less than1: a direct new
ideal unit, a direct finite unit, or two sealed halves resolves it. No cyclic
chain of ideal transfers is assumed to terminate without proof.

## Finite and ideal budgets

The [finite8 coupon budget](fresh_finite8_a17_d11_coupon_budget_20261017.md)
already covers exactly the stated coupons together with distance≤8,
distance9 at a whole-component leaf, and whole finite components of order
at most17, for n≤80. Its total is at most8n.

To restore all finite distance9 pairs, nothing is added in ranks at most17.
At rank r≥18 the only irreducible finite diagrams are A, B and D. Each has
r−9 pairs at distance9. Exactly two of these pairs contain a leaf for A/B,
and exactly three for D; there are no leaf-to-leaf distance9 exceptions at
these ranks. The additional costs are r−11 and r−12, respectively, each at
most r. Summing over components adds at most n. The new finite total is
therefore at most9n with all enlarged coupons included. The certificate
checks the literal pair sets at every relevant rank18–80, alongside this
rank-independent path count.

The small-factor costs and the all-rank ideal5/R5 spine proof give the factor
rank envelope

| Rank r | Normalized factor cost bound b(r) |
|---:|---:|
|1–12|r/2|
|13|87/14|
|14|20/3|
|15|111/16|
|16|127/17|
|at least17|163/23|

For an actual cusp whose link is a product of simplices of ranks r_i, let
V be its number of incident edges. Each pair in a varying factor of order
r_i+1 repeats V/(r_i+1) times. Hence the ideal-vertex payment is at most
V times the sum of b(r_i). The actual link and factor-fiber lemmas justify
this multiplicity even at nonsimple ideal vertices; no global simple-cusp
or arithmetic assumption has entered.

## Verification scope

The checker binds the complete inherited core and low-core proofs, the1025
triangle cases, the571-pattern target sets, both universal D transfers,
both nine-case sealed replacements, the cycle/weak-profile theorems and the
finite coupon budget. It checks all16 triangle-pair inclusions, both literal
order12/order13 immediate alternatives, terminal finite payment totals,
the small costs and every finite9 incremental pair set. The all-label
classification, local geometric lemmas and arbitrary-rank counting proofs
remain explicit inherited mathematical inputs. Generation and a fresh
separate-process saved replay both terminated exit0/PASS. A separate full
source/mathematical read found no gap in the short-path reduction before
the949-case classification, triangle-pair retention, immediate order12/13
alternatives, literal finite replacements, finite9 incremental cost or
actual cusp multiplicities. That review did not rerun the571-map census,
the949 core census or the saved row arithmetic. Syntax, all38 direct
bindings and report links pass. Source, data and report are final/frozen.
No external refereeing or proof-assistant formalization is claimed.
