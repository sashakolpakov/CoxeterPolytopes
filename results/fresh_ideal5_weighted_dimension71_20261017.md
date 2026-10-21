# Unrestricted finite-volume dimension at most71

The order13 transfers and the weighted order17 rule, together with the exact
global comparison, give the following computer-assisted conclusion:

**Every finite-sided, finite-volume real hyperbolic Coxeter polytope has
dimension at most71.**

There is no arithmetic, compactness, simple-cusp, cusp-factor or facet-count
restriction. The [source](../scripts/fresh_ideal5_weighted_dimension71_20261017.py)
and [certificate](fresh_ideal5_weighted_dimension71_20261017.json) bind the new
local interfaces and verify all4,913 reduced cusp profiles in dimensions72,
73 and74. The completed [dimension74 theorem](fresh_ideal5_piecewise_dimension74_20261017.md)
handles every higher dimension. Verification status is recorded below.

## The complete modified face rule

Retain the complete finite8/ideal5 face theorem and its validated successors,
with these changes. At finite vertices use the
[joint coupon budget](fresh_finite8_d12_b10_coupon_budget_20261017.md): the old
half between whole A1 and a short leaf of whole D12 becomes one, and a new
half joins whole A2 to the doubled terminal of whole B10. All earlier rules,
including A4/D13-neighbor halves, remain. The total finite cost is at most8n
through dimension75. This dimension range covers all new comparisons here.

At ideal factors of order at most12, all pairs remain paid. At order13 the
four classical affine costs are now A65, B75, C73 and D74. At orders14,15,16,
the prior rules and costs are unchanged. At order17, twelve affine-D16 pair
weights fall from one to one half, reducing its cost126 to120. Other order17
types are unchanged. The larger-order ideal5/R5 rule is unchanged.

The essential coverage statement remains: every triangle receives at least2
and every quadrilateral at least1. Here is the full interface argument,
including the short-path and immediate-factor obligations.

For order13 cycles, the [universal cycle13 transfer](fresh_cycle13_ideal5_transfer_20261017.md)
removes exactly the13 distance-six pairs. Its remaining65 pairs have distance
at most5. The missing pairs occur in the short8-first stage of the original
proof, so absence from selected overlap maps would not suffice. The new
universal transfer supplies a finite short-path payment or an ideal short
path. Its final isolated profile requires energy10/13, excluded by the
complete strict spectator spectrum. It needs no new finite coupon.

For affine B12 of order13, remove the singleton(0,10) and the two distance11
pairs(10,11),(10,12), keeping the other75 pairs. The
[B12 singleton proof](fresh_order13_affine_b_singleton_transfer_20261017.md)
seals the two singleton cores and supplies the two finite halves, one old
and one new A2/B10 half. Its four saved conditional distance11 contexts
are promoted here:166 and1744 have an existing whole finite B13 unit;
170 and2137 have whole A1 against a short leaf of whole D12 and receive
the newly budgeted full coupon. The source checks the actual whole types
and marked roles in those saved finite corners. Their sealing is the same
positive-affine-kernel and full-core radical argument as in the singleton
proof. The B12 proof also supplies a strict positive determinant gap for
all universal D14 `inner4` alternatives, so none of those formerly immediate
payments is silently lost.

For affine D12 of order13, the
[universal opposite-fork-leaf transfer](fresh_affine_d12_leaf_pair_transfer_20261017.md)
removes all four distance10 leaf pairs, leaving74. This covers all127
affected selected contexts and all universal D14 `inner_next3` leaf-pair
alternatives. It uses the new full A1/D12 coupon where needed; its two-sided
boundary instead gives a short path or two existing halves.

For affine C12, retain the73-pair basic-plus-selected set in the
[complete order13 support certificate](fresh_order13_payment_support_20261017.md).
Its five omitted pairs have distance at least10. Every elementary triangle
pair and every pair of distance at most8 remains a unit. Every selected
marked image is retained, and there is no C-type immediate order13 demand.
Thus this change needs no new geometric local transfer. The primary checker
checks all these inclusions for the other three types, the complete801
eligible-pattern domain including larger component bounds, and all six
literal D14 immediate records.

These exhaust the new order13 obligations in the complete face theorem.
Its low-core alternatives have whole component order at most12. The D13
universal transfer has immediate factor order12. Its D14 transfer has the
order13 B/D alternatives just discharged. The other sealed order16/17
transfers end in finite payments or retained short paths. The general
feature-free tree replacement ends in finite or retained feature-path
payments. There is no additional whole-order13 assumption being imported.

At order17 use the [twelve half-weight transfers](fresh_d17_partial_weight_transfer_20261017.md).
All affected pairs have distance greater than9, so short8 and triangle units
are unchanged. The complete saved marked-image support consists of twelve
sealed shared-L25 contexts. In each, the selected ideal corner keeps one
half and a different actual finite corner has an existing half. Thus the
face still receives at least1. This is a sum of weights, not a claim that
a half is a full payment. The prior zero-weight D17 singleton retains its
two-finite-half replacement. No finite cost is added by these twelve changes.

The original complete low/disjoint/overlapping core classification and
induced-map censuses remain explicit frozen inputs. The modifications above
check every place their previous ideal payments change. All finite coupons
apply to whole actual finite components, not merely finite subdiagrams
inside ideal corners.

## The new ideal envelope and concentration domain

The resulting normalized ideal-factor costs are bounded as follows:

|Affine rank r|Normalized cost|
|---|---:|
|1 through11|r/2|
|12|75/13|
|13|87/14|
|14|97/15|
|15|111/16|
|16|120/17|
|17 and larger|163/23|

The full small-order maxima are taken over all four types. In particular
order17 has A85, B111, C108 and D120. The all-rank tail is unchanged.
For an actual nonsimple cusp with factor ranks r_i and V link vertices,
each factor pair-weight sum is multiplied by V/(r_i+1), exactly as in the
inherited factor-fiber theorem. Fractional weights preserve this linear
counting formula.

For total cusp rank d=n-1 write

    H(t)=product_i(1+t+...+t^r_i), V=H(1),
    G=sum_j [t^j]H(t)*min(j,d-j).

The normalized curvature term is (5G-dV)/(2V). Subtracting the sum of the
costs in the table gives the ideal margin.

The inherited convex-order lemma decreases normalized curvature when
(a,b) with1<=a<=b is concentrated to(a-1,b+1). We use only moves that also
preserve or increase total cost. Low ranks concentrate up to11, leaving at
most one rank from1 through10. Two ranks12 may become11 and13: their cost
increases by exactly16/91. Two ranks at least18 concentrate until one is17,
preserving the constant tail cost. Each move strictly increases the bounded
integer sum of squared parts at fixed total rank, so the process terminates.

Every partition therefore reduces to repeated11,13,14,15,16,17; at most
one12; at most one low rank1 through10; and at most one rank at least18.
The older move eliminating repeated15 is not used: its cost difference
has become negative under the new rank16 coefficient. Arbitrarily many
rank15 factors are explicitly included in the new enumeration.

## Exact global comparison

Every reduced profile in each of the three dimensions has positive ideal
margin. The checker evaluates exact rational values throughout:

|Ambient dimension|Reduced profiles|Minimum ideal margin|Finite slack beta_n-8n|
|---|---:|---:|---:|
|72|1530|6875/49152|18|
|73|1635|41401/145152|28|
|74|1748|135455/304128|37|

The respective minimizing profiles are(11,15,15,15,15),
(11,11,11,11,13,15), and(10,11,11,11,15,15). Here beta_n is n(n-6)/8 for
even n and(n-1)(n-5)/8 for odd n. The inherited global shape inequality also
has the strictly positive constant5(n-1). Summing the retained face-payment
inequalities and comparing with the finite and ideal vertex budgets gives
a contradiction in each of dimensions72,73,74. The completed74 theorem
excludes all dimensions at least75. Together they prove the ceiling71.

At ambient dimension71, the profile(11,11,11,11,11,15) has margin
-5597/248832. Thus this envelope does not prove a bound70. This is a scalar
obstruction, not a polytope construction or a claim that such a cusp is
realized.

The separate compact arithmetic bound21 and arbitrary compact bound29 are
unchanged. The effective ambient ceiling for the compact-facet alternative
becomes71; its older auxiliary proof remains79. No new finite-volume
polytope or new example dimension is claimed.

## Verification status

Primary generation and a fresh separate-process saved replay both terminated
exit0/PASS on all4,913 profiles and the new local interfaces. The separate
[standard-library reconstruction](fresh_ideal5_weighted_dimension71_20261017_review.md)
agrees with every exact profile row, minimum and boundary. It uses a different
partition enumeration, sliding polynomial convolution, cumulative-G formula
and absolute-deviation identity. It also independently reconstructs all13,034
finite-budget endpoints and the order13/weighted17 interfaces.

Two full internal source and mathematical reads of the assembled argument
found no identified gap in the complete local provenance, whole-component
coupons, weighted face sums, actual cusp multiplicities, new concentration
domain, finite budget or inherited higher-dimensional tail. Every local input
has generation, fresh replay and internal mathematical/source reviews.
All31 direct primary source bindings, Python syntax and local links pass.
The complete inherited overlap and induced-map censuses retain their frozen
validation records and are not rerun here. Source, data and report are
final/frozen. All reviews are internal to the agent workflow. External
mathematician refereeing and Lean or other proof-assistant formalization
are not claimed.

    python3 scripts/fresh_ideal5_weighted_dimension71_20261017.py
    python3 scripts/fresh_ideal5_weighted_dimension71_20261017.py --verify
