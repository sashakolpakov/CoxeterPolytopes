# Unrestricted finite-volume dimension at most68

**Every finite-sided, finite-volume real hyperbolic Coxeter polytope has
dimension at most68.** The proof excludes dimension69 by combining the
complete order12 affine transfers with an exact finite-corner budget.
The completed [dimension69 theorem](fresh_ideal5_weighted_dimension69_20261017.md)
excludes every dimension at least70. There is no arithmetic, compactness,
simple-cusp, cusp-factor or facet-count hypothesis.

The [source](../scripts/fresh_ideal5_weighted_dimension68_20261017.py) and
[certificate](fresh_ideal5_weighted_dimension68_20261017.json) bind the local
proofs and calculate the complete new domain of2,060 cusp profiles.
Verification scope is recorded below.

## One simultaneous rule for all faces

Start with the complete finite8/ideal5 rule of the preceding theorem.
The following are the only changes to its weights. Finite type names
denote rank; affine order denotes the number of nodes.

|Factor|Change|New same-factor cost|
|---|---|---:|
|Whole finite A17|Delete its endpoint pair|135|
|Affine A11, order12|Delete the six antipodal pairs|60|
|Affine B11, order12|Give the two distance10 pairs weight1/2|65|
|Affine C11, order12|Delete the endpoint pair at distance11|65|
|Affine D11, order12|Delete the four opposite short-leaf pairs|62|

The finite cross-component rule has three specified changes: upgrade the
existing A1/D11-short-leaf half to one; add a unit joining whole A1 to the
middle node of whole A11; and add halves joining either node of whole A2
to that A11 middle. All other weights and coupons remain those of the
completed predecessor. An upgrade replaces the previous half rather than
adding a second copy.

The [A17 deletion](fresh_finite_a17_endpoint_payment_support_20261017.md)
audits every selected and indirect use of its endpoint unit. Its sole
debited old literal configuration, case662, still has two finite D17
units. All short finite paths and the protected leaf9 pairs remain, as
do the old cross-component coupons. Its later-transfer audit includes the
D13/D14 equality cases and the order16/order17 half-weight rules.

The [cycle12 theorem](fresh_cycle12_ideal5_transfer_20261017.md) is universal:
it covers every quadrilateral with an antipodal marked pair in a whole
affine A11 factor, including short-path and selected-core uses. Its exact
end-profile inequalities and energy-one spectator classification replace
the removed pair by a short ideal path, a finite unit, or the stated A11
coupons. All its finite unit outputs have rank below17 or lie on a path
of length at most6, so the A17 endpoint deletion cannot remove them.

The [D11 theorem](fresh_affine_d11_leaf_pair_transfer_20261017.md) likewise
covers selected, low-core and immediate-factor occurrences. Its all-rank
[spectator lemma](fresh_spectator_five_quarters_transfer_20261017.md)
resolves the energy5/4 boundary. The surviving alternatives are a short
path of length at most5, whole finite units, the full A1/D11-short-leaf
coupon, or two inherited finite halves. The literal A1+A3 spectator cases
have finite corners of rank16; the arbitrary one-sided corners use whole
A1/D11 or D12. None requires the deleted finite A17 endpoint.

The [B11 half-transfer](fresh_order12_affine_b_half_transfer_20261017.md)
covers all14 affected high contexts and28 induced maps. Radical cores
seal the finite alternative; nonsingular cores admit only a uniquely
forced zero-residual simple extension, or already have a sealed finite
B11 component. A retained ideal half plus a finite half supplies one.
The complete69 low-core variants do not use either changed pair. Every
immediate D13 occurrence of these B pairs is excluded by a strict Schur
gap of1/9 or1/11, including arbitrary spectators. No finite coupon is added.

The [C11 endpoint deletion](fresh_order12_affine_c_endpoint_reduction_20261017.md)
independently reconstructs all induced weighted embeddings of the758
eligible high patterns. Its15 successful contexts and30 maps avoid the
endpoint pair. The complete low domain and immediate D13 alternatives
are B or D, with no C occurrence. The report checks each later transfer's
actual terminal outputs, including the new cycle12/D11 interfaces and
A17 case662. No finite coupon is added.

These changes compose without a circular replacement. First discharge any
path of length at most8 by the complete short-path branch, using the new
universal cycle12 theorem if it is an antipodal pair. This need not leave
a unit at that particular corner: the cycle theorem pays the whole face,
returning finite coupons or ideal paths of length at most4. It therefore
does not recurse into another changed pair. New D alternatives use finite
payments or ideal paths of length at most5, retained in every family.
After this short-path branch, new B alternatives use retained finite
halves/units and the same retained B half. The C and A17 deletions preserve
those outputs. The older
transfers have only the finite, short-path, same-order D half or explicit
small B/D outputs checked in the local provenance reports. The latter
are covered by the universal D transfer and the complete B immediate
analysis. The complete2,614 high overlaps,66 chosen low closures with
three optional extensions, and disjoint-core exclusions remain the
classification inputs; no sampling replaces them.

Triangles retain both required payments. The finite elementary rule uses
distance at most8, whereas the only deleted finite pair has distance16.
For ideal factors the elementary rule uses distance at most5, or distance
at most8 with a nearby feature; exceptional whole factors have order at
most9. The cycle's deleted antipodal distance6 pairs have no feature.
The B/C/D changed pairs have distances10/11/9. The source reconstructs
the distances and feature sets and checks every elementary triangle pair
at weight one. Thus triangles receive at least2 and quadrilaterals at
least1 under this single rule, including at nonsimple cusps.

## The exact finite budget includes every new coupon

The [revised forest theorem](fresh_a17_endpoint_d11_upgrade_forest_budget_20261017.md)
and its [independent checker](fresh_a17_endpoint_d11_upgrade_forest_budget_20261017_review.md)
prove the exact finite maximum544 at rank69 after the A17 deletion and
D11 upgrade. They keep A17 and D17 separate: the former has base135 and
coupon roles, while the latter has base136 without those roles. The
complete finite-type domination and1,608 large-component vectors reduce
the remaining problem to exact Bellman bounds with attaining witnesses.
The independent implementation counts pairs directly and verifies every
upper-bound and attainment certificate.

The two cycle12 coupons also fit544, by the
[A11 domination theorem](fresh_a11_cycle_coupon_domination_20261017.md).
Replace each whole A11 component by D11. Both have same-component cost55.
The new A11 total cross-cost against A1 or A2 is1; upgraded D11 has totals
5/2 and3 respectively. Against every other component D11 costs at least
as much. Large/large interactions vanish, so successive substitutions
preserve rank and never decrease the total cost. The old maximizing
forests omit A11 and remain attaining witnesses. In particular A1+4D17
attains544. This establishes the budget for the actual simultaneous rule,
not a smaller rule omitting its new coupons.

## Ideal envelope and the changed concentration domain

All affine factors of order12 are classical A/B/C/D, so their maximum
cost is65. The older factors are unchanged. A cusp whose Euclidean link
has V vertices and factor ranks r_i therefore has total weight at most
V times the sum of c(r_i), with

|Affine rank r|c(r)|
|---|---:|
|1 through10|r/2|
|11|65/12|
|12|75/13|
|13|43/7|
|14|32/5|
|15|107/16|
|16|120/17|
|17 and larger|163/23|

The factor multiplier is V/(r+1), from the inherited
[actual factor-fiber count](fresh_cusp_factor_payment_framework_20261017.md);
it does not require a simple cusp. For total rank d=n-1 set

    H(t)=product_i(1+t+...+t^r_i), V=H(1),
    G=sum_j [t^j]H(t)*min(j,d-j).

The normalized ideal curvature is(5G-dV)/(2V). Subtracting the sum of
c(r_i) defines the ideal margin.

The inherited [all-rank convex-order lemma](fresh_distance14_unrestricted_dimension156_20261017.md)
decreases this curvature when two ranks are concentrated. The old low-rank
cap11 cannot be reused: rank11 now
has a cheaper cost. Instead concentrate pairs of ranks below10 only until
one disappears or reaches10. Their total linear cost is preserved. Keep
arbitrary rank11 parts. Two rank12 parts can be replaced by11 and13,
increasing cost by

    65/12 + 43/7 - 2*(75/13) = 23/1092 > 0.

Two ranks at least18 concentrate until one reaches17, preserving the
constant high-rank cost. Every nontrivial move increases the bounded
integer sum of squared ranks. Hence it suffices to retain arbitrary
copies of10,11,13,14,15,16,17; at most one12; at most one rank1 through9;
and at most one rank at least18. Curvature decreases and cost does not,
so the terminal margin is no larger than the original margin.

At d=68 this complete terminal domain has2,060 profiles. Direct exact
polynomial multiplication gives positive margin for every profile. The
minimum is3641/13608, attained at(11,11,11,11,11,13). No old terminal
profile list is reused.

## The strict contradiction in dimension69

Let D5 be the sum of5 minus the number of sides over all two-faces, and
let v_f count finite vertices. Face coverage and the finite/ideal budgets
give

    D5 <= 544*v_f + sum_(ideal v) V_v*sum_i c(r_vi).

The inherited [shape inequality](fresh_cusp_link_shape_curvature_20261017.md)
gives at n=69, d=68

    D5 >= 544*v_f + sum_(ideal v) (5G_v-dV_v)/2 + 340.

Here beta69=(69-1)(69-5)/8=544. Every ideal difference is strictly positive,
the finite difference is zero, and the additional constant340 is positive
even when there are no ideal vertices. The two inequalities contradict
each other. This excludes dimension69; the completed ceiling69 supplies
the higher-dimensional tail, proving ceiling68.

The present finite rule does not prove a further reduction: at ambient68
its maximum is still544, attained by4D17, while beta68=527. The resulting
finite deficit17 remains even though all1,906 terminal ideal profiles
there were positive in the [conditional comparison](fresh_rank11_cost65_scalar_frontier_20261017.md).
No ceiling67 is asserted.

Compact arithmetic remains at most21 and arbitrary compact Coxeter
polytopes at most29. The effective compact-facet ambient ceiling is68;
the older auxiliary proof remains79. There is no new finite-volume
example or new construction dimension.

## Verification status

Primary generation and a fresh separate-process replay both terminated
exit0/PASS on the complete2,060 profile domain and all direct local
interfaces. The [independent standard-library reconstruction](fresh_ideal5_weighted_dimension68_20261017_review.md)
agrees with every profile polynomial, vertex count, curvature, cost and
margin using a separately enumerated optional-first domain, sliding
convolution, cumulative-composition and absolute-deviation identities.
It independently reconstructs all264 canonical order12 pair weights and
elementary triangle sets,228 A11 partner comparisons and76 finite attaining
witness costs. Its own report records the final bound generation and replay.

Separate full internal source/mathematical reads passed the corrected
concentration domain, actual cusp multiplicities, simultaneous finite
budget, zero finite slack, global additive strictness and inherited tail.
A further composition review checked the short8-first ordering and every
new local interface; small independent saved-set checks found elementary
triangle-set sizes60/61/60/62 for A/B/C/D, all retained at unit weight.
No lost payment or circular replacement was identified.

These global checks do not repeat the inherited all-label geometric
censuses, supporting-root closure proofs or complete Bellman upper-bound
calculation. Those frozen proofs and their separate checkers retain their
stated scopes. Source, certificate and report are final/frozen. All reviews
are internal; no external refereeing or proof-assistant formalization is
claimed.
