# Complete finite8 / ideal5 face payments with the existing finite coupons

This package completes the finite8 quadrilateral coverage that was explicitly left open in the77 proof. It proves a local payment theorem and its vertex budgets; it does **not** lower the current unrestricted dimension bound77. The exact unchanged-envelope obstruction is recorded below.

Let P be a finite-sided finite-volume hyperbolic Coxeter polytope of dimension at most80. Finite vertices are simple; ideal vertices may have arbitrarily many simplex factors. At a finite corner give a unit when its two marks have graph distance at most8, lie in one whole connected component of order at most17, or have distance9 with at least one mark a leaf of that whole component. Add the frozen half-coupons between a whole small component of rank at most5 and a whole D11–D16 component at a fork center/short leaf, a whole A13/A15 component at its unique middle node, or a whole A17 component at positions5/11. These are exactly the [finite8 budgeted rules](fresh_finite8_a17_d11_coupon_budget_20261017.md).

At an ideal corner, use the actual varying affine cusp factor. Pay all pairs in factors of order at most13. For orders14 through17 use the exact frozen ideal5 pair sets, as retained by the [new complete map census](fresh_distance8_ideal5_small_factor_20261017.md). At larger orders pay distance at most5, or distance6 through9 when an endpoint is within5 of a feature (a branch vertex or endpoint of a nonsimple bond).

**Every triangle receives at least2 and every quadrilateral receives at least1.** The finite-vertex sum is at most8n. The normalized ideal-factor envelope, for affine rank r, is unchanged:

|r|factor cost per link-vertex incidence|
|---|---:|
|1–12|r/2|
|13|87/14|
|14|20/3|
|15|111/16|
|16|127/17|
|at least17|163/23|

The finite and ideal bounds use their established actual-corner multiplicities. An ideal vertex with link factors of ranks r_i and V link vertices contributes the sum of its factor pair costs multiplied by `V/(r_i+1)`. No port-transitivity statement at a nonsimple cusp is assumed.

## Complete face coverage

The [ideal5/finite8 triangle theorem](fresh_ideal5_finite8_triangle_payments_20261017.md) already gives two payments per triangle. The new checker verifies that every elementary ideal triangle pair in each of the16 small target diagrams remains paid. The whole-factor alternatives in that proof have order at most9 and are retained. Thus changing finite9 to finite8 creates no triangle gap.

For a quadrilateral, first handle an existing adjacent marked path of length at most8. If the corner is finite it is paid directly. If ideal and not already paid, the universal [finite8 cycle](fresh_finite8_cycle_replacement_20261017.md), [ideal5 local](fresh_finite8_ideal5_local_replacement_20261017.md) and [cycle14](fresh_cycle14_ideal5_transfer_20261017.md) replacements give another retained payment. Factors of order at most13 are wholly paid. For classical tree factors of orders14 through17 the exact retained pair sets contain every pair of distance at most8; the checker verifies this explicitly. The only missing small short pairs are cyclic and lie in the universal cycle domains. For larger tree factors, an unpaid path is outside feature radius5, exactly the hypothesis of the tree replacement. These replacements already use finite8, protected-leaf9 or the retained coupons; finite distance9 with two nonleaf marks is never imported here.

We may therefore assume the star has **no** marked path of length at most8. Choose minimal indefinite opposite-side cores. The actual edge diagrams make both marked deletions positive definite; the closed opposite-side argument applies even when the quadrilateral has ideal vertices. Minimal cores of order at most5, including the unbounded opposite-black product in order3, are fully handled by the frozen [distance8 low-core package](fresh_distance8_low_core_20261017.md). Its66 surviving low/high overlaps have full marked-component bound at most12, which is retained on both finite and ideal sides.

For two high cores, disjoint cores are excluded by the complete [all-label disjoint certificate](fresh_distance8_high_disjoint_20261017.md). Overlapping cores are exactly the new [82,197-consistent-identification domain](fresh_distance8_high_overlap_20261017.md): strict edge tests leave2,614 contexts, and every possible additional cross-core edge is either a short path or excluded by a nonnegative quadratic witness. All2,614 have a complete supporting-white certificate. The finite outcomes pay by whole component at most17, a protected marked leaf at distance9, or the new adaptive finite lemma. Other outcomes are basic ideal5 payments or full varying affine components of order at most17.

Every selected pattern with possible factor order14–17 is included in the new complete1,153-pattern induced-map census. Its images require no new ideal pair beyond the retained sets and their validated transfer domains. The omitted cycle pairs have universal replacements. The nine omitted D13/D14 fork pairs have universal transfers whose immediate whole-affine alternatives have orders12/13. These alternatives remain wholly paid.

The removed D pairs at orders16 and17 are different: their transfer is literal, not an arbitrary embedding argument. The new census matches all18 affected whole marked-core keys exactly to the frozen sealed cases, and the new full overlap certificate supplies the absence of extra cross-core edges. The local replacements then give finite whole-component units or finite half-coupons, or a length-two marked path. Their old global no-short9 reduction is **not** imported as coverage of the new domain. The source checks the saved finite outcomes and confirms that every old finite unit is still paid by distance at most8 or whole component at most17. It likewise rechecks all144 D13/D14 equality outcomes under the finite whole-component rule. This closes the change from the older finite9 statement without treating its weaker finite hypothesis as interchangeable.

## Vertex costs and exact remaining scalar obstruction

The finite8 coupon budget proves the sum at each finite vertex is at most8n for n at most80. There is no longer an extra n charge for general finite distance-nine pairs. The exact ideal pair sets and all-rank ideal5/R5 spine formulas are unchanged, so the envelope in the table is the same as in the77 theorem.

The fixed-envelope shape comparison still fails at dimension77 on the formal rank profile `(12,16,16,16,16)`. Its normalized curvature-minus-budget is

`-16664/1085773`.

The link has1,085,773 vertices, so its unnormalized cusp deficit is16,664. The improved finite term has positive slack68 at n77, but no universal finite-vertex multiplicity relation has been proved that compensates this cusp deficit. The global strict shape constant is380; it does not by itself cancel an unbounded number of such formal cusp contributions. This profile is a scalar obstruction, not a claimed geometric realization.

For that single profile, eliminating its negative margin would require lowering the order13 factor's total pair cost by `16664/83521`, or lowering each of its four order17 factor costs by `4166/63869`. One whole pair would exceed either threshold, but such a reduction must be geometrically proved for the corresponding worst factor family and checked against all profiles. The source records the exact reduced-profile minima at dimensions70 through77; all remain negative with the unchanged envelope. It therefore makes no new dimension claim.

## Verification scope

The [binder](../scripts/fresh_finite8_ideal5_combined_payments_20261017.py) verifies all frozen input digests, the complete domain and selected-pattern interfaces, exact small pair-set equality, triangle and short8 retention, immediate ideal alternatives, and finite-only literal replacement payments. It also recomputes the unchanged-envelope diagnostic profiles. It does not repeat the82,197 overlap construction or1,153-pattern embedding search; those are separate freshly generated and replayed certificates.

Final generation83039 and separate fresh saved replay85271 both terminated exit0/PASS. All51 direct source hashes, syntax and report links passed. Root completed a full source/mathematical composition read, including the preceding binder and triangle proof, and found no gap in the short8-first reduction, complete domain split,1,153 maps, literal finite8-safe replacements,144 equality checks or the distinction between the improved finite budget and unchanged ideal envelope. This review did not duplicate the reduced-profile arithmetic or the separate overlap/embedding censuses. Source, data and report are final/frozen. Published finite/affine and minimal-indefinite classifications remain explicit mathematical inputs. No external refereeing or proof-assistant formalization is claimed.
