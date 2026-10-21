# Actual payment reallocation reveals a second retained-rule obstruction

This model allows the three previously obstructing source/context bonus coefficients to decrease and adds payments at the other corners of their three constraining quadrilaterals. Its exact optimum is **2579221/60000 ≈ 42.98701667**, above the marked-facet target 35.

The earlier forest 7A1+2I2(5)+2A5 is released. A different forest, **6A1+H3+4A3**, with an isolated A1 marking the facet, gives the new exact lower bound. Neither long-component discount affects it. Its remaining fixed old bonuses cost 2579221/60000, and all compensation coefficients on this forest are nonnegative. The exact LP dual consists of this single vertex constraint with multiplier 1. Thus the result exposes the cost of leaving most old dimension-22 bonuses fixed; it does not show that payment reallocation or the new same-hub family cannot help in a broader model.

The changes allowed by this model are concrete:

- The A5 I2(5)-cross, I2(5) A1-cross, and A5 connected bonus families can decrease to zero in the actual context (A1,I2(5),H3)=(7,2,0).
- Eight alternative corner orbits from prior quadrilateral rows 10,80,94 receive independent nonnegative payments in their literal contexts. These include connected B6/E7 and cross corners A1×A4, I2(5)×B4/E6/H3/A4, and H3×A3.
- A new family pays between roots in two separate I2(5) components in context (6,3,0). Its payment is γ/[2b(b−1)], where b is the number of I2(5) components. Its whole-vertex total is γ; if the facet root lies in an I2(5) component, its total is γ(b−1)/b.

The last normalization was separately checked by the construction agent: there are 4*C(b,2)=2b(b−1) unordered cross-root pairs, and a marked root removes 2(b−1) of them. That narrow source read was not a full DP or face-coverage review. The producer also recounts individual root pairs for all 2,272 exact context maxima, independently of the component attribution used by the dynamic program.

All 1,979 inherited baseline-zero quadrilateral rows remain explicit constraints in 71,036 common spectator contexts, giving 4,278 distinct affine profiles. The selected 163 triangle families give 122,416 marked-context cases and 1,712 affine profiles. The exact rational primal satisfies every retained face obligation and the complete marked rank 21 finite-forest DP. The new alternative-corner and same-hub variables are zero at this optimum because the second fixed forest already blocks further improvement.

The [source](../scripts/fresh_sqrt5_marked_reallocation_proposal_20261016.py) generates and freshly replays the [certificate](fresh_sqrt5_marked_reallocation_proposal_20261016.json), which records the rational primal, exact dual, every vertex context maximum, and all retained affine face profiles. The [preceding context-scaling diagnostic](fresh_sqrt5_marked_old_bonus_context_scaling_20261016.md) identifies the three quadrilaterals used to choose the reallocation directions.

The next model should free the existing old orbit weights themselves, rather than retain almost all coefficients from an dimension-22 rule optimized for a different budget.

**Scope:** only the inherited zero-corner quadrilaterals and the selected triangle families are covered. Arbitrary connected-white extensions and newly affected paid-corner quadrilaterals still require complete coverage. This is a restricted-model diagnostic; no compact arithmetic example or dimension bound is obtained.
