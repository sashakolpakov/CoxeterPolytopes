# A retained-bonus obstruction to the marked compensation proposal

The exact optimum of this restricted proposal is **1516843/35000 ≈ 43.33837143**, above the required marked-facet budget 35. This is a failure of the specified retained-bonus model, not an obstruction to compensated weight rules in general.

The model keeps the independently verified dimension 22 H3 rule, indexing its bonuses with one additional common external A1 while dividing cross payments by the actual dimension 21 hub multiplicities. It adds two discounts on components of rank at least 6: ordinary distance-two pairs with both bonds3, and adjacent bond3 pairs. Nonnegative compensation variables cover H4, H3, I2(5), and I2(10) pair orbits. H4 and H3 channels distinguish whether the retained facet root lies in the same component. H4 compensation is disabled on forests containing five H4 components.

The one-row exact lower bound comes from **7A1+2I2(5)+2A5**, with one A1 root marking the facet. Neither discount changes this forest. Its retained cost is 1516843/35000, and every compensation coefficient is nonnegative. A rational feasible primal attains that cost, so the optimum is certified without relying on a floating-point dual.

The complete finite-component dynamic program checks every marked rank 21 forest, using exact counts of A1, I2(5), H3, and H4. It produces 2,272 context/marked-hub maxima, all independently recounted by literal individual root pairs to verify the hub normalization and facet-root debit. The underlying catalogue contains 68 finite component types and 3,617,699 marked types up to the symmetries recorded in the earlier census.

The selected triangle proposal comprises 163 literal families: the 76 minimal cores, four connected-white extension families of the H core, and the I2(10) tail family. All three corner forests are checked to be finite. Common spectator contexts and retained facet-root positions produce 122,416 cases and 1,685 distinct affine profiles. These selected diagram obligations are conservative proposal constraints; no completeness claim is made for arbitrary connected white extensions.

All 1,979 inherited baseline-zero quadrilateral rows pass in 71,036 common contexts, with minimum payment 140001/140000. Every corner has zero baseline, so neither proposed discount changes those rows. The additional A1 used in the old-rule index is a definition of the new payment lookup; it is not asserted to be an actual wall. It corresponds to a valid n22 spectator context, and using the smaller actual A1 denominator cannot reduce a corner payment.

The dominant cost of the obstruction is the old cross payment into A5. For each A5 component the retained contributions are baseline 7, connected bonus 20467/10000, and I2(5)-cross bonus 48899/5000. Each I2(5) component contributes baseline 1 plus A1-cross bonus 32247/17500. Thus further progress requires changing some old bonuses as well as compensating the new discounts.

The exact generation and fresh saved replay are performed by the [source](../scripts/fresh_sqrt5_marked_compensation_proposal_20261016.py); the [certificate](fresh_sqrt5_marked_compensation_proposal_20261016.json) records the rational primal, the lower-bound forest, every context maximum, all selected triangle profiles, and the inherited quadrilateral checks. The [discount census](fresh_sqrt5_marked_discount_census_20261016.md) and [actual triangle counterexamples](fresh_sqrt5_marked_discount_triangles_20261016.md) motivate the proposal.

**Scope:** no complete coverage of newly affected triangles or paid-corner quadrilaterals has been established. No compact arithmetic polytope is constructed or excluded, and no dimension bound changes.
