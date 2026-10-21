# An exact two-term barrier at 36

Every rule in this specified payment family has marked-vertex maximum **at least 36**, above the required 35. The proof uses one quadrilateral and one vertex; it does not depend on any H4 compensation choice.

Let w be the inherited connected D5 payment on the pair joining the long-arm endpoint to a short-arm leaf. Quadrilateral 64 in the inherited face catalogue has corners A4+A4, A3+D5, A3+D5, and D4+D4. The first and last corners receive no payment in this family; each middle corner receives w. With common external components A1+A12, the quadrilateral obligation is therefore **2w ≥ 1**.

At the rank-21 forest **A1+4D5**, mark the isolated A1 as the facet root. Each D5 has eight baseline-paid pairs and two pairs paid w. All pairs involving the marked A1 are omitted, and pairs between D5 components receive zero. Neither long-component discount applies to rank-five D5. The marked total is consequently **32+8w ≥ 36**.

The [certificate](fresh_sqrt5_marked_h4_remainder_free_bonus_20261016.json) records the literal quadrilateral labels, its four corner types, all ten D5 root pairs, the marked forest, and the two exact dual multipliers: 4 for the face obligation and 1 for the vertex inequality. The [source](../scripts/fresh_sqrt5_marked_h4_remainder_free_bonus_20261016.py) reconstructs these data and verifies exact weak duality. This barrier also applies to the preceding free-bonus and H4 compensation models, because it uses no newly split H4 coefficient. Their previously reported near-35 numbers were lower bounds from incomplete vertex relaxations, not optima.

The run added one targeted distinction to the [H4 facet-context model](fresh_sqrt5_marked_h4_facet_context_free_bonus_20261016.md): at actual H4 count four, compensation distinguishes whether the remaining rank-five forest contains D4 or D5. The only positive cases are D5 and A1+D4. A separate narrow source read by the construction agent checked this partition and the common external-context accounting; it was not an independent full LP review.

The model has 30,392 variables and 48,018 distinct face profiles. It covers 163 selected triangle families in 166,711 marked contexts and all 1,979 inherited baseline-zero quadrilateral rows in 71,036 contexts. Complete marked finite-forest separation is retained. Generation, final literal-witness refresh, and a fresh saved replay verify the exact certificate.

The next enlargement should let the A4+A4 and D4+D4 corners carry some of this quadrilateral's cost, using the literal endpoint/leaf cross orbits. Changing H4 compensation alone cannot remove the barrier.

**Scope:** an exact obstruction to the specified payment family, not a dimension theorem. Complete coverage of newly affected triangles and paid-corner quadrilaterals remains open. The argument does not exclude a compact arithmetic polytope in dimension 21.
