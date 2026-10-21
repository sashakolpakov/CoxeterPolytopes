# Exact repair after the cross-coface obstruction

The 241 new cross-coface cap capacities invalidate the previous fixed witness. An exact monotone repair, with every inherited matrix and local/joint obligation retained, has ratio **10.2470841112**, below the target 12. This is a repaired witness, not an optimized value or a feasible weight rule.

The [script](../scripts/fresh_sqrt5_cross_coface_cap_repair_20261017.py) and [saved certificate](fresh_sqrt5_cross_coface_cap_repair_20261017.json) invoke no optimizer. Generation and a fresh saved replay pass. All inputs are bound by hashes, including the full model's recursive `prepared.sources` dependencies.

## Changed masses and complete inherited matrix

The starting point is the [five-port diagnostic](fresh_sqrt5_joint_five_port_diagnostic_20261017.md), with 2,219 positive face masses: 1,864 simplex rows and 355 fixed actual quadrilateral provenances. Its ratio was 12.3449447329. The [cross-coface cap capacities](fresh_sqrt5_cross_coface_cap_capacity_20261017.md) require 53 simplex masses to increase; all were already positive and have dimension 1 or 4. The quadrilateral masses and their common provenance assignments remain unchanged.

For each new capacity, the script sets its simplex mass to at least its positive quad load, taking the maximum when several capacities share that simplex. It then rebuilds all 82,806 inherited face obligations and all 6,558 vertex budgets, and recounts all **333,829 positive-cost columns and 94,481 signed columns** at the changed masses. It checks every new cap capacity separately.

Every inherited signed transfer has positive terms only on quad rows and a negative term on its containing simplex. Consequently the unchanged quad masses and increased simplex masses preserve those transfers. The same monotonicity preserves the previously checked complete coface-base packing cones and their labelled subset capacities. No new independent provenance selection is used.

## Rechecked complete local covers and joint charges

The old local cover certificates cannot simply be reused at changed masses. The script revalidates all **2,830 literal disjoint-packing terms** for all **593 forests and 3,485 block orbits**. It recomputes the new orbit marginals and appends one-block packings only for uncovered marginal deficits. There are **240 added terms**. The saved additions and final marginal coverage give a fresh certificate for every nonnegative orbit-local simplex-block charge, without reoptimizing any local cone.

The raises change 241 block marginals in 177 forests. The script recomputes the complete affected orbit histograms for the joint three-, four-, and five-port families, retaining the literal unchanged quad footprints. It checks the full frozen domains:

| Family | Total inequalities checked | Orbits with recomputed simplex increments |
|---|---:|---:|
| Strengthened three-port | 44,260 | 11,448 |
| Four-port | 90,619 | 22,389 |
| Fifteen five-port facets | 2,153,595 | 34,256 |

For unchanged forests the bound inputs remain identical. Weighted binomial sum identities check every affected forest's new intersection multiplicities. The five-port calculation uses exact intersection sizes, so nonmonotone simplex coefficients cause no loss of information. Finally, the script takes the maximum required vertex mass across the inherited matrix, updated covers, and every joint charge, and rechecks all matrix and cap columns.

## Result and limits

The final dual objective is `79262181156857/5775000000000`; its normalizer is approximately 1.3394105938, giving the exact rational ratio stored in the certificate, approximately **10.247084111225892**. There are 176 increased vertex masses.

A simpler independently bounded repair adds one unit of mass for every additional simplex-corner occurrence, producing ratio 10.2314338178. The refined construction uses actual marginal deficits and joint loads and slightly improves that bound. The updated old-matrix costs alone, while retaining the previous vertex masses, give ratio about 11.39331; the explicit added packing covers dominate the final repair. Neither calculation is an optimality statement.

This package establishes that the new constraints have been incorporated into an actual exact repaired witness, rather than merely invalidating an old one. A ratio below 12 does **not** prove primal feasibility, rule out another dual above 12, or establish a dimension bound. It justifies testing a genuinely enlarged model. The unrestricted dimension-163 proof and the compact arithmetic dimension-21 proof are unchanged.
