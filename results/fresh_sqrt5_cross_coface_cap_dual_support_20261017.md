# The new cross-coface dual: exact support and restrictions

The new full-model dual has ratio **12.6927766928**, with 3,447 positive face masses and 770 positive vertex masses. Removing all 923 unselected merged quad rows still leaves an exact witness of ratio **12.3573995146**, using only 336 explicitly mapped actual quad provenances. No simplex row has to be restored in that restriction.

The [script](../scripts/fresh_sqrt5_cross_coface_cap_dual_support_20261017.py) and [exact data](fresh_sqrt5_cross_coface_cap_dual_support_20261017.json) reconstruct and check the full frozen [cross-coface model](fresh_sqrt5_cross_coface_cap_full_model_20261017.md). Generation and a fresh saved replay pass; no optimizer is invoked.

## Actual support

| Positive row category | Count |
|---|---:|
| Simplex rows in the preceding complete local/joint feature domain | 1,334 |
| Simplex rows needing fresh complete local features | 854 |
| Actual quad rows previously checked in the joint scan | 242 |
| Other selected actual quad rows | 94 |
| Unselected merged quad rows | 923 |

The 241 cross-coface capacities have 73 positive loads; 46 are tight and 27 have positive slack. Every capacity is rechecked exactly at the new masses. The certificate records all loads, slacks, containing-simplex rows, and actual quad coefficients.

Fresh full simplex-corner reconstruction covers all **2,188 positive simplex rows, 4,169 block orbits, and 6,763 corner occurrences**. It checks the literal pair and one-port multiplicities against the inherited row coefficients. This is a feature reconstruction, not a packing-cover certificate at the new masses.

## Exact restrictions

Each restriction starts with a named subset of the new exact masses, then uses the full model's monotone signed-column repair. Every resulting rational dual is checked against all **94,722 signed columns**, all saved positive columns, and all 82,806 face obligations.

| Seed domain | Final positive faces | Restored simplex rows outside seed | Exact bound, approximately |
|---|---:|---:|---:|
| All simplex rows and selected actual quad rows | 2,524 | 0 | 12.3573995146 |
| Old complete simplex domain and selected actual quad rows | 1,861 | 191 | 12.1996287336 |
| Old complete joint domain | 1,576 | 0 | 12.1891681974 |
| Simplex rows only | 2,188 | 0 | 10.7842619302 |

In the first restriction the 2,188 simplex masses and 336 selected actual quad masses are unchanged; only the unselected merged quad masses are removed, and the positive-cost vertex masses are recomputed. The second row must not be described as an old-feature-only witness: its 191 forced simplex rows are explicitly exported. No merged quad row is assigned an independently convenient provenance.

## Next obligation and limits

The clean 2,524-face restriction isolates an obstruction supported on actual quad provenances. Its 854 new simplex rows need freshly reconstructed local covers, and 94 actual quad rows need new joint feature checks. Old fixed-mass packing mixtures have not been reused without checking marginal requirements. A successor can impose additional cross-coface geometry on this consistent provenance domain before undertaking those checks.

These are exact witnesses for the declared finite model and its support restrictions. They do not establish a dimension bound or assert that all remaining local/geometric charge families are satisfied. The simplex-only ratio below 12 is a particular restricted witness, not an optimality or feasibility result.
