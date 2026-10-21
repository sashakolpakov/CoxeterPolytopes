# The sparse label dependence of the T8 rank obstruction

A small literal witness check removes the unused leaf-label restrictions from the two frozen T8 polynomial theorems. In particular, it excludes the [last saved mandatory state](fresh_h9_t8_middle_20261014_search_case_08.json), although that state lies outside the narrower written hypotheses of the earlier reports. No further search or polynomial bank was run, and no later surviving assignment is asserted.

The [checker](../scripts/fresh_h9_t8_sparse_frontier_20261014.py) and [certificate](fresh_h9_t8_sparse_frontier_20261014.json) pass generation and saved replay. They use a direct two-coefficient Fraction calculation for the new null vectors and bind the frozen incidence and polynomial certificates. They do not repeat the earlier coefficient banks.

The only old guards involving roots outside

```text
S + {11,14,17},  S=[1,2,4,5,7,8,12,13],
```

were the full-vertex witnesses bounding the contacts of 14 with the two anchors. They can be replaced by these actual four-wall faces:

| Anchor | Four-wall support | Containing vertex index |
|---|---|---:|
| 11 | `[11,14,2,4]` | 72 |
| 17 | `[17,14,2,4]` | 322 |

For either support in the displayed order, test label four on the anchor–14 edge, retain labels three on `14–2` and `14–4`, and lower all three outer contacts to right angles. The resulting star has positive null vector

```text
(sqrt(2)/2, 1, 1/2, 1/2).
```

The exact Gram times this vector is zero. A higher tested label or any nonright outer contact only decreases its quadratic value. Since the four walls belong to an actual compact vertex, their Gram must be positive definite. Therefore both anchor–14 labels are at most three, independently of every other leaf and every anchor-own label.

All remaining cross-block, anchor-budget, third-wall-budget and shared-root guards already use only these eleven selected roots. The checker confirms their literal support containment and verifies the two missing triples that force nonempty anchor-own supports. Their proofs therefore impose no labels on leaves 0, 3, 6, 15 or 16 beyond whatever the incidence and compact Coxeter conditions themselves require.

The exact remaining label hypotheses are:

```text
(1,2) >= 7,   (4,5) = 4,   (7,8) >= 4,   (12,13) >= 4;
(1,14) = 2,  (2,14) = (4,14) = (5,14) = 3.
```

All inequalities here refer to finite integer labels, with no upper cutoff. The four own contacts of anchors 11 and 17 with roots 12 and 13 are arbitrary finite labels at least two; each anchor has at least one nonright own contact because its three-wall set is a missing face. All other mandatory and omitted meeting contacts may vary.

Under these sparse hypotheses, the domains, residual formulas and exact positive coefficient banks of the [complete own-pattern theorem](fresh_h9_t8_asymmetric_20261014.md) and its [middle-label extension](fresh_h9_t8_middle_20261014.md) apply without change. Their rank contradiction therefore excludes this larger family. This is a logical dependency reduction, not a new arithmetic parameter sweep.

The saved state with `(6,7)=(7,15)=2` satisfies every sparse hypothesis, so those decreased leaf contacts do not provide an escape from the obstruction. A future continuation should skip this entire sparse family, then seek an assignment changing one of the eight listed core contacts. This note does not run that continuation and does not exclude the whole case-8 incidence.

```sh
python scripts/fresh_h9_t8_sparse_frontier_20261014.py --verify
```
