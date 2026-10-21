# Independent review of the rank-22 three-hub weights

**Generation and a fresh saved replay pass.** The independent exact maximum is

    21977/500 = 43.954 < 44,

with margin **23/500**. The [reviewer](../scripts/fresh_sqrt5_dimension22_h3_cross_20261016_review.py) imports neither an optimizer nor a producer. The [certificate](fresh_sqrt5_dimension22_h3_cross_20261016_review.json) binds the [primary weight data](fresh_sqrt5_dimension22_h3_cross_20261016_cap2.json), the separate graph engine, and the literal face inputs. This verifies a weight rule; a global dimension conclusion requires the separate complete geometric coverage and field composition.

The exact checks cover:

| Obligation | Count |
|---|---:|
| Original literal face maps rechecked | 2,397 |
| Original literals fitting dimension 22 | 1,746 |
| All finite-queue singular stops accounted for before filtering | 2,070 |
| Finite stops fitting dimension 22 | 1,113 |
| Prior catalogue rows fitting dimension 22 | 1,038 |
| Total checked face rows | 2,151 |
| Common external-count cases | 89,255 |
| Exact vertex count contexts | 410 |
| Finite irreducible component types through rank 22 | 71 |
| Nonnegative bonus variables | 30,169 |

Every checked quadrilateral sum is at least one. The smallest slack is exactly **3/560000**. The four orbit families have 232 connected pair orbits, 178 A1 destination-root orbits, 97 I2(5) destination-root orbits, and six H3 destination-root orbits. The last six occur in A3, A4, D4 and D5 components. All orbit types, labels, marked roots/pairs and complete variable domains are reconstructed with independent canonical keys; the producer's encoded keys serve only as opaque variable names.

## The three cross families are disjoint

Let a,b,c be the **actual** numbers of entire A1, I2(5), and H3 components in a finite vertex diagram. The coefficients select a context

    (a, min(b,2), min(c,1)).

The cross payments themselves use the actual denominators:

| Source component | Excluded destination component types | Denominator |
|---|---|---:|
| A1 | A1 | a |
| I2(5) | A1, I2(5) | 2b |
| H3 | A1, I2(5), H3 | 3c |

Thus A1-to-H3 pairs remain in the A1 family, and I2(5)-to-H3 pairs remain in the I2(5) family. Neither is counted again by the new H3 family. Two roots in the same connected component use the connected rule instead. The reviewer asserts directly at each marked corner that at most one of the four families applies.

For any eligible destination root there are exactly a, 2b, or 3c source roots. Summing over source roots cancels the corresponding denominator. Capping b and recording only H3 presence affects the selected coefficient, never the denominator. This is why component costs remain additive once exact a,b,c have been fixed.

## All forests and all common spectators

The independent graph engine constructs finite A/B/D/E/F/H/I2 diagrams and identifies all rooted and marked-pair orbits. It includes the rank-two aliases A2=I2(3) and B2=I2(4) exactly once. The 71 component models are the complete finite alphabet inherited from the reviewed Q(sqrt5) classification, with ranks at most 22.

For each exact vertex triple a,b,c, the remaining rank is

    22-a-2b-3c.

The residual forest excludes exactly A1, I2(5), and H3 because their numbers are already fixed. Its possible total ranks are zero and every integer at least two: A2 and A3 suffice. A type-first unbounded-coin recurrence, with ascending residual rank and reversed component order, independently computes the maximum over this entire finite-forest domain. This differs from the producer's rank-first recurrence. All 410 context maxima agree exactly, and each producer maximizing forest is checked directly for its rank, all three hub counts, and its rational total. H3 components themselves retain their internal weight and can receive A1/I2(5) rooted contributions; removing them from the residual DP does not discard those contributions.

For a k-node zero-touch face core, the external rank is 24-k. The reviewer independently enumerates every common triple of external A1, I2(5), H3 counts whose remaining rank is zero or at least two. The **same** external triple is added to the internal counts at all four corners. It reconstructs the four literal finite corner diagrams and evaluates their payments exactly. The dimension filter is checked on both sides of every old literal map, including the exact dotted parameter.

All 2,070 finite-queue singular stops are accounted for before filtering, including expanded stops. The old catalogue plus those fitting stops supplies the 2,151 rows here. The complete geometric coverage binder must still account for its other branch reductions and verify that every actual face reaches an input row or an unchanged baseline payment. This review does not replace that obligation.

The independent triangular-face certificate supplies all 76 markings with total baseline at least two. Every new coefficient is nonnegative, so the triangular proof and the short-path/entire-H4 payments used by the attachment queues remain valid. No attachment is pruned merely because an H3 bonus occurs in an intermediate principal diagram.

## Validation and interface

The saved primary explicitly declares dimension 22, target 44, capped I2(5) context 2, H3-presence context 1, and the feasible rational rule. The reviewer asserts each of these metadata fields and the full variable index set. It checks the exact source bindings of its input certificates. Two separate source reads found no normalization, domain or dynamic-programming inconsistency.

The module exposes `Rule(primary).payment(row)` with records containing `external_A1`, `external_I25`, `external_H3`, `corner_weights`, and `payment`. It also exposes `component_bounds(rule)`, returning models, totals keyed by exact `(a,b,c)`, sorted context maxima, and orbit multiplicities. These are available to the separate coverage binder.

This is an exact feasible-rule check, not a proof of numerical optimality. It makes no standalone global dimension, new-example, external-validation or priority claim. All frozen earlier reviews remain unchanged.

```sh
python3 scripts/fresh_sqrt5_dimension22_h3_cross_20261016_review.py
python3 scripts/fresh_sqrt5_dimension22_h3_cross_20261016_review.py --verify
```
