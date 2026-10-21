# Independent review of the sparse T8 hypotheses

**PASS.** The [independent checker](../scripts/fresh_h9_t8_sparse_frontier_20261014_review.py) and [certificate](fresh_h9_t8_sparse_frontier_20261014_review.json) validate the [sparse dependency extension](fresh_h9_t8_sparse_frontier_20261014.md). Generation and saved replay pass. The conclusion applies to compact H9 realizations of the literal case-8 incidence under exactly these eight conditions:

```text
(1,2) >= 7,   (4,5) = 4,   (7,8) >= 4,   (12,13) >= 4;
(1,14) = 2,  (2,14) = (4,14) = (5,14) = 3.
```

The inequalities refer to finite integer labels, without an upper cutoff. All four own contacts of anchors 11 and 17 with roots 12 and 13 may vary. Their nonempty supports follow from the two literal missing triples. Other mandatory and meeting labels need not retain their earlier saved values.

The review checks more than containment of guard supports in eleven roots. It rebuilds every retained or relaxed label in the local argument from the eight hypotheses or a preceding proved domain. It imports no producer arithmetic. The previously independent polynomial certificates are bound and reused; their coefficient banks are not regenerated here.

## Complete dependency audit

Put `S=[1,2,4,5,7,8,12,13]` and add roots 11, 14 and 17. The checker verifies that S occurs in actual vertex 172 and that every one of the 55 pairs among these eleven roots is a meeting pair. Thus all selected labels are finite. The following steps account for every selected Gram entry.

| Guard | Exact labels or previously proved facts it uses | Consequence |
|---|---|---|
| Twelve cross guards involving `(1,2)` | Internal label at least 7; test another nonright contact at label 3 | Roots of different S-blocks are right |
| Twelve other cross guards | Both internal labels at least 4; test their cross contact at label 3 | Remaining S-block cross contacts are right |
| Four retained-own-root anchor budgets | Cross-right S, internal lower bounds, one own contact at least 3 | Anchors are right to 1 and 2, have at most one other foreign contact, and that contact is label 3 |
| The same four budgets | A foreign contact exists; choose either active own root | Every active own contact is label 3 in that case |
| Two third-wall budgets | Cross-right S and fixed `14–2=14–4/5=3` | Root 14 is right to 7, 8, 12 and 13 |
| Two new affine stars | Fixed `14–2=14–4=3`; test `14–anchor >=4` | Both third-anchor labels are 2 or 3 |
| Shared-root triangles | Both active own contacts at least 3; test mutual anchor label at least 3 | Anchors with a shared own root are mutually right |

All actual support containments are checked literally; no incidence automorphism is assumed. The alternative own-support patterns have their own containing vertices. No row of this audit uses a label involving leaves 0, 3, 6, 15 or 16.

For the first twelve cross guards, the checker uses the three-node path with relaxed edge labels `(6,3)` and positive null vector `(sqrt(3)/2,1,1/2)`. The actual first label is larger, so the quadratic value is negative. For the other twelve it independently reconstructs the affine `(4,3,4)` path and null vector `(1,sqrt(2),sqrt(2),1)`. Every unmentioned edge is relaxed to right. Increasing any actual cosine decreases the quadratic value on these positive vectors, so unknown contacts cannot invalidate either exclusion.

For the anchor budgets, an active own contact costs at least `1/4`. A contact into either middle block costs at least `1/2`; a contact into the first block costs more than 1. Two foreign contacts therefore violate positive definiteness. A foreign label at least 4 costs at least 1 by itself. An own label at least 4 together with a foreign contact costs at least `1/2+1/2=1`, also forbidden. Both choices of retained own root are checked for both anchors, so this conclusion applies to every active own contact, including double-support patterns.

Each third-wall budget retains the fixed contacts to root 2 and one of roots 4 or 5, costing `1/4+1/4`. Any additional nonright contact to the full `(7,8)` or `(12,13)` block costs at least `1/2`. Their sum reaches the forbidden Schur value 1. This argument uses no outside leaf and no own-anchor label.

## The replacement stars remove the outside-leaf dependencies

The literal supports are `[11,14,2,4]` in actual vertex 72 and `[17,14,2,4]` in actual vertex 322. In either displayed order, set the anchor–14 label to 4, keep the two fixed label-3 contacts from 14 to roots 2 and 4, and relax all outer contacts to right. The independent exact calculation verifies

`G (sqrt(2)/2, 1, 1/2, 1/2)^T = 0`.

Every vector entry is positive. Larger tested labels or any nonright outer contact only decrease its quadratic value. Hence positive definiteness at the actual vertex excludes every anchor–14 label at least 4. This proof is independent of the anchors' own supports and all other leaf labels. It replaces all eighteen former full-vertex upper-contact witnesses; those witnesses are not silently inherited under the weaker assumptions.

The two missing triples `{11,12,13}` and `{12,13,17}` are checked against the literal incidence. An empty own support would make the corresponding triple positive definite, contradicting its being a missing face. For a shared own root, the actual three-wall support has an affine all-label-3 relaxation, forcing the mutual anchor contact right. For disjoint own supports the mutual contact remains arbitrary finite, covered by the whole cosine interval `[0,1]` in the existing polynomial relaxation.

## Why the existing polynomial banks now apply

After these deductions, all 55 selected Gram entries have exactly the domains used by the [independently checked complete own-pattern theorem](fresh_h9_t8_asymmetric_20261014_review.md) and its [middle-label extension](fresh_h9_t8_middle_20261014_review.md). There are no selected entries left dependent on outside labels.

The six continuous boxes still cover all nine nonempty own-support patterns. At last-block labels 4, 5 or 6, a single own label 3 would make its missing triple positive semidefinite, contrary to compactness. Thus the low-label single-support boxes correctly start at own label 4 and have no foreign contacts. The high-label threshold follows from the same exact Chebyshev calculation; the small review independently regenerates that polynomial and its rational threshold evaluations.

The middle-label split is exhaustive: label 4 uses the frozen 408-case bank; label 5 uses 256 changed cases and 152 unchanged cases; labels at least 6 prohibit any foreign contact to that middle block by Schur cost at least 1 plus positive own cost. These deductions use only the sparse domains established above.

The prior independent certificates reconstructed the determinant identities and every coefficient in those banks. Their source manifests are checked recursively here. Their positive residual determinant therefore contradicts the rank at most two left after the positive rank-eight support in ambient dimension ten. The additional freedom outside the selected eleven roots cannot repair that rank contradiction.

The latest saved mandatory state satisfies all eight hypotheses and is excluded. This audit runs no continuation and asserts no later survivor. It does not exclude the whole case-8 incidence or all compact H9 polytopes.

```sh
python3 scripts/fresh_h9_t8_sparse_frontier_20261014_review.py --verify
```

The scope is a fresh exact local-witness and logical dependency audit plus reuse of frozen independently verified polynomial banks. It is not another polynomial enumeration or a new metric search.
