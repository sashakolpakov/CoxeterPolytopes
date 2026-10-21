# Twelve literal core placements in star case 0

The two core-only single-contact obstructions transfer to **12 literal placements**, giving 48 label predicates. The [transfer checker](../scripts/fresh_h9_star_core_transfers_20261014.py) and [certificate](fresh_h9_star_core_transfers_20261014.json) pass generation and saved replay. They verify the exact incidence and hypotheses of each transfer; they do not repeat the determinant or continuous-parameter proofs.

Write the four block pairs as
\[
B_0=(3,4),\ B_1=(6,7),\ B_2=(9,10),\ B_3=(12,13),
\]
with positive leaves \(P=(2,5,8,11)\) and cut leaves \(C=(14,15,16,17)\). For any anchor block \(i\), choose another owner block \(j\). Take the anchors \(P_i,C_i\), and the third leaf \(C_j\) when \(j<i\), or \(P_j\) when \(j>i\).

The source core consists of all eight block roots, anchors 11 and 17, and third leaf 14. The transfer maps its first block to \(B_j\), its last block to \(B_i\), the two remaining blocks to the other two block pairs, and \((11,17,14)\) to the selected anchors and third leaf. Each map is injective. The checker proves equality of the complete maximal induced face lists and the complete induced minimal-nonface lists, not merely containment of a few convenient witnesses. There are nine maximal faces and three minimal nonfaces in every core. All 108 mapped maximal faces have explicit containing actual vertices.

The source theorems are the [shared-single core obstruction](fresh_h9_star_mixed_singles_20261014.md) and [opposite-single core obstruction](fresh_h9_star_opposite_core_20261014.md). In every placement the anchor-block internal label is at least 7; the other three internal labels are at least 4; the third leaf has its two own labels equal to 3; and the anchors have either shared or opposite single own labels equal to 3, with their other own contacts right. Both root orientations are included. All unlisted contacts and dotted distances remain free.

All domain witnesses and rank calculations in those theorems use only their stated rank core. Preserving the complete induced core face list therefore transports the actual-vertex constraints, and preserving the block/contact hypotheses transports the algebra unchanged. No global automorphism of the full incidence is assumed or needed. The two middle-block labels are interchangeable in the source theorems, so choosing their increasing order loses no hypothesis.

The certificate records four predicates for each placement. It is conditional on the bound source theorems; a search composition additionally requires independent reviews of both source theorems and these transfers. This is not a whole-incidence exclusion.

```sh
python scripts/fresh_h9_star_core_transfers_20261014.py --verify
```

The [independent review](fresh_h9_star_core_transfers_20261014_review.md) passes generation and saved replay. Its separate standard-library implementation checks all 26,624 induced-subset decisions for the source and 12 targets, all 48 predicates, 108 maximal-face witnesses and 4,123 source-proof witness occurrences. It also binds both source-theorem independent reviews. No source algebra is re-enumerated by this incidence review.

The prior 1,271-visit cursor is preserved as the input to the separate continuation below.

The [composition checker](../scripts/fresh_h9_star_core_transfers_20261014_gates.py) requires all three independent reviews before appending predicates. Four of the 48 images are already among the previous 209 predicates; the other 44 give [253 proved predicates](fresh_h9_star_core_transfers_20261014_gates.json). Generation and saved replay pass. The [continuation wrapper](../scripts/fresh_h9_star_core_transfers_20261014_search.py) retains the frozen DFS implementation, changing only the namespace and bound certificate references.

## Bounded continuation and preserved frontier

[Batch 011](fresh_h9_star_core_transfers_20261014_search_batch_011.json) resumes the saved 1,271-visit stack with a budget of 100,000 additional visits. It stops naturally at the first full mandatory-label/bridge survivor after **19,566 new visits**, giving **20,837 cumulative visits and nine pending frames**. The local-transfer predicates prune 3,150 partial states. The remaining 80,434 visits of the configured budget are unspent; this is a survivor stop, not exhaustion or a node-limit failure.

The new state's own-contact data are:

| Block roots | Internal label | Positive leaf / own labels | Cut leaf / own labels |
|---|---:|---|---|
| 3, 4 | 4 | 2 / (3,3) | 14 / (3,4) |
| 6, 7 | ≥7 | 5 / (2,3) | 15 / (2,3) |
| 9, 10 | ≥7 | 8 / (2,3) | 16 / (2,3) |
| 12, 13 | ≥7 | 11 / (2,3) | 17 / (2,3) |

The third leaf's \((3,4)\) profile lies outside the new core theorems, which require \((3,3)\). This is only a necessary assignment: omitted meeting contacts and dotted distances still need a consistent Lorentzian metric. No compact or finite-volume realization is claimed.

The engine saves all pending paths, chosen domain indices and sibling cursors every 1,000 visits. Its terminal file retains the complete stack and solution. A continuation must validate those paths and supply a proved gate before discarding this survivor; it cannot silently restart or skip the saved state. No prior proof, checkpoint or search prefix was edited. The final saved replay passed, reproducing the complete terminal JSON and trace; its measured wall time was **89.812 seconds**, including startup and provenance checks.

```sh
python scripts/fresh_h9_star_core_transfers_20261014_gates.py --verify
python scripts/fresh_h9_star_core_transfers_20261014_search.py --batch 11 --verify
```

## Next structural task

The same 12 induced-core maps can potentially transfer the [finite-core theorem](fresh_h9_star_single_right_20261014_finite_core.md), whose four internal labels are each 4 or 5 and whose three core leaves have own contacts \((3,3)\). There are 16 source predicates, so this suggests at most 192 mapped predicates before deduplication. That extension has **not** been implemented or used here. It requires a separate literal-hypothesis check and independently reviewed binder; it would reuse the already reviewed determinant bank rather than create another arithmetic bank.

That finite-core extension does **not** cover the newly saved state, which has three high internal labels. The immediate mathematical question is instead whether the core argument extends to the third leaf's \((3,4)\) profile, with fresh actual-vertex angle bounds and a consistent treatment of all three unbounded internal labels. No such extension or further search is started in this phase.
