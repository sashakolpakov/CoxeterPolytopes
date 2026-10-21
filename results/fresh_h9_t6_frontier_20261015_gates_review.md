# Independent literal and source review of the T6 binder

**PASS.** The binder contains 155 inherited conditions and three valid literal transfers of the reviewed three-large-block theorem, for 158 distinct predicates. The separate [checker](../scripts/fresh_h9_t6_frontier_20261015_gates_review.py) and [saved output](fresh_h9_t6_frontier_20261015_gates_review.json) reconstruct the premises without importing the producer or its arithmetic. This is a source/literal audit, not another exhaustive search.

The checker rebuilds all 488 vertices from the complete minimal-nonface clutter, including its 11 missing triangles and 10 disjoint pairs. It checks every role assignment and threshold in the 23 dense-high, 88 strong-contact and 44 two-dense conditions, including all 2,288 recorded finite-face and threshold-pair witnesses. The analytic inequalities behind these three template families remain inherited from their hash-bound independent reviews.

For the three-large-block theorem, the four doubled blocks are `(1,2)`, `(6,7)`, `(9,10)`, `(12,13)`. Independent enumeration of the missing-triangle leaves and actual path faces gives:

| Blocks | All eligible ordered leaf choices | Selected choice |
|---|---|---|
| 012 | (14,5,8), (14,5,17) | (14,5,8) |
| 013 | (14,5,11), (14,5,17) | (14,5,11) |
| 023 | none | none |
| 123 | (15,16,11) | (15,16,11) |

Each selected triple has the three literal missing Lannér triples, an actual spherical leaf triangle, twelve actual cross-block path guards and twelve actual foreign-leaf path guards. A block label at least six forces each tested contact right: otherwise the path has a nonpositive three-root principal minor. The separately reviewed three-timelike-leaf argument therefore transfers to precisely these three placements. No whole-nerve symmetry or relabeling is assumed. This transfer by itself does **not** exclude the 023 placement.

The empty-start independent traversal's arithmetic and recursive function bodies are unchanged from the [audited star traversal](fresh_h9_star_closure_hostile_20261015.md). The changed run inputs are case 6, its incidence, the 158-predicate binder, and corresponding counters/status strings. The recursion still starts with no labels and has no node or time cutoff. This review separately verifies all sixteen simultaneous twin swaps on the full T6 vertex set and clutter; every chosen missing-triangle leaf lies outside all twin pairs, so the four ordering guards can be imposed simultaneously. The symbolic category `7` still represents every finite label at least seven. No fixed-equality predicate mistakes that category for the literal label seven.

The producer's exhaustive-search conclusion and counters belong to its own saved traversal and independent replay. This review does not rerun either search or claim a further catalogue-wide exclusion. It validates the binder and source adaptation supporting that calculation.

Generation and saved replay of the small review, syntax, manifest and local-link checks pass. Frozen inputs are bound by SHA-256; no prior artifact was changed.

```sh
python scripts/fresh_h9_t6_frontier_20261015_gates_review.py
python scripts/fresh_h9_t6_frontier_20261015_gates_review.py --verify
```
