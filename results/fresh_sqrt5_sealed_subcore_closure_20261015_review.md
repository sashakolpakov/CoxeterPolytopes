# Independent sealed-subcore containment review

**PASS.** The [separate checker](../scripts/fresh_sqrt5_sealed_subcore_closure_20261015_review.py) and [result](fresh_sqrt5_sealed_subcore_closure_20261015_review.json) independently reproduce all 458 first extensions and all 15,922 eligible white subsets in the [primary certificate](fresh_sqrt5_sealed_subcore_closure_20261015.md). Generation and saved replay pass.

The checker uses explicit complete-graph isomorphism backtracking, not the producer's canonicalization. It tries all eight permutations preserving the two black opposite pairs and every compatible white bijection. Coarse degree signatures only discard provably nonisomorphic candidates. Every reported map is checked against all literal edge labels and the black/white partition. The complete match sets agree, including unsuccessful subsets.

Independently recomputed full physical Gram inertias give the same split:

| First extension | Count |
|---|---:|
| Contains a reviewed sealed subcore | 86 |
| Singular, with no such subcore | 61 |
| Nonsingular, with no such subcore | 311 |

All 86 covered extensions are singular. The containment implication is sound: a subcore retains all four black facets and their opposite pairing; every omitted white node is an external white node for that same marked quadrilateral. The independently reviewed forest-sensitive rule already covers arbitrary such extensions of the sealed subcore, either through a baseline payment or through orthogonal external white components. It therefore pays the full face, not merely its displayed principal diagram.

The arithmetic uses the frozen independent Fraction field helper from the earlier contact review, read-only. The accepted first-extension input and the sealed-core arbitrary-extension theorem remain explicit reviewed dependencies. This result does not assert closure of the remaining 61/311 extensions, arbitrary-depth attachment trees, or the 24 dotted families, and gives no dimension bound.
