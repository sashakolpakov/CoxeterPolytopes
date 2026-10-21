# Independent review of local transfers to four other tree incidences

The [separate checker](../scripts/fresh_h9_other_trees_20261014_transfer_review.py) and [certificate](fresh_h9_other_trees_20261014_transfer_review.json) pass generation and saved replay. They verify every analytic-theorem condition set used for cases 17, 22, 25 and 26 in the [transfer certificate](fresh_h9_other_trees_20261014_transfers.json). No isomorphism between these incidences and the old T8 incidence is assumed.

The review imports no producer code. It reconstructs each complete rank-nine vertex set from its minimal-nonface clutter and compares it with the frozen primal geometry. It then checks injectivity of every local map, all mapped missing triples, and all 2,524 actual-face premises. These include the positive eight-wall supports, all 45 or 55 core meeting pairs, every cross-block and anchor Schur budget, the replacement affine-star guards, and the nine own-support pattern guards. A pattern with unrestricted finite anchor contact needs only the literal meeting pair; every pattern forcing a right contact retains its stronger three-wall guard.

There are seventeen valid two-high-block maps and twelve sparse all-block maps in the supplied certificate. The new checker verifies all twenty-nine maps, not the exhaustiveness of their discovery. Every deduplicated condition set points to a verified map and is exactly the image of the reviewed theorem's fixed labels and unrestricted finite-label lower bounds. Cases 25 and 26 use only the two-high theorem; cases 17 and 22 also use the sparse polynomial theorem.

The numerical identities and all-parameter polynomial/negative-plane arguments are inherited from the separately reviewed frozen T8 certificates. They depend only on the listed local incidences and label hypotheses; translating those premises preserves the proofs. This review checks that translation, not a new arithmetic derivation of the old theorems.

Whole-incidence exclusion still requires a completed mandatory-label search using these valid gates. A surviving assignment or a node-limit stop would not be an exclusion. The separate primary and independent traversal certificates record those outcomes.

```sh
python scripts/fresh_h9_other_trees_20261014_transfer_review.py --verify
```
