# Complete-face composition for compact arithmetic Q(√5)

The [binder](../scripts/fresh_sqrt5_dimension24_coverage_20261015.py) and [complete certificate](fresh_sqrt5_dimension24_coverage_20261015.json) pass generation and saved replay. They bind the completed dimension-27 coverage certificate, the final independently reviewed arithmetic queues and the three separate exact weight reviews.

The certified conclusion is: every compact arithmetic hyperbolic Coxeter polytope whose ground field is exactly Q(√5) has dimension at most 24. The [preceding coverage binder](fresh_sqrt5_dimension27_coverage_20261015.md) excludes dimensions 28 and 29; this composition adds dimensions 25, 26 and 27. The previously established compact upper bound 29 covers all higher dimensions. Other fields, noncompact polytopes and arbitrary nonarithmetic polytopes are outside this statement.

## Rules on whole vertex diagrams

Use the three [complete-queue weight rules](fresh_sqrt5_complete_face_weights_20261015.md), with their [separate exact implementation review](fresh_sqrt5_complete_face_weights_20261015_review.md). Their common baseline puts one on pairs at graph distance one or two and the reviewed extra units on the 5-edge and endpoint pair of an entire H4. All additional weights are nonnegative.

Connected-pair bonuses depend on the entire connected spherical component, its marked pair and the exact number a of entire A1 components in the vertex. A disconnected pair joining an A1 to a root of another entire component receives u_a/a. No A1–A1 pair is charged. For each root of the other component there are exactly a such pairs, so its total contribution is u_a. This justifies additive component totals for every fixed a; the separate review verifies the multiplicities and this normalization literally.

The reviewer recognizes labelled finite diagrams intrinsically, including the A2/I2(3) and B2/I2(4) aliases. Its exact type-first rank dynamic program is independent of the optimizer and of the producer's rank-first recurrence. It checks all finite component types through the target dimension and all realizable exact A1 counts. The sole impossible count is n−1, since the remaining rank one cannot be a non-A1 component. No numerical optimization status or optimality assertion is used.

| Dimension | Exact maximum vertex weight | Strict lower threshold |
|---|---:|---:|
| 25 | 59 | 60 |
| 26 | 31,017/500 | 65 |
| 27 | 67,009/1,000 | 143/2 |

The current binder independently recalculates these maxima from the literal rational rules and the reviewed finite-component implementation. It also rechecks all 76 minimal triangle markings using the weaker H4-only baseline. Every triangle gets at least two, and this persists under allowed spherical extensions. Nonnegative bonuses preserve all existing short-path and entire-H4 quadrilateral payments.

## Complete coverage with the current weights

The geometric and arithmetic reductions are those of the preceding complete coverage proof. An unpaid compact quadrilateral yields two minimal arithmetic Lannér pieces containing its opposite black pairs. All shared-white and finite cross contacts occur in the reviewed 5,440-support/1,180-class reduction. After field, signature and persistent-payment filters, the remaining cases are the 112 sealed singular finite cores, 44 sealed nonsingular rank exclusions, 130 finite cores with 458 possible first attachments, and 24 dotted parameter families.

The new weights are not assumed to preserve any old long-pair payment. The binder reconstructs all 458 literal first attachments, every sealed-subcore and representative map, all 147 singular first attachments, all 51 marked singular frontier states and all seventeen closure states. The seven touching continuations are explicitly mapped to their closure representatives. Each row is paid again under each current dimension-specific rule whenever its corner rank fits that dimension. The same direct check is made for all 264 singular dotted-queue states and all 2,070 singular finite-queue states.

The finite queue has 5,172 complete marked states from all 311 nonsingular literal starts, with no pending state. Every singular stopping choice is included, even if it admits further attachments or was paid by an earlier rule. Nonsingular no-touch stars are rank-excluded; every admissible next touching root is followed unless a persistent baseline payment or a necessary arithmetic obstruction ends that branch. The separate full queue review is a required source dependency. The dotted queue's 438 states cover all 32 exact parameter choices in the 24 families and likewise retain every singular stopping choice.

A full quadrilateral star in dimension n has n+2 roots. Thus the queue's 31-node ceiling covers every target n≤29 and is a dimension bound, not a search cutoff. For a displayed k-node singular core, its corner rank is k−2. It contributes no case at dimensions below that rank. Otherwise the external rank is n−k+2: every possible exact number e of external A1 components is checked, excluding only a leftover non-A1 rank of one. No assumption of sealing is made for an intermediate singular stopping case: it is paid only when the remaining white roots are orthogonal spectators, while all touching alternatives continue through the queue.

Conditional sealing is used exactly as in the reviewed branches. If every nonzero next-white contact is excluded or creates a persistent payment, all uncharged outside roots are orthogonal to the core. Its marked components are then entire. Orthogonal spectators may be arbitrary spherical forests, and their exact A1 contribution is included. A sealed nonsingular star with only positive-definite orthogonal spectators cannot have the required one-dimensional full-star relation.

## Strict contradiction

Once triangles receive two and quadrilaterals one, every k-gonal two-face has total weight at least 5−k; nonnegativity suffices for k≥5. Summing over faces gives total vertex weight at least `5 f₂ − binom(n,2) f₀`.

For a simple compact polytope, the symmetric nonnegative h-vector expresses f₂/f₀ as a weighted average of `(binom(n−i,2)+binom(i,2))/2`. Its minimum produces the thresholds in the table. The inequality is strict because h₀=hₙ=1 and the extreme paired coefficient is larger than the minimum. Thus the average vertex weight must exceed 60, 65 or 143/2 in dimensions 25, 26 or 27, respectively. The independently verified upper bounds contradict each case.

This composition is a source-bound finite verification plus the stated all-size geometric reductions, not an enumeration of global polytopes. It makes no priority claim. Complete-mode output exposes the precise ground field, excluded dimensions, upper bound and completion flag so that a separate all-field theorem cannot mistake a partial rule for complete coverage.

```sh
python3 scripts/fresh_sqrt5_dimension24_coverage_20261015.py
python3 scripts/fresh_sqrt5_dimension24_coverage_20261015.py --verify
```

The final certificate directly checks 2661 zero-touch rows, with dimension-feasible payment counts 25: 18330, 26: 20372, 27: 22565. It verifies 471 literal maps and 70 recursive manifests containing 568 source-hash references. Generation, saved replay, syntax and local-link checks pass.
