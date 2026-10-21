# Closing the two-chain stars and three new weight obligations

The 112 principal stars from the [all-length two-chain calculation](fresh_sqrt5_two_chain_completion_20261015.md) reduce to 70 marked diagram classes. Their complete unpaid attachment analysis leaves no open contact queue. Of the 112 original records, 109 are covered by the existing weight rule and three yield new exact local quadrilateral constraints. This closes the stated attachment construction, not the full Q(sqrt5) quadrilateral classification or a global dimension bound.

The [source](../scripts/fresh_sqrt5_two_chain_star_closure_20261015.py) and [certificate](fresh_sqrt5_two_chain_star_closure_20261015.json) preserve the all-length producer unchanged. Every quotient member has a literal full-label bijection preserving the opposite marked pairs. Known containment uses the reviewed 169 sealed cores and five closed starting configurations from the [attachment closure](fresh_sqrt5_singular_attachment_closure_20261015.md).

| Initial outcome | Marked classes | Original records |
|---|---:|---:|
| Contains a known sealed core | 32 | 54 |
| Contains a known closed starting configuration | 7 | 12 |
| New sealed, covered by current weights | 26 | 41 |
| New sealed, unpaid | 3 | 3 |
| One further unpaid contact | 2 | 2 |
| Total | 70 | 112 |

The containment search checks 9,758 white subsets, stopping at a sufficient literal witness in matched cases and checking every relevant size in unmatched cases. Every sealed-subcore match additionally verifies that the remaining white nodes are orthogonal to that subcore. The seven closed-start matches instead inherit its complete attachment proof; they are not assumed orthogonal.

The 31 unmatched classes have 303 complete next-white contact profiles over the full finite alphabet `{2,3,4,5,6,10}`. The exact outcomes are 230 kernel-relation failures, 40 cyclic-field failures, 29 certified affine principal witnesses, two persistent H4 payments and two surviving contacts. The only surviving contacts occur at classes 64 and 65. Their two explicit children each have five candidate next-white profiles; all ten violate the exact old kernel relation. Both children are covered by the current weight rule in either external-A1 context. The zero-addition configurations at classes 64 and 65 are also covered.

Here “sealed” is conditional on the absence of a persistent short-path or H4 payment. Any external white node directly touching the core supplies one of the completely enumerated profiles. A path through other external nodes cannot avoid a first direct contact. Thus a sealed unpaid configuration has its marked corner components fixed under arbitrary unpaid extensions. Other external components are orthogonal spectators; their possible effect on the presence of an A1 component is checked in both contexts. If a tested addition is already paid by a short path or H4, that payment persists under every valid finite-corner extension. These facts justify exhausting the two one-step children rather than imposing an attachment-depth cutoff.

## Three genuine local unpaid quadrilaterals

Classes 50, 55 and 62 are sealed but receive weight zero under the frozen rule. Each has exactly one connected marked corner; its other three marked pairs lie in different connected components. The certificate reconstructs the residual four-normal Gram matrix in a Lorentz three-space using the separate frozen Fraction arithmetic helper. It tests all six pair intersections, obtaining exactly the four prescribed timelike corner rays, rejects both opposite pairs, and gives an exact strict timelike center in the same future component as all four rays.

Consequently these are actual compact quadrilateral faces of local arithmetic-compatible Coxeter chambers. They are not merely singular Gram matrices. The full ambient chambers are not claimed compact or finite volume.

| Class | Original core / profiles | Sole connected marked component | Marked positions along its long arm | A1 context |
|---|---|---|---|---|
| 50 | 986 / `(1,1)` | D10 | distances 3 and 6 from the trivalent node | both absent and present |
| 55 | 1036 / `(2,2)` | D9 | distances 3 and 6 from the trivalent node | present |
| 62 | 1110 / `(1,1)` | D10 | distances 3 and 7 from the trivalent node | both absent and present |

All three arise from two length-one chains, with repeated first-contact profiles allowed. For a component-local weight rule assigning zero to disconnected marked pairs, the displayed orbit must receive weight at least one in the indicated context. The exact orbit keys and literal geometry are saved in `new_weight_obligations` and `unpaid_geometry`. This necessary charge statement is not asserted for rules that can place weight on disconnected pairs or depend on more of the ambient forest.

The [immutable star export](fresh_sqrt5_two_chain_star_closure_20261015_export.json) contains all 31 newly checked sealed marked classes: 29 initial classes and the two children, of which 28 are already paid and three are unpaid. Its [generator](../scripts/fresh_sqrt5_two_chain_star_closure_20261015_export.py) checks the bound source manifests and all marked-class keys; generation and saved replay pass. It supplies the next weight optimization without mutating earlier weights.

These three new obligations are the remaining output of this branch. They have not been incorporated into the frozen global weight optimization. Other initial attachment topologies, including branching components and several nodes of one component touching the core, remain outside the two-chain calculation.

Generation and saved replay pass. The source manifests bind all inherited contact, containment, arithmetic and weight dependencies. Independent review of the new all-length producer is being coordinated separately; this closure is not presented as independently recounted.

```sh
python3 scripts/fresh_sqrt5_two_chain_star_closure_20261015.py --verify
```
