# Complete last-block label coverage with three fixed dense blocks

**Star case0 has no compact H9 realization when its first three blocks have all own contacts `(3,3)` and the last internal label is at least4.** Every last own label may be arbitrary, and all omitted meeting labels and dotted distances remain free. This is a contact-family theorem, not a whole-incidence exclusion.

The [logical checker](../scripts/fresh_h9_star_lastblock_complete_20261014.py) and [certificate](fresh_h9_star_lastblock_complete_20261014.json) bind the three exact arithmetic families and their separate independent reviews:

- [Last internal4 or5](fresh_h9_star_lastblock_20261014.md), with [independent review](fresh_h9_star_lastblock_20261014_review.md): arbitrary own labels, reduced to290 profile pairs and77,952 residual cases.
- [Shared single π/3 contacts, internal at least7](fresh_h9_star_shared_high_20261014.md), with [independent review](fresh_h9_star_shared_high_20261014_review.md):7,840 continuous affine cases.
- [Opposite single π/3 contacts, internal at least7](fresh_h9_star_opposite_high_20261014.md), with [independent review](fresh_h9_star_opposite_high_20261014_review.md):47,040 continuous quadratic cases.

For last internal labels at least6, the remaining reduction is elementary and exhaustive. If any of the four own labels is at least4, one of four literal reviewed strong-contact placements excludes it. Otherwise all own labels lie in `{2,3}`. If either anchor has dense profile `(3,3)`, one of two literal reviewed two-dense placements excludes it.

At internal6 the only Lannér profile in `{2,3}²` is `(3,3)`, so that label is completely covered. At every internal label at least7 the only profiles are `(2,3)`, `(3,2)` and `(3,3)`: `(2,2)` is never Lannér, while the reciprocal sums for the two single profiles are strictly below1 precisely beyond6. After the dense profiles are removed, the four ordered single-profile pairs are exactly the shared/opposite families and their checked root swaps. The certificate verifies each literal match, rather than assuming an unproved incidence symmetry.

The first three internal labels are automatically at least4 from their `(3,3)` own missing triples. Their possible larger values are already included in the component theorems and reviewed high-block gates. The logical checker introduces no new arithmetic engine or cutoff. It requires each independent review to bind the exact current primary certificate before producing its result.

## Resumable search

The [search](../scripts/fresh_h9_star_lastblock_complete_20261014_search.py) appends one complete-family predicate, giving159 bound conditions, and resumes the105-visit cursor. [Batch006](fresh_h9_star_lastblock_complete_20261014_search_batch_006.json) reaches another mandatory-label/bridge survivor after81 additional visits:186 cumulative, with nine pending frames. Every earlier sibling cursor is retained; the explored prefix is not restarted.

The new state changes cut leaf16's own profile on `(9,10)` to `(3,4)`. All other own profiles are `(3,3)`, and the four internal labels are `(4,4,4,5)`. This lies outside the fixed-first-three-own hypothesis. It remains a necessary label assignment: omitted meeting contacts and dotted distances have not been completed to a Lorentzian Gram matrix. No compact polytope or whole-incidence closure is claimed.

A concrete next step is to audit which fixed own labels enter the residual arithmetic and which enter only monotone vertex failures. The core arithmetic uses `S8 + {11,17,14}`; this may allow some noncore contacts to increase without a new polynomial bank. Such a widening is not part of the present theorem.

Generation and saved replay pass using:

```sh
python scripts/fresh_h9_star_lastblock_complete_20261014.py --verify
python scripts/fresh_h9_star_lastblock_complete_20261014_search.py --batch 6 --verify
```

The independent reviews reconstruct the three arithmetic families separately. The new composition receives a source/proof coverage check; it is not a fourth duplicate of those large determinant censuses.
