# Independent exact review of the complete-queue weight rules

The separate [checker](../scripts/fresh_sqrt5_complete_face_weights_20261015_review.py) reproduces the rational face payments and vertex maxima for all three rules in the [producer report](fresh_sqrt5_complete_face_weights_20261015.md). Generation and saved replay pass for each dimension.

| Dimension | Face rows, with repetitions | Exact spectator cases | Maximum vertex weight | Strict face threshold |
|---|---:|---:|---:|---:|
| 25 | 2,709 | 21,897 | 59 | 60 |
| 26 | 2,905 | 24,609 | 31,017/500 | 65 |
| 27 | 3,110 | 27,523 | 67,009/1,000 | 143/2 |

The saved independent certificates are [dimension 25](fresh_sqrt5_complete_face_weights_20261015_review_n25.json), [dimension 26](fresh_sqrt5_complete_face_weights_20261015_review_n26.json), and [dimension 27](fresh_sqrt5_complete_face_weights_20261015_review_n27.json). These certify weights on the stated inputs; a complete dimension theorem additionally requires the reviewed arithmetic attachment queues and the global face-coverage composition.

## Independent reconstruction and scope

The checker imports no optimizer or weight-producer arithmetic. It reuses the frozen [separate labelled-tree/Fraction implementation](fresh_sqrt5_all_a1_cross_corrected_20261015_review.md), reconstructs whole spherical components and marked or rooted orbits directly from each literal graph, and checks all variable indices and nonnegative grid coefficients. Intrinsic component recognition treats A2 as I2(3) and B2 as I2(4); the earlier naming defect is not inherited.

Every one of the completed finite queue's 2,070 singular states is accounted for before filtering by the target dimension. The included counts are 1,488, 1,623 and 1,767. Expanded singular states, previously paid states and boundary states are all retained. There is no sealed-only or old-payment pruning. The earlier 2,397 literal inputs are separately mapped to their reviewed parameter-marked representatives; every label, marked-node bijection, dotted ratio and dimension filter is checked. Extra conservative inputs need not be realizable to make the verified inequalities sufficient.

For a displayed corner of rank r in dimension n, the external spherical rank is n−r. Every possible exact number e of external A1 components is tested: 0≤e≤n−r, excluding only n−r−e=1, which cannot be a forest of non-A1 components. The entire corner's A1 count is its internal count plus e. The checker recomputes each of the four angle weights, including the rational division in disconnected-pair weights, and requires their sum to be at least one.

If a full vertex has a>0 A1 components, each pair joining one of them to a specified root of a non-A1 component receives u_a/a. There are exactly a such pairs for that root. Their total is u_a, which is checked explicitly before using additive component costs. No A1–A1 pair is charged. Thus a fixed exact A1 count reduces the vertex maximum to an ordinary finite-forest rank problem. The reviewer builds the finite Dynkin types independently, sums every root-pair and rooted-node multiplicity, and uses a type-first unbounded rank recurrence, different from the producer's rank-first separation. It also checks every saved maximizing forest. All feasible exact A1 counts are included, including the all-A1 vertex; n−1 is the sole impossible count.

The three exact maxima are feasible upper bounds, not claims of optimality. The 1/1000 upward rounding is checked coefficient by coefficient, followed by all face and vertex inequalities. No numerical solver status is used as evidence. The H4-only triangle baseline is nonnegative and agrees with the independently reviewed 76 triangle markings; extra bonuses cannot lower those payments.

This is independent implementation and proof-scope review within the same collaborative research effort. It does not repeat the finite queue's field arithmetic or certify its geometric completeness. Nor does it assume that older closure payments transfer to new weights: those branches must be evaluated directly by the separate coverage binder.

```sh
python3 scripts/fresh_sqrt5_complete_face_weights_20261015_review.py --dimension 25 --verify
python3 scripts/fresh_sqrt5_complete_face_weights_20261015_review.py --dimension 26 --verify
python3 scripts/fresh_sqrt5_complete_face_weights_20261015_review.py --dimension 27 --verify
```
