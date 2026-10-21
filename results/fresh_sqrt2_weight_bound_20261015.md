# Independent weights for the Q(√2) dimension-bound argument

The [checker](../scripts/fresh_sqrt2_weight_bound_20261015.py) and
[certificate](fresh_sqrt2_weight_bound_20261015.json) independently verify
the weight accounting for the proposed compact arithmetic Q(√2) bound 20.
The implication requires the separate exhaustive quadrilateral-star and
external-attachment proof. This package alone does not establish that
classification or a global theorem.

Give an angle base weight 1 when its two marked nodes have distance at most
2 in the spherical vertex diagram, and 0 otherwise. Add the following
weights, always using the **entire connected component** of that diagram:

| Component | Additional marked pair | Weight |
|---|---|---|
| F4 | the two endpoints of its 3–4–3 path | 1 |
| D5, arms 1,1,2 | the long leaf and either short leaf | 1/2 each |
| E6, arms 1,2,2 | the near node on one long arm and far node on the other | 1/2 each |
| E7, arms 1,2,3 | near node on the length-3 arm, far node on the length-2 arm | 1/4 |
| E8, arms 1,2,4 | near node on the length-4 arm, far node on the length-2 arm | 1/4 |

All these marked pairs have distance 3, so no base weight is double-counted.
The checker tests invariance under every automorphism of each model diagram.
It reconstructs the four finite corners of each of the seven sealed
finite-label cores directly from their labelled graphs; their quadrilateral
weight is at least 1. It separately reconstructs the three survivors of the
[pendant-tree lemma](fresh_sqrt2_dotted_tail_20261015.md), obtaining totals
1,1,2. It imports no producer arithmetic or graph routines.

A triangular face has a minimal indefinite diagram containing all three
marked nodes. Otherwise that diagram would occur in a finite corner.
By the Q(√2) Lannér classification, it is a connected order-3 diagram or
a cycle of order 4 or 5. The first has at least two nonright edges. In the
cycle cases, the three positive integer gaps between marked nodes sum to
4 or 5, so at least two gaps have length at most 2. Those paths belong to
their respective finite corners. Thus every triangle receives at least 2.
The checker enumerates all cycle markings as a separate check of this step.

For a finite path of rank r≥2, the base weight is 2r−3. For a finite D/E
tree it is 2r−2, by counting edges and pairs of neighbors. Rank 1 contributes
zero. The complete component totals after correction are:

    F4: 6;  D5: 9;  E6: 11;  E7: 49/4;  E8: 57/4.

Thus the deficit from 2r is at least 2 except for D5/E6, where it is 1,
and E7/E8, where it is 7/4. At total vertex rank n≥9, a single component
cannot be one of these exceptions. If there are two or more components,
their deficits sum to at least 2. Consequently every vertex total is at
most 2n−2. This is an all-rank argument; the checker also verifies ranks
through 60 and the maximum over every rank partition in that range.

If the exhaustive star classification proves that every quadrilateral
receives at least 1 under these entire-component rules, then every k-gonal
two-face receives at least 5−k. Double-counting gives

    sum(vertex weights) ≥ 5 f₂ − binom(n,2) f₀.

The symmetric nonnegative h-vector gives the strict average lower bound
n(n−6)/8 for even n and (n−1)(n−5)/8 for odd n, as derived in the frozen
[independent cubic-field review](fresh_compact_global_bound_20261015_review.md).
At n=21 it is strictly greater than 40, whereas 2n−2=40. For n=21+2u and
n=22+2u, the non-strict lower-minus-upper gaps are respectively

    u²/2 + 5u,       u²/2 + 11u/2 + 2.

They are nonnegative for every integer u≥0, and strictness excludes equality.
The conditional dimension consequence is therefore n≤20.

The original face correspondence, arithmetic criterion, finite diagrams,
and weighted-face framework are the published inputs from
[Vinberg's paper](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).
The remaining obligation belongs to the global classification package:
all initially possible minimal opposite pairs, all cross contacts and
shared whites, and all arbitrarily large external white components must
reduce to the sealed cores and pendant survivors used here. A bounded
enumeration of only the initial union would not suffice.

Replay with `python3 scripts/fresh_sqrt2_weight_bound_20261015.py --verify`.
