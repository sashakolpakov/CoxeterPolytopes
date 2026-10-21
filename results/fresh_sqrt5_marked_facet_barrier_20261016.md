# The unchanged marked-facet baseline stops at 38

**The proposed facet model needs another change before an LP is worthwhile.** Its required dimension-20 threshold is **35**, while the existing nonnegative baseline already has exact maximum **38**. The new central-E8 incidence inequality leaves both maximizing marked vertex types unrestricted.

The [exact checker](../scripts/fresh_sqrt5_marked_facet_barrier_20261016.py) and [certificate](fresh_sqrt5_marked_facet_barrier_20261016.json) use the frozen finite component catalogue and baseline, enumerate all 724 rooted component cases, and perform exact finite-forest optimization through rank 21. No LP is run.

For a simple polytope of dimension `m`, the symmetric nonnegative h-vector gives

\[
f_0=\sum_jh_j,\qquad
f_2=\sum_j\binom{m-j}{2}h_j.
\]

After pairing `j` and `m−j`, the minimum coefficient is 45 when `m=20`, and 50 when `m=21`. The ratio is strictly larger than this minimum because `h₀=hₘ=1` and their coefficient is larger. A corner-weight rule paying at least `5−k` on each k-gonal two-face therefore has strict vertex-average lower bounds

\[
5\cdot45-\binom{20}{2}=35
\quad\text{and}\quad
5\cdot50-\binom{21}{2}=40,
\]

respectively. Thus the facet threshold is 35; the original whole-polytope threshold 40 is recovered by the same normalization.

The inherited baseline pays one unit for a connected root pair at distance at most two. H4 has its existing extra unit on the 5-bond and on the pair at distance three. Disconnected pairs receive zero. Marking a facet root `e` removes every pair containing `e`. For a marked component, its contribution is its full baseline total minus the sum of all pairs incident to `e`; the remaining components contribute their ordinary totals. Exact finite-forest dynamic programming therefore covers every marked rank-21 vertex diagram.

The maximum is **38**, attained precisely by these two marked types:

| Ambient vertex diagram | Facet root | Facet baseline total |
|---|---|---|
| A1 + D20 | The isolated A1 | 38 |
| D21 | The far endpoint of its long arm | 38 |

Neither type has a root deletion yielding the central-E8 support `A1+A2+A4+D5+E8`. Both sides of the new [narrow-facet incidence inequality](fresh_sqrt5_narrow_e8_ports_20261016.md) are therefore zero at a frequency concentrated on either marked type. Such a frequency is feasible in that incidence relaxation. The exact maximum remains 38 after imposing the entire weighted inequality, not just its coarser 4/5 consequence.

Any added nonnegative corner bonuses can only increase these two costs. Consequently, retaining this baseline and adding only those bonuses and the central-E8 incidence restriction cannot force the required upper bound 35. This is a model obstruction, not a realized-polytope claim. The baseline itself is not asserted to cover every required face payment.

For comparison, D5 + 2E8 has marked-facet baseline maximum **34**. Its problematic role in the earlier models comes from the additional payments, while this facet obstruction already occurs in the baseline.

Useful next changes are a different baseline or constraints on the actual marked D20/D21 transitions. Preserving the facet root through an edge can impose more information than matching the unmarked endpoint type.

Generation and saved replay passed. Source, certificate and report are frozen.

```sh
python3 scripts/fresh_sqrt5_marked_facet_barrier_20261016.py
python3 scripts/fresh_sqrt5_marked_facet_barrier_20261016.py --verify
```
