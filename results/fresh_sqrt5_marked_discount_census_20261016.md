# Exhaustive marked-vertex census for baseline discounts

The [exact census](fresh_sqrt5_marked_discount_census_20261016.json) and [checker](../scripts/fresh_sqrt5_marked_discount_census_20261016.py) cover all **3,617,699 marked rank-21 finite-forest types**, up to automorphisms preserving the distinguished facet root. Exactly **374** have baseline at least 35:

| Baseline | Marked types |
|---:|---:|
|38|2|
|37|9|
|36|68|
|35|295|

Generation and fresh saved replay pass. The component catalogue and original baseline are inherited from the [marked-facet baseline review](fresh_sqrt5_marked_facet_barrier_20261016.md). The threshold for the selected 20-dimensional facet remains 35. This is a vertex-budget diagnostic, not a face-covering theorem or dimension bound.

## Completeness and corner data

There are 68 allowed connected finite types through rank 21 and 311,574 unmarked rank-21 forests. For a component define its integer deficit as `2*rank-baseline`. Every deficit is positive. A marked forest has baseline at least 35 precisely when the marked component's deficit plus its removed incident payment, together with all other component deficits, is at most 7. Exact rank/deficit multiset enumeration therefore isolates every high-cost marked type. A separate unbounded coin-change count covers all forests, including those below the threshold.

Every record identifies the marked component, its facet-root orbit, and the complete remaining component multiset. Shared corner catalogues retain the literal unordered root pairs and distinguish the facet root from the two corner roots. For each paid orbit they record the component rank, distance between the corner roots, intervening bond labels, distances from the facet root, distance to the nearest branch, distance to the nearest multiple bond, and the central node's degree when applicable.

The two baseline-38 types are `A1+D20` with the isolated root marked, and `D21` with its long-arm endpoint marked. The nine baseline-37 types are:

- `A1+A20`, `A1+B20`, and `A1+D16+H4`, with the isolated root marked;
- `A21` with an endpoint marked;
- `B21` with either endpoint marked;
- `D17+H4` with the long-arm endpoint of `D17` marked;
- `D21` with either a fork leaf or the node next to the long-arm endpoint marked.

The `B` cases matter: ordinary simply-laced paths occur inside their long tails. A discount restricted by the component names `A` and `D` leaves eight types above 35 and has maximum 37, regardless of how much its eligible payments are reduced.

## Exact discounted budgets

Set `delta=3/19`. A distance-two pair is ordinary here when both intervening bonds have label 3. The following rules subtract `delta` from each eligible unit payment; all remaining payments stay at the inherited baseline.

| Discounted pairs | Exact marked maximum | Violated minimal triangle cores |
|---|---:|---:|
| Distance-two pairs in `A` or `D` components only |37|3|
| All ordinary distance-two pairs |35|8|
| Ordinary distance-two pairs in components of rank at least 6 |35|0|
| All adjacent label-3 pairs |35|35|
| Adjacent label-3 pairs in components of rank at least 6 |35|0|

Each maximum is recomputed by an exact unbounded finite-forest dynamic program for all ranks through 21, then by every marked-root case. The certificate saves every rank value and all extremizing forests. For the two ordinary distance-two rules, `3/19` is the smallest uniform discount that reaches 35 when there is no unchanged obstruction: both baseline-38 maximizers contain exactly 19 eligible corners. Their discounted totals are exactly 35.

The long-component rules leave all 76 minimal triangle cores unchanged because their corner components have rank at most four. This does **not** establish triangle coverage: added common white roots can turn a small corner component into a long one. The separate [actual triangle witnesses](fresh_sqrt5_marked_discount_triangles_20261016.md) show that both long-component rules fail after such extensions, despite satisfying every vertex budget.

The useful proposal is therefore a discounted rule with carefully chosen compensation. The catalogue makes the marked contexts and available vertex slack explicit. Arbitrary connected white extensions and quadrilaterals must be covered anew; the old pruning of already-paid short corners is not assumed valid after any discount.
