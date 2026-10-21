# Frontier and bounded contact scout from the nine new vertices

This stage retains the exact [19-wall quartet chamber](fresh_h9_h8_disjoint_walk_20261014.md). It introduces no additional walls. The nine newly created finite vertices have **71 open germs**, certified exactly; the 56 still-open germs at older vertices are not rescanned.

The [frontier certificate](fresh_h9_h8_disjoint_walk_20261014_next_frontier.json) uses the correct field `Q(√5,√K)`, `K=6+16√5/5`. For every germ it reconstructs the eight-wall positive-definite support and its forward direction `D`, checks `P·D=0`, and gives a strictly feasible spacelike continuation beyond the opposite null limit. At that limit precisely the eight support roots remain active, so the link is unbounded. This is a complete queue at these nine finite vertices, not a complete projective-ray census.

The local spherical profile census gives **1,065 finite profiles and 50 symbolic isolated dihedral tails**. The normalized support Grams still lie in `Q(√5)` and have only labels 2, 3 and 5. The symbolic tail parameters remain unrestricted in the exact profile statement.

The [numerical scout](fresh_h9_h8_disjoint_walk_20261014_next_scout.json) samples the same explicitly listed finite-anchor and dihedral-tail labels as the previous stage. Its 1,515 sampled profiles produce 532,152 real anchor roots and **21,024 deduplicated accepted forward-endpoint proposals**. None cuts the saved escape. The best numerical escape pairing is approximately `−0.44951922`; this is a diagnostic, not an exact lower bound or an all-label exclusion.

The final implementation computes physical pairings and Lorentz-plane normalization at 100 decimal digits before solving the numerical quadratic events. It evaluates a cap through its wall-pairing formula `Az+B+C/z`, avoiding cancellation in a directly boosted physical normal. For `z>1`, its endpoint is `P+tD` with

`t=((z²−1)/(z²+1))√(−P²/D²)`.

The candidate checks require forward `t`, a strictly negative endpoint norm and strict nonsupport wall signs. A guard rejects numerically unresolved near-null endpoints; no such skip occurred in the final run. Remaining wall angles are compared against labels 2 through 256 with tolerance `10⁻⁸`, or a dotted gap of `10⁻⁷`. The unbounded tail labels are sampled only at `6,7,8,10,12,16,24,32,64`; anchor labels are `2,3,4,5,6,7,8,9,10,11,12,16,20,24,32,64`.

The full proposal bank is [deterministic gzip](../data/enumerations/fresh_h9_h8_disjoint_walk_20261014_next_scout_bank.json.gz). Its hash is recorded in the numerical manifest. There is no exact new-cap certificate and no deeper wall-count extension in this stage. Continuous parameters and unsampled labels are not excluded by this scout; proper ideal endpoint profiles are also outside this finite-profile stage.

Generation and saved replay pass for both files, as do syntax, source hashes and report links:

```sh
python scripts/fresh_h9_h8_disjoint_walk_20261014_next_frontier.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h9_h8_disjoint_walk_20261014_next_scout.py --verify
```
