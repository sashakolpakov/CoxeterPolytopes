# Zero baseline discounts still exceed budget 40

The existing global dimension-21 payment family cannot reach budget 40 when both baseline discounts are fixed to zero. An exact nonnegative combination of 1,709 selected face, vertex and model constraints proves

`T ≥ 43232664297439 / 1055499993667 > 40`

(approximately 40.9594169179). This is a lower bound for this restricted payment family, not a geometric non-existence theorem or a claim that the exact optimum was recovered.

The test preserves all 32,554 variables, their literal normalization, and the 4,762 previously reconstructed global vertex cuts. It fixes the distance-two and adjacent deductions to zero. The 1,979 inherited baseline-zero quadrilateral rows are evaluated in all 227,408 common spectator contexts. The remaining quadrilaterals already have baseline payment at least one. Every graph in the complete 996-state triangle catalogue has baseline at least two. Since every additional payment is nonnegative, these latter obligations need no new compensation when both discounts vanish.

The fixed-budget solve found its finite constraint set infeasible. A subsequent minimum-budget solve returned approximately 40.9703439891; rational weak-duality recovery produced the smaller exact bound above. Both numerical solves finished within the single five-minute optimization budget. No numerical optimum is used as a proof.

The [source](../scripts/fresh_sqrt5_global_zero_discount_20261017.py) rechecks the exact multipliers, their nonnegative combined column coefficients and the resulting rational bound. Every vertex appearing in the dual is independently reassembled from component costs and recounted as literal unordered root pairs. Every face term is matched to the reconstructed complete-context obligations. The [saved certificate](fresh_sqrt5_global_zero_discount_20261017.json) records all terms and their face or forest metadata. Generation and a fresh saved replay pass:

```sh
python scripts/fresh_sqrt5_global_zero_discount_20261017.py
python scripts/fresh_sqrt5_global_zero_discount_20261017.py --verify
```

Thus simply removing both discounts from the expanded bonus family does not repair the attempted bound-20 proof. This does not exclude different baseline reallocation, other payment families, or stronger global incidence arguments. The frozen compact arithmetic bound-21 candidate is unchanged.
