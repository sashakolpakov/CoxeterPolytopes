# Center/short-arm omission scout: a rejected refinement

The [exact scout](../scripts/fresh_nikulin_center_omission_scout_20261016.py) and [saved output](fresh_nikulin_center_omission_scout_20261016.json) test replacing the previous E7/E8 omission pair by the center and short-arm node, numbered `(0,1)`. This change does **not** improve the unrestricted finite-volume dimension bound. It introduces a cheap coincident-row obstruction, and the proposed far-height inequality fails by `4/3`.

For either ordinary E7 or E8, put both endpoint contact rows at node 1, with contact value one. Deleting the center leaves that node isolated with Gram entry two. The endpoint costs are therefore `(1/2,0)`, while the two normalized contact distributions coincide and `W=0`. In the saved output these are E7 model 972, case 17, and E8 model 981, case 15; all indices are zero-based. Their exact convex certificates have weights `(1,1)`, zero multipliers, and zero slacks.

Consequently, the proposed inequality

\[
W+\frac43(C_0+C_1)\ \geq\ \frac{2}{M^2}
\]

already fails at equal endpoint heights, `M=1`: its left side is `2/3`, versus right side `2`. It also fails in the intended far range: scale the first row by `20/19`, obtaining `M=20/19`, `C0=200/361`, `C1=0`, and still `W=0`. The left side is `800/1083`, whereas the right side is `361/200`. No global extension or finite-volume realization is asserted. Avoid-node variants are not needed for this failure.

The scout uses exact rational convex certificates to cover 1,002 affine models and 62,552 continuous support domains. Classical affine types are included through affine rank 162, together with the exceptional types and the modified/avoid cases used by the inherited proof. For E7/E8, `(0,1)` is used when it avoids the forbidden node; otherwise the previous pair is retained. The target is `2r/9` for an ordinary factor and `2(r−1)/9` for a modified factor. Each endpoint has budget two, and the near-height ratio cap is `20/19`.

There are 8,625 domains with positive loss upper bound. Exact ordered-cost domination reduces them to 15 records. Packing those records under endpoint budgets `(2,2)`, with at most one modified factor, gives the relaxed maximum `80/9`. Its witness consists of four ordinary B2 records with costs `(0,1/2)` and loss `4/9`, and four ordinary E8 records with costs `(1/2,0)` and loss `16/9`. The packing forgets total rank and compatibility between factors. The value `80/9` is therefore a maximum of this numerical relaxation, not a certified affine-rank-162 Gram configuration or a finite-volume polytope. It is already much worse than the frozen loss-four rule.

The far-height candidate was checked on 525 changed E7/E8 domains with retained total cost below `3/2`; its exact minimum margin is `−4/3`. The explicit ordinary short-arm rows above attain that margin, so failure does not depend on an unattained continuous lower bound.

The finite rank range is deliberate: rank 162 is the next possible exclusion case under the previously reviewed unrestricted bound 163. This scout provides no all-rank theorem and no improved dimension bound. It imports the frozen exact convex helper read-only; it has generation and saved-output replay, but no separate independent arithmetic review. The existing bound 163 is unchanged.

Replay:

```sh
python3 scripts/fresh_nikulin_center_omission_scout_20261016.py --verify
```

Together with the [full rank-162 Gram obstruction and strict-null escape](fresh_nikulin_rank162_local_obstruction_20261016.md), this leaves incidence and finite-volume completion as useful additional constraints. The Gram obstruction has all present wall pairings below 14 but an explicitly infinite-volume chamber; the present scout shows why simply moving the E7/E8 omissions to the center is insufficient. Adaptive choices after seeing all contact rows remain untested.
