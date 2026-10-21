# All-real affine endpoint gate for the original 62 edge germs

Every single Coxeter wall that turns one of the original 15-wall chamber's 62 open finite-corner germs into a **bounded-link ideal endpoint** strictly misses the saved spacelike escape. This covers all real parameters and all integer labels on the remaining old-wall contacts. It does not exclude adding intermediate walls or completing different edge supports.

The [certificate](fresh_h9_h8_disjoint_walk_20261014_affine_gate.json) and [producer](../scripts/fresh_h9_h8_disjoint_walk_20261014_affine_gate.py) cover all **14 proper-affine profiles**, using 50 interval boxes of depth at most three. Generation and saved replay both pass.

## Local profiles and the complete parameter family

The [frozen local affine census](fresh_h9_h8_disjoint_walk_20261014_affine.json) supplies the 14 profiles. Its scope is one new wall added to the eight positive-definite normals of an old edge. A bounded Euclidean eight-dimensional link with these nine facets is a simplex, so its Gram diagram must be irreducible affine. The rank-eight affine classification permits only labels 3 and 4, new-node degree at most four, and a tree or the affine `A8` cycle. The census checks 6,056 assignments on 22 eligible supports; the other 40 supports already contain an incompatible label-5 edge. This completeness argument is inherited, not independently enumerated here.

Normalize the eight support roots to square two and write their Gram and inverse as `G,J`. A contact vector `c` has entries `1` or `√2` at its nonright contacts. The producer reconstructs and checks

`w=Jc>0`, `cᵀJc=2`, and `p=−Σ w_i N_i`, hence `p²=2`.

It also checks the normalized Gram against physical old-root products. The projected weights lie in `Q(√5)+√2 Q(√5)`; the coefficient of `√2` in the Schur cost vanishes exactly.

Let `X` be the old finite endpoint and `Y` its saved outward edge direction. Set `U=X/√(−X²)` and `V=Y/√(Y²)`. Then `U²=−1`, `V²=1`, `U·V=0`, and both vectors are perpendicular to the support. The complete outgoing ideal family is

`r=p+h(U+V)`, `h>0`.

Indeed, the positive affine relation is `r+Σw_i N_i=r−p`. Pairing this null vector with any strict future center of the proposed extension is negative, so it is future null. Of the two future null directions in the complementary Lorentz plane, `U−V` has positive pairing with the released old wall. Since `p` has nonnegative pairing with that wall, it cannot give an acute extension. Thus the outgoing direction is forced. The limit `h=0` has no strict interior because its positive relation vanishes, but the checker includes it conservatively.

## Exact all-parameter exclusion

For each old nonsupport root, normalized to square two, and for the literal escape `E`, the pairing takes the affine form

`r·N_j=b_j+h(a_j+c_j)`.

The coefficient names in the JSON refer to this formula, not to the local contact vector above. All coefficients are reconstructed from the physical vectors. Exact outward 128-bit dyadic intervals establish `a_j≤0` and `a_j+c_j<0`; for the escape they also establish `a_E<0<c_E`.

Consequently every pairing decreases with `h`. A certified integer upper bound `H` makes `r·E<0` for all `h≥H`. The remaining closed interval `[0,H]` is exhaustively split until each box either has `r·E<0`, or an old-wall pairing lies outside the allowed Coxeter set. No floating-point value selects an accepted sign or a stopping condition.

The old-wall gate isolates labels `2,…,15` and retains the entire interval `(-∞,−2cos(π/16)]`. This interval includes every label at least 16, all parallel contacts and all dotted contacts. It is an overapproximation, so the result has no finite-label cutoff. The 50 boxes have 32 terminal boxes: ten are rejected by the escape inequality and 22 by an old-wall angle gap. Every unbounded parameter tail is separately rejected by its certified escape bound.

The exact scalar engine is reused from the [earlier universal endpoint checker](../scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py), specialized to residual norm zero. This is a new coefficient reconstruction and replay, not an independent implementation of interval arithmetic.

## Scope

The literal escape has positive square `−6837/256+6543√5/512` and is strict against all 15 original walls. This result proves it also remains strict against any admissible single-wall proper ideal completion of these 62 supports. It strengthens the earlier bounded affine-event scout for this purpose.

The theorem concerns the original 15-wall seed. It does not classify ideal completions on the later 71 supports, ends with unbounded links, walls that first create intermediate finite vertices, or other normal directions. The separate finite-endpoint gates and the exact quartet remain distinct artifacts. No finite-volume nonexistence theorem for arbitrary completions is claimed.

```sh
python scripts/fresh_h9_h8_disjoint_walk_20261014_affine_gate.py --verify
```

Source hashes bind the seed, exact frontier, local affine census, arithmetic helper and interval engine. No frozen predecessor was changed.
