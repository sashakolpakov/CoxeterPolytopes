All 144 fixed-profile Lorentz dual witnesses pass a separate reconstruction in the physical 21-dimensional ambient space. The [review checker](../scripts/fresh_h20_e7_residual_gate_20261014_physical_review.py) and [certificate](fresh_h20_e7_residual_gate_20261014_physical_review.json) do not import the primary coordinate builder, optimization model or dual checker. They reuse the existing exact multiquadratic arithmetic engine; this is not an independent arithmetic implementation.

For each record in the [primary certificate](fresh_h20_e7_residual_gate_20261014_dual.json), the review constructs

`X = E − Σ λᵢNᵢ − Σ μⱼAⱼ`

directly from the physical normals and positive rational multipliers. Here Aⱼ denotes a physical corner wall of the completed H13 face. If r is acute with these retained walls, then `r·E ≤ r·X`.

The review subtracts X's E7 projection to obtain U. It independently verifies that U is orthogonal to all seven E7 roots, future timelike, and has square equal to the negative of the primary dual quantity. With `r=−w+y`, the complete finite-volume H13 face cone forces the residual y to be future causal. For a norm-two r, `y²=−m` where `m=w²−2`. Therefore

`r·E ≤ −w·X − sqrt(m(−U²))`.

The direct physical computations reproduce all 144 primary constants and masses, then prove the strict upper bound `−1/100000` by exact sign comparisons and squaring. The 129 bounds with no wall multipliers and the 15 bounds using retained-wall inequalities are both covered.

This validates the coordinate signs and multiplier directions by another representation. The completeness and causal-cone property of the H13 face remain inherited from its existing exact construction and independent review. The conclusion is limited to the listed profiles; it does not exclude the 100 positive numerical relaxations or unrestricted Coxeter labels.

Generation and saved replay both pass:

```sh
python3 scripts/fresh_h20_e7_residual_gate_20261014_physical_review.py
python3 scripts/fresh_h20_e7_residual_gate_20261014_physical_review.py --verify
```
