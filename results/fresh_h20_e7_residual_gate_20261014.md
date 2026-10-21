The corrected all-residual scout has **144 negative and 100 positive numerical relaxations** among 244 tested fixed E7 contact profiles. All 144 negative cases now have exact Lorentz-dual certificates, each proving `r·E<−1/100000` for every allowed intrinsic residual. The 100 positive relaxations remain numerical proposals; no Coxeter cap or finite-volume construction is certified by them.

The [scout](../scripts/fresh_h20_e7_residual_gate_20261014_scout.py) previously formed Gram inverses and physical pairings in floating point. This left small nonzero coefficients in constraints that should be exactly constant, producing unreliable optimization bounds. The new [exact coefficient helper](../scripts/fresh_h20_e7_residual_gate_20261014_exact_data.py) forms both inverse matrices over the rationals and every physical pairing in Q(√2,√79,√11). It converts only the final field values to floats, after exact-zero detection and a 90-digit evaluation with a sign check. There are exactly 396 zero face-wall coefficients, 20 zero time-wall coefficients and 147 zero E7-weight coefficients. The helper also verifies dual-weight identities, orthogonality and time orientation.

The underlying [87-wall branch](fresh_h20_e7_intrinsic_caps_20261014_cusps.json) and [complete H13 face](fresh_h20_e7_intrinsic_caps_20261014.json) are unchanged. Write an intrinsic residual as `y=−Σ d_i λ_i+tT`, with `T²=−2`, `d_i≥0`, and let J be the rational A11+A2 inverse Gram matrix. For a fixed E7 coefficient vector c, put `w=Σ c_j ω_j` and `m²=cᵀIc−2`, where I is the E7 inverse Gram matrix. A norm-two candidate `r=−w+y` satisfies

`dᵀJd+m²=2t²`.

Acuteness against the complete H13 source forces y into its future causal cone: the source's full projective census consists of 44 finite and two null rays. Thus t is nonnegative. The scout replaces the displayed equality by `dᵀJd+m²≤2t²`, a convex relaxation. It also imposes all 87 old-wall acuteness inequalities, but not their discrete angle conditions.

The fixed profiles comprise one or two nonzero E7 contacts with labels 3, 4 or 6, and three nonzero contacts all with label 3. Of the 245 formal choices, the single label-3 contact at E7 position 5/current wall 31 has `w²=3/2`, requiring a spacelike residual of square 1/2. The complete face cone rules that out. The mass-zero profile at position 3/current wall 15 is retained in the scout. All other 244 programs terminated with numerical status `optimal`; their values and residual coordinates are in the [corrected scout output](fresh_h20_e7_residual_gate_20261014_scout.json).

The [exact dual checker](../scripts/fresh_h20_e7_residual_gate_20261014_dual.py) uses no optimizer or floating-point acceptance. Its [certificate](fresh_h20_e7_residual_gate_20261014_dual.json) lists all 144 proved profiles. It uses rational multipliers in the [recipe file](fresh_h20_e7_residual_gate_20261014_dual_recipes.json) and exact arithmetic in Q(√2,√79,√11,√3). Of these proofs, 129 need only the unconstrained Lorentz bound; 15 use old-wall inequalities as well. Sixty-nine proved profiles use only labels 3 and 4.

Here is the checked dual inequality. Set `C_ij=(λ_i,N_j)`, `τ_j=(T,N_j)`, `P_ij=(ω_i,N_j)`, `a_i=−(λ_i,E)`, `B=−(T,E)` and `A=−(w,E)`. Old-wall acuteness reads `−Cᵀd+τt−Pᵀc≤0`. For rational nonnegative old-wall multipliers λ and nonnegative-coordinate multipliers μ, define

`ℓ=a+Cλ+μ`, `h=B+τ·λ`, `K=A+(Pλ)·c`.

Subtracting λ times the old-wall inequalities and adding `μ·d≥0` gives

`r·E ≤ K+ℓ·d−ht`.

Every saved recipe satisfies `h>0` and `D=h²/2−ℓᵀJ⁻¹ℓ>0` exactly. The reverse Lorentz Cauchy inequality therefore gives

`r·E ≤ K−sqrt(m²D)`.

For each profile the checker proves either `K+1/100000<0`, or `m²D−(K+1/100000)²>0` with a nonnegative shifted constant. This establishes the stated strict margin for all intrinsic residuals, including those not lying on any enumerated face edge or vertex ray. The multipliers were suggested by the solver, then moved toward a strict dual interior and rationalized; their validity is entirely checked by the exact inequalities.

Useful **numerical** positive proposals include single contact at current wall 65. For labels 3, 4 and 6, respectively, the optimizer has only face-contact positions 0 and 12 appreciably nonzero, approximately `(d0,d12,t)=(1.43151,0.65216,1.51834)`, `(1.65563,0.96238,2.13660)` and `(1.68631,1.21577,2.60692)`. Their relaxed escape objectives are about 1.20, 1.43 and 1.50. Mixed contacts at current walls 18 and 65 also survive the relaxation. These coordinates have been passed to the constructive lane; their nonintegral contact magnitudes still need full Coxeter angle gates.

The 144 certificates exclude only their listed fixed E7 profiles. They are not an all-label or all-profile completion obstruction. The exact checker also does not interpret a positive solver optimum as an actual cap. No prior frozen source, certificate or report was modified; the explicitly unfrozen numerical scout was corrected.

Exact generation and saved replay both pass:

```sh
python scripts/fresh_h20_e7_residual_gate_20261014_dual.py
python scripts/fresh_h20_e7_residual_gate_20261014_dual.py --verify
```
