The 42 saved mixed-contact cap proposals all retain one common exact spacelike escape, even if they are imposed simultaneously. The selected one-parameter refinement is also excluded for every real parameter on its positive-time branch. These are scoped infinite-volume certificates, not a nonexistence theorem for dimension 20.

The [exact checker](../scripts/fresh_h20_e7_replacement_escape_20261014.py) and [certificate](fresh_h20_e7_replacement_escape_20261014.json) rebuild the recipes in the [interior scout](fresh_h20_e7_causal_residuals_20261014_interior_scout.json) from the physical normals of the [87-wall chamber](fresh_h20_e7_intrinsic_caps_20261014_cusps.json). Both dual bases use physical norm-two roots; the scaled A2 vectors in the parent's face-coordinate basis are not substituted for those roots. Each positive time coefficient is reconstructed from the norm-two equation. The proof uses exact arithmetic in Q(√2,√79,√11,√3), with an additional square root handled by sign comparison and squaring.

The new rational point E has

`E² = 4374521734010421633 / 1000000000000000000 > 0`.

It pairs strictly negatively with all 85 retained walls after deleting current walls 61 and 75, and with every one of the 42 proposed cap normals. The inherited strict timelike center also satisfies all these inequalities. Consequently every subset of the proposed caps with this retained core still has infinite volume. This does not require the proposed caps to be mutually Coxeter, and makes no claim that they are. The [continuous checker](../scripts/fresh_h20_e7_replacement_escape_20261014_continuous.py) additionally verifies that E satisfies even the two deleted old walls strictly.

The particularly simple proposal 34 is

`r(2) = −ω₂−ω₄−ω₆−√2 λ₀−2λ₁₂+5T`,

where the indices refer to positions in the seven E7 weights and thirteen physical A11+A2 weights, and `T²=−2`. This normal has square two and cuts the former escape. It nevertheless retains the new E. Cutting one known escape is therefore insufficient evidence of finite volume.

More strongly, consider the full continuous family

`r(x) = V−xλ₁₂+sqrt((71+x²)/3)T`,

`V = −ω₂−ω₄−ω₆−√2 λ₀`, for any real x.

The [continuous certificate](fresh_h20_e7_replacement_escape_20261014_continuous.json) checks that V, λ₁₂ and T are mutually orthogonal, with squares `148/3`, `2/3` and `−2`. Thus every displayed r(x) has square two. With `C=V·E`, `a=−λ₁₂·E` and `B=−T·E`, it proves exactly

`C>0`, `B>0`, `D=B²−3a²>0`, and `71D/3−C²>0`.

The reverse Lorentz Cauchy inequality then gives

`E·r(x) = C+ax−B sqrt((71+x²)/3) ≤ C−sqrt(71D/3) < 0`

for every real x. The upper bound is approximately −0.350289; that decimal is explanatory and is not used by the checker. In particular, the refined contact at `x²=86/11`, `t²=289/11` retains this escape. No finite coherent collection consisting only of this family and any subset of the old 87 walls can close the chamber. Other E7 profiles, independent source coefficients and different retained cores remain outside this statement.

Both exact checkers pass generation and saved replay with `--verify`. The numerical escape scouts are exploratory artifacts; validity rests on the rational point and exact inequalities in these certificates.
