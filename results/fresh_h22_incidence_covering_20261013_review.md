# Review of the covering reformulation and missed-cusp certificates

**PASS.** The parent agent reviewed the [covering argument](fresh_h22_incidence_mixed_height_20261013_covering.md), its [stdlib integer checker](../scripts/fresh_h22_incidence_mixed_height_20261013_covering.py), and the [separate π/7 checker](../scripts/fresh_h22_incidence_prism_20261013_uncapped.py). Both frozen `--verify` commands were then run afresh and passed. The latter imports the separate Fraction audit's arithmetic; this review does not provide a third implementation of that arithmetic.

The all-point covering equivalence is correct under the stated retained-source and two-sided-height hypotheses. Normalizing the old timelike coordinate reduces upper and lower feasibility to scalar bounds on z. A bound exceeding one produces a feasible spacelike point, which can be made strict. Bounds at most one for every old timelike point also control old-null limits by homogeneity and continuity. Consequently the finite polyhedral cone lies in the future Lorentz cone, and its projective polytope has finitely many finite or ideal vertices and finite volume. A truncated or sampled covering check does not establish the all-point hypothesis.

For future timelike old projections, the sublevel condition is exactly a hyperbolic ball; for null projections, a horoball. The explicit missed future-null ray is feasible for every old wall and strictly negative on every added cap. The vertical vector gives a nonzero horospherical recession direction. The direct strict spacelike witnesses avoid any reliance on a complete cusp census.

The mixed-bank checker verifies all 322 null directions and 720 doubled timelike directions, as well as the final √22 projection, with integer products and exact squared radical inequalities. The scale bound a>2 is explicitly justified. The π/7 checker reconstructs all 23 recipe normals using independently rounded Fraction intervals. It also checks that the literal 301 parent null projections are a subset of the 322-direction list and binds that bank by hash. Absolute height bounds cover either height sign of every π/7 cap, including the two reflected cusp-repair additions.

Thus the independent witnesses prove infinite volume for all three mixed 1,336-facet chambers and the 835/837-facet π/7 chambers. The finite-family obstruction is scoped to caps that remain strict at this missed cusp; a new null projection proportional to its ray can evade it. The earlier repository result already identified the general timelike-cap mechanism, and is credited in the report.

```sh
python3 scripts/fresh_h22_incidence_mixed_height_20261013_covering.py --verify
python3 scripts/fresh_h22_incidence_prism_20261013_uncapped.py --verify
```
