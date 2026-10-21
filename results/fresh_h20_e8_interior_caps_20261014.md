# Coupled interior E8 caps: 87 walls and four finite edges

**2026-10-14 — exact construction and independent replay PASS; volume remains infinite.** The [primary checker](../scripts/fresh_h20_e8_interior_caps_20261014_exact.py) and [exact data](fresh_h20_e8_interior_caps_20261014_exact.json) retain all [79 previous walls](fresh_h20_e8_cusp_completion_20261014.md) and add eight compatible norm-two facets. They cut the previous strict escape, create eight actual finite vertices and four finite edges, and preserve the completed source-face cusps. A different exact strict spacelike escape survives. The [independent arithmetic and geometry review](fresh_h20_e8_interior_caps_20261014_exact_review.md) imports no producer arithmetic and repeats all new pair and vertex tests.

The construction uses interior directions of the complete H12 face, with nonzero E8 components. It is outside the previously excluded orthogonal E8 additions and the fixed paired curves.

## Exact coupled pair and its four height images

Order the retained E8 roots as `[12,13,14,15,18,26,33,67]` and use zero-based weight positions. Let w1=ω3+ω6 and w2=ω4+ω5; their squares are 66 and 68. The two lifted source vertices V0,V1 have squares −1 and −13/121, with mutual pairing −1. Define

\[
r_i=-w_i+a_iV_0+b_iV_1,
\]

with

\[
\begin{aligned}
a_1&=17/8-\sqrt{869}/24,\\
b_1&=-2057/104-99\sqrt{11}/52+121\sqrt{869}/312+187\sqrt{79}/52,\\
a_2&=33/16-\sqrt{869}/48,\\
b_2&=-3993/208-99\sqrt{11}/104+121\sqrt{869}/624+363\sqrt{79}/104.
\end{aligned}
\]

All four coefficients are positive. Both roots have norm two and are mutually right. Their E8 contacts are π/3 at the specified two positions and right at the other six. Both are right to old walls 0 through 10; their pairings with old wall 19 are −1 and zero. Both cut the saved 79-wall escape.

Taking the four independent sign choices of ambient height coordinates 19 and 20 gives walls 79–82 from r1 and walls 83–86 from r2, ordered `++,+−,−+,−−`. The checker does not assume that the old 79-wall chamber shares these symmetries: it directly verifies every new/old and new/new relation. The exact field is Q(√2,√79,√11).

There are 660 new unordered pairs: 144 right, 20 at π/3, and 496 strictly dotted. A common strict timelike center and a strict point on every supporting wall certify all 87 actual facets.

## Finite geometry and the remaining escape

Every new cap has exactly 21 finite neighbors. Testing all 210 choices of nineteen companions for each cap gives 1,680 exact tests. Positive-definite Gram matrices are followed by physical vertex reconstruction and all 87 inequalities. The result is a complete finite-vertex census incident to the new caps: four A11+A7+A2 vertices and four A11+E7+A2 vertices. It is not a complete ambient vertex census.

The four corresponding pairs share positive-definite rank-19 supports, giving four actual finite edges. At the unflipped pair, the two supports differ by exchanging old wall 19 and new wall 79. The E8 profile `{4,5}` already supplies its finite corner individually; the other cap closes the adjacent finite edge. All 35 previously certified H12-face finite vertices and both completed cusp links survive, since every new wall is strict there.

The exact future spacelike escape saved in the certificate has norm

`1336028658157091827 / 250000000000000000`.

It is strict at all 87 walls, proving infinite volume. This is finite-edge progress, not a finite-volume example or a complete corner/cusp census. The retained irrational Cartan cycle `12+8√2` and rational indefinite two-wall block still give the previous conditional non-quasi-arithmetic obstruction to any cofinite completion retaining those facets.

## Bounded numerical discovery and next structural route

The [source-edge scout](fresh_h20_e8_interior_caps_20261014_scout.json) tested the 36 one/two-node π/3 E8 profiles, omitting the single null-residual case. It used positive combinations on 234 edges of the 37-vertex H12 source, with one additional old-wall anchor at labels 2,3,4,6. Of 212,454 quadratic branches, 105,723 were positive-cone branches, 9,377 cut the target and 818 were acute. Twenty-two passed the numerical Coxeter screen; their only compatible pair supplied these recipes. These are bounded numerical counts, not an exhaustive all-label or all-direction theorem.

The root's [87-wall volume scout](fresh_h20_e8_orbit_caps_20261014_volume_scout.json) used 288 LP proposals and found 101 strict numerical escapes. Its strongest proposal was rounded and strictified to obtain the exact escape above. Acceptance and volume conclusions depend on the exact checks, not those numerical predicates.

The new E7 corner yields a more useful next problem than another source-edge bank. The [exact H13 face certificate](fresh_h20_e7_face_20261014.md) identifies its current A11+A2 cone, rational intrinsic form `diag(A11,A2,−2)`, and the unique allowed nonright E7 meeting profile. A complete intrinsic source and its ambient-compatible lifts remain to be found.

```sh
python3 scripts/fresh_h20_e8_interior_caps_20261014_exact.py --verify
python3 scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py --verify
python3 scripts/fresh_h20_e7_face_20261014.py --verify
```
