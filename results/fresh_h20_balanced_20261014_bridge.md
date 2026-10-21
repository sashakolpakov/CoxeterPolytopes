# Balanced H20: central contacts and the remaining radial gap

An exact three-cap configuration now joins central-contact and vertex-direction caps in the unit-inradius `A6 × A6 × A7` cusp model. It is coherent with the complete radius-`sqrt(28)` vertex layer, preserves its 392 finite endpoints, and still has a strict spacelike escape. Separately, an all-label geometric bound shows that **every cap with strictly interior center and finite or parallel contact with the radius-two central cap lies in the Euclidean ball of radius five**. A cap reaching farther out while retaining that central facet must have its center on the alcove boundary or be dotted to the central facet.

The exact [checker](../scripts/fresh_h20_balanced_20261014_bridge_exact.py) and [certificate](fresh_h20_balanced_20261014_bridge_exact.json) replay with:

```sh
python3 scripts/fresh_h20_balanced_20261014_bridge_exact.py --verify
```

The arithmetic is `Fraction`. The shared [coordinate module](../scripts/fresh_h20_balanced_20261014.py) is hash-bound; no numerical candidate is trusted by this checker.

A separate [FLINT rational-matrix checker](../scripts/fresh_h20_balanced_bridge_review_20261014.py) and [saved audit](fresh_h20_balanced_bridge_review_20261014.json) independently replayed all 87,153 pairs, all 418 actual-facet projection witnesses, the 392 preserved finite endpoints, and the strict escape. Generation and saved replay both passed. Its total pair counts are 884 at π/3, 7,659 right, 397 parallel, and 78,213 dotted.

## Exact central-contact configuration

Write `vA_i`, `vB_j`, and `vC_k` for the three factor vertices from the shared module, and `z_ijk = vA_i + vB_j + vC_k`. Every product vertex has square `196`; the least nonzero product-vertex distance squared is `84`.

Start with the 22 unit seed normals, the central sphere of radius two, and all 392 vertex spheres of radius `sqrt(28)`. Their pair gates are exact: nearest A6 vertex pairs have angle π/3, distance-squared `112` gives parallel mirrors, and larger distances give dotted pairs. Each original vertex has 19 right seed contacts and its own spherical facet, hence an actual finite endpoint.

Add these three spheres:

| sphere | Euclidean center | radius |
|---|---|---|
| face bridge | `(3/11) [vA_0 + vB_0 + (vC_0+vC_2+vC_5)/3]` | `16/11` |
| first radial cap | `z_000/4` | `3/2` |
| second radial cap | `z_002/4` | `3/2` |

The bridge center has square `1092/121`. It meets the central sphere at π/3 and both radial caps at right angles. The two radial caps are mutually dotted and parallel to the central sphere. The bridge additionally has three parallel seed relations. No boundedness assertion is made about these new parallel-pair cusps.

The 1,248 new pairs are exactly 56 at π/3, two right, five parallel, and 1,185 strictly dotted. Thus the literal configuration has 418 actual facets. In coordinates with spatial Gram `B` and `D·U=-1`, the common strict center is `(0,200,1)`, of square `-400`. For any wall, orthogonally projecting this point to that wall remains timelike and strictly satisfies every other inequality, proving facet actuality.

All 392 prior finite endpoints `(z,(z²+28)/2,1)` remain feasible, with all nonincident inequalities strict. Their 19 seed roots form `A6+A6+A7`, and their own vertex cap is orthogonal to those roots; this gives a positive-definite rank-20 incidence Gram. This is not a complete vertex census.

The exact escape is

```text
x = (2/5) z_000,
E = (x, (x²-1)/2, 1),   E² = 1.
```

It is future-oriented and strict on all 418 inequalities. Its smallest spherical ball power is `54/25 > 1`, so the configuration has infinite volume.

The bridge and seed wall zero have Cartan two-cycle `841/64`, which is not an algebraic integer. Consequently any finite-volume Coxeter completion retaining these facets would be nonarithmetic, by Vinberg's integrality condition (V3), stated in [Bogachev–Yorov, Theorem 2.2](https://nvbogachev.github.io/Quasi_arithmetic_Coxeter_prisms-2024.pdf). This conditional observation does not assert existence, and does not exclude quasi-arithmeticity.

## All-label radius-five gate

Let `C` be the bounded unit alcove product. Its vertices have norm 14, hence `|b|≤14` throughout `C`. Let a new hemisphere have center `c` strictly inside `C`, radius `R>0`, and Coxeter-compatible relations to all 22 old seed facets. Put

```text
d = min_i (1 - v_i·c) > 0.
```

Each affine factor has unit offsets and marks one, so its distances sum to the number of its walls. Therefore `d≤1`. If `d=1`, then `c=0`, which cannot have acute contact with the central sphere. Otherwise set `g=1-d`. Every seed inequality gives `c/g ∈ C`, so `|c|≤14g`.

At a seed wall realizing `d`, the normalized cosine is `d/R`. A positive finite Coxeter cosine is at least `1/2`; parallel or dotted pairs have cosine at least one. Thus `R≤2d=2(1-g)` for **every finite label**, with no field assumption. Finite or parallel contact with the central radius-two sphere also gives `|c|≤R+2`. Hence

```text
|c|+R ≤ min(12g+2, 6-4g) ≤ 5,
```

where the two linear bounds meet at `g=1/4`.

For `x=(2/5)z_000`, `|x|=28/5`, so every ball in this entire family has power at least `(3/5)²=9/25` at `x`. The Lorentz lift of this same spatial point with square `1/4` is therefore a strict escape for any finite selection of such caps together with the central sphere and the uniform vertex layer. The argument includes all scaled proper-face centers with scale strictly below one; it excludes arbitrary boundary centers, for which zero seed distances are right contacts and the bound `R≤2d` cannot be used.

## Bounded numerical searches

These proposals are kept separate from the exact theorem and construction.

| scout | tested scope | outcome |
|---|---|---|
| [literal face centers](fresh_h20_balanced_20261014_bridge_scout.json) | 4,437 product-face representatives; central/seed labels 2–12 plus parallel | 1,542 representative candidates, 513,930 expanded caps; all inner probes at fractions `3/20` and `1/4` remain uncovered |
| [scaled face centers](fresh_h20_balanced_20261014_bridge_scaled_scout.json) | seed labels 3–6; central labels 2–6 or parallel | 99,279 accepted representatives; nine have proposed finite contacts to the 24 radial families; the exact motif above comes from one of these |
| [joint bank and boundary scout](fresh_h20_balanced_20261014_bridge_joint_scout.json) | 9,408 radial caps plus the root's 284,592 endpoint caps | all 392 probes at fraction `2/5` remain uncovered, with worst minimum power 2.16 |
| same joint scout, boundary part | 0–3 prescribed π/3 contacts per factor, other retained seed contacts right, one free dotted anchor per factor; central labels 2–16 or parallel | 1,806,336 profile/label trials; 2,639 accepted proposals; none cuts any of the 392 transported `2/5` targets |

The boundary search permits zero or negative endpoint Schur residual, so it genuinely goes beyond finite-endpoint profiles. It is nevertheless bounded in contact labels and support patterns. The 24 radial families have no numerical new–new finite contact at labels 2–60; that observation is not an all-label theorem. No joint compatible selection or full covering certificate is inferred from probe coverage.

Reproducible numerical producers are [literal faces](../scripts/fresh_h20_balanced_20261014_bridge_scout.py), [scaled faces](../scripts/fresh_h20_balanced_20261014_bridge_scaled_scout.py), and [joint/boundary testing](../scripts/fresh_h20_balanced_20261014_bridge_joint_scout.py). The literal-face expansion is stored compactly in [the NPZ](../data/enumerations/fresh_h20_balanced_20261014_bridge_scout.npz); its JSON retains orbit recipes rather than half a million duplicated rows.

The remaining constructive choices are boundary-center caps with changed finite-contact profiles, a changed central radius, or replacement of the central facet. Larger interior-center banks cannot cross the exact radius-five barrier while that radius-two facet is retained.
