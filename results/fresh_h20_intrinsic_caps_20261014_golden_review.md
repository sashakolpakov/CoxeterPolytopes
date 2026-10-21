# Scoped review of the golden 275- and 527-wall constructions

**PASS.** The sphere formulas, actual-facet and finite-endpoint arguments, and strict escapes in the [275-wall construction](fresh_h20_intrinsic_caps_20261014_golden_seed.md) and [527-wall extension](fresh_h20_intrinsic_caps_20261014_golden_extension.md) have no issue in this review. The [review certificate](fresh_h20_intrinsic_caps_20261014_golden_review.json) binds their frozen sources and exact data. This is source/proof review plus independent symbolic identities; it is not a second implementation or full replay of their 138,601 pair checks.

In the stated null coordinates, a sphere of center `a` and radius `R` has normal `(a,(a²−R²)/2,1)`, of square `R²`. Two such normals pair to `(R²+T²−|a−b|²)/2`. A point above `a` with height squared `h²` is `(a,(a²+h²)/2,1)`, of square `−h²`; its pairing with a sphere centered at `b` is `(T²−h²−|a−b|²)/2`. The [small separate checker](../scripts/fresh_h20_intrinsic_caps_20261014_golden_review.py) derives these identities symbolically without importing either producer. They justify the producers' distance-based pair checks and the preservation of the old height-one endpoints.

Each of the 252 old endpoints has nineteen independent spherical seed facets and its own orthogonal cap. Its physical rank-twenty Gram is positive definite, and every other inequality is strict. Thus it is an actual finite vertex. Every vertical seed facet occurs at one of these vertices, so a local orthant perturbation gives its relative interior. The central facet has its own strict timelike facet point. In the extension, all these endpoint and central-facet points survive, and every new facet has a timelike point with only its own inequality zero. This proves the stated 275 and 527 actual facet counts, not merely an inequality count. The extension's new facet points are not being counted as new vertices.

The saved common centers are strict and future timelike. Each saved escape is future spacelike and strict on every wall of its respective chamber. Such a strict spacelike direction, together with a strict timelike interior point, supplies an open set of ideal escape directions and hence infinite volume. The first new half-vertex cap cuts the preceding escape with pairing one, but the quarter-vertex escape remains strict on all 527 walls. Closing the original 252 vertical edges does not imply global volume closure.

The source gates check the original 37,675 pairs and all 100,926 added pairs exactly. The independent symbolic review checks their formulas, count accounting, source manifests, and the geometric implications; it does not repeat that scalar census or establish a complete corner census. The conditional arithmetic-field discussion in the seed artifact is outside this geometry review.

```sh
python scripts/fresh_h20_intrinsic_caps_20261014_golden_review.py --verify
```
