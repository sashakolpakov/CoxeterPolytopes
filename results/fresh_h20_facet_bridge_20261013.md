# H20 facet bridges and a polar-reflected replacement

No finite-volume example was found. This round changes the construction from the previous dotted symmetry orbit: it searches for walls that meet the last facet, and for a coherent polar-reflected replacement of that facet. Two complete, narrowly scoped obstructions are now exact. The current geometric checkpoint remains the [176-wall chamber](fresh_h20_coupled_completion_20261013.md), which has a strict spacelike escape.

The [exact certificate](../data/enumerations/fresh_h20_facet_bridge_20261013_exact.json) is rebuilt by the [verifier](../scripts/fresh_h20_facet_bridge_20261013_exact.py). Its deterministic `--verify` replay passed. A separate agent also completed a [source review and fresh replay](fresh_h20_facet_bridge_20261013_review.md), finding no issue. It uses the frozen 176-wall data and the prior separate rational/radical arithmetic. Its interval calculations use integer outward rounding at scale 10^60; the trigonometric enclosures use the previously checked rational Machin identity and cosine Taylor remainder. This is a new proof using existing arithmetic helpers, not an independent implementation of the 176-wall checkpoint.

## Complete perpendicular-bridge obstruction

Facet 175 has exactly 18 finite neighbors, indexed

```
0,2,3,4,5,6,8,9,11,12,13,14,15,16,17,20,22,40.
```

Consider every positive-norm proposed wall perpendicular both to facet 175 and to those 18 normals. The verifier first computes the exact rank 18 kernel over Q(√6), imposes the additional facet 175 equation, and certifies a positive-definite two-dimensional Gram matrix. Thus every proposed wall can be scaled to a unit vector on this plane's unit circle.

The normalized projection of old wall 18 onto the plane has square approximately 0.8813428425610371, rigorously below cos²(π/9) and above cos²(π/8). Consequently its relation to a unit candidate cannot be dotted or parallel, and any compatible finite label must be one of 2,3,4,5,6,7,8. Each label gives exactly two directions. All 14 directions are enclosed rigorously, and each has a positive pairing or a cosine strictly between adjacent permitted Coxeter cosines with one of old walls 10,30,50.

This exhausts all real directions and all integer finite labels for this particular perpendicular support. It does **not** exclude nonright relations to facet 175, changes to one or more of its 18 neighbor relations, or replacements that make old inequalities redundant. The assertion is about complete pair compatibility with the retained normals; it does not infer actual forbidden ridges merely from positive-definite pair Grams.

## Complete obstruction to the polar-reflected bank

Let P be the positive-norm hyperideal point saved with the 175-wall predecessor. For each predecessor normal r, take its reflection

```
s_P(r) = r − 2(P,r)P/P².
```

This is a coherent source of new walls: reflection preserves all mutual pairings inside the reflected bank. A numerical cross-pair scout found 152 compatible nontrivial reflected profiles and 20 fixed profiles, with three rejected reflected profiles. That 327-wall numerical selection was only a proposal; no exact327-wall Coxeter claim is needed for the following stronger rejection.

The exact verifier checks `(P,r)≤0` for all 175 predecessor normals. It also rechecks the current176-wall strict future spacelike escape X, including every one of its 176 strict inequalities, and proves `(P,X)<0`. Hence

```
(s_P(r),X) = (r,X) − 2(P,r)(P,X)/P² < 0
```

for every predecessor normal r. The same strict escape therefore survives the **entire 175-entry reflected bank**, together with all 176 existing inequalities. Every subset of this combined dictionary is infinite-volume, including a replacement that deletes facet 175 or any other current walls. This argument does not depend on the numerical compatibility screen. It concerns this one polar reflection, not arbitrary chamber assemblies or successive reflections about different hyperplanes.

## Bounded numerical construction attempts

These scouts are proposals and rejections in floating-point arithmetic; their counts do not constitute exact classification results.

| Family | Scope and outcome |
|---|---|
| [One released neighbor, right to175](../data/enumerations/fresh_h20_facet_bridge_20261013_release_scout.json) |3456 angle profiles; no fully compatible candidate.|
| [One released neighbor, nonright/parallel to175](../data/enumerations/fresh_h20_facet_bridge_20261013_nonright_scout.json) |13824 profiles;255 compatible numerical candidates, all missing both chosen exposed facet ends.|
| [Continuous remaining dotted anchors](../data/enumerations/fresh_h20_facet_bridge_20261013_continuous_scout.json) |1458 prescribed-label profiles;177 feasible angular-arc families, none cutting either exposed end. This includes continuous parameter intervals rather than only sampled points.|
| [Five-direction replacement](../data/enumerations/fresh_h20_facet_bridge_20261013_replacement_scout.json) |3125 coefficient profiles;2912 distinct rounded numerical normals passed the old 175 pair gate. At most 19 old hyperplanes meet any one of them; no volume conclusion follows from this count.|
| [Coherent polar reflection](../data/enumerations/fresh_h20_facet_bridge_20261013_reflection_scout.json) |152 accepted nontrivial reflected profiles; the exact universal escape argument above rejects the entire larger reflected bank.|

The corresponding scripts have the same basenames under `scripts/`. The [initial plane diagnosis](../data/enumerations/fresh_h20_facet_bridge_20261013_diagnosis.json) and [14-branch scout](../data/enumerations/fresh_h20_facet_bridge_20261013_right_scout.json) are superseded by the exact perpendicular-bridge proof.

The remaining constructive problem requires changing more of the neighbor relations or using a different polar hyperplane/assembly. The proved obstructions do not establish nonexistence of finite-volume Coxeter polytopes in H20, nor noncompletability of the 176-wall chamber by unrestricted compatible additions.

Replay from the repository root:

```sh
python scripts/fresh_h20_facet_bridge_20261013_exact.py --verify
```
