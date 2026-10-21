# Independent review of the enlarged H5 bank obstruction

The four-witness argument is sound for the frozen dictionary of **38,453 entries**: 125 optional old walls, 38,209 original formula branches, and 119 opposite branches. No pairwise Coxeter-compatible selection from this dictionary has finite volume. The result allows deleting optional old walls. It does not exclude other normals, other seeds, or arbitrary completions.

Reviewed sources are the [enlarged verifier](../scripts/fresh_h9_global_h5_20261013_enlarged_verify.py), [recipes](../data/enumerations/fresh_h9_global_h5_20261013_enlarged_recipes.json), and [exact certificate](../data/enumerations/fresh_h9_global_h5_20261013_enlarged_verify.json), together with the [original-bank review](fresh_h9_global_h5_20261013_review.md). A separate [standard-library Fraction checker](../scripts/fresh_h9_global_h5_20261013_enlarged_review.py) records its own [audit artifact and hashes](fresh_h9_global_h5_20261013_enlarged_review.json).

The logical cover is exhaustive within this dictionary. If a selection uses none of the 119 opposite branches, the original common strict escape survives. Otherwise select any such branch. Its assigned witness has a complete cutter list, and every member of that list is incompatible with the selected branch. All remaining entries pair strictly negatively with the witness. The four lists contain 15, 54, 1, and 2 cutters, respectively; they include possible cutters among all 125 old walls and all candidate entries. Optional deletions therefore do not evade the argument.

The primary replay checks all **153,812 candidate–witness signs** and **1,198 actual cross-pairings**. The latter consist of 449 positive pairings and 749 pairings strictly between consecutive finite-angle values. Such an open interval excludes every integer Coxeter label by monotonicity of \(-2\cos(\pi/m)\). Its endpoints are strictly above \(-2\), so parallel or dotted pairs are not discarded accidentally. The largest gap used is between labels 28 and 29; this is an observed certificate value, not a label cutoff. The projection formula used to compute candidate–candidate pairings is the exact orthogonal decomposition against one candidate's nine-contact span.

All four witnesses have positive ambient norm and the correct future orientation. Their H5 projections also have positive norm and the correct orientation. Keeping the four named walls as the H5 face makes every permitted selected normal nonpositive on those walls, so orthogonal projection can only decrease its pairing with a witness. The H9 conclusion itself does not require those four walls to remain selected.

The separate checker independently recomputes the four exact ambient norms and orientations, their projected sign bounds, the coverage logic for all 119 branches, and every conflict inequality. It recomputes the cosine enclosures with its own rational rounding and Taylor routine. It does **not** independently recompute all 153,812 root signs or the 1,198 root pairings; those are covered by mathematical source review and a fresh full primary-verifier replay. Dependencies on previously certified normal enclosures are explicit and hash-bound.

Validation: primary saved replay **PASS**; separate Fraction checker generation and saved replay **PASS**. No mathematical issue found in this scoped audit.

```sh
PYTHONPATH=src:. python scripts/fresh_h9_global_h5_20261013_enlarged_verify.py --verify
python scripts/fresh_h9_global_h5_20261013_enlarged_review.py --verify
```
