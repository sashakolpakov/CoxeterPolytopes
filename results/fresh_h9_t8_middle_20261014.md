# Allowing every middle-block label in the T8 contact family

The [complete own-contact theorem](fresh_h9_t8_asymmetric_20261014.md) extends when the middle internal label on `(7,8)` is allowed to be **any finite integer at least four**. Together with the first internal label `m >= 7`, the last internal label `n >= 4`, and all four anchor-own labels free, this is a family exclusion, not an exclusion of the whole case-8 incidence.

The [checker](../scripts/fresh_h9_t8_middle_20261014.py) and [certificate](fresh_h9_t8_middle_20261014.json) preserve all earlier files. They bind the preceding theorem and the [next necessary assignment](fresh_h9_t8_asymmetric_20261014_search_case_08.json), which first introduced label five on `(7,8)`.

The [independent review](fresh_h9_t8_middle_20261014_review.md) reconstructs the changed block inverse and all 256 new polynomial cases with separate Fraction arithmetic, matches every coefficient, verifies the inherited 152-case complement and checks the domain argument. Its generation and saved replays pass; no producer arithmetic or SymPy is imported.

## Domain completeness

All literal actual-vertex witnesses in the preceding theorem remain valid when this middle label increases from four: increasing a nonright cosine decreases the corresponding Gram entry. A nonpositive relaxed Gram witness therefore remains nonpositive. The inverse of the isolated two-node block has entrywise nonnegative entries increasing with its cosine. Hence all Schur-budget lower bounds persist as well.

In particular, the four blocks of the positive eight-wall support remain orthogonal; both anchors still have at most one foreign contact, of label three; any foreign contact forces every own nonright contact of that anchor to have label three; wall 14's support projection is unchanged; and its two anchor contacts are still restricted to labels two or three. Each of the nine own-support patterns inherits its own literal witnesses. No new incidence symmetry is assumed.

Let `p` denote the middle label on `(7,8)`. The complete split is:

- `p=4`: the preceding 408-case certificate applies literally.
- `p=5`: only cases with at least one foreign contact at 7 or 8 change. There are 256 such cases; the other 152 are unchanged.
- `p>=6`: a label-three contact into that block has projection cost `1/[4 sin(pi/p)^2] >= 1`, before adding the positive own-contact cost. The actual vertex must be positive definite, so both anchors are right to 7 and 8. All remaining polynomial cases are among the 152 inherited cases.

This split includes every finite integer label, including the unbounded tail.

## The changed coefficient bank

For `p=5`, the self cost of a foreign label-three contact is

```text
(5+sqrt(5))/10,
```

and the overlap for contacts to opposite roots of that block is

```text
(5+3sqrt(5))/20.
```

These follow directly from the inverse of `[[1,-cos(pi/5)],[-cos(pi/5),1]]`. All own projections, the third-wall residual pairings, the six continuous parameter boxes and the lower value `Q=25/19+sqrt(2)/2` are unchanged. Substituting the two changed costs into the same general determinant identities gives 256 exact cases over `Q(sqrt(2),sqrt(5))`.

The certificate stores their power and Bernstein coefficients in the basis `1,sqrt(2),sqrt(5),sqrt(10)`. Positivity uses exact nested quadratic-field comparisons. All 4,320 coefficients are strictly positive; their minimum is `7/50-31sqrt(5)/2000`. Thus both the residual determinant numerator and the negative anchor determinant numerator are positive, and the same monotonic extension in the first internal label and the rank-eight Schur contradiction apply. The arithmetic does not substitute values for dotted distances.

## Replay and scope

```sh
python scripts/fresh_h9_t8_middle_20261014.py --verify
python scripts/fresh_h9_t8_middle_20261014_search.py --case 8 --limit 100000 --verify
```

Generation and saved replay pass for both commands. The separate [continuation](../scripts/fresh_h9_t8_middle_20261014_search.py) stops at its next necessary assignment after 1,450 nodes: 12,491 spherical-vertex failures, 972 forced-orthogonal negative-block failures, and 27 completed mandatory assignments excluded by the whole new family. Its 100,000-node limit is not reached. It is a first-survivor search, not a completed search of this incidence.

The [preserved next state](fresh_h9_t8_middle_20261014_search_case_08.json) changes the own contacts of leaf 6 with `(7,8)` from `(3,3)` to `(2,3)`, and those of leaf 15 likewise to `(2,3)`. The internal label `(7,8)` now denotes every finite integer at least seven. Anchors 11 and 17 both have own labels `(3,3)`, and the last block `(12,13)` has label four. The first block remains the unbounded `m>=7` category. This state is outside the proved fixed-label scope and is not a metric realization.

A useful next task is to isolate which literal vertex witnesses actually need the old leaf-6 and leaf-15 contacts. The selected residual formulas use neither leaf, but the stated domain theorem binds all other mandatory labels. Removing that restriction requires a fresh source-bound guard argument, not merely observing that the formulas are unchanged. That work is left open here; no further parameter extension or search is included.

Only the explicitly listed mandatory family is excluded. Other mandatory contacts remain fixed. No compact H9 example, and no exclusion of the whole incidence, is claimed here.
