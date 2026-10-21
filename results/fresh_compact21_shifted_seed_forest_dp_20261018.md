# A complete face-feasible seed, but no improved rank-21 budget

The completed rank-22 H3-cross rule has an exact retained primal. A guaranteed
lift to rank 21 uses its bonuses at A1 index a+1, where a is the actual entire
A1 count. The complete exact finite-forest test of this lifted seed gives

    min over 0<=t<=1 of the new forest maximum = 21977/500 = 43.954.

An unchanged attaining forest is A1+A4+2A6+2I2(5). Its budget is independent of
the new H4 parameter. Thus this fully face-feasible seed does not reach the
target 40, and it cannot simply replace the sparse below-40 support candidate.
This does not refute the funded rule or rule out a better full-model vector.

## Why the index shift is necessary for the inherited guarantee

The source is
`fresh_sqrt5_dimension22_h3_cross_20261016_cap2.json`, whose certificate stores
30,169 bonuses and its exact upper coordinate. Its independent completed proof
is inherited here. Setting the newly introduced A4-specific cross bonuses to
zero and replicating across A4 contexts preserves the old three-hub rule.

A same-index restriction does not automatically inherit every face inequality:
a rank-21 case with zero untracked residual would have residual one in rank 22,
which was excluded from the old common-spectator list. We therefore add one
whole orthogonal A1 spectator, using source index a+1 everywhere. This maps every
old rank-21 native case to an admitted rank-22 case with the same untracked
residual. Connected and I2(5)/H3-cross payments agree. At an actual A1-cross
angle, a>=1 and its denominator changes from a+1 to a, so its nonnegative
payment can only increase. When a=0 there is no such marked A1 angle. All
shifted variable indices lie inside the rank-22 domain. This establishes face
feasibility conditional on the previously completed source coverage.

## Exact complete finite-forest calculation

The script streams the 30,170 scalar primal entries and the 30,169 variable
records. It builds exact marked-pair and rooted profiles for all 68 allowed
finite component types through rank 21. Entire A1, I2(5) and H3 counts are fixed;
their capped contexts are then used in the coefficient lookup. The remaining
rank is filled by the full list of nonhub finite components. All 362 admissible
exact hub contexts are checked.

Each rank/context retains the three best **distinct** old forests. Only two
physical forest types have decreased budgets: H4+A1+2E8 and D5+2E8. Their
active/inactive states do not create another physical forest. Three retained
forests therefore suffice to recover a best unchanged one after excluding
both decreases. All four positively adjusted source forests are evaluated
explicitly even if absent from these top-three lists.

The D5 uniform saving is min(1,2 beta), with beta the extreme E8 endpoint-pair
bonus of this lifted rule. The H4 transfer gives seven affine upper lines:
the unchanged maximum, H4 with slope -1, D5 with slope zero, and the four
source forests with slopes 1,2,1,1. Checking endpoints and their pairwise
intersections in [0,1] gives the exact best envelope. No optimizer is imported.

Selected old forest totals are 43.9527 at H4+A1+2E8, 42.5308 at A4+A1+2E8,
42.383 at A5+2E8, 41.1623 at B4+A1+2E8 and 43.8831 at B5+2E8. The exact
certificate retains all seven lines and compact attaining frontiers.

## Reproducibility and resource record

```sh
python3 -B scripts/fresh_compact21_shifted_seed_forest_dp_20261018.py --verify
```

Generation passed in chunk `ee4588`, and fresh separate-process replay passed
in `7f557a`. The guard recorded 2.481 seconds and 67,862,528 bytes peak combined
RSS, with 2,901 bytes worker output. Worker peak was 36,159,488 bytes. Initial
free disk was 198.34 GiB and available memory was 7.30 GiB. Both runs held the
global lock, enforced 59 seconds and 1 GiB RSS, exited and removed temporary
output. No numerical package, Lean process or full constraint matrix was used.

Verification status: exact internal finite DP plus written lifting argument;
completed rank-22 mathematical coverage is an inherited input. No new dimension
bound, new full-model primal at 40, external review or Lean verification is
claimed.
