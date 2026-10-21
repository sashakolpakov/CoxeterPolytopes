# Independent arithmetic review of the Prokhorov refinement

**PASS.** The [separate checker](../scripts/dimension_bound_prokhorov_20261015_review.py) imports no producer code. It reconstructs graph distances by breadth-first search on160 connected spherical graph shapes through rank80 and evaluates the unchanged Prokhorov weights using exact fractions. Its [saved certificate](dimension_bound_prokhorov_20261015_review.json) passes generation and replay.

For rank at least17, paths have deficit `29N−3W=176` and D diagrams have deficit173. Independently of the finite small-graph census, a component of rank `r≤16` has deficit at least `r(61−3r)/2≥13r/2`. Consequently every spherical forest of total rank `N≥27` has `W≤(29N−173)/3`.

Substituting `C=29/3−173/[3(n−1)]` into Prokhorov's published `n<96C+68` gives the necessary inequality `n²−997n+6532<0`. For `n=991+t`, its left side is `t²+985t+586>0`. Thus the unchanged published geometric lemmas imply **n≤990**. Separate rational calculations reproduce the parity-sensitive incidence inequalities and confirm the endpoint distinction:990 remains allowed by this calculation.

This review checks arithmetic and the disconnected-graph argument. The [independent primary-source audit](dimension_bound_prokhorov_20261015_proof_review.md) checks the geometric imports, including ideal endpoints and nonsimple cusps. The [main report](dimension_bound_prokhorov_20261015.md) gives the complete derivation and the historical limitation: Nikulin previously announced the stronger general bound `n<300`. We make no new-best-bound or priority claim.

```sh
python scripts/dimension_bound_prokhorov_20261015_review.py --verify
```
