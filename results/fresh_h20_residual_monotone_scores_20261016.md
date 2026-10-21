# Numerical residual scores redirect the search toward the largest open region

The [diagnostic](../scripts/fresh_h20_residual_monotone_scores_20261016.py)
and [saved data](fresh_h20_residual_monotone_scores_20261016.json) compare
the 151-wall base, its pi/8 extension and the 153-wall triangle extension
in the same affine C+four-caps coordinates. These are numerical search
scores, not exact cone or volume certificates.

| Chamber | Affine Euclidean volume | Largest vertex squared radius | Sampled exterior-direction fraction |
|---|---:|---:|---:|
| 151 | 1.597973574 | 16.33457530 | 0.962845 |
| 152, pi/8 branch | 1.597971859 | 16.33457530 | 0.962402 |
| 153, compact triangle | 1.597970887 | 16.33457530 | 0.962296 |

The polytope volumes come from floating Qhull convex hulls of the saved
numerical vertices. The direction fractions use 65,536 shared, seeded
uniform directions in the positive orthant of S3. The JSON includes
sampling standard errors, which do not cover numerical geometry errors.
It also records a radial Monte Carlo estimate of the Euclidean volume
outside the unit ball. This is not hyperbolic volume.

Adding a halfspace shrinks the region. Its exterior vertex count can
nevertheless increase, so that count alone is a poor progress criterion.
Both new caps reduce the numerical affine volume, but together remove
only about 0.000168 percent of the starting residual region and leave
the largest exterior radius unchanged. Their exact new finite vertices
and compact triangle remain valid local results; they do not close the
large exterior region.

The three largest retained exterior vertices have active sets
{133,141,144,149}, {133,143,144,149} and {141,143,144,149}. All have
z142 approximately 4.02259. These old133/149 boundaries are the next
target. The earlier theorem excluding straight-edge E4 cutters did not
exclude walls cutting this different region.

A decrease in these numerical scores is useful search guidance, not a
finite-volume proof. Even an apparently zero sampled exterior fraction
would not establish closure. An exact cone or other complete geometric
certificate is still required, and closing this residual face alone
would not establish finite volume of the whole H20 chamber.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_residual_monotone_scores_20261016.py --verify
```
