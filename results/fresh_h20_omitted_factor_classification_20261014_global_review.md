# Independent global-witness review: one omitted affine-factor facet

The global obstruction is sound within its stated rank-at-least-four, single-factor, one-omitted-facet scope. No issue was found in the actual-ridge argument, endpoint treatment, two exceptional parameters, or transfer of the 112 previously known directions. The result does not concern arbitrary H20 sections or classify the separate rank-three cases.

A [separate checker](../scripts/fresh_h20_omitted_factor_classification_20261014_global_review.py) uses only the Python standard library and imports no producer or project arithmetic. Its [saved certificate](fresh_h20_omitted_factor_classification_20261014_global_review.json) binds the frozen source and data. Both generation and saved replay pass. This is separate witness arithmetic plus mathematical source review; it does not independently enumerate the local directions or repeat the prior all-angle theorem.

| Separate reconstruction | Count |
|---|---:|
| Continuous new-plane families | 63 |
| Rational ambient face points, with all210 signs and exact zero sets | 126 |
| Exact interval-polynomial positivity certificates | 208 |
| Exceptional actual finite vertices | 2 |

The arithmetic uses rational coefficients in the squarefree radical basis `1,sqrt(2),sqrt(3),sqrt(6),sqrt(5),sqrt(10),sqrt(15),sqrt(30)`. Multiplication is reconstructed from gcds of squarefree radicands; signs use rational bounds from integer square roots. This is separate from the producer's bit-mask multiplication. It also reconstructs the source Gram matrix and strict timelike center directly.

## Why the angle witnesses are genuine ridges

The checked points x0 and x1 are future timelike, satisfy all source inequalities, and have exactly the advertised ambient face zero sets. Their normalization gives `(D,x0)=(D,x1)=-1`, `(e0,x0)=0`, and `(e0,x1)=-L`. Thus `(1-u)x0+u*x1` lies in the section with normal `e0-u*L*D` for every `0<u<1`, stays future timelike, and has precisely the common zero set K.

The old K Gram is positive definite. Its projected normals remain independent: the section normal has nonzero pairing `uL` with x0, whereas every old K root is orthogonal to x0. The projected roots lie in the positive-definite tangent space at the constructed timelike point. Their constraints consequently form a local orthant. Perturbing within that orthant while keeping either selected pair zero makes all other K inequalities strict, proving an actual ridge. A positive-definite pair Gram by itself would not have sufficed; the strict face-point and independence checks supply the missing geometry here.

The review recomputes every projected pairing and norm polynomial. It independently checks the Bernstein data by transforming each polynomial to the archived subinterval, reconstructing its Bernstein coefficients, and verifying their strict signs. Endpoint zeros are restored as explicit positive factors `u^a*(1-u)^b`. The subintervals cover [0,1] without gaps. This validates all 208 strict open-interval claims, rather than sampling parameter values.

At u=1, the only additional active source inequality is the omitted wall. The reconstructed ambient identity expresses it, on the section, as a nonnegative combination of K inequalities, with no coefficient outside K. It is redundant. The projected K independence argument still applies. The separate checker verifies the stated endpoint angle failures whenever this is an interior section. If every retained coefficient vanishes, the section is the old omitted facet hyperplane and does not meet the strict interior.

## The two exceptional parameters

For row110, the actual ridge has squared cosine `u²/(4-u²)`, so its only possible Coxeter value for `0<u<1` is 1/4. This leaves `u=2/sqrt(5)`. Row83 also has one endpoint requiring a separate witness.

For both exceptions the review checks the full normal, its unit square, the exact vertex, all210 inequalities, the exact twenty-wall zero set, and strict timelike/future orientation. It solves the rational old-Gram projection independently and verifies that the residual normal has positive square. The rank-one Schur criterion then makes the projected twenty-wall Gram positive definite. The designated pairs are actual ridges at a simple finite vertex and have squared cosine strictly between zero and 1/4. No facet-redundancy assumption is needed.

## Dependencies and scope

The 112 old-direction matches use equality of the full 22-coordinate normal planes modulo the cusp ray; matching only active-factor pairings would have required an additional premise. The producer now checks that full identity explicitly. Their global exclusion depends on the frozen [59-plane result](fresh_h20_cusp_sections_20261013.md), whose source manifest is bound.

One new family uses an angle interval between 1/2 and 3/4. That interval alone would contain the label-five value; its exclusion correctly invokes the separately proved [uniform crystallographic-angle theorem](fresh_h20_crystallographic_angle_gate_20261014.md) for hypothetical finite-volume Coxeter sections of this source. The other gap and obtuse witnesses need no arithmeticity assumption. The [independent review of that theorem](fresh_h20_crystallographic_angle_gate_20261014_review.md) remains a distinct dependency, not something re-proved by this checker.

```sh
python scripts/fresh_h20_omitted_factor_classification_20261014_global_review.py --verify
```

The [main report](fresh_h20_omitted_factor_classification_20261014.md) accurately separates local classification, global actual-ridge exclusion, the 48 irrational directions, and rank-three work. The present review supports that scope and makes no unrestricted dimension bound or nonexistence claim.
