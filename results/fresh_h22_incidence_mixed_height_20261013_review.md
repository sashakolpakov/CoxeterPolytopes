# Independent source review of the mixed-height H22 branch

**PASS.** Both frozen exact scripts replayed successfully. I reviewed their algebra, geometric conclusions and the [report](fresh_h22_incidence_mixed_height_20261013.md). This is source/proof review and fresh execution of the supplied checker, not a separate arithmetic implementation or an additional search.

The parity solve and exact orbit closures correctly recover322 null projections and720 norm-minus-one projections. The old-wall relations are checked directly: the scaled-null projections have only zero or sufficiently negative old pairings, and all720 timelike projections occur in the directly checked binary records. Thus the conclusions do not depend on assuming that a form-preserving generator necessarily permutes the old walls.

The changed scale `a=(sqrt(6)+sqrt(3))/2` repairs the cross relation exactly. For `b=-(delta,x)` in `{2,3,4}`, equal-height pairing is `-a*b+sqrt(6)`: at b=2 this equals `-sqrt(3)`, giving pi/6 between norm-two roots; at b>=3 it is strictly dotted. Opposite signs decrease the pairing. Distinct null projections are dotted and opposite caps over one null projection are parallel. The exact timelike-cap graph is three K240 components. Its omitted negative pairings have rational nonintegral `2*cos(theta)`, so cannot be finite Coxeter angles of any label. Choosing both height signs from any one component gives the three stated1334-facet chambers; no exhaustive classification of all possible signed selections is claimed.

The actual-facet and finite-vertex arguments are valid. All selected roots have positive norm, nonpositive pairings and the common strict timelike center `(h,0)`. Projection onto any wall gives a timelike relative-interior facet point. At each archived representative `(x,+/-1/sqrt(3))`, the old21-wall Gram is positive definite, the point has square `-2/3`, and every nonincident inequality is strict. The saved common spacelike vector is strict even on the incompatible union of all2294 inequalities, so every subset of that finite dictionary has infinite volume.

For the [two-cap extension](../scripts/fresh_h22_incidence_mixed_height_20261013_gap.py), `u²=-20` makes `(u,+/-sqrt(22))` norm two. All old relations and the uniform dotted bounds against both new banks are exact; the caps have mutual pairing `-42` and cut the preceding escape. Their points `(u,+/-20/sqrt(22))` have square `-20/11`. The21 old right contacts are positive definite, and all other inequalities are strict: their height magnitude is smaller than `sqrt(22)`, so the previously established uniform negative bounds also suffice at these points. These are two actual new finite vertices per branch.

The new rational escape is checked with Fraction arithmetic for the form, orientation and old pairings, and with rigorously outward square-root bounds for the new banks. It is strictly feasible on all2296 inequalities in the relaxed union and has square `1023438164699189/5000000000000000>0`. Hence all three1336-facet branches remain infinite-volume. No finite-volume or unrestricted impossibility claim follows from the dictionary check.

The cyclic-product observation is appropriately conditional: a rational noninteger length-two cyclic product prevents an arithmetic finite-volume reflection completion retaining that pair of facets. It is not used to claim finite volume. The separate integrality statement is restricted to the unscaled sqrt-two bank.

Fresh replays:

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_incidence_mixed_height_20261013.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_incidence_mixed_height_20261013_gap.py --verify
```

The matrix products use explicit int64 bounds. I also independently bounded the remaining direct bilinear operations using the frozen arrays: a common absolute bound is20,210,602,896, and a conservative bound for the pair squares is20,076,613,795,713,024; both are below2^63. No issue was found.
