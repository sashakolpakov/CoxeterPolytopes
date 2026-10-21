# An exact barrier to additive distance-two weights over Q(sqrt5)

The [rational LP certificate](fresh_sqrt5_sealed_weights_20261015.json) gives a concrete obstruction to the proposed **fixed distance-two baseline plus nonnegative long-pair corrections that depend only on the entire connected component and marked-pair orbit, with zero weight on disconnected marked pairs**. On the 112 sealed singular finite cores in the current necessary list, the best common component bound is \(2r+7\), not \(2r\). More sharply, six local faces force the total weight of a D7 component to be at least 21, and nine force an E8 component to be at least 23. Thus the rank-29 forest D7+D7+D7+E8 has forced weight at least 86 under this rule, exceeding the desired upper bound 84.

This excludes this particular local weighting strategy. It is **not** a dimension lower bound, a compactness claim for these full chambers, or an obstruction to using smaller baseline weights, payments on disconnected pairs, dependence on the full vertex forest, joint-face constraints, or genuinely global arithmetic information.

The [checker](../scripts/fresh_sqrt5_sealed_weights_20261015.py) takes the 112 singular finite cores with no accepted unpaid extra-white profile from the separate contact census. Its baseline is one on graph-distance-at-most-two pairs, plus the established marked-5-edge and endpoint corrections in an entire H4 component. The 94 variables assign nonnegative corrections to long marked-pair orbits in 16 finite component types. A quadrilateral must receive total weight at least one. A component budget counts each orbit with its exact number of pairs; no incidence multiplicity is omitted. Canonical rooted-tree encodings determine component types and marked-pair orbits exactly.

The optimizer only proposes coefficients. Every primal inequality, dual inequality and strong-duality equality is verified using `Fraction`. For the common-deficit optimization \(w(C)\le2\operatorname{rank}(C)-\delta\), the exact optimum is \(\delta=-7\). The upper certificate for this optimum uses just nine face inequalities and the E8 component inequality, all with multiplier minus one.

The separate component optimizations give the following sparse lower certificates:

| Component | Baseline | Forced additional weight | Forced total | Face-class indices in its dual |
|---|---:|---:|---:|---|
| D7 | 12 | 9 | 21 | 535, 563, 843, 967, 1032, 1052 |
| E8 | 14 | 9 | 23 | 940, 958, 982, 984, 993, 995, 1107, 1143, 1145 |

Each of these 15 faces has only one connected marked corner; all other marked pairs are disconnected. In E8 the nine faces require nine different long-pair weights to be at least one. In D7 the dual also accounts for the two-element pair orbits induced by interchanging its short arms. The exact LP data include these multiplicities and all dual multipliers.

These are actual local face geometries, rather than artifacts of a principal-inertia relaxation. The [separate geometry checker](../scripts/fresh_sqrt5_sealed_weights_20261015_geometry.py) and [certificate](fresh_sqrt5_sealed_weights_20261015_geometry.json) project the full white support of each witness and reconstruct its four black normals in a Lorentzian three-space. Across all 15 witnesses, the complete 90 pair-intersection calculations give exactly 60 feasible timelike rays—the four expected adjacent corners per witness—and reject all 30 opposite pairs. A strict timelike center is verified in each case, and all four vertices lie on its future sheet. The resulting two-faces are compact quadrilaterals. Their full principal Grams have one negative direction and positive-semidefinite conjugates; the finite contacts are the prescribed Coxeter labels.

A strict center of the local chamber follows by perturbing a face-interior point a sufficiently small distance in the negative white-dual direction. The local diagrams can also be augmented by orthogonal finite components to give the displayed rank-29 vertex forest. Neither operation proves that the resulting entire chamber is compact or extends to a compact polytope. Excluding such global extensions, or coupling different two-faces at one vertex, remains a possible route around this barrier. Likewise, the one-white census's sealing conclusion is conditional on the face remaining unpaid by an H4 correction; that inherited condition is not strengthened here.

Generation and saved replays pass. The contact census is an explicit input, and its entire arbitrary-extension proof is not independently repeated by this package. The next finite diagnostic should allow reductions of distance-two weights while retaining all justified triangle constraints. A successful relaxed LP would still require complete control of white extensions and the remaining nonsingular/dotted classes.

```sh
python3 scripts/fresh_sqrt5_sealed_weights_20261015.py --verify
python3 scripts/fresh_sqrt5_sealed_weights_20261015_geometry.py --verify
```
