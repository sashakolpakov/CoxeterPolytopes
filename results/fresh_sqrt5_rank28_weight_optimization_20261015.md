# Exact rank-28 obstruction for Boolean A1-sensitive weights

The optimum of the specified partial weight system is **78**, above the target 77. This is an exact obstruction to this weighting scheme, not a dimension or polytope nonexistence theorem. The [checker](../scripts/fresh_sqrt5_rank28_weight_optimization_20261015.py) and [certificate](fresh_sqrt5_rank28_weight_optimization_20261015.json) pass generation and saved replay.

The input consists of all 174 currently covered face graphs: 112 old sealed cores, 46 newly sealed singular stars, eleven further sealed attachment states, and the five earlier nonsealed stars in their zero-further-contact case. Both possible external-A1 states are imposed, giving 348 literal inequalities and 174 distinct inequalities. The reviewed attachment closure supplies the other cases of those five nonsealed stars.

The distance-two/H4 baseline remains fixed. There are 72 connected long marked-pair orbits appearing in these graphs and an independent nonnegative bonus for each orbit in each full-vertex A1 context: 144 variables. Every other pair bonus is zero. This retains the established triangle, short-path, and H4 payments. The optimization does not assume unrestricted monotonicity under enlarging a spherical component; its face coverage is inherited from the reviewed sealing arguments.

The optimizer only proposes rational coefficients. Exact Fraction arithmetic checks every primal inequality, every dual inequality, and equality of primal and dual objective. A complete rank-partition dynamic program over all 89 allowed irreducible finite types through rank 28 supplies separating forest constraints. Two LP iterations suffice; the final saved replay reconstructs the finite system and verifies the certificate without invoking the optimizer.

The returned rule has maximum 76 without an A1, attained at \(3E_8+H_4\), and maximum 78 with an A1, attained at \(A_1+D_6+3D_7\). All finite forests of rank 28 are included, not just forests observed in the face data.

There is also a short exact lower-bound certificate, independent of optimization. In the A1-present context:

* Classes 1018, 543 and 545 each force a single D6 bonus to be at least one. Their orbit multiplicities are 2, 2 and 1. The D6 baseline is 10, so its total is at least 15.
* Classes 535, 563, 843, 967, 1032 and 1052 force the six D7 bonus orbits, of total multiplicity 9. The D7 baseline is 12, so its total is at least 21.

Consequently the rank-28 forest \(A_1+D_6+3D_7\) has weight at least \(15+3\cdot21=78\). These nine forcing constraints already belong to the old 112-core input; adding the later graphs does not cause the obstruction.

The precise limitation suggests a next test. Class 1018 has forcing corner forest \(D_6+2A_1\), whereas the maximizing forest above contains only one A1. Classes 543 and 545 need only \(D_6+A_1\). The Boolean flag identifies these different contexts. A rule distinguishing zero, one, and at least two A1 components is therefore not excluded by this certificate.

The remaining finite and dotted quadrilateral fronts are still outside the covered face collection. Even an optimum at most 77 for a richer rule would need complete face coverage before yielding a dimension conclusion. This package introduces no new geometric realization, rank cutoff, or universal impossibility claim.
