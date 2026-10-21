# Separation of the packing rule on every saved simplex star

The exact packing rule passes all 10,185 flexible stars that constrained its construction. Applied unchanged to the other 2,235 saved stars, it fails 1,417 fixed-discriminant cases. This is an explicit coverage failure of the generic rule outside its construction domain, not a failure of the valid simplex-packing theorem or of its partial feasibility certificate.

The [source](../scripts/fresh_sqrt5_rank13_packing_all_star_separation_20261017.py) and [result](fresh_sqrt5_rank13_packing_all_star_separation_20261017.json) bind the frozen [packing rule](fresh_sqrt5_rank13_simplex_packing_repair_20261017.md), [packing-charge theorem and export](fresh_simplex_port_incidence_20261017.md), and [complete saved selected-endpoint catalogue](fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017.md).

| Exact evaluation | Count |
| --- | ---: |
| Saved simplex stars | 12,420 |
| Flexible stars, all paid | 10,185 |
| Fixed-class stars paid | 818 |
| Fixed-class stars unpaid | 1,417 |
| Fixed classes with an unpaid saved star | 128 of 158 |
| Fixed classes whose saved stars all pass this rule | 30 |
| Saved stars using an absent forest | 198 |
| Additional forest types appearing in those stars | 97 |

Every old pair-orbit coefficient is evaluated exactly. Each packing charge is recomputed from the literal full port block and its complete diagram-automorphism orbit; it is not extrapolated from the old 45 exported nonzero coefficients. The original flexible-domain inequalities are independently recovered as passing cases.

The rule is explicitly extended to all 6,789 finite rank-13 forests by assigning zero pair payments outside its original 1,870 forests. All packing charges are supported on their original forests. Therefore each of the remaining 4,919 vertex forests has budget zero, and the complete vertex universe has maximum at most `239780393396/19999999999 < 12`. This completes the definition of vertex payments only; it supplies no missing face coverage.

## Small literal failure

Star 1237, in fixed discriminant class 27 of the saved catalogue, has a connected four-root tetrahedral core. Its central root is joined to three leaves by bonds 3, 3, and 5. Its orthogonal spectator forest is

`2 A2 + D4 + I2(5)`.

Deleting a core root gives the four finite vertex forests

- `2 A2 + D4 + H3 + I2(5)` twice;
- `3 A1 + 2 A2 + D4 + I2(5)`;
- `2 A2 + A3 + D4 + I2(5)`.

Every one of these forests already belongs to the optimized model. The rule pays exactly zero on all 12 triangle-corner occurrences of this tetrahedron and has no active packing contribution there. Its aggregate demand is 8. Thus this failure is independent of the zero extension to absent forests.

The core is an explicit arithmetic compact Coxeter tetrahedron locally. Write `phi=(1+sqrt(5))/2`, normalize all four roots to norm 2, and use central products `-1,-1,-phi`, with orthogonal leaves. Eliminating the three leaves gives the scalar Schur complement `(1-phi)/2`, negative physically and positive under conjugation. The core determinant is `4(1-phi)=2-2sqrt(5)`. Its one-root deletions are `3A1`, `A3`, and two copies of `H3`, all positive definite. The inverse Gram has negative entries, so its four dual rays are timelike on one sheet and span the compact tetrahedron. Its Gram entries are integral in `Z[phi]`.

The spectator determinant is `9 * 4 * (3-phi)`. Multiplying by the core determinant gives the saved full determinant `360-216sqrt(5)`, with the correct opposite signs at the two real places. Orthogonal finite spectators preserve the compact local tetrahedral face and arithmetic admissibility. This constructs local face data, not a compact ambient 13-dimensional polytope.

The JSON also retains the first zero-payment star by index, star 144 in class 18, and the worst payment and number of failures for every fixed class. No optimization was run in this package.

## Limits

Fixed-class failures must be interpreted using the common ambient form: a rule intended only for a class that excludes a given star need not pay that star. Conversely, the 30 classes passing this saved-star test do not acquire a dimension theorem. The saved catalogue still omits other connected cores and face families, including quadrilaterals.

Generation and a fresh saved reconstruction replay passed. The unchanged generic rule has a concrete unpaid fixed-class tetrahedron, while its exact flexible-domain success remains valid.
