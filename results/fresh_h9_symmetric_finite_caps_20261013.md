# Complete finite golden leaf-angle census on the symmetric H9 seed

No additional cap was found. Every square-two real normal whose sixteen leaf pairings lie in `{0,-1,-phi}` and which is Coxeter-compatible with the two center walls belongs to the already known orbit of 70 balanced caps. This statement allows arbitrary real center pairings, including dotted ones. A separate, broader numerical scout also found no new orbit, but its result is not an exact exclusion.

The [18-wall seed and 70 caps](../data/enumerations/fresh_h9_symmetric_balanced_caps_20261013.json) are described in the [symmetric-filling report](fresh_h9_symmetric_filling_20261013.md). Set `b = (1+2phi)/2` and `p = 1+phi`. Exact row identities give

```
u_i + v_i = u_j + v_j,
c_0 + c_1 = (p/b)(u_i + v_i).
```

The ten normals `c_0,c_1,u_0,...,u_7` form a basis. The archived inverse Gram matrix certifies this, while the eight displayed independent relations account for the seed's full nullspace.

Write the nonnegative leaf-pairing magnitudes as `l_i,r_i`. Their sums must be constant. The six sums of two elements from `{0,1,phi}` are distinct, so every pair uses one fixed unordered pair of magnitudes. If those magnitudes differ, permutation symmetry reduces the assignments to the nine possible counts of the first orientation. For equal magnitudes there is one assignment. This gives exactly 30 representative orbits.

For each orbit, the sum of the two center pairings is fixed at `-p(l+r)/b`. The first center pairing is the sole remaining coordinate. Substitution into the inverse Gram matrix gives one exact quadratic equation imposing normal square two. Both real roots are checked, including the repeated-root case. No center-angle alphabet or bound on dotted distance is imposed.

The exact result is:

| Outcome | Real orbit solutions |
|---|---:|
| At least one positive center pairing | 52 |
| A center pairing strictly between -1 and 0 | 6 |
| Compatible balanced orbit | 1 |
| Total | 59 |

For square-two normals, a Coxeter pairing is zero, at most -1 at a nonright finite angle, or at most -2 for parallel or disjoint walls. Thus the first two rows reject every possible center label. The final orbit has four leaf pairs `(0,-phi)`, four `(-phi,0)`, and center pairings `(-1,-1)`. Its 70 placements coincide exactly with the old bank. It supplies no new mirror cutting the saved 20-wall escape.

The [certificate producer](../scripts/fresh_h9_symmetric_finite_caps_20261013.py) uses the seed verifier's Fraction arithmetic. Its [exact archive](../data/enumerations/fresh_h9_symmetric_finite_caps_20261013.json) contains the inverse basis Gram matrix, every quadratic and discriminant, both center-pairing expressions, and exact sign decisions. A [separate replay](../scripts/fresh_h9_symmetric_finite_caps_verify_20261013.py) uses the repository's independent Golden arithmetic implementation, reconstructs the seed and all orbit equations, verifies the inverse certificate, and matches all 70 accepted placements. The [audit passes](fresh_h9_symmetric_finite_caps_verify_20261013.json) and rejects a corrupted norm equation. A second agent source-reviewed the completeness argument and independently reran the replay successfully.

The [broader numerical scout](../scripts/fresh_h9_symmetric_finite_caps_scout_20261013.py) uses finite labels `2,3,4,5,6,8,10,12,18` for leaves and finite center angles, also allowing arbitrary dotted or parallel center pairings. It checks 333 representative assignment orbits. Twenty-one norm equations have no real solution; among the remaining roots, 600 have a positive center pairing and 22 fail the declared finite-center angle gate. The only survivor is the already certified balanced orbit. The [numerical archive](../data/enumerations/fresh_h9_symmetric_finite_caps_scout_20261013.json) records its tolerances through the source, pair-sum separation, and all surviving proposals. This is a numerical scout, not a proof that no additional cap exists for those labels. Its representative pairings with the selected caps also do not constitute an exhaustive placement search.

The exact theorem is limited to leaf pairings `{0,-1,-phi}`. It does not exclude dotted leaf pairs, other finite labels, deleting old facets, or different seeds. The broader numerical scout does not enlarge the exact theorem's scope.

```sh
PYTHONPATH=scripts python scripts/fresh_h9_symmetric_finite_caps_20261013.py --verify
PYTHONPATH=src python scripts/fresh_h9_symmetric_finite_caps_verify_20261013.py
python scripts/fresh_h9_symmetric_finite_caps_scout_20261013.py
```
