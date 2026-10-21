# Independent review of the shared high-label star branch

The [separate reviewer](../scripts/fresh_h9_star_shared_high_20261014_review.py) reconstructs all **7,840 cases** of the [shared high-label family](fresh_h9_star_shared_high_20261014.md) using rational arithmetic in a fixed algebraic basis. It does not import the family's producer or symbolic algebra code. The [saved certificate](fresh_h9_star_shared_high_20261014_review.json) records a passing reconstruction.

For every case it builds the projected three-by-three Gram matrix, evaluates its determinant by permutations at parameter values 0, 4 and 5, and recovers the affine coefficients. The dependence is affine because the variable part has rank one. All slopes are positive; 7,392 cases are already strictly positive at 4. The remaining 448 cases fail positive definiteness at a literal vertex. Thus the argument covers the entire real interval `D>4`, including every integer last-block label at least seven.

The review also checks the complete Cartesian product of cases, all new vertex guards, projection-budget placements, the literal twin symmetry and the placement of earlier high-label gates. Those earlier gate proofs remain explicit, hash-bound dependencies. This excludes the specified contact family; the full eighteen-facet incidence remains open.

```sh
python3 scripts/fresh_h9_star_shared_high_20261014_review.py --verify
```
