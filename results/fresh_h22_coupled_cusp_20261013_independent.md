# Independent review of the 32 coupled-cusp bridge obstructions

**PASS.** The [independent checker](../scripts/fresh_h22_coupled_cusp_20261013_independent.py) imports no producer code or SymPy. It rebuilds the ambient rational form, its signature, the fixed null vector, the five finite chains, all 32 endpoint pairs, and the positive rank-17 projection space. Separate Fraction and quadratic-field arithmetic reproduce every coefficient in the [original certificate](fresh_h22_coupled_cusp_20261013.json). The [saved independent audit](fresh_h22_coupled_cusp_20261013_independent.json) also verifies the strict target inequalities against all 210 old walls.

For each endpoint pair, it confirms that every multiplier is positive, that

```text
z = X − sum(lambda_j*w_j)
```

is perpendicular to the specified null vector delta, and that the independently computed projection quantities satisfy

```text
rho > 0,    kappa < 0,    zeta >= 0,
kappa² − rho*zeta > 0.
```

The mathematical implication is sound. The span U of the seventeen prescribed old normals is positive definite and lies in `delta^perp`. The quotient `delta^perp / R*delta` is Euclidean. After removing U, four Euclidean directions remain. Hence every real norm-two candidate r with the prescribed pairings obeys

```text
(r,z) <= kappa + sqrt(rho*zeta) < 0.
```

The null component is unrestricted and contributes nothing to this bound. If r is compatible with all old walls, positivity of the multipliers gives

```text
(r,X) = (r,z) + sum(lambda_j*(r,w_j)) <= (r,z) < 0.
```

Thus the target remains strict for every such bridge, without any field, coordinate-size, null-shift or direction cutoff. Any finite collection of these bridges retains the strict spacelike escape and cannot complete a finite-volume chamber while retaining the old walls.

The scope is correctly restricted to the fixed `delta4`, orthogonality to old wall 26, and the stated pairings with endpoints of two different finite chains. The argument includes the proposed coherent affine-A20 completions assembled from these bridges. It does not exclude different cusp orbits, different incidence profiles, caps not incident to this null vector, or removal of old walls.

```sh
python scripts/fresh_h22_coupled_cusp_20261013_independent.py --verify
```
