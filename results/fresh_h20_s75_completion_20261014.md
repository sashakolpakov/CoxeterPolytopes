# A continuous escape gate for the untouched S75 face

The [exact certificate](fresh_h20_s75_completion_20261014_gate.json) proves that **every acute Coxeter cap inducing a boundary geodesic on the retained S75 hyperbolic plane strictly misses the current 98-wall escape**. It covers arbitrary real normals, lengths and labels, without a coefficient lattice or a height cutoff. Consequently no finite extension using only such caps can complete the current chamber to finite volume.

This is a restriction on one route. It does not exclude a mixed extension that also adds walls whose S75 projections are causal, or a construction that changes retained walls. No new wall is added in this package.

## Why this face was tested

The [98-wall construction](fresh_h20_bridge_edge_closure_20261014.md) closes one edge on S61 while preserving the entire alternative S75 plane. The [independent geometric review](fresh_h20_bridge_edge_closure_20261014_review.md) establishes the latter assertion exactly. Its support is

```
S75 = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,17,31,75,87],
```

with norm-two Gram of type `A10+D7+A1`.

The natural constructive test allowed nonright contacts with low inverse-diagonal support nodes, while retaining the necessary right contact to wall 75. The separate [complete profile gate](fresh_h20_s75_projection_gate_20261014.md) reduces this to 17 explicit profiles and the single-node-17 family with arbitrary finite label. This includes possible noncrystallographic labels; there is no restriction inherited from an arithmetic root algorithm.

The [bounded numerical scout](fresh_h20_s75_completion_20261014_scout.json) sampled 12 labels on that one-parameter family, giving 29 profile cases. It used nine specified finite anchors and labels `2,3,4,5,6,8,10,12`, plus parallel contact. Among 84,636 normalized proposals, 25,184 cut the current escape, but none passed all acuteness gates. These counts are numerical diagnostics only. Physical Gram entries, the residual basis, and projected wall/escape pairings were evaluated at 120 decimal digits before conversion to floating point; exact support zeros and known root norms were imposed. The later eigensolver, linear programs and quadratic gates remain numerical and do not prove nonexistence.

The exact result below is stronger and independent of that failed search. The scout only suggested a small basis for the coefficient identity.

## Mandatory right contacts

Let `J` be the inverse of the positive-definite support Gram. It is entrywise nonnegative. For a proposed new normal `r` with `r²=q>0`, put `b_i=r·N_i≤0`. Its orthogonal projection to the face has square

`q − bᵀJb`.

A boundary geodesic requires this square to be positive. At any nonright Coxeter contact with a norm-two support wall, including parallel or dotted contacts,

`b_i² ≥ q/2`.

Since all products `b_i b_j` are nonnegative, a diagonal entry `J_ii≥2` would force `bᵀJb≥q`. Therefore every cap inducing a geodesic must be right to

```
F = [3,4,5,6,7,8,13,14,15,75].
```

The checker verifies the rational inverse and these inequalities exactly. This argument alone already handles every label; enumeration of the remaining profiles is not needed for the escape proof.

## One exact identity covers all remaining profiles

Let `E` be the rational strict escape retained by the frozen 98-wall chamber. The certificate gives the exact physical-vector identity

`E = Σ_{i∈F} μ_i N_i + Σ_{j∈P} λ_j N_j`,

where

```
P = [0,2,18,31,65,79,84,85,88,89,92],
```

every `λ_j` is strictly positive, and the 21 normals indexed by `F∪P` are linearly independent. The ten coefficients `μ_i` are unrestricted. All coefficients lie in the 16-element multiquadratic field already used by the predecessor; neither the later bridge square root nor the large new boost radicals are needed.

The producer solves the identity using the field's 336-dimensional rational regular representation, then independently recomputes every physical coordinate using Fraction radical multiplication. It checks exact matrix rank and the signs of all eleven positive coefficients. A numerical residual is not used as a certificate.

For any nonzero normal `r` right to `F` and acute to the eleven retained walls in `P`,

`r·E = Σ_{j∈P} λ_j(r·N_j) ≤ 0`.

Equality would make all eleven pairings zero. Together with the ten right contacts and full rank, that would imply `r=0`. Hence **`r·E<0`**. The argument imposes no norm equation, field restriction or discrete label assumption after the mandatory-right reduction.

The saved signs give a sharper geometric interpretation. Only the coefficients at walls **8 and 15** are negative; all other 19 coefficients are strictly positive. Thus right contacts to those two walls already force strict failure to cut `E`. Every acute escape-cutting Coxeter cap must have a nonright contact to at least one of them. Their inverse diagonals are `24/11` and `4`, so such a cap's S75 projection has square at most `−q/11<0`. In any coherent extension it is future timelike and therefore leaves the entire S75 plane strict. Cutting this escape and adding a boundary to S75 require different caps under the retained-wall conditions.

The current witness has exact square

`E² = 20181174579467302741 / 10^18 > 0`.

It is already strict against every old wall. Thus any finite extension made only from S75 boundary-producing Coxeter caps retains this same strict spacelike escape and remains infinite-volume. No perturbation or uniform numerical margin is required.

## Scope and reproduction

The result does not say that S75 cannot acquire a finite-area boundary: such walls may exist while missing this particular escape. It says that completing that face and cutting the saved global escape cannot be accomplished by this boundary-producing family alone. A successful mixed construction must also use a wall outside that family, or change the retained conditions used in the identity. This distinction prevents promotion of a support-local gate into a global nonexistence claim.

The [exact producer/verifier](../scripts/fresh_h20_s75_completion_20261014_gate.py) binds the current chamber, its escape, the independently verified support Gram and all arithmetic dependencies.

```
PYTHONPATH=src:. python scripts/fresh_h20_s75_completion_20261014_gate.py --verify
```

Generation and saved replay pass for both the exact gate and the numerical scout. The [independent coefficient and proof review](fresh_h20_s75_completion_20261014_gate_review.md) reconstructs all physical coordinate identities and signs with separate arithmetic and verifies full rank by a finite-field determinant. Frozen predecessor files were not modified.
