# Independent 810-wall review and local cusp diagnosis

**PASS.** The frozen [810-wall producer](../scripts/fresh_h20_balanced_integer_exact_20261014.py) passes a fresh saved replay. The separate [review checker](../scripts/fresh_h20_balanced_integer_exact_20261014_review.py) and [review certificate](fresh_h20_balanced_integer_exact_20261014_review.json) reconstruct the literal chamber using standard-library integers and `Fraction`, without importing the producer, FLINT, or project geometry. No source or mathematical issue was found.

The independent replay reconstructs all 810 normals from the type-A inverse formula and saved recipes. It checks all 327,645 pairs: 10,405 right, 6,318 at `pi/3`, 410 parallel, and 310,512 dotted. It verifies a strict timelike center, all 392 old finite endpoints and every inequality there, all 392 designated `2z/5` cuts, and the norm-one strict future spacelike escape at `3z_0/5`. Strict center and nonpositive pairings justify all 810 actual facets by orthogonal projection. The old endpoints have positive-definite incidence Gram `A6 + A6 + A7 + A1`. These are retained endpoints, not a complete new-vertex census.

The escape proves that the **present 810-wall chamber has infinite volume**. Neither the escape nor a deficient cusp link proves that all further compatible completions have infinite volume.

## One complete cusp-link diagnosis

Use zero-based wall indices. With unit seed roots, the future null vector

`delta = r_3 + r_4 + r_418/2`

has exactly the incident set

`{3,4,9,12,15,17,20,418}`.

The normalized first triple is affine `A2`, with positive kernel `(1,1,1)` and sum `delta`. The other five roots, indexed `9,12,15,17,20`, are mutually orthogonal unit roots and orthogonal to that triple. Thus the incident diagram is affine `A2` plus five spherical `A1` components. Its Gram rank is seven; its normal-vector rank is eight. Every other one of the 802 walls satisfies `(delta,r) <= -1/2`.

The horospherical tangent space `delta^perp / R delta` has dimension 19. This link has 12 missing tangent directions and five one-sided spherical factors. The vector `-r_9` is already an explicit positive-norm recession direction: it is orthogonal to `delta` and nonpositive on every incident root. Consequently this link is unbounded. This is an actual ideal boundary point; no claim that it is a combinatorial ideal vertex is needed.

## Local completion lemma and its scope

Let a finite, full-dimensional acute Coxeter chamber have strict timelike center `H`, and let `delta` be a nonzero feasible future null vector. Suppose its incident diagram decomposes into positive affine components and spherical components which admit Euclidean affine Coxeter extensions. Then, allowing arbitrary real new normals, the link at this **one** null direction can be made bounded by finitely many additional Coxeter walls, retaining every old wall and `H`.

To see this, work in the positive-definite quotient `delta^perp / R delta`. Complete each spherical component with its negative-highest-root normal, and choose orthogonal opposite pairs spanning any remaining quotient directions. Lift these normals to `delta^perp`, then replace every new lift `b` by `b+t delta`. Norms and all pairings between incident or new normals are independent of `t`. All required affine Coxeter angles are therefore unchanged. Positive shifts can make each new affine kernel sum a positive multiple of `delta`.

For an old nonincident wall `r`, one has `(delta,r)<0`. Hence, for all sufficiently large positive shifts, `(b+t delta,r)` is strictly below the parallel threshold, so the relation is dotted. The same shift makes `(H,b+t delta)<0`. There are only finitely many inequalities, so all can be satisfied simultaneously. The resulting incident quotient spans the full tangent space, with a positive affine kernel in every component. A recession vector nonpositive on all its normals must pair to zero with every normal, and therefore vanish. This proves boundedness of the link; rank alone would not suffice.

The hypotheses hold for the representative cusp above. The five `A1` factors are completed by `-r_i+t delta`; the checker verifies that **the explicit common choice `t=100`** gives norm-one mates, preserves all 810 old Coxeter relations and the strict center, and makes the five mates mutually orthogonal. Choose 12 orthogonal quotient directions perpendicular to the old rank-seven quotient span and add their opposite pairs with sufficiently large positive null shifts. This gives an existence proof using 29 additional walls, whose final link is affine `A2 + 17 A1`, of total Euclidean rank 19. No 839-wall construction or global geometry assertion is being archived here.

The same existence argument works simultaneously for any **finite set** of distinct feasible future null directions satisfying the hypotheses. Choose all new shifts equal to a sufficiently large `T`. For roots assigned to different directions, their cross pairing has leading term `T²(delta_i,delta_j)<0`; the remaining terms are at most linear in `T`. Thus all cross-cusp relations become dotted. Each new wall is also strict at the other null directions. Same-cusp relations remain unchanged. This avoids conflicts that particular smaller, uncoordinated shifts could create.

Even this finite-set observation does not produce finite volume, repair the saved escape, or imply a compact completion: it does not control other ideal boundary directions or further unbounded regions of the resulting chamber. The representative `A2 + 5 A1` cusp supplies **no structural obstruction to a parabolic repair** compatible with the full seed, central cap, and all remaining old walls.

```sh
python scripts/fresh_h20_balanced_integer_exact_20261014.py --verify
python scripts/fresh_h20_balanced_integer_exact_20261014_review.py --verify
```
