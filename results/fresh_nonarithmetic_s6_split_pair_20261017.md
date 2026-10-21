# Two separate endpoint caps joined by a finite corner

This exact construction adds walls 63 and 64 to the frozen 63-wall H9 chamber. They terminate the old edges `[30,44]` and `[30,52]` separately and meet at angle π/4. Their common corner gives a connected three-corner patch. All 65 walls are actual facets, the complete finite-vertex count is 196, and no affine circuit is introduced. The previous strict escape survives; this is not a finite-volume completion or a closure of the whole S6 face.

The source is [the exact checker](../scripts/fresh_nonarithmetic_s6_split_pair_20261017.py), with [saved data](fresh_nonarithmetic_s6_split_pair_20261017.json). Its inputs are the [63-wall cap](fresh_nonarithmetic_escape62_20261017_cap.md), [complete cone](fresh_nonarithmetic_escape62_20261017_cone.md), and the frozen endpoint-domain label table. The exploratory numerical scout is not a proof dependency.

## Exact recipe

Put `S6=[0,1,2,3,5,6]`. The two positive definite rank-eight frames are

```
F_A = S6 + [30,44]
F_B = S6 + [30,52].
```

For each frame, project the retained strict center H into its Lorentzian orthogonal plane and normalize the result to a future unit timelike vector k0. Project old wall 52 for the first frame, or old wall 44 for the second, and remove its k0 component. Normalize this remaining spacelike vector to k1. Thus `k0²=-1`, `k1²=1`, and `k0·k1=0`.

The perpendicular norm-two family is

```
r(t) = t (k0+s k1)/sqrt(2) + (-k0+s k1)/(sqrt(2)t),  s in {-1,1}.
```

Set the first parameter to exactly 7. Exactly one of its two signs satisfies all the required strict dotted gates. Write the second family as `u a+b/u`; its contact with the first wall is π/4 precisely when

```
(r(7)·a) u² + sqrt(2) u + (r(7)·b) = 0.
```

Both signs and both quadratic roots are checked. Exactly one branch with `u>1` passes all old gates; its parameter is approximately `5.528687752868377`. All normal coordinates are defined by this algebraic recipe, not by rounded numerical fitting.

The implementation uses literal finite frame Grams and the identities

```
q0 = H² - hᵀ G⁻¹ h
w·k0 = (seed·H - cᵀ G⁻¹ h)/sqrt(-q0)
q1 = 2 - cᵀ G⁻¹ c + (w·k0)²
```

to avoid dependency growth from expanded coordinate products. Exact identities establish the prescribed zero contacts and norms; rational outward bounds certify all strict signs and discriminants. The arithmetic engine is inherited from the frozen geometric packages; this is a primary calculation, not an independent arithmetic review.

## Contacts and actual vertices

Wall 63 is right to `F_A`. Wall 64 is right to `F_B` and meets wall 63 at π/4. All 110 other new contacts are strictly dotted. The largest dotted products are approximately −2.04172047 and −2.04500611, both with old wall 8. The checker records all 127 new contact decisions.

The old strict center remains strictly inside all 65 halfspaces. Projecting it onto each wall gives a future timelike point satisfying every other wall inequality strictly, proving that every wall is an actual facet.

The new finite vertices have exactly these supports:

```
S6 + [30,44,63]
S6 + [30,52,64]
S6 + [30,63,64].
```

Their positive definite supporting Grams, future timelike intersection rays, and all nonincident strict inequalities are checked. Wall 63 has exactly eight finite preceding neighbors, giving one possible vertex support. Wall 64 has nine, giving nine supports: two are the displayed vertices and the other seven contain the dotted pair `[52,63]`. This exhausts all new finite vertices, including potential nonsimple ones. The checker also verifies that all 193 prior finite vertices remain strict against both additions. Consequently the complete finite census is 196, without requiring a full cone enumeration.

The preceding chamber has no affine circuit. The only nonright finite bond involving either new root is the isolated π/4 bond `[63,64]`; its component is spherical I2(4). Hence neither addition creates an affine circuit.

## Remaining open geometry and export

The previous exact strict escape remains strict against both new walls, with products approximately −1476.61628 and −712.64479. It therefore still proves infinite volume of this particular 65-wall chamber. This pair avoids the earlier obstruction to one wall terminating both old edges by assigning the two endpoint duties to different walls. It does not settle whether the resulting patch admits a compact completion.

`geometry()` returns all 65 computational normals, their bilinear form, the basis, and the common center. The basis is unchanged:

```
[0,1,2,3,51,55,7,30,5,6].
```

The JSON exports `new_normals[].computational`, `new_normals[].physical`, `physical_common_center`, `new_finite_vertices`, and the complete sequential finite-support census. A full 65-wall cone and its boundary ledger are separate work.

Validation: generation and fresh saved replay both passed. Direct source hashes, syntax, and report links were checked before freezing. The source, saved data, and report are final; any continuation uses a new prefix.
