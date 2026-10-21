# A forced noncrystallographic cusp closes the fixed-core completion route

**The fixed 44-, 45-, and 46-wall H9 cores cannot be completed to finite volume by adding acute-compatible walls while retaining their inherited facets.** Four already-present walls suffice. This is stronger than the explicit escapes of the saved chambers, but is not an obstruction to replacing inherited walls or completing the ten-wall seed by a different construction.

The [exact certificate](../data/enumerations/fresh_h9_fixed_core_cusp_obstruction_20261013.json) and [standalone checker](../scripts/fresh_h9_fixed_core_cusp_obstruction_20261013.py) reconstruct, from the authoritative seed and old-cap archive, the four normals

```
a = old31, b = seed8, r = seed6, s = seed7,

Gram(a,b,r,s) =
 [ 2 -2    0    0 ]
 [-2  2    0    0 ]
 [ 0  0    2 -phi ]
 [ 0  0 -phi    2 ],       phi^2=phi+1, phi>0.
```

All four have norm two. The sum `v=a+b` is nonzero and null, is future-directed relative to the saved strict timelike center, and is orthogonal to all four normals. The reflection product in `r,s` has trace `7+phi` on the ten-dimensional Lorentz space and trace `5+phi` on the eight-dimensional Euclidean quotient `v-perp / Rv`.

## Why further acute walls cannot repair this cusp

Suppose a full-dimensional finite-facet Coxeter completion retains these four walls and has nonpositive pairings between every pair of distinct outward normals. For any additional normal `w`,

```
(v,w)=(a,w)+(b,w)<=0.
```

Thus no such wall removes the ideal direction of `v`. Choose a strict timelike point `x` of the hypothetical completion. Then `v+epsilon*x` is strictly feasible and timelike for sufficiently small positive `epsilon`, so `v` is still in its ideal boundary. **The saved center need not survive the new cuts**; it is used only to certify the original orientation and that `v` is nonzero.

Send `v` to infinity in the upper half-space model. Incident facets become vertical planes. Above all nonincident hemispheres, the chamber is the product of a Euclidean polyhedron `C` in dimension eight with a vertical half-line. It has a strict interior because the completion does. If `C` were unbounded, its Euclidean volume, and hence this hyperbolic cusp volume, would be infinite. Finite volume therefore forces `C` to be bounded.

The inherited walls `r,s` remain facets of this link. Indeed, project a strict point of `C` orthogonally onto either wall. Every other inequality remains strict because its outward normal has nonpositive pairing with the projected wall's normal. The positive-definite two-wall Gram also permits simultaneous projection, giving their actual pi/5 ridge. Their reflections therefore belong to the cocompact Euclidean reflection group of the bounded Coxeter link.

Bieberbach's theorem gives this group a full-rank translation lattice. Its full translation subgroup is invariant under conjugation, so the linear part of every group element preserves a lattice and has integer trace. This is the uniform case of Theorem 2.9 in Tsachik Gelander's author-hosted [*Lectures on Lattices*, page 5](https://www.math.utah.edu/pcmi12/lecture_notes/gelander.pdf).

But the product of the `r,s` reflections fixes six Euclidean directions and rotates their two-plane through `2pi/5`. Its linear trace is

```
6 + 2cos(2pi/5) = 5 + phi,
```

which is irrational. This contradiction proves the claimed obstruction. Equivalently, the exact ambient trace loses two eigenvalues equal to one when passing to `v-perp / Rv`, giving the same result without numerical trigonometry.

## Replay and scope

```sh
python -m scripts.fresh_h9_fixed_core_cusp_obstruction_20261013
python -m scripts.fresh_h9_fixed_core_cusp_obstruction_20261013 --verify
```

The checker uses Fraction arithmetic, reconstructs the Lorentz signature, Gram inverse, four normals, null ray, strict center and reflection matrices, and rejects forged Gram and trace entries. It does not perform a new search or a 46-by-46 pair census. Separate root-agent calculation also confirmed the four-wall Gram, null ray and quotient trace.

The old cap `old31` is present in all three literal cores, so this closes their finite-volume completion route by acute-compatible addition. The 43-wall checkpoint that deletes old caps 31, 51 and 197 is outside this four-wall theorem. Its separately verified positive affine circuits obstruct compactness by acute additions, but this particular argument does not exclude its finite-volume completion. The ten-wall seed alone is likewise not excluded.
