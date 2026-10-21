# Complete ambient cone of the58-wall H9 chamber

The [three-wall tetrahedral construction](fresh_nonarithmetic_pi8_tetrahedron_completion_20261017.md)
has exactly **184 finite vertices, no ideal rays, and2,614 exterior extreme
rays** in its complete58-wall cone. All179 old finite vertices persist; the
five new vertices are the endpoint on the old[7,30] edge and the four vertices
of the new compact tetrahedral face.

The chamber still has infinite volume. The certificate constructs an exact
future-oriented spacelike vector strict against all58 walls, so an open patch
of the ideal boundary remains. This is a complete cone calculation, not merely
a check of previously known exterior rays.

Starting from the [independently verified55-wall cone](fresh_nonarithmetic_nonright_s6_20261016_cone_review.md),
the successive exact clips give:

| Added wall | Removed exterior rays | New crossing rays | Finite vertices | Exterior rays |
|---|---:|---:|---:|---:|
|55, the pi8 endpoint cap|15|45|180|2,538|
|56, first tetrahedron side|15|45|181|2,567|
|57, second tetrahedron side|35|85|184|2,614|

The count of exterior extreme rays is not a measure of the size of the open
region: clipping can replace a removed extreme ray by several new ones.

Every old ray has a strict sign against each newly added halfspace. All
opposite-sign pairs with at least eight common active walls are examined.
A certified rank-eight minor proves that the common face is an edge: the
two distinct parent rays already give the opposite rank bound. The positive
crossing combination therefore supplies exactly the new extreme ray.

Repeated interval clipping initially obscured the sign of one small positive
norm. The final algorithm reconstructs each crossing ray from nine independent
active equations, fixing a certified nonzero coordinate to its sign. The
invertible minor proves that this is exactly the same projective ray after
positive normalization. It reduces interval dependency without changing the
geometry. Every final norm has a strict certified sign. All162,284 final
wall inequalities are checked, including strictness at every nonincident wall;
active equalities follow from the clipping and reconstruction identities.

The complete original H4+A2 three-face now has14 finite and8 exterior rays.
The two added tetrahedron sides do not change that face. The distinct
H4+I2(8) three-face has exactly the four finite tetrahedron vertices and no
exterior rays, confirming whole-face closure inside the ambient cone.

For the remaining escape, take the exterior ray with active set
`[0,1,2,3,5,6,7,33,49]` and add `2^-34` times the saved strict center.
The resulting exact vector has square approximately `7.19447e-11`, negative
product with the center, and strictly negative product with every wall. A
segment from the center reaches the future null cone with all wall
inequalities strict, proving infinite volume.

The [source](../scripts/fresh_nonarithmetic_pi8_tetrahedron_cone_20261017.py)
and [certificate](fresh_nonarithmetic_pi8_tetrahedron_cone_20261017.json)
are a primary exact clipping calculation using the independently reviewed
old cone and the new exact root definitions. A separate agent read the
clipping completeness, coordinate conventions and direct reconstruction
argument; this is a conceptual source audit, not independent new-root
arithmetic or external review. Generation and a fresh saved replay pass.

```sh
python3 scripts/fresh_nonarithmetic_pi8_tetrahedron_cone_20261017.py
python3 scripts/fresh_nonarithmetic_pi8_tetrahedron_cone_20261017.py --verify
```
