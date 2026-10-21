# Independent review of the actual compact triangle face

**PASS.** The [independent reviewer](../scripts/fresh_h20_pi8_triangle_face_20261016_review.py)
and [certificate](fresh_h20_pi8_triangle_face_20261016_review.json)
reconstruct the [153-wall chamber](fresh_h20_pi8_triangle_face_20261016.md)
and its whole compact triangle face. Generation and a fresh saved replay
both pass. The primary report's pending-review statement is historical.
A separate conceptual review checks the face-closure and finite-support
arguments below, without repeating the physical computation.

The preceding 152-wall physical geometry is independently rebuilt. The
new root is the sum of three orthogonal vectors with squares 1,3/2,-1/2,
giving norm 2 exactly. Its prescribed contacts are eighteen right angles,
pi/4 to wall 142 and pi/6 to wall 145. Every remaining old contact is
strictly dotted. Exact recipe identities justify the equalities; rational
outward intervals certify the 132 strict inequalities, all 153 actual
facets, and preservation of five proper cusps.

Two general four-by-four Gram inversions, each with positive pivots,
independently construct the new finite corners. The reviewer does not
reuse the producer's intermediate H18/two-root projection formula. Both
corners are future timelike, have exactly the stated twenty incidences,
and satisfy every other wall inequality strictly. The new root has
twenty finite old neighbors. Of the twenty possible nineteen-element
supports, eighteen contain the indefinite triangle on 142,145,152;
the remaining two are exactly the checked corners. Thus there are two
new finite vertices, bringing the known total to at least 174.

The retained support S16+{150,151} is positive definite of rank 18.
The three side normals 142,145,152 lie in its Lorentz three-dimensional
orthogonal complement. Their normalized Gram determinant is
8-2*2-2*3=-2 and their three proper pair determinants are 4,2,1.
Together with the checked future orientation, these give the compact
hyperbolic triangle with angles pi/2,pi/4,pi/6 and area pi/12.

Its three vertices are the preceding pi/8 corner and the two new corners.
Every additional wall inequality holds at these vertices, hence throughout
their positive cone and hyperbolic convex hull. Conversely the three side
inequalities bound that compact triangle. Consequently this is the whole
actual codimension-18 face, not a triangle inside a larger or truncated
face. The argument uses explicit perpendicular side roots and does not
assume arbitrary spherical-support faces are Coxeter faces.

The old exact E6 escape still pairs strictly negatively with all 153
normals and has positive square. The ambient chamber therefore retains
infinite volume. This review proves neither a compact H3/H4 face nor a
finite-volume H20 example. The ambient finite-vertex census remains
partial, and the old isolated parameter and chamber acuteness are named
inherited inputs. Separate floating residual-cone scores are not proof inputs.

```sh
python3 scripts/fresh_h20_pi8_triangle_face_20261016_review.py --verify
```
