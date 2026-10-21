# An exact compact triangle face in the 153-wall chamber

One exact finite-contact cap extends the
[pi/8 alternative 152-wall chamber](fresh_h20_second_cusp_pi8_cap_20261016.md)
to **153 actual facets**. It completes the codimension-18 face supported
by `S16+{150,151}` to a compact hyperbolic Coxeter triangle with angles
`pi/2, pi/4, pi/6`. Two of that triangle's three vertices are new finite
H20 vertices. The certified total is now **at least 174 finite
vertices**, and the five proper cusps persist.

The ambient chamber still has infinite volume: the old exact strict
escape E6pi8 remains. This closes a whole two-dimensional face, not the
four-dimensional residual face or the whole H20 chamber. The optional
all-dotted 152-wall cutoff is not retained.

The [source](../scripts/fresh_h20_pi8_triangle_face_20261016.py) and
[certificate](fresh_h20_pi8_triangle_face_20261016.json) give the exact
recipe, every new contact, two complete new finite-corner checks and
the remaining escape. Generation and a fresh saved replay pass; the
separate numerical diagnostic also passes its tolerance replay.
Independent review of this new package is pending.
All predecessor files remain frozen.

## Exact root and all old contacts

Let C8 be the finite corner of the pi/8 package, perpendicular to
`S16+{142,145,150,151}`, and set `k8=-C8^2>0`. Define

```
r152 = -n142/sqrt(q142) - (sqrt(3)/2)n145 + C8/sqrt(2*k8),
q142 = 309/200.
```

The three displayed summands are pairwise orthogonal. Their squared
norms are `1, 3/2, -1/2`, so `r152^2=2` exactly. The new root is right
to S16 and both roots 150 and 151 of the retained I2(8) component.
Its products with 142 and 145 are respectively `-sqrt(q142)` and
`-sqrt(3)`, giving pi/4 and pi/6. The full 152-wall contact check gives
**18 right, one pi/4, one pi/6 and 132 strictly dotted** contacts.

The retained strict future center remains inside. Acute projections
from it certify every one of the 153 walls is an actual facet. Each
of the five old proper cusp rays has strictly negative pairing with
r152, so its bounded link persists.

## Complete finite-corner census and the actual compact face

The new root has exactly 20 finite old neighbors:

```
S16 + {142,145,150,151}.
```

There are 20 ways to choose 19 of these old roots for a finite vertex
containing r152. Exactly two work: omit 142 or omit 145. Their diagrams
are respectively

```
A10 + A4 + A2 + I2(8) + I2(6),
A10 + A4 + A2 + I2(8) + I2(4).
```

Both Gram matrices are positive definite. Projecting the common center
onto each 20-wall intersection gives a future timelike point. Every
one of its 153 pairings is checked, including strict negativity away
from its incident walls.

Each of the other 18 candidate supports contains the indefinite
three-node diagram on `{142,145,152}`. After normalizing roots to
squared norm 2, its Gram matrix is

```
[[2,       0, -sqrt(2)],
 [0,       2, -sqrt(3)],
 [-sqrt(2), -sqrt(3), 2]].
```

Its determinant is `-2`; all proper two-node determinants are positive.
Thus this is a complete census of new finite vertices, not a bounded
corner search. Every old finite vertex persists, including C8.

The support `S16+{150,151}` is positive definite of rank 18. Its
orthogonal complement is a Lorentz 3-space, and the three side roots
142, 145 and 152 all lie in it. Their three proper two-node diagrams
are elliptic and their full diagram is the displayed Lanner triangle.
The three future finite vertices are C8 and the two newly checked
corners. The new wall is strictly negative at C8, since
`r152.C8=-sqrt(k8/2)`.

All other retained inequalities hold at these three vertices and
therefore throughout their hyperbolic convex hull. Conversely the
three side inequalities bound precisely their compact triangle.
Consequently this is the **whole actual codimension-18 face**, with
area `pi/12` at curvature -1. No general assumption that arbitrary
spherical faces are Coxeter faces is used: these side roots are
explicitly perpendicular to the support, and all vertices and retained
inequalities are checked.

## Remaining escape and residual diagnostic

The exact E6pi8 from the previous package has positive squared norm
and strictly negative pairing with all 153 walls, including the new
root (approximately `-0.4174481876`). Its future-chart sign remains
negative. This proves that the ambient chamber has infinite volume.

The separate [numerical source](../scripts/fresh_h20_pi8_triangle_face_20261016_residual_scout.py)
and [saved comparison](fresh_h20_pi8_triangle_face_20261016_residual_scout.json)
track all projected inequalities in the residual four-dimensional
affine polytope. Floating-point Qhull reports:

| Residual polytope | Vertices | Clearly inside | Within 1e-6 of sphere | Clearly outside |
|---|---:|---:|---:|---:|
| pi/8 152-wall base | 154 | 14 | 19 | 121 |
| Triangle-face 153 walls | 163 | 16 | 19 | 128 |

The cap cuts seven exterior vertices and introduces two finite and
fourteen exterior vertices. These counts are numerical diagnostics,
with nearly degenerate projected incidences left subject to exact
verification. They do not constitute a complete cone or ideal-vertex
census. The exact compact triangle and two new finite vertices above
are independent of this diagnostic. The remaining exterior frontier
has not closed.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_pi8_triangle_face_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_pi8_triangle_face_20261016_residual_scout.py --verify
```
