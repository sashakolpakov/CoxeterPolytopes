# A distinct local quadratic space, but an infinite reflection chamber

The lattice

\[
L=\langle-\phi\rangle\perp I_7\perp(4-\phi)I_2,
\qquad \phi^2=\phi+1,
\]

is **nonreflective**. Two explicit integral future isometries preserve its full reflection chamber; their product has trace \(50+76\phi>10\), hence infinite order. The chamber has infinite volume. In particular, no finite-volume Coxeter polytope can have all its wall reflections preserving this literal lattice. This is an internally checked exact obstruction, with independent verification of the starting lattice, bounded root census and full-chamber symmetry. It is not a finite-volume example or a classification of other lattices in its quadratic space.

The [setup](fresh_h9_witt_twist_20261016.md) proves that this space is not even rationally similar to any \(E_8\)-plus-binary space: its local Witt index at \((4-\phi)\) is three, whereas the latter spaces contain a split eight-dimensional subspace there. Thus this experiment is structurally different from the earlier discriminant11/55/131 runs. That distinction does not prevent the different symmetry found here. The [independent setup review](fresh_h9_witt_twist_20261016_review.md) checks the local invariant, admissibility, complete four-norm alphabet and complete zero sector.

## Complete height20 prefix

Write \(\beta=4-\phi\), use coordinates \((x,y_1,\ldots,y_7,z_1,z_2)\), and take \(H=e_0\). The four primitive reflection-norm representatives are \(1,2,\beta,2\beta\). All 106 unit-normalized oriented zero roots form \(B_7+B_2\), whose nine simple walls impose

\[
y_1\ge\cdots\ge y_7\ge0,\qquad z_1\ge z_2\ge0.
\]

For a positive-height root of norm \(q\), the norm equation and height are

\[
\sum y_i^2+\beta\sum z_j^2=q+\phi x^2,
\qquad h=\phi x^2/q,\qquad x>0.
\]

The [producer](../scripts/fresh_h9_witt_twist_pilot_20261016.py) and [saved certificate](fresh_h9_witt_twist_pilot_20261016.json) enumerate every root with \(h\le20\) in the zero chamber, using explicit coefficient boxes, both real embeddings, all four norms and the long-root divisibility filter. Chronological selection yields 13 control levels, 1,537 recursion nodes and **16 walls**, consisting of the nine zero walls and seven positive-height walls. All 120 wall pairs pass exactly: 86 right, eleven label3, two label4, five label5 and sixteen dotted. A common strict timelike center and sixteen actual facet points are checked.

The [independent prefix review](fresh_h9_witt_twist_pilot_20261016_review.md) uses a different arithmetic implementation, larger boxes and a 3+4+2 spatial-block join. It reconstructs 972 block candidates and 32 complete dominant norm vectors before chronological filtering, reproducing every control level and all sixteen walls. Its generation and fresh saved replay pass.

## Geometry and the two maps

The [cone checker](../scripts/fresh_h9_witt_twist_pilot_20261016_geometry.py) and [certificate](fresh_h9_witt_twist_pilot_20261016_geometry.json) give a complete cone census of **37 finite vertices, no ideal vertices and 165 hyperideal rays**, with no rays on the opposite sheet. The prefix has an exact strict spacelike escape of square

\[
\frac{34117}{968}-\frac{656957}{30976}\phi>0.
\]

All sixteen escape inequalities are strict. This proves that the bounded prefix has infinite volume; the symmetry argument below proves the stronger statement for the entire integral reflection chamber.

The [map screen](../scripts/fresh_h9_witt_twist_pilot_20261016_maps.py) and [saved data](fresh_h9_witt_twist_pilot_20261016_maps.json) tested 149 norm-colored corner-frame maps, retaining 44 integral maps. Two of them give the required product. The obstruction itself is isolated in the [minimal symmetry checker](../scripts/fresh_h9_witt_twist_pilot_20261016_symmetry.py) and [certificate](fresh_h9_witt_twist_pilot_20261016_symmetry.json); it requires neither completeness of the corner census nor completeness of the map screen.

| Generator | Ordered source wall IDs | Ordered target wall IDs | Positive apex scale | Determinant | Trace |
|---|---|---|---|---|---|
| \(T_0\) | 0,1,2,4,6,7,8,9,11 | 4,3,2,0,6,7,8,11,9 | \(3+4\phi\) | −1 | 4 |
| \(T_1\) | 0,1,2,3,4,5,6,8,10 | 10,14,2,3,4,5,6,8,0 | \(2+2\phi\) | 1 | 6 |

The source/target vertex indices are respectively \((11,24)\) and \((2,22)\). For each pair, the checker reconstructs both physical apex vectors and their complete nine-wall incidence, verifies their positive-definite spherical Gram matrices and nonnegative inverse matrices, and checks all sixteen wall inequalities. It reconstructs \(T_i\) from the ten-vector frames. Each matrix and inverse is integral, preserves the literal form, maps the nine outward normals exactly in the displayed order, and sends the future apex to the stated positive multiple of the future target apex.

The product \(T=T_0T_1\) is checked directly:

\[
T^tGT=G,\qquad T,T^{-1}\in\mathrm{GL}_{10}(\mathbb Z[\phi]),
\qquad \det T=-1,\qquad \operatorname{tr}T=50+76\phi>10.
\]

The certificate stores the literal generators, product and inverses. A finite-order real matrix has eigenvalues of modulus one and trace at most its dimension, so this product has infinite order.

## Why the maps stabilize the full chamber

The complete primitive norm alphabet and zero sector make the height20 walls a genuine prefix of the full Vinberg chamber. At any of the four checked vertices \(v\), write its active normals as \(n_1,\ldots,n_9\) and their Gram matrix as \(S\). Exact arithmetic verifies

\[
H=a v+\sum_i b_i n_i,
\qquad a>0,\qquad b=S^{-1}(H\cdot n_i)_i\le0.
\]

Every later simple wall \(r\) satisfies \(H\cdot r<0\) and \(n_i\cdot r\le0\). Consequently

\[
0>H\cdot r=a(v\cdot r)+\sum_i b_i(n_i\cdot r)\ge a(v\cdot r),
\]

so \(v\cdot r<0\). The corner remains in the full chamber and acquires no further incident walls. The same conclusion applies at each source and target corner.

An integral future isometry normalizes the full lattice-reflection group, because it carries integral reflection directions to integral reflection directions. Here it also carries the complete outward local chamber sector at one full corner to that at another. Local finiteness of the discrete mirror arrangement therefore identifies the image chamber with the original chamber. Both \(T_0\) and \(T_1\), hence their product, stabilize it.

The conjugate form is positive definite. Integral matrices bounded at the physical embedding are bounded at the conjugate embedding as well, so the admissible integral isometry group is discrete. An infinite subgroup preserving the chamber produces infinitely many disjoint translates of a sufficiently small interior ball, chosen about a point with trivial stabilizer. Thus the chamber has infinite volume. Any reflection subgroup generated by walls preserving this lattice lies inside the full lattice-reflection group; a finite-covolume such subgroup would force finite covolume for the full group, a contradiction.

This leaves nonintegral geometric reflections and other literal lattices open. It makes no claim that determinant squareclass alone distinguishes hosts, that all lattices in this rational quadratic space are nonreflective, or that the full chamber stabilizer has been classified.

## Replay and provenance

All four producer scripts below have passed generation and saved replay. Every saved output binds its source and exact input artifacts; the symmetry checker binds only the setup, bounded prefix, literal corner data and arithmetic helper. The full cone and map search are useful discovery records, while the nonreflectivity proof uses only the four verified corners and two literal maps.

```sh
python3 scripts/fresh_h9_witt_twist_pilot_20261016.py --verify
python3 scripts/fresh_h9_witt_twist_pilot_20261016_geometry.py --verify
python3 scripts/fresh_h9_witt_twist_pilot_20261016_maps.py --verify
python3 scripts/fresh_h9_witt_twist_pilot_20261016_symmetry.py --verify
```

The independent setup and prefix reviews are linked above. The separate [independent symmetry review](fresh_h9_witt_twist_symmetry_20261016_review.md) passes generation and fresh saved replay. It reconstructs the four apices by projection of \(H\), uses separate project-field arithmetic and linear algebra, and checks both literal maps, their inverses, outward sectors and product. This producer's cone census is not an independence premise. No higher-height enumeration is needed for this obstruction.
