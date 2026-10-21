# An exact compact prism on the π9 face

Two new coherent Coxeter walls close the entire actual H3 face supported by

`S'=[0,1,2,3,58,59]=H4+I2(9)`

as a compact triangular prism. The resulting ambient chamber has62 genuine facets. The construction creates five finite vertices: two when wall60 is added, then three more when wall61 is added. The old60-wall strict escape remains strictly feasible against both new walls, so the ambient chamber still has infinite volume.

The [exact checker](../scripts/fresh_nonarithmetic_pi9_prism_completion_20261017.py) and [certificate](fresh_nonarithmetic_pi9_prism_completion_20261017.json) bind the frozen [60-wall construction](fresh_nonarithmetic_pi8_quad_completion_20261017.md) and [complete60 cone](fresh_nonarithmetic_pi8_quad_cone_20261017.md). Generation and a fresh saved replay pass. These are primary exact checks using the shared physical interval engine; a separate reconstruction or full62-cone computation is not claimed here.

## The two normals

Write `φ=(1+sqrt(5))/2`. The prism retains old sides5,6,55. New wall60, denoted `v`, is right to `S'` and5, has angleπ/5 with6, and angleπ/4 with55. All51 other old contacts are strictly dotted.

The old side Gram on `[5,6,55]` is `A2+A1`. Let `p` be its exact dual-frame projection with products `[0,−φ,−sqrt(2)]`, and let `k` be the future unit timelike vector at the old finite corner `S'+[5,6,55]`. Then

`v=p+sqrt(p²−2)*k`.

Here `p²−2>0`; the chosen sign puts the new wall strictly beyond the persistent old finite vertex. The opposite sign would exclude that vertex.

New wall61, denoted `u`, is the oriented norm2 polar of the hyperideal corner `S'+[6,55,v]`. Project old5 to the one-dimensional positive orthogonal complement of that nine-root Lorentzian span, obtaining `w`; then

`u=−sqrt(2/w²)*w`.

Thus `u` is right to `S',6,55,v`, and all52 remaining contacts are strictly dotted. In particular, it is dotted to5. The exact prism Gram, in order `[5,6,55,v,u]`, is

```
[ 2  -1    0       0       -d ]
[-1   2    0      -φ        0 ]
[ 0   0    2     -sqrt(2)   0 ]
[ 0  -φ  -sqrt(2)  2       0 ]
[-d   0    0       0        2 ]
```

where `d²=4+2φ` and `d≈2.6899940479`. The latter identity also follows directly by eliminating the `[6,55,v]` block: its inverse entry at6 is `−φ`, hence `w²=2+φ`.

The121 new pair relations are checked exactly. A common strict future center is constructed by perturbing the retained finite corner toward the old center. Projection of this center to each facet gives62 strict future facet witnesses. These checks certify a genuine acute Coxeter chamber, including both new facets.

## The whole face is compact

The triangular bases are old5 and new61. The rectangular sides are6,55,60. In the four-dimensional Lorentzian complement of `S'`, the five side normals span the whole space and have a strict future common interior.

All ten triples of side hyperplanes are examined. Exactly six projective intersection rays are feasible and future timelike:

- `[5,6,55]`, `[5,6,60]`, `[5,55,60]`;
- `[6,55,61]`, `[6,60,61]`, `[55,60,61]`.

Both orientations of each of the other four intersection lines violate a side inequality. Every one of the six feasible vertices satisfies all other native wall inequalities strictly. Consequently the six-ray cone is exactly the full five-side cone, and its compact projectivized convex hull is the entire actual `S'` face. This proves compact prism closure, rather than only six compatible local corners.

The finite-vertex census is sequential. Wall60 has nine finite old neighbors; of its nine possible eight-neighbor supports, two give actual finite vertices and seven contain the indefinite `[6,55,60]` triple. Wall61 again has nine finite preceding neighbors; three supports give actual finite vertices and six retain the same indefinite triple. Old finite vertices persist under the coherent acute extension. Two of the vertices created by61 also lie on60, so the final wall60 is incident to four of the displayed prism vertices.

A complete connected affine-circuit screen on the62 native normals finds673 positive-definite connected supports, two indefinite supports, and no affine circuit. This is a persistent-cusp check, not a volume criterion.

## Remaining ambient region and exports

The old60-wall strict spacelike escape has new products approximately

`E·v=−0.0006433840308493`, `E·u=−0.0013288291549334`.

Its norm is strictly positive, its future-chart product is negative, and all62 wall products are strictly negative. Interpolating with the common interior gives a strict ideal point and an open ideal region. Thus this prefix still has infinite volume; the completed compact prism is a local construction milestone.

For a subsequent independent cone replay, the JSON exports `new_normals[].physical`, `physical_common_center`, `face_support`, and `complete_triple_intersections[].incident`. The source `geometry()` returns all62 computational normals, their bilinear product, and the original physical basis; physical conversion is `P.combine(vector,z['basis'])`. The old sources and60-wall certificate remain unchanged.
