# Independent review of the 158-wall H20 cap and compact quadrilateral

**PASS.** Independent generation and a fresh saved replay reproduce the new
root, all 157 old contacts, 158 actual facets, exactly two new finite vertices,
and the entire compact quadrilateral on `S16+{151,156}`. The old strict escape
survives, so this is local geometric progress in an infinite-volume chamber.

The [review source](../scripts/fresh_h20_compact_quad_cap_20261016_review.py)
and [saved certificate](fresh_h20_compact_quad_cap_20261016_review.json) bind
the frozen [primary source](../scripts/fresh_h20_compact_quad_cap_20261016.py)
and [primary certificate](fresh_h20_compact_quad_cap_20261016.json).

## Independence and inherited arithmetic

The checker does not import the new producer. It replays the assertions of
the frozen [independent 157-wall reconstruction](fresh_h20_release144_cap_20261016_review.md),
using that review's physical coordinates, algebraic parameter isolation,
and outward rational interval engine. It obtains the reconstructed geometry
by replacing only the predecessor's final packaging step in an in-memory
syntax tree; no predecessor file changes.

The new calculation uses old wall142 as its kernel seed, reverses the four
contact equations to `(156,151,145,95)`, and implements its own interval
Gaussian inverse. It checks both norm-two branches; only one has all required
dotted products and a strict common center. Projection identities give the
prescribed norm and contacts exactly. Interval overlap checks on these
identities are consistency checks; they are not used alone to infer equality.
All strict signs are certified by rational interval separation.

This is an independent reconstruction of the new cap and face, sharing the
stated physical arithmetic and old geometry. It is not a second independent
implementation of the original chain or its algebraic parameter isolation.

## New cap and ambient incidences

The reconstructed norm-two wall157 has 18 right contacts, angle pi/6 with95,
angle pi/3 with145, and 137 strictly dotted contacts. The retained center is
strict for all158 walls. Its projection to each wall is timelike; acuteness
makes every other inequality strict there, proving all158 facets actual.

The finite old neighbors of157 are precisely `S16+{95,145,151,156}`. Of the
twenty possible nineteen-root supports, eighteen contain the old dotted pair
95/145. The other two give the two new actual finite vertices constructed
below. This is a complete census on the new facet, not a complete ambient
finite-vertex census. Existing finite corners persist by acuteness and the
common strict interior point; the inherited total is at least184.

All five previously certified proper cusp vectors remain strictly inside the
new halfspace. The old compact triangle vertices are also retained strictly.

## Independent whole-face proof

The eighteen-root support `S16+{151,156}` is positive definite, so its
orthogonal complement has Lorentz signature `(2,1)`. The checker projects
the four sides95,142,145,157 directly to this three-space. The first three
projected normals have strictly negative Gram determinant and form a basis.

Each consecutive pair has positive definite Gram. Independent two-by-two
projection formulas produce four future timelike vertices. At each vertex,
all twenty active equations hold by the projection identities and all other
wall products are strictly negative. This checks all632 vertex-wall products.

For completeness, the checker also enumerates all six pairs of side lines in
the three-dimensional basis and tests both orientations of each intersection
ray. Exactly four rays satisfy the four side inequalities, and all four are
future timelike. The two opposite-side pairs have no feasible orientation.
The side cone is pointed because its normals span the three-space; its four
extreme rays therefore generate the cone. Its hyperbolic section is exactly
the compact quadrilateral. Every other retained wall is strict at the four
vertices, hence contains this entire quadrilateral. This proves the whole
actual face is compact, rather than only proving four isolated corners.

The adjacent intrinsic cosine squares are `0,1/9,1/3,3/4`. The face is not a
Coxeter polygon; the ambient wall contacts are Coxeter contacts.

## Remaining escape and replay

The exact old escape has positive norm

```
6401659894965851604205604657941 / 206232955413654050000000000000000000
```

and all158 products are strictly negative. Its future chart is inherited from
the replayed predecessor. Thus the full chamber still has infinite volume.
No finite-volume construction or new dimension is claimed.

Generation and fresh saved replay both passed. The saved status is
`PASS_INDEPENDENT_158_WALL_CAP_AND_COMPLETE_COMPACT_QUADRILATERAL`.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_compact_quad_cap_20261016_review.py --verify
```
