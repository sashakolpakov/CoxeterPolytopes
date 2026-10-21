# Five affine closures produce a second completed cusp

Five exact walls extend the [146-wall chamber](fresh_h20_h4_parallel_partner_20261016.md)
to **151 actual facets**. The ideal vertex `nu2` now has a bounded
19-dimensional Euclidean link. Its diagram is

```
affine A10 + affine A4 + affine A2 + affine A1 + affine A1 + affine A1.
```

The four previously certified proper cusps and every old finite vertex
persist. There are no new finite vertices. The strict spacelike escape
E5 remains, so this chamber still has infinite volume.

The [source](../scripts/fresh_h20_second_cusp_completion_20261016.py) and
[certificate](fresh_h20_second_cusp_completion_20261016.json) give all
exact recipes, rational interval bounds, height choices and diagram checks.
Generation and a fresh saved replay both pass.
This is a primary certificate; independent review of this new package is
pending. Frozen predecessors are unchanged.

## Exact construction

Retain the notation of the 146-wall certificate:

```
nu2 = C/sqrt(-C^2) - n144/sqrt(n144^2),
n145/sqrt(2) + n141/sqrt(n141^2) = z*nu2,  z>0.
```

The complete old incidence consists of the parallel pair `{141,145}` and
five mutually orthogonal spherical components:

| Component | Old walls | Added wall | Height | New squared norm |
|---|---|---:|---:|---:|
| A10 | 1 through 10 | 146 | 3 | 2 |
| A4 | 12 through 15 | 147 | 3 | 2 |
| A2 | 61, 87 | 148 | 3 | 2 |
| A1 | 142 | 149 | 5 | 309/200 |
| A1 | 143 | 150 | 6 | 19717/22500 |

For each row with old component roots `R_i`, define

```
p_i = -sum(R_i),
n_(146+i) = p_i + h_i*nu2.
```

All old component roots are perpendicular to `nu2`, and `nu2` is null.
For an Ar path with r at least 2, its root sum has squared norm 2,
product 1 with each endpoint and product 0 with its other roots.
Thus its added root has pi/3 contacts at the two endpoints. For each A1
we retain the original normal length: the product with its added root
is minus that squared norm, giving a parallel pair. No assertion that
walls 142 and 143 have norm 2 is used. The five added roots are mutually
orthogonal.

Every old wall outside the 20-node incidence has strictly negative
pairing with `nu2`. Consequently increasing a height decreases every
such new-old product. The certificate chooses the displayed positive
integers above strict interval bounds that simultaneously ensure:
all 126 nonincident old contacts are dotted, the common future center
remains strict, E5 remains strict, and the four old proper cusp rays
have negative pairing with every added root. The height choices are
canonical for this recipe, with no minimal-height claim.

All **740** new pair checks pass: **102 right, 6 pi/3, 2 parallel and
630 dotted**. Exact zero and finite-angle identities follow from the
component construction; intervals certify the strict comparisons.

## Bounded cusp link and actual facets

Normalize each incident root to squared norm 2 only for the diagram
calculation. Each Ar closure has a positive all-ones null relation.
Its physical root sum is the positive multiple `h_i*nu2`.
For the two A1 closures the normalized root sum is
`sqrt(2/q_i)*h_i*nu2`. The old parallel pair has normalized root sum
`sqrt(2)*z*nu2`. All six null relations therefore have the same future
null direction, with strictly positive multipliers.

The complete 25-node incident Gram matrix is the orthogonal sum of
these six affine diagrams. Exact rational elimination gives 19 positive
and 6 zero pivots. Their Euclidean spans have total dimension 19, and
each factor is a bounded affine Coxeter simplex. The horospherical link
is their bounded product. The old 20 incident normals already span
`nu2` perpendicular, so this is the same actual ideal vertex, now a
proper cusp. All other wall pairings with `nu2` are strictly negative.

The retained future center is strictly inside all 151 halfspaces.
For any wall n, its projection from this center onto n perpendicular is
timelike and has strictly negative pairing with every other wall:
the coefficient removed is negative and all mutual products are
nonpositive. This verifies all **151 facets are actual**.

The four old proper cusp points pair strictly negatively with every new
wall, preserving their links. Thus at least five proper cusps are now
certified. This is not a complete cusp census.

## Finite vertices and remaining escape

A finite vertex containing a new wall cannot contain any of the 126
old nonincident walls, because each such pair is dotted. Its incident
roots would therefore lie inside the 25-node cusp star. That star has
Gram rank 19, and cannot contain the 20 independent positive-definite
normals required at a finite H20 Coxeter vertex. Hence there are no new
finite vertices. Old finite vertices persist by the acute-wall
persistence lemma, leaving at least 171 known finite vertices; the
ambient finite census remains incomplete.

The exact escape inherited from the 146-wall certificate is

```
E5 = x - e/100 - (49/50)f - (3/10)omega17 - (1/25)omega31
     - (1/100000)sum(other16 support dual roots),
E5^2 = 31500061/5000000000 > 0.
```

All 151 pairings are strictly negative, and E5 lies in the future chart
of the common center. It therefore certifies infinite volume.

One further exact restriction follows from the completed link. An acute
additional spacelike normal r cannot pass through `nu2`: if `r.nu2=0`,
each positive affine relation forces r to be orthogonal to every root
in its component. Then r is orthogonal to all 25 incident roots, whose
physical span is `nu2` perpendicular, forcing r to be proportional to
the null vector `nu2`. Every additional coherent wall must instead
have strictly negative pairing with `nu2`. This is a consequence of
the completed cusp and does not exclude caps away from it.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_second_cusp_completion_20261016.py --verify
```
