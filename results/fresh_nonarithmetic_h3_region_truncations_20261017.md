# Three complete single-wall truncation obstructions in the retained H3 face

For the actual face supported by `S6=[0,1,2,3,5,6]` of the retained 55-wall
chamber, a single added Coxeter wall cannot terminate all the selected open
edge ends around any of these three patches:

* triangle `[7,30,49]`;
* quadrilateral `[8,30,44,52]`;
* quadrilateral `[30,34,47,52]`.

This covers all finite Coxeter labels and arbitrary coefficient fields. It
does not exclude a coordinated repair by several new facets. No wall is
added, and no finite-volume claim is made.

The [source](../scripts/fresh_nonarithmetic_h3_region_truncations_20261017.py)
and [certificate](fresh_nonarithmetic_h3_region_truncations_20261017.json)
reuse the exact physical normals and outward rational interval engine. The
input [complete endpoint domains](fresh_nonarithmetic_h3_endpoint_domains_20261017.md)
provide the all-label finite-endpoint necessity; no sampled angle cutoff is
used for an exclusion here.

## The five patches and complete profile intersections

The complete old S6 cone has thirteen finite and seven exterior vertices.
Three exterior vertices give triangular patches. The other four are joined
in pairs across the nongeodesic raw cone edges `[44,52]` and `[47,52]`, giving
two quadrilateral patches. Their listed geodesic boundary ends account for
all seventeen old open edge ends, with multiplicity checked against the
complete cone inventory.

| Old side walls | Boundary edges | Finite profiles | High-label tails |
|---|---|---:|---:|
| 7,30,49 | 7–30, 7–49, 30–49 | 2 | 0 |
| 7,33,49 | 7–33, 7–49, 33–49 | 5 | 0 |
| 30,33,49 | 30–33, 30–49, 33–49 | 8 | 0 |
| 8,30,44,52 | 8–44, 8–52, 30–44, 30–52 | 84 | 4 |
| 30,34,47,52 | 30–47, 30–52, 34–47, 34–52 | 84 | 4 |

A single wall terminating all these selected ends must give a positive
definite rank-nine Gram at each endpoint. The frozen endpoint theorem
therefore applies at each boundary edge. Intersecting its complete profile
lists gives the table. In the certificate, code `6` means the entire tail
`c=2 cos(pi/m)`, integer `m>=6`, relaxed to `sqrt(3)<=c<2`; it never means
only the label 6.

Only a complete necessary-profile intersection is recorded for the second
and third triangles. Those two patches are not excluded or certified
closed by this package.

## Exact rejection of the three patches

All physical normals are positively normalized to square 2. For an
independent frame `F` of prescribed wall contacts, write `G=Gram(F)` and
`p=F G^{-1}c`, where `c` is the negative contact vector.

For `[7,30,49]`, the nine-element frame `S6+[7,30,49]` leaves a
one-dimensional positive kernel. Projecting the old center off the frame
gives a nonzero vector `k`; all possible normalized roots are exactly
`p +/- sqrt((2-p²)/k²) k`. The two permitted profiles both have positive
discriminant, hence exactly four real branches. Two have positive product
with old wall 4. The other two have products with old wall 50 strictly in
`(-1,0)` and `(-sqrt(2),-1)`, respectively. These are the gaps between
consecutive Coxeter labels 2/3 and 3/4. No other finite label can fill either
gap, and neither product is dotted or parallel.

Each quadrilateral frame has ten independent normals. Every finite contact
profile therefore specifies a unique vector. All 168 such vectors have
square different from 2, certified by disjoint outward rational intervals.
The smallest absolute norm discrepancy is greater than 0.0043, so none is
an undecided near-zero computation.

Each of the eight tail recipes gives an affine vector `p+c*d`. Its norm
minus 2 is `a*c²+b*c+e`. The certificate proves `a>0`, `b>=0`, and a strictly
positive value at `c=sqrt(3)`. Thus its norm is above 2 throughout the whole
continuous interval `sqrt(3)<=c<2`. This also covers every integer label in
the tail.

All branches are rejected before imposing a retained particular center or
mutual new-wall conditions. The result therefore does not rely on those
extra restrictions. In a coherent finite-volume completion the retained
pi/5 support makes this H3 face compact, so its open edges really do need
finite endpoints. The obstruction established here is specifically to one
wall supplying all selected endpoints of a patch; different endpoint walls
remain possible.

## Validation and next scope

Generation and fresh saved replay pass. The checker verifies the profile
counts, all seventeen boundary-end incidences, frame invertibility, four
triangle branches, 168 nonzero norm discrepancies, and eight continuous
tail inequalities. It shares the frozen physical interval arithmetic; no
independent arithmetic reconstruction is claimed.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_nonarithmetic_h3_region_truncations_20261017.py --verify
```

The unresolved constructive step is to connect multiple endpoint walls
across the three failed patches while retaining every old native Coxeter
gate. The five-patch inventory supplies finite targets for that search; an
increase in facet or finite-vertex counts alone would not close the H3 face.
