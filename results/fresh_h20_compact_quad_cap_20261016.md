# A compact quadrilateral face in the 158-wall H20 chamber

The [157-wall chamber](fresh_h20_release144_cap_20261016.md) extends to
**158 actual facets and at least 184 finite vertices**. Its previously open
two-face with support `S16+{151,156}` is now a **whole compact quadrilateral**.
The old compact triangle and all five proper cusps survive. A strict spacelike
escape still proves that the full chamber has **infinite volume**.

The [exact construction](../scripts/fresh_h20_compact_quad_cap_20261016.py)
and [certificate](fresh_h20_compact_quad_cap_20261016.json) passed generation
and a fresh saved replay. The source replays the frozen independent
157-wall geometry through all its assertions before constructing the new
root. It returns that geometry by replacing only the final report-packaging
part of an in-memory syntax tree; no predecessor file is modified.

The [independent review](fresh_h20_compact_quad_cap_20261016_review.md)
also passes generation and fresh saved replay. It uses a different kernel
seed and reversed contact frame, then projects directly to the Lorentz
three-space of the face and exhausts all six side pairs. Exactly four
feasible future timelike rays remain, independently confirming compactness.

## The new wall

The new norm-two root, numbered157, is right to `S16+{151,156}`, has
angle pi/6 with95 and pi/3 with145, and is strictly dotted to every other
retained wall. These are **18 right contacts, two finite nonright contacts,
and 137 dotted contacts**, all checked with rational interval arithmetic.

Here is a coordinate-independent definition. Let

```
U = (n95,n145,n151,n156),
a = n95.n145 < -2,
G = [[2,a,0,0], [a,2,0,-1], [0,0,2,0], [0,-1,0,2]],
t = (-sqrt(3),-1,0,0).
```

Let `Hres` be the common future center projected perpendicular to S16. Put

```
p = sum_i (G^-1 t)_i U_i,
k = Hres - proj_U(Hres),
r157 = p + sqrt((2-p^2)/k^2) k.
```

The checks give `p^2≈-4.4053296511`, `k^2≈245.7654265236`, and a positive
radicand. Both signs of the square root are tested; the displayed branch
is the unique one passing the common-center and all old contact conditions.
The projection identities prove the norm and prescribed contacts exactly.
Strict interval inequalities prove every dotted relation.

The common center stays strict. Its orthogonal projection to each wall
is timelike and strict on every other wall, proving all158 facets actual.
Five strict negative cusp products preserve the already completed links.

## Why the whole face is compact

The four sides of the face are, in cyclic order, **95,142,145,157**.
All four intersections of consecutive sides are finite vertices:

| Side pair | Status |
|---|---|
| 95,142 | Previously certified finite vertex |
| 142,145 | Previously certified finite vertex |
| 145,157 | New finite vertex |
| 157,95 | New finite vertex |

Each full vertex support is S16 together with151,156 and the indicated
pair. The checker constructs each point by orthogonal projection, proves
the incident Gram positive definite, and verifies all158 wall products.
Every nonincident product is strictly negative.

In the projective face plane, the four side lines therefore support the
convex quadrilateral formed by these four future timelike vertices. Their
four inward halfplanes intersect in precisely that quadrilateral. Every
remaining wall contains all four vertices and hence their convex hull.
This certifies the entire face, not merely four isolated finite corners.

The only finite old neighbors of157 are `S16+{95,145,151,156}`. Of their
nineteen-root subsets, only those omitting95 or145 are spherical; the other
eighteen contain the dotted pair95/145. Thus exactly two finite vertices
are added, giving the inherited lower bound184.

The face itself is not a Coxeter quadrilateral. Its intrinsic angle cosines
are `0,1/3,1/sqrt(3),sqrt(3)/2`. Coxeter dihedral angles of the ambient
chamber do not force every lower-dimensional face to be a Coxeter polytope.

## What remains open

The previous exact escape E7 has new-wall product
approximately **-0.156378643898990**. Its positive rational norm and all158
strict negative wall products are replayed, so the ambient chamber still
has infinite volume. No new finite-volume example or dimension is claimed.

The [bounded numerical scout](../scripts/fresh_h20_pentagon_closure_scout_20261016.py)
and [data](fresh_h20_pentagon_closure_scout_20261016.json) tested31,460
prescribed four-contact patterns with labels2 through12 and both real
branches. It found4,914 distinct numerical candidates, including24 that
make this face compact and nine with near-sphere corners. Candidate3309,
with labels `(2,2,6,3)` on `(151,156,95,145)`, was selected for the exact
construction above. The finite search is not a completeness claim over
all Coxeter extensions. Generation and fresh numerical replay both pass.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_compact_quad_cap_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_pentagon_closure_scout_20261016.py --verify
```
