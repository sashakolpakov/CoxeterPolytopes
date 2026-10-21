# Independent audit of the H20 section-gluing construction

**Outcome:** the five-core-seam obstruction is valid under its stated
assembly restrictions. Allowing additional seams escapes that obstruction:
the proposed six-tile patch is an actual convex, finite-volume H20 polytope
with 580 facets. It is **not Coxeter**. An exact interior point witnesses a
surviving non-Coxeter ridge; the complete boundary census also finds obtuse
dihedral angles.

The [independent checker](../scripts/h20_gluing_independent_audit_20261013.py)
reconstructs the geometry from the frozen Borcherds lattice, checks all
210 section defect covers, recomputes the complete A2 seam group, reverses
the holonomy product convention as a cross-check, and certifies the patch
using rational cone identities. The [certificate](../data/enumerations/h20_gluing_independent_audit_20261013.json)
contains all 580 boundary normals and their support proofs.

## The projected walls and defect edges are actual faces

The frozen Gram has 42 square-two roots and 168 square-four roots, all
off-diagonal pairings nonpositive. The sum `h` of all 210 roots has square
`-137214` and strictly negative pairing with every root. Exact leading
principal minors of the ambient form have 21 positive values followed by
`-4`, verifying Lorentz signature `(21,1)`.

Projection of `h` onto an old facet is timelike and strict for every other
old inequality. The projected wall Gram remains nonpositive off the
diagonal. Every positive-norm projected wall is a true facet: project this
strict point onto that wall; all other inequalities stay strict. Projected
normals of nonpositive square have no hyperplane meeting the open
hyperbolic section and impose no additional facet.

The same argument applies to a pair of projected walls with a positive
definite two-by-two Gram. Its inverse has nonnegative entries, so the
orthogonal projection of the strict point onto their common subspace is
timelike and strictly satisfies every other inequality. Thus the pair
really is a ridge. This justifies the defect graph geometrically; it is
not merely a list of potentially intersecting hyperplanes.

The independent replay checks the explicit matching and covering witnesses
for all 210 sections. Short sections have 209 facets and minimum defect
cover 5; long sections have 82 facets and minimum cover 40. The five short
core facets each have defect degree 104. Therefore their set is the unique
cover of size five: omitting any one would require all 104 of its defect
neighbors in the cover. This is a statement about which facets must be
made internal in this gluing setup, not a lower bound on arbitrary H20
Coxeter polytopes.

All these squared cosines are rational. The only rational squared cosines
of Coxeter angles `pi/m`, `m>=2`, are `0,1/4,1/2,3/4`, so the defect gate
does not inadvertently omit a larger finite label.

## Two A2 seam placements and the restricted obstruction

The common A2 seam of short facets 1 and 22 has 208 actual facets. An
independent normalized-Gram graph computation gives its complete metric
automorphism group order 2304. The archived generators generate all 2304
permutations; their integral ambient lifts preserve the form and all 210
original roots. The four remaining core facets on this seam have an orbit
of size two, so their stabilizer has order 1152. Core preservation is
intrinsic, since these facets are distinguished by their defect pattern.
The lifts realize this subgroup as parent-section seam symmetries.

Consequently the two geometric neighbor placements, modulo parent-section
symmetries, are reflection in the seam itself and odd rolling. This uses
the entire seam isometry group, not just a selected list of twists.

Write

\[
 \gamma=\arccos\sqrt{2/3},\qquad
 \beta=\pi/2-\gamma,\qquad
 \alpha=2\gamma=\arccos(1/3).
\]

Then `alpha+2 beta=pi`. In a construction where **every tile glues its
five core facets and all its noncore facets remain exposed**, a ridge
sector of angle alpha has two core sides. A self-reflection seam joins
alpha to alpha. Rolling joins alpha to beta, and that beta sector ends at
an exposed noncore side. Thus a finite chain of `a` alpha sectors has two
beta ends and angle

\[
 a\alpha+2\beta=\pi+(a-1)\alpha.
\]

Convexity permits only `a=1`. A cycle made solely of alpha sectors cannot
close to `2 pi`: `alpha/pi` is irrational, since `cos(alpha)=1/3` and the
rational-cosine criterion excludes this value at a rational multiple of
pi. Hence a finite convex Coxeter assembly with these restrictions must
use rolling, not self-reflection, across each core seam.

The connected short-root odd graph has 42 vertices and 105 edges, hence
cycle rank 64. The checker independently develops the representative
six-cycle with **right** products of the root reflections, obtaining the
inverse of the archived left-product matrix. Its square is `I+N`, with
`N^3=0`, `rank(N)=2`, and `N` nonzero on the H20 section. It therefore has
infinite order there.

If a finite developed assembly closed every core seam using rolling,
following this same graph loop repeatedly would visit only finitely many
tiles. Some positive holonomy power would stabilize a tile. Its finite
facet permutation group has finite order, whereas the computed unipotent
holonomy has infinite order. This is a contradiction.

**Scope:** this excludes finite developed, full-facet assemblies of the
specified short sections in which only the five core seam types are
internal and every such seam is closed. It does not exclude additional
noncore seams, mixed section types, cut tiles, more general gluing
identifications, or H20 Coxeter polytopes. The six-tile patch below uses
exactly the omitted freedom.

## Exact six-tile patch

Use one-based roots `f=1`, `j=22`, `k=23`; their Gram is

\[
 \begin{pmatrix}2&-1&-1\\-1&2&0\\-1&0&2\end{pmatrix}.
\]

Let `Q_f` and `Q_k` denote the intrinsic facet sections. Put

\[
 u=2j+f,\qquad v=f+k,\qquad n=f+j+k.
\]

Their squares are respectively 6, 2, 2. Let `D`, `U`, and `S` be reflection
in these three normals. Inside `f^perp`, the six tiles are

\[
 A=Q_f,\quad DA=DQ_f,\quad B=UQ_k,\quad DB=DUQ_k,
 \quad C=SUQ_k,\quad DC=DSUQ_k.
\]

The exact checker compares complete projected facet-wall sets on each
shared plane. It confirms precisely the six intended full-seam matches:

| tile pair | original wall labels in the two tiles |
|---|---|
| A, DA | 22, 22 |
| A, B | 23, 1 |
| DA, DB | 23, 1 |
| B, C | 22, 22 |
| DB, DC | 22, 22 |
| C, DC | 1, 1 |

The B–C and DB–DC matches use the additional noncore seams. Two further
opposite-normal coincidences, A–DC and DA–C, have different facet sections;
they are explicitly rejected as full-seam matches. They involve facets
already internal through the six matches above.

### Exact nonoverlap and support certificates

All six tiles have strict timelike interior points. For each of their
15 pairs, the certificate chooses a facet normal of one tile whose
negative is a nonnegative rational combination of facet normals of the
other. The tile interiors therefore lie on opposite sides of that plane
and cannot overlap.

Cancel the twelve facet occurrences in the six full-seam matches and
merge positively proportional remaining normals. Exactly 580 distinct
oriented boundary normals remain. For every such normal `r` and every
tile `T`, the certificate gives an identity

\[
 r=\sum_i c_i r_i^{(T)},\qquad c_i\geq0.
\]

There are 3,480 identities. Every coefficient is a rational number; the
checker verifies every coordinate equality exactly. Numerical linear
programming was used only to find coefficients when generating this
certificate. Replay does not invoke the LP solver.

It follows that every one of the 580 halfspaces contains **all** six
tiles. Every such plane is a true facet of their union: it contains a
full, 19-dimensional facet of at least one owner tile, and is supporting
for the whole union. Therefore none of these 580 planes is merely a
redundant inequality or a tangency outside the hyperbolic slice.

### Why the union equals the convex halfspace intersection

Let `K` be the intersection of these 580 halfspaces inside the open Klein
ball and `P` the six-tile union. The support certificates give `P subset K`.
Suppose a point of `K` lay outside `P`. Choose a generic segment from an
interior tile point to such a point. A small perturbation avoids the
finitely many faces of codimension at least two. Its first exit from the
closed finite union must then cross the relative interior of a tile
facet. A fully matched seam is locally interior to the union, since the
two tiles lie on opposite sides and have identical full facet sections.
Thus an exit must cross one of the remaining facets. But its supporting
halfspace contains both endpoints of the segment, which cannot cross it
outward. This is a contradiction.

Hence `P=K`, proving global convexity. Finite volume follows from the
finite volume of the six Borcherds facet sections; disjoint interiors
also make the volumes additive. No assumption of final acuteness enters
this argument.

The three-tile union `A+B+C` is a convex half of this patch. The normal
`u` supports A and C directly. In B it is exactly the first wall normal
plus twice the wall normal with original label 22. The reflected tiles
are in the opposite halfspace. These additional exact identities prove
`A+B+C = P intersect {x:(x,u)<=0}`.

## Why the patch is not Coxeter

The two base walls with original labels 24 and 49 remain supporting
facets. Their primitive normals `r,s` have Gram

\[
 \begin{pmatrix}6&-2\\-2&6\end{pmatrix},
 \qquad \cos^2\theta=1/9.
\]

Let `h` be the archived strict center of tile A and write
`a=(h,r)`, `b=(h,s)`. The checker verifies the explicit point

\[
 p=32h-(6a+2b)r-(2a+6b)s
\]

is timelike, lies on these two walls, and is strict for all 207 other
walls of A. The full original tile ridge is therefore contained in both
global supporting planes, proving this is an actual surviving H18 ridge.
Its angle is alpha, which is not a Coxeter angle.

The companion [boundary census](../data/enumerations/new_dimension_h20_six_patch_boundary_20261013.json)
uses actual tile-ridge normal planes to recognize final ridges. I reviewed
its exact plane-key construction and the recognition argument. A global
supporting pair is an actual ridge when its normal plane agrees with that
of a genuine tile ridge. Conversely, a generic point of a global ridge
lies on a tile codimension-two face: finitely many higher-codimension
faces cannot cover the ridge, and the normal cone at an interior tile
facet is only one-dimensional. This gives completeness without assuming
the final polytope is acute.

The census finds 19,685 actual ridges in the 363-facet three-tile half and
33,222 actual ridges in the 580-facet six-tile patch. For the latter,
6,474 positive-definite global normal pairs are **not** realized ridges.
Counting intersecting hyperplanes alone would therefore be incorrect.
The six-tile patch has 132 actual obtuse ridges, including seven with
squared cosine `49/81`, angle `2 alpha`. These are additional failures
of the Coxeter condition, not evidence against more selective repairs.

## Replay and scope

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/h20_gluing_independent_audit_20261013.py
python3 -m py_compile scripts/h20_gluing_independent_audit_20261013.py
```

To regenerate all rational cone witnesses from LP selections, use
`--rebuild`. Every selected witness is still checked exactly before it
is saved. The default replay consumes only the saved rational identities.

This audit certifies a concrete convex non-Coxeter H20 assembly and a
restricted obstruction. It does not certify a new Coxeter example, a
new highest dimension, or worldwide novelty of the auxiliary construction.
