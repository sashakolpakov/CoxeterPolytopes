# H20 facet gluing: exact seam choices, holonomy, and a six-tile angle repair

This construction search starts with the **non-Coxeter** H20 facets of
Borcherds' 210-facet H21 chamber. It constructs an exact six-tile,
convex, finite-volume H20 polytope with 580 facets. **It is not Coxeter.**
An independent exact certificate also excludes every attempt to obtain
a finite-volume Coxeter polytope by deleting walls from its 580-wall set.
It also rules out a specific complete gluing family: finite assemblies
using only the five mandatory A2 seams of each short facet. Gluing
additional seam types remains a distinct construction problem.

The executable is
[`new_dimension_h20_facet_gluing_20261013.py`](../scripts/new_dimension_h20_facet_gluing_20261013.py).
The full exact archive is
[`new_dimension_h20_facet_gluing_20261013.json`](../data/enumerations/new_dimension_h20_facet_gluing_20261013.json).
It includes both representative section Grams, all 210 defect-cover
certificates, complete seam-group generators and ambient lifts, ten
neighbor-placement matrices, 1,120 natural six-cycles, holonomy, and every
wall of the six-tile patch. Wall labels below are one-based ambient labels.

## Sections and unavoidable glued faces

The input is the frozen integral `U+D20(-1)` lattice and its 210 wall
vectors, reconstructed from Alexeev's transcription. The original chamber
is Borcherds' classical example, not a new result. See
[Borcherds, Example 5](https://math.berkeley.edu/~reb/papers/al/al.pdf) and
[Alexeev's matrix source](https://arxiv.org/abs/2209.09110).

The calculation uses the sign convention in which wall normals have
positive norm. For a facet normal `f`, it computes

`H_ij = G_ij - G_if G_fj / G_ff`.

The positive-norm projections give the intrinsic facet walls. Elliptic
pairs give genuine ridges. In particular the angle defects below concern
actual geometry, not just arbitrary pairs of supporting hyperplanes.

| Ambient facet type | Number | Intrinsic H20 walls | Bad ridges | Minimum number of glued faces needed |
|---|---:|---:|---:|---:|
| short, norm 2 | 42 | 209 | 510 | 5 |
| long, norm 4 | 168 | 82 | 480 | 40 |

Here a bad ridge has squared cosine outside `0,1/4,1/2,3/4`.
Every bad ridge must touch at least one glued facet of its tile, or its
angle survives unchanged on the final boundary. The script supplies an
exact vertex cover and an equally large disjoint-edge matching for every
one of the 210 defect graphs. No numerical optimization is needed to
verify the minimum.

Consequently any face-to-face convex assembly has at least six tiles;
an assembly containing a long tile has at least 41. Every assembly whose
tile-adjacency graph is a finite tree is impossible, since a leaf cannot
hide its defects with one glued facet. This statement does not cover all
successive whole-cluster doublings: such doublings can introduce cycles.

For the short facet at wall 1, the minimum cover consists of walls
`22,23,24,49,50`. These are its five A2 seams. Each has defect degree 104,
so every vertex cover of size below 104 must contain all five.

## Exact angle bookkeeping

Write

`gamma = acos(sqrt(2/3))`, `beta = pi/2-gamma`, `alpha = 2gamma`.

The non-Coxeter angles have squared cosines respectively `2/3`, `1/3`,
and `1/9`. The other angles are `pi/2`, `pi/3`, and `pi/4`.
Because `cos(2gamma)=1/3`, `gamma/pi` is irrational: if `2gamma` were a
rational multiple of pi, its doubled cosine would be a rational algebraic
integer, whereas `2/3` is not an integer.

Thus local angle repairs really are possible:

- `beta+gamma=pi/2`;
- `alpha+2beta=pi`;
- `2alpha+4beta=2pi`.

But two adjacent alpha sectors cannot lie at a convex Coxeter boundary,
even after adding further tile sectors. Cancelling their irrational
contribution requires at least four beta sectors, already giving a full
`2pi`. They must become an interior ridge.

## All twists of the mandatory A2 seam

The seam at ambient walls 1 and 22 is an H19 polytope with 208 facets.
The complete normalized weighted Gram automorphism group has order 2304.
The parent short-section seam stabilizer has order 1152 and index two.
The archive includes all five group generators, exact integral extensions
to the ambient lattice, and their action on the two four-element alpha
ridge sets. Four generators preserve the ambient pair separately; the
fifth interchanges it.

After quotienting by symmetries of the attached tile, there are exactly
two actual neighbors across each A2 seam:

1. **Self-reflection.** Reflect the short section across its own seam.
   If the ambient seam root is `j`, its intrinsic normal is proportional
   to `u=2j+f`, with `u^2=6`. The reflection matrix has denominator 3.
2. **Odd rolling.** Attach the adjacent ambient facet at `j`, transported
   to `f^perp` by reflection in the norm-two root `f+j`.

The first choice pairs alpha with alpha. The second pairs alpha with
beta. Consequently ordinary rolling does not exhaust the possible
intrinsic gluings; the self-reflection branch must be retained when
seeking constructions outside the ambient arithmetic reflection group.

All ten neighbors, two across each of five seams, were explicitly built.
The search compared exact oriented full wall sets for all two-step and
three-step placements. It found **no cycle closure of length 3, 4, 5, or
6** with different initial seams. There are 90 retained two-step
placements and 720 distinct three-step wall sets. This test is restricted
to A2 seam moves; the six-tile repair below uses two additional seams.

## Why using only A2 seams cannot close finitely

Every beta or gamma sector of a short tile has a noncore side. If only
the five A2 seams are glued, that side remains on the boundary. Hence
beta and gamma sectors cannot occur around an interior ridge in this
restricted family.

Self-reflection produces two adjacent alpha sectors, which cannot be
repaired at a convex boundary and cannot close an interior ridge using
only alpha sectors. So this family must use odd rolling everywhere.

The short-root odd diagram is the point-line incidence graph of
`PG(2,4)`: 42 vertices, 105 edges, and cycle rank 64. Standard rolling
develops its universal cover. The free part of the reflection centralizer
is the fundamental group of the odd diagram, as described by
[Allcock's reflection-centralizer theorem](https://arxiv.org/pdf/1106.3822).
It has rank 64 here, so no finite closed rolling assembly results.

There is also an exact concrete holonomy witness. The script enumerates
all 1,120 incidence hexagons and verifies their common affine A5 Gram
pattern. For the cycle `[2,5,4,51,8,47]`, it multiplies the six integral
transport reflections. Its restriction to H20 is nonidentity; its square
is a nontrivial unipotent matrix with `(T^2-I)^3=0` and
`rank(T^2-I)=2`. Thus this cycle has infinite-order holonomy, rather than
closing in H20. The matrix itself is archived.

This is an obstruction for the stated A2-only construction. It does not
exclude the use of additional seams, mixed short/long sections, or
arbitrary hyperplane sections.

## Concrete extra seams

The short section has the following five wall profiles. The seam counts
are computed for a representative of each profile; the archive retains
all labels and the full normalized Gram histograms. Equal profiles alone
are not used as a proof of isometry.

| Profile | Number | Irrational ridge sectors per wall | Representative seam |
|---|---:|---|---:|
| short orthogonal, no defects | 16 | none | H19/190 |
| short orthogonal, same point/line family | 20 | one beta | H19/196 |
| long orthogonal | 120 | two gamma | H19/81 |
| short A2 core | 5 | four alpha, 52 beta, 48 gamma | H19/208 |
| long B2 | 48 | five beta | H19/50 |

Every seam admits gluing to a reflected copy of itself. The five displayed
representatives are pairwise nonisometric, already distinguished by their
facet counts. Thus one cannot simply identify a beta seam with a gamma
seam on the basis of complementary angle values.

The H19/81 and H19/50 seams also admit **mixed short/long attachments**:
they are literally the same intersection of the two relevant ambient
facets, viewed from either parent. The 50-facet B2 intersection itself
has only Coxeter angles `pi/2,pi/3,pi/4`. These give concrete interfaces
for a larger mixed construction, although two-piece gluing does not
remove the parents' other defects.

## Six tiles that close the minimal interior angle

Take the short roots `f=1`, `j=22`, `k=23`, so that `f` is the central
node of their A3 diagram. Let `Q_f` and `Q_k` be their H20 facets and set

`D = reflection in 2j+f`,
`U = reflection in f+k`,
`S = reflection in n`, where `n=f+j+k`.

The vector `n` has norm two and is orthogonal to `f`. All six pieces lie
in `f^perp`:

`A=Q_f`, `DA`, `B=U Q_k`, `DB`, `C=S B`, `DC`.

Their cycle around the common H18 ridge is

`A -- DA -- DB -- DC -- C -- B -- A`.

The sectors are `alpha,alpha,beta,beta,beta,beta`, summing exactly to
`2pi`. All six full seam matches are verified from the exact projected
wall sets. Four interfaces have 208 facets; the additional interfaces
`B--C` and `DB--DC` have 196 facets and use the noncore short-short
orthogonal seam. The final core seam closes because the relevant normal
of `C` is again `2j+f`.

The independent [geometry audit](h20_gluing_independent_audit_20261013.md)
now certifies global convexity, disjoint tile interiors and finite volume.
There are exactly 580 true boundary facets. For every boundary normal,
the audit supplies an exact nonnegative rational combination of the
wall normals of each of the six tiles: 3,480 exact support certificates.
It also supplies exact separators for all 15 tile pairs. The six full
facet identifications, supporting halfspaces and disjoint interiors show
that the union equals their convex intersection. This bypasses the
A2-only holonomy obstruction by using the additional noncore seams.

The three-tile half `P=A union B union C` is convex too: it is the six-tile
union intersected with the negative halfspace of `u=2j+f`. This normal is
a wall of A and C, and is the sum of one B wall and twice another B wall.
The reflected triple lies on its other side.

These facts do not make either polytope Coxeter. An untouched ridge
between the original walls 24 and 49 still has angle alpha.

## Genuine boundary ridges and the cost of this repair

The [boundary script](../scripts/new_dimension_h20_six_patch_boundary_20261013.py)
and [exact archive](../data/enumerations/new_dimension_h20_six_patch_boundary_20261013.json)
recognize a global ridge only when its two-normal plane equals that of
a genuine codimension-two face of an original acute tile. Every global
ridge has such a witness: a generic point of its relative interior avoids
the finite collection of tile faces of codimension at least three, and
cannot lie inside a tile or in the relative interior of just one tile
facet. Every recognized plane contains exactly two global facet normals.
This matters because the glued polytope has obtuse angles; positive
definiteness of an arbitrary pair of its supporting normals alone does
not prove incidence. For the six-tile union, 6,474 such raw pairs are
excluded because they do not give actual ridges.

| Boundary angle | Two tiles A,DA | Three tiles A,B,C | Six tiles |
|---|---:|---:|---:|
| right | 8,412 | 10,155 | 16,824 |
| acute, squared cosine 1/2 | 3,744 | 4,536 | 8,136 |
| acute, squared cosine 1/4 | 3,072 | 3,948 | 6,376 |
| acute, squared cosine 1/3 | 416 | 515 | 824 |
| acute, squared cosine 1/9 | 60 | 16 | 114 |
| acute, squared cosine 2/3 | 384 | 504 | 816 |
| obtuse, squared cosine 1/3 | 0 | 7 | 14 |
| obtuse, squared cosine 1/4 | 0 | 4 | 8 |
| obtuse, squared cosine 1/9 | 52 | 0 | 103 |
| obtuse, squared cosine 49/81 | 4 | 0 | 7 |
| **All actual ridges** | **16,144** | **19,685** | **33,222** |
| **Facets** | **312** | **363** | **580** |
| **Non-Coxeter ridges** | **916** | **1,046** | **1,886** |
| **Exact minimum defect vertex cover** | **108** | **12** | **217** |

The minimum covers have explicit combinatorial certificates: disjoint
edges, supplemented for the six-tile case by four disjoint triangle
components, match the displayed cover sizes. Numerical optimization
only proposes the covers; these certificates prove the optima.

The angle with squared cosine 49/81 is the obtuse angle `2alpha`.
The initial double has four such ridges; the six-tile repair has seven.
Thus this isolated move repairs one compulsory interior ridge while
creating four others, a net increase of three. It does not prove that all
future assemblies grow forever: overlapping repairs might consume
several existing defects at once.

## Why full symmetric collars grow indefinitely

For the five core neighbors `j_i`, the ten roots
`n_ij=f+j_i+j_j` all have norm two and lie in `f^perp`.
Two have pairing zero when their index pairs meet, and pairing minus two
when those pairs are disjoint. The latter relation gives the Petersen
graph on the ten mirrors. The archive includes their complete Gram and
reflection matrices.

The product of the reflections in `n_(22,23)` and `n_(24,49)` has
`rank(T-I)=2`, `rank((T-I)^2)=1`, and `(T-I)^3=0`. It is a nonidentity
unipotent matrix, so its order is infinite. Closing a collar under all
ten reflections therefore cannot produce a finite symmetric tile
assembly. This does not exclude a selective finite patch or a finite
fundamental domain for the infinite reflection group.

## Exact codimension-two angle balance

The irrational angle contribution is balanced before any construction.
For one short tile, the 510 irrational ridges divide as follows:

- 10 ambient A3 faces with the short facet in the central position,
  contributing alpha;
- 20 ambient A3 faces with it at an endpoint, contributing beta;
- 240 relevant ambient B3 faces in the short central position,
  contributing beta;
- 240 relevant ambient B3 faces at the short endpoint, contributing gamma.

The full chamber has 420 such A3 faces and 10,080 relevant B3 faces.
At each identical ambient H18 face, the short-facet incidences are
`alpha,beta,beta` for A3 and `beta,gamma` for B3. The irrational
coefficients cancel exactly, with no comparison of unknown volumes:
`2-1-1=0` and `-1+1=0` respectively. All ambient triple labels are archived.

Consequently the sum of the codimension-two Dehn components of all
42 congruent short facets is zero. In the rational, torsion-free target
this gives zero for each short facet, including the projection using
H18 ridge volumes. This is **only the codimension-two component**;
no claim is made about other generalized Dehn invariants or the existence
of a compatible finite assembly.

## Exhausted facet-deletion construction

The [selection script](../scripts/new_dimension_h20_compatible_subsets_20261013.py)
uses the complete pair-conflict graph on all 580 outward normals.
A pair conflicts if its normal pairing is positive, or if its hyperplanes
meet and its squared cosine is outside `0,1/4,1/2,3/4`. There are 7,072
conflict edges. This is the complete Gram condition, not just the
actual-ridge census of the original six-tile union.

The bounded round used 32 fixed greedy seeds, extensions of the original
204 compatible noncore walls, and six MILP objectives with 15-second
limits, including three symmetry-constrained objectives. It produced
17 distinct maximal compatible subsets of 360--363 walls, plus a
168-wall subset preserved by D. The latter is retained before maximal
completion, since adding single walls can destroy that symmetry.
The [candidate archive](../data/enumerations/new_dimension_h20_compatible_subsets_20261013.json)
records their exact indices, compatibility and symmetry checks.

Thirteen of the eighteen candidates have an explicit nonintegral
length-two Vinberg cyclic product `4(r,s)^2/(r^2 s^2)`; for example
candidate 001 contains walls 43 and 387, with norms 4 and 18, pairing -14,
and cyclic product `98/9`. This is an intrinsic obstruction to arithmetic
integrality, unlike failure to preserve the particular original lattice.
It would be relevant if these chambers had finite volume.

All eighteen instead have strict exact spacelike escape certificates in
the [independent escape archive](../data/enumerations/h20_gluing_independent_escape_20261013.json).
Each certificate is an integer vector in `f^perp`, with positive square
and strictly negative pairing against every retained normal. Together
with a common strict timelike center, it gives an open set of ideal
directions in the chamber and therefore proves infinite volume.

More strongly, **every possible subset of the 580-wall dictionary is
excluded**, using just two deletion escapes. In the zero-based canonical
boundary ordering, walls 0 and 74 have Gram

`[[4,-4],[-4,6]]`,

so their actual ridge has squared cosine `2/3`. One exact spacelike
vector satisfies all 579 other wall inequalities strictly when wall 0 is
removed; another does so when wall 74 is removed. Thus every finite-volume
subset must retain both walls. Their genuine original ridge survives
every deletion retaining them and is not a Coxeter ridge.

The independent
[deletion-obstruction report](h20_gluing_independent_deletion_obstruction_20261013.md)
and [primary certificate](../data/enumerations/h20_gluing_independent_deletion_obstruction_20261013.json)
also give an explicit timelike point lying on exactly these two walls,
strictly inside all other halfspaces. The auditor simplifies the first
escape to `c+10r_0`, of square 22, and independently verifies
[all 580 single-wall deletion escapes](../data/enumerations/h20_gluing_independent_all_deletions_20261013.json).
Thus every wall of this particular finite-volume patch is individually
necessary for finite volume.

The [minimal integer replay](../scripts/new_dimension_h20_counterexample_sat_20261013.py)
imports no other project verifier. Its
[certificate](../data/enumerations/new_dimension_h20_counterexample_sat_20261013.json)
contains the two vectors, their positive squares
`21024997679999942` and `6718042333490830`, and the exact Gram.
Both vectors pair at most -2 with every other wall and positively with
the omitted wall. The necessary Boolean clauses are simply
`x_0`, `x_74`, and `not x_0 or not x_74`.
An exhaustive four-case truth table and an exported Z3 proof both confirm
UNSAT. A forged timelike replacement witness is rejected on replay.

This closes the finite-dictionary deletion route, not the H20 problem.
A successful continuation needs genuinely new supporting planes or a
different assembly, for example through additional noncore or mixed
short/long seams. The six-tile move and its exact geometry remain usable
building blocks, but no new finite-volume Coxeter polytope is claimed.

```sh
python3 scripts/new_dimension_h20_facet_gluing_20261013.py
python3 scripts/new_dimension_h20_six_patch_boundary_20261013.py
python3 scripts/new_dimension_h20_counterexample_sat_20261013.py --verify
```
