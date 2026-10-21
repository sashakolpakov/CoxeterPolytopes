# Polar cuts, new root blocks, and H20 gluing — 13 October 2026

**No new compact or finite-volume Coxeter dimension has been established.**
This continues the [fresh dimension search](fresh_dimension_search_2026_10_13.md).
The targets remain compact H9 and finite-volume H20. The current calculations
resolve a gap in the earlier polar-cut screen, test different positive root
blocks, identify three disguised copies of previously searched H9 lattices,
and construct a six-piece H20 angle repair.

The constructive gains are a 372-facet H9 Coxeter cone obtained using
targeted arithmetic cuts, and an exactly verified convex 580-facet H20
assembly. The former still has infinite volume; the latter still has
non-Coxeter angles. Two exact witnesses rule out every attempt to fix the
H20 assembly merely by deleting some of its 580 facets.

## Removing invisible polar cuts

The previous numerical screen applied its angle test before checking whether
each proposed polar truncation actually survived as a facet. This matters:
many proposed cuts are redundant inside hyperbolic space.

For a simple hyperideal vertex of an acute connected Lorentzian old-wall
diagram, its polar cut is visible as a hyperbolic facet exactly when its
incident old-wall diagram is connected. The
[polar-reduction report](polar_reduction_h9_20261013.md) gives the proof,
including why all invisible cuts can be removed simultaneously and why the
old facets remain visible.

Replaying the same bounded family removes 1,036,376 invisible cuts and retains
268,980. All 67,344 noncusped systems still fail a retained polar/old-wall
angle test. The 2,465 degenerate and 43,290 cusp-bearing skeletons remain
outside this test. The geometry criterion is proved exactly, but the large
screen still uses floating-point rank, signature, and angle decisions; its
negative result is not a general nonexistence theorem.

An exact algebraic example verifies the disconnected case. An independent
convex-optimization check of its 24 hyperideal cuts agrees with the criterion:
five visible, nineteen invisible. The report links the reproducible scripts
and archives.

## E7+A2 hosts, with the missing G2 reflection included

Over `O=Z[phi]`, `phi^2=phi+1`, the hosts
`<-alpha>+E7+A2`, for `alpha=phi,2phi`, have a complete height-zero chamber
of type **E7+G2**. The extra A2-vector `(1,-1)` has norm six and must be
included; using only the A2 Weyl chamber would omit a reflection.

| alpha | Complete candidate norm alphabet | Height | Walls | Exact escape norm |
|---|---|---:|---:|---|
| phi | 1,2,3,4,6,12 | 30 | 29 | `(44-27phi)/5` |
| 2phi | 2,4,6,12 | 30 | 30 | `(1514-774phi)/1083` |

Both escape norms are strictly positive. The saved vectors satisfy every
wall inequality, and each cone has a strict timelike interior. Thus these
finite prefixes have infinite volume. Neither full reflection group has
been decided. A standalone rational-arithmetic verifier checks the saved
certificates; six structural and certificate tests pass. See the
[root-block report](fresh_root_block_h9_report.md).

These are new search presentations in this repository, not established new
reflection groups. In particular the second host shares a determinant square
class with an earlier E8 host; an overclosure relation remains possible.

## Field-native H4 blocks and exact duplicate detection

The H4 Gram matrix has diagonal two and consecutive pairings `-phi,-1,-1`.
Its determinant `5-3phi=phi^-4` is a unit. It is positive definite at both
real embeddings. We implemented exact dominant-weight enumeration for sums
of these blocks, with an optional scalar coordinate.

Two calibrations close: `<-phi>+H4` gives the known compact H4 simplex
`[5,3,3,4]`, and `<-2phi>+H4+H4` gives the known compact eleven-facet
Bugaenko H8 chamber, with 41 finite vertices. An explicit integral isometry,
with integral inverse and exact transport of all eleven walls, proves the
second identification.

Extending and conjugating that isometry proves that **all three tested H9
H4 hosts are integral copies of earlier E8 hosts**. They must not be counted
as new lattices or new reflection groups. The matrices and exact Gram
congruences are recorded by the
[isometry script](../scripts/fresh_h4_blocks_isometry_20261013.py).

| Host | Height | Walls | Exact escape norm |
|---|---:|---:|---|
| H9: `<-phi>+H4^2+<1>` | 30 | 22 | `(246-152phi)/5` |
| H9: `<-2phi>+H4^2+<1>` | 30 | 23 | `(466-287phi)/5` |
| H9: `<-sqrt5>+H4^2+<1>` | 30 | 21 | `56-69phi/2` |
| H12: `<-2phi>+H4^3` | 30 | 36 | `(7576-4682phi)/5` |
| H13: `<-2phi>+H4^3+<1>` | 30 | 65 | `290-179phi` |

Every displayed norm is positive. All five finite prefixes therefore remain
open. The H12 and H13 presentations also admit an E8+H4 decomposition using
the H8 isometry; no worldwide novelty is claimed for those hosts.

The [independent audit](fresh_h4_blocks_independent_report_20261013.md) uses
standalone Fraction arithmetic. It verifies all five escapes, root
integrality, pairwise acuteness, chronology, strict interiors, active ranks,
and all wall inequalities. It independently enumerates the control
coordinates and height-zero roots, checks the 120-root H4 orbit, and
replays the two compact calibrations and all four lattice identifications.
Thirteen tests pass, including rejection of forged certificate data.
Positive-height exhaustion is audited in the source and calibrated; it was
not independently reimplemented.

The duplicate identification nevertheless gives a useful search method:
transport walls from the different Vinberg basepoints into a common E8
model. After identifying duplicate normals, all cross-pairings are
compatible. Adding those walls and their
images under the transported H4-block swap enlarges the old `2phi` prefix
from 184 to **364 walls**, and the `sqrt5` prefix from 43 to **85 walls**.
Both old exact escape vectors survive those added inequalities. These are
larger open Coxeter cones; the search gain comes from changing basepoints
rather than increasing a uniform height cutoff. Every finite angle is
checked explicitly. A common interior and nonpositive pairings alone
would not establish that the two complete lattice chambers coincide, nor
would they exclude an acute but non-Coxeter angle such as `2pi/5`.

An escape-directed search subsequently finds a compatible cut for the
`2phi` host. Saturating the spatial direction of its escape reduces the
root equation to

`61 t^2 - 2phi x^2 = 1`.

The solution `t=11+18phi`, `x=48+78phi` gives the norm-one root

`r=(x, t*(14,28,40,32,24,16,8,20,sqrt5))`.

It pairs nonpositively with all 364 walls, negatively with the common
interior, and positively with the former escape, with pairing
`(6+5phi)/7`. Its transported block-swap image is compatible as well.
Together they give 366 walls and remove both members of that escape
orbit. The displayed root's original-basepoint height is `27144+43920phi`, roughly
98,000, far beyond the earlier uniform cutoff of 300. Three further
bounded adaptive rounds add another six compatible mirrors. The resulting
**372-wall Coxeter cone still has an exact escape**, with norm
`(424-262phi)/175 > 0`. Every finite angle in the full wall list is
verified; the construction does not assert chronological completeness
or that this is the full lattice reflection chamber.
The [multi-basepoint and targeted-cut report](fresh_multibasepoint_h9_report_20261013.md)
records all four cut orbits, the failed bounded `sqrt5` slice search,
standalone certificate verifiers, and six passing regression tests. One
regression explicitly rejects an integral H4 root pair at angle `2pi/5`.

## H20: extra seams permit a six-piece repair

The two types of H20 facets of Borcherds' H21 chamber have 209 and 82 walls.
Neither is Coxeter. Their forbidden angles can be written

`gamma=acos(sqrt(2/3)), beta=pi/2-gamma, alpha=2gamma`.

The identities `beta+gamma=pi/2`, `alpha+2beta=pi`, and
`2alpha+4beta=2pi` give exact ways to cancel the irrational angle parts.
The [gluing report](new_dimension_h20_facet_gluing_20261013.md) records full
seam isometries, a restricted obstruction for using only five A2 seams per
short tile, and the explicit six-tile construction that escapes that
restriction by using two further seams.

For ambient short roots `f,j,k` forming A3 with central root `f`, let
`D` reflect in `2j+f`, `U` in `f+k`, and `S` in `f+j+k`. The pieces

`Q_f, DQ_f, UQ_k, DUQ_k, SUQ_k, DSUQ_k`

all lie in `f^perp`. Six exact full-facet identifications close an interior
ridge with angles `alpha,alpha,beta,beta,beta,beta`. This is a concrete
construction with archived wall vectors, not just a numerical angle match.
The [independent global audit](h20_gluing_independent_audit_20261013.md)
certifies a **convex finite-volume H20
polytope with 580 facets**: 3,480 exact nonnegative-combination certificates
show that each boundary plane supports all six pieces, and exact separating
certificates establish nonoverlap of every pair of pieces.

Its outer boundary is **not Coxeter**. The
[actual-ridge census](../data/enumerations/new_dimension_h20_six_patch_boundary_20261013.json)
finds 33,222 ridges, including 1,754 irrational acute defects and 132 obtuse
angles. Seven of the latter have angle `2alpha`, which cannot be repaired
into a convex rational-angle boundary by adding the available sectors:
cancelling their irrational contribution already requires a full `2pi`.
They must become interior ridges in any such extension. Thus repairing one
local cycle leaves substantial new global obligations. The census uses
actual tile-ridge incidence; it rejects 6,474 positive-definite pairs of
boundary hyperplanes that are not genuine ridges of the union.

### Removing facets also fails, for the entire 580-plane dictionary

A separate construction enlarges the six-tile polytope by dropping facets.
An exact compatibility graph tests every pair of the 580 planes, including
pairs that were not ridges before enlargement. Seventeen selected maximal
compatible subsets have 360--363 walls; a symmetric eighteenth subset has
168. Every one has an independently certified spacelike escape. Thirteen
also fail Vinberg's integrality test through an explicit nonintegral
length-two cyclic product. None has finite volume.

The stronger final obstruction covers **every subset of these 580 oriented
planes**, not just those eighteen selections. Walls 0 and 74, in the saved
zero-based ordering, are each individually necessary for finite volume:
after deleting either one from the full dictionary, an exact spacelike
vector strictly satisfies all 579 remaining inequalities. Any smaller
subset omitting that wall has the same escape.

Thus a finite-volume subset must retain both walls. Their Gram matrix is

```
[[4, -4],
 [-4, 6]],
```

giving forbidden squared cosine `2/3`. The boundary census verifies that
they already form an actual ridge of the six-tile polytope. That ridge
survives every deletion retaining both facets. Selected walls also remain
true facets because each already supports a full facet of the original
polytope. They therefore cannot both occur in a Coxeter result. This
closes the facet-removal construction completely for this finite dictionary.
It does not exclude adding new planes, gluing different pieces, or H20
Coxeter polytopes in general.
The [deletion-obstruction report](h20_gluing_independent_deletion_obstruction_20261013.md)
gives the two witnesses and an explicit timelike point on their surviving
ridge. Its standalone checker uses only the Python standard library and
replays the proof without an optimizer, network access, or project imports.
The optional larger archive independently checks a deletion escape for
every one of the 580 facets.

## What remains open

The exact H9 escape witnesses apply to the saved finite wall lists. They
do not prove that the full arithmetic reflection groups have infinite
covolume. The H20 obstruction closes the displayed finite plane dictionary,
not the dimension. Further work needs additional H20 planes or different
assemblies, or further compatible H9 cuts from genuinely different
basepoints and arithmetic slices.
