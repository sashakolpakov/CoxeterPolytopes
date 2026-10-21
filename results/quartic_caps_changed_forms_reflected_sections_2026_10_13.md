# Quartic caps, changed arithmetic forms, and reflected sections

This round has **not produced a new compact or finite-volume Coxeter
polytope in a new dimension**. It develops four different constructions:
coherent H9 caps over larger number fields, a changed arithmetic H9 space,
a new H9 seed with nonsymmetric caps, and H20 sections after changing the
ambient H21 chamber. Exact checks distinguish successful local repairs
from the remaining volume or angle failures.

## Coherent H9 completion over larger fields

The [full construction report](fresh_h9_quartic_completion_20261013.md)
records a limitation exposed by the final dependency audit: the retained
core contains four parallel pairs, between old cap 31 and seed 8, old cap
51 and seeds 8 and 9, and old cap 197 and seed 1. Their positive null
vectors remain feasible under every acute-compatible added wall.
**Keeping this core therefore excludes compact completion, regardless of
how many further walls are added.** A stronger cusp obstruction also
excludes finite-volume completion of this fixed core by acute-compatible
wall addition, as explained below. The ten-wall seed alone remains open.

The cap search retains a repaired 32-wall core, including a quartic cap
that replaces one incompatible earlier wall. Solving the full relation to
that quartic anchor allows new normals beyond coefficientwise solutions
over the golden field. Target-guided distance choices produce 80 caps;
an adaptive pass preserving a 35-wall configuration produces another 80;
a pass preserving a 43-wall configuration produces five more.

An exact mirror comparison establishes **165 distinct mirrors outside the
preceding 2,865 saved entries**. These are new relative to the archived
search, not a worldwide novelty claim. Some close gaps that the preceding
compatible selections left open. Compatible selections reach 44
walls while preserving the repaired core. The tested configurations still
have exact strict spacelike escapes, so the number of walls is not evidence
of finite volume.

A subsequent construction releases the prescribed distance to the quartic
anchor. It instead imposes angle π/3 with an older polar cap and solves for
two seed distances. With `phi²=phi+1`, the unique positive root of

```
(−20+25phi/2)s² + (−25+15phi)s + (−42/5+24phi/5) = 0
```

gives nonzero seed-pairing magnitudes `sqrt(2+phi)s` and
`sqrt(2+phi)(s+1)`. Exact checks give norm two, all 44 required Coxeter
relations, a common strict interior, and a cut of a previously persistent
escape. This is a **166th distinct new mirror**, outside all 3,030 entries
preceding it. The resulting **45-wall H9 configuration still has an exact
strict spacelike escape**. Both the useful local repair and the remaining
infinite volume are certified.

That escape simplifies to a feasible spacelike ray with nonzero dual
coordinates `y6=2−phi`, `y9=1`. A second targeted cap has seed-pairing
magnitudes `sqrt(2)` and `sqrt(12+20phi)` in those positions, with all
other seed pairings zero. It meets seed wall 6 at π/4 and is disjoint
from all 35 retained nonseed walls. This **167th new mirror** closes the
opposite gap and extends the compatible configuration to **46 walls**.
A different exact strict spacelike escape still remains; its dominant
coordinates are 3, 8, and 9. This is another local repair, not a
finite-volume completion.

The literal cap archives are
[targeted](../data/enumerations/fresh_h9_quartic_completion_20261013_targeted.json),
[adaptive](../data/enumerations/fresh_h9_quartic_completion_20261013_adaptive.json),
and [second adaptive](../data/enumerations/fresh_h9_quartic_completion_20261013_second.json).
The [mirror comparison](../data/enumerations/fresh_h9_quartic_completion_20261013_mirrors.json)
and [44-wall geometry](../data/enumerations/fresh_h9_quartic_completion_20261013_second_geometry.json)
record the exact distinctions and remaining escape. The
[two-distance polar construction](../data/enumerations/fresh_h9_quartic_completion_20261013_polar.json)
contains the new quadratic equation, cap, all retained-wall relations,
and the 45-wall escape. A separate
[independent reconstruction and review](fresh_h9_quartic_completion_20261013_polar_independent.md)
checks the unique positive root, all 44 pairings, distinctness from the
preceding entries, and the strict future spacelike escape; forged cap
and escape certificates are rejected.
The [opposite cap and 46-wall escape](../data/enumerations/fresh_h9_quartic_completion_20261013_opposite.json)
record the second targeted extension. Its
[independent reconstruction](fresh_h9_quartic_completion_20261013_opposite_independent.md)
checks all 45 relations, the cut of the old ray, distinctness from all
3,031 preceding entries, and the new strict escape. Corrupted cap and
escape certificates are rejected here as well.

## A genuinely changed arithmetic H9 space

The [arithmetic construction report](fresh_h9_changed_form_20261013.md)
uses the cubic field `T^3−3T−1=0`, at its largest real root, and the form

```
<-2T> + [[2,−T],[−T,2]] + A6 + <2>.
```

Its determinant has odd valuation at the inert prime 7. This proves it
is not field-similar to six specified cubic-field H9 spaces already in
the repository. Complete bounded root censuses, with independently
reconstructed time and spatial bounds, produce 23 H9 walls through height
100. They satisfy all Coxeter pair gates and have a strict timelike
interior, but an exact spacelike escape remains.

A local congruence argument excludes every primitive reflection root in
the binary plane generated by the time direction and the relevant A6 end
weight, at every height. The follow-up then leaves that plane: it permits
two or three A6 covectors, both binary-block coordinates, and the scalar.
Three searches with 81 time/norm controls each, through height 1000, find
no cutter of the displayed escape. The independent replays pass.

The bound-independent obstruction concerns only that binary plane. The
larger searches exhaust only their declared slices and height bounds;
the full arithmetic host remains undecided.

## A new symmetric seed and seventy nonsymmetric caps

The [symmetric-filling report](fresh_h9_symmetric_filling_20261013.md)
starts from an 18-wall noncompact H9 chamber with eight pairs of leaves.
It allows simultaneous changes to 112 old right angles. Rank ten forces
two algebraic branches, and explicit escape formulas exclude compactness
for the entire family when extra walls pair constantly with each leaf
species. Every changed branch has infinite volume. The unchanged parallel
case is excluded only from compactness.

To leave that restriction constructively, a new golden-field seed admits
**70 exact caps with nonconstant leaf pairings**, at angles π/2, π/3,
and π/5. Its Gram signature is `(9,1,8 zero)`, with a positive semidefinite
rank-ten conjugate. The latter does not establish arithmeticity: an
explicit Cartan cyclic product is nonintegral.

Each cap cuts its associated balanced escape. Between caps, exact pairings
permit only 35 complementary pairs and no triples. Every such pair gives
a 20-wall Coxeter cone with a separately checked strict spacelike escape.
Thus this particular 70-cap dictionary cannot complete the seed; arbitrary
other nonsymmetric caps are not excluded. Independent Fraction replays
check the polynomial obstruction, all cap formulas, all 2,415 cap pairs,
and all 35 escapes. A separate agent review found no issue.

The subsequent [complete finite leaf-angle census](fresh_h9_symmetric_finite_caps_20261013.md)
shows these 70 caps exhaust every norm-two real normal with leaf pairings
in `{0,−1,−phi}`, allowing arbitrary real center pairings. Thirty
assignment orbits give 59 real solutions: 52 have a positive center
pairing, six have a forbidden pairing strictly between −1 and zero,
and one orbit gives exactly the existing caps. An independent arithmetic
implementation replays the complete census. A separate numerical scout
with more finite leaf labels examines 333 orbits and also finds only that
orbit; it does not enlarge the exact exclusion.

## H20 sections after repairing the ambient chamber

The [reflected-section report](fresh_h20_reflected_sections_20261013.md)
doubles the original Borcherds H21 chamber across wall 169. The union is
an exactly certified convex finite-volume chamber with 388 distinct
boundary facets. For one earlier section plane, this really replaces an
actual forbidden ridge by an actual right angle. Another forbidden ridge
survives, so that section is still not Coxeter.

Six new irrational section planes retain a complete Coxeter cusp and
produce actual right angles involving new reflected ambient walls. All
six are finite-volume H20 polytopes. Each nevertheless has another exact
actual ridge whose angle is not a Coxeter angle.

A separate finite-star obstruction covers all twelve fixed planes from
the preceding round. Any finite spherical reflection star generated by
a subset of the 210 original facet reflections and containing the
original chamber must cover certain actual bad ridges by its reflecting
walls. Exact Gram obstructions make those covering requirements
inconsistent. The twelve certificates use 56 actual-ridge clauses,
33 nonelliptic-pair clauses, and 11 larger null-vector clauses. Independent
geometry and Boolean replays pass and reject corrupted certificates.

This excludes only these finite spherical stars and their relevant
subunions containing the original chamber, for those fixed planes. It leaves open
changed planes and more general chamber assemblies. The six constructive
sections above explicitly change the plane.

## Status

The new roots and angle repairs are usable construction data. None is a
new finite-volume Coxeter polytope, and the scoped obstructions do not
settle a dimension-wide existence question. The research goal remains
open. Detailed reports separate finite bounds, exact continuous
obstructions, geometry certificates, and independent replays.

For a compact H9 continuation, the inherited parallel pairs must first be
broken by replacing walls. Removing old caps 31, 51, and 197 removes these
four pairs, but the resulting 43-wall subset still contains the exact
affine triangles `{seed3,seed4,old496}` and
`{seed4,seed5,old496}`. Each has diagonal two, off-diagonal minus one,
and positive null vector `(1,1,1)`. They also prohibit compact completion
by acute-compatible wall additions. Compactness therefore requires
further replacements, such as removing old cap 496 or seed wall 4, and
checking any other remaining affine circuits. The ten seed walls alone
have no affine subset: any positive semidefinite subset omits seed 8 or
seed 9, and the corresponding maximal diagrams are positive-definite
`A7 + I2(5)` and `E7 + I2(5)`.

The [exact compactness audit](../data/enumerations/fresh_h9_quartic_completion_20261013_compact_scope.json)
records this 43-wall deletion checkpoint, its inherited strict escape,
and five verified affine circuits: the two triangles and three affine
`A3` diagrams. With the seed fixed, breaking these five circuits requires
changing at least four further nonseed walls: adaptive cap 19, old caps
1028 and 496, and either old cap 310 or 837. This minimum concerns only
the displayed circuits; it is not an exhaustive affine census or a
certificate that those changes suffice for compact completion.

The [four-wall cusp certificate and standalone replay](fresh_h9_fixed_core_cusp_obstruction_20261013.md)
give a stronger reason to replace the original retained core. The parallel
pair `{old31,seed8}` is orthogonal to both seed 6 and
seed 7, which meet at π/5. Its positive null ray survives every
acute-compatible added wall. If the completed chamber had finite volume,
this ray would have a bounded Euclidean Coxeter cusp link containing that
π/5 pair. This is impossible: the product of the two reflections has trace
`5+phi` on the eight-dimensional horosphere, whereas preservation of a
full cusp translation lattice would make its trace an integer, by
[Bieberbach's theorem (Gelander, Theorem 2.9)](https://www.math.utah.edu/pcmi12/lecture_notes/gelander.pdf).
The same
obstruction applies to each of the four inherited parallel pairs.

Thus the **fixed 44-, 45-, and 46-wall cores have no finite-volume completion by adding
acute-compatible walls**, and further construction must change retained
walls. This statement does not exclude other completions of the ten-wall
seed. The parallel-free 43-wall deletion checkpoint avoids these four
particular cusp obstructions, but its five displayed affine circuits
still exclude compact completion without further changes.

The [37-entry target archive](../data/enumerations/fresh_h9_quartic_completion_20261013_next_targets.json)
preserves the last spacelike escape at index 36, with dominant dual
coordinates 3, 8, and 9. It is not proportional to any of the earlier 36
targets. It remains useful geometric data after changing the core;
cutting that ray while retaining all 46 walls cannot solve the
finite-volume problem.
