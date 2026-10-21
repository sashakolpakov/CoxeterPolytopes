# A9 nonreflectivity, a compact prism face, and complete cap orbits

**No new compact or finite-volume Coxeter example has been established.** The
strongest result of this continuation is a full-lattice obstruction: the
quartic A9 arithmetic host is nonreflective. Separate geometric constructions
produce a compact H4 face in H9 and exact compatible extensions in H20 and
H22, but every resulting whole chamber has infinite volume.

| Search | Exact result | Consequence |
|---|---|---|
| Quartic arithmetic H9, A9 host | Complete edge paths and an infinite-order integral chamber symmetry | The full reflection chamber has infinite volume; this lattice is closed |
| Geometric H9 | 96 actual facets, including an actual compact H4 prism face | A strict spacelike escape survives in the whole chamber |
| Geometric H20 | 176 actual facets, a repaired cusp, and two successive gap cuts | A new strict escape survives; further walls all dotted to the last facet cannot complete volume |
| Geometric H22 | 2,736,428 actual facets, for every new label k≥6 and old center label m≥3 | One strict escape proves infinite volume throughout the family |

## A full arithmetic obstruction from complete edge walks

The host is `L=⟨−(T−1)⟩⊕A9`, over the quartic field with
`T=2cos(π/8)` and `T⁴−4T²+2=0`. Eight further complete first-mirror
calculations inspect 1,592 exact relative-orbit interval points. With reversal,
the resulting prefix contains 35 coherent walls, 26 known finite vertices and
25 certified connecting edges.

The decisive step leaves the previously explored H8 facet. Its endpoint and
another reached corner have complete spherical stabilizers of type
`A4+A4+A1`. An explicitly integral Lorentz isometry maps the entire outward
corner sector of one to the other. Because complete first-mirror paths place
both corners in the same full reflection chamber, this map preserves that
full chamber, not merely the saved wall prefix.

Its trace is exactly greater than ten at the distinguished real embedding,
proving infinite order. Discreteness of the arithmetic isometry group then
forces infinite chamber volume. This closes the third quartic host examined
by this method. The unenumerated finite-prefix directions do not weaken the
global conclusion, and further wall enumeration in this lattice is unnecessary.

The [full report](fresh_h9_a9_outer_20261013.md) includes the explicit matrix,
integral inverse, oriented paths, complete corner stabilizers, and trace.
The separate Fraction matrix/path replay and the complete supporting edge
replays pass. Certified PARI ideal-class and unit computations are an explicit
dependency; no independent implementation of those algorithms is claimed.

## A compact H4 face inside an infinite-volume H9 chamber

Two coupled exact caps extend the earlier 94 facets to 96. One joins four
new endpoints to the old compact tetrahedron, creating an actual compact
Coxeter H4 prism with eight vertices. The other cuts the preceding chamber
escape. All 189 new pair relations pass, and all 1,024 possible new
positive-semidefinite subdiagrams are classified, with no affine circuit.

The reached finite component has ten vertices, nineteen finite edges and
52 open rays. A new exact strict escape remains. Thus compactness of this
face is a concrete geometric improvement, but it gives neither compactness
of the ambient H9 chamber nor a new lower-dimensional Coxeter type.

The [construction and frontier](fresh_h9_face_bridge_20261013.md) and
[independent symbolic/source review](fresh_h9_face_bridge_20261013_review.md)
include the full prism boundedness argument. The final precision replay
checks 887 strict enclosures at 160 digits, nested inside the original
enclosures; it reuses the interval arithmetic implementation.

A follow-up [all-label prism-lift obstruction](fresh_h9_prism_lift_20261013.md)
rules out using this exact H4 prism as the base of a compact H5 prism while
keeping all six side hyperplanes. The proposed top must have finite Coxeter
angles with all six sides. Their linear dependence contains the nonabelian
quartic number `√((9+3φ)/11)`, whereas finite-angle cosines belong to abelian
extensions of Q. This forces one pairing to vanish. Positive definiteness at
an H4 corner forces four further pairings to vanish, leaving only a normal
to the already existing base hyperplane. There is no proper new top.

This rules out the specified one-top prism extension for every finite angle
label. It does not exclude multiple caps, changed sides, or arbitrary H5
faces inside another H9 construction. The standalone Fraction certificate and
[independent field/geometry review](fresh_h9_prism_lift_20261013_review.md)
both pass; the same endpoint argument also excludes ideal top vertices in
this ansatz.

## H20: repaired cusp, coherent orbit, and an obstruction to further copies

Six new walls first complete the deficient cusp to affine components of
ranks `8+8+1+1+1=19`, with strictly positive kernels. This proves the
horospherical link bounded, rather than inferring boundedness from rank alone.
A coupled cap then frees a previously constrained anchor to be dotted and
cuts the old escape. Seventy additional coherent copies under a 72-element
isometry group give 175 walls; a second coupled cap gives 176.

The final cap cuts the preceding exact escape but exposes another. Its only
intersecting neighbors are 18 retained walls. Projecting a proposed
completion into this facet leaves a spacelike tangent direction perpendicular
to all 18 constraints. Every wall dotted to the facet projects to a timelike
normal and stays strict at the resulting null escape. A small perturbation
proves that no finite addition consisting only of walls dotted to this facet
can produce finite volume. Intersecting or parallel additions and wall
replacements remain outside this obstruction.

The [construction and general obstruction](fresh_h20_coupled_completion_20261013.md)
and [independent source review and replays](fresh_h20_coupled_completion_20261013_review.md)
pass through the exact 176-wall stage. A numerical 247-wall scout is explicitly
not an exact chamber certificate. No arbitrary H20 completion is excluded.

## H22: a complete orbit and an unbounded-label family

A rational old-space vector is determined by the 19 facets of an old H2 face
and three additional pairing equations. Its complete orbit has 241,920
elements. Taking both height signs adds 483,840 walls to the previous
2,252,588 facets.

The new norm and height are `q_k=1/(2cos²(π/k))` and
`H_k=√(31/3+q_k)`. Exact orbit extrema prove all pair conditions for every
integer `k≥6`; the inherited center parameter remains any `m≥3`.
The same full-orbit construction fails the Coxeter angle test for `k=3,4,5`.
There is no finite angle cutoff hidden in this family calculation.

At `k=6`, a new parallel pair has a deficient cusp link with two missing
horospherical directions. An exact strict spacelike perturbation near that
cusp remains feasible for every `k≥6` and `m≥3`. Its square is
`647099999951/35000000000000000>0`. Hence every member of this entire
family has infinite volume, despite cutting the previous escape.

The [family certificate](fresh_h22_gap_face_family_20261013.md) and
[independent source/proof review and replay](fresh_h22_gap_face_family_20261013_review.md)
check the orbit, actual facets, small-label failures, and uniform escape.
They retain explicit hash-anchored dependencies on the preceding exact
chamber certificates. Changing the integer label alone cannot repair volume;
new walls or a different orbit could still change the conclusion.

## What this changes in the search

The arithmetic route now needs a genuinely different lattice or quadratic
form. The geometric route needs compatible connections between the remaining
ends or coherent replacement of obstructing choices. More copies of the
same caps, more closed edges at one corner, and more completed inherited
cusps are insufficient when an exact global escape survives.

The [Vinberg–Allcock audit](vinberg_allcock_search_strategy_20261013.md)
remains the algorithmic guide. Vinberg orders mirrors by distance, while
Allcock traverses chamber edges; the published edge-walking algorithm is
also a lattice algorithm. Diagram extensions can guide a broader geometric
search, but require separate metric, firstness where applicable, and global
closure proofs.
