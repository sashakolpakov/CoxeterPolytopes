# Full-host obstruction, a new H20 construction lane, and remote H9 vertices

No new compact or finite-volume Coxeter dimension has been established.
This round does change the search substantially: the tested quartic
arithmetic host is now proved globally nonreflective, a new lift into H20
gives an exact 85-facet infinite-volume chamber, and a different geometric
construction cuts escapes that survived the earlier H9 cap banks.

## A complete stopping theorem for the quartic host

The [quartic report](fresh_h9_quartic_next_frontier_20261013.md) proves that
both lattices

\[
\langle-(T-1)\rangle+F_4+\langle1\rangle,
\qquad
\langle-(T-1)\rangle+F_4+F_4+\langle1\rangle,
\quad T=2\cos(\pi/8),
\]

have infinite-volume full reflection chambers. Two integral involutions
fix complete local chamber sectors at corners connected by the first two
certified nearest-mirror steps. Their product preserves that same chamber
and has H5 trace

\[
-14-8T+26T^2+14T^3\approx148.3085>6.
\]

It therefore has infinite order. Discreteness gives disjoint positive-volume
interior-ball translates, excluding finite chamber volume. Extension by
the identity on the other F4 block proves the H9 statement. The independent
audit re-enumerates both connecting intervals and verifies the exact
oriented-sector connection, integral inverses, form preservation, and traces.
This excludes every finite-volume reflection subgroup preserving these
specific lattices, not merely the saved wall prefixes.

Before this conclusion, a tighter census using shifted integral pairings
near the current corner found two additional first mirrors of norm
2T²−T³. Each complete interval needs 54 quadratic equations. Their orbits
produce a coherent 194-wall H9 cone; two much longer proposals, not certified
as nearest mirrors, also give a 202-wall subgroup cone. These constructions and all
independent replays are retained, but no further search in this host is
justified. The exact interval machinery remains useful for other hosts.

## A different route into dimension 20

The [H20 report](fresh_dimension_bridge_20261013.md) starts from the known
36-wall Vinberg H19 chamber. It retains 24 short roots and replaces each of
twelve long roots by two coupled lifts. Twenty-two additional walls complete
the eleven remaining deficient inherited cusp links. The resulting 70-facet
Coxeter chamber has all 107 inherited cusp links of full Euclidean rank 19.

An exact escape led to twelve compatible polar truncations, producing 82
facets. A separate Fraction implementation checks all 3,321 wall pairs and
the preservation of every inherited cusp link. The final exact strict
escape has norm 3849/32768 and center pairing −14839/128. Closing the inherited
cusps therefore does not close this chamber.

The search then replaces failed right-angle polar truncations with caps
having nineteen right angles and one prescribed finite angle. All 107
accepted entries are exactly reconstructed in real radical towers and
independently checked against the 82 retained walls. The same strict escape
survives all 107 simultaneously, so no subset of that defined bank can give
finite volume. See [the independent bank certificate](fresh_dimension_bridge_angle_repairs_20261013_independent.md).

Allowing eighteen right angles and two finite angles gives 933 numerical
proposals on three further target rays. Three selected caps are reconstructed
exactly, each with norm 2/27 and two prescribed π/6 angles. They retain all 82
old walls and are pairwise dotted, giving the coherent 85-facet H20 chamber
in [the exact extension](../data/enumerations/fresh_dimension_bridge_two_angle_20261013_exact.json).
They cut three specified hyperideal rays, but the same strict escape
survives. A [separate arithmetic audit](fresh_dimension_bridge_two_angle_20261013_independent.md)
checks all 249 new pair relations, target cuts, and the surviving escape.
Only these three two-angle proposals receive a positive exact construction
claim; the full 933-entry scout is numerical.

Any future finite-volume completion retaining the 70-wall stage would be
nonarithmetic: an existing normalized cyclic product is (11+4√6)/3,
which is not an algebraic integer. The H20 report also gives a complete
classification of 8,195 labelled one-coordinate lifts with the original
projected walls, and a separate obstruction to finite-volume completion of
its two golden branches. These results do not obstruct arbitrary completions
of the rational lifted branch used for the 85-wall construction.

## Cutting the old H9 remote escapes

The [remote-insertion report](fresh_h9_remote_caps_20261013.md) extends the
previous 79-facet geometric chamber to 85 facets. An exact orthogonal frame
in the complement of five old walls produces four mutually orthogonal caps.
They cut the remote witness surviving the earlier 3,857-cap bank and create
an actual 9A1 finite vertex.

A second construction uses two coupled caps and an exact π/3 anchor
condition. It cuts a different escape, closes an old open edge and creates
a finite path through H3+6A1 and 9A1 vertices. Every old wall is retained.
The three new vertices have an exact local frontier census: two finite
edges and 24 directed open ends, with the shared finite edge counted from
both of its new endpoints. A different strict escape proves that the final
85-facet chamber still has infinite volume.

The frame identities, finite-angle equation, actual vertices, facets and
strict escape signs have producer-free exact replays. A separate [source
review and higher-precision interval execution](fresh_h9_remote_caps_independent_20261013.md)
corroborate the 83- and
85-wall constructions; this second execution is explicitly not a second
interval-arithmetic implementation.

## New starting metrics

The [unequal-block and parameter-intersection scout](fresh_h9_asymmetric_blocks_20261013.md)
tests 388 fixed slices, 3,415 three-label curve intersections and 41,934
four-label intersections. Fifty-four of the latter full cap banks
numerically bound a region inside the Klein sphere, but mutually compatible
selections do not: all 47,195 maximal cliques in those 54 banks retain
outside vertices. This is bounded numerical reconnaissance, not an exact
parameter-space exclusion. It motivates using the remaining asymmetric
core metric parameter to enforce compatibility between cap families.

The continuation should use a different arithmetic host, new metric data,
or the explicit remaining H20/H9 gaps. The completed full-host theorem and
finite-bank exclusions prevent further work from merely extending known
dead branches. None of these infinite-volume constructions is asserted to
be a globally new abstract polytope.
