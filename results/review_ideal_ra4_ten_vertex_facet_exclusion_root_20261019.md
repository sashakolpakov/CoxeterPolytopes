# Figure7 exclusion: completed written review and exact model replays

19 October2026. **PASS for the restricted structural theorem**, with the
published and inherited inputs below. This note supplements the frozen-by-
hash preparation report; it does not change the current16-page manuscript
or its older conditional Lean interface. No new dimension ceiling follows.

The [full proof](ideal_ra4_ten_vertex_facet_exclusion_research_20261019.md)
has SHA256
`35a8bddbffd215251f815ee02de4964feec0eeee03e1a177d042afd2dcd60763`.
The [independent full mathematical review](review_ideal_ra4_ten_vertex_facet_exclusion_resource_20261019.md)
passes. Root checked the metric derivation, all three cusp classes and
the unbounded multiplicative-ratio argument directly while writing it.

The complete facet spectrum O,D,B,A_k permits no B or A_k. The new B
model has cusp ratios s,s³ with s²=phi, and its AEJ triangle is long
at A and short at E,J. The possible three-distinct boxes involving B
would have exactly one opposite pair of A5 faces and only B neighbors;
the required pentagonal ridge cannot occur. The other three-distinct
box and equal-long boxes are excluded by the retained side capacities.
The remaining single-long boxes force the same AEJ neighboring facet
to be B at one cusp and O/D at the other two, impossible. Only then
is the prior all-antiprism theorem applied.

The independent exact checker was read in full. Its quartic arithmetic
uses the irreducible relation s^4−s²−1=0. Its physical sign test first
compares in Q(phi) using sqrt5, then compares a+s*b by the sign of
a²−phi*b² when the summands have opposite signs. No floating point,
numerical fitting or general symbolic-algebra library is used.
The wall normals and null vertices are built independently from the
explicit centers/radii; inversion computes every cusp, not only its
representative orbit.

The checker source SHA256 is
`9c3c9b478f515255754dea9bb407054b481a122babe611d40c43fd0c008f3451`.
Both resource-gated single jobs passed:

* [Generation](ideal_ra4_ten_vertex_model_generation_20261019.json):
  0.094seconds,56,754,176-byte process-tree peak.
* [Separate replay](ideal_ra4_ten_vertex_model_replay_20261019.json):
  0.103seconds,56,803,328-byte peak.
* [Exact model certificate](ideal_ra4_ten_vertex_model_20261019.json):
  12 unit wall normals,10 feasible ideal null vertices,20 edges,
  120 vertex/wall products,66 wall pairs, the complete12 face cycles,
  both automorphisms and all10 inversion rectangles, including the
  AEJ long/short/short signature.

These finite checks support the explicit model calculations; they do
not enumerate arbitrary polytopes or prove the universal box argument.
The written proof derives the metric inside an existing Figure7 facet
using its combinatorial symmetries and standard three-dimensional
right-angled uniqueness. The numerical consequence below additionally
uses the published complete small-facet classification through ten
vertices. The earlier O/D cutting/chain theorem, all-antiprism exclusion
and ideal5 kappa identity are separately inherited reviewed proofs.

Consequently every entirely ideal RA4 polytope whose facets all have
at most ten vertices is an O/D chain. Every positive-kappa ideal4 facet
of an entirely ideal RA5 polytope contains a three-face with at least
eleven vertices. The universal ideal dimension ceiling remains6.

    python3 -B scripts/run_bounded_research_check_20261018.py -- scripts/check_ideal_ra4_ten_vertex_model_20261019.py --verify
