# Fifteen-facet H9 search and geometric H20 constructions

Date: 2026-10-14. **No new finite-volume polytope or dimension is established.**
The preceding [direct-Gale and two-coordinate checkpoint](direct_gale_and_two_axis_constructions_2026_10_14.md)
remains frozen. This checkpoint changes both the combinatorial search
family and the priorities for geometric completion.

## Literature restrictions applied before further search

The user's arithmetic warning is correct for noncompact examples.
Esselmann's theorem gives `n<=19` or `n=21`; there is no noncompact
arithmetic example in H20 or H22 and above. The publisher's
[original abstract](https://www.sciencedirect.com/science/article/pii/S0022314X96901419)
states the maximal reflective integral Lorentzian lattice rank as 22,
meaning hyperbolic dimension 21. Borcherds' updated
[original construction](https://math.berkeley.edu/~reb/papers/al/al.pdf)
records the uniqueness in lattice ranks at least 21, including the
missing rank 21 / hyperbolic dimension 20. The range is explicitly
restated in [Belolipetsky, §1](https://arxiv.org/pdf/1506.03111).

Compact arithmetic examples are a different category: Vinberg's
[original theorem](https://www.mathnet.ru/eng/mmo444)
excludes dimensions at least 30, not all dimensions above 21.
The same compact upper bound applies without arithmeticity. Thus H9
is a compact-dimension target; a nonarithmetic H20 example would fill
the missing finite-volume dimension, while H22 would raise the maximum
known dimension. The [dimension review](arithmetic_dimension_bounds_and_effective_high_dimension_search.md)
records these distinctions and the recent literature check.

Allcock's [published edgewalking algorithm](https://web.ma.utexas.edu/users/allcock/research/edgewalking3.pdf)
explicitly uses Z-lattices with a rational Lorentzian form. Both it and
Vinberg's lattice implementation address arithmetic hosts. An unrestricted
geometric search can use edge traversal, but must independently construct
compatible walls and establish global finite volume. It cannot invoke the
published lattice algorithm as a complete nonarithmetic search.

## A stronger stopping rule for the golden-angle branch

The [face obstruction](fresh_noncrystallographic_face_gate_20261014.md),
with a [separate review](fresh_noncrystallographic_face_gate_20261014_review.md),
combines Borcherds' Coxeter-face theorem and the Euclidean cusp trace
condition: an actual ridge with finite label outside `{2,3,4,6}` in a
finite-volume Hn Coxeter polytope is itself a compact Coxeter H(n−2)
polytope. In particular, retaining a pi/5 ridge in H20 requires compact
H18. This is a consequence of established theorems, not a new unrestricted
dimension bound or a proof that H20 is impossible.

We constructed actual finite corners on all sixteen central pi/5 ridges
of the weighted 39-wall core: 288 small caps give 327 actual facets and
exactly sixteen finite vertices. One additional wall closes one finite
edge, giving 328 facets and seventeen finite vertices. The
[construction report](fresh_h20_central_ridge_20261014.md) and
[independent different-field-basis audit](fresh_h20_central_ridge_20261014_review.md)
verify both stages. Both retain a strict norm-one escape and have infinite
volume. The complete edge census also includes 524,288 seed-only edges
missed by traversal from the finite vertices.

More decisively, any completion in the specified orthogonal-contact
scheme would induce a compact right-angled H18 face. The compact
right-angled dimension bound is four. This excludes the entire stated
completion scheme, regardless of additional cap count. Stop this scheme;
nonright contacts or replacement of the ridge are outside the exclusion.

## Nonarithmetic hybrid lift

The [coupled-lift report](fresh_h20_coupled_lift_20261014.md) starts from
the previous 60-wall H20 lift of the H18 hybrid. Keeping all old projected
directions and replacing the short-wall heights by arbitrary finite
planar sets does not repair the deficient cusp under the stated
positive-convex-recovery assumption. New projected normals are necessary
within that construction.

The new construction adds 32 jointly compatible walls to the 60-wall
chamber. Exact checks certify 2,416 new pairs: 592 right and 1,824
strictly dotted. It creates sixteen positive-definite rank-20 corners,
cuts all four sign images of the saved escape, and preserves the three
previously repaired bounded cusp links. A new strict rational escape
still proves infinite volume, and all 320 outgoing germs at these
sixteen corners remain open. These are local geometric constructions,
not a finite-volume certificate or a complete vertex census.

The [separate geometry/source review](fresh_h20_coupled_lift_20261014_review.md)
passes. Its no-import Fraction implementation checks all 92 norms and
strict center/escape signs, four support frames and sign isometries,
1,472 finite-vertex inequalities and 96 strict new-cap/cusp signs.
The full wall-pair interval checks are source-reviewed and replayed in
the producer, not independently reimplemented by this review.

The retained hybrid relations force any finite-volume completion to be
non-quasi-arithmetic. Its finite ridge labels 2, 3 and 4 avoid the compact
H18 requirement above. Further work must change the missing endpoint
contacts or the core; simply adding more isolated small corners is not
a completion strategy.

## Fifteen-facet compact H9 search

The [fifteen-facet report](fresh_h9_fifteen_facet_20261014.md) considers
one codimension-two-through-five cut of each of the 770 saved realizable
fourteen-facet Gale nerves. It is a new bounded family, not an enumeration
of all fifteen-facet H9 polytopes.

| Stage | Count |
|---|---:|
| Literal source faces considered | 2,124,486 |
| Cuts potentially repairing inherited obstructions | 75,311 |
| Verified operation orbits | 45,888 |
| Distinct fifteen-facet nerves | 45,679 |
| Excluded by inherited local or block constraints | 39,931 |
| Excluded by two new exact local patterns | 4,240 |
| Excluded by complete necessary-label searches | 174 |
| Still unsearched | 1,334 |

The label searches visit 2,847 states. They allow missing triangles and
all finite labels: labels 2 through 6 are literal, and category 7
represents every integer label at least 7, with the necessary isolated
I2(m) restriction at spherical vertices. No failed numerical metric fit
is used as an exclusion. No complete necessary-label survivor was found
in the 174 selected cases.

The [separate incidence replay](fresh_h9_fifteen_facet_20261014_incidence_independent.json)
reconstructs all 9,817,894 vertex incidences and checks the 209 duplicate
operation identifications. Separate arithmetic replays check the
[first 150 label searches](fresh_h9_fifteen_facet_20261014_search_independent.json)
and [final 24 searches](fresh_h9_fifteen_facet_20261014_metric_priority_independent.json).
They reuse DFS propagation and the classical Lannér tables, while using
a separate finite-diagram classifier and finite-field determinants.
The [peer source/proof review](fresh_h9_fifteen_facet_20261014_review.md)
has the narrower scope stated in its report.
The remaining 1,334 cases are genuine untested work, not candidates
known to admit a Coxeter metric.

## Archive and next search conditions

The [archive audit](finite_ridge_corners_and_fifteen_facet_search_2026_10_14_audit.json)
checks report links and recursively declared source hashes. Mathematical
replays and independent reviews are separate from this provenance check.

Continue only outside closed families. For compact H9, test the saved
untested incidence types before expanding the family. For noncompact
H20, retain the nonarithmetic gate and require an actual endpoint/cusp
closure mechanism. The golden orthogonal-face family is closed by a
theorem and should not be enlarged. None of these exclusions settles
the unrestricted existence questions.
