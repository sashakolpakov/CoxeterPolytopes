# Quadratic-five face constraints and the residual contact budget

This checkpoint develops the remaining high-dimensional arithmetic field
Q(√5) and sharpens the unrestricted finite-volume dimension argument.
The preceding [Q(√2) bound20 and norm31 obstruction checkpoint](sqrt2_dimension_bound_and_norm31_obstruction_2026_10_15.md)
remains frozen. No new finite-volume polytope or dimension is established.

## Dimension bounds and their status

The [dimension165 binder](fresh_nikulin_height_charge_20261015_dimension165.md)
combines an independently reviewed aggregate estimate with the frozen
geometric case split. Expensive exceptional affine factors consume at least
3/2 of the shared contact budget four, so at most two can fail unit charging.
This lowers the unrestricted candidate bound166 to165.

A further [dimension164 proof candidate](fresh_nikulin_defect_frontier_20261015_dimension164.md)
uses the fact that two ordinary E8 exceptions leave total retained cost
strictly below one elsewhere. This allows at most one retained coordinate;
an all-rank estimate bounds the remaining loss by2/3. The common-height window
is widened to21/20, and the aggregate loss drops from40/9 to38/9.
The [independent review](fresh_nikulin_defect_frontier_20261015_review.md)
passes: it reconstructs all40 exceptional models,71 single-contact and16
multiple-contact cases,65 subunit coordinates,580 convex certificates and
the unbounded classical-family arguments. The final binder independently
checks the numerical threshold and complete geometric branch split.

These unrestricted bounds impose no arithmeticity, compactness or simple-cusp
assumption. They inherit the previously reviewed narrow-facet geometry and
compact-facet bound143. They are research proof candidates, with external
validation, novelty and current-best literature status unestablished.
The new unrestricted estimates do not lower the separate arithmetic ceiling29.
The earlier internally reviewed compact arithmetic bounds20 over Q(√2) and
Q(cos(2π/7)) still leave Q(√5) as the possible ground field in dimensions22–29.

## Complete finite reductions over Q(√5)

The [Lannér catalogue](fresh_sqrt5_lanner_catalogue_20261015.md) and its
[independent review](fresh_sqrt5_lanner_catalogue_20261015_review.md) establish
the full finite-label alphabet2,3,4,5,6,10 and the allowed minimal indefinite
pieces:12 triangles,6 order-four types and4 order-five types, besides the
separately parameterized dotted pair. No integer label cutoff is used.

The [union enumeration](fresh_sqrt5_quad_exhaustion_20261015.md) and
[independent enumeration](fresh_sqrt5_quad_exhaustion_20261015_incidence_review.md)
cover every pair of opposite minimal indefinite pieces in a quadrilateral
with no marked corner path of length at most two. They include shared white
nodes and all cross contacts:153 marked piece representatives and127,875
identifications produce5,440 labelled unions and1,180 marked classes.
This is a complete union census under its explicit no-short-path hypothesis,
not a classification of all external white attachments.

The [principal arithmetic filters](fresh_sqrt5_quad_filters_20261015.md),
independently checked in the [separate review](fresh_sqrt5_quad_filters_20261015_review.md),
remove784 classes:617 force two physical negative directions,89 violate
the cyclic ground field, and78 fail conjugate positivity. There remain362
finite-label classes and34 parameterized dotted classes. Exact positive
radical rescalings place accepted finite Grams in Q(√5); both embeddings
are tested without rounded eigenvalues.

The [triangle rule](fresh_sqrt5_triangle_weights_20261015.md) and
[independent review](fresh_sqrt5_triangle_weights_20261015_review.md) cover
all76 three-black markings of the catalogue, including arbitrary external
white extensions. Distance-at-most-two weights plus the specified H4 bonuses
give every triangle weight at least2. The finite-component totals are at most
2n−2 for n≥5. Quadrilateral coverage is a separate obligation.

The [complete one-white contact domains](fresh_sqrt5_quad_contacts_20261015.md)
and their [independent reconstruction](fresh_sqrt5_quad_contacts_20261015_review.md)
reduce the396 surviving classes further:

- 84 already receive sufficient weight from an entire H4 component.
- 2 dotted classes contain explicitly certified affine principal diagrams and
  are incompatible with compactness.
- 112 singular finite cores admit no further touching white node while the
  quadrilateral remains unpaid. Their relevant corner components are entire.
- 44 nonsingular finite cores likewise cannot receive unpaid touching nodes;
  the full-star rank condition therefore excludes them from the unpaid case.
- 130 nonsingular finite cores retain458 necessary one-node profiles, and24
  dotted families retain369 unresolved profiles.

The contact search exhausts4,811 structurally possible one-node profiles.
The assertions about arbitrary external whites use an explicit first-contact
argument. “Sealed” is conditional on remaining unpaid:19 singular cores can
instead acquire an H4 payment. No multi-node closure or arithmetic distance
classification is asserted for the remaining130/24 classes.

## A precise weighting obstruction, and a possible way past it

The [additive-weight diagnostic](fresh_sqrt5_sealed_weights_20261015.md) has
an exact small dual obstruction. Fix weight1 on every distance-one or
distance-two marked pair, allow only nonnegative additional weights depending
on the entire connected component and marked pair, and give disconnected
marked pairs weight0. Nine local quadrilaterals force nine additional E8
pair weights, and six force six additional D7 pair weights. Hence their
component totals are at least23 and21 respectively. The rank29 forest
D7³+E8 would then cost at least86, above the84 needed by the intended face-count
argument. This excludes that specified weighting scheme, not general
component weights, cross-component payments or dimension29 polytopes.

The fifteen witness stars are actual compact local quadrilateral faces:
exact arithmetic checks90 pair intersections,60 feasible future timelike
corner rays,30 rejected opposite intersections and strict future interior
points. They are not asserted to coexist in one ambient finite-volume
polytope. The [independent review](fresh_sqrt5_sealed_weights_20261015_review.md)
passes and explicitly verifies the scope of the obstruction.

Allowing existing weights to decrease gives a feasible
[free-weight relaxation](fresh_sqrt5_free_weights_20261015.md):210 included
face constraints on28 represented finite component types have exact optimal
common deficit21/10 from twice their rank. Exact rational primal and dual
certificates verify this finite optimization. It is not a global weight rule.
Missing obligations are explicit: quadrilaterals formerly paid by short paths,
arbitrary triangle extensions, the nonsingular/dotted frontier, and all
spherical components and induced embeddings up to rank29.

## Continuation and archive scope

The next arithmetic enumeration removes the no-short-path hypothesis, since
reducing those weights invalidates its automatic quadrilateral payment.
Its bounded pilot completed six left-piece blocks and passed an exact
mid-DFS pause-and-resume check; the full enumeration continues separately.
Further nonexistence work retains the two individual endpoint budgets, each
at most2, instead of only their sum. Those ongoing investigations use fresh
prefixes and are outside this frozen checkpoint.

The geometric H9/H20 constructions, compact-H9 incidence cursor and closed
arithmetic norm31 host retain their preceding status. No additional search
of a lattice already proved nonreflective is justified by these face results.

The [archive auditor](../scripts/audit_quadratic_five_residual_checkpoint_20261015.py)
and [manifest](quadratic_five_faces_and_residual_charge_2026_10_15_audit.json)
bind new reports, recursive source hashes and current entrypoints while
checking preservation of the preceding frozen archive. The completed producer
checkpoint `fresh_sqrt5_quad_exhaustion_20261015_running.json` is explicitly
excluded; the authoritative completed enumeration is bound. Archive integrity
is distinct from mathematical verification and from the stated scope of each
finite diagnostic.
