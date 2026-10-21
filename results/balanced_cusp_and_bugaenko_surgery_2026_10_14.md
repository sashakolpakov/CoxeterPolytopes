# Balanced cusp construction and compact-H9 incidence surgery

**No new compact or finite-volume hyperbolic Coxeter polytope is established.**
The new H20 construction has 418 actual facets and 392 finite vertices, with
an independently verified strict infinite-volume escape. The compact-H9
branch supplies 2,953 incidence types and additional exact obstructions.

The dimension restrictions were checked again against
[Esselmann's original paper](https://www.sciencedirect.com/science/article/pii/S0022314X96901419),
[Borcherds' updated thesis remark](https://math.berkeley.edu/~reb/papers/thesis/thesis.pdf),
and [Belolipetsky's precise statement](https://arxiv.org/pdf/1506.03111).
Noncompact arithmetic examples occur only in dimensions at most 19 or in
dimension 21. The general arithmetic and compact upper bound is 29.
The targets remain nonarithmetic finite-volume H20 and compact H9.

Vinberg's algorithm enumerates mirrors by distance from a chosen point.
[Allcock's published edge-walking algorithm](https://web.ma.utexas.edu/users/allcock/research/edgewalking3.pdf)
finds neighboring vertices through rank-two Lorentzian lattice calculations;
it explicitly works over Z. Its geometric idea motivates endpoint and
incidence searches, but using it outside that lattice setting requires fresh
metric, completeness, and finite-volume checks. It is not an existing
arithmetic-free enumeration theorem.

The balanced cusp has affine diagram `A6 × A6 × A7`, 22 seed walls and
392 edge germs. At unit inradii every cross-section vertex has squared norm
196. The [complete endpoint census](fresh_h20_balanced_endpoint_20261014.md)
has 37 positive-definite and four semidefinite profiles. It permits every
real endpoint normal and every omitted-wall Coxeter label. The four PSD
profiles are not being treated as certified bounded cusp links.

The [all-radius central-gap theorem](fresh_h20_balanced_endpoint_20261014_central_gap.md)
shows that every endpoint cap is strictly separated from the radius-two
central cap and leaves a common strict escape, at both unit and golden
third-factor scales. A [separate implementation](fresh_h20_balanced_independent_20261014.md)
rederives all 723 reduced trials, all 392 transports, all 82 radius bounds,
and the final escape inequalities. Generation and saved replay pass.

The [bridge construction](fresh_h20_balanced_20261014_bridge.md) adds a
mixed-profile contact motif to the central cap and complete radius-sqrt(28)
vertex layer. Its face-center cap has radius 16/11 and scaling 3/11. It
meets two radius-3/2 radial caps orthogonally and the central cap at pi/3.
All 418 walls are actual facets. The
[independent matrix checker](../scripts/fresh_h20_balanced_bridge_review_20261014.py)
and [saved result](fresh_h20_balanced_bridge_review_20261014.json) confirm all
87,153 pairs, every actual-facet witness, all 392 finite endpoints, and a
strict norm-one escape. Both generation and saved replay pass.

The same bridge report proves a useful unrestricted placement bound:
every cap with center strictly inside the unit cusp cross-section and finite
or parallel contact with the central cap lies within the Euclidean ball of
radius five. The point `2z/5`, of norm 28/5, lies outside every such ball.
Caps whose centers lie on cross-section faces and caps separated from the
central cap remain outside this theorem.

The [sharper endpoint theorem](fresh_h20_balanced_half_20261014.md) proves
that `2z/5` has ball power at least `164/25` against every endpoint cap.
It improves the radius bounds by excluding impossible finite angles through
the algebraic integrality of `2cos(pi/m)`, without assuming arithmeticity
of the polytope. Therefore endpoint caps together with **all** interior-center
central-contact bridges still leave the same strict norm-1/4 escape.
The [separate radial checker](../scripts/fresh_h20_balanced_radial_review_20261014.py)
and [record](fresh_h20_balanced_radial_review_20261014.json) recompute all
16,072 power bounds, assuming the sharper radius bounds; the all-label
proof itself receives source review rather than a second classification.
The secondary half-probe compatibility theorem is not needed for this
volume obstruction.

The [endpoint and central-radius scout report](fresh_h20_balanced_20261014_endpoint_scout.md)
archives the two numerical endpoint banks and a changed-central-radius
search. None supplies a completion. The bridge report separately archives
literal face-center, scaled-center and boundary-profile scouts, preserving
the distinction between a bounded numerical search and an exact exclusion.

For compact H9, wedging the known Bugaenko H8 incidence type and then making
two codimension-2 through codimension-5 face truncations yields 73,513 ordered
recipe orbits. Of these, 5,986 pass the missing-face-size gate and collapse
to 2,953 distinct final nerves. These are polytopal incidence candidates;
a Lorentzian Coxeter realization remains a separate condition. The
[compact search record](fresh_h9_bugaenko_surgery_20261014.md) distinguishes
complete local exclusions, bounded rank searches, unresolved limits, and
numerical Gram fitting. No valid metric survives in the tested proposals.
Eleven full incidence types are excluded; 2,942 remain unresolved. Independent
incidence reconstruction checks 432,736 vertex occurrences. The independent
rank checker verifies the first 217 completed metric-label assignments and
2,358 archived minor records accounting for 30,162 rejected category
extensions. The one modular-rank survivor fails the negative-index condition.
These are exclusions of the displayed assignments, not of their entire
incidence types.

The [archive audit](balanced_cusp_and_bugaenko_surgery_2026_10_14_audit.json)
checks source hashes and local links. Its scope is provenance, not a further
mathematical proof.
