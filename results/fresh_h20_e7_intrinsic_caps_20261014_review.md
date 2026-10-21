# Independent review of the 85-wall H20 replacement

The [reviewer](../scripts/fresh_h20_e7_intrinsic_caps_20261014_review.py) and
[certificate](fresh_h20_e7_intrinsic_caps_20261014_review.json) independently
reconstruct the two transferred caps, the four removed old walls, and all
85 resulting normals. Sparse square-class multiplication with rational
interval signs supplies a separate arithmetic implementation from the
producer's recursive quadratic arithmetic.

All 3,570 wall pairs pass: 1,336 right, 151 at pi/3, four at pi/4, 73
parallel, and 2,006 divergent. The checker verifies a strict interior point
and a strict timelike point on each of the 85 facets.

An explicit isometry transfers all 46 vertices of the independently
verified H13 source. Exactly 44 are actual finite ambient vertices. The
other two are ideal rays with bounded intrinsic links but an explicit
spacelike recession in each ambient link. The other 63 old nonincident
walls have future timelike projections and impose no further boundaries
on this complete face. Intrinsic completeness inherits the
[separate sparse-support source census](fresh_h13_a11a2_completion_20261014_review.md).

The original exact spacelike escape is strict against all 85 walls, proving
infinite ambient volume. This review does not claim a census of every
ambient vertex or repeat all rejected fixed-87 source orientations.

Replay: `python3 scripts/fresh_h20_e7_intrinsic_caps_20261014_review.py --verify`.
