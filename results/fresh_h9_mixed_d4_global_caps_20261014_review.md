# Separate review of the dotted-cap example and coverage argument

PASS. The [review checker](../scripts/fresh_h9_mixed_d4_global_caps_20261014_review.py) uses the original construction's quartic scalar implementation, separately from the producer's review-engine arithmetic. Its [certificate](fresh_h9_mixed_d4_global_caps_20261014_review.json) reconstructs the literal interpolation at `t=65323/65536`, verifies all twelve strict dotted inequalities, and rechecks the escape cut, common interior, old facet actuality and timelike projected old normals. The old twelve-wall chamber is inherited.

The [general argument](fresh_h9_mixed_d4_global_caps_20261014.md) was reviewed directly. Future-null interpolation gives a strict escape cutter whose normalized old pairings tend to negative infinity. Its entire H8 facet has a uniform inward collar, proving infinite volume. This explains why counting escape-cutting walls alone gives no completion criterion.

The inward normal-coordinate sign and the resulting ball-radius formula agree with the old inequalities. A strict open-ball cover of a compact perpendicular base gives a uniform finite height bound. Conversely, an open base patch outside all closed balls has unbounded inward height and infinite volume. Closed-ball coverage by itself is not asserted sufficient.

The parallel-contact boundary case also passes: when every old/new side pairing is parallel or dotted, the closed balls lie in the base, with at most finitely many tangency points on each H7 facet. A facet point outside those tangencies has a nearby uncovered interior patch. Thus finite contacts are necessary for this H8 attachment mechanism. Acuteness also puts every old projected center in the base, making containment of their hyperbolic convex hull in the closed-ball union necessary.

This review does not produce a lower-dimensional source placement or independently certify the bounded numerical coverage scout. Generation and saved replay pass for the exact example review.
