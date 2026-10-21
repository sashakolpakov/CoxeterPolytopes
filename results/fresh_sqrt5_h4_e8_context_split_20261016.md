# H4 compensation distinguished by E8 presence

The [source](../scripts/fresh_sqrt5_h4_e8_context_split_20261016.py) and
[certificate](fresh_sqrt5_h4_e8_context_split_20261016.json) add 108 independent
H4 compensation coefficients, for 32,554 variables total. They distinguish
whether the full finite vertex forest contains E8, at every feasible H4 count
1–3. All six ordinary H3 cross families retain their actual-count split.

The four-triangle [285/8 barrier](fresh_sqrt5_h4_triangle_compensation_barrier_20261016.md)
is explicitly separated: its triangles retain the E8-absent coefficients,
whereas the marked A1+3H4+E8 vertex uses the new coefficients. All 18,035
inherited vertex cuts were rebuilt literally and compared with additive
component sums; 63 profiles changed.

Common external E8 counts are enumerated across all corners of each face.
The residual-forest oracle tracks E8 presence, including the marked component.
All 996 triangle states and 37 paid-quadrilateral stopping cases were rebuilt,
giving 2,013,264 triangle and 140,193 paid-quad marked contexts. The 1,979
inherited baseline-zero quadrilaterals have no H4 compensation and remain
unchanged. There are 109,311 distinct face profiles.

The fixed-35 finite-cut problem is still infeasible. Exact weak duality proves

\[
T\ge\frac{552073685102341}{14608999970782}\approx37.7899709909.
\]

This is a lower bound, not an exact optimum. The numerical finite-cut minimum
is approximately 37.8181060763, with adjacent discount approximately
0.2398617037. The exact dual contains 317 vertex, 204 triangle, 995 quad, and
one distance-two-zero terms. Generation and fresh saved replay both pass.
No full maximizing-forest separation was needed to establish infeasibility.

A separate preflight tied every new coefficient to its old value and compared
the old and split exact forest oracles at H4 counts 1, 2, and 3, all available
facet modes, and three deterministic weight vectors. All 27 old maxima agreed
with the corresponding 54 split maxima; every split winner was recounted by
literal root pairs. This tests the revised state accounting, not full rule
feasibility.

The dominant dual vertices now include marked 2A1+3H4+D7,
A1+2H4+D4+E8, and 2A1+H4+E8+E7. Splitting E8 presence resolves one concrete
coupling but does not address all compensation costs. A numerical pruning
scout did not produce another small transparent barrier; it is not part of
this certificate.

This result excludes the stated payment family at 35. Complete newly affected
paid-quadrilateral extension coverage remains open, and no geometric
dimension-bound claim follows. The frozen compact arithmetic bound-21 proof
candidate remains unchanged.
