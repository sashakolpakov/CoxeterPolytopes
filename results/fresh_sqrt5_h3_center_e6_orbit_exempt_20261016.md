# Exempting the E6 edge orbit does not suffice

The [source](../scripts/fresh_sqrt5_h3_center_e6_orbit_exempt_20261016.py) and
[certificate](fresh_sqrt5_h3_center_e6_orbit_exempt_20261016.json) exempt exactly
the two E6 edges joining the trivalent root to its degree-two neighbors from
the adjacent discount. The distance-two discount remains zero; one uniform
adjacent discount, bounded by 2/3, applies elsewhere in components of rank≥6.

All 18,035 inherited vertex cuts were rebuilt by literal unordered root pairs;
1,299 profiles changed, and only their adjacent-discount coefficient changed.
All component/marked costs and all 996 triangle states, 37 paid-quad stopping
rows, and inherited baseline-zero quads were recomputed. Both former killer
quads 4266 and 7746 now pay exactly 1 at every zero-reserved-spectator white
facet mark.

The enlarged family nevertheless has the exact restricted-model lower bound

\[
T\ge\frac{1859938611880960}{48400999854797}
\approx38.42768987130>35.
\]

This is a lower bound, not a certified optimum. The finite-cut numerical
relaxation is approximately 38.43412022101 with adjacent discount approximately
0.11617106569. The exact dual contains 256 vertex terms, 127 triangle terms,
925 quad terms, and the zero-distance-two restriction. Generation and fresh
saved replay pass. There is no geometric nonexistence or dimension-bound claim.

The largest new quad demands come from classes 7661 (E8 paid corner), 7667
(D7 paid corner), and 6156 (E7 paid corner). A concrete remaining defect is that
an inherited H3×A4-end coefficient depends only on H3 presence, while its
individual pair payment divides by three times the actual H3 count. A face
with five H3 copies therefore imposes a large coefficient that is also paid
at vertices with only one H3 copy. The
[small multiplicity barrier](fresh_sqrt5_h3_multiplicity_barrier_20261016.md)
isolates an exact lower bound 4997/139>35 from one such face and two vertices.

Complete affected paid-quadrilateral extension coverage remains open. The
reduced eligibility makes the inherited complete triangle domain a safe
overapproximation; no old paid-face pruning is silently reused.
