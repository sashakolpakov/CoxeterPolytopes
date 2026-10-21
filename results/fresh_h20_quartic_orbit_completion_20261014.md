# Changed-chamber construction with quartic cusp sections

This was a positive construction attempt outside the now-excluded cusp sections of the original chamber. It retains a known valid projected cusp link but allows its other facets to come from a larger ambient reflection orbit. No finite-volume polytope was obtained. The results here are numerical necessary-condition tests, not an exact exclusion of all these normals or all reflected mirrors.

The [producer](../scripts/fresh_h20_quartic_orbit_completion_20261014.py) uses each of the 44 quartic records in the [exact local census](fresh_h20_rank3_omitted_sections_20261014.md), at shift parameters `u=1/10, 1/4, 1/2, 3/4, 9/10`. For each of these 220 sections, it keeps the cusp-incident original walls except the designated omitted inequality. Their projected Euclidean link is already certified by the local census. The other 210 source inequalities are not imposed.

The larger mirror pool contains the original roots and their reflections in every original simple mirror. After orientation and exact integral-vector deduplication it has **29,043 roots**: 2,163 of norm two and 26,880 of norm four. The [compressed bank](../data/enumerations/fresh_h20_quartic_orbit_completion_20261014_bank.npz) records each root and its original-wall/reflection recipe. This is one reflection from the original wall set, not the full ambient reflection orbit.

Only roots strictly negative on the chosen cusp are possible closing caps. The scout projects them into the section, requires a positive projected square, and tests compatibility with every retained cusp wall. It allows divergent pairs as well as Coxeter angles. For the fixed quartic normal, every squared cosine belongs to its quartic field. A finite angle `pi/m` therefore has `phi(m)/2` dividing four, giving the necessary alphabet

```
2, 3, 4, 5, 6, 8, 10, 12, 15, 16, 20, 24, 30.
```

This is a superset of the individual field's actual cyclotomic subfields. No arithmeticity condition is imposed on the proposed H20 chamber. In particular, the original-chamber theorem restricting labels to `{2,3,4,6}` does not apply after this change of walls.

The [numerical compatibility archive](fresh_h20_quartic_orbit_completion_20261014.json) records 89,242 acute cap proposals, of which 8,136 pass the numerical Coxeter-angle test; another 96 are retained as near-angle remainders. These counts have multiplicity across sections and do not count inequivalent mirrors. Compatibility with the cusp walls is only a local necessary condition.

The [edge-coverage test](../scripts/fresh_h20_quartic_orbit_completion_20261014_edges.py) then enumerates every edge leaving the fixed cusp: omit one facet from each irreducible affine component of its link. The nineteen remaining projected normals have a positive-definite Gram A. If a candidate cap has projected norm N and pairings b with these walls, a necessary condition for a finite or ideal endpoint is

```
N - b^T A^(-1) b >= 0.
```

For a finite endpoint the inequality is strict; the zero case is retained as a possible ideal endpoint, without claiming that its full cusp link is bounded. This Schur-complement test avoids mistaking pairwise angle compatibility for global closure. It is applied to every proposed cap, including the 96 near-angle remainders.

Across all 220 sections there are **42,420 tested cusp edges, and none has a proposed endpoint in this bank**. The [coverage archive](fresh_h20_quartic_orbit_completion_20261014_edges.json) records every omitted-facet tuple. Consequently this attempt produces no chamber to submit to an exact finite-volume verifier.

Both scouts use double-precision geometry after evaluating the algebraic normals at high precision. Their thresholds, finite mirror pool, and five parameter choices prevent an exact or continuous exclusion claim. The next useful question is whether the integral norm-two/norm-four condition itself forbids endpoint caps for a particular quartic normal; merely enlarging this mirror pool without examining that obstruction would be poorly targeted.

The [separate review](fresh_h20_quartic_orbit_completion_20261014_review.md) passes. Its independent integer checker reconstructs all bank vectors and reflection recipes, verifies the quartic degree alphabet, and exhausts the archived edge tuples. It does not independently certify the floating-point endpoint signs.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_quartic_orbit_completion_20261014.py
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_quartic_orbit_completion_20261014_edges.py
```

These commands regenerate the numerical archives. They are not exact saved-certificate replays.
