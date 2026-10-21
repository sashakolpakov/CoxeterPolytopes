# Independent all-label classification of the 1,068 exterior polars

The [independent checker](../scripts/fresh_nonarithmetic_joint_regions_20261016_polar_gate_review.py) and [certificate](fresh_nonarithmetic_joint_regions_20261016_polar_gate_review.json) pass generation and fresh saved replay. Exactly **16** polar walls of the fixed complete 36-wall cone are Coxeter-compatible; the other **1,052** are excluded for every integer finite angle label. No angle-label cutoff or numerical tolerance enters the proof.

The complete old ray list and rational coordinate enclosures come from the [independent 36-wall cone review](fresh_nonarithmetic_escape35_20261016_cone_review.md). The sixteen survivors, including the recovered label-10 candidate at exterior-ray index 8, are tied to that list and accepted by the [independent batch review](fresh_nonarithmetic_joint_regions_20261016_batch_review.md). Survivor indices and rejected indices form a disjoint exhaustive partition of all 1,068 old exterior rays.

For each rejected ray `V`, the review independently decodes the proposed blocking old wall `n`, recomputes intervals for `V²`, `n²`, and `(V,n)`, and encloses

```
4(V,n)^2 / (V^2 n^2).
```

The producer's suggested wall and consecutive labels are treated only as witnesses to check. None of its numerical ray products, threshold intervals, or exact arithmetic helpers is used. The review proves a strict gap between the independently enclosed product and the two adjacent label values `4cos²(pi/m)` and `4cos²(pi/(m+1))`. These values increase strictly for all integers `m>=2` and approach 4, so a gap between consecutive values excludes every possible finite label.

The trigonometric enclosures use a different identity from the producer:

```
pi = 4*(atan(1/2)+atan(1/3)).
```

Each arctangent has a 240-term exact rational alternating-series bracket. Cosine is enclosed with 90 exact rational alternating terms at the opposite argument endpoints, using its monotonicity on `[0,pi/2]`. The certificate saves independent thresholds and positive rational margins for every rejection.

All sixteen accepted polar walls preserve the independently checked strict spacelike escape. Therefore every coherent subset of this one-generation polar bank remains infinite-volume. This does not exclude nonpolar walls or new polar directions arising after further clipping.
