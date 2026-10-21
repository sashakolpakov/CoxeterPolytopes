# A deficient cusp of the radius-four replacement

The [radius-four replacement](fresh_h20_intrinsic_caps_20261014_radius4_replacement.md) creates an explicitly unbounded cusp link, in addition to its saved strict escape. The [small exact checker](../scripts/fresh_h20_intrinsic_caps_20261014_parallel_cusp_gap.py) and [certificate](fresh_h20_intrinsic_caps_20261014_parallel_cusp_gap.json) bind the frozen construction and inspect every one of its 587 walls.

Let r be the new norm16 cap at global index23, and let s be old source wall22, global index15, of norm1. Their pairing is−4. Thus

```text
delta=r+4s,
delta²=0,
(delta,strict_center)=-709/6-2phi < 0.
```

This is a nonzero future null vector. Every wall pairs nonpositively with it. Exactly seventeen inequalities are incident, at indices

```text
2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,19,23.
```

All other 570 pairings are strictly negative. The incident Gram is rational and has exact rank16. This is the rank of the Euclidean projected normals in `delta-perp / R*delta`, whose dimension is19; it is not a claim about ambient coordinate rank.

Consequently the common orthogonal complement of the incident projected normals has dimension3. Both signs of each such vector are directions of recession of the horospherical link. The saved strict timelike point supplies a strict interior point of that link, and the 570 nonincident walls impose no inequalities sufficiently deep in this cusp. Therefore the actual link has **three-dimensional lineality and is unbounded**.

This diagnoses one cusp only. A finite-volume continuation retaining this parallel pair must add incident constraints that remove this defect, or replace the responsible walls; a complete cusp census is not claimed.

```sh
python3 scripts/fresh_h20_intrinsic_caps_20261014_parallel_cusp_gap.py --verify
```

Generation and saved replay pass. The frozen replacement source and certificate were not changed.
