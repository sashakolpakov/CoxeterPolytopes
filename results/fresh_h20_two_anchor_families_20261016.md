# Continuous two-anchor H20 cap diagnostics

The saved 158-wall chamber admits new numerical cap families missed by the older four-contact banks. These diagnostics recover 16 clusters with numerical Coxeter contacts, but none closes the whole tracked H3 face. No new exact wall, finite-volume chamber, or all-parameter exclusion is claimed.

All six source/data pairs below passed generation and a fresh saved replay. Their direct source manifests and Python syntax pass. The entire prefix is frozen after this report. The replays reproduce the same numerical implementation; they are not independent exact arithmetic reviews.

## Parameter domain and its limits

The [base source](../scripts/fresh_h20_two_anchor_families_20261016.py) and [data](fresh_h20_two_anchor_families_20261016.json) use the frozen 158-wall geometry and keep S16 orthogonal. They prescribe contacts to wall 151 and two distinct walls from `{31,95,132,142,145,150,152,156,157}`. The 151 label is in `{2,3,4,5,6}` and the other labels are in `{2,3,4,5,6,8,10}`: 8,820 prescribed patterns.

For an independent three-root frame, the prescribed projection is fixed and the remaining norm equation is a conic in its two-dimensional orthogonal complement. Circle, spacelike-hyperbola, timelike-hyperbola, and null-line cases are represented by a quadratic numerator divided by a quadratic denominator. Every variable old-wall inequality becomes a quadratic inequality in that parameter. The implementation splits at all computed real roots on the full real line or positive half-line, with no artificial boost cutoff. However, its coefficients, root isolation and sign tests are floating point. The output is **numerical interval evidence**, not a certified exact interval theorem.

Every unprescribed variable old contact is required to be dotted or parallel. Constant contacts are checked separately. Isolated extra finite-contact parameter values therefore lie outside this family domain. Samples are chosen inside the surviving intervals; their rational parameter strings refer to numerically normalized kernel bases and do not yet define exact algebraic roots.

The result is 2,730 one-interval families: 1,503 circles, 1,133 timelike hyperbolas, 71 spacelike hyperbolas and 23 null lines. There are 41,984 distinct sampled candidates, each rechecked against all 158 numerical angle gates, its norm and the center orientation. The inherited center is timelike and weak: its zero pairings are exactly S16. Every cap is strict against it, so a small perturbation into the S16 chamber supplies a common strict center for a fixed finite candidate set. The tracked initial H3 face has seven finite and eight exterior vertices, with maximum squared radius about 1.00780360535 in the fixed chart.

## Pair and third-cap tests

The [pair source](../scripts/fresh_h20_two_anchor_families_20261016_pairs.py) and [data](fresh_h20_two_anchor_families_20261016_pairs.json) select 160 first samples and solve prescribed mutual labels `{2,3,4,5,6,8,10}` against every family. Of 5,577,796 quadratic roots, 5,650 pass the strong pair filter; 219 coverage/label representatives are retained. Sixty-four full numerical face censuses were performed. The best has seven finite and nine exterior vertices and maximum squared radius about 1.00273667200.

The [all-dotted third-cap test](../scripts/fresh_h20_two_anchor_families_20261016_clusters.py), with [data](fresh_h20_two_anchor_families_20261016_clusters.json), intersects each family with mutual dotted inequalities and the uncovered saved sample inequalities for 32 selected pairs. It finds no covering third interval in that specified domain. The [finite-angle third-cap test](../scripts/fresh_h20_two_anchor_families_20261016_finite_third.py), with [data](fresh_h20_two_anchor_families_20261016_finite_third.json), checks 2,273,510 quadratic roots; 270,012 pass the new-pair gates, but none covers the remaining saved samples. These statements concern selected pairs and the listed family/contact domain only.

## Acute relaxation and coordinated angles

The [acute relaxation](../scripts/fresh_h20_two_anchor_families_20261016_acute_cover.py), with [data](fresh_h20_two_anchor_families_20261016_acute_cover.json), uses all 42,203 samples and solved partners. Old contacts remain Coxeter, while new-new contacts need only be acute. A sampled clique cover is recomputed after adding an ideal sample from each remaining exterior region. It reaches its explicit **100-iteration guard**; it neither proves infeasibility nor finds a closed face. The best step, 62, uses 13 caps and has maximum squared radius about 1.00000477444, but 22 new-new pairs are not Coxeter. Such a cluster is a geometric relaxation, not a Coxeter construction.

The [coordinated-angle source](../scripts/fresh_h20_two_anchor_families_20261016_coordinated_angles.py) and [data](fresh_h20_two_anchor_families_20261016_coordinated_angles.json) vary every family parameter simultaneously in the first four acute covers. They try the two closest target angles for each initially non-Coxeter pair, preserve initially finite contacts, and retain dotted inequalities and all 21 initial sample cuts. Each least-squares run is limited to 1,200 evaluations. All 136 target assignments are saved; 16 pass the numerical old/new angle and norm gates. The best recovered cluster has maximum squared radius about 1.00046022782, nine finite, fourteen exterior and one near-sphere vertex. No recovered cluster closes the face. The near-sphere class is deliberately unresolved numerically.

The fresh saved replays of the acute and coordinated stages ended successfully in sessions 19417 and 55479. Their respective JSON hashes are `3c8aa4406f8c467bf15d873750aae2727347184cd2bca4b8ab0704cbc530f6aa` and `ee26a0aabf852f240e8d0f852391b94a5c247ca1a12acd499a7513d6ea43d328`.

## Next finite target

The promising unresolved target is incidence-preserving angle recovery in a nearly closed acute cluster. Intrinsic finite adjacency must use normals projected off wall 151, not their ambient norms. A compact polygonal new facet requires at least three finite neighboring side planes. Another improvement is to split each allowed family interval at every current ideal-sample pairing zero, preserving narrow coverage windows rather than selecting a few high-coverage samples. Both tasks belong to fresh prefixes; this bounded diagnostic is frozen.
