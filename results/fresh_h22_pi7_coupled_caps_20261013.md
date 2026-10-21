# Coupled H22 endpoint caps and an actual compact quadrilateral

Two new exactly verified caps retain the [837-wall branch](fresh_h22_pi7_cusp_pair_20261013.json), meet one of its pi/7 caps at right angles, and close two previously certified open edges. Together they give **839 actual facets**, three new finite vertices and an actual compact H2 quadrilateral with angles **pi/4, pi/2, pi/2, pi/2**. The two repaired cusp links remain bounded.

The chamber nevertheless has **infinite volume**. Both the previous rational escape and a simpler integer escape of square1586 remain strict against every839 wall. A specified old cusp is also untouched. This is a coherent endpoint/face construction, not a new finite-volume H22 example or an unrestricted obstruction to further completion.

The [exact certificate](fresh_h22_pi7_coupled_caps_20261013.json) and [checker](../scripts/fresh_h22_pi7_coupled_caps_20261013.py) contain the formulas, all1,675 changed pair gates, four actual vertex points, two preserved cusp links and both global escapes.

## Construction from two edges at one finite vertex

Use vertex15 of the saved835-wall geometry. Its21 rational incident walls are

`S=[4,6,9,11,14,16,18,21,22,23,24,25,27,28,40,46,48,49,50,63,196]`.

The22nd incident wall is829, namely pi/7 cap17, candidate246. It is right to every wall in S. The rational Gram on S is positive definite. The [exact frontier certificate](fresh_h22_pi7_edgewalk_20261013.json) proves that releasing27 or9 gives unclosed edges; their saved row indices are342 and332.

Let p be the strict center projected perpendicular to all22 incident normals. It is future timelike. Define

`Z=sqrt(2/(-p²))*p`, so `Z²=-2`.

Let u and v be the rational Gram-dual vectors for released walls27 and9 respectively. Exact inversion gives

`u²=2`, `v²=1`, `(u,v)=1`, and `(Z,u)=(Z,v)=0`.

The new normals are

`rA=(209Z-241u)/120`,

`rB=(241Z-418v)/sqrt29281`.

Their norms are exactly two and their mutual pairing is exactly zero. These statements follow from the displayed frame Gram, not from numerical interval overlap. The first cap has hyperbolic boost `e^t=15/4`; the second boost is solved from its prescribed right angle with the first.

Each cap is right to the21 retained normals of its old edge, including cap829. Every other pairing against the837 retained walls is strictly dotted. The total changed-pair count is43 right and1,632 dotted. In particular, both new caps are strongly dotted to the two reflected pi/7 caps that completed the prior cusp links.

The old null edge limits can be scaled to `Z-u` and `Z-sqrt2*v`. The respective new-cap pairings are

`8/15 > 0`, and `(-482+418sqrt2)/sqrt29281 > 0`.

Thus these are actual cutters of the old unclosed rays, rather than additions dotted to every pi/7 cap or caps that merely increase the wall count.

## The compact face and exact geometry

Fix the common20 walls in `S union {829}` after deleting27 and9. In their H2 intersection, the four vertices have frame coordinates `(Z,u,v)`

| Vertex | Coordinates | Square |
|---|---|---:|
|Old corner|(1,0,0)|−2|
|First endpoint|(241,−209,0)|−28800|
|Joint endpoint|(50369,−29281,−28800)|−843292800|
|Second endpoint|(209,0,−241)|−29281|

Every point is future timelike, has exactly its stated22 zero wall inequalities, and is strict on all other839 walls. The corresponding incident Grams are positive definite: the edge blocks are principal spherical blocks of the old corner, and each new normal is perpendicular to its block. At the joint vertex, both new normals are mutually orthogonal and perpendicular to the common20-wall block.

The inverse dual minor is `[[1,-1],[-1,2]]`, giving intrinsic angle pi/4 between the old side walls. All three remaining angles are right. The four points and sides form a compact quadrilateral. Every other wall is strict at all its vertices, hence cannot truncate its convex hull. This proves an actual compact H2 face in the retained chamber, rather than merely a formal diagram or a numerically suggested cycle.

The common strict timelike center and global nonpositive pairings also prove all839 inequalities are actual facets: projection to any wall remains timelike and strict on every other wall.

No new positive affine subdiagram is introduced. Any subset containing only one new cap and no dotted pair consists of that cap and a principal positive-definite edge block, all orthogonal to the cap. A subset containing both new caps must use only their common20-wall block; its Gram again splits into positive-definite blocks. This statement concerns newly introduced diagrams; inherited ideal vertices and bounded cusp links remain present, so it is not a compactness claim.

Both new caps pair strictly negatively with the two previously repaired cusp rays. Their incident sets and bounded-link certificates are therefore unchanged. The report does not claim a complete global839-wall vertex or cusp census.

## Remaining volume obstruction

The previous strict rational escape still has square

`3956990142533/12500000000000 > 0`.

All839 of its strict inequalities are replayed. The additional [uncapped-cusp certificate](fresh_h22_incidence_prism_20261013_uncapped.json) supplies the old null ray

`delta=(20,-47,-32,-63,-93,-76,-58,-40,-20,80,52,106,160,135,111,87,64,-20,-40,-29,-16,42)`.

Both new normals pair strictly negatively with `(delta,0)`, approximately−48.1830 and−41.3789, so this old cusp receives no new incident wall. In the coordinate convention `diag(Q,2)`, the simple vector `(100delta+h,40)` has square1586 and center pairing−814. The checker independently verifies that it is strict against every839 wall; in particular its pairings with the two new caps are approximately−4573.19 and−3945.31. Thus the global failure of finite volume is exact and persists despite the new compact face.

## Search and verification scope

A bounded numerical scout first tested all-right endpoint insertions on the475 certified strict open edges that retain a pi/7 cap. Powers-of-two boosts with exponents between1 and20 supplied compatible proposals on each. None cut the old remote escape. A second scout coupled two such edge families at one vertex, using a rational first boost and labels2,3,4,5,6,7,8,10,12,16; it saved11,601 proposals across4,670 edge pairs. These are proposal counts, not an exact census of distinct mirrors or all possible couplings. Only the two selected normals above are claim-bearing.

The exact checker imports no scout, optimizer or floating-point calculation. Rational projection identities prove every required equality. Outward-rounded integer radical intervals prove strict signs, using the frozen predecessor certificates as explicit dependencies. The [independent source/formula review and fresh replay](fresh_h22_pi7_coupled_caps_20261013_review.md) passed, including both escape certificates and the direct edge-limit cuts. That review uses the same rigorous interval engine; it is not a second implementation of all839 roots.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_pi7_coupled_caps_20261013.py --verify
```

Further constructive work must address another remaining gap or the unchanged cusp; no additional cap search was run after freezing this839-wall checkpoint.
