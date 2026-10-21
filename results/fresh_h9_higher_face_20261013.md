# A coupled roof over the compact H4 prism

This construction retains the previous106 H9 facets and adds12 exactly certified facets. Its distinguished H5 face now has19 facets and25 finite vertices in one connected component. All eight original edges leaving the compact H4 prism acquire finite endpoints. The H5 face has41 closed finite edges and43 remaining open rays. **Neither this H5 face nor the118-wall H9 chamber has finite volume:** the final certificate contains strict spacelike escapes for both.

The positive result is a linked roof and subsequent truncations, with finite metric couplings between newly constructed walls. This advances beyond isolated endpoint insertions and beyond the previously excluded single-top prism construction. It does not give a new finite-volume example or a completion obstruction for arbitrary further caps.

## Construction

The retained compact H4 prism has tetrahedral sides

`T=(seed6, crosspair0, nextH4, compact0)`

and bases `wall17, facePrism`. Its five common orthogonal walls are `seed1, seed2, compact1, compact2, compact3`. Freeing `seed1` leaves the actual H5 face fixed by the other four walls. Every new normal is perpendicular to those four; consequently the resulting19-facet H5 face is itself Coxeter, with the same physical angle labels.

Number the original prism corners by `2k+b`, where the kth side of T is omitted and `b=0,1` chooses `wall17,facePrism`. Seven roof normals are determined as follows. An unlisted pairing with one of the eight edge walls is zero.

| New wall | Original corner edge | Nonright edge pairing | Additional metric contact |
| --- | ---: | --- | --- |
| roof0 | 1 | none | π/16 to wall17 |
| roof1 | 2 | none | π/12 to roof0 |
| roof2 | 3 | none | π/2 to roof0 |
| roof3 | 5 | none | π/3 to roof0 |
| roof4 | 4 | none | π/3 to roof1 |
| roof5 | 7 | π/3 to seed6 | π/5 to roof0 |
| roof6 | 6 | π/3 to seed6 | π/4 to roof5 |

The first wall also closes the edge at corner0. It has a dotted pairing with seed6, so it does not satisfy the all-finite-side hypothesis of the [single-top obstruction](fresh_h9_prism_lift_20261013.md). Several caps and finite new–new contacts are essential here. A trial π/3 contact between roof5 and roof6 was rejected because it creates the affine3–3–3 triangle with seed6; the retained π/4 choice avoids that circuit.

The roof formulas use the unit future timelike vertex direction `v`, the unit spacelike freed-wall direction `u`, and a fixed projection `p` to the eight edge normals. Put `q=2-p²>0`. For an anchor a define

`A=−(v,a), B=−(u,a), Δ=A²−B²>0, z=((p,a)+2cos(π/m))/sqrt(q)`.

Then

`s=(Az+B sqrt(Δ+z²))/Δ`,

`c=(Bz+A sqrt(Δ+z²))/Δ`,

`r=p+sqrt(q)(s v−c u)`.

The identities `c²−s²=1` and `As−Bc=z` prove the norm and prescribed finite angle exactly. Bounds on Δ select a real, oriented branch. The π/12 and π/16 values are explicit nested radicals; there is no assumed small-label bound. The [roof certificate](../data/enumerations/fresh_h9_higher_face_20261013_verify.json) establishes12 new actual finite vertices and all eight old vertical edge closures.

Three further caps prescribe complete finite pairing profiles at remote hyperideal vertices of this roof. In each row the cap is also right to the four fixed H5 walls:

| Cap | Five other prescribed walls | Angle labels |
| --- | --- | --- |
| adaptive0 | nextH4, crosspair0, roof0, roof1, roof5 | 2,3,2,3,4 |
| adaptive1 | wall17, crosspair0, roof0, roof1, adaptive0 | 5,2,5,3,2 |
| secondAdaptive0 | seed6, crosspair0, roof5, roof6, adaptive0 | 3,2,5,8,8 |

These caps solve nine linear pairing conditions and a one-dimensional norm equation. If E is the matrix of incident normals, `M=EGEᵀ`, and b is the desired covector, the particular vector is `p=bM⁻¹E`. Project the common center perpendicular to E to obtain z. The chosen root is

`r=p−sqrt((2-p²)/z²) z`.

Certified nonzero elimination pivots justify the inverse. Orthogonality supplies the exact norm identity, rather than an interval overlap being treated as equality. The [full-face certificate](../data/enumerations/fresh_h9_higher_face_20261013_truncation_verify.json) verifies these three caps. They cut remote gaps while preserving the roof vertices, but create no new finite vertices themselves.

A polar truncation, `roofPolar`, is perpendicular to the nine normals

`seed2, seed6, nextH4, crosspair0, compact1, compact2, compact3, roof5, roof6`.

Their orthogonal complement is positive and one-dimensional. Projecting the common center gives z and the outward root is `−sqrt(2/z²)z`. This creates three finite vertices.

Finally `roofMetric` is right to

`seed2, wall17, crosspair0, compact1, compact2, compact3, roof1, roof4`

and meets compact0 at π/3. The positive two-dimensional complement reduces this to one linear condition and a quadratic norm equation. The selected branch creates two more finite vertices. Both final normals retain all preceding walls.

## Exact geometry and its limits

The [final certificate](../data/enumerations/fresh_h9_higher_face_20261013_final_verify.json) checks all1,338 changed pairs. Exact projection identities supply norm2 and the prescribed finite equalities; outward-rounded integer intervals certify strict dotted comparisons and orientation. A common strict timelike interior together with nonpositive pairings makes every retained inequality a true facet: projecting that interior onto any one wall remains timelike and makes every other inequality strict.

All12 additions have nine finite neighbors when inserted. Their complete6,144 finite-neighbor subset tests give1,632 exclusions by a dotted pair,3,952 positive-definite subsets,560 indefinite subsets and no affine subset. Together with the inherited106-wall certificate, this excludes positive affine subdiagrams in this literal118-wall configuration. It does not establish compactness.

For the H5 face, the verifier tests every old/new wall against the four fixed normals. Exactly19 walls meet it as facets; the others have an indefinite Gram obstruction. All19 are orthogonal to the fixed normals. The finite-vertex census checks all `binomial(19,5)=11,628` subsets:11,594 contain a dotted pair,25 are positive definite and9 are indefinite. Every positive case has an explicit future timelike point satisfying all other118-wall inequalities strictly. Thus these are all finite vertices of this face, including the eight unchanged prism vertices. They are simple.

Sharing four H5 walls gives the full finite-edge graph. It has one25-vertex component,41 finite edges and43 open rays; `2·41+43=5·25`. No additional retained facet can supply an unrecorded finite endpoint. The affine census also rules out a ninth-facet-supported ideal endpoint. The open rays still reach the light cone; this does not claim that their ideal limits can never be repaired by new walls.

The exact final strict escape has norm

`(21111078902767917419386970717 − 13047074150472097078498634521 φ) / 4000000000000000000000000`

`≈0.1173683909761>0`.

Its future orientation and every one of118 strict inequalities are checked. Orthogonally projecting it onto the four fixed walls gives another strictly spacelike point, strict on every nonfixed inequality, so the H5 face itself also has infinite volume. This is stronger than a numerical projective-vertex failure and avoids a false compactness conclusion from the finite-vertex census.

The [rational witness](../data/enumerations/fresh_h9_higher_face_20261013_final_escape.json) is derived from the last numerical gap, but its verification does not trust numerical geometry. Its leading unsolved support is `seed6, nextH4, facePrism, roof3, adaptive0` inside the H5 face. This is a concrete different place to seek a metric coupling or a coherent replacement. The current report gives no all-parameter obstruction there.

## Reproduction and review

```sh
python scripts/fresh_h9_higher_face_20261013_verify.py --verify
python scripts/fresh_h9_higher_face_20261013_truncation_verify.py --verify
python scripts/fresh_h9_higher_face_20261013_final_verify.py --verify
```

All three use the existing exact Q(φ) arithmetic and outward-rounded radical intervals, with source-bound predecessor certificates. They import no numerical scout or optimizer. The numerical scripts under the same prefix propose labels and branches only. Their bounded profile searches and incomplete labels are not exclusion theorems.

The normalized-hyperbola formula and the finite-vertex/affine arguments received an independent source review from the arithmetic agent, with no issue found. That review did not repeat all3,584 roof subset tests. The parent agent's [completed review](fresh_h9_higher_face_20261013_review.json) independently checks the symbolic norm, anchor and projection identities and replays the construction at160 decimal places:18,236 rigorous enclosures nest inside the archived80-digit bounds, with12 predecessor-certificate comparisons. The [review script](../scripts/fresh_h9_higher_face_20261013_review.py) uses the same interval implementation; this is source review and higher-precision replay, not a separate implementation of all118 normals.
