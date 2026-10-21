# H9: a compact tetrahedral face in a 94-wall chamber

Four exactly compatible caps extend the [90-wall chamber](fresh_h9_crossplane_repair_20261013.md) to **94 actual facets** and cut its saved gap. They create an actual compact Coxeter **H3 tetrahedral face**, with diagram `3–5–3`, four finite vertices and six finite edges. This supplies finite corners involving `crosspair0`, which previously had none.

The entire H9 chamber still has **infinite volume**. An exact strict future spacelike point survives all 94 walls. The complete frontier of the four new vertices has 24 further open rays. No compact H9 example or new finite-volume dimension is claimed; the lower-dimensional tetrahedron is not claimed to be a new Coxeter type.

## Construction and metric constraints

All normals have square two. The old spherical block

`S = [crosspair0, seed1, seed2, seed6, wall17]`

has type `A2+3A1` and rank five. Project `nextH4` orthogonally to `S-perp`, obtaining `a`. Its exact square is

`q = (4−2φ)/3 > 0`.

In `S-perp ∩ a-perp`, use the projected common center, projected old escape, and coordinate vectors 0 and 3 to obtain an orthonormal frame `e0,e1,e2,e3` with Gram `diag(−1,1,1,1)`. Gram-Schmidt divisions and square roots are certified by strict rational interval bounds.

Choose the rational boost `k=13/4`, so `ct=(k+1/k)/2=185/104` and `st=(k−1/k)/2=153/104`. Put

```text
P = ct e0 + st e1,
B = st e0 + ct e1,
nu = 1/q − 2 = (3φ−1)/2,
r0 = −a/q + sqrt(nu) P,
r1 = sqrt(2/3) B + 2/sqrt(3) e2,
r2 = sqrt(2/3) B − 1/sqrt(3) e2 + e3,
r3 = sqrt(2/3) B − 1/sqrt(3) e2 − e3.
```

The exact frame identities give four norm-two, mutually orthogonal normals. All are perpendicular to `S`. The first has pairing −1 with `nextH4`, hence angle `pi/3`; the other three are perpendicular to it. All other new/old relations are strictly dotted. The 366 new pair checks consist of 29 right angles, one `pi/3` angle and 336 dotted pairs. The first new cap cuts the old strict escape; the other three supply the finite face.

The common center is strict for every new normal. Orthogonal projection of that center onto any wall remains timelike and strictly satisfies every other inequality, because all off-diagonal wall pairings are nonpositive. Thus all 94 inequalities are actual facets.

The [exact recipe](../data/enumerations/fresh_h9_compact_next_20261013_caps.json) fixes every source, frame choice and parameter. The [producer-free replay](../data/enumerations/fresh_h9_compact_next_20261013_verify.json) checks identities through the construction and strict signs through outward integer intervals. A small residual containing zero is not used as a proof of equality.

## The compact face and its complete new frontier

The six common facets

`[seed1, seed2, wall17, r1, r2, r3]`

have type `6A1`. They are perpendicular to all four bounding normals

`[seed6, crosspair0, nextH4, r0]`.

Those four normals have the exact norm-two Gram diagram `3–5–3`. Its determinant is `5−4φ<0`; each three-node principal subdiagram is positive definite, and every entry of its inverse is strictly negative. Its four future timelike dual rays generate its entire feasible cone. Therefore its projectivization is a compact tetrahedron, not just a collection of sampled finite vertices.

Each tetrahedral vertex is checked by its full rank-nine physical Gram, timelike norm, future orientation, and strict inequalities against every nonincident wall of the 94-wall chamber. Thus no other wall truncates the tetrahedron. Two vertex diagrams have type `H3+6A1`; the other two have type `A2+7A1`.

Every positive-semidefinite subdiagram meeting a new wall must use old walls from `S ∪ {nextH4}`, because every other old relation is dotted. All 960 such subsets are checked: 896 are positive definite and 64 are indefinite. None is positive affine. The earlier absence of affine circuits is inherited from the replayed 90-wall certificate. This census also proves that exactly four finite vertices involve the new walls.

The four new vertices have 30 distinct outgoing rank-eight spherical edges. Testing every possible ninth facet proves that six edges close inside the tetrahedron and 24 remain open to the light cone, with no ideal vertex at their other ends. This is the complete frontier of these four vertices, not a census of every inherited vertex. The 24 archived edge supports are concrete endpoint constraints for a future attempt to join this component to the rest of the chamber.

## Whole-chamber volume gate

The exact saved dual vector `y` has square

`(4764512044018601985931865 − 2484892311721629081203273 φ) / 16000000000000000000000000 > 0`,

approximately `0.0464919891`, and future-center pairing `−558547622989/250000000000`. All 94 dual wall inequalities are strictly positive. This certifies an open spacelike patch in the feasible projective cone and proves infinite volume.

Numerical whole-cone enumeration proposed this witness. It found 68,511 projective vertices, including 68,485 spacelike ones, and maximum square approximately `0.0464946261`. Those counts are diagnostics, not exact volume claims. The unperturbed proposal has support `[0,2,3,4,5,8,9]` and active walls `[seed1,seed6,seed7,wall16,wall17,wall28,wall67,wall70,coupled0]`; the exact archived witness is strictly inside every inequality.

## Search scope and validation

The initial [single-anchor corner scout](../data/enumerations/fresh_h9_compact_next_20261013_scout.json) tested 204 coordinate-derived frames and 251,648 real boost solutions, producing 458 numerical four-cap entries. These were proposals, not an exact cap census. The subsequent [finite-edge scout](../data/enumerations/fresh_h9_compact_next_20261013_edge_scout.json) used two spherical rank-five blocks, their old spherical extensions, specified angle labels, 600 boosts and 180 triangle phases. It produced 555 numerical entries. The selected `nextH4/pi3` construction was exactified separately.

A preliminary `seed4/pi6` profile was discarded because its prescribed `3–6` three-node path is positive affine. Larger-label alternatives and the unselected numerical entries receive no exact validity or finite-volume claim here. Neither bounded scout excludes arbitrary caps or coherent replacements.

```sh
python scripts/fresh_h9_compact_next_20261013_verify.py --verify
python scripts/fresh_h9_compact_next_20261013_regressions.py
```

The [rejection tests](fresh_h9_compact_next_20261013_regressions.json) reject a timelike replacement escape and the previous escape that the new cap cuts. The verifier uses the shared exact `Q(φ)` and rigorous interval machinery; it imports no numerical producer and runs no numerical search. Both the [separate symbolic checks and nested 160-digit replay](fresh_h9_compact_next_20261013_review.md) and a [second source/formula review](fresh_h9_compact_next_20261013_source_review.md) pass; both explicitly reuse the existing interval implementation.
