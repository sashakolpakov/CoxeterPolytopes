# A mixed-height H22 construction with an exact remaining gap

This branch first produces three exact Coxeter chambers in H22, each with **1,334 actual facets**, then adds two compatible caps that cut their common gap, giving **1,336 facets** per branch. The final chambers still have infinite volume, certified by a new strict spacelike vector. This is a construction within this search, with no claim of literature novelty or a new finite-volume example.

The reproducible [exact checker](../scripts/fresh_h22_incidence_mixed_height_20261013.py) and [certificate](../data/enumerations/fresh_h22_incidence_mixed_height_20261013.json) passed a fresh deterministic replay. A separate [source/proof review and replay of both construction stages](fresh_h22_incidence_mixed_height_20261013_review.md) also passed. The earlier [numerical scout](../scripts/fresh_h22_incidence_mixed_height_20261013_scout.py) and its [output](../data/enumerations/fresh_h22_incidence_mixed_height_20261013_scout.json) are only exploratory evidence; the final conclusion uses the exact checker.

## The changed construction

Write the ambient form as `Q + z²`, where Q is the original H21 Borcherds form, and let h be the old strict center, with h²=−14. The binary incidence formula is

`x=(A v′−Σ u_i P_i)/2`, with `x²=(weight−A²)/2`.

An independent binary linear solve enumerates all 2,048 parity solutions for each parity and then checks every old-wall angle. It reproduces these complete low-type counts:

| Type | Binary profiles | Full orbit under the archived generators | Projection norm |
|---|---:|---:|---:|
| A=3, weight=9 | 280 | 280 | 0 |
| A=4, weight=16 | 21 | 42 | 0 |
| A=3, weight=7 | 360 | 720 together with the next row | −1 |
| A=4, weight=14 | 360 | same 720 orbit | −1 |

The missing half of the 42-element null orbit lies outside the chosen binary coordinate chart. Both null orbits together give 322 distinct integral projections δ. The timelike orbit gives 720 projections x.

The unscaled height-√2 null caps and height-√3 timelike caps have forbidden same-sign cross angles. A changed scale repairs all these relations simultaneously:

`cδ±=(aδ, ±√2)`, where `a=(√6+√3)/2`, and `rx±=(x, ±√3)`.

Every cap has norm two. The exact old pairings satisfy `−(δ,x)∈{2,3,4}`. For equal height signs, the cross pairing is `−ab+√6`: it equals `−√3` when b=2, hence gives π/6, and is strictly dotted for b=3 or4. Opposite signs are dotted. Distinct scaled null caps are dotted; the two signs of one null projection are parallel. All relations with the old 210 walls are Coxeter relations.

The 720 same-sign √3 caps have a compatibility graph consisting exactly of **three disjoint complete graphs K240**. Opposite-height pairs are all dotted. Thus choosing both signs of any one component gives

`210 old + 644 scaled-null + 480 timelike = 1,334 facets`.

The strict common timelike point and nonpositive pairings certify every retained inequality as an actual facet by orthogonal projection. The archive also verifies a representative actual finite vertex for each component: `(x,±1/√3)`, of norm −2/3, with 21 spherical old walls and its new cap. No complete vertex or cusp census is claimed.

## Arithmetic distinction

The pure **unscaled** √2 bank preserves an integral lattice after rescaling the added coordinate: its form is `diag(Q,2)`, all its new roots have integral coordinates and norm two, and `Qw` is even for every old norm-four wall. Rationality alone would not establish this; these divisibility checks do.

The √3 caps behave differently. Every one has an old-wall length-two cyclic product `4(r,w)²/(r²w²)` that is a rational noninteger; the archive records one witness per cap. Thus a finite-volume reflection completion retaining these facets would violate the algebraic-integrality part of [Vinberg's arithmeticity criterion, Theorem 2 (1967)](https://doi.org/10.1070/SM1967v001n03ABEH001992). This is a conditional nonarithmeticity test, not a finite-volume claim. The reflection-integrality convention is also stated directly in [Vinberg's primary paper, §4](https://sfb701.math.uni-bielefeld.de/preprints/sfb14047.pdf).

For the scaled-null branch there is a stronger conditional distinction. A scaled-null cap and an old short wall with `(δ,w)=−1` have cycle `a²=(9+6√2)/4`. Thus the cyclic field contains Q(√2). The normalized Gram of the 42 old short walls is rational and indefinite: their sum is 3h, so its value on the all-ones coefficient vector is −63. Any embedding moving √2 fixes this indefinite principal submatrix. Vinberg's quasi-arithmeticity condition requires positive-semidefinite conjugate Gram matrices under embeddings nontrivial on the cyclic field, which is impossible here. If a finite-volume completion retains these facets, it must therefore be non-quasi-arithmetic. This reasoning assumes the usual cofinite criterion; it is not a claim that the existing infinite-volume chambers are new lattices.

## Exact universal escape

Set

`u=(4,0,0,0,0,0,0,0,0,15,10,20,30,25,20,16,12,0,0,0,0,8)`

and `v=(100u+h,2400)`. The exact checks give

`v²=5,527,986 >0`, `(v,(h,0))=−16,014`, and every old-wall pairing is at most −1.

Across all 322 null projections, `(100u+h,δ)≤−3,206`. Across all 720 timelike projections, `(100u+h,x)≤−4,207`. Consequently both signed new pairings are strictly negative:

`−3206a+2400√2 <0`, and `−4207+2400√3 <0`.

The checker proves these signs by integer squaring. This strict future spacelike vector certifies infinite volume even for the full, mutually incompatible dictionary of 2,294 inequalities. It therefore excludes every subset of that dictionary. It does not exclude changed projection directions, changed scales, or other new facets.

## A constructive step outside that dictionary

The old projection u of the escape satisfies u²=−20. This gives two new caps

`s±=(u,±√22)`.

They have norm two, 21 right old-wall contacts, one π/3 old-wall contact, and 188 dotted old-wall relations. Both are strictly dotted to every scaled-null and √3 cap, and their mutual pairing is −42. Thus both can be added to each of the three chambers without deleting a wall. They strictly cut the preceding escape because `−2160+2400√22>0`.

They also supply actual finite vertices `(u,±20/√22)`, of square −20/11, with the 21 old orthogonal walls. The old21 Gram is positive definite, and every nonincident inequality is strict there. The [two-cap checker](../scripts/fresh_h22_incidence_mixed_height_20261013_gap.py) and [certificate](../data/enumerations/fresh_h22_incidence_mixed_height_20261013_gap.json) verify the formulas, all pair classes, vertices and target cut.

The final result remains infinite-volume: a new rational strict escape has square

`1023438164699189/5000000000000000 >0`.

It survives all 2,296 inequalities of the relaxed union, so all three coherent 1,336-facet branches fail. Its exact vector and field-aware sign bounds are archived. This is the final frozen construction checkpoint; further work would need a genuinely new cap direction or changed retained walls.

A subsequent [independent covering and cusp diagnostic](fresh_h22_incidence_mixed_height_20261013_covering.md) identifies a literal untouched old cusp and an even simpler strict escape of square 8,386. It explains which old cusp direction must acquire new incident caps and why finite collections of timelike-projection caps cannot repair that particular end. It explicitly cites the repository's prior use of this cusp mechanism.

## A concrete local frontier for the next construction

The [small exact facet diagnostic](../scripts/fresh_h22_incidence_mixed_height_20261013_frontier.py) and [output](../data/enumerations/fresh_h22_incidence_mixed_height_20261013_frontier.json) inspect the 21 edges leaving each new cap's saved finite vertex. Every other retained wall is dotted to that cap except its 21 right old neighbors and old wall 14 at π/3. Consequently the intrinsic facet geometry is determined by those 22 old walls; dotted walls cannot cut its hyperbolic interior.

Let `d_i` be the Gram-dual vector to its ith old right neighbor and write `a_i=(d_i,w14)≤0`. The edge is

`v(t)=(u−t d_i,±20/√22)`, with `v(t)²=−20/11+t² d_i²` and `(v(t),w14)=−1−t a_i`.

All 21 cases are checked exactly. Releasing old walls **12, 13 or15** gives a new actual finite endpoint; releasing any of the other **18** walls gives an open edge reaching the light cone before another facet. No case closes at an ideal vertex. These counts apply separately to each of the two new caps in every branch. They are a local frontier, not a global vertex census.

For a particularly simple target, release old wall 1. Its dual d has `d²=8`, `(d,w14)=0` and `(d,u)=0`. The explicit endpoint

`η=(u−√(5/22)d,±20/√22)`

is future null. Exactly the other 20 old right walls and the same-sign new cap vanish there. Their Gram is positive definite, so this is an open ideal edge end, not a completed affine cusp link. All other walls are strict. In particular, the projection of a dotted wall into the cap hyperplane is future timelike and is strict at every future null point of that hyperplane.

This yields a scoped obstruction to repeating only dotted additions. Moving a little farther along the same edge gives a spacelike point with the retained zero inequalities unchanged; after any finite set of dotted additions, those added inequalities are still strict. A small move toward the proposed completion's own strict interior point produces a strict spacelike escape. A successful continuation must therefore introduce a finite-angle or parallel neighbor of the new cap, or replace it.

Replay:

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_incidence_mixed_height_20261013.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_incidence_mixed_height_20261013_gap.py --verify
python scripts/fresh_h22_incidence_mixed_height_20261013_frontier.py --verify
```
