# A proper-cusp completion works locally, but increases unresolved boundary

The exact conditional affine calculation and the numerical realization are separate. The former proves a completion mechanism for any exact coherent base satisfying the stated hypotheses. The latter evaluates that mechanism on the unpromoted175-wall numerical branch and recounts its complete target-H3 boundary. Generation and a fresh saved replay pass; no exact promotion of this base is claimed.

The base uses triangular replacement20, continuation14, and continuation112 from the [replacement](fresh_h20_replacement150_triangles_20261017.md) and [continuation](fresh_h20_replacement150_continuation_20261017.md) packages. It retains the original walls except31,150,158. Its worst exterior plane intersection is160/174/176, with squared radius1.0354125853001426. That exterior intersection is **not** asserted to be a cusp.

Instead take the forward ideal endpoint ν of the genuine open edge174/176. The full incident support is

`W = T + {174,176} = A10 + A4 + A2 + A1 + A2`,

of spherical rank19. Numerically, ν has positive time coordinate, zero norm and exactly these19 incident retained walls. Every other retained wall has negative product withν; the largest is about−0.338814. No noncrystallographic component obstructs an affine completion.

For each of the five components C, define

`r_C = −Σ_{j∈C} n_j + h_C ν`, with `h_C > 0`.

All components are normalized A-types. Their sums have square2. Consequently the new roots have square2, distinct new roots are orthogonal, and each completed component is affine A of the same rank. For A1 the new pair is parallel. The exact positive relation is

`Σ_{j∈C} n_j + r_C = h_C ν`.

The saved integer matrices verify Euclidean Gram rank19 and normal-span rank20, equal to `ν⊥`. The five affine components therefore give a bounded product of Euclidean alcoves in the19-dimensional cusp link. This is a proper cusp, not a deficient affine subdiagram.

For any old nonincident wall N, write `p_C = −Σ_C n_j`. Since `ν·N < 0`, the sufficient dotted-contact bound is

`h_C > (p_C·N + √(2N²))/(-ν·N)`.

Together with the analogous strict-center inequality, finitely many such bounds give a finite sufficient height. This proves existence conditionally on an exact base with the indicated feasible ray. The selected numerical heights are8,7,7,7,24; the corresponding largest dotted-height bounds are saved with their blocking walls. A common future point on the T-face is strict for all non-T old and new roots. A sufficiently small transverse perturbation gives an ambient strict point. These physical margins and the displayed heights have not been independently certified with outward intervals.

Four new roots project to positive multiples of the future null rayν in `T⊥`; their inequalities cannot cut the interior of that H3. They are therefore omitted from its side ledger, with this analytic justification. The fifth, completing174/176 to affine A2, has positive projected square and is an actual target side.

All counts below concern the **target H3 only**, not an ambient vertex or cusp census. An open end is an ideal edge endpoint not certified as the completed proper cusp. Thus edges ending atν are counted as closed at that end.

| Target branch | Finite vertices | Proper ideal vertices | One-open edges | Two-open edges | Open ends | Exterior skeleton components | Max radius² |
|---|---:|---:|---:|---:|---:|---:|---:|
| Before completion |25|0|9|1|11|1|1.035413|
| Five-root completion, A2 height24 |25|1|10|1|12|1|1.031289|
| A2 height giving160 label2,3,or4 |25|1|10|2|14|2|1.031289|
| A2 height giving160 label5 through12 |25|1|10|1|12|1|1.031289|

The exterior skeleton component count is a numerical polyhedral graph diagnostic, not a topological classification of the ideal boundary. At height24, the raw positive-side polytope has25 finite vertices,10 exterior vertices and one point on the unit sphere; that last point is identified withν by the conditional affine relations.

There is a useful one-parameter description of the A2 family. With `p = −n174−n176`, prescribing the native160 cosine c gives

`h(c) = (p·n160 + 2c)/(-ν·n160)`.

The right-contact height is approximately18.085050 and the parallel height22.658132. At the right-contact endpoint all other unprescribed native products are already strictly dotted, with minimum numerical margin greater than0.7058; their slopes in h are negative. Thus the exact linear implication covers all c in[0,1] whenever those endpoint margins are rigorously verified. The stored bounded realization checks all finite labels2 through12; all eleven pass. No finite-label cutoff is used in the conditional linear implication. Heights above the parallel threshold are also permitted by the same old-wall dotted inequalities.

The outcome is unfavorable for whole-face closure: a proper cusp is created, but the number of unresolved ends increases. No exact promotion or larger cusp batch is selected on this evidence. The [source](../scripts/fresh_h20_replacement150_parabolic_cusp_20261017.py) and [data](fresh_h20_replacement150_parabolic_cusp_20261017.json) preserve the exact affine matrices, positive relations, native height bounds, center checks, every bounded height event and full numerical boundary ledgers.
