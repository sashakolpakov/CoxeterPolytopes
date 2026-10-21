# A targeted finite-vertex walk from the 73-wall chamber

Four coherent wall replacements cut the original leading gap and two subsequent gaps while retaining all ten seed facets and polar wall `old348`. The final chamber still has **73 genuine Coxeter facets**. Four new actual finite vertices survive simultaneously, and the full exact affine-subdiagram check passes. An explicit strict spacelike escape still proves **infinite volume**; no compact or finite-volume example is claimed.

The change from the [preceding construction](fresh_h9_gap_completion_20261013.md) is that each new cap may meet its eight-wall spherical block at non-right finite angles. The new cap and those eight walls form a positive-definite rank-nine diagram. A finite angle to another retained wall then fixes the remaining parameter by a quadratic norm equation.

## Construction

The starting exact escape has seed-coordinate support `[2,3,6,8,9]` before strict perturbation. The numerical proposal stage considers 21 rank-eight spherical blocks made from retained square-two normals over `Q(phi)`. For each block it prescribes up to three nonzero finite couplings, keeps only positive-definite rank-nine extensions, and imposes a finite angle to an outside anchor. Proposals must cut the particular saved escape, retain every seed wall, and require at most two optional-wall removals. These numerical filters are not an exhaustive nonexistence theorem.

The exact construction solves nine linear pairing equations in the anchor's quartic field, leaving `p=u+t v`. The norm equation is

`A t²+B t+C=0`,

where `A=vᵀIv`, `B=2uᵀIv`, `C=uᵀIu−2`, and `I` is the inverse seed Gram. The chosen real quadratic branch is checked exactly. Every accepted cap has norm two, the common strict timelike interior, the intended seed and anchor angles, and a strictly negative pairing with the target in the archive's dual-coordinate inequality convention.

The selected progression is:

| Stage | Optional walls replaced | New walls | Facets | Largest numerical normalized vertex norm |
|---|---|---|---:|---:|
| Starting chamber | — | — | 73 | 0.0567845 |
| First edge | `vertex18` | `walk0` | 73 | 0.0548449 |
| Next edges | `old456, coupled6` | `bridge0, bridge2` | 73 | 0.0501912 |
| Changed four-coordinate gap | `vertex22` | `step0` | 73 | 0.0493593 |

The finite-vertex totals from the numerical projective census are 43, 43, 43, and 46 respectively. These totals are numerical; the four selected new vertices have independent exact point certificates.

The first cap uses `old717`, with field `Q(phi,sqrt(31556/605+46628phi/605))`, as a `pi/5` anchor. It meets `seed3` at `pi/3` and is orthogonal to the other seven walls of its spherical block. Its new finite vertex has diagram `A2+D5+I2(5)`.

The next caps use the quartic `required0` anchor and the existing `vertex12` cap. They cut the next exact escape while retaining `walk0`. The third step targets the newly exposed support `[1,6,8,9]`. Its cap is orthogonal to `old469`, whose field is `Q(phi,sqrt(16phi))`, and meets `seed6` at `pi/3`. The new vertex has diagram `A1+D5+H3`.

All four caps remain together in the final chamber. Each saved vertex point is future timelike, lies on exactly nine selected walls, is strict on every other inequality, and has an exactly positive-definite physical incident Gram matrix. This verifies actual vertices rather than merely spherical abstract subdiagrams.

## Exact gates and rejected alternatives

The first family has eleven certified candidate caps. They all cut the original escape and require `vertex18` removed, but their accepted compatibility graph has no two-cap edge. Seven singleton replacements pass the affine gate and have exact escapes; four have positive affine subdiagrams and are rejected for compactness.

The second stage certifies three targeted alternatives. Five tested coherent selections pass the pair and affine checks, and each has a strict escape. The third stage certifies three further alternatives. Four of its seven selections are rejected by the affine circuit

`{seed1,seed2,seed3,step1}`,

an affine `A3` four-cycle. Three selections remain affine-free but still have exact escapes. Removing additional walls can reopen a larger gap: one third-stage alternative raises the numerical maximum norm to about 0.0562. The final checkpoint therefore uses only `step0` from that stage.

The [producer-free walk replay](../scripts/fresh_h9_vertex_walk_20261013_verify.py) passes:

- seventeen exact cap constructions and their target cuts;
- twenty-three coherent selections, including eight affine rejections;
- 59,870 selected-chamber pair checks, with 3,783 distinct pairs recomputed;
- 18,342 exact affine-frontier checks;
- forty-seven actual finite-vertex checks across individual and selected configurations;
- 180 strict future spacelike escapes.

The affine census is complete for each selected chamber: it extends every connected positive-definite subset, and every proper principal subset of an irreducible affine diagram is positive definite. The finite-angle pool is an acceptance alphabet for the general composita; no exclusion of every possible label outside that pool is claimed.

All facets are genuine. Projecting the common strict timelike interior onto any one wall remains timelike and preserves every other strict inequality because the selected pairings are nonpositive.

## Final checkpoint and remaining gap

The [final exact checkpoint](../data/enumerations/fresh_h9_vertex_walk_20261013_final.json) retains all four new caps and checks all four new finite vertices against all 73 walls. Its first strict escape has norm

`78091742702949107/5000000000000000000
+ (207283610259155447/10000000000000000000) phi > 0`.

The unperturbed proposal now has support `[0,3,7,9]`. Both coefficients of the exact norm are positive. The saved point is strict on every wall and has negative pairing with the future center, so this is a direct infinite-volume certificate. It does not exclude further coherent replacements or additions.

The [final checkpoint replay](../scripts/fresh_h9_vertex_walk_20261013_final_verify.py) additionally verifies the simultaneous four-vertex certificates and anchors them to the full walk audit. The [exact cap archives](../data/enumerations/fresh_h9_vertex_walk_20261013_caps.json), [second-stage archive](../data/enumerations/fresh_h9_vertex_walk_20261013_second_caps.json), and [third-stage archive](../data/enumerations/fresh_h9_vertex_walk_20261013_third_caps.json) save the full line equations, quadratic branches, Gram data, retained pairings, and target cuts.

A [separate arithmetic implementation](fresh_h9_vertex_walk_independent_geometry_20261013.md) also passes the final 73-normal Gram/norm and strict-center checks, the first strict spacelike escape, and the `step0` actual finite vertex with every wall inequality and its physical rank-nine positive Gram. Its narrower geometry scope leaves the complete angle and affine census to the producer-free walk replay above.

```sh
PYTHONPATH=src:. python scripts/fresh_h9_vertex_walk_20261013_verify.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_vertex_walk_20261013_final_verify.py --verify
```

For reuse, `chamber(3)` in [the archive-only loader](../scripts/fresh_h9_vertex_walk_20261013_base.py) returns the final Gram, inverse, strict-center covector, all 73 exact normal profiles, and its geometry case. The underlying case is index 4 of the third geometry archive. The bounded numerical proposals and the finite set of tested replacements do not exhaust these continuous cap families. No comparison with all earlier mirror banks or claim of global novelty is made.
