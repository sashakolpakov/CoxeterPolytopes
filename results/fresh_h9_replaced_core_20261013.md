# Replacing the obstructed H9 core

This round constructs two changed 39-wall cores and extends one to **45 actual Coxeter facets**, with no parallel pair and no affine principal subdiagram. It does **not** produce a finite-volume polytope: an exact strict future spacelike escape survives. The additions are tested against a coherent retained configuration throughout.

The earlier [four-wall cusp obstruction](fresh_h9_fixed_core_cusp_obstruction_20261013.md) applies to the old fixed 44/45/46-wall configurations. The chambers here remove the offending old walls before adding anything. They are outside that obstruction's fixed-core scope. This is not a nonexistence result for completions of the ten-wall seed.

## Rebuilding the core

All indices below retain the earlier archive's conventions. Start with the exact 46 walls in [the opposite-cap archive](../data/enumerations/fresh_h9_quartic_completion_20261013_opposite.json). Remove `old31, old51, old197, old496, old1028, new4`, where `new4` is the earlier adaptive cap 19. Remove either `old310` or `old837` as well. This hits the four old parallel pairs and every previously detected positive affine circuit while retaining all ten seed walls.

The [core certificate](../data/enumerations/fresh_h9_replaced_core_20261013_cores.json) checks every retained pair exactly. Both resulting 39-wall configurations have a common strict timelike interior. Their complete connected affine censuses test 426 and 508 principal subsets respectively and find no affine diagram.

The census extends every connected positive-definite principal subset by every finite-paired vertex that joins it. An irreducible affine Coxeter Gram matrix has all proper principal subsets positive definite, so every connected affine diagram is reached. Indefinite subsets and already-affine subsets need no further extension. Disconnected parabolic diagrams would contain a connected affine component. Parallel pairs are checked separately and are absent. This is a complete affine census of these specific chambers.

Absence of an affine diagram does not prove compactness or finite volume. Both starting cores have exact strict spacelike escapes.

## Coupled additions

Let `K=Q(phi)`, with `phi²=phi+1`, and let `G` be the ten-wall seed Gram matrix. Every new normal is represented by its seed pairings `p=-w`, coefficients `c=G⁻¹p`, and norm equation `wᵀG⁻¹w=2`. The retained polar wall `old348` has norm `(6−2phi)/5`.

The two-distance construction fixes all but two entries of `w` to zero, one, or phi, imposes a finite angle with the polar wall, and solves the resulting quadratic norm equation exactly. Its fields are represented by explicit positive-radical towers, including degenerate formal radicals when the actual field is smaller.

The first two accepted additions produce a coherent 41-wall chamber. A further two-distance scan did not cut its leading saved gap. A compatible combination proposed from that scan instead created two affine diagrams and was rejected before any volume claim.

The next construction changes the coupling: it uses the **actual first new cap as a quartic anchor**, with field `K(sqrt(132−60phi))`. Three seed distances are free. A finite angle with that anchor and a second finite angle with the polar wall give two linear equations over the quartic field; the norm condition is then quadratic over it. This yields 27 exactly accepted, leading-gap-cutting proposals. Two coexist with the retained 41 walls, giving an affine-free 43-wall chamber.

The final two-distance search is constrained by all 43 retained walls and targets their actual escapes. Selection is performed separately for each core variant. Every added prefix receives an exact affine census. In particular, one otherwise pair-compatible prefix is rejected because it contains both an affine triangle and a six-node affine diagram. Two accepted final configurations have 45 facets.

The preferred final configuration is case 0 in [the final selection archive](../data/enumerations/fresh_h9_replaced_core_20261013_round4_selection.json). It uses the variant that removes `old837`, and these six additions:

| Source cap | Free seed distances | Fixed nonzero seed magnitudes | Polar label | Additional anchor |
|---|---|---|---|---|
| `polar_caps[0]` | 1, 8 | `w3=1, w9=phi` | 10 | — |
| `round2_caps[0]` | 3, 9 | `w2=w5=phi` | 10 | — |
| `quartic_caps[0]` | 3, 1, 8 | `w5=w9=phi` | 10 | first cap, label 5 |
| `quartic_caps[26]` | 3, 1, 9 | `w5=1, w8=phi` | 10 | first cap, label 5 |
| `round4_caps[2]` | 2, 9 | `w5=w6=phi` | 2 | — |
| `round4_caps[5]` | 0, 9 | `w4=w6=phi` | 4 | — |

All unlisted fixed seed magnitudes are zero. The source prefixes are `data/enumerations/fresh_h9_replaced_core_20261013_`. The archives contain the exact affine lines, norm polynomials, discriminants, selected roots, all seed coefficients, and pair certificates. The other final configuration replaces the last cap by `round4_caps[10]`, with polar label 6.

These additions are distinct from one another and every retained wall: an equal outward mirror would give a positive pairing and fail the exact pair gate. They also cannot reproduce any of the seven deleted mirrors: those have a parallel pair or one of the saved affine circuits with the still-retained seed/core walls. No comparison with every older repository cap bank, and no claim of global mirror novelty, is made.

## Geometry and exact rejection

The selected progression is:

| Facets | Numerical maximum normalized vertex norm | Exact result |
|---:|---:|---|
| 39 | 0.12408086 | strict spacelike escape |
| 40 | 0.08080948 | strict spacelike escape |
| 41 | 0.07357411 | strict spacelike escape |
| 43 | 0.06636894 | strict spacelike escape |
| 45 | 0.06339322 | strict spacelike escape |

The numerical column is only a search heuristic in the fixed positive dual section `sum(y)=1`. It is not a compactness certificate or an exact optimization claim. Each conclusion in the last column is instead justified by an independently replayable strict escape.

For the final chamber, the first saved escape has rational dual coordinates `y`, all 45 inequalities strict, and exact norm

`(684182955592023711 − 325126963008839837 phi) / 2500000000000000000 > 0`.

Its pairing with the saved future timelike center is

`−1598251911/1000000000 − (2077989/8000000) phi < 0`.

The underlying unperturbed numerical vertex has support `[1,3,6,8,9]`; the certified rational perturbation has all ten dual coordinates strictly positive. Thus this is an actual open end of the changed 45-wall chamber, checked against all six additions. There are twelve saved strict escapes for each final configuration.

Every reported wall is a true facet. For a strict timelike interior point `x` and a retained spacelike outward normal `r`, its orthogonal projection `x_r=x−(x,r)r/(r,r)` remains future timelike. Because every other normal has nonpositive pairing with `r`, every other inequality is still strict at `x_r`. This supplies a relative-interior point on each facet. The ten independent seed normals also give the full required Gram rank and signature.

The final 45-wall affine census tests 651 connected principal subsets and finds none. Its remaining infinite volume is established by the strict spacelike escape, not by an inherited cusp. Neither this escape nor the bounded cap searches excludes other completions.

## Replay and scope

[The producer-free checker](../scripts/fresh_h9_replaced_core_20261013_verify.py) is independent of this round's numerical proposal, selection, and volume modules. It reconstructs each selected cap's anchor equations and norm quadratic, verifies all exact pair relations, checks the affine-census frontier for completeness, and replays every strict escape. The four accepted-stage archives give:

| Archive suffix | Configurations | Exact pair checks | Affine inertia checks | Strict escapes |
|---|---:|---:|---:|---:|
| `selection` | 6 | 4,680 | 2,958 | 72 |
| `round2_selection` | 4 | 3,280 | 2,060 | 48 |
| `quartic_selection` | 6 | 5,334 | 3,268 | 72 |
| `round4_selection` | 2 | 1,980 | 1,302 | 24 |

The separate `round3_selection` record is an exact affine rejection, not a volume candidate. The final check rejects a forged cap coefficient and a timelike replacement for the claimed spacelike escape.

The [second-agent geometry audit](fresh_h9_replaced_core_independent_geometry_20261013.md) uses a separate Golden-based recursive radical implementation and a frozen input bundle. It checks the seed signature/inverse, all 45 root norms and seed Gram pairings, the strict center, and all 45 inequalities of the first strict spacelike escape. It relies on the separate pair-angle and affine-census replay for those hypotheses; it does not claim to duplicate them.

```sh
PYTHONPATH=. python scripts/fresh_h9_replaced_core_20261013_cores.py --verify
PYTHONPATH=. python scripts/fresh_h9_replaced_core_20261013_verify.py --source data/enumerations/fresh_h9_replaced_core_20261013_round4_selection.json --verify --self-test
PYTHONPATH=src python scripts/fresh_h9_replaced_core_independent_geometry_20261013.py --verify
```

The proposal families are bounded, numerically screened searches. Their accepted roots have exact producer certificates; the replay table concerns selected roots and chambers. Numerical filtering and rounding-based proposal deduplication are not completeness proofs. Finite angle recognition is exact for each accepted label. The saved label pool is an acceptance alphabet here: composita of different cap fields can be larger than those in the preceding round, so an unrecognized pair is not asserted to be impossible at every Coxeter label. There is no exhaustive cap-family, subset-selection, or whole-seed nonexistence claim.
