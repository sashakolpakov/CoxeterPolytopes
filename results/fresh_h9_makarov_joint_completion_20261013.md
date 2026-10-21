# A new finite corner beyond the old edge-cap banks

The frozen 27-facet H9 chamber has been extended to **29 actual Coxeter facets**, with a new actual finite vertex on a previously unused codimension-seven face. Both added walls cut the exact escape that defeats the earlier 25,101-entry cap bank. Nevertheless, the 29-wall chamber has another exact strict spacelike escape and therefore **infinite volume**. This is a local construction, not a compact H9 example.

The primary certificate is [coupled_exact.json](fresh_h9_makarov_joint_completion_20261013_coupled_exact.json), replayed by [coupled_exact.py](../scripts/fresh_h9_makarov_joint_completion_20261013_coupled_exact.py). The [root source/proof review](fresh_h9_makarov_joint_completion_20261013_coupled_exact_review.md) also replayed it. A [separate arithmetic audit](fresh_h9_joint_frontier_20261013_coupled_independent.md) independently checks the new geometry and gives the complete census of **66 finite vertices and 455 edges**: 144 finite edges, 306 edges with one finite endpoint, and five with neither endpoint finite. The construction retains all 27 facets from the [previous exact chamber](fresh_h9_makarov_joint_20261013.md); its source files are unchanged.

## The exact new corner

Let `K = Q(sqrt(5))`. Both new mirrors have projective normals in K. They are orthogonal to old walls

    S = [0, 2, 3, 5, 7, 8, 9].

The first projective vector `v` is uniquely determined by these seven right equations and the free coordinates

    (v_1, v_4, v_9) = (207/4, 10, -129/4).

Its exact square is

    q = 24967/12 - 29745*sqrt(5)/32 > 0,

and its norm-two normal is `sqrt(2/q)*v`. The second projective K-vector is orthogonal to `S` and to `v`; its exact rational free coordinates and normalization are stored in the certificate. The choice of those rational coordinates is a proposal, not a minimality claim.

All **55 new pairings** are checked exactly: 15 are right, and 40 are strictly dotted. In particular, no parallel pair is introduced. The common old strict timelike center remains strict. Projection of that center to each wall gives a timelike point strict on all other inequalities, so both new walls and all retained walls are actual facets.

The nine walls `S + [27,28]` have a positive-definite Gram matrix: the old PD7 block is orthogonal to both new positive norm vectors, which are mutually right. The certificate projects the center to their common line and checks a future timelike point with zeros **exactly** at those nine walls and strict negative pairings with every other wall. This is an actual finite vertex. The new caps have no finite contact with any old wall outside S, so this new corner initially lies in a different component of the finite-edge graph from the old corner graph; it is a new place from which to seek bridges.

Both new caps cut the [old-bank exact escape](fresh_h9_makarov_joint_completion_20261013_old_bank_gate.json). The following different integer point is strict for all 29 walls:

    (15272584, 19332841, 13954294, 8575997, 3198152,
     706462, -27462556, -33326515, -22217679, -11108841).

Its exact square is `2712511767010814 - 1210493722603484*sqrt(5) > 0`; its center pairing is negative. The certificate checks every wall sign. Thus it gives an open ideal end and proves infinite volume.

## What the preceding searches establish

The first row below is an exact exclusion for a fixed dictionary and fixed retained chamber. The other rows are numerical scouts and do not prove that the displayed families are exhausted.

| Search | Scope and outcome |
| --- | --- |
| Earlier cap dictionary | [Exact old-bank gate](fresh_h9_makarov_joint_completion_20261013_old_bank_gate.md): all 25,101 archived entries either miss one common strict escape or conflict with a retained wall. No compatible selection from that dictionary completes the retained 27 walls. |
| New edges at finite corners | [Finite-anchor scout](fresh_h9_makarov_joint_completion_20261013_profiles_all.json): 99,207 rounded proposals from the 297 open directions. Their whole-bank relaxation, together with the older compatible bank, still has a numerical escape. |
| Free nonedge contacts | [Continuous scout](fresh_h9_makarov_joint_completion_20261013_continuous.json): 13,819 sampled-profile instances, 41,567 rational-boost proposals. All nonedge contacts are strictly dotted. None cuts the later relaxed-bank gap. |
| Five edges without finite endpoints | The [exact complete edge census](fresh_h9_joint_frontier_20261013_all_edges.json) supplies five previously missed supports. The [remote-edge scout](fresh_h9_makarov_joint_completion_20261013_remote_edges.json) tested 269 profile instances, 102,220 finite-anchor quadratic roots and 538 continuous intervals, yielding 6,413 rounded proposals. None cuts either the later gap or the exact old-bank escape. |
| PD7 face mechanism | [All-right residual scout](fresh_h9_makarov_joint_completion_20261013_faces.json): 2,234 finite-pair seven-cliques, 1,339 numerically PD supports, and optimization of at most 150 promising residual cones. Five proposals cut the exact old-bank escape. The simplest yielded the two exact caps above. |
| New corner's outgoing profiles | [Nine-edge scout](fresh_h9_makarov_joint_completion_20261013_coupled_frontier.json): 1,094 profile instances and 457,120 quadratic roots yielded 8,523 proposals. Many cut the new 29-wall escape, but their [whole-bank relaxation](fresh_h9_makarov_joint_completion_20261013_coupled_relaxed_volume.json) still escapes. All also miss the later 99,207-bank gap. |
| Changed PD7 support at that later gap | [Bounded all-right scout](fresh_h9_makarov_joint_completion_20261013_faces_latest.json): 449 supports permit a cut under acuteness alone; optimization of the first 150 with every other relation dotted found no cutter. This is not an all-parameter obstruction. |
| Nonright PD7 profiles at that gap | A [small changed-profile pass](fresh_h9_makarov_joint_completion_20261013_nonright_faces.json) on the 16 nearest exact supports tests 137 single-contact profiles, with labels 3/4/5. It finds 24 numerical cutters, subject to the remaining old relations being dotted or parallel. These are proposals only. |

The final [unconditioned union test](fresh_h9_makarov_joint_completion_20261013_nonright_union_volume.json) includes the old compatible entries, the 99,207 new finite-corner entries, remote-edge proposals, both all-right PD7 passes, the 131 feasible nonright proposals, and the new-corner profile bank: 117,228 archived cap entries, not deduplicated, plus the 11 seed walls. After 100 linear-programming objectives, a numerical feasible point still has square about `4.67096e-5 > 0`; it is also strict for both exact new caps. This point has not been upgraded to an exact certificate for that large union. The reproducible scout is [union_volume.py](../scripts/fresh_h9_makarov_joint_completion_20261013_union_volume.py). Thus the small nonright-profile pass cuts the previous gap but does not yet justify compatible selection or another exact cap insertion. The next search can change the contact pattern or retained-wall choice; this numerical result does not exclude an unsampled continuous cap.

Finite spherical profiles do not describe all possible added walls. In particular, an isolated wall near a strict ideal gap can have no finite contact with an existing wall; coupled constructions can first create a new finite corner, as the PD7 construction demonstrates. Conversely, cutting a saved escape, closing old edges, or finding more finite vertices is insufficient for finite volume. Every positive proposal here is followed by a global escape test.

## Reproduction

The exact saved cap/vertex/escape replay requires no new optimization:

    python scripts/fresh_h9_makarov_joint_completion_20261013_coupled_exact.py --verify
    python scripts/fresh_h9_makarov_joint_completion_20261013_old_bank_gate.py --verify

The scripts with suffixes `_profiles.py`, `_continuous.py`, `_remote_edges.py`, `_faces.py` and `_coupled_frontier.py` reproduce numerical scouting only. Their tolerances, sampled label lists and scope are stored in their outputs. The exact construction depends on the hash-bound 27-wall certificate and its earlier checks; it does not repeat that predecessor's entire global geometry census.
