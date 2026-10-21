# L59 fixed-seed checker: source preparation only

20 October 2026. **NOT EXECUTED.** No checker import, generation, replay,
root search, cone census, Lean job or build was run in this lane.

Prepared source:
`scripts/check_compact_h9_golden_discriminant59_edges_20261020.py`,
231 lines / 11,806 bytes, SHA256
`01a569b66b77c38056255a94dce765e4784d809dfe86735ac5bef869755edf3d`.
Its written geometric/profile-completeness premise is the corrected route
note SHA256
`725d8775cebec51d8c18b4b0a0d3be0255e8b1472621adb4d3cefa656adb429e`.

The source-only peer review corrected the arbitrary-arrangement finite
angle list by adding 2-phi and 3-phi. The five nonzero fundamental-wall
values and their native labels are unchanged. The worker now checks
prime-59 divisibility for all seven arbitrary-arrangement values, while
classifying the sixteen seed walls with the fundamental list. The route
also explicitly derives X=x-y in O and the E8-dual integrality of Z.
These changes precede every generation or replay; no earlier certificate
is being replaced or reinterpreted.

The source imports only the standard library and the pinned pure
`field`, `inverse`, `inertia`, `bilinear` definitions of the existing
local-model checker (SHA256 `185d2064...e1bae9`). It calls none of that
checker's model, producer or main functions. The other two input files
are the small calibrated Bugaenko example and its complete transcript.
Every input, the helper and the new worker are hashed before and after.

The fixed normal list has sixteen entries: the eleven calibrated
Bugaenko walls, r, c, u, Sc and Su. Thus it is closed under the proved
chamber involution S; Sv is already the opposite native Bugaenko wall.
The checker reconstructs S from its E8-plus-apex columns, verifies S^2=1,
its full Gram identity, integral entries and the complete native-wall
permutation. It records all 120 distinct pairs, all root norms and full
reflection-matrix integrality, and an exact S-fixed strict center.

In addition, it verifies 41 inherited base-facet vertices and three new
vertices by positive Gram projection and every one of the sixteen wall
inequalities. This is a 44-vertex lower bound, deliberately not an
exhaustive vertex census. The E7 pentagon and both first-simple-edge
arguments remain the written arguments; no list of saved vertices is
used as a substitute for those proofs.

For all five finite endpoint profiles it reconstructs the new root
directly from the full Gram matrix and extracts ALL coefficients of the
quadratic norm identity in X,Z. The six exact interpolation values form
a basis for degree-at-most-two polynomials, rather than a sampling test.
It checks the exact contact polynomial and all numerical inequalities
used by the norm/unit elimination, including the unique A2 root and the
rational-norm-45 nonsquare. The primitive norm inventory, finite-diagram
classification, unit group and the resulting universal exclusion are
explicit written premises, not claims of finite enumeration by this
program.

The natural sixteen-wall seed is still noncompact. On c, the long r,
Sc, Su and all three positive-height Bugaenko walls are ultraparallel.
Only the eight E8 walls and u bound the intrinsic facet; u is orthogonal
to seven of those walls and ultraparallel to native node 7. Each initial
c-vertex edge other than the node-7 direction therefore lacks a finite
endpoint among these walls. Adding Su does not repair this. No full
cone enumeration is proposed merely to rediscover noncompactness.

Root-only prospective argv, each in the usual shared-lock/resource
wrapper (59 s parent, 1 GiB process-tree RSS, 6 GiB start-memory gate):

    ["scripts/check_compact_h9_golden_discriminant59_edges_20261020.py", "--write"]
    ["scripts/check_compact_h9_golden_discriminant59_edges_20261020.py", "--verify"]

The worker has its own 45-second alarm, takes no inner lock, refuses to
overwrite its certificate, and uses canonical JSON comparison on a
separate replay. Expected work is sixteen fixed vectors, 120 pairs,
five quadratic identities and forty-four 9-by-9 positive Gram systems;
there is no unbounded variable, parameter sweep or scientific import.
Execution and final resource measurements are pending root review and
scheduling. Full-host reflectivity and compactness remain open.
