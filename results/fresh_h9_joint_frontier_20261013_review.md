# Review of the exact 27-wall frontier census

**PASS.** The [frontier source](../scripts/fresh_h9_joint_frontier_20261013.py)
and [frozen certificate](fresh_h9_joint_frontier_20261013.json) correctly
certify 65 finite vertices and all 441 edges incident to them: 144 have
two finite ends and 297 continue to the ideal boundary with no existing
wall terminating them. The finite-edge graph is connected. Fifteen of
these vertices are additional to the 10 inherited vertices and the 40
endpoints in the preceding construction.

This review checks the mathematical coverage and source, and reruns the
primary saved checker. It is not a second arithmetic implementation or a
duplicate independent enumeration. The separate producer-free replay
described in the [main report](fresh_h9_joint_frontier_20261013.md) has a
different scope from this review.

All normals are normalized to square two without changing their oriented
hyperplanes. The checker recomputes every pair: finite contacts have Gram
entries `0,-1,-phi` over `Q(sqrt(5))`, while every dotted comparison is
strict. Thus any positive-definite nine-normal vertex support must be a
nine-clique of the finite-pair graph. All 351 such cliques are enumerated;
286 fail positive definiteness. For each remaining support, exact center
projection produces a future timelike feasible point whose zero set is
exactly that support. A finite Coxeter vertex is simple, and these tests
therefore give the full finite-vertex census.

At each finite vertex, deleting an incident wall gives a rank-eight
positive-definite edge support. Every existing wall is tested as a
possible ninth wall. A dotted pair rules out positive semidefiniteness;
otherwise the exact Schur complement tests the complete finite/affine
possibilities. Positive complements reconnect to the complete finite
vertex list. Zero complements produce an explicitly checked nonzero
future null residual and test all wall inequalities. This excludes
unexamined existing-wall endpoints on the recorded edges.

Every one of the 297 remaining edges has an exact feasible spacelike
continuation beyond its null endpoint, strict on every wall outside the
eight-wall support. The proposed rational step is numerical, but its
norm, orientation and all inequalities are checked exactly. A small
perturbation toward the common strict center also yields a strict
spacelike escape. The count is consistent with all finite-vertex
incidences: `65*9 = 2*144 + 297`.

The spherical endpoint-profile argument is complete for the stated local
problem. A finite Coxeter diagram is a forest. A new node consequently
meets at most one node of each old connected component. Degree four is
impossible, and degree three requires all new labels to be three, by the
positive-definite star Schur complement. Any label at least six must
belong to an isolated dihedral component; all these `I2(m)`, `m>=3`,
possibilities are retained as tails at isolated old A1 nodes. All other
labels 3,4,5 and positive Schur complements are explicitly enumerated.
The resulting 10,057 finite profiles and 418 dihedral tails are necessary
local endpoint data, not assertions that corresponding global new walls
exist.

The source explicitly limits completeness to finite vertices and their
incident edges. It does not census edges with no finite endpoint, all
exterior projective rays, affine new-wall profiles, or arbitrary new
mirrors. No finite-volume or compactness claim is inferred.

Fresh saved replay passed:

```sh
python3 scripts/fresh_h9_joint_frontier_20261013.py --verify
```

Frozen SHA-256 bindings:

| File | SHA-256 |
| --- | --- |
| `scripts/fresh_h9_joint_frontier_20261013.py` | `f8ff35d26b29a42107df6d9a53662d74bebdca12674809e2a91f68915848ac8b` |
| `results/fresh_h9_joint_frontier_20261013.json` | `9bc257ed0393ffd125d088c6ca81da210b9110b809bb88e35173bc41d5de1d23` |

No unresolved mathematical issue was found in this scope.
