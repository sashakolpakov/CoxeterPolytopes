# Review of the H22 replacement orbit and its global limitations

**PASS, scoped source review and fresh exact replay.** The complete stabilizer
orbit is a coherent infinite-volume chamber, and the escape and endpoint
obstructions have the stated restricted scopes. This review does not certify
a new finite-volume example or a complete census of the changed chamber.

The review read the complete orbit checker, the complete 22-edge checker and
the separate endpoint Schur-bound checker. Fresh saved-output replays of all
three passed. The accompanying [manifest](fresh_h22_replacement_global_20261013_review.json)
binds the reviewed source and certificate versions.

The [orbit calculation](fresh_h22_replacement_orbit_20261013.md) uses the full
5,760-element stabilizer within the separately certified old group. It checks
preservation of both short-wall parts and all 301 unsigned parent recipes;
this justifies reducing orbit pairs to a representative. Integer transport,
strict pair signs, the common parameter tail and actual finite vertices are
checked directly. Positive coefficients and an invertible three-functional
matrix prove the escape obstruction on each entire residual plane, with no
height or angle cutoff. The final inherited-cusp diagnostic verifies all 24
incident walls, positive affine relations of total quotient rank 20 and an
explicit vertical direction in a 21-dimensional link. The new construction
uses integer array arithmetic with pre-evaluation size bounds. Its old group
and source reconstruction remain explicit dependencies; this review is not a
second independent implementation of those dependencies.

The [edge calculation](fresh_h22_replacement_edges_20261013.md) enumerates all
22 directions at the specified finite vertex. Exact dual vectors give a null
point and a feasible spacelike continuation for each direction. These are
open ends, not bounded cusps. Every direction retains one new cap, to which
all nonincident old walls are strictly dotted; this also explains why no
existing wall can terminate any of these edges.

The endpoint checker independently replays both literal rational decompositions
of the escape from the [dual certificate](fresh_h22_replacement_escape_dual_20261013.md).
They force any acute escape-cutting wall to pair negatively with both old
walls 7 and 168. Each 21-wall edge Gram has entrywise nonnegative inverse.
Normalizing an added wall to square two, the smallest allowed nonright
Coxeter pairings already give projected square at least three. Its extended
Gram therefore cannot be positive semidefinite. This excludes an escape-cutting
finite or ideal endpoint on every one of these 22 edges, with no finite-label
bound. Parallel relations in an affine support are included: their magnitudes
are larger than the label-three lower bound. Strictly dotted relations are
excluded by positive semidefiniteness.

The only report correction was to include parallel, as well as finite-angle,
relations in that affine-support explanation. No certificate change was needed.
Intermediate additions that leave the escape feasible, subsequent vertices,
normals outside the residual planes and further old-wall replacements all
remain possible. None of the checks proves that this retained core has no
finite-volume completion.

```sh
python3 scripts/fresh_h22_replacement_orbit_20261013.py --verify
python3 scripts/fresh_h22_replacement_edges_20261013.py --verify
python3 scripts/fresh_h22_replacement_edges_20261013_endpoint_gate.py --verify
```
