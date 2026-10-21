# Root review of the two trivalent H9 endpoints

The frozen exact checker passed a separate root replay:

`python3 scripts/fresh_h9_component_bridge_20261013_verify.py --verify`

The source audit checked that the nine finite pairings and norm two follow
from exact projection identities, rather than from intervals containing
zero. All other new pairs have strictly dotted upper bounds. The affine
subset enumeration covers every subset involving a new wall: each new wall
has exactly nine finite neighbors, and their mutual dotted pair excludes
any positive-semidefinite subset containing both new walls. Each of the two
new finite vertices is tested against every nonincident wall. The strict
center and nonpositive normal pairs justify actual facets.

The strict remaining escape is checked against all 104 walls and has positive
norm and the correct time orientation. Thus this is constructive progress,
with infinite volume still proved. The cut of the entire preceding cap-bank
escape has an explicit dependency on that bank's certificate and hashes.

A second reconstruction uses 110-digit arithmetic, a nine-equation coordinate
linear system, its one-dimensional nullspace, and a scalar quadratic. It
does not use the new normals' successive orthogonal-projection recipe.
It agrees with the rigorous normal enclosures to within `1.4e-72`, reproduces
the finite pairings, and checks all 205 new pairs. The two largest old dotted
pairings are approximately `-2.30115541426` and `-2.12244566254`; the new
cross-pairing is approximately `-70.77975700289`. The remaining escape's
smallest inward inequality is greater than `0.0000093109`.

This high-precision reconstruction is supplemental numerical evidence; the
exact producer-free replay supplies the certificate. The prior 102 normals
remain a shared, explicitly replayed dependency.

The [separate script](../scripts/fresh_h9_component_bridge_20261013_root_review.py)
and [saved result](fresh_h9_component_bridge_20261013_root_review.json) record
the reconstruction and source hashes.
