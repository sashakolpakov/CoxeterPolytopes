# Fixed-chart ideal-boundary scores for the replacement hexagon

This [numerical diagnostic](../scripts/fresh_h20_replacement_hexagon_scores_20261017.py) compares the three final hexagon alternatives and their precursor branches in the same Lorentz chart of the target H3 face. The [saved output](fresh_h20_replacement_hexagon_scores_20261017.json) records all sample counts and the SHA256 of every Boolean survival mask. The separately certified [171-wall construction](fresh_h20_replacement_hexagon_20261017.md) uses final label2.

The deterministic rule maps 524,288 unsmoothed two-dimensional Sobol points to the uniform unit sphere using height `2u0-1` and longitude `2*pi*u1`. Every surviving sphere point satisfies the numerical projected halfspaces. Support equalities and future timelike projected normals do not restrict that sphere. No error estimate or exact-area assertion is attached to these counts.

| Branch | Actual walls | Surviving samples | Sample fraction |
|---|---:|---:|---:|
| Original chamber |168|95,833|0.182787|
| Delete158, add168 |168|95,833|0.182787|
| Delete31 and158, add168 and169 |168|81,967|0.156340|
| Same deletions, add168 through171 |170|60,290|0.114994|
| Final hexagon, label2 |171|56,712|0.108170|
| Alternative final label3 |171|57,693|0.110041|
| Alternative final label4 |171|57,979|0.110586|

The final alternatives are separate branches, not simultaneous additions. Label2 leaves the fewest sampled ideal directions and is the simplest final contact. Removing old31 does not remove an actual target-face side: its projection to the target H3 is future timelike. Native ambient compatibility still requires its explicit deletion in the construction.

These are bounded numerical comparisons, not certified areas, hyperbolic volumes, complete cone censuses, or continuous-search exclusions. The exact171 certificate independently proves that a strict spacelike escape remains, so its ambient chamber has infinite volume despite its compact hexagonal two-face.

Status: generation and fresh saved replay passed; source syntax, direct source hashes, and report links were checked. This diagnostic is frozen separately from the exact construction.
