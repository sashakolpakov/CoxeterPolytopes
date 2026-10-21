# Balanced cusp numerical endpoint and central-radius scouts

No finite-volume candidate was found. The banks below propose individual
caps; pairwise compatibility and sampled point coverage do not establish
finite volume.

The [endpoint scout](../scripts/fresh_h20_balanced_20261014_endpoint_scout.py)
uses the complete 41-profile census and all 392 transports. It anchors radii
at omitted-wall or central contacts with finite labels 2 through 16, 20, 24,
30, and 60, and at parallel contacts. Each cap passes its seed and central
pair gates numerically. Mutual cap compatibility is not imposed.

| Scale choice | Caps | Finite endpoint profiles | PSD endpoint profiles | Uncovered probes |
|---|---:|---:|---:|---:|
| `(1,1,1)` | 284,592 | 255,080 | 29,512 | 392 of 1,961 |
| `(1,1,phi)` | 284,200 | 256,648 | 27,552 | 784 of 1,961 |

The probes are the origin and every old vertex scaled by 1/8, 1/4, 1/2,
3/4, and 1. The unit bank misses all quarter-vertex probes. The PSD records
have no bounded cusp-link certificate. Details and array hashes are in the
[unit record](fresh_h20_balanced_20261014_unit_endpoint_scout.json) and
[golden record](fresh_h20_balanced_20261014_golden_endpoint_scout.json).
The [all-radius theorem](fresh_h20_balanced_endpoint_20261014_central_gap.md)
is stronger than this sampled failure.

A separate [central-radius scout](../scripts/fresh_h20_balanced_central_radius_20261014.py)
changes the unit-seed central radius to `R0=1/cos(pi/k)` and tests identical
scaled-vertex caps `w=gamma*z`, `R=(1-gamma)/cos(pi/n)` with a prescribed
central contact `pi/m`. It tests 980,100 profiles: `k,n=3,...,100` or
parallel, and `m=2,...,100` or parallel, across all 41 old vertex-distance
values. There is no surviving finite pair contact. Seven numerical near
hits, all at large inferred labels, are rejected by separate 192-bit Arb
intervals. The [saved record](fresh_h20_balanced_central_radius_20261014.json)
contains those nonzero difference bounds.

The large scout remains numerical: its seven near-hit rejections are
certified, but this is not an exact exhaustive exclusion of the parameter
box. It does not test pairs with different cap profiles. The mixed-profile
connector found by the bridge search is therefore outside this negative
result.
