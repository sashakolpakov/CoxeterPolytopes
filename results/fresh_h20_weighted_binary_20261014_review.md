# Independent compressed review of the binary cap construction

**PASS.** The frozen [producer](../scripts/fresh_h20_weighted_binary_20261014.py) and [certificate](fresh_h20_weighted_binary_20261014.json) give a coherent 1,048,615-facet chamber, strictly cut every half-vertex probe, and preserve 524,288 actual finite endpoints. The saved strict spacelike escape proves that the chamber still has infinite volume.

The separate [review checker](../scripts/fresh_h20_weighted_binary_20261014_review.py) and [review record](fresh_h20_weighted_binary_20261014_review.json) use only standard-library `Fraction` arithmetic in `Q(sqrt(5))`, with no producer or project imports. The primary arithmetic uses coefficients in `Q(phi)`. Independent generation and saved replay pass. This is a compressed all-pair proof review and independent reconstruction, not a million-wall Gram-matrix calculation or a finite-vertex census. The separate local frontier work is outside this review.

## Why the compressed pairs cover the whole orbit

The full-coordinate triples are `A={0,1,2}` and `B={0,3,4}`. Their intersection has size one, symmetric difference size four, and complement among the ten unit coordinates size five. The type depends on unit-plus-golden sign parity; the heavy sign does not enter it.

For equal types the squared center distance is

`4*h_full + h_other + 2*h_heavy + (1+phi)*h_golden`.

The unit-plus-golden difference parity is even. For different types the squared distance is

`1 + 4*h_intersection + 2*h_symmetric_difference + h_neither + 2*h_heavy + (1+phi)*h_golden`,

with odd difference parity. The constant one is the sum of four squared half-unit offsets. A sign reversal in such a coordinate changes its squared contribution from `1/4` to `9/4`, giving coefficient two. These formulas apply to either reference type; they use only the indicated coordinate-class cardinalities.

The review independently reconstructs a literal center pair for every one of the **828** compressed rows, checks its squared distance and Coxeter relation, and recomputes the binomial multiplicities. Each type half has exactly `2^18` sign choices. The reference finite contacts are 27 right, four at `pi/3`, and 64 at `pi/5`, giving exactly 95 finite binary neighbors. Their explicit masks are also reconstructed and checked, with no duplicates. All remaining pairs are parallel or dotted.

The explicit seed relations are independently recalculated: three right, seven at `pi/3`, one at `pi/4`, eight at `pi/5`, and nineteen dotted. The central contact is dotted. Every matched half-probe has power `-1/4`, so its cap cuts it strictly.

## Endpoint layer, facets and finite endpoints

The radius-`sqrt(2)` endpoint orbit is compatible in full. Its normalized pair cosine is `H-1`, where `H` is the weighted Hamming difference; the smallest positive case is a right contact, and every other case is parallel or dotted. The independent checker recomputes these weighted counts.

Every binary–endpoint squared center distance is at least `17/4+2phi>6>(1+sqrt(2))²`. Hence all such pairs are dotted. This minimum follows coordinate by coordinate and holds for both binary types and every sign pattern.

All sphere normals are `(c,(|c|²-R²)/2,1)`. The timelike point `(0,100,1)` is strict for the seed, central, binary and endpoint walls. Since all pairings are nonpositive, projection onto any one wall yields a timelike point strict on every other wall. Distinct centers and radii give distinct walls. This proves that all counted walls are actual facets.

At a box vertex `z`, the retained endpoint point is

`P_z=(z,(N+2)/2,1)`, with `P_z²=-2`.

It lies on nineteen mutually orthogonal seed walls and its own endpoint cap, which form a positive-definite twenty-wall Gram matrix. Its pairing with a binary cap of center `c` is `-(|z-c|²+1)/2<0`; with another endpoint cap it is `-|z-w|²/2<0`. The central and remaining seed inequalities are also strict. Thus all 524,288 endpoint vertices are actual. No assertion about completeness of the finite-vertex set is made.

## Uniform strict escape

The saved spatial point has ten unit coordinates zero and all heavy/golden coordinates `99/100`. At this point, every binary cap has unit-coordinate squared-distance contribution `19/4`, independently of type and signs. Matching the other signs minimizes its remaining contribution, giving the recorded uniform power. The endpoint minimum is obtained in the same way.

The review reconstructs the norm-one Lorentz lift, all three minimum powers, the minimum normalized seed gap `1/100`, and the future orientation. Every minimum power exceeds one. The escape therefore survives the entire orbit and endpoint layer, despite strict coverage of all half-vertex probes.

```sh
python scripts/fresh_h20_weighted_binary_20261014_review.py --verify
```
