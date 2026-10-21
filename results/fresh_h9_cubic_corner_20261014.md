# A nearby fixed point extends a compact-H9 arithmetic search

**The cubic-field host now has 32 certified walls of its full reflection chamber, but this finite prefix still has infinite volume.** There is no new compact polytope and no proof that the full lattice is nonreflective. The useful advance is a chamber-symmetry certificate at a nearby fixed face, which validates a distant corner without enumerating every mirror out to that corner.

Let `θ=2cos(2π/7)`, with `θ³+θ²−2θ−1=0`. The lattice is `⟨−2θ⟩+E8+⟨1⟩`, over the integers of the real cubic field. Its nonidentity real conjugates are positive definite, so a finite-volume arithmetic chamber here would be compact. This dimension-nine route is unaffected by the noncompact arithmetic dimension restriction. The [previous bounded-search scope audit](fresh_cubic7_h9_independent_scope_audit_20261013.md) supplies the complete norm alphabet, initial chamber and height-100 enumeration premise.

## Fixed-point certificate

The [exact corner search](../scripts/fresh_h9_cubic_corner_20261014.py) checks 31 numerical corner proposals in three saved cubic hosts, and constructs integral, form-preserving maps from norm-preserving diagram isomorphisms. Its [saved output](fresh_h9_cubic_corner_20261014.json) distinguishes proposed prefix corners from corners of the full chamber. Only the middle host supplies a nonidentity map T; it is an involution.

T sends the initial corner to a proposed corner at distance height about 1009.44, beyond the complete prefix. Direct transport from that corner alone would be unjustified. Instead, set `H=u+Tu`, where `u=(1,0,…,0)`. The [extension certificate](fresh_h9_cubic_corner_20261014_extension.json) proves:

- T is an integral lattice isometry with integral inverse, `T²=I`, and preserves the future cone.
- H is timelike, fixed by T, satisfies every old wall inequality, and has radial distance height `2+7θ+3θ²<16<100`.
- T positively permutes precisely the seven incident simple walls at H.

For a point X, the radial height is `sinh²d(u,X)=2θ X₀²/(−X²)−1`. For a root r, the mirror height is `2θ r₀²/r²`. Every omitted wall meeting the segment from u to H, or containing H, would therefore have height below 100. Completeness of the inherited prefix rules it out. Thus H's local chamber sector is the full chamber's sector. T preserves the reflection arrangement and this sector, so it preserves the entire chamber. All images of the original seventeen walls are consequently actual walls of that same chamber.

Their union has 22 walls and cuts the original escape. The [separate Cubic7 arithmetic review](fresh_h9_cubic_corner_20261014_review.md) reconstructs the isometry, fixed point, oriented incident-wall permutation and all 231 pairs, without importing the Fraction-polynomial producer. It verifies the use of the bounded-distance premise but does not repeat the inherited root enumeration.

## Targeted continuation and full prefix census

The 22-wall [exact elliptic-face census](fresh_h9_cubic_corner_20261014_edges.json) exposes proposed corners at distance height about 140.3776. This selects the next complete distance bound 141. The [continuation](../scripts/fresh_h9_cubic_corner_20261014_h141.py) matches every previous control level, enumerates the two new levels, and adds six walls. The [raw prefix](fresh_h9_cubic_corner_20261014_h141.json) has 23 walls; its [union with its T-image](fresh_h9_cubic_corner_20261014_h141_orbit.json) has 32. Every pair passes exact Coxeter conditions, a common strict timelike center remains, and T permutes all 32 walls.

The [complete elliptic rank-eight/rank-nine census](fresh_h9_cubic_corner_20261014_h141_edges.json) checks all eligible subsets of the finite-pair graph, with exact positive-definiteness tests. It finds 24 finite prefix vertices; 16 are certified full-chamber corners by the complete distance-141 neighborhoods of u or Tu. Among 200 elliptic rank-eight edges, 28 have two finite endpoints, 160 have one, and twelve have none. This is a census of the saved prefix, not a completed walk of the full chamber.

The [map search over these corners](fresh_h9_cubic_corner_20261014_h141_symmetries.json) finds only the identity and T among integral maps induced by literal norm-preserving corner-diagram isomorphisms. No infinite-order stabilizer is established. Other lattice symmetries and undiscovered corners remain possible.

Finally, a [strict spacelike escape](fresh_h9_cubic_corner_20261014_escape.json) satisfies all 32 wall inequalities strictly. Its positive square and all inequality signs are checked with both Cubic7 and the separate Fraction-polynomial implementation. Together with the strict timelike center, it proves infinite volume of this prefix. An undiscovered full-chamber mirror may cut it, so this is not a nonreflectivity theorem for the host.

## Replays and next step

The new two-level enumeration was rerun completely and matched the saved certificate. The orbit, strict-escape, full census and symmetry saved replays passed, as did the separate fixed-point audit. A [narrow height-141 review](fresh_h9_cubic_corner_20261014_h141_review.md) independently reconstructs the saved new root packets and all 496 pairs of the 32-wall orbit, and reviews the continuation's control-matching logic. It does not independently enumerate the new control or candidate census. These scopes should not be conflated: exact saved geometry, independent arithmetic, and independent exhaustive enumeration are different checks.

```sh
PYTHONPATH=src python3 scripts/fresh_h9_cubic_corner_20261014_h141.py --verify
python3 scripts/fresh_h9_cubic_corner_20261014_h141_orbit.py --verify
python3 scripts/fresh_h9_cubic_corner_20261014_edges.py --source results/fresh_h9_cubic_corner_20261014_h141_orbit.json --output results/fresh_h9_cubic_corner_20261014_h141_edges.json --verify
python3 scripts/fresh_h9_cubic_corner_20261014_symmetries.py --source results/fresh_h9_cubic_corner_20261014_h141_orbit.json --edges results/fresh_h9_cubic_corner_20261014_h141_edges.json --output results/fresh_h9_cubic_corner_20261014_h141_symmetries.json --verify
python3 scripts/fresh_h9_cubic_corner_20261014_escape.py --verify
python3 scripts/fresh_h9_cubic_corner_20261014_review.py --verify
```

Further work should resolve selected open edges or recognize another full-chamber symmetry. The saved escape and finite-edge census provide concrete targets. Increasing a radius alone, or finding a high prefix corner without validating its full sector, is insufficient to decide reflectivity.
