# A complete local packing check after 815 cross-coface capacities

An exact witness of ratio **12.311548327014249 > 12** survives the expanded 815 cross-coface capacities and the **entire nonnegative local simplex-block packing family** on its positive support. Its quadrilateral masses use 336 fixed actual provenances. This is an obstruction for the specified incidence model, not an ambient dimension theorem.

The [script](../scripts/fresh_sqrt5_orbit_cap_restricted_local_cones_20261017.py) and [exact certificates](fresh_sqrt5_orbit_cap_restricted_local_cones_20261017.json) use fresh changed-mass local computations. Generation and fresh saved replay pass; replay invokes no optimizer. No global LP was run.

## Common provenance and repaired face masses

The starting point is the clean `all_S_and_selected_actual_Q` restriction in the [support certificate](fresh_sqrt5_cross_coface_cap_dual_support_20261017.md). It has 2,188 simplex rows, 336 actual quad rows, and ratio 12.3573995146. It has no unselected merged quad rows or hidden restored simplex rows.

The [expanded orbit/integrality cap package](fresh_sqrt5_orbit_integral_cap_capacity_20261017.md) supplies 815 signed capacities. The script checks each positive coefficient against the inherited common actual-provenance map, fixes every quad mass, and raises the required containing-simplex masses. Nineteen simplex masses increase; the positive face support remains 2,524 rows. Recounting the entire frozen model gives ratio 12.3453216395 before local packing separation.

The subsequent check retains all 82,806 inherited face obligations and all **94,722 inherited signed columns**, as well as every supplied cap. Increasing vertex masses preserves those exact checks. The 815-column family includes the earlier 241 cap constraints, so it is not 815 additional independent capacities on top of those 241.

## Fresh complete local covers

Full literal simplex-corner reconstruction gives **6,763 corner occurrences, 4,169 block orbits, and 648 contributing forests**. Every block orbit is enumerated in the actual weighted rank-13 forest. Pair and one-port multiplicities are checked against the existing face rows.

For each forest, the checker enumerates the complete disjoint-block count-vector domain and solves one small local packing problem. Resource guards are explicit: at most one million cached vector entries, 30 seconds of enumeration, and 20 seconds for its local solve. There are **no open or resource-limited cases**.

The output contains **558 exact local optima and 90 exact cover brackets**, with **3,507 explicit disjoint-packing mixture terms**. A bracket still has a rigorously verified cover at its reported upper mass; no exact optimality claim is made for those 90 cases. All 648 covers dominate the newly computed orbit marginals. No old fixed-mass cover was silently reused.

The separation identifies **84 exact new charge columns** violated before the local repair. Their orbit charges, finite face coefficients, and exact excesses are exported. After increasing vertex masses to cover the fresh mixtures, the resulting ratio is exactly

\[
\frac{44535347825512544000000000000000}
{3617363685101424635575857650731}
=12.311548327014249\ldots.
\]

## Why this covers the full local family

Each saved mixture consists of literal disjoint port blocks. Averaging an entire packing under the weighted-forest automorphism group preserves disjointness and makes its marginals uniform on each complete orbit. The saved mixture dominates every required orbit marginal. Thus any nonnegative charge on those blocks, normalized by its maximum total on one disjoint packing, has dual load at most the repaired vertex mass. This establishes the full nonnegative local-family obstruction, even where the saved cover is only an upper bracket rather than an exact local optimum.

Additional unoccupied block orbits cannot invalidate the certificate: their required marginals are zero, and the displayed packings remain admissible in the larger abstract domain. Actual geometric compatibility of different simplex or prism configurations is a separate restriction, not asserted by these abstract mixtures.

## Remaining scope

The new face masses are checked against all inherited finite joint cuts. The complete three-, four-, and five-port joint families have **not yet** been rescanned at these changed masses and expanded actual-quad support. A further full joint separation can therefore expose additional constraints. This package proves neither an optimal global model value nor a stronger compact arithmetic dimension bound.
