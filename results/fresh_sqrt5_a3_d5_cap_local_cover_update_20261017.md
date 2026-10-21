# Exact witness repair after the A3/D5 cap constraints

The repaired partial-model witness has exact ratio **12.128700376798244>12**, down from12.12966051709919. Both newly identified A3/D5 inequalities are satisfied, together with all inherited cap and matrix constraints. This is an exact obstruction for the stated incidence model, not an optimized bound or dimension theorem.

The [source](../scripts/fresh_sqrt5_a3_d5_cap_local_cover_update_20261017.py) and [certificate](fresh_sqrt5_a3_d5_cap_local_cover_update_20261017.json) bind the [two exact cap witnesses](fresh_sqrt5_a3_d5_cap_witnesses_20261017.md). They raise containing-simplex rows81801 and78789 to meet capacities2 and1 respectively. Both containing simplices are maximal edges, so their curvature demands are zero. Every quadrilateral mass is unchanged.

The [previous complete local covers](fresh_sqrt5_triple_groundfield_cap_local_cover_update_20261017.md) supply3,710 literal packing terms. All their orbit keys, literal-mask membership, disjointness and old marginal coverage are rechecked. Exactly two singleton terms cover the changed marginals, producing3,712 valid terms over the same4,169 block orbits and648 forests. Two vertex masses increase above the [previous joint-port baseline](fresh_sqrt5_triple_groundfield_cap_joint_ports_20261017.md).

The support remains2,188 simplex stars and336 actual quadrilateral rows. All94,722 inherited signed matrix columns pass directly, as do2,056 recorded cap constraints:851 original profile/clique columns,555 original multispectator columns,93 triple-A1 columns,555 fixed-ground-field successor columns and the two A3/D5 constraints. Keeping the stronger555 successors while omitting their weaker predecessors would leave1,501 family columns; this certificate retains both for explicit preservation.

Simultaneous automorphism averaging of the saved whole disjoint packings gives complete uniform marginal coverage, so every nonnegative local simplex-block charge remains covered at the new face masses. This is a statement about the entire local packing family, not just an inherited finite collection of weights. It does not establish geometric realizability of the incidence witness.

Generation and a fresh saved reconstruction replay both terminate with PASS. Exact objective, normalizer, face and vertex masses, all cap loads and complete packing mixtures are saved. No optimizer is invoked. The complete joint three-/four-/five-port check at these changed masses is a separate successor.
