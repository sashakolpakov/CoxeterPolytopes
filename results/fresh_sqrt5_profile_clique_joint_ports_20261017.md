# Complete joint three-, four-, and five-port checks after 851 cap constraints

An exact witness of ratio **12.25714671503906 > 12** survives all 851 cross-coface cap constraints, the complete nonnegative local simplex-block packing family, and the 17 specified joint-port families. This remains an obstruction for a partial incidence model, not a dimension theorem or a global optimum claim.

The [source](../scripts/fresh_sqrt5_profile_clique_joint_ports_20261017.py) and [saved certificate](fresh_sqrt5_profile_clique_joint_ports_20261017.json) retain exactly the face masses from the [final 851-cap local-cover update](fresh_sqrt5_profile_clique_local_cover_update_20261017.md). There are 2,524 positive face rows: 2,188 simplex stars and 336 literal actual quadrilateral provenances. No representative is reselected, no merged quad mass is substituted, and no optimizer is invoked.

## Complete common-orbit checks

Every simplex occurrence contributes its actual port block. The 4,169 complete weighted-forest block orbits are averaged against the same marked subset orbits used for quadrilateral corners. Each of the 2,294 actual quad-aggregate corners has a square-pair star with one hub and up to four leaves. For a single edge either hub choice gives exactly the same subsets containing that edge. The script reconstructs every literal quad corner from its saved full core and orthogonal spectators.

For a subset of size m, simplex features record the exact intersection-size histogram h=2,3,4,5. Quad features record exactly z=1,2,3,4 occupied star edges. Nonmonotone simplex coefficients are applied to these joint histograms, without separate minimization. Every subset orbit is enumerated in each contributing rank-13 forest.

| Subset size | Literal subsets | Subset orbits | Exact inequalities | New violated columns |
|---|---:|---:|---:|---:|
| 3 | 185,328 | 51,166 | 51,166 | 21 |
| 4 | 463,320 | 106,351 | 106,351 | 17 |
| 5 | 833,976 | 170,261 | 2,553,915 | 60 |

The families are the strengthened three-port inequality 2s+t+q<=3, the four-port inequality r+q_ge2<=3, and all 15 exact facets from the [five-port necessary domain](fresh_sqrt5_joint_five_port_domain_20261017.md). Their hypotheses restrict the selected quadrilateral aggregates to the actual maximal dotted-prism family. Arbitrary all-finite product aggregates are not silently included.

For each simplex block of size k, summing its h-intersection coefficient over all m-subsets gives binom(k,h)binom(13-k,m-h). For a quad corner with e square edges, the z-edge sum is binom(e,z)binom(12-e,m-1-z). Every such identity is checked exactly. Full domain and histogram streams are hashed; a fresh replay reconstructs those streams and all identities. Every newly violated column is exported with exact rational face coefficients, and its load is independently recounted from those coefficients and the same alpha vector.

## Union and scope

The scan checks 2,711,432 inequalities on 648 contributing forests. It exports 98 newly violated columns and increases 18 vertex masses. The union uses the per-forest maximum across every family. Outside the contributing support all new loads are zero.

Only vertex masses increase. Therefore every inherited matrix column, all 94,722 signed constraints, the 851 cap constraints, and the revalidated complete local packing covers remain valid. The exact final ratio is

\[
\frac{283748127467729879563707541448224505846026541755517058125904000000000000000}{23149606842804745538681805678996675190943173077937679324802469205844370541}.
\]

Generation and fresh saved replay both pass. Parent and frontier agents separately read the source and mathematical accounting and found no issue; those narrow reviews did not duplicate the full arithmetic scan. The cap families have a separate [independent Fraction-arithmetic and mapping review](fresh_sqrt5_orbit_and_orientation_caps_20261017_review.md).

A separate standard-library Fraction check also recomputes the per-vertex maximum union, all 98 exported loads, the unchanged face masses, and the final objective/normalizer ratio.

These 17 families do not exhaust all possible joint geometric compatibility constraints. No new global LP, feasibility claim, or improved ambient dimension bound follows from this certificate.
