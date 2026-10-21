# Incremental complete local covers after cross-orientation cap tightening

The cross-orientation cap family leaves an exact incidence-model witness of ratio **12.26680604413596 > 12**. This does not establish a dimension bound or an optimal model value.

The [source](../scripts/fresh_sqrt5_tightened_cap_local_cover_update_20261017.py) and [saved certificate](fresh_sqrt5_tightened_cap_local_cover_update_20261017.json) fix all 336 actual quadrilateral masses and minimally raise 24 containing-simplex masses. The positive support stays at 2,524 faces: 2,188 simplex rows and 336 actual quad rows. The strengthened [815-cap family](fresh_sqrt5_cross_orientation_cap_capacity_20261017.md) is checked through the same inherited provenance map.

The starting [complete local covers](fresh_sqrt5_orbit_cap_restricted_local_cones_20261017.md) contain 3,507 explicit literal packing terms on 648 forests and 4,169 block orbits. Every term is transported using its exact marked-forest/block key, checked for literal disjointness, and recounted at the old required marginals. The new face masses require 46 additional singleton packing terms. No new block orbit is needed. Forty-five vertex masses increase.

These mixtures dominate every new orbit marginal. Averaging each entire packing under the finite-forest automorphism group preserves disjointness and makes the required within-orbit marginals uniform. Thus the certificate still covers **all nonnegative local simplex-block charges**, not just a saved finite collection. No local or global optimizer is invoked; the new cover is a rigorous upper cover, with no optimality claim.

All inherited positive-cost columns and 94,722 signed columns are rechecked at the final masses. Both the old and strengthened 815-column packages are checked, giving 1,630 cap checks. The exact ratio is

\[
\frac{113909244873926154974846443172151213124664770186272000000000000000}{9285974235190540396580922599007305055491739131005554509601078103}.
\]

Generation and fresh saved replay both pass. Replay rebuilds all coefficient maps and literal cover checks without optimization. Complete joint three-, four-, and five-port scans at these changed masses remain a separate obligation. The following distinct-profile clique successor may strengthen the masses again; this intermediate certificate is preserved unchanged.
