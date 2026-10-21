# The arithmetic witness passes every disjoint-base packing capacity

The current compact-arithmetic witness passes the **complete abstract
disjoint-base packing cone** in all34,520 represented coface slices.
There are510,656 inequalities in total;65,069 have to be evaluated on the
870 slices with positive demand. All remaining slices have zero demand
and nonnegative capacity. There are no violations, and the maximum
capacity utilization is exactly1.

The exact ratio remains

    66908102699388132000000000000000 /
    5419878674792788039126804769429
      = 12.344944732908832... > 12.

No face or vertex masses change, and no optimizer is called. This result
exhausts a specified family of necessary inequalities on the saved
averaged slices; it is neither an ambient dimension theorem nor a proof
that the witness describes an actual chamber.

## Complete domain and exact reconstruction

The [labelled coface calculation](fresh_sqrt5_labelled_coface_base_packing_20261017.md)
reconstructs all6,016 bases of1,858 actual prism provenances. It preserves
the black/white coloring of the simplex star, simultaneous coface/subset
actions and every identical finite spectator copy. Its164,219 subset-floor
inequalities have no violation, but those inequalities alone do not
exhaust disjointness.

The [complete packing theorem](fresh_simplex_base_packing_cone_20261017.md)
supplies the missing full domain. A simplex has at most five vertices,
so there are at most26 nontrivial base subsets and52 disjoint packings.
The five-vertex covering dual reduces to a ten-variable Petersen-graph
polytope. Its294 vertices are exhaustively found by an exact
half-integrality argument and enumeration. Their lifted inequalities,
and the smaller domains on2--4 vertices, are necessary and sufficient
for an abstract packing mixture at the stated capacity.

The [scan](../scripts/fresh_sqrt5_complete_coface_base_cones_20261017.py)
and [saved data](fresh_sqrt5_complete_coface_base_cones_20261017.json)
import no producer arithmetic. They reconstruct every nonzero block
marginal from the saved actual-provenance coefficients and the common
face masses. They check its total against the masses of the source bases
and reconstruct capacity as coface-orbit multiplicity times simplex mass.
They then evaluate every full-cone inequality with exact integer scaling.

| Simplex dimension | Positive coface slices |
|---|---:|
|1|382|
|2|101|
|3|172|
|4|215|

For each positive slice the saved result includes a maximizing inequality,
the exact minimum abstract packing-cover mass, the capacity and the slack.
The zero-demand slices pass all nonnegative inequalities identically.
The total source-base mass is checked globally as well as slice by slice.

## What this rules out, and what remains

Every nonnegative linear capacity that follows solely from disjointness
of bases in one of these averaged coface slices is implied by the complete
domain. Adding more such capacities cannot invalidate this fixed witness.
That statement includes all intersecting-family and odd-cycle constraints,
not only the earlier subset-floor inequalities.

It does not test whether disjoint bases have compatible Coxeter diagrams,
whether separately feasible coface slices arise simultaneously from one
simplex star, or whether their other faces satisfy additional global
incidence restrictions. Those are genuinely additional geometric inputs.
It also does not complete any face classes omitted by the larger model.
The next useful strengthening must add one of these inputs; another scan
or optimization over the same disjointness domain cannot lower this witness.

Generation and a fresh saved replay pass for the domain and this scan.
The [domain's independent checker](review_simplex_base_packing_cone_20261017.md)
reconstructs all packings from set
partitions and identifies its vertices by active-matrix ranks, separately
from the producer's odd-cycle test. The [accompanying independent review](fresh_sqrt5_complete_coface_base_cones_20261017_review.md)
also checks the scan's normalization and scope; it does not repeat the
geometric base-to-coface reconstruction. These are internal checks, not
external refereeing or proof-assistant formalization.
