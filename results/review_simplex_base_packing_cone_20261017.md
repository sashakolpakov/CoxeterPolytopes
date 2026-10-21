# Independent exact review of the complete small base-packing cones

**PASS.** The stated inequalities completely describe the nonnegative marginal vectors covered by disjoint-subset packing mixtures on two through five vertices. No gap was found in the LP-duality reduction, half-integrality proof or exhaustive enumeration.

The independent [review checker](../scripts/review_simplex_base_packing_cone_20261017.py) and [saved review](review_simplex_base_packing_cone_20261017.json) use standard-library Fraction arithmetic and no producer imports. Generation and a fresh `--verify` reconstruction both passed. They bind:

* Producer source SHA-256 `eb2835e57d52cb3a87a52dc9382d97b9eb90360376f30a7461c455329f451fa2`.
* Producer data SHA-256 `0fb80cf8f7fb3fd3718cd4192ec5a6f776110382cb894bd9f11dc6f463ceebe3`.

See the [producer report](fresh_simplex_base_packing_cone_20261017.md) for the domain and geometric motivation. This review does not assert that abstract packings are geometrically realizable.

## Completeness of the reduction

Every counted base contains at least two vertices, so at most two bases can be disjoint when the ground set has at most five vertices. The minimum mixture mass covering nonnegative demands has a finite feasible primal: singleton-base packings cover every demand. Its dual has nonnegative charges bounded by one on each base and charge sum at most one on each disjoint pair. Ordinary finite-dimensional LP duality therefore applies without numerical assumptions.

The downward-closed packing domain makes covering and exact marginals equivalent: excess block mass can be removed by deleting that block from an appropriate part of the mixture. This retains a packing and does not increase total mass. Thus all listed dual inequalities are sufficient as well as necessary, including at capacity zero.

For five vertices, each triple is compatible only with its complementary pair, and every four- or five-set is isolated in the compatibility graph. Because demands are nonnegative, an optimum can set the latter charges to one and each triple charge to one minus its complementary pair charge. This leaves the ten pair charges with bounds between zero and one and pairwise-sum constraints on the Petersen graph. The reduced objective coefficients may be negative; the argument does not assume otherwise.

At an extreme point, examine the graph of tight equalities between strictly fractional coordinates. No tight equality joins a strictly fractional coordinate to a coordinate fixed at zero or one. A bipartite fractional component permits a small alternating two-sided perturbation while maintaining every tight equality and every formerly strict inequality. Hence it cannot occur at a vertex. A nonbipartite component contains an odd cycle; its equations force all its coordinates to one half. This proves half-integrality. Conversely, integral coordinates fixed by their bounds and half-coordinate components containing odd cycles provide full active rank. The entire vertex set is therefore found on the three-value grid.

For at most three ground vertices all bases mutually intersect. For four vertices the only compatibility edges are the three complementary pair matchings. The one, one and eight corresponding inequalities follow directly. These reductions, together with the complete five-vertex enumeration, establish the full nonnegative marginal domain; a finite grid search by itself would not do so.

## Separate arithmetic reconstruction

The reviewer independently enumerated all \(3^{10}\) pair-charge assignments. For each of the 1,959 feasible assignments it formed the full active constraint matrix and computed its rank by rational elimination. It did **not** use the producer's odd-cycle vertex predicate. Exactly 294 assignments have rank ten, with half-coordinate counts

\[
0:76,\quad5:12,\quad6:60,\quad7:90,\quad8:45,\quad9:10,\quad10:1.
\]

Their pair-coordinate set agrees exactly with the producer's lifted rows. The complement-triple and large-subset lift identities were checked separately.

The integer packings were independently reconstructed as set partitions with singleton parts discarded. This gives 2, 5, 15 and 52 packings on two, three, four and five vertices. The resulting literal packing sets match the saved domain. Every one of the 304 inequalities was evaluated on every relevant packing. The reviewer computed affine rank directly by subtracting one tight packing from the others, obtaining full facet dimension for every row; this independently confirms the claimed rank25 for all 294 five-vertex facets.

Finally, the four pairwise-intersecting blocks of mass one third were checked against all 80 subset-floor conditions, including the single-vertex conditions. They satisfy these conditions while their intersecting-family mass is four thirds. Thus the full domain is strictly stronger than that family.

## Scope

This is an independent exact arithmetic and mathematical review of the abstract packing domain. It supplies no ambient dimension conclusion. Its application to prism incidences additionally needs the reviewed disjoint-base theorem and correctly symmetrized coface demands; those are separate hypotheses, rather than conclusions of this small enumeration.

**FINAL/FROZEN.** Generation, fresh saved replay, syntax, source hashes and local links passed. This is an internal review, not external refereeing or formal proof-assistant verification.
