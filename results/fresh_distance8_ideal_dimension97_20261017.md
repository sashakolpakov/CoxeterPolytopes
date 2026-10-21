# Finite distance nine and ideal feature payments: dimension97

**Complete internally checked computer-assisted proof.** Primary generation
and fresh saved replay pass. Separate mathematical/source review and
independent exact scalar reconstruction pass, with no identified gap.
The argument below gives dimension at most **97** for every finite-sided,
finite-volume real hyperbolic Coxeter polytope. There is no arithmetic,
compactness, simple-cusp, cusp-factor-rank or facet-count restriction.
The previously completed compact-facet ambient ceiling **79** is unchanged.

This reduction inherits the complete [dimension108 proof](fresh_selective_distance9_dimension108_20261017.md)
and its [independent internal review](fresh_selective_distance9_dimension108_20261017_review.md).
It changes only the ideal-vertex payments. The [source](../scripts/fresh_distance8_ideal_dimension97_20261017.py)
and [certificate](fresh_distance8_ideal_dimension97_20261017.json) bind the
new local geometric and all-rank counting inputs, then check the scalar transfer.
The separate distance-eight core explorations are not prerequisites.

## Corner rule and payment of every face

At finite vertices retain the completed rule: a pair receives weight1
when its graph distance is at most9 or its connected component has order
at most17. The total finite-vertex cost remains at most9n.

At an ideal vertex, use its actual irreducible affine factors. A feature
is a vertex of degree at least3 or either endpoint of a nonsimple bond.
Give a pair in one factor weight1 in any of these cases:

* its distance is at most8;
* the factor has at most17 vertices;
* its distance is9 and at least one mark is within14 graph steps of a feature.

Every other pair has weight0. This is intrinsic under diagram
automorphisms. The [factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md)
identifies each actual ideal face corner with a pair in its varying factor,
including at nonsimple cusps.

The inherited triangle proof supplies two corner paths of length at most8.
Both remain paid. For a quadrilateral, the dimension108 proof supplies a
corner of distance at most9 or with marked component of order at most17.
A finite corner, a small component, or a distance-at-most8 corner remains
paid immediately. In the sole remaining case it is an ideal distance9
corner in a factor of order greater than17. If it is near a feature, the
new rule pays it. Otherwise the [universal ideal-corner replacement](fresh_distance8_ideal_corner_replacement_20261017.md)
supplies an adjacent corner path of length at most8, which is paid.
That lemma has a separate [independent exact and mathematical review](fresh_distance8_ideal_corner_replacement_20261017_review.md).

Thus every triangle receives at least2 and every quadrilateral at least1.
Nonnegative weights suffice for all larger polygons. The total assigned
weight is at least the sum D5 of5 minus the number of sides over all two-faces.

The replacement lemma uses complete one-root finite/affine extension
lists, weak positive cusp heights and an exact four-root determinant.
It needs no exclusion of small affine diagrams. Its cycle branch has
uniform normalized determinant gap5/81; its feature-free tree branch
has gap1/9. Positive-semidefinite spectator elimination preserves these
contradictions. Its geometric scope is essential; a numerical height or
rank cutoff would not justify this transfer.

## Ideal cost and all-rank shape reduction

The [affine feature-budget theorem](fresh_affine_distance8_feature_payment_20261017.md)
proves that a connected factor with m vertices has at most

    min(8m, m(m−1)/2)

paid unordered pairs. Its proof covers every rank by explicit path, cycle
and fork counts. In particular affine D reaches8m in all orders m>=43;
the estimate does not rely on a finite rank sample.

At a cusp with affine ranks r_i, write d=n−1 and V=product(r_i+1).
Each pair in factor i occurs in exactly V/(r_i+1) face-corner fibers.
The total cost divided by V is therefore at most

    sum_i min(8,r_i/2).

For the link H-polynomial h(t)=product_i(1+t+...+t^(r_i)), put
G=sum_j h_j min(j,d−j). The inherited [shape inequality](fresh_cusp_link_shape_curvature_20261017.md)
gives normalized ideal curvature (5G−dV)/(2V), together with a strictly
positive global slack. The finite curvature coefficient is

    beta_n=n(n−6)/8                   if n is even,
    beta_n=(n−1)(n−5)/8               if n is odd.

The inherited all-rank convex-order transfer applies to cost
sum min(8,r_i/2), with cap16. It reduces all integer rank partitions of
d to copies of16, at most one smaller positive part and at most one
larger part. The cost is preserved and the curvature can only decrease
under this reduction. Therefore it suffices to check those shapes.

## Exact dimension conclusion

For link dimensions d=97 through107, the checker enumerates all957
reduced shapes. Each has strictly positive curvature minus cost. The
minimum is112101/1419857, at ranks(1,16,16,16,16,16,16) and ambient
dimension98. It also checks beta_n>9n in every ambient dimension98
through108. The resulting lower bound for D5 exceeds the upper bound
from corner payments, a contradiction.

The completed108 theorem excludes every larger ambient dimension.
Consequently the new ceiling is97. At ambient97, the profile consisting
of six ranks16 has margin−185352/1419857. This proves that this particular
per-cusp scalar comparison cannot lower the ceiling to96; it is not a
realized finite-volume polytope.

Finite weights were retained, so the completed compact-facet ceiling79,
compact-arithmetic ceiling21 and unrestricted compact ceiling29 are
unchanged. No new finite-volume construction is asserted.

## Verification scope

The producer computes the link polynomial by exact integer convolution,
then checks the curvature again by a separate inclusion-exclusion moment
formula. It verifies all957 shapes, the finite coefficients and the97
boundary profile, and binds the completed108 proof and both new local
lemmas. Final generation and fresh saved replay both terminated exit0
with `PASS_FINITE9_IDEAL8_FEATURE_DIMENSION97`. The separate independent
checker reconstructs all957 shapes, minima and the boundary using a
different shape enumeration, sliding-window convolution and cumulative-tail
inclusion-exclusion identity. Its exact values agree throughout. Two separate
internal mathematical/source reads confirm the face-payment transfer,
finite9n cost, ideal factor multiplicity, cap16 reduction and inherited108
tail. The linked independent review records its final bindings and saved
replay; it does not rerun the inherited108 high-core census.
External mathematician refereeing, publication priority and proof-assistant
formalization are separate unestablished matters.
