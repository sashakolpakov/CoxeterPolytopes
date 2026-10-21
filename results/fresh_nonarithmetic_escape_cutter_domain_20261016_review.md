# Independent all-real escape-cutter obstruction in H9

**Every acute spacelike wall normal that cuts the specified escape is strictly
dotted to each of the three new caps, 49, 50 and 51.** Thus no single additional
wall can both cut that escape and meet one of these caps. This does not exclude
several walls that perform the two tasks separately.

The [independent reviewer](../scripts/fresh_nonarithmetic_escape_cutter_domain_20261016_review.py)
and [certificate](fresh_nonarithmetic_escape_cutter_domain_20261016_review.json)
use the [independently reconstructed complete 52-wall cone](fresh_nonarithmetic_affine_filtered_20261016_review.md),
not the producer's clipping implementation. The
[primary certificate](fresh_nonarithmetic_escape_cutter_domain_20261016.md)
is used for comparison and source binding.

Let the retained outward normals be n_i, the new caps a_j have square two,
and the strict spacelike escape be E. The cone

    K = {r : r.n_i <= 0 for every i, r.E >= 0}

contains every possible acute escape-cutting normal, including the relaxed
boundary r.E=0. A complete independent clip leaves 146 old generators and
adds 196 crossing generators. Every crossing is certified by a rank-eight
common-covector minor. There are 342 generators in total: one timelike and
341 spacelike. All 18,126 final inequalities are checked; exact crossing
identities establish the boundary zeros.

For every generator v and each cap a_j, the reviewer checks strictly

    (v - (v.a_j)/2 * a_j)^2 < 0,
    (v - (v.a_j)/2 * a_j).(H - (H.a_j)/2 * a_j) < 0.

Here H is the independently reconstructed future timelike interior point.
Its projection into a_j-perpendicular is timelike. Therefore all 342
projections lie in the same open future timelike cone of that Lorentz
hyperplane. Convexity makes every nonzero conic combination strictly
timelike there as well; zero projections and opposite-time cancellation
are excluded by the strict orientation checks.

Consequently any spacelike r in K satisfies

    r^2 - (r.a_j)^2/2 < 0,   r.a_j <= 0,
    hence r.a_j < -sqrt(2*r^2).

This is precisely a strict dotted contact. The proof imposes **no field,
angle-label, root-normalization, or old-center-preservation restriction**.
It requires retention of the specified 52 acute inequalities. It is an
obstruction to a single combined cutter and facet-closing wall, not a
nonexistence theorem for every finite-volume completion.

```sh
python3 scripts/fresh_nonarithmetic_escape_cutter_domain_20261016_review.py --verify
```
