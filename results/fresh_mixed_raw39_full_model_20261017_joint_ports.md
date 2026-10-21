# Selected joint-port constraints at the new mixed-raw39 masses

The complete represented local-simplex covers and the17 selected joint
simplex/prism families still admit an exact witness above12. After the
joint repair its ratio is

    562393570581318125000000000000 / 45153865815029704235917954789
      = 12.45504809898519...

This follows the [new local-cover repair](fresh_mixed_raw39_full_model_20261017_local_cones.md),
whose ratio was12.46109928685114. All5059 face masses remain unchanged;
only21 vertex allowances increase. The [source](../scripts/fresh_mixed_raw39_full_model_20261017_joint_ports.py)
and [saved certificate](fresh_mixed_raw39_full_model_20261017_joint_ports.json)
perform no optimization.

## Actual geometric domain

The scan reconstructs all2886 positive maximal-simplex/edge rows and their
5373 complete literal block orbits. It also reconstructs519 actual
quadrilateral classes, contributing3504 actual corner occurrences:
461 have the inherited explicit provenance map, and58 additional
raw23/raw25/raw39 rows have exactly one saved source context. Each such
context is checked against its canonical full star, demand and spectator
forest. No representative is chosen from a merged row with multiple
possible contexts.

The remaining1654 positive merged quadrilateral rows receive zero new
coefficients. This is valid underpayment for the nonnegative charge
families; it is not complete geometric coverage of those rows. The
explicit provenance map and all selected literal corner descriptions
are saved, so the separate orbit averages cannot silently choose
incompatible representatives.

The represented corner forests total821. For every forest, all286
three-subsets,715 four-subsets and1287 five-subsets of its13 ports are
enumerated and grouped under the same weighted-forest automorphism
invariant. Simplex intersection sizes and the numbers of occupied
quadrilateral star pairs are averaged on these common orbits. This
retains exact intersection sizes, including where a five-port facet's
coefficient is not monotone in that size.

## Complete scan of the17 selected families

The inherited proved families are the strengthened three-port inequality,
the four-port inequality, and all15 nonnegative facets of the five-port
necessary domain. Their original capacities and exact rational
normalizations are preserved.

| Subset size | Literal subsets | Orbits | Inequalities | Violations |
|---|---:|---:|---:|---:|
|3|234806|63482|63482|7|
|4|587015|130860|130860|10|
|5|1056627|208252|3123780|51|

Thus3,318,122 exact inequalities are tested and68 genuinely violated
columns are exported. Each exported load is independently recounted
from its face coefficients and the unchanged face-mass vector. Taking
the maximum required allowance at each forest raises21 allowances.
This is a simultaneous repair of all17 families, not separate
normalizers or independently selected provenance distributions.

For a simplex block of size k, summing its exact-intersection histogram
over all subsets gives `binom(k,h) binom(13-k,s-h)`. For a quadrilateral
star with e occupied pairs, the corresponding positive-edge histogram
is `binom(e,z) binom(12-e,s-1-z)`. Both identities are asserted for every
literal feature. The common subset-orbit and exact histogram digests
are saved. Singleton maximal edges contribute zero to all these
intersection-at-least-two features.

The predecessor directly checked all433753 finite-model columns,
including96176 zero-budget signed transfers and52 signed incidence
columns with positive vertex cost. Since alpha is fixed and every
vertex allowance increases or stays unchanged, all those checks and
the complete local covering mixtures remain valid. No stale-alpha
cover is reused: the predecessor rebuilt them at these exact masses.

## Verification and next missing geometry

Primary generation (session80491) and fresh saved replay (session89243)
both finish with exit0 and PASS. The replay recomputes the entire exact
scan without an optimizer. The new result does not prove feasibility at target12 or
an arithmetic dimension theorem. It covers the17 named joint families,
not every possible joint local condition, and leaves the1654 unexpanded
quadrilateral rows explicitly outside the new geometric coefficients.

A separate read-only support ranking identifies a concrete further
cohort: raw24 has15 explicitly represented and30 uniquely represented
positive classes, hence360 actual prism corners. The30 new unique
classes have total demand-weighted mass0.0983888376, the largest such
unexpanded single-family cohort. Its seven-node star has only finite
labels2,3,5 and dotted square `(23+9*sqrt(5))/2`, of norm31, so it is
compatible with the existing exact quadratic-extension projection
method. No raw24 pair census or new capacity conclusion is claimed here.
