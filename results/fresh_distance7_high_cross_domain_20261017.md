# Additional cross edges in the distance-seven overlap domain

All6,691 surviving base contexts from the
[complete identification domain](fresh_distance7_high_overlap_domain_20261017.md)
now have every possible extra cross-core edge tested. In6,647 contexts
the base Gram is forced: every additional edge either creates a marked
path of length at most7 or contradicts an actual edge's strict positive
definiteness. The remaining44 contexts have just48 individually possible
positions:40 have one position and four have two.

The [source](../scripts/fresh_distance7_high_cross_domain_20261017.py)
and [certificate](fresh_distance7_high_cross_domain_20261017.json) retain
all base contexts, every tested position and every exclusion witness.
They do not yet classify the labels or simultaneous subsets of the48
remaining positions, arbitrary supporting whites, low/disjoint cores,
or retained finite and ideal payments. No new dimension bound follows.

## Completeness and monotone exclusions

The two minimal indefinite cores are induced diagrams. Every entry with
both endpoints in either copy, including an entry involving a shared
white root, is therefore already fixed. The only unspecified entries
join a node exclusive to the left copy to a node exclusive to the right
copy. The source enumerates their full Cartesian product in every base
context. There is no restriction to white endpoints or connected overlaps.

Each nonzero native contact between norm-two roots has magnitude at
least1. Insert a simple bond at the proposed position. A shortest path
between the opposite marked sets has no internal marked node, after
truncating at the nearest marked endpoints. A path of length at most7
therefore gives the short-path alternative required by this reduction.
Additional bonds cannot remove that path.

Otherwise inspect every actual edge principal diagram, comprising all
supporting whites and one black root. Where possible, a positive diagonal
congruence makes the inserted native Gram rational. If a doubled-cycle
parity prevents that congruence, replace every nonzero bond by a unit
bond. This latter matrix bounds the native quadratic form from above on
nonnegative vectors. A saved nonnegative vector of square at most zero
contradicts the strict positive definiteness of an actual edge. Increasing
any bond or adding any other nonpositive contact only decreases its square.
Thus the witness excludes every higher label and every simultaneous
additional-contact assignment containing that position.

Actual corners are also tested with strictly negative nonnegative-vector
witnesses. Null vectors alone are never used to reject a PSD corner.
No context or remaining position is eliminated by this extra corner
test in the present data.

|Extra position classification|Occurrences|
|---|---:|
|Marked path of length at most7|356860|
|Strict actual-edge obstruction|29613|
|Not excluded by these tests|48|
|Total|386521|

The29,613 edge exclusions use24,313 distinct saved witnesses. The source
checks each exact witness when constructing it. The repeated-context
cache only shares identical literal matrices and marked sets; it does
not identify merely isomorphic configurations or discard orientations.

## Scope and verification

Generation79915 and separate saved replay70298 both terminated exit0/PASS.
A separate internal source/mathematical review checked the complete exclusive
cross-position domain, positive gauges and unit-bond comparison, shortest-path
truncation, strict corner witnesses and the stated limited scope. It reused the
inherited helpers and did not independently enumerate all positions. The frozen identification
census supplies the completeness of the6,691 base contexts. This successor
does not repeat that matching enumeration, but it reconstructs every union,
every possible extra position, every shortest path and every actual-edge
or corner test in its own domain.

The48 surviving positions are necessary possibilities, not asserted
geometric realizations. Their all-label/subset closure is a separate
successor. Even closing that domain would still leave the supporting-root
and full quadrilateral-payment obligations stated above.

    python3 scripts/fresh_distance7_high_cross_domain_20261017.py
    python3 scripts/fresh_distance7_high_cross_domain_20261017.py --verify
