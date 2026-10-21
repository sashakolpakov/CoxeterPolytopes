# Complete labels and subsets of the residual distance-seven cross positions

The48 residual positions in44 contexts from the
[cross-position domain](fresh_distance7_high_cross_domain_20261017.md)
have now been closed under all native labels and all simultaneous subsets.
There are54 nonempty assignments. Exactly12 pass the necessary native-core
tests; all twelve add a single simple bond and have signature(17,1,0).
Thirty-eight assignments have at least two negative directions and four
fail strict positive definiteness of an actual edge.

The [source](../scripts/fresh_distance7_high_cross_completions_20261017.py)
and [certificate](fresh_distance7_high_cross_completions_20261017.json)
retain the label exclusions and each fully specified assignment. The6,691
zero-extra-edge base contexts remain separate. Arbitrary supporting whites,
low/disjoint cores and quadrilateral payments are still open obligations.
This is a complete necessary cross-contact domain, not a finite7 payment
theorem or an improved dimension bound.

## Every label and simultaneous subset

At46 positions an inserted label4 already supplies a nonnegative witness
of square at most zero in an actual edge. Increasing its native magnitude
or adding other nonpositive contacts preserves the strict-edge obstruction.
Only label3 can therefore remain at those positions.

The other two positions belong to contexts4545 and6003. Their base graphs
are simply laced. At each position, inserting the rational magnitude3/2
gives a nonnegative strict negative-square witness in an actual corner:
respectively -1/40 and -1/16. For every native label at least5, the physical
norm-two contact magnitude is at least(1+sqrt(5))/2>3/2. Monotonicity therefore
excludes all those labels, including parallel and dotted contacts. Only
labels3 and4 require enumeration at these two positions.

Each context has at most two residual positions. The source enumerates
every nonempty subset, and every combination of its surviving labels.
All other extra entries are zero, as forced by the preceding monotone
exclusions in an otherwise unpaid case. This produces all54 assignments.
The empty subsets belong to the predecessor's base contexts.

For each fully specified Gram, the checker tests the four actual edge
principal matrices for positive definiteness, the four corner matrices
for positive semidefiniteness, and the full core for Lorentz index at most1,
stopping as soon as a contradiction is found. Every survivor passes them all.
The full-core exclusion is used only at this stage, when every entry is
fixed; it is not assumed monotone under unspecified additional contacts.

## Exact inertia

Where possible, a positive diagonal congruence converts the native Gram
to a rational matrix. Otherwise the source works exactly in Q(sqrt(2)),
with label4 represented by -sqrt(2). The sign of a+b*sqrt(2) is decided
from rational signs and a^2-2*b^2. Symmetric elimination uses a nonzero
one-dimensional pivot or, when all diagonal entries vanish, a hyperbolic
two-dimensional pivot. The certificate saves the permutations, exact
pivots and signature. No floating-point eigenvalue test is used.

The twelve surviving base-context indices are:

    1678,1963,1994,2205,2997,3093,
    3104,4460,5262,5275,5281,5304.

The certificate saves their native graphs and literal marks. Passing
these necessary tests does not assert a geometric realization or a
finite-volume chamber.

## Verification

Initial generation79286 and final generation6951 terminated exit0/PASS.
The final generation binds the completed cross-position report. Fresh
saved replay91164 also terminated exit0/PASS. A separate internal full
source/report review checked the exact field signs and congruence steps,
monotone all-label gates, subset completeness and limited scope. Its
saved-schema/count check confirmed all54 outcomes and the12 listed
survivors; it did not independently recompute witness or inertia arithmetic.

    python3 scripts/fresh_distance7_high_cross_completions_20261017.py
    python3 scripts/fresh_distance7_high_cross_completions_20261017.py --verify
