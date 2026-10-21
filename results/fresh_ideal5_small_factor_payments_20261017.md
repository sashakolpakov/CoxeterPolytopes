# Complete small-factor census for the ideal5 rule

This package reconstructs every inherited preferred-component pattern that can
occur in an affine factor of order14 or larger, and records every labelled
embedding into the four classical targets of orders14–17. It also supplies
the all-rank cost of the basic ideal5 rule. Local face coverage and the global
dimension comparison are separate composition steps; no dimension claim is
made by this census alone.

The [source](../scripts/fresh_ideal5_small_factor_payments_20261017.py) and
[certificate](fresh_ideal5_small_factor_payments_20261017.json) preserve the
complete frozen [distance9 classification](fresh_selective_distance9_dimension108_20261017.md).
They do not substitute the still separate distance8 overlap domain.

## Exact domain and pair sets

At a cusp, a feature means a vertex of degree at least3 or an endpoint of a
bond with label at least4. The basic rule pays marked distance at most5, and
distances6–9 when a mark is within5 steps of a feature. This rule is intrinsic
under weighted-diagram automorphisms.

The complete949 high-core overlap cases have preferred marked components
with proved extension bounds at most17. Exactly571 have bound at least14:
134 of bound14,177 of bound15,180 of bound16, and80 of bound17. Each pattern
is reconstructed from its literal two cores, overlap map and preferred
corner. For a target of order m, only patterns with bound at least m and
order at most m are considered. The inherited explicit backtracking tries
every unused target vertex consistent with all previous edges **and
nonedges**. It therefore enumerates all induced labelled injections, not
only one graph isomorphism or one orbit representative.

The paid set is the union of the basic rule and every marked image, followed
by precisely these validated removals:

* In cycles, remove distances greater than5. Order14 uses the
  [new cycle14 transfer](fresh_cycle14_ideal5_transfer_20261017.md); orders15–17
  use the [ideal5 local theorem](fresh_finite8_ideal5_local_replacement_20261017.md),
  including its sealed cycle16 coupon.
* In affine D, remove the nine pairs across the two fork triples. Orders14
  and15 use the universal [D13](fresh_affine_d13_long_pair_transfer_20261017.md)
  and [D14](fresh_affine_d14_long_pair_transfer_20261017.md) transfers. Orders16
  and17 use the previously proved literal sealed-core replacements in
  [order16](fresh_order16_finite_cross_payments_20261017.md) and
  [order17](fresh_ideal7_finite_cross_payments_20261017.md). The checker matches
  the complete affected case-index sets to those saved proofs, not merely
  the number of removed pairs.

The exact retained costs are:

| Factor order | A cycle | B | C | D |
|---:|---:|---:|---:|---:|
|14|70|87|84|82|
|15|75|100|95|96|
|16|80|107|103|111|
|17|85|111|108|127|

All emitted maps, marked images, basic pairs, additional pairs, removed pairs
and final pairs are saved. Thus the table is accompanied by literal pair
sets rather than only scalar totals. Exceptional affine types have order at
most9 and remain wholly paid in the intended composition.

## All-rank basic cost

For every classical affine factor of order m≥18, the exact basic costs are

| Family | Cost |
|---|---|
|A|5m|
|C|5m−15 + sum of min(m−k,14), k=6,7,8,9|
|B|5m−6 + sum of min(m−k,13), k=8,9,10,11|
|D|5m+3 + sum of min(m−k,12), k=10,11,12,13|

Here is an arbitrary-rank derivation. A cycle has no features and exactly m
pairs at each distance1–5. For a tree, use a simple spine with M vertices.
At distance at most5 the spine contributes5M−15. A spine pair at distance k
has M−k possible starting positions. Paying the first a and last b positions
of a list of length L selects exactly min(L,a+b) positions, including when
the endpoint intervals overlap.

For C, M=m and the doubled terminal bonds make the first two and last two
vertices features. The radius5 neighborhoods reach the first seven and last
seven spine vertices. At each distance k=6–9, the paid starting positions
are the first seven and last seven, giving min(m−k,14).

For B, M=m−2, with two short leaves at the left end and a doubled bond at the
right end. The left branch neighborhood contributes the first six starting
positions; the right doubled end contributes the last seven. Thus the
spine extras are min(m−k−2,13), k=6–9. The two leaves contribute ten pairs
at distances1–5 and their mutual pair, changing the short count to5m−14.
They also contribute two pairs at each of distances6–9, all within one
step of the branch. Adding eight gives the displayed5m−6 constant.

For D, M=m−4, with two leaves at each end. Each branch gives six starting
positions, so the spine extras are min(m−k−4,12), k=6–9. The four leaves
contribute twenty short leaf-spine pairs and two same-end leaf pairs,
giving5m−13. Their sixteen extra leaf-spine pairs change this to5m+3.
Opposite-end leaf pairs have distance M+1≥15 and do not contribute. These
arguments cover every m≥18 and do not extrapolate from a finite sample.

The four minima change one at a time, so each tree family has five linear
branches, of slopes9,8,7,6,5. C changes at orders20–23; B at21–24; D at22–25.
The certificate saves all15 branches and their exact endpoint comparison
with163m/23. The final slope5 tails have increasing slack. The respective
normalized maxima are75/11 for C,160/23 for B, and163/23 for D; the common
bound163/23 is attained by D at order23. The cycle bound5 is smaller.

Consequently the intended factor-rank envelope, writing r=m−1, is

    r/2 for 1≤r≤12;
    87/14 for r=13;
    20/3 for r=14;
    111/16 for r=15;
    127/17 for r=16;
    163/23 for r≥17.

This is a counting envelope. Its use for actual cusp incidences additionally
requires the factor multiplicity theorem and complete face-payment transfer.

## Verification scope

Generation passed the complete571-pattern embedding census, all16 small
target pair sets, exact correspondence of the order16/17 removed contexts,
248 breadth-first graph checks and all15 infinite-tail branch comparisons.
A fresh separate-process saved replay also terminated exit0/PASS. The finite
graph checks corroborate the written spine proof. A separate complete
source/mathematical read found no gap in the inherited domain filter,
all-edge/nonedge backtracking, pair unions, exact deletion-case correspondence
or all-rank spine counts. That review independently derived the15 branches
and their maxima, but did not duplicate the571-pattern map census. Syntax,
direct source bindings and report links pass. Source, data and report are
final/frozen; no claim of external refereeing or formalization is made.
