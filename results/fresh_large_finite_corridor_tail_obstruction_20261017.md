# Deep finite corridors force short-path alternatives at64 corners

The [source](../scripts/fresh_large_finite_corridor_tail_obstruction_20261017.py)
and [exact certificate](fresh_large_finite_corridor_tail_obstruction_20261017.json)
prove a short-path alternative for a particular large-component completion at64 corners in the
remaining [251-context supporting frontier](fresh_distance7_direct_contact_features_20261017.md).
They do not delete or halve any payment yet.

Suppose the entire displayed corner is one connected simple path, its two
marked roots are distance8 apart, and its complete old core is nonsingular
Lorentzian. In each of the64 cases selected below, if its actual whole finite
component has both marks at graph distance at least20 from every endpoint,
branch, and endpoint of a nonsimple bond, the full actual supporting star
has an opposite-black path of length at most7. Equivalently, such a deep
corridor is impossible in the branch with no such short-path alternative
anywhere in the full supporting star. Absence of a short path in the saved
old core alone is not this additional branch hypothesis. The result holds
in all ranks for finite A, B and D, with no bound on the number or labels
of additional supporting walls.

## Complete subdomain and its limitations

The frozen joint-domain theorem left548 contexts. The direct-contact theorem
resolved297, leaving251. This checker inventories all four corners of each
remaining context, so there are1,004 literal corners. Classification is in
the following order; these rows are disjoint.

|First applicable classification|Corners|
|---|---:|
|Core has a radical|52|
|Marks disconnected in the displayed corner|483|
|Displayed corner has other component blocks|367|
|Marked component is not a simple path|16|
|Marked distance is not8|22|
|Eligible single-path corner|64|

The64 eligible corners occur in64 distinct still-open contexts. They were
also all still open before the direct-contact refinement. None is counted
as a new fully resolved supporting context: the conclusion concerns only
the proposed deep-corridor completion at that particular corner, and gives
a short7 alternative rather than unconditional geometric impossibility.

The single-block requirement is essential to this packet. The source asserts
that the path contains every old core node except the two omitted black
roots. It does not infer this from connectedness of the two marks. Thus
there is no other displayed corner component for a new endpoint root to
join. Corners with such spectator blocks remain outside this argument.
The native matrices, corner maps, complete contact/census provenance and
earlier structural classifications are frozen inputs. Their parent census
is not repeated here.

## Forced white paths exclude arbitrary extra contacts

Normalize every wall root to squared norm2. All64 selected literal cores
already have this gauge. A native simple bond has product−1. Define a
feature of a finite A/B/D diagram to be a vertex of degree other than2,
or either endpoint of a bond labeled more than3.

Assume both marks have feature distance at least20. The path joining them
has length8. Starting at either mark in the outward direction, away from
the other mark, there are20 successive vertices of the same simple
corridor. They are all supporting whites: the two omitted black roots do
not belong to this actual corner, and there are no further selected black
roots along either outward arm. The induced diagram of these20 whites is
A20. A feature at distance exactly20 is harmless: the segment through it
still has20 vertices and simple internal bonds. A doubled incoming bond
would place its other endpoint at distance19, contrary to the hypothesis.

An omitted black root occurs with all these supporting whites in its actual
edge diagram. That diagram is strictly positive definite. The inverse of
the norm2 A20 Cartan matrix is entrywise positive, with diagonal

    J[d,d] = d(21−d)/21,       1 ≤ d ≤ 20,

where d is outward distance from the marked black. Consequently:

* At every distance3 through18, J[d,d]>2. Even one unit contact from the
  omitted black would have projection square greater than2, impossible.
  Arbitrary other nonpositive contacts only increase the square, by
  positivity of the inverse.
* At distance2, J[2,2]=38/21>1. A label at least4 has squared contact
  magnitude at least2 and therefore also has projection square greater
  than2. A nonzero contact there must be the unit contact of label3.

These exclusions hold for either omitted black, with every other contact
allowed. No label cutoff or sampled support assumption is used.

The old path's outward endpoints have distances a and b from their nearest
marks, with1≤a,b≤7 in every selected case. Continue each endpoint until
outward distance17. The resulting new simple tails have lengths

    m_left = 17−a,       m_right = 17−b,

at least10 each. A first tail root contacts its old endpoint with product−1.
It can have an additional old-core contact only if a=1 (or b=1), putting
it at outward distance2, and then only a simple contact to one omitted
black. Contacting both omitted blacks would create an opposite-black path
of length2, already in the retained short-path branch. Every later tail
root is orthogonal to the entire old core. The whole actual finite diagram
is a tree, so roots of the left and right new tails are also orthogonal to
one another.

## Exact Lorentz-complement obstruction

Let G be the old core Gram matrix and J=G^{-1}. The source reconstructs
each selected native Gram, checks signature(one negative, no zero), and
recomputes its saved inverse. Its actual span is a nondegenerate Lorentz
subspace; the ambient orthogonal complement is positive definite.

Write x and y for the positive unit contact vectors of the first left
and right tail roots. Put

    q_left = x^T J x,     q_right = y^T J y,     t = x^T J y.

All other tail roots have zero old-core projection. The Gram matrix of the
tail residuals in the positive complement must be positive semidefinite.
Eliminating all but the first root of each simple tail leaves the necessary
two-by-two Schur matrix

    [ (m_left+1)/m_left − q_left          −t                  ]
    [             −t            (m_right+1)/m_right − q_right ].

The diagonal correction follows directly from the endpoint inverse
(m−1)/m of A_{m−1}. All denominators are positive. The source lists every
possible endpoint/omitted-black support, checks any retained short7
alternative, and evaluates this matrix exactly over the rationals.

All112 joint entry profiles either give a short7 alternative in the full
star or fail a necessary positive-semidefinite condition:

|Reason|Profiles|
|---|---:|
|Negative Schur diagonal|79|
|An entry already makes a marked path of length at most7|25|
|Nonnegative diagonals but negative determinant|8|

Thus every eligible deep corridor supplies the short7 alternative.
Only under the full-star no-short7 branch hypothesis does the result
contradict its existence. Some profiles routed to a shortcut have positive
Schur matrices; those profiles are not claimed geometrically impossible.
This proof controls long tails without imposing a bound on the total
actual component order or on additional walls elsewhere.

## Conditional candidate counts, not an implemented payment rule

The motivating invariant rule would halve the weight of a distance8 pair
whose two marks both have feature distance at least20. Direct graph counts
give the following numbers of such pairs in a whole component of order n:

    A_n: max(0,n−48),    B_n: max(0,n−49),    D_n: max(0,n−50).

The source independently counts these pairs for18≤n≤75 and checks the
formulas. For all n, the canonical path indices give the formula directly:
in A_n the pair is(i,i+8) with20≤i≤n−29; in B_n, with its doubled bond
at0–1, it has21≤i≤n−29; and in D_n, with spine0 throughn−3 and two
short leavesn−2,n−1, it has20≤i≤n−31. Empty ranges contribute zero.
The current cost formulas are cross-checked against the frozen
[current finite-budget rule](fresh_d16_d17_d14_coupon_forest_budget_20261017.md);
its dynamic program is not rerun. At order65, the hypothetical effects are:

|Whole component|Candidate pairs|Current same-component cost|Hypothetical cost|
|---|---:|---:|---:|
|A65|17|486|477.5|
|B65|16|486|478|
|D65|15|487|479.5|

Including A and B matters because their current costs also exceed480.
These numbers are only a target. This packet has not proved coverage for
multi-block or disconnected old corners, the already-certified supporting
outputs and their actual retained pair maps, low/disjoint cases, or later
universal payment interfaces. In particular, a whole-component bound at
most17 is not automatically a paid corner under the current refined
weights. No new coupon, finite-forest budget, payment system or dimension
bound is asserted. The unrestricted dimension ceiling remains65.

## Verification

Final corrected generation60502 terminated exit0/PASS (terminal chunkcd998f),
and separate fresh replay43615 terminated exit0/PASS (terminal chunk3ac6be).
Earlier runs preceded the explicit full-star branch statement and are not
used as the final replay. The new calculation reconstructs selected native cores, inverses,
corner maps, the A20 inverse and all112 rational tail-profile exclusions;
it does not rerun the inherited all-label overlap census.

The construction agent completed an independent full source/report
mathematical read, including the A20 all-label argument, single-block
first-entry completeness, positive complement, two-tail Schur formula,
all-rank A/B/D counts and frozen cost provenance. That review identified
and corrected the earlier statement that treated a short7 alternative as
an absolute geometric contradiction. Its saved-only Fraction check found
that14 of the25 shortcut-routed matrices pass the two-by-two PSD test;
the final theorem explicitly allows these short7 alternatives. The review
did not repeat the64-core/112-profile enumeration or any parent census.

Root completed a full source/report mathematical read, followed by a
targeted read of the corrected source/data statement and a reread of the
corrected report through the Schur proof. Its final review confirms only
the stated full-star short7 alternative and the exclusion under the
explicit full-star no-short7 branch. This scope includes the frozen
current-cost comparison and all-rank index ranges, but no independent
112-profile arithmetic or upstream census replay. Direct metadata audit
e1ff6f passed the source AST, all18 bindings and all4 local report links.

    python3 scripts/fresh_large_finite_corridor_tail_obstruction_20261017.py
    python3 scripts/fresh_large_finite_corridor_tail_obstruction_20261017.py --verify
