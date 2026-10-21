# Endpoint marks give116 further deep-corridor alternatives

The [source](../scripts/fresh_large_finite_corridor_endpoint_obstruction_20261017.py)
and [certificate](fresh_large_finite_corridor_endpoint_obstruction_20261017.json)
inventory all188 endpoint-mark corners left outside the
[preceding principal-core lemma](fresh_large_finite_corridor_principal_obstruction_20261017.md).
Exactly116 pass the present principal-core filters. At each, a whole finite
A/B/D corridor with both marked roots at feature distance at least20 forces
an opposite-black path of length at most7 in the full actual supporting star.

Thus the deep corridor is excluded only under the explicit full-star
no-short7 hypothesis. A short path is an alternative payment route, not
an unconditional geometric contradiction. No pair weight is changed here.

## Complete endpoint inventory

The frozen prior inventory contains188 corners, in120 parent contexts,
whose displayed marked component is a simple path, whose marks are
distance8 apart, and where at least one mark is an endpoint of that path.
The source reconstructs each path and actual corner, then restricts its
Gram to the path and the two omitted black roots. All other old components
are dropped from this principal matrix.

|Principal-core outcome|Corners|
|---|---:|
|Physical norm2 gauge, nonsingular Lorentzian|116|
|Positive definite, signature(11,0,0)|30|
|Positive semidefinite, signature(10,0,1)|10|
|Positive semidefinite, signature(11,0,1)|12|
|Positive semidefinite, signature(12,0,1)|4|
|Nonunit rational gauge, outside the current filter|16|

Signatures list positive, negative and zero counts. The72 unselected cases
are not geometrically excluded. A positive principal core has a Lorentzian
orthogonal complement, so the residual-PSD argument below cannot be used
unchanged. The16 gauge cases may be approachable after rescaling; they
are not asserted to lack a physical norm2 normalization.

The116 selected corners lie in96 parent contexts. Each has exactly one
endpoint mark, so the final Schur matrices happen all to have order3,
although the source handles the complete at-most-four-root setup. Their
old corner diagrams have two component blocks in72 cases and three in44.
No old whole-component isolation is assumed.

## Complete contacts at the first two outward whites

Use physical wall roots of squared norm2. A feature is a vertex of degree
other than2, or either endpoint of a bond labeled more than3. Feature
distance at least20 supplies an induced outward A20 path of supporting
whites from each mark, as in the preceding lemma. Its inverse is positive,
with diagonal

    J[d,d] = d(21−d)/21,       1≤d≤20.

Every omitted black occurs with these whites in its actual edge diagram,
which is strictly positive definite. Any other nonpositive contacts only
increase the relevant projection square. The resulting all-label bounds
are:

* At distance1, only labels3 and4 are possible. For label5 or larger the
  squared contact magnitude is greater than5/2, and (20/21)(5/2)>2.
* At distance2, only label3 is possible, since (38/21)·2>2 excludes
  the squared magnitude2 of label4 and every larger label.
* At distances3 through18, even a unit contact has projection square
  greater than2; no omitted-black contact is possible.

For an endpoint mark, the first two new whites may therefore carry omitted
contacts. In the full-star no-short7 branch, at most one such contact occurs
across those two whites. Indeed, the same omitted black contacting both
would have projection square at least2 against their A2 block: its positive
inverse is (1/3)[[2,1],[1,2]], and two unit contacts already attain2. Strict
edge positivity forbids this, including equality. Different omitted blacks
contacting the two whites give an opposite-black path of length3; both
contacting one white give length2.

This proves the complete seven-profile domain for a distance0 end:

* Neither white has an omitted-black contact.
* The first white contacts either omitted black with label3 or4: four
  profiles.
* The second white contacts either omitted black with label3: two profiles.

At displayed endpoint distance1, only the first new white can carry an
omitted contact, giving three profiles: none or either omitted black with
label3. At distance at least2, the endpoint's unit attachment is the only
old-principal-core contact. The native attachment to the displayed path
itself is always simple, since it lies inside the assumed deep corridor.

These domains are overdomains of actual completions: no claim is made that
every listed contact pattern is realizable. Further short routes through
the first/second white bond are checked explicitly on the augmented graph.

## Exact residual Gram obstruction

The selected principal matrix G has one negative eigenvalue and no zero
eigenvalue. Its actual span is a nondegenerate Lorentz subspace, with
positive definite ambient orthogonal complement. Let J=G^{-1}. Continue
each outward tail through distance17; for displayed endpoint distance a,
its length is m=17−a.

At a=0 retain the first two tail roots and eliminate the later A_{m−2}
segment. Before subtracting projection to the principal core, the retained
block is

    [  2       −1       ]
    [ −1     m/(m−1)    ].

At a≥1 retain only the first root, giving the scalar (m+1)/m. The opposite
tails are orthogonal in the actual finite tree. Subtract X^T J X, where
the columns of X contain the positive contact magnitudes of retained roots
against the principal core. Each magnitude is1 or sqrt(2). Every remaining
tail root is orthogonal to that principal core by the A20 restrictions and
the simple-path incidence. The resulting matrix, of order at most4, must
be positive semidefinite.

Contacts to dropped spectator components remain unrestricted. Later tail
vertices may be old spectator whites. The immediate first root cannot
already be an old spectator adjacent to the displayed path, because it
would then lie in the same displayed component. None of these observations
changes the necessary Gram calculation in the smaller Lorentz span.

The source first inserts the retained whites and their forced bonds into
the old-core graph and finds an explicit shortest opposite-black route.
If a retained tail root is actually an old spectator root, this abstract
graph may display it twice; mapping the duplicate back to that actual root
still gives a walk of the recorded length, hence a path no longer than it.
An opposite-black route with other black vertices can be shortened to a
segment between consecutive black vertices of different colors, with only
whites in its interior. Thus a recorded route of length at most7 supplies
the required full-star alternative. Such profiles need not fail a Gram test.

For every other profile, a strictly negative principal minor of the residual
Schur matrix is saved. Exact arithmetic is provided by the
[Q(sqrt(2)) helper](../scripts/fresh_quadratic_sqrt2_exact_20261017.py): elements
are Fraction pairs(a,b) representing a+b·sqrt(2), signs are decided by
rational comparisons with a²−2b², and determinants use exact row elimination.
No floating tolerance or sampled label bound is involved.

## Exhaustive outcomes and remaining scope

All1,428 joint seven/three-profile choices are covered:

|Outcome|Profiles|
|---|---:|
|Explicit full-star short7 alternative|788|
|Strictly negative Schur principal minor|640|

The negative witnesses have order1 in216 cases, order2 in320, and order3
in104. There are no remaining PSD profiles in the selected subdomain.
The640 Gram obstructions and788 alternatives must not be conflated into
an unconditional nonexistence statement.

Together with the earlier164-corner result, this gives280 conditional corner
alternatives across174 parent contexts. It does not resolve those contexts
as a whole: an actual marked pair might occur in an untreated corner or
have a different component shape. The72 excluded endpoint cases, disconnected
marked components, earlier structural outputs with actual retained pair
maps, and all low/disjoint/later payment interfaces remain separate work.
There is no new finite-forest budget, pair-weight rule or dimension bound.
The unrestricted ceiling remains65.

## Verification

The bounded pilot found the same1,428 profile totals using a recursive
determinant implementation. Authoritative generation14274 terminated
exit0/PASS (terminal chunk4aa5d2). Separate fresh saved replay58737 also
terminated exit0/PASS (terminal chunk59ee69). Direct audit1c78fd passed
both source ASTs, all15 direct bindings and all4 report links. Root and the
construction agent each completed a full source/report/helper review.

The shared exact arithmetic helper independently compared row elimination
against recursive Laplace expansion on72 matrices of orders1 through6 and
checked exact near-zero Pell signs and PSD boundary examples. These checks
run in the producer too. They complement the algebraic source review; they
do not replace it. The inherited all-label supporting census and the
preceding164-corner enumeration are not rerun.

The construction agent independently read the full endpoint producer,
report and Q(sqrt(2)) helper. It checked the complete seven/three contact
domains from strict actual-edge A20/A2 tests, the necessary restriction
after dropping spectators, augmented routes mapping to actual walks,
the positive-complement residual condition, both tail Schur blocks and
exact field signs, row swaps and principal-minor tests. It confirmed the
explicit72 unresolved cases and the short7-alternative-only conclusion.
This review did not rerun the188-case/1,428-profile calculation, the72
helper determinant comparisons or the parent census.

Root read the full endpoint producer, report and exact sqrt(2) helper.
Its review checked seven-profile completeness, the A2 strict-equality
exclusion, the A20 actual-edge argument, old spectator identities mapping
augmented routes to actual walks, the consecutive opposite-color black
segment, both Schur corrections, and exact field signs, determinants and
negative-minor tests. It also checked that all72 filtered cases remain
unresolved by this packet. This was a source/mathematical read, without
independent116-case/1,428-profile arithmetic or parent-census replay.

    python3 scripts/fresh_large_finite_corridor_endpoint_obstruction_20261017.py
    python3 scripts/fresh_large_finite_corridor_endpoint_obstruction_20261017.py --verify
