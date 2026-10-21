# Direct attachment restrictions resolve297 supporting cases

The [joint-contact frontier](fresh_distance7_joint_contact_domains_20261017.md)
left548 high-overlap contexts open. The present
[source](../scripts/fresh_distance7_direct_contact_features_20261017.py)
and [certificate](fresh_distance7_direct_contact_features_20261017.json)
resolve297 of them, leaving251. The combined frontier now has6,452
structural certificates:5,858 whole-component bounds at most17 and594
protected-feature outcomes. It is still not a complete finite7 theorem.

The improvement uses already proved necessary contact restrictions
directly. Once no further supporting wall can touch a given core node,
that node's degree is fixed. An extra positive-inverse/end-set condition
is unnecessary for preserving its leaf status.

## Radical sign constraints remove impossible contacts

Each reconstructed full core has one negative direction. Its root span
therefore contains a timelike vector and is nondegenerate as an actual
Lorentz subspace. Every Gram radical vector z consequently gives an
actual zero relation among its displayed wall normals n_i:

    sum_i z_i*n_i = 0.

Let e be an additional supporting wall. The inherited joint domain
contains every possible core node with e*n_i<0. All these products
are nonpositive, so if e contacts a node i with z_i nonzero, it must
also contact some node k whose coefficient has the opposite sign.
Otherwise dotting the relation with e gives a strict nonzero sum.

Write d_L and d_R for graph distances to the two opposite black pairs.
A contact with both i and k creates a marked route of length at most

    min(d_L(i)+d_R(k), d_R(i)+d_L(k)) + 2.

If every opposite-sign candidate gives a route of length at most7,
node i cannot be contacted in the branch under consideration. The
empty opposite-sign set also forbids contact. The source records every
such node and every opposite-sign path bound, removes these contacts,
and repeats until stable. Deletions within one round all use the domain
before that round, so there is no circular exclusion. This argument
does not bound labels or the number of additional supporting walls.

There is also an amplitude constraint. Every core wall occurs with the
additional supporting white in an actual edge diagram, whose Gram is
strictly positive definite. For physical norm2 roots, a nonzero contact
therefore has magnitude theta with1<=theta<2. When all remaining contact
nodes have gauge1, contact at i contributes at least abs(z_i) to its
sign side of the radical balance. The opposite sign side is strictly
less than2*sum_k abs(z_k), with k ranging over every allowed opposite-sign
node. If the former bound is at least the latter, contact at i is
impossible. The strict upper endpoint excludes equality too. This does
not discretize or impose a finite upper bound on Coxeter labels.

The source records this second reason separately and repeats both
pruning tests to stability. All162 pending contexts with a nontrivial
radical have their domains reduced;28 use an amplitude exclusion.
Zero-coefficient nodes remain permitted unless a separate necessity
excludes them. Surviving contacts are not assumed to occur in an actual
chamber.

## Features preserved without an inverse splitting condition

Fix one actual corner and its displayed connected component containing
both marks. If no node of that component lies in the final contact
domain, no new supporting wall can attach to it. An indirect attached
component would need a first direct attachment. The displayed marked
component is therefore whole.

Otherwise, any displayed leaf outside the domain stays a leaf. A
displayed simple leaf whose degree-two neighbor is also outside the
domain retains both its simple bond and the degree of its neighbor.
Existing branches and nonsimple bonds also persist. Paths already
present provide upper bounds for the distances in any completion.

The source records a protected marked leaf at distance at most9 from
the other mark. When there is a protected simple leaf with a fixed
degree-two neighbor and the displayed order is at least10, the whole
component must be finite: the irreducible affine classification has no
such configuration at order at least10. Otherwise an existing branch
or nonsimple bond within5 of a mark supplies the inherited ideal-feature
alternative. These are exactly the protected finite-leaf9/ideal-feature5
interfaces of the frozen supporting theorem, proved here from direct
contact exclusions rather than an inverse estimate.

No positive definiteness is imposed on an ideal corner, and no bound
on the order of a protected-feature completion is asserted.

## Complete outcomes and limits

All548 pending contexts, on106 distinct literal core matrices, are tested.
Direct contact restrictions alone certify226. Radical pruning supplies
another71, giving the following new certificates:

|Outcome|Contexts|
|---|---:|
|Whole marked component, order at most16|27|
|Protected leaf with forced finite completion|84|
|Protected leaf with finite/ideal-feature alternative|186|
|Still unresolved|251|

The whole-component orders are9 (5 cases),11 (2),12 (8),13 (2),15 (2)
and16 (8). Thus the new bounded cases add27 to the previous5,831,
while the270 new protected outcomes add to the previous324. All12 nonempty
cross assignments remain covered by earlier certificates. The251 unresolved
cases are explicit and have not been declared impossible.

This result controls arbitrary supporting attachments. It does not
yet choose final finite7 pair weights, verify their selected embeddings,
close low/disjoint cases, or count a new global budget. The dimension
ceiling remains65; no new polytope is constructed.

## Verification

The final generation16147 and fresh saved replay91542 both terminated
exit0/PASS. The earlier path-only generation18228 and replay88664 also
passed, but the final pair verifies the complete amplitude refinement.
A separate full source/report mathematical review passed the actual
radical relations, simultaneous pruning, strict finite-contact amplitude
endpoint, direct degree preservation and all-rank affine-classification
step. The reviewer did not independently replay the arithmetic. The calculation rebuilds
native Grams, radicals and actual corner data, then checks every pruning
step and every protected-feature output. The joint corner/edge domains
and the upstream6,703-context census are frozen inputs, not rerun here.

    python3 scripts/fresh_distance7_direct_contact_features_20261017.py
    python3 scripts/fresh_distance7_direct_contact_features_20261017.py --verify
