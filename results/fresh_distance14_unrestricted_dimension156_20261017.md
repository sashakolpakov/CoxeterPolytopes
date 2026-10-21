# An unrestricted dimension156 bound from distance14 payments

**Complete internally checked computer-assisted proof.** Every
finite-sided, finite-volume real hyperbolic Coxeter polytope has dimension
at most156. There is no arithmetic, compactness, cusp-simplicity,
minimum-factor-rank, factor-count or facet-count hypothesis. The proof
combines the inherited bound163 with the new exclusion of dimensions
157 through163. No new finite-volume example is constructed.

Two ingredients are new: a distance14 quadrilateral payment and a
convex-order reduction of arbitrary cusp-rank partitions to a small,
explicit family. The argument uses no unproved fine-payment extension
at small cusp factors.

## The radius-seven obstruction cannot occur on both sides

Use the actual two-face star W plus its side roots. As in the
[distance15 proof](fresh_simple_cusp_distance_bound_20261017.md), W plus
one side is positive definite, and W plus an adjacent pair is positive
semidefinite, including ideal corners with nonsimple ambient cusps.
Opposite closed sides are disjoint, so W plus either opposite pair is
indefinite. These are actual-face facts, not assumptions about an
arbitrary abstract diagram.

Let U={u1,u2}, V={v1,v2} be opposite pairs. Choose inclusion-minimal
indefinite diagrams L in W+U and M in W+V. Both marked deletions of
each core are positive definite because they lie in actual edge
supports. The frozen complete minimal-diagram classification and its
384 finite-marked-pair calculations show that every core node has
distance at most7 from its marked pair. Orders at most5 have the
stronger elementary radius3 bound. At radius7 the only possibilities
are the following, with labels from the frozen catalogue:

| Core | Black pair | Unique farthest node |
|---|---|---|
| L23 |0,1|8|
| L25 |1,2 or1,3 or2,3|9|

Suppose no adjacent black pair in the quadrilateral can be joined by
a white path of length at most14. Then L and M are disjoint: a shared
white node would give a path of length at most7+7. They cannot be
orthogonal because both are indefinite and the ambient negative index
is one. Any bond between them must join nodes of radius7 on both
sides; otherwise it gives a white path of length at most14. There is
only one radius7 node in each listed case, so the unique possible cross
bond joins these two farthest nodes. All other cross products vanish.

The L23 case is impossible even before considering the other core.
In its finite deletion of black0, the inverse Gram diagonal at node8
is9/2 (and after deleting black1 it is4). An exterior root of square2
with a nonright contact there has magnitude at least1. Entrywise
nonnegativity of the finite inverse makes its Schur cost at least9/2,
contradicting positive semidefiniteness. The exterior farthest root is
white, so this extension is in fact contained in an actual finite edge
diagram. Additional nonpositive contacts can only increase this cost.

It remains to join two literal L25 diagrams at their respective node9.
Each has a simple A3 tail7--8--9. These two tails belong to an actual
corner and have no other cross bonds. Eliminating the two interior
nodes on each side leaves diagonal4/3 at each proposed bridge endpoint.
The bridge magnitude must therefore be at most4/3. A nonright Coxeter
product is1 or at least sqrt2, so the bridge has label3, magnitude1.
This includes all finite labels, parallel and dotted possibilities;
there is no angle cutoff.

The resulting20-node diagram has negative index two. Each L25 has a
fork center0 with short arms of lengths1 and2 and a long arm of
length6. Eliminating the short arms leaves effective diagonal5/6 at
each fork center. The two centers are joined by a13-edge simple path.
Eliminating its12 internal nodes leaves

    [ -7/78  -1/13 ]
    [ -1/13  -7/78 ],

whose eigenvalues are−1/6 and−1/78. Every eliminated block is positive
definite. Thus the actual star would have two negative directions,
which is impossible. The checker also derives this Schur complement
directly from the literal20-root Gram.

Consequently every quadrilateral has an adjacent black pair connected
by a white path of length at most14. The triangle argument is unchanged:
a minimal indefinite core has order at most10 and contains all three
black roots, giving two marked paths of length at most8 in the
corresponding actual corners.

## Unrestricted vertex costs

Give weight1 to a corner whose marked roots have graph distance at
most14 in its full finite or affine vertex diagram; give weight0
otherwise. The paths just proved occur inside the actual corner
diagrams, so triangles receive at least2 and quadrilaterals at least1.
Larger polygons require only nonnegative payment. Extra factors at an
ideal vertex do not change the marked pair's varying factor.

The [all-rank classical graph formulas](fresh_simple_cusp_distance_bound_20261017.md)
bound pairs at each distance by the component order; small exceptional
types satisfy the bound directly. Thus a finite vertex costs at most14n.
For an affine factor of rank r and order r+1, both the distance bound
14(r+1) and the total-pair bound r(r+1)/2 apply. The
[actual cusp-fiber multiplicity](fresh_cusp_factor_payment_framework_20261017.md)
therefore gives, for a cusp link with V=product(r_i+1),

    cost/V <= sum_i min(14,r_i/2).

This holds for every irreducible affine type and any number of factors.

## Concentrating rank partitions without changing this cost

For centered discrete uniforms X_r, the
[shape identity](fresh_cusp_link_shape_curvature_20261017.md) is

    kappa/V = 3d/4-(5/2)E|sum_i X_(r_i)|,    d=sum_i r_i.

For integers1<=a<=b, replacing (a,b) by(a−1,b+1) increases their sum
in convex order. Here is a direct all-rank proof. Before centering,
the probability at j of U_a+U_b, with each U_r uniform on0,...,r, is

    min(j+1,a+1,b+1,a+b-j+1)/[(a+1)(b+1)].

After the replacement the denominator is a(b+2), smaller by b−a+1.
The new probability is strictly larger on the tails j<a or j>b,
and strictly smaller on the central interval a<=j<=b. On that
interval the probabilities are respectively1/(b+2) and1/(b+1).
Both distributions are symmetric about the same mean(a+b)/2 and
have total mass one. Symmetrizing any convex test function makes it
nondecreasing with distance from that mean. Subtracting its value at
the crossing radius proves that the change in expectation is
nonnegative. This proves convex order, including a=1, where X_0 is
the constant zero.

Adding the other independent uniforms preserves the comparison,
by applying it to the convex function x -> E|x+Y|. Consequently
the replacement decreases normalized curvature at fixed d.

Apply this transfer only to two ranks strictly below28, stopping when
one becomes0 or the other28, or to two ranks strictly above28,
stopping when the smaller becomes28. In each case
sum min(14,r_i/2) is unchanged. Ranks0 can be discarded. Repeating
terminates, for example because the integer sum of squared ranks
strictly increases while d stays fixed. Every initial rank partition
therefore has curvature-minus-cost at least that of a profile

    (28,...,28), optionally one a with1<=a<=27,
    optionally one b>=29.

There is at most one rank of each optional kind. This is a complete
comparison for all partitions, not a census of sampled cusp shapes.

## Exact exclusion of dimensions157 through163

For each d=156,...,162 the checker enumerates every profile of the
displayed reduced form with sum d. It computes the exact curvature
from polynomial convolution. A separate inclusion-exclusion formula
for the truncated first moment independently verifies each result.
Every reduced profile satisfies

    kappa/V > sum_i min(14,r_i/2).

The smallest margin over the entire checked range is
1476616/20511149, attained at d=156 with ranks(28,28,28,28,28,16).
Every other rank partition inherits the inequality from concentration.
Also beta_n>14n throughout n=157,...,163.

The actual face payments give

    D5 <= 14n v_f + sum_ideal_v V_v sum_i min(14,r_vi/2).

The unrestricted shape theorem gives

    D5 >= beta_n v_f + sum_ideal_v kappa(L_v) + 5(n-1).

These inequalities contradict each other for every n in that range.
The already completed [unrestricted163 proof](dimension163_proof_status_20261017.md)
excludes all higher dimensions. Hence n<=156.

The distance14 comparison alone does not exclude n=156: at d=155,
the local link profile(28,28,28,28,28,15) has margin−4172075/82044596.
This is a limitation of the payment comparison, not an ambient example
or a claim that the dimension156 bound is sharp.

## Verification scope

The [checker](../scripts/fresh_distance14_unrestricted_dimension156_20261017.py)
and [certificate](fresh_distance14_unrestricted_dimension156_20261017.json)
bind the frozen catalogue, the radius proof, actual factor geometry,
unrestricted shape inequality and inherited163 proof. New arithmetic
checks the forbidden Gram configurations, the all-rank probability
identities and every concentrated profile in the necessary seven
dimensions. No previous large catalogue is re-enumerated.

Generation and a fresh saved replay pass all924 concentrated profiles,
the physical Gram obstructions and direct dependency hashes. Separate
construction-agent and frontier-agent mathematical/source reads pass
the actual-face geometry, unique farthest-node bridge, all-label bridge
gate, convex-order transfer and complete partition reduction. Neither
read claims to repeat the inherited classification census.

The [independent review](fresh_distance14_unrestricted_dimension156_20261017_review.md)
records a separately implemented rational-arithmetic check of the
literal Gram obstructions and all924 reduced profiles, together with
its saved replay and written mathematical review. Its arithmetic does
not import this producer's implementation. The written geometric steps
and convex-order proof are not formally verified by source-integrity
checks. No unresolved step is currently identified. External refereeing,
proof-assistant formalization and publication-priority claims remain
unestablished.
