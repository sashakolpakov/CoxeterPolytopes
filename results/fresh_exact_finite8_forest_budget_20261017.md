# Exact maxima of the current finite-corner rule

The unchanged finite-corner payment rule has maximum cost **548 in rank69**
and **553 in rank70**, over all abstract finite Coxeter forests of those
ranks. The first is attained by `A1+A17+A17+A17+A17`; the second by
`A2+A17+A17+A17+A17`. These are universal upper bounds at actual finite
vertices, with equality witnesses in the larger class of abstract vertex
diagrams. No ambient realization of a maximizing diagram is asserted.

The [source](../scripts/fresh_exact_finite8_forest_budget_20261017.py) and
[certificate](fresh_exact_finite8_forest_budget_20261017.json) give the exact
maximum in every rank1 through75. This is a refinement of the existing8n
budget. It changes no face payment, adds no finite coupon, and proves no
new ambient dimension bound by itself.

## The complete unchanged rule

The inherited [base rule](fresh_finite8_leaf9_half_coupon_budget_20261017.md)
pays one on a same-component pair if its distance is at most8, or its
distance is9 with a whole-component leaf mark, or its whole component
has rank at most17. Its exact component cost is binom(r,2) for r<=17,
8r-34 for A/B at r>=18, and8r-32 for D at r>=18. These are all-rank graph
identities using the finite Coxeter classification, not sampled formulas.

Retain all [broad half-payments](fresh_finite8_a17_d11_coupon_budget_20261017.md):
between distinct whole components, one of rank at most5 with arbitrary
mark and the other D11 through D16 at its three fork marks, A13/A15 at
its middle, or A17 at its two specified path positions. Also retain the
[A4/D13 neighbor halves](fresh_finite8_a4_d13_neighbor_coupon_budget_20261017.md),
and the [D12 short-leaf upgrade and A2/B10 halves](fresh_finite8_d12_b10_coupon_budget_20261017.md).
Every unordered marked pair is counted once.

## Domination by one representative per rank

For the maximization problem each connected finite type may be replaced
by the following type of the same rank without decreasing its total cost:

|Rank|Representative|Reason|
|---|---|---|
|1 through5|A_r|Same complete base cost and broad small-component eligibility; A1/A2/A4 retain every additional small-type coupon|
|6 through9|A_r|Same complete base cost; no cross-component coupon|
|10|B10|Same complete base cost, plus its A2 coupon|
|11 through16|D_r|Same complete base cost; three broad large marks dominate the one A13/A15 middle mark; D12/D13 retain their additions|
|17|A17|Same complete base cost, plus its two broad large marks|
|18 and above|D_r|Base8r-32 dominates the A/B base8r-34; no cross-component coupon|

The exceptional finite types and arbitrary I2(m) occur only among the
small ranks already covered. The domination argument uses counts of all
eligible marked pairs, not a graph embedding or a geometric deformation.
Every coupon coefficient is nonnegative. It is therefore valid to replace
one component at a time: no other component's contribution decreases.
Conversely, every listed representative is itself a finite Coxeter type,
so a maximum attained among these representatives is also attained in the
full abstract finite-forest class.

Let a_i count representative A_i components for1<=i<=5, let
`a=sum_i i*a_i`, let t count D11 through D16, and let u,b10,b12,b13 count
A17,B10,D12,D13 respectively. For a representative forest the cost is exactly

    sum_components base(rank)
      + a*(3t/2+u) + a1*b12 + a2*b10 + a4*b13.

There are three broad eligible marks on each D and two on each A17,
each paid at one half against all a small marked nodes. The D12 upgrade
increases each of its two short-leaf weights by one half against each A1,
giving a1*b12. A2/B10 has two pairs of weight one half, giving a2*b10.
The two A4 endpoints against the distinct D13 neighbor give a4*b13.
These additions involve disjoint pair roles and count no pair twice.

## Complete exact maximization

Enumerate all count vectors of the eight large coupon ranks10 through17
whose total rank is at most75. There are1,608 such vectors. For fixed large
counts the cost of each remaining connected component is additive: its
rank is1 through9 or at least18, and its entire cross-payment contribution
is fixed by the large counts.

Only the tuple `(t,u,b10,b12,b13)` matters for this residual additive
problem. There are377 distinct tuples. Double every cost to use integers.
For each tuple k let w_k(r) be the doubled residual component cost, including
all its cross-payments, and define

    F_k(0)=0,
    F_k(s)=max_(allowed r<=s) [F_k(s-r)+w_k(r)].

Allowed ranks are1,...,9 and18,...,75. A rank greater than s cannot occur.
This Bellman identity is an induction proof of completeness: removing any
one component of a residual forest leaves a lower-rank residual forest;
conversely every maximizing predecessor may be extended by that component.
The positive component ranks make the recursion finite. The certificate
stores all76 values and an attaining last rank for every tuple.

For each ambient rank n, maximize the doubled large-component base cost
plus F_k(n minus their rank) over every large vector of rank at most n.
This covers every representative forest. The source also reconstructs an
attaining forest for each n and checks its full cost directly using the
displayed coupon formula.

## Maxima near the current frontier

|Rank n|Exact maximum|Slack below8n|One maximizing forest|
|---:|---:|---:|---|
|66|513|15|D15+3 A17|
|67|528|8|D16+3 A17|
|68|544|0|4 A17|
|69|548|4|A1+4 A17|
|70|553|7|A2+4 A17|
|71|559|9|A3+4 A17|
|72|566|10|A4+4 A17|
|73|574|10|A5+4 A17|
|74|578|14|A1+A5+4 A17|
|75|583|17|A2+A5+4 A17|

At ambient69 the finite curvature coefficient is544. The unchanged payment
rule can cost548, so even this exact abstract-forest budget exceeds the
available coefficient by4. Its obstruction is explicit: the four A17
components already cost544, and their half-coupons against A1 add4.
Thus simply sharpening the estimate of the same rule cannot make every
finite-vertex coefficient fit at69. A smaller geometric payment rule,
restrictions on actual vertex diagrams, or a global compensation argument
would be needed. This does not prove that the abstract obstruction occurs
in an actual finite-volume chamber.

At ambient70 the exact553 cost is7 below the finite curvature coefficient560,
whereas the older8n bound only gives equality. The existing dimension
argument may retain that older bound and its positive global additive
constant; this refinement is a separate theorem with no new ideal costs.

## Verification status

Primary generation and a fresh separate-process saved replay both terminated
exit0/PASS on all1,608 large vectors,377 residual tables and75 positive-rank
maxima. The [independent certificate review](fresh_exact_finite8_forest_budget_20261017_review.md)
uses no producer arithmetic. It verifies all885,950 Bellman upper-bound
alternatives and28,275 attaining predecessors,22,504 global upper bounds,
and all maximizing forests by direct marked-pair counting. It reconstructs
the large-vector domain by component number and nondecreasing multisets,
and separately checks type domination on228 finite-type profiles and51,984
ordered partner comparisons, including the symbolic arbitrary I2 family.
Its report records final generation and fresh replay.

A separate full internal source/report mathematical read passed the type
domination, exact coupon counts, complete Bellman and large-vector domains,
and the distinction between abstract attainment and actual realizability.
Finite-type classification and the unchanged geometric payment rule are
inherited inputs; the present computation does not establish any new face
coverage. No external refereeing, formalization, ambient realization or
new dimension theorem is claimed. Source, certificate and report are final/frozen.
