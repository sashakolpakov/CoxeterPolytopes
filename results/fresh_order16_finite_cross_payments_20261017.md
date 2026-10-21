# Reduced order16 ideal payments with finite cross-component half-payments

Every triangle receives at least2 and every quadrilateral at least1 under
the rule below. Finite-vertex cost is at most9n for n<=108. The normalized
ideal cost per affine factor of rank r is bounded by

    b(r) = r/2       for 1<=r<=14,
           7        for r=15,
           127/17   for r=16,
           217/29   for r>=17.

This is a complete local payment and cost statement; dimension transfer
is a separate argument. The [source](../scripts/fresh_order16_finite_cross_payments_20261017.py)
and [certificate](fresh_order16_finite_cross_payments_20261017.json) bind
the [completed preceding payment theorem](fresh_ideal7_finite_cross_payments_20261017.md),
the [order16 embedding census](fresh_order16_preferred_embeddings_20261017.md)
and the [universal antipodal cycle replacement](fresh_cycle16_half_coupon_transfer_20261017.md).

## The rule and reduction to preferred components

Retain the preceding finite rule: unit payment for connected marked
pairs at distance at most9 or in whole components of order at most17,
and half-payment between A1/D16 or A2/D14 at the stated D fork marks.
Add these finite half-payments:

* Between whole A1 and D15 components, at the D fork center or either short leaf.
* Between whole A2 and D13 components, at the D fork center or either short leaf.
* Between whole A1 and A15 components, at the unique middle node of A15.

Either A2 mark is allowed. No payment is assigned merely from a subdiagram
inside a larger component, and these disconnected rules are disjoint
from the connected unit rule and from one another.

Retain the preceding ideal rules at every order except16. At order16,
the A, B and C paid sets are the complete basic-plus-preferred-image
unions in the new census, of sizes112,107,103. For D, start with all120
pairs and remove the eight pairs at distance at least12 and the pair of
fork centers, which has distance11. The remaining set has111 pairs.
All rules are intrinsic under automorphisms of the whole affine factor.

Every triangle retains its two distance7 unit payments. Suppose a
quadrilateral has total payment less than1. It can have no adjacent
white-interior path of length at most9. At a finite corner such a path
pays directly. The inherited ideal replacement handles all larger
factors. In smaller factors all pairs are retained. At order16, all
B/C/D distances at most9 are retained. The sole A exception is the
antipodal distance8 pair; the new cycle lemma supplies another unit
payment or two sealed finite A1/A15 half-payments. Thus this exception
cannot occur in an unpaid face either.

The complete inherited no-short9 classification now applies unchanged,
including its cross-bond exclusions and arbitrary-support extension
bounds. Its preferred connected marked component is either finite,
or in an actual affine factor of order at most17. Finite corners and
orders at most15 pay directly; order17 retains the completed preceding
transfer. For order16, every preferred component must occur among the
260 induced embedding domains. All its images are paid except for the
nine D pairs handled next.

## The nine D15 exceptions and arbitrary additional support

The nine deleted D pairs occur only in inherited contexts
384,394,401,648,663,676,814,829,842. Each literal core consists of two
L25 diagrams with the common edge matched by8->9 and9->8. It has18
roots and Gram signature(17,1,0). The corner omitting roots3 and13
contains the entire16-node affine D15 factor F.

Every extra supporting white e is orthogonal to every root of F. Indeed,
at that actual ideal corner the positive affine kernel pairs to zero
with e and every individual Coxeter product is nonpositive. The only
remaining possible contacts with the core are at3 and13.

If both contacts are nonzero, the actual corner marked by3 and13 has
the white-interior path3-e-13 of length2 and is already paid. If exactly
one contact is nonzero, write its magnitude as t>=1 in the norm-two
normalization. The inverse of the18-root core Gram has contact block

    [ 4 -4 ]
    [-4  4 ]

on3,13. Therefore the orthogonal projection of e into the Lorentzian
core span has square4t^2>=4. The orthogonal complement of that span is
positive definite, so e itself would have square at least4, contradicting
its square2. This rules out a single contact without a label cutoff.

Consequently, in an otherwise unpaid face every additional supporting
white is orthogonal to the entire18-root core. The remaining white
diagram is positive definite, as it supports an actual hyperbolic face.
It cannot enlarge or join any of the finite core components. All
three other corners remain finite with their displayed whole components.
This proves support control in every ambient dimension; it is not an
observation restricted to a particular small chamber.

Exact rational corner reconstruction gives the same finite-payment
pattern as for the preceding order17 exceptions. Five cases already
have a connected finite unit payment. Cases384,401,814,842 instead
have two half-payments of the new A1/D15 or A2/D13 kinds. Every affected
quadrilateral therefore receives at least1. The certificate checks all
nine cases, the positive affine kernel, the Lorentz signature and the
inverse contact block exactly.

The separate cycle lemma treats all possible antipodal order16 A stars.
Its complete finite/affine extension analysis and determinant equality
force the two simple-leaf equality core, unless another unit corner is
already paid. Spectator equality then either yields a path2 payment or
seals the A1/A15 components. This universal step is essential to the
new112-pair A rule.

## Combined finite budget

Let a,b be the numbers of whole A1,A2 components; let t count all
D15,D16,A15 components and v count all D13,D14 components. Their ranks
are at least1,2,15,13 respectively. Compared with9 times their ranks,
the retained unit costs save at least9,17,24,35 per component.
All finite half-payments together cost at most(3/2)at+3bv; this
overestimates the A1/A15-center contribution, whose coefficient is1/2.
Thus

    cost <=9n+a(3t/2-9)+b(3v-17)-24t-35v,
    a+2b+15t+13v <=n.

For n<=108, t<=7 and v<=8. If t<=6 and v<=5, every term after9n is
nonpositive. If t=7, then v=0 and a<=3; the possible positive term is
at most9/2, while24t=168. If v=6,7,8, respectively, b<=15,8,2;
the b contribution is at most15,32,14, while35v is210,245,280.
In these cases t<=2, so the a term is nonpositive. This exhausts the
integer domain and proves the finite9n bound. Arbitrary other finite
components retain cost at most9 times their ranks.

The source additionally checks26,817 endpoints for ambient dimensions
87 through108. Its expression is affine in a, so checking the two
rank-constrained a endpoints for every b,t,v covers all integer counts.
These checks corroborate the analytic bound; no bound on ambient
dimension is inferred from them alone.

## Every-rank ideal envelope

Factors of order at most15 retain every pair, giving rank cost r/2.
At order16 the maximum of112,107,103,111 is112, giving112/16=7.
At order17 the unchanged maximum is127, giving127/17.

For order m>=18, the retained all-rank A cost is7m. The B/C/D formulas
are the previously proved piecewise linear counts, with successive
slopes9,8,7. Dividing by m, their largest values occur at their finite
breakpoints. The maximum is217/29, attained by D at m=29; B is bounded
by7+11/28, C by7+8/27, and A by7. The source checks every linear
branch endpoint and the increasing slack of each infinite final branch
against217m/29. This is an all-rank argument, not extrapolation from
the earlier finite BFS range.

The actual cusp factor-fiber theorem multiplies each factor's pair count
by V/(r+1), where V is the number of link vertices. Summing gives the
normalized ideal budget sum_i b(r_i), including at nonsimple cusps.

## Verification scope

The source uses exact rational linear algebra for the exceptional cores,
complete inherited embeddings for the preferred components, the separate
universal cycle theorem, finite-budget endpoints and all-rank linear
branches. Primary generation and a separate fresh saved replay both
terminate exit0 with the identical certificate. A separate full
mathematical/source read confirms the260-pattern domain, no-short9
reduction, nine D15 exceptions and inverse4 projection obstruction,
whole-component support argument, all new coupon types, combined finite
budget and every-rank envelope. That review does not duplicate the
260 embeddings, nine-core arithmetic or26,817 endpoint calculations.
The independent cycle package has its own generation and fresh replay.
The inherited complete
classification remains an explicit input; this package does not claim
an independent rerun of that entire census or external refereeing.
