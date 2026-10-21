# Order10 D opposite-fork weights27/32 with explicit finite compensation

This local theorem lowers each of the four opposite-fork leaf pairs of
whole affine D9, of order10, from1 to27/32. Its total falls45→355/8.
Together with the frozen [cycle10](fresh_order10_cycle_half_transfer_20261018.md),
[B10](fresh_order10_b_terminal_half_transfer_20261018.md) and
[C10](fresh_order10_c_endpoint_zero_transfer_20261018.md) rules, the
order10 family totals are85/2,44,44,355/8. The all-family envelope is71/16
per factor, a saving1/16 from9/2.

Write δ=5/32. The new finite payments are precisely:

* Weightδ for every whole A1..A8 endpoint/whole D9 short-leaf pair.
  A1 has one endpoint and D9 has two short leaves, giving increment2δ=5/16.
  A2..A8 have two endpoints, giving increment4δ=5/8 for each component pair.
* Restore the unique endpoint pair of whole finite A17, canonically(0,16),
  from0 toδ. Its total rises451/4→3613/32.

All these ordinary payments are intrinsic whole-component roles. There
is no spectator condition or multiplicity of witnesses. No other finite
same-component weight decreases. The source checks simultaneous actual
and representative base/cross domination, but does **not** compute a new
forest upper bound. Lowering rank9 also requires a new concentration
argument before any cusp or global consequence; this report proves no
new dimension bound.

The [producer](../scripts/fresh_order10_d_opposite_fork_fractional_transfer_20261018.py)
and [exact local certificate](fresh_order10_d_opposite_fork_fractional_transfer_20261018.json)
bind the current finite8/depth6 rule, the three preceding order10 proofs,
current finite rank15–17 role certificates and the inherited all-rank
spectator inputs. Frozen predecessor files are unchanged.

## Actual-star setup and complete end profiles

Use spine0,...,5, leaves6,7 at0 and8,9 at5. Select a=6,b=8;
fork swaps cover the other three changed pairs. The positive affine
kernel ν has coefficient1 on the spine and1/2 at every leaf. Its marked
leaf distance is7.

In an actual quadrilateral the opposite black pairs are{a,c},{b,d}.
The selected whole varying affine factor F and the positive definite
remaining supporting forest Z satisfy the inherited
[factor-fiber framework](fresh_cusp_factor_payment_framework_20261017.md).
Every actual edge is strictly positive definite, and each actual corner
is positive semidefinite. Other inert factors are treated by that
framework. The actual full Gram has at most one negative direction.

An adjacent black contact or a common nonzero c/d coordinate in Z gives
a distinct actual length1 or2 unit. Thus we may continue in the branch
where adjacent black products vanish and the nonnegative native
spectator contact vectors x,y have disjoint supports. This is a proved
alternative, not an artificial assumption that contacts vanish.

Delete a or b from F. The remaining connected finite D9 has inverse
entrywise nonnegative, diagonal1 at the old mate, diagonal2 at its fork
center and all other diagonals greater than2. The other marked short
leaf has diagonal9/4. Two nonzero contacts already exceed2 in principal
projection; the source records the exact bound. The only possible PSD
one-root profiles against this entire old D9 are therefore:

|Profile|Old supporting contact|Height h|Moment l|Residual c0|
|---|---|---|---|---|
|isolated|none|t/2|0|2|
|end3|mate, magnitude1|(t+1)/2|1|1|
|end4|mate, magnitude√2|(t+√2)/2|√2|0|
|next3|center, magnitude1|1+t/2|2|0|

Here t is the nonnegative magnitude of the contact to the omitted old
mark. It may be dotted or parallel; no upper bound is imposed. In the
isolated case t>0 because h>0. A nonzero native magnitude is at least1,
and the next possible value after1 is at least√2. At the old mate the
PSD condition permits only1 or√2, and at the center only1. This covers
all native labels, not a finite label search. Added supporting roots
cannot repair a non-PSD principal extension.

The end4/next3 extensions are whole affine B10/D10. Their positive
kernels force every supporting spectator contact to that extension to
vanish; hence the resulting factor is whole. The marked B10 pair has
weight1/2, and the marked D10 opposite-fork pair has weight27/32. Either
is a distinct actual corner. Together with the original27/32 they pay
at least1. Their native t-values do not change this sealing argument.

## Physical Schur identity and the extreme reduction

Eliminate the eight old supporting roots. Scaling a,b by their kernel
coefficients leaves (1/9)[[1,−1],[−1,1]]. The source verifies both end
identities, all16 cross-profile identities and the symbolic determinant.
Writing g=2l/h+c0/h², the rank relation is

    0 = 1 − (g_c+g_d)/9
          + (x/h−y/k)^T G_Z^(-1)(x/h−y/k)/9.                 (1)

The energy is nonnegative. After the immediate affine alternatives,
the only end cost above4 is the extreme profile h=1/2,g=8. It occurs
in exactly two native ways: isolated t=1 or end3 t=0. Thus the extreme
root has one unit old contact, either at the selected old leaf or its
mate. All other admissible isolated/end3 profiles have g≤4, including
arbitrarily large opposite contact t. Two nonextreme sides leave a
strict gap1/9 in(1), so at least one side is extreme.

Assume c is extreme. If d is nonextreme, then

    E = (2x−y/k)^T G_Z^(-1)(2x−y/k) = g_d−1 ≤3.             (2)

We first justify a completely unmixed branch, or obtain a short paid
path. Each of c,d contacts at most one node per connected component of
Z: two contacts would create a cycle in the strictly finite augmented
edge. If a mixed component has c/d contact distance≤2, the distinct
actual c–d path has length≤4 and retains a unit.

Otherwise choose mixed nodes i,j at distance≥3, with native magnitudes
p,q≥1. Their Gram block is2I, so coordinate Cauchy in the **full**
positive definite spectator metric gives

    E ≥ 2p² + q²/(2k²) ≥ 2 + 1/(2k²).                    (3)

For the isolated d profile, k≥1/√2 and g_d=2/k². The lower bound in(3)
minus the required energy in(2) is3−3/(2k²), strictly positive unless
k=1/√2. In that boundary case E=3 and equality forces p=q=1 and

    2x−√2y = G_Z(e_i−e_j/√2).                            (4)

Let u be the first node after i toward j. It is not adjacent to j.
Because c contacts only i in this component, the u-coordinate in(4)
forces a nonzero y_u, in addition to y_j. This violates the finite
one-contact property for d. The argument uses the full equality vector,
so no unexamined spectator can absorb the defect.

For a nonextreme end3 d profile, t≥1, k≥1 and
 g_d=2/k+1/k²≤3. Then E≤2, contradicting(3) outright. Thus the remaining
branch has no mixed component. Its c contribution obeys x^T J_Z x≤3/4.
The complete one-sided theorem in
[energy below one](fresh_spectator_energy_below_one_20261017.md) makes it
either empty or a single unit endpoint contact on a whole A_r,r≤3.
This conclusion concerns the entire x vector, not a chosen subpath.

## Both extremes: complete energy7/4 alternatives

If both sides are extreme, divide(1)'s energy by4:

    (x−y)^T G_Z^(-1)(x−y)=7/4.                            (5)

A mixed contact distance≤3 gives a distinct actual path of length≤5.
Its current weight is at least1/2: the cycle10 and cycle11 distance5
halves are the only relevant reduced short-five cases. Thus27/32 plus
that distinct payment exceeds1.

For a mixed component whose contacts are farther apart, write its
inverse block [[a,c],[c,b]]. Separate strict edge positivity requires
a,b<2, and its contribution is ap²+bq²−2cpq. We apply the all-rank inverse
identities and positive quadratic remainder from the
[five-quarters spectator proof](fresh_spectator_five_quarters_transfer_20261017.md),
using their intermediate identities, not its smaller-energy conclusion.

For A_r put N=r+1 and let d≥4 be the contact distance. At unit contacts

    E_unit=d(N−d)/N ≥ 2d/(d+2).

Consequently d≤14, since30/17>7/4. For each d=4,...,14,
N lies between d+2 and4d²/(4d−7). The source exhausts exactly these
finite ranges and every position satisfying the strict diagonal bound.
This includes the A6 distance4 pairs of energy12/7; it does not silently
restrict all A cases to endpoints.

The remaining classical and exceptional possibilities follow directly
from the same complete inverse formulas: B has no far pair of allowed
nodes; D has only long/short pairs in D6,D7, of energy3/2,7/4; E6 has
its eligible endpoint block of energy4/3; E7 has only one allowable node
and E8 none. F4,H3,H4,I2 have diameter≤3. This covers arbitrary ranks and
native labels.

For every remaining inverse block, the source checks nonnegative
quadratic gradients at(√2,1) and(1,√2); its positive quadratic remainder
then makes those the respective quadrant minima. Both minima exceed7/4.
Therefore every surviving far contact is unit. Every other active
component contributes at least1/2. For one-sided components this is
J_ii≥1/2; for a mixed component metric Cauchy gives

    E_component ≥ (p+q)²/(4−2G_ij) >1/2.

Each listed far unit energy plus1/2 exceeds7/4. A lower unit energy
therefore cannot be supplemented to(5), and all other components must
be uncontacted. The exact remaining possibilities are **whole A15 with
endpoint contacts**, or **whole D7 with long/short-leaf contacts**.

The source checks all four old attachment choices for A15 and both
orientations/all four choices for D7. It checks every actual edge,
all four corners, whole-component roles and exact inertia:

* If A15 attaches through both selected old leaves, the c,d corner is
  whole finite A17 with its endpoints marked. Its restored weightδ
  combines with27/32 to give1.
* In every other A15 attachment choice, a distinct adjacent corner is
  whole D25 with a marked short leaf at distance8, a retained full unit.
* Six D7 configurations fail strict actual-edge positivity or corner
  PSD. The two remaining configurations have a whole affine E8 of
  order9 at the c,d corner, with marked distance7 and positive kernel.
  The whole order≤9 ideal unit is unchanged.

Untouched spectator components are orthogonal direct summands, so the
literal checks here apply to arbitrary supporting forests after the
preceding full energy argument. These are not selected-core embeddings
being mistaken for universal classification.

Finally, if there is no mixed component in(5), one side has total
one-sided energy≤7/8. The one-sided theorem makes it empty or one unit
endpoint contact on a whole A_r,r≤7. There cannot be two such components,
since each costs at least1/2. Select that side for the following corner.

## Whole finite compensation and exact costs

Choose the extreme root just obtained, say c, and the adjacent corner
c,b. If its old unit contact is to the omitted selected leaf a, its
whole factor is A_(r+1) with c an endpoint, while b is a short leaf of
a separate whole D9. The new weightδ pays the remaining deficit exactly.
This includes empty support r=0 and every r≤7.

If c instead contacts the old mate, c,b belong to one whole finite
D_(10+r), with b a short leaf at distance8. This is a retained unit for
r≤7. The source checks all16 literal alternatives and current short8
unit tables for D15,D16,D17. Smaller finite orders retain these units;
the generic finite6/7/8 deletions require both marks deep, which fails
at the marked short leaf. The D25 output above is protected for the same
reason. No long same-component D15–17 weight is being spent.

The A17 endpoint restoration is an increase, so it preserves all earlier
A17 obligations. The source binds the complete current A17 weighted table,
checks the old endpoint weight0, and checks its new total3613/32. It leaves
the finite A14 proof's A15(1,13) unit, finite D14 proof's D15(1,12) half,
and every current D16/D17 alternative untouched.

The actual A1..A8/D9 increments are5/16,5/8,...,5/8. Existing B10/cycle10
coupons are included simultaneously, not substituted for these new roles.
For representative counting, the prior B10 allowances cover ranks1–4;
the additional rank5–8/B9 allowances are5/8 each. The rank6 representative
is E6, not A6. Competing actual types with the same rank are dominated
simultaneously, rather than adding all their increments to one component.
The source checks every ordered actual-type pair and actual/upper base
inequality at scale32. A new complete forest certificate must still use
these changed items and the A17 base restoration.

## Triangles, current priorities and scope

A third root in a triangle containing a changed D10 pair must give PSD
extensions of both deleted D9 diagrams. After deleting a its supporting
contact set is contained in{0,7}; after deleting b it is contained in{5,9}.
These sets are disjoint. Contact to either old marked leaf is also
forbidden by the opposite deleted inverse diagonal9/4. The third root
would be orthogonal to F, contradicting positive height. Thus no actual
triangle contains a newly changed D10 pair.

Give the independent cycle10/C10/B10 arguments priority before this D10
argument. Cycle10 and C10 outputs are unchanged. The B10 heavy-next3
immediate alternative formerly supplied a D10 unit; its weight is now
27/32, still more than the1/2 needed alongside B10's original half.
Thus the B10 proof remains independently valid. The present D10 immediate
B/D alternatives retain enough weight without invoking another selected
witness or the same theorem recursively.

When none of the changed order10 corners is present, inherit the current
face theorem pointwise, with only nonnegative finite additions. When an
older finite6/7/8 or ideal theorem would use a changed D10 short-seven
payment, that actual corner activates the independent order10 proof
instead. The new D10 argument itself uses only independently retained
short paths≤5, protected marked-leaf units, whole affine order≤9 units,
or its explicitly distinct finite coupons/restored A17 endpoint.
The triangle and quadrilateral lower bounds remain2 and1.

The full rank9 envelope now changes. The former low9 concentration cannot
be cited unchanged. A possible separate route is to stop the inherited
linear-cost concentration at rank8, where all factor costs remain r/2,
and to enumerate every terminal partition with arbitrary ranks≥8 and at
most one factor1..7. Neither that proof nor a new forest count is claimed
by this local packet. In particular, no dimension57 theorem follows here.

Final generation86861/fb02de and separate fresh replay56316/3944b6 both
terminated exit0/PASS. The root completed a full source/report
mathematical read and reread the below-one and five-quarters proofs.
The construction reviewer independently read the full source/report,
the complete five-quarters inverse formulas, the below-one proof,
K's PD/inertia helpers and the large-D distance9 leaf rule. Both reviews
checked arbitrary-support profiles and equalities, the complete7/4
spectrum, all twelve literal alternatives, current role retention,
coupon/base accounting, triangles and simultaneous priority. Neither
reviewer executed the producer or replayed spectrum, inertia, coupon,
type-pair or historical census arithmetic. Exact source checks do not formally verify the inherited
actual-star framework or replay the historical all-rank finite/affine
classification. No external refereeing, proof-assistant verification,
new global bound or finite-volume example is claimed.

    python3 scripts/fresh_order10_d_opposite_fork_fractional_transfer_20261018.py
    python3 scripts/fresh_order10_d_opposite_fork_fractional_transfer_20261018.py --verify
