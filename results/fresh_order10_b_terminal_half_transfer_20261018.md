# Two order10 affine B terminal payments can be halved

In a whole varying affine B9 factor of order10, halve the two pairs
joining its doubled terminal to either short fork leaf. The family
cost falls from45 to44. The distinct-corner proof requires two new
ordinary half-coupons: whole A1 versus the doubled terminal of whole
B9, and an endpoint of whole A3 versus that same B9 role. Their actual
component-pair increments are respectively1/2 and1. Existing whole
A2/B9 and A4-endpoint/B9 halves remain unchanged.

This composes with the independent
[C10 endpoint zero theorem](fresh_order10_c_endpoint_zero_transfer_20261018.md)
and the [completed58 rule](fresh_a15_d15_dimension58_20261018.md).
The order10 family totals become A45, B44, C44, D45, so the all-family
rank9 envelope is still9/2. The [source](../scripts/fresh_order10_b_terminal_half_transfer_20261018.py)
and [certificate](fresh_order10_b_terminal_half_transfer_20261018.json)
check the physical identities, whole finite roles and exact ordinary
component costs. No new forest upper or dimension bound is asserted.

## Complete physical profiles and immediate whole factors

Use an actual quadrilateral star W+a+b+c+d with opposite pairs
{a,c},{b,d}. The selected whole affine factor has spine0,...,7, doubled
bond6--7 and fork leaves8,9 at0. Set a=8,b=7. The changed pairs are
(7,8),(7,9), both at distance8. They form an invariant orbit under the
fork interchange. Every other pair retains its current unit.

By the [factor-fiber framework](fresh_cusp_factor_payment_framework_20261017.md),
the other supporting roots Z form a positive-definite diagram orthogonal
to the selected factor. Positive kernel coefficients are1 on0,...,6,
1/sqrt2 at7 and1/2 at8,9. The opposite heights h,k are positive. An
adjacent-black product or shared spectator contact immediately gives
a distinct retained unit at distance1 or2. Hence assume those products
vanish and the nonnegative spectator vectors x,y have disjoint supports.

Deleting a gives a B9 path; deleting b gives D9. In both cases the
only inverse diagonals at most2 are1 at the simple remote endpoint
and2 at the next node. Their inverses are entrywise nonnegative.
Thus two nonzero native old contacts have projection greater than2;
a label at least5 also fails the single-contact endpoint bound.
The complete possibilities are the four familiar profiles below.
The producer checks both full inverse diagonals, inverse signs and
the two-contact exclusion, rather than inferring completeness from a
bounded label catalogue.

The omitted opposite-mark contact has magnitude t, either0 or at least1
with no upper bound. Write g=2l/h+c0/h².

|End|Profile|Height|Moment l|c0|
|---|---|---|---:|---:|
|fork|isolated|t/2|0|2|
|fork|end3|(t+1)/2|1|1|
|fork|end4|(t+sqrt2)/2|sqrt2|0|
|fork|next3|1+t/2|2|0|
|heavy|isolated|t/sqrt2|0|2|
|heavy|end3|t/sqrt2+1|1|1|
|heavy|end4|t/sqrt2+sqrt2|sqrt2|0|
|heavy|next3|t/sqrt2+1|2|0|

Here fork end3/end4 contacts the old mate9; fork next3 contacts the
center0. Heavy contacts occur at6 or5. Positive height requires t≥1
for an isolated profile; the other profiles allow t=0. Polynomial
numerator/denominator signs verify the cost bounds on each whole ray.

End4 and next3 give actual whole affine order10 factors at the adjacent
corner. Their positive kernels force every additional supporting white
orthogonal to the displayed factor, so the roles cannot be spoiled by
spectators. The four possibilities are

|Original end|end4 target|next3 target|
|---|---|---|
|fork|C10 terminal zero|B10 changed half|
|heavy|B10 changed half|D10 opposite-fork unit|

If the C10 target occurs, its independent theorem pays the face; give
that actual-corner case priority. Otherwise every immediate target
contributes at least1/2 at a corner distinct from the original B half.
The two payments total at least1. This is a direct pair of retained
weights, not a recursive invocation of the B transfer. All four native
Grams, positive kernels and marked distances are checked in the source.

## An extreme fork is necessary

Only isolated and end3 profiles now remain. The fork maximum8 occurs
exactly at isolated t=1 or end3 t=0. Both have h=1/2 and exactly one
simple old contact, to a or its mate9. All other fork profiles have
g≤4: a nonextreme isolated magnitude is at least sqrt2; a nonextreme
end3 magnitude is at least1. The heavy cost is at most4, and its end3
cost is at most3.

Eliminating the eight old supporting roots gives resistance R=9.
The fresh physical identities and all16 cross identities yield the same
symbolic determinant formula checked in the C10 packet:

    0=1−(g_c+g_d)/9
          +(x/h−y/k)^T G_Z^(-1)(x/h−y/k)/9.

Two nonextreme ends give a strict lower bound1/9. Therefore the fork
end is extreme. Taking it as c, the exact equation becomes

    E=(2x−y/k)^T G_Z^(-1)(2x−y/k)=g_d−1≤3.       (1)

For any nonzero x-coordinate i, coordinate Cauchy gives E≥2x_i²≥2.
Except for heavy end3 t=0 and heavy isolated t=1, the native cost
bounds give E<2, forcing the entire x vector to vanish. Arbitrary
contacts of d are allowed: that root is absent from the paying corner.

For heavy end3 t=0, k=1 and E=2. A nonzero x-coordinate forces equality,
x_i=1, and

    2x−y=G_Z e_i.

All other coordinates of x vanish. A neighbor of i gives an adjacent
c–i–j–d path of length3, which pays a retained full unit. Otherwise i
is a whole isolated A1 and every remaining contact vanishes. Thus this
case either pays independently or leaves x empty or a unit contact to
a whole A1.

## The new energy-three case has no far mixed component

The last case is heavy isolated t=1, hence k=1/sqrt2 and

    E=(2x−sqrt2*y)^T G_Z^(-1)(2x−sqrt2*y)=3.       (2)

Both Z+c and Z+d are strictly positive-definite principal diagrams of
actual edges. Thus each root contacts at most one vertex in any one
connected spectator component: two contacts would create a cycle in a
finite diagram. If a spectator component has both contacts at distance
at most2, the corresponding c–d path has length at most4 and pays an
unchanged full unit.

Suppose a mixed component instead has contact nodes i,j at distance
at least3, with x_i=p≥1 and y_j=q≥1. They are nonadjacent, so the
two-coordinate Gram-metric Cauchy bound gives

    E≥(2p)²/2+(sqrt2*q)²/2=2p²+q²≥3.

Equality in(2) forces p=q=1 and, in the full spectator space,

    2x−sqrt2*y=G_Z(e_i−e_j/sqrt2).               (3)

Let z be the first neighbor of i along the spectator path to j.
Since the distance is at least3, z is not adjacent to j. The z-coordinate
of(3) is strictly negative, forcing y_z>0. But d already contacts j
in that same connected component, contradicting the preceding finite
tree condition. Thus no far mixed component is possible. This argument
uses the entire spectator metric and handles other arbitrary contacts;
it does not assume an empty complement or a fixed spectator rank.

In the remaining case every contacted component is one-sided. Put
X=x^T G_Z^(-1)x and Y=y^T G_Z^(-1)y. Equation(2) is4X+2Y=3, so X≤3/4.
The one-sided part of the complete
[below-one spectrum theorem](fresh_spectator_energy_below_one_20261017.md)
then says x is empty or consists of one unit contact at an endpoint of
a whole A_r, r≤3. Indeed every active component costs at least1/2,
nonunit contacts cost at least1, and a single non-A inverse diagonal
is at least1. On A_r the endpoint energy is r/(r+1), which is at most
3/4 exactly for r≤3. Every other spectator component is orthogonal
to c. The opposite root d may still have unrestricted contacts on
its own components; it is absent from the chosen paying corner.

## Literal whole payments and exact coupon costs

All remaining branches have x empty or one unit endpoint contact to
whole A_r, r≤3. Consider the distinct actual corner(c,b).

If c's extreme old contact was to the omitted a, then c lies in whole
A_{r+1}, marked at an endpoint, while b lies at the doubled terminal
of a separate whole B9. No other spectator attaches to c or to the
old factor, so these are whole components. The needed halves are

|Small whole component|Existing or new|Actual component-pair total|
|---|---|---:|
|A1|new unique-vertex half|1/2|
|A2|existing two endpoint halves|1|
|A3|new two endpoint halves|1|
|A4|existing two endpoint halves|1|

The existing roles are bound to the
[order12 B theorem](fresh_order12_b_center_terminal_transfer_20261018.md)
and [order13 sparse B theorem](fresh_order13_bc_sparse_terminal_transfer_20261018.md).
If c instead contacted its old mate9, the selected component is a whole
finite B_{10+r} with marked distance9 and b its doubled terminal leaf.
For r≤3 every such role retains a full unit in the completed rule.
The source reconstructs all eight whole finite corners and their marks.
Every coupon branch contributes1/2 in addition to the original ideal
half, and every joined branch contributes a full unit.

The two new actual increments are1/2 for A1/B9 and1 for A3/B9.
Their corresponding representative allowances must increase by those
same amounts; no hidden existing slack is claimed. The source checks
simultaneous domination for every inherited ordered actual-type pair.
It does not run a finite-forest optimization. A separate conditional
cost preflight also included prospective cycle10 coupons and left the
old ambient58 upper384.5 unchanged, but that private diagnostic is not
an authoritative combined certificate or a geometric cycle theorem.

## Current priority, triangles and scope

First handle any actual C10 endpoint corner by its independent zero
theorem. In its absence, an actual newly halved B10 corner has the
retained half plus a distinct retained half/unit just proved, or an
independent short unit. If neither changed kind occurs, the old actual
weights agree pointwise with the completed58 rule, with only nonnegative
coupon additions. Its face theorem applies unchanged. This covers older
uses of the B10 unit without assuming their selected witness survives.

All short outputs have length at most4. The only same-component finite
outputs are the explicit B10 through B13 distance9 roles. Current finite
rank14–17 transfers and all generic deep finite6/7/8 zeros miss them.
No current D16, A17 or D17 same-component obligation is changed. The
changed ideal distance8 pairs are outside the triangle theorem's ideal
distance-at-most7 and whole-affine-order-at-most9 outputs. Triangles
therefore retain total2 and quadrilaterals total1.

Because A and D remain45, the rank9 envelope and its existing low9
concentration argument are unchanged. A future all-family improvement
needs its own concentration justification. The completed unrestricted
bound remains58 by its frozen proof; this packet supplies no bound57,
no exact new forest maximum and no finite-volume example.

Final source-bound generationea64de and separate fresh replaydc33bb
both terminated exit0/PASS. The parent and a peer completed full
mathematical producer/report reads, including the final B10--B13 cost
guard. The parent reread the C10 determinant/profile source and the
complete below-one spectrum report. The peer also reread the inherited
Gram/graph/distance/positive-pivot helpers, the complete one-sided
spectator argument and the existing A2/B9 and A4/B9 role provenance.
Both passed the full native inverse/profile domain, immediate whole
affine alternatives, extreme-fork and E2/E3 equality arguments, whole
finite roles, exact coupon/representative increments and current
short-output/triangle priority. Neither reviewer executed any inverse,
scalar, corner, type-pair, forest or historical census calculation.
Direct metadata audit passes24 source bindings,8 report links and
source syntax. Source, certificate and report are final/frozen.

Exact scalar and literal checks do not
constitute independent proofs of the inherited actual-star framework
or the all-rank below-one spectrum. No external refereeing or proof
assistant verification is claimed.

    python3 scripts/fresh_order10_b_terminal_half_transfer_20261018.py
    python3 scripts/fresh_order10_b_terminal_half_transfer_20261018.py --verify
