# Half payments on the five antipodal pairs of the affine10-cycle

Assign weight1/2 to each distance-five pair in a whole varying affine
10-cycle. The five antipodal pairs reduce its total cost from45 to85/2.
Every affected quadrilateral retains total at least1, and an affected
ideal corner cannot occur in a triangle.

Add a half at an endpoint of each whole A_r, 1≤r≤4, against the unique
middle node of a separate whole A9. The actual component-pair increments
are1/2,1,1,1 respectively. They fit the representative allowances already
required by the [B10 half theorem](fresh_order10_b_terminal_half_transfer_20261018.md);
there is no further representative increment. All other weights are
unchanged. The resulting order10 family totals are A85/2, B44, C44, D45.
The all-family rank9 envelope remains9/2, so no dimension improvement
or changed low9 concentration result follows from this packet alone.

The [producer](../scripts/fresh_order10_cycle_half_transfer_20261018.py)
and [certificate](fresh_order10_cycle_half_transfer_20261018.json) check
the new physical profiles, height identity, literal whole-corner roles,
triangle domains and component costs. The arbitrary-rank mixed spectrum
used below is the bound theorem in the
[order11 B/C/D transfer](fresh_order11_bcd_half_pair_transfer_20261018.md),
not a new bounded spectator search.

## Complete actual-corner profiles

Use the actual quadrilateral star W+a+b+c+d with opposite pairs
{a,c},{b,d}. Let the selected whole affine factor F be the cycle0,...,9,
with a=0,b=5 and both arcs of length5. Its positive kernel is all ones.
The [factor-fiber framework](fresh_cusp_factor_payment_framework_20261017.md)
gives a positive-definite spectator diagram Z orthogonal to F, positive
opposite heights h,k and strict positive definiteness of W plus each
single black root. Adjacent-black contact or a common spectator contact
already gives a distinct full unit at distance1 or2. Hence assume those
products vanish and the nonnegative vectors x,y on Z have disjoint supports.

Deleting either mark leaves A9. Its only inverse diagonals at most2
occur at the two endpoints and their next nodes. At an endpoint the
allowed single contacts are1 and sqrt2; at a second node only1 is
allowed. All other native magnitudes fail principal semidefiniteness.
The inverse is positive. Two unit contacts are admissible only at the
two endpoints, where the projection is exactly2; every other pair has
projection greater than2. Increasing either endpoint magnitude makes
that equality strict. Three contacts cost at least3·9/10>2. The producer
checks this complete principal domain directly.

Thus a new root has no old contact, one endpoint3/4 or second-node3
contact on either side, or both endpoint3 contacts. The opposite-mark
contact is unrestricted by this principal test. Each of the six
single-contact profiles yields an adjacent marked path of length at
most5. Its opposite old endpoint is a simple leaf with a fixed degree-two
neighbor after arbitrary spectators are restored, since Z is orthogonal
to F. No classical affine component has that leaf configuration;
exceptional affine factors are smaller than the displayed ten roots.
The whole marked component is therefore finite. Its distance-at-most-five
unit remains unchanged. The source checks the six principal finite
Grams, distances and protected leaf/neighbor positions.

Both endpoint3 contacts instead give a whole affine10-cycle at the
adjacent corner. Its positive kernel forces every additional white
orthogonal to the displayed cycle. Its marked distance is again5, so
that different corner contributes a half, as does the original corner.
These two actual weights already sum to1. No recursive cycle theorem
is used for this branch.

## The isolated height cases

Only the two isolated profiles remain. They touch F solely at their
opposite old marks, so h,k are native positive magnitudes:1 or at least
sqrt2, including all higher labels and dotted products. The two equal
arcs give alpha=1/5+1/5=2/5. Direct Schur elimination yields

    0=1−2alpha/h²−2alpha/k²
           +alpha*(x/h−y/k)^T G_Z^(-1)(x/h−y/k).      (1)

The source expands the symbolic determinant to verify this identity.
If both heights are at least sqrt2, the first three terms are at least
1/5, impossible. Up to reflection take h=1. The required normalized
spectator energy is

    E=(x−y/k)^T G_Z^(-1)(x−y/k)=−1/2+2/k².           (2)

If k>sqrt2, then E<1/2. Every nonzero x-coordinate gives E≥x_i²/2≥1/2
by coordinate Cauchy, so the entire x vector vanishes. At k=sqrt2,
E=1/2. Either x is empty or a nonzero coordinate forces x_i=1 and
equality in the full-space bound:

    x−y/sqrt2=(1/2)G_Z e_i.

All other x-coordinates vanish. A neighbor of i then gives an adjacent
c–i–j–d path of length3 and an independent full unit. With no neighbor,
i is a whole A1 and all remaining contacts vanish. Thus these unequal
height cases either pay immediately or select a root whose spectator
support is empty or one unit contact to a whole A1. The other root's
arbitrary contacts do not affect the chosen corner because it is absent
there. Values making the right side of(2) negative are simply impossible.

At k=1 the required unscaled difference energy is3/2. This is the only
case needing the following mixed-component argument.

## Complete spectator energy3/2 and its whole alternatives

Strict positivity of Z+c and Z+d makes each augmented diagram a finite
forest. Each new root therefore contacts at most one node in any one
connected spectator component. If a mixed component has its two contact
nodes at distance at most3, the c–d path has length at most5. At its
actual corner it contributes at least1/2: the only current short-five
halves are the cycle11 orbit and the new cycle10 antipodal orbit.
Together with the original half this suffices.

For distance at least4, apply the complete all-rank inverse calculation
in the energy3/2 section of the
[coupled order11 theorem](fresh_order11_bcd_half_pair_transfer_20261018.md).
Its strict separate-extension hypotheses are exactly the present ones.
It includes every classical rank, exceptional finite type and native
magnitude. The possible whole components with energy at most3/2 are

|Whole component and contacts|Unit-contact energy|
|---|---:|
|A5 endpoints|4/3|
|A6 endpoints|10/7|
|A7 endpoints|3/2|
|D6 long/short leaves|3/2|
|E6 eligible endpoints|4/3|

Every nonunit magnitude in these blocks costs strictly more than3/2.
Every additional active component costs at least1/2. Thus A5,A6,E6
cannot reach exactly3/2 alone or with another component; only whole A7
endpoints or whole D6 long/short leaves remain, with every other
spectator contact zero. The producer rechecks the inherited scalar
block certificate and matches its frozen record. The all-rank
completeness proof remains an explicit inherited input.

Here the new roots attach to the omitted selected old marks, unlike
the four fork-attachment choices in the older theorem. At the(c,d)
corner, a whole A7 spectator with the two new simple ends becomes a
whole A9, marked at endpoints with distance8. A whole D6 spectator
becomes whole E8, marked at simple leaves with distance6; reflection
handles either orientation. Both retain full units in the current rule.
Their supporting complements are orthogonal direct summands. The source
checks all three literal cases, every strict edge, the full Lorentz
signature, the paying finite corner, its whole component and leaf roles.
No larger same-component payment is consumed.

If no component is mixed, write X=x^T G_Z^(-1)x and Y=y^T G_Z^(-1)y.
Then X+Y=3/2. Choose an empty side if one exists; otherwise choose the
side of energy at most3/4. The one-sided part of the
[complete spectrum below one](fresh_spectator_energy_below_one_20261017.md)
forces a single unit contact at an endpoint of a whole A_r, r≤3.
This controls the entire chosen spectator vector, rather than just one
coordinate. All its other components are orthogonal to the selected
paying factor.

## Four exact whole-component coupons

In every remaining case the chosen root, say c, has no spectator
contact or one unit endpoint contact to a whole A_r, r≤3. At the distinct
corner(c,b), the old mark a and opposite root d are absent. The component
through c is consequently whole A_{r+1}, with c an endpoint. The other
marked component is the whole path1,...,9, and b=5 is its unique middle.
No additional supporting root joins these components: F was a whole
affine factor and the chosen spectator vector has just been completely
determined. The four new halves therefore apply literally.

Their exact component-pair increments are1/2 for A1/A9 and1 each for
A2/A9, A3/A9 and A4/A9. The producer checks all four whole finite
corners and the middle-node role. All four old actual costs were zero.
Under the current simultaneous representative map A9 maps to B9.
The B10 theorem has already provided upper allowances1/2,1,1,1 for
these four small ranks, so no additional representative increment is
needed. The complete ordered type-pair domination is checked with both
sets of actual coupons present; competing actual A9 and B9 types are
not incorrectly added as if one component had both types.

Every remaining branch has the original ideal half and a different
actual corner half, totaling1. This completes the arbitrary-support
quadrilateral proof.

## Triangles and simultaneous priorities

For an actual triangle with selected marks0,5, the third root must give
PSD extensions of both deleted A9 paths. Deleting0 permits contacts only
at{1,2,8,9}; deleting5 permits only{3,4,6,7}. These sets are disjoint,
and the old marked nodes are interior forbidden positions in the
opposite deleted paths. Thus the third root is orthogonal to F, contrary
to its strictly positive height. No actual triangle contains a changed
cycle10 corner. Rotation and reflection cover all five antipodal pairs.

Give any actual cycle10 antipodal corner the independent priority above.
It supplies either a retained unit or the original half plus a distinct
weight at least1/2. The C10 and B10 proofs use short paths of length at
most4, their explicit whole finite roles and their retained affine
halves; the new cycle10 changes none of those outputs. In the absence
of all changed actual order10 corners, the prior face proof is inherited
pointwise, with nonnegative added coupons. This also resolves an older
finite6 or ideal argument that had used a now-halved short-five cycle
payment: its actual corner activates the independent cycle priority.

All finite output roles are whole small components or the specified
coupons. They survive every current finite6/7/8 deletion and finite
rank14–17 transfer. The all-family order10 envelope remains9/2 because
D still has45. Hence the old low9 concentration argument can be kept
without alteration. A future reduction of that envelope requires a new
concentration argument; this theorem does not supply one.

Generation `f93d65` and separate fresh replay `fc9aaf` both exited0/PASS.
The root completed a full producer/report mathematical read (`68f39e`),
including the inherited energy3/2 classification/roles (`e64a01`) and
the below-one argument (`da70f3`). The construction reviewer independently
read the complete producer/report, `K.mixed_spectrum`, the all-rank
energy3/2 section, the one-sided below-one proof and current B10 roles.
Both reviews checked the arbitrary-support cases, actual whole-corner
outputs, coupon accounting, triangle domains and simultaneous priorities;
neither reviewer executed the producer or replayed the profile, inertia,
type-pair or historical census arithmetic. Source checks are exact but do not formally verify the
actual-star framework or repeat the historical all-rank classification.
No new finite-forest maximum, dimension bound, external refereeing,
proof-assistant verification or finite-volume example is claimed.

    python3 scripts/fresh_order10_cycle_half_transfer_20261018.py
    python3 scripts/fresh_order10_cycle_half_transfer_20261018.py --verify
