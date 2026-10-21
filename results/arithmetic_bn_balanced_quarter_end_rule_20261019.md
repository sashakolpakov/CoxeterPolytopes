# A balanced quarter-unit rule at the simple end of B19/B20

19 October 2026. **New written consumer theorem, awaiting independent
review.** No scientific job, proof script, optimizer, Lean process or build
was run. Existing pinned reports and charge sources are unchanged.

The new rule preserves all face payments in the eight high-label fields.
Combined with the separately reviewed simple-end distance-two deletion,
it brings the full B20 vertex cost to 35.25. The same cost occurs at
B19+A1 after its compensating quarter is added. This does not prove
allowance 35 or lower a dimension bound; the A20 and D20 obstructions
also remain.

## 1. Exact field and whole-forest scope

Let P be a compact arithmetic simplest-type Coxeter n-polytope,
n=19 or20, with ground field in

    Q(sqrt(3)), Q(sqrt(6)), Q(sqrt(2),sqrt(3)),
    K9,K11,K15,K16,K20,       Km=Q(cos(2pi/m)).

Use the inherited charge from
`results/fresh_arithmetic_high_label_20261016.md`: one unit for
marked distance at most two, and a further unit for an entire marked
I2(m) component with m>=6. Every triangle has an unchanged entire
high-label I2 corner paying two, and every quadrilateral originally
receives at least one. These are the previously proved eight-field
inputs, not new assumptions asserting this note's conclusion.

Make the following three modifications. Each context is the ENTIRE
rank-n finite vertex forest; there are no unlisted spectators.

1. At a whole B_n vertex, replace the simple-end ADJACENT baseline
   unit by 1/4. Its simple end is the end opposite the terminal
   4-bond.
2. At a whole B_(n-1)+A1 vertex, add 1/4 to the mark consisting of
   the isolated root and the simple endpoint of the B_(n-1) component.
3. At a whole B_(n-2)+A1+A1 vertex, add 1/4 to the pair of isolated
   roots.

There is exactly one marked role of each kind at a forest of its
specified type. The two new recipients have original baseline zero.
No unit is deleted at the I2(m) components, and no other weight is
changed by these three modifications.

**Balanced-rule theorem.** These modifications preserve payment at
least two on every triangle and at least one on every quadrilateral.
All weights are nonnegative. They can also be combined with the
whole-B_n simple-end distance-two deletion and the whole-D_n
whisker-pair deletion specified in Section 5 below.

The theorem concerns a fractional charge, not a new Coxeter angle or
a count of physical polygons. Its proof does not require that an
arbitrary local Gram pattern extend to a compact global chamber.

## 2. The affected quad has a sealed large-B white support

Only the first modification decreases a weight. Consider any quad
containing this whole-B_n simple-end adjacent role. Call that pair
a,b, with a the simple endpoint. Its fixed white support is the
remaining connected B_(n-2) path W, including the remote terminal
4-bond. At every corner of the quad the ENTIRE finite diagram is
W plus its two adjacent black roots: a compact Coxeter polytope is
simple, and the quad has codimension n-2. Thus no additional outside
factor or vertex wall is being suppressed in what follows.

Every individual black root occurring on this quad has a finite
one-root continuation of W, because it occurs in a positive definite
corner. Since W has rank at least 17 and already contains a 4-bond,
the connected finite component containing W can only be a B path.
The exceptional finite diagrams are too small, and A,D,E diagrams
have no 4-bond. The old terminal 4-bond must remain terminal.

Consequently each black root has exactly one of two profiles:

    I: no contact with W;
    N: one simple contact at W's simple endpoint, no other W contact.

This is a structural finite-diagram classification, not a contact
restriction imposed in order to make the calculation close. A
contact at the short endpoint would move the old 4-bond into the
interior, another 4-bond would give two double bonds, and any other
attachment would branch or cycle. None is a finite large B path.

## 3. Every finite corner of this quad pays at least one quarter

Examine the two black profiles at any corner. The following list is
complete, because W plus those two roots is its entire finite diagram.

### I,N

If the two black roots have a nonright product, their component is
the whole connected B_n path. Its remote 4-bond is already in W,
so the new black-black bond must be simple. The marked pair is
precisely the simple-end adjacent role in modification 1, paying 1/4.

If their product is zero, the entire forest is B_(n-1)+A1, with
the N root as the simple endpoint of the connected component and
the I root isolated. It is precisely the recipient in modification 2,
and pays 1/4.

No other finite label is possible in the connected case: a second
4-bond, a label at least five, or any exceptional connected type is
incompatible with this large finite B path.

### I,I

If the black product is zero, the forest is B_(n-2)+2A1. The marked
pair is precisely the recipient in modification 3, paying 1/4.

If it is nonzero, the entire forest is B_(n-2)+I2(m) with m>=3.
The marked pair retains its original baseline one (or two including
the high-label addition). It is not any reduced whole-B_n role.
Thus this corner pays at least one.

### N,N is impossible

Both black roots would attach at W's same simple endpoint. If their
mutual product is zero, that endpoint has three neighbors, and the
connected component containing the old 4-bond branches. If the
product is nonzero there is already a triangle. Neither is a finite
large B path. No finite exceptional component can contain W's rank
at least 17. Hence this profile cannot be an actual finite corner.

Every corner of the quad therefore pays at least 1/4, with all
simultaneous occurrences of the reduced role included. Its total
is at least one. No assumption about how many corners were originally
paid, which opposite sides are dotted, or which cyclic attachment
pattern occurs is required.

This explains why adding only an isolated-pair coupon would be
insufficient: a quad can alternate I,N around all four sides and
have no I,I corner. The cross recipient in modification 2 covers
the right I,N corners in exactly that case.

## 4. Other faces and unchanged lower bounds

Triangles retain their entire high-label I2 corner paying two. Such
a marked component cannot be a whole connected B_n vertex, so the
only reduction in modification 1 does not change that payment.
The new cross roles only increase weights elsewhere.

A quad not containing the reduced adjacent role has no decreased
weight from these three modifications, so its old lower bound survives.
Section 3 covers every quad containing the reduced role. Polygons
with k>=5 require only a nonpositive amount 5-k, and all new weights
are nonnegative. Thus these modifications form a complete consumer
argument within the stated scope.

In particular positive new cross roles do not invalidate any old
face lower bound: they only add charge. Their cost at a vertex must,
however, be retained in the global upper bound; it is not an uncharged
transfer.

## 5. Composition with the two earlier deletions

The source
`results/arithmetic_bn_simple_end_distance_two_refund_20261019.md`
proves deletion of the unit at the first and third roots of a whole
B_n simple end. Its two-face has white graph

    A1+B_(n-3),

which is disconnected. The adjacent role reduced here has connected
white graph B_(n-2). These two roles cannot occur at different corners
of the same two-face: the actual white roots of that face are fixed.
Therefore no quad treated in Section 3 loses any additional charge
from the distance-two deletion, and none of that theorem's affected
quads contains the reduced adjacent role. Nonnegative coupons can
only help either proof.

The whole-D_n whisker-pair deletion is covered by
`results/arithmetic_dn_whisker_pair_triangle_retention_20261019.md`
and `results/arithmetic_dn_whisker_pair_quad_retention_20261019.md`.
A face carrying either of the whole-B_n reduced roles has a white
4-bond. Every corner contains that same bond and therefore cannot
be a whole D_n, which is simply laced. No quad can mix the D_n
deletion with either B_n reduction.

All triangles remain paid at their untouched high-label I2 corner
after every listed modification. Larger polygons retain nonnegative
weights. Thus the composition follows without asserting arbitrary
spectator extensions of any local theorem.

## 6. Exact vertex costs and the remaining gap

The ordinary baseline at a large B_r component is 2r-3. Its entire
I2 high-label additions do not occur at any of the following recipient
or source forests. At dimension n=20, the changed costs are therefore

| Whole forest | Baseline | Net change | Final cost |
|---|---:|---:|---:|
|B20|37|-1 distance-two, -3/4 adjacent|141/4 = 35.25|
|B19+A1|35|+1/4 cross coupon|141/4 = 35.25|
|B18+2A1|33|+1/4 isolated-pair coupon|133/4 = 33.25|

Each coupon is once at the specified vertex, since its marked
recipient pair is unique. These are actual per-vertex costs, not
costs per component occurrence that can be omitted on a forest.

The full A20 cost is still 37. The full D20 cost is still 37 after
its separate one-unit whisker deletion; large spectator cases are
not newly covered. Hence no rank-20 forest allowance 35 follows.
Even within the B branch one more quarter must be removed or funded
at BOTH B20 and B19+A1 without losing these face payments. This is
a precise remaining interface, not a claim that the fractional gap
may be ignored in the strict global averaging inequality.

The earlier exact Q(sqrt(3)) survivor with four full B_n corners
shows that 1/4 is the smallest uniform adjacent weight compatible
with that local model if no other payment is added at its same four
corners. The present consumer proof meets that constraint and also
covers the mixed corner patterns absent from that model.

Independent review of this new written proof is pending. No code
implementing a new global charge, enumeration certificate, manuscript
claim, formal theorem or dimension bound was changed.
