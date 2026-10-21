# Finite distance six at feature depth six

Delete the finite distance-six unit whenever both marked nodes have graph
distance at least six from every endpoint, branch vertex and endpoint of a
nonsimple bond in their whole classical component. Add a unit between a
whole A1 and either middle node of a whole A12. Relative to the
[current combined rule](fresh_order15_d16_d17_ambient59_frontier_20261018.md),
every triangle still receives at least two and every quadrilateral at least
one. All ideal weights and all representative cross-component allowances
remain unchanged. This is a local theorem; no new dimension bound or
complete finite-forest maximum is asserted here.

The [producer](../scripts/fresh_depth6_finite6_zero_transfer_20261018.py)
and [certificate](fresh_depth6_finite6_zero_transfer_20261018.json) promote
the [complete adaptive diagnostic](fresh_depth6_finite6_adaptive_pilot_20261018.json)
by resolving its universal supporting-white obligations and all seven
remaining branches. The diagnostic itself remains frozen with its more
limited scope. In particular, its unsuccessful single-mark projection
test is not treated as a seal.

## Actual star and the supporting-white restriction

Use the actual quadrilateral star W+a+b+c+d, with opposite pairs a,c and
b,d and all roots of square two. W plus any one black node is positive
definite; actual corners are positive semidefinite; W+a+c and W+b+d are
indefinite. The full Gram has at most one negative direction. These are
the inherited actual-star hypotheses, with arbitrary additional supporting
whites allowed.

The changed finite corner contains a simple A19 path, numbered 0 through18,
with a=6 and b=12. Work first on the branch where no adjacent black pair
has a white-interior path of length at most five. Then c.b=d.a=c.d=0.
An additional supporting white can contact the retained path only at 0 or18:
every interior position already has two simple path neighbors and is within
five steps of a or b, so another contact would create a forbidden branch
feature. No additional white contacts both endpoints, since that would
create a cycle in the selected finite forest. These restrictions hold
for all supporting whites, regardless of their labels or the ambient rank.

The retained W is A6+A5+A6. Its inverse entries are nonnegative and the
smallest diagonal is5/6. Three nonzero contacts of c give projection at
least5/2>2. A label at least five has squared magnitude greater than5/2,
giving projection greater than25/12>2. Thus at most two retained white
contacts occur, each of magnitude one or sqrt(2).

The no-short-five condition permits positions0,1,2,3,4,5,7,17,18. All163
contact vectors with at most two nonzero entries are checked. Strict W+c
positivity and semidefiniteness of (A19 minus a)+c leave exactly19 profiles,
in certificate order:

    0 empty;
    1 {18}; 2 {18:sqrt(2)}; 3 {17};
    4 {7}; 5 {7,18}; 6 {7:sqrt(2)};
    7 {5}; 8 {5,18}; 9 {5,7}; 10 {5:sqrt(2)};
    11 {4}; 12 {3}; 13 {2}; 14 {1};
    15 {0}; 16 {0,18}; 17 {0,7}; 18 {0:sqrt(2)}.

Unwritten magnitudes are one. Reflect i to18-i for d. The opposite
products t=-a.c and u=-b.d are independently zero or at least one;
the calculation includes arbitrarily large dotted products.

## Low projections force one of two actual extra whites

Let J be the A19 inverse and q the projection of c onto that path,
including its a-contact t. If q is at most two, t is zero. The only
such profiles have indices0,1,2,3,14,15,16,18. In each, c has no retained
contact within four white steps of a.

Take a minimal indefinite principal diagram K in W+a+c. It contains
both a and c; deleting either is positive definite. K is connected,
has order at most ten, and its component through a after deleting c
inherits every branch or nonsimple-edge feature from the selected finite
component. Such a feature is at distance at least six from a. Also the
distance from a to c in K is at least six: a shorter contact would lie
in the retained path where the profile excludes it, and reaching an
additional white through an endpoint already takes at least seven steps.
Thus K cannot have order at most five.

The inherited [complete higher-order catalogue](fresh_high_rank_quasilanner_census_20261017.md)
has26 diagrams and768 ordered distinct finite-deletion pairs. The exact
screen allows disconnected K minus c. Under the two distance conditions
only diagram25, with (a,c)=(8,1) or(8,2), survives. It is the E10 tree
with simple edges

    0-1, 0-2, 2-3, 0-4, 4-5, 5-6, 6-7, 7-8, 8-9.

Its c-contact is five steps from a. It cannot point toward b, since
that would give a short adjacent c-b path. Therefore the complete
retained profile is {1}. There is an additional actual supporting white z
of one of two types:

* z extends the path simply at0 and is orthogonal to c and the rest
  of the retained path;
* z is joined simply to c and is orthogonal to the retained path.

The explicit E10 maps in the [preceding adaptive theorem](fresh_depth6_finite7_adaptive_transfer_20261018.md)
relabel to the following maps from its nodes0 through9, with c=19 and z=21:

    extension: 1,19,0,21,2,3,4,5,6,7;
    separate:  1,0,19,21,2,3,4,5,6,7.

For the extension, a second retained contact of z would make a selected
finite cycle. For the separate white, the map fixes orthogonality to0
through7 and the feature restriction excludes every other interior
position. A possible contact at18 is excluded by strict W+c positivity:
the left A6 contributes10/7 to the c-projection, and the far A6 together
with z contributes at least7/8. Their sum129/56 exceeds two. For a native
far-end magnitude p at least one, the latter diagonal is
1/(2-p^2*6/7), at least7/8 whenever the denominator is positive; a
nonpositive denominator already contradicts the white-block positivity.

In both cases d.z=0, since otherwise c and d have an adjacent path of
length at most four. If reflected witnesses are needed on both sides,
they are distinct and orthogonal. Equality or a mutual contact gives a
c-d path of length at most five whenever a separate white occurs; two
extensions would instead make a cycle in the selected finite forest.
Thus the positive block is the retained path with zero, one or two
simple extensions and zero, one or two separate A1 roots. No unspecified
supporting contact is assigned zero without these arguments.

## Unbounded Schur branches and native gap

After the forced enlargement, both black projections q,r exceed two.
Write s for their cross projection. Split each opposite magnitude into
zero or1+T, respectively1+U, with T,U nonnegative. The Schur block is

    [[2-q,-s],[-s,2-r]].

The complete829 adaptive branches have822 strictly positive determinant
polynomials with nonnegative exact Q(sqrt(2)) coefficients. Both
eigenvalues are negative, contradicting the Lorentz index. Branches
with a low projection other than the forced singleton have already
been excluded by the complete minimal-core argument. A profile pair
with a short c-d path is handled independently below.

Exactly seven branches remain, given as (c-profile,d-profile,t,u):

    (0,4,1,0), (0,5,1,0), (4,0,0,1),
    (4,4,0,0), (4,5,0,0), (5,0,0,1), (5,4,0,0).

In a remaining branch with a nonzero opposite product, all nonconstant
determinant coefficients are nonnegative and substitution sqrt(2) for
that product gives a strictly positive determinant. The native gap
between one and sqrt(2) therefore forces that product to be exactly one.
Up to reflection the seven branches give four exact21-root Grams.
Each has inertia(20,1,0), and every retained W plus one black is
positive definite. The following arguments address arbitrary additional
supporting whites in all four cases.

## Four residual configurations

Write c=19,d=20. All unspecified core contacts vanish.

**Mixed case:** c joins6 and d joins11, both simply. An additional white
touching d can contact only d,0,18: a c-contact gives a short path, and
interior contacts violate the selected feature restriction. With H=G^-1,
the H-block on0,18 is positive definite. Minimizing the projection
quadratic form over both endpoint coordinates without sign restrictions
gives at least(9/4)*v_d^2, hence greater than two. Since the core already
contains a timelike direction, its orthogonal complement is positive
definite; a square-two root cannot have this projection. Thus d is a
whole-component leaf. Its neighbor11 has fixed degree two at the actual
corner (a,d), whose marked distance is six. A finite classical
completion pays a retained unit because a mark is a whole leaf. Finite
exceptional completions have rank at most eight and also retain their
basic distance-six units. An affine
classical leaf adjacent to a degree-two node along a simple edge is
impossible; the exceptional affine alternatives have order at most nine
and retain their units. This pays at a distinct actual corner.

**Mixed cycle:** c joins6; d joins11 and0. The corner (a,d) contains the
affine13-cycle on0 through11 and d. Its strictly positive kernel vector
is all ones. In a positive-semidefinite actual corner, every additional
white is orthogonal to this whole cycle: its nonpositive contacts have
zero sum against the kernel. A white touching c can therefore contact
only c and18 in the core. Free minimization over18 gives projection at
least(78/7)*v_c^2>2. Thus c is a whole A1 at corner (c,d). The component
through d is the whole A12 path

    5,4,3,2,1,0,d,11,10,9,8,7.

It is sealed by the affine kernel and the c-seal. Here d is middle
node6, or5 in reverse order. The new whole A1/A12 middle-node unit
pays this distinct corner. This is the only new coupon.

**Central case:** c joins7 and d joins11. It suffices to prove that at
least one of c,d is a whole leaf. Suppose an additional white z touches
c. It cannot touch d, and its only possible retained contacts are0,18,
at most one of them. In the positive-definite W+c principal, the three
relevant blocks are A6: the left path0 through5, the middle path c,7
through11, and the right path13 through18. Strict positivity and native
magnitudes leave exactly

    {c:1}, {c:sqrt(2)}, {c:1,0:1}, {c:1,18:1}.

Indeed a label at least five already exceeds the single-contact budget;
with two contacts, any sqrt(2) exceeds the budget18/7. Exact projection
onto G excludes {c:sqrt(2)} and {c,0}, whose projections are7/3 and29/12.
The survivors are {c} and {c,18}, with projections7/6 and17/12.
Reflection leaves {d} and {d,0} for a white w touching d.

The whites z,w must be distinct and orthogonal, since otherwise they
give a c-d path of length two or three. The two-white Schur determinants
for the four surviving contact pairs are, in order,

    -13/12, -7/8, -7/8, -5/6.

Each gives an additional negative direction beyond the one in G.
Therefore both marks cannot have additional neighbors. At least one
is a whole simple leaf with fixed degree-two neighbor7 or11. The actual
corner (c,d) has marked distance six and pays a retained unit by the
same finite/affine classification as the mixed case. The diagnostic's
failed one-mark lower bound13/18 is superseded by this two-neighbor
argument, not silently strengthened.

**Central cycle:** c joins7; d joins11 and0. The same affine13-cycle
seals all additional whites against0 through11 and d. A possible white
touching c contacts only c and18; its minimized projection is at least
(78/31)*v_c^2>2. Thus c is sealed. Corner (c,d) contains the whole A13
path obtained by appending c to the mixed-cycle A12 path. Its marked
distance is six and c is a whole leaf, so its unit remains. Reflection
handles the other three residual branches.

## Independent short paths and triangles

Every current finite marked path of length at most five retains its
payment. Every current ideal short-five unit is also unchanged. The
exception is an affine11-cycle pair at distance five, whose half is
handled by its independent [replacement theorem](fresh_cycle11_half_pair_transfer_20261018.md).
That proof uses shorter retained paths, a distinct cycle11 half, or
finite paths of length six. Its only two length-six single profiles,
end3_1 and end4_1, have a fixed whole leaf1 at distance four from the
old marked node5. They therefore cannot trigger the present deletion.
Consequently the short-five branch still gives a full quadrilateral
payment without invoking the present theorem recursively.

The producer reconstructs all1,025 high-core triangle triples. Every
one of their54 elementary distance-six payments has a persistent
branch or nonsimple-edge feature within five steps. The three
exceptional alternate paths of length six have a protected whole leaf
within two, one and zero steps of a mark, respectively. The fourth
exceptional alternate has length seven and is unchanged. Lower-order
cores retain the two paths of length at most five supplied by the
inherited all-label Steiner argument. Every triangle still receives two.

For simultaneous composition, inspect actual corners. If a newly changed
finite distance-six corner occurs, use this universal proof directly.
If none occurs, all actual corner weights equal the predecessor weights
plus the nonnegative new coupon, so the completed predecessor whole-face
payment applies. This priority avoids requiring every historical
preferred witness to survive in a different hypothetical face. The
independent short-five interface above is the only replacement theorem
invoked within the new proof.

## Actual costs and verification

For A_n,B_n,D_n the new deleted-pair counts are respectively
max(0,n-18), max(0,n-19), max(0,n-20). The admissible marks form a
simple spine interval; counting distance-six starting positions gives
these formulas for every rank. Direct graph enumeration checks236
diagrams through rank80. These are additional distance-six deletions,
not replacements for the preceding distance-seven deletion counts.

The whole A1/A12 coupon adds exactly two to that actual component-pair
cost, previously zero. The unchanged representative A1/D12 allowance is
5/2. All32,400 ordered actual type pairs with each type's rank at most59
satisfy simultaneous domination, including every inherited coupon;
every actual same-component cost is dominated as well. Actual A12 and
D12 are different types and their costs are not added to one another.
All216 current literal replacement configurations are reconstructed;
none of their corner payments changes and every total remains at least
one. This finite check supplements the universal geometry above.

Generation session22306/chunk9f27a2 passed with
`PASS_DEPTH6_FINITE6_ZERO_LOCAL_TRANSFER`. Separate fresh replay
session29495/chunkea8eac also passed with that status. The construction
peer completed a full producer, adaptive-diagnostic and report mathematical
composition read: PASS. Its scope covered the complete contact domain,
disconnected minimal-core screen and relabelled E10 witnesses, mutual and
opposite witness orthogonality, unbounded/native branch reduction, arbitrary
supporting-white seals, central two-neighbor argument, independent short-five
routes, triangle protections, actual-corner priority and full actual versus
representative cost accounting. That review did not execute the profile,
Schur, inertia, literal, type-pair or catalogue calculations. The root also
read the inherited short-five cycle11 and triangle proofs in full.

The source, data and report are final and frozen. These are internal
mathematical reviews, not external refereeing or formal verification.
The producer freshly reconstructs the adaptive profiles, catalogue
screen and Schur branches, the four exact residual arguments, triangle
interfaces, graph counts, literal configurations and actual type-pair
domination. It does not rerun the historical catalogue census, prove
the written arbitrary-support reasoning by a proof assistant, or
perform a finite-forest dynamic program. The explicit relabelled E10
maps above are inherited from the preceding exact map check.

    python3 scripts/fresh_depth6_finite6_zero_transfer_20261018.py --verify
