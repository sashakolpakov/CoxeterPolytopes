# Finite distance seven at feature depth six

In the current simultaneous face-payment rule, delete the unit for a
finite marked pair at distance seven whenever both marks are at least
six steps from every endpoint, branch vertex and endpoint of a nonsimple
bond in their whole A, B or D component. Every triangle still receives
at least two and every quadrilateral at least one. No new coupon is
required. This strengthens the [depth-eight theorem](fresh_depth8_finite7_profile_transfer_20261018.md).

The current baseline is the [sparse D16 and fractional order17 rule](fresh_small_d16_coupon_and_ideal17_rebalance_20261018.md).
The [producer](../scripts/fresh_depth6_finite7_adaptive_transfer_20261018.py)
and [exact certificate](fresh_depth6_finite7_adaptive_transfer_20261018.json)
address local geometry and all-rank deleted-pair formulas. A complete
finite-forest count and global dimension conclusion are separate obligations.
Generation, separate fresh replay and two full internal mathematical
composition reviews pass. This local packet is frozen.

## Actual star and complete retained profiles

Use the inherited actual quadrilateral star W+a+b+c+d, with opposite
pairs a,c and b,d, and normalize each root to square two. Each W plus
one black is positive definite; adjacent corners are positive
semidefinite and opposite pairs with W are indefinite. The full Gram
has at most one negative direction. Work on the branch with no adjacent
black pair joined by a white-interior path of length at most six.
Consequently c.b=d.a=c.d=0.

The changed corner supplies a simple path A20, numbered 0 through19,
with a=6 and b=13. Features can occur at its endpoints but not in its
interior. Additional supporting whites remain arbitrary. The retained
W splits as A6+A6+A6. Its inverse entries are nonnegative and its
smallest inverse diagonal is6/7. Three nonzero contacts from c already
give projection at least18/7>2. A contact label at least five has square
greater than5/2, likewise giving projection greater than15/7>2.
Thus c has at most two retained white contacts, each1 or sqrt(2).

The no-short-six assumption leaves positions0,1,2,3,4,5,7,19. Exhaust
all129 contact vectors on these positions. Strict positivity of W+c
and semidefiniteness of the adjacent retained corner (A20 minus a)+c
leave exactly18 profiles:

    empty;
    {19}, {19:sqrt(2)};
    {7}, {7,19}, {7:sqrt(2)};
    {5}, {5,19}, {5,7}, {5:sqrt(2)};
    {4}, {3}, {2}, {1};
    {0}, {0,19}, {0,7}, {0:sqrt(2)}.

Unwritten coefficients are one. The d profiles are reflected by
i→19−i. The opposite magnitudes t=−a.c and u=−b.d are independently
zero or at least one; dotted and parallel values are included.

## The two forced supporting-white alternatives

Let J be the A20 inverse, x the c profile including x_a=t, and
q=x^T J x. If q>2 no enlargement is necessary. Suppose q≤2. Then t=0,
since J_aa=14/3. Every retained contact of c is at distance at least
five from a: the only profiles with q≤2 are empty, {19},
{19:sqrt(2)}, {1}, {0}, {0,19}, {0:sqrt(2)}.

Choose a minimal indefinite principal diagram K in W+a+c. Both a and
c belong to K; deleting either leaves a positive-definite diagram.
K is connected and has order at most ten. The order bound and the
compact Lannér bound of five are explicit published classification
inputs, already used by the preceding theorem. The
[complete classical higher-order catalogue](fresh_high_rank_quasilanner_census_20261017.md)
reproduces all26 minimal diagrams of orders six through ten.

Let B0 be the component through a of K minus c. It is an induced subtree
of the actual selected finite component. Any branch or endpoint of a
nonsimple edge in B0 is also such a feature in that actual component;
it must be at distance at least six from a. Leaves of B0 are deliberately
not treated as features, since the actual component may extend them.
The distance from a to c in K is at least six. Indeed, a contact at
distance at most four in B0 would be inside the retained path, where
the profile has none. Reaching outside that path needs at least seven
white steps. A diagram with at most five nodes is therefore impossible.

Enumerate all768 ordered pairs of distinct finite-deletion nodes (a,c)
in all26 higher-order diagrams, allowing K minus c to be disconnected.
Require the two necessary distance conditions just stated. Exactly two
ordered pairs survive: diagram25, with (a,c)=(8,1) or(8,2). Its edges are

    0-1, 0-2, 2-3, 0-4, 4-5, 5-6, 6-7, 7-8, 8-9,

all simple. These are the two relevant positions in the classical E10
tree. The c contact in B0 is five steps from a. It cannot be on the
b side, since that would give an adjacent c-to-b path of length three.
It is therefore retained position1, and the complete profile is {1}.
The exact labelled maps give one of the following necessary witnesses:

* **Path extension:** an additional supporting white z is joined simply
  to retained position0. Its products with c and positions1 through7
  vanish. As a member of the selected finite tree it cannot contact any
  other retained path node: that would make a cycle. Thus the retained
  positive block extends to A21.
* **Separate white:** an additional supporting white z has c.z=−1 and
  is orthogonal to retained positions0 through7. It is in a different
  component of K minus c. Any contact to an interior position8 through18
  would create a feature within six steps of a or b in the selected
  actual finite component. The only other possible retained contact is
  position19. If present, strict W+c positivity fails: the left A6
  contributes10/7 to the projection of c; the right A6 together with z
  contributes at least7/8. Their sum129/56 exceeds two. Hence z is
  orthogonal to the entire A20 and supplies an independent A1 block.

In the second argument, the native contact to19 is at least one.
The inverse diagonal at z in A6+z is
1/(2−p²·6/7)≥7/8; a nonpositive denominator already contradicts the
positive-definite supporting-white block. Extra contacts cannot restore
positivity. The c contacts in the displayed principals are specified
by the exact E10 map and the retained profile.

In either alternative d.z=0: otherwise there is an adjacent c-to-d
white path of length four or two, respectively. Reflect the entire
argument for d when its original projection r≤2. If both sides need
witnesses, the two extra whites are distinct and orthogonal. Equality
or a nonzero mutual contact gives a short adjacent path of length at
most five whenever a separate white is involved. If both are path
extensions, equality or a mutual contact makes a cycle in the selected
finite component. These arguments allow every unspecified supporting
root; there is no assumption that K minus c is connected or that all
whites lie in the displayed core.

## Complete adaptive Schur certificate

Use the positive block consisting of A20, any forced simple endpoint
extensions, and any forced separate A1 roots. Both c and d now have
projection greater than two. A separate A1 contributes exactly1/2 to
its own projection and zero to the cross projection. Path extensions
use the exact A21 or A22 inverse. All new contacts have just been
forced, rather than assigned for convenience.

For every pair of retained profiles, first discard a pair already
giving a c-to-d white path of length at most six. Split t into0 or1+T
and u into0 or1+U, with T,U≥0. If a projection is at most two, exclude
the branch unless its profile is the forced singleton above; in that
case test both possible witness types. There are766 adaptive branches.

Writing q,r,s for the two new projections and their cross projection,
the Schur block is

    [[2-q, -s], [-s, 2-r]].

In765 branches the polynomial (q−2)(r−2)−s² has nonnegative exact
Q(sqrt(2)) coefficients and strictly positive constant, with minimum5/84. This makes
both eigenvalues negative, contradicting Lorentz signature. The
certificate handles arbitrary nonnegative T,U, not a finite angle grid.
The sole remaining branch is t=u=0 with c joined simply to7 and d
joined simply to12; it needs no witness. Its determinant is−4/7.
This branch is Lorentz compatible and is not discarded by the sign test.

## The remaining branch forces both replacement marks to be leaves

Let G be the complete22-root Gram A20+c+d in that remaining branch.
Its signature is(21,1,0). Suppose an additional actual supporting
white z touches c. No-short-six forces its contacts with b,d and retained
positions9 through17 to vanish. Its contacts with any middle white
7 through12 also vanish: W+c already contains the path c−7−...−12,
and an extra z touching c and that path would create a cycle in the
positive-definite actual edge diagram W+c.

Consequently the magnitude vector v=−(z.root) can be supported only on

    L={0,1,2,3,4,5,6,c},       R={18,19},

with every entry nonnegative and v_c≥1. Put H=G^−1. The block H_RR is
positive definite. Minimize v^T H v over the R coordinates with no sign
restriction; this only weakens the lower bound. The remaining matrix

    H_LL − H_LR (H_RR)^−1 H_RL

is entrywise positive and has c,c entry5/2. Thus the squared norm of
the projection of z onto the core span is at least5/2. The core span
is nondegenerate and contains a timelike direction, so its orthogonal
complement in actual Lorentz space is positive definite. A square-two
root cannot have that projection. Therefore no additional supporting
white touches c. Reflection proves the same for d.

In the actual opposite corner containing c,d, both are whole-component
leaves joined by the simple length-seven path c−7−...−12−d. The finite
classification permits only A8 or D9. The affine classification permits
only affine D9 of order ten or affine E8 of order nine. Finite exceptional
diagrams have smaller diameter; doubled endpoint bonds and affine cycles
cannot contain this whole-leaf simple path. Each named marked pair retains
its unit under the current rule. Additional supporting roots may attach
inside the path, but the classification includes precisely those allowed
completions. No claim that the entire A20 is sealed is needed.

The paying corner is distinct from the deleted selected corner. This
completes the no-short-six branch without any new coupon.

## Current-rule composition and all-rank counts

On the complementary short-six branch use the already independent
short-path routes from the [depth-eleven theorem](fresh_depth11_finite7_zero_transfer_20261018.md)
and their current composition in the [depth-eight theorem](fresh_depth8_finite7_profile_transfer_20261018.md).
Their finite distance-seven/eight outputs have a feature within five
steps or a fixed whole-component marked leaf; their remaining outputs
are retained shorter paths, separate coupons or small ideal factors.
Threshold six preserves all those outputs. Every elementary triangle
distance-seven payment likewise has a feature within five steps, and
the exceptional triangle alternate has a fixed whole-component leaf.
The later finite rank15–17 changes concern distances at least ten;
the changed ideal order17 roles have distance eleven. None is used
recursively as a short-six output. The current sparse-D16 removal keeps
the entire existing output inventory. Triangle total two and quadrilateral
total one are consequently preserved simultaneously.

The numbers of newly deleted distance-seven pairs, relative to paying
all such pairs, are max(0,n−19), max(0,n−20), max(0,n−21) in A_n,B_n,D_n.
Both marks must lie in the unique simple spine interval at feature
depth at least six. Counting its possible starting positions gives
these formulas for every rank. Direct graph enumeration through rank80
checks them in236 finite diagrams. Relative to the old depth-eight
rule, subtract the corresponding old counts n−23,n−24,n−25, each
truncated at zero. All other payments remain unchanged.

## Verification status

Generation session89253, terminal chunk965e1a, and separate fresh replay
session27237, terminal chunk80313a, both pass exit0. Each reconstructs
all129 candidate profiles,768 ordered classical finite-deletion pairs,
the two exact E10 maps, all766 adaptive polynomial branches, the central
inverse/Schur seal and236 graph counts. They do not rerun the historical
26-diagram attachment census or the preceding geometric theorems.

The construction and order16 peers each completed a full mathematical
source/report composition read, including the inherited short-six and
triangle interfaces. Both checked disconnected minimal-core deletion,
genuine new supporting whites, unspecified contacts and their exclusion,
mutual patch orthogonality, unbounded polynomial substitutions, the
Lorentz projection inequality, whole-leaf completion classification and
all-rank counts. Neither executed the package or repeated the catalogue,
profile, Schur, inertia or historical convolution computations. Neither
found a correction. These are internal agent reviews, not external
mathematician refereeing.

Metadata check9cd315 verifies16 direct source bindings, source syntax and
local report links. It is an integrity check, not another mathematical
replay. The classical catalogue, current short-six and triangle geometry,
and finite/affine classification are explicit inherited inputs. No
proof-assistant verification is claimed.

    python3 scripts/fresh_depth6_finite7_adaptive_transfer_20261018.py --verify
