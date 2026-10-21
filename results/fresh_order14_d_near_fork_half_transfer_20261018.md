# Four universal affine D13 near-fork halves

For an actual ideal factor of type affine D13, of order14, the four
near-fork/far-leaf pairs may be reduced from1 to1/2. In the labeling with
spine0--1--...--9, leaves10,11 at0 and leaves12,13 at9, they are
(1,12),(1,13),(8,10),(8,11). They form one intrinsic distance9 orbit.
The total falls from82 to80. The fork-center pair(0,9) was already zero
and is unchanged.

There is exactly one new finite coupon: weight1/2 between a whole A1 and
the unique long-arm neighbor of the fork center of a separate whole D13.
For finite D13 with spine0--...--10 and short leaves11,12, the large mark
is9. There is one eligible marked pair, so the increment is exactly1/2
per component pair, raising its current actual total from3/2 to2. This
coupon is distinct from the D13 fork-center and short-leaf roles.

Together with the frozen B14 and C14 reductions, the ideal order14 family
totals are A70, B80, C80, D80, with envelope80/14. The theorem covers
arbitrary native labels, arbitrary supporting spectators and every actual
occurrence of the changed ideal orbit. A complete finite-forest count and
any new ambient bound are separate obligations. The
[producer](../scripts/fresh_order14_d_near_fork_half_transfer_20261018.py)
and [certificate](fresh_order14_d_near_fork_half_transfer_20261018.json)
check the exact physical Gram identities and marked component roles.

## The actual star and complete principal profiles

Write the quadrilateral star as W+a+b+c+d, with opposite pairs{a,c},{b,d}.
The actual edge diagrams W plus one black root are strictly positive
definite; adjacent corners are positive semidefinite. Roots have square2.
At W+a+b, the inherited factor-fiber theorem supplies the whole varying
affine factor F and a positive-definite spectator diagram Z orthogonal to
F. This applies without a simple-cusp assumption. Use a=1 and b=12;
automorphisms give the other three changed pairs. The positive kernel of F
has coefficient1 on the spine and1/2 at each leaf.

If b.c, a.d or c.d is negative, a distinct actual corner has a retained
distance1 unit. If c and d share a nonzero spectator contact, a retained
distance2 unit pays. Thus assume those black products vanish and the
nonnegative spectator contact vectors x,y have disjoint supports. The
opposite products a.c and b.d remain unrestricted. Every nonzero native
magnitude is at least1; a magnitude greater than1 is at least sqrt2,
including the parallel and dotted possibilities. Kernel heights are
positive.

Deleting a from F leaves A3 on10--0--11 and finite D10 on2,...,9,12,13.
Deleting b leaves finite D13, with long endpoint13. In either large finite
block, inverse diagonals at its first two long-arm vertices are1,2, and
all other diagonals exceed2. A single nonzero contact in a positive-
semidefinite one-root extension can therefore only be simple or doubled
at the first vertex, or simple at the second. Two contacts would create a
cycle in a connected diagram still containing the old fork; finite and
affine Coxeter classification excludes this, and exceptional diagrams
are too small to contain D10. This proves the complete profile list for
all native labels, rather than assuming a finite label cutoff.

The doubled-first and simple-second profiles create an affine component
with positive kernel. In the actual positive-semidefinite corner that
kernel forces every other corner root to be orthogonal to that whole
component. In particular the opposite new root's old A3 and spectator
contacts vanish when such a near profile occurs. Its product with the
omitted opposite black root remains unrestricted.

Put t=-a.c, and u=-10.c, v=-0.c, w=-11.c. Define

    K = [[3/4,1/2,1/4],[1/2,1,1/2],[1/4,1/2,3/4]],
    e=(u,v,w) K (u,v,w)^T,     h0=t+(u+w)/2+v.

Eliminating F minus a,b gives resistance R=10. The near profiles have
height h, residual constant c0 and moment l as follows:

|Near profile|h|c0|l|
|---|---|---|---|
|Only outer A3/opposite contacts|h0|2-e|0|
|Simple contact at2|h0+1|1-e|1|
|Doubled contact at2|t+sqrt2|0|sqrt2|
|Simple contact at3|t+1|0|2|

For the far root put t=-b.d. Its complete profiles are:

|Far profile|k|d0|m|
|---|---|---|---|
|No contact to F minus b|t/2|2|0|
|Simple contact at13|(t+1)/2|1|1|
|Doubled contact at13|(t+sqrt2)/2|0|sqrt2|
|Simple contact at9|t/2+1|0|2|

The residual new-root squares before eliminating Z are c0+l²/10 and
d0+m²/10. Every cross projection is lm/10. These statements, the two
complete inverse-diagonal lists and all16 cross identities are checked
directly in the norm-two Gram, not after an unrecorded gauge change.

The last two **near** profiles already have independent alternatives.
A doubled contact at2 gives a distinct whole affine B10 factor of order11,
with doubled-terminal/fork-leaf marks. Its weight is currently zero, so
apply the independent BC11-zero theorem at that actual corner first.
It supplies a full retained payment. A simple contact at3 gives a whole
affine D10 factor of order11 with opposite fork-leaf marks, of current
weight1/2. That distinct half plus the selected D14 half is1. The source
checks the positive kernels, native canonical embeddings and current
weights of both outputs. The far doubled/second profiles need not be
spent as alternative payments; the inequalities below exclude them.

## Complete native cost bounds

For an end define g=2l/h+c0/h², with the corresponding far variables.
After the two immediate near alternatives have been handled, every near
profile has g<=5, with a unique extreme type: a simple contact to one of
the old outer leaves10,11, no opposite contact, height1/2 and g=5. Every
other near profile has g<=3. Here is a complete argument.

For the outer-only profile, A3+c is a principal subdiagram of the strict
actual edge. Thus e<2. Any two nonzero contacts give
e>=Kii+Kjj+2Kij>=2, so at most one A3 contact is possible. With none,
positive height implies t>=1 and g<=2. A contact at the middle0 is simple
and gives g<=1. A single leaf contact q with t=0 gives
g=8/q²-3: this is5 for q=1 and at most1 for q>=sqrt2. If t>=1 instead,
c0<=5/4 and h>=3/2, giving g<=5/9. These bounds include all allowed native
opposite products, not just simple ones. For the simple-at2 profile,
positive semidefiniteness of D10+c plus A3 gives0<=c0<=1 and h>=1,
hence g<=3. Additional spectators only strengthen these restrictions.

At the far end the isolated profile with t=1 and the simple-at13 profile
with t=0 both have k=1/2 and g=8. Every other far profile has g<=4:
the isolated nonunit case has t>=sqrt2; simple-at13 with t>=1 has g<=3;
and the doubled-at13 or simple-at9 cases have g<=4 for every t>=0.
Thus these are the complete two far extreme types.

Eliminating Z, with J=G_Z^(-1), gives the exact determinant identity

    0 = 1-(g_c+g_d)/10 + E/10,
    E=(x/h-y/k)^T J (x/h-y/k) >= 0.

The source verifies this identity formally, including the spectator
off-diagonal sign. If the far end is not extreme, g_c+g_d<=5+4=9, giving
a strict positive gap at least1/10. Hence the far end is extreme.

## Whole-factor payments with arbitrary spectators

First suppose the near end is not extreme. Then

    E=g_c+8-10=g_c-2<=1,     k=1/2.

If any y_i is nonzero, disjoint supports give x_i=0. Gram-metric coordinate
Cauchy, using G_ii=2, implies

    E >= (x_i/h-y_i/k)²/G_ii = 2*y_i² >= 2,

a contradiction. Thus the entire far vector y vanishes, without any
restriction on x or on arbitrary connected spectator tails. In the
chosen paying corner W+a+d, the opposite root c is absent.

* For the isolated far profile, d contacts only the omitted b. In the
  chosen corner it is a whole A1; F minus b is a whole finite D13, with
  a at the unique long-arm neighbor of its fork center. The new half
  coupon plus the selected ideal half gives1.
* For the simple-at13 far profile, d extends that D13 to a whole finite
  D14. The marked canonical pair is(0,10), of retained weight1. The
  current finite D14 refinement changes(0,11), so this unit survives.

The whole-factor assertions follow because Z was orthogonal to F and
y=0. Contacts belonging only to c do not enlarge this paying corner.
The source reconstructs both actual finite Grams and marked roles.

If both ends are extreme, h=k=1/2 and the determinant gives

    (x-y)^T J (x-y) = (5+8-10)/4 = 3/4.

Both separate spectator extensions are strictly positive definite, as
principal subdiagrams of actual edges. The complete
[below-one spectator theorem](fresh_spectator_energy_below_one_20261017.md)
therefore applies. At rational energy3/4 its only possibilities are a
single simple contact at an endpoint of one whole A3, with the other
support empty, or two simple contacts at adjacent nodes of one whole A3.
Every other component has zero contacts. In the two-nonempty case there
is an actual path of length3 between c and d, giving a distinct unit.

If the far support is empty, use the previous far-corner payment. If the
near support is empty, c extends the old A3 to a whole A4 with c an
endpoint. The old far component is a whole D10 with b a short leaf.
The existing
[small-component/D10 short-leaf coupon](fresh_large_finite_d_leaf_half_transfer_20261017.md)
has weight1/2 for every mark of a whole component of order at most5.
It applies to this whole A4 endpoint and combines with the selected half.
The other root d is absent from the chosen corner, so its arbitrary
spectator contacts cannot spoil either whole factor. The source checks
the literal A4/D10 roles; no A7 or additional cross-component coupon is
being assumed.

## Current rule, priorities and complete scope

The canonical current D14 table is bound to the completed dimension64
interface and its frozen successors: total82, center pair(0,9) zero, and
all four proposed pairs literal units. The native automorphism orbit is
checked exactly. The current
[B14 reduction](fresh_order14_b_basic_distance9_transfer_20261018.md)
and C14 reduction both have total80; their recorded family provenance
supplies the unchanged A70. Thus all four current family totals are
70,80,80,80 after this theorem.

For any actual quadrilateral, first apply the independent BC11-zero
theorem if a changed BC11 corner occurs. Its short1/2, small finite and
retained D11-half or typed-coupon outputs avoid every new D14 role.
Otherwise apply this theorem whenever a changed D14 corner occurs.
The two immediate near profiles, the far-empty cases and the complete
energy3/4 alternatives above give an independent full payment. If neither
occurs, the actual ideal weights agree pointwise with the preceding
proved rule and the finite coupon only increases, so the whole-face
proof is inherited. This priority concerns actual corners and does not
assume that an earlier selected witness avoids the new orbit.

The proof is universal and consequently covers occurrences in every
selected, protected, adaptive or low/disjoint catalogue branch, including
all74 protected bases. No absence from saved embeddings is substituted
for this coverage. In contrast to the B/C basic-distance9 deletions,
there is no new selected-map or protected-base census in this packet.

All ideal distances at most8 are unchanged literal units for this D14
factor. The independent current short8-to-short7-to-short6/whole-A8
route is therefore available to the other B14/C14 and finite14 proofs.
The new theorem itself uses only short1/2/3, retained affine D11 halves,
whole finite D14(0,10), A4/D10 short-leaf halves and the new A1/D13-neighbor
half. It uses no reduced finite D15--D17 same-component payment, no
A2/D13-center coupon, and no A1/A13-middle coupon. Its outputs therefore
do not depend on the separate proposed rebalances of those roles.

The elementary D14 triangle target set is explicitly unchanged. The
inherited finite7/ideal5 triangle theorem uses short or protected
small-factor outputs, and does not spend these distance9 near-fork pairs.
Triangle minimum2 and quadrilateral minimum1 persist. A fresh audit of
all216 current literal configurations, starting from the frozen finite
B14/D14 table, finds zero occurrences of either new ideal debit or new
coupon. This is an interface check; universal coverage rests on the
arbitrary-contact proof above, not on the216-row observation.

The only finite increment is exactly1/2 for the one whole A1/D13 role.
The existing A4/D10 coupon is already budgeted. No finite maximum is
asserted here, and the family envelope80/14 alone is not a dimension58
proof.

## Verification

Final generation session29119/chunka50bf9 and separate fresh replay
session47272/chunk0e7165 terminated exit0/PASS. The
checks verify the two inverse lists, four near and four far profiles,
all16 cross identities, determinant formula and native-bound endpoint
values, two immediate affine corners, three whole finite corner roles,
native orbit/triangle retention, current family provenance and216 literal
configurations. The all-native principal-extension classification and
arbitrary-rank below-one theorem are written mathematical inputs; these
are not inferred from finite checks.

Parent completed a full producer/report mathematical read and checked
the bound below-one theorem at3/4. The all-native principal profiles,
R10 physical identities, strict cost gap, entire-vector Cauchy argument,
whole finite roles, independent BC11 priority, retained D11 half and
triangle/current216 interfaces pass. A second reviewer independently
read the full producer/report, complete below-one report, analogous D13
center-profile and sealing proof, factor-fiber/determinant interfaces and
BC11-zero outputs. That read also passes the native bounds, complete
energy3/4 alternatives and unique new D13-node9 coupon. Neither reviewer
executed inverse, profile, literal, census or finite-count calculations.
Source, certificate and report are final/frozen after direct binding,
local-link and syntax checks. No external refereeing or proof-assistant
formalization is claimed.

    python3 scripts/fresh_order14_d_near_fork_half_transfer_20261018.py --verify
