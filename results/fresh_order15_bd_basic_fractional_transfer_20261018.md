# An economical order15 B/D transfer

The current ideal order15 family totals A75/B92/C91/D92 may be changed
simultaneously to **A75/B91/C91/D367/4**. Thus the all-family envelope
falls from92/15 to367/60, a total saving1/4. The B change is a unit
deletion. The D change reduces four units to15/16 and adds exactly one
finite coupon of weight1/16. This is a local payment theorem; it does
not assert a new finite maximum or dimension58 bound.

The new coupon joins a whole A1 to a separate whole finite D14, marked
at its unique long-arm vertex two edges from the fork center. In the
canonical finite D14 with spine0--...--11 and short leaves12,13, the
large mark is9. There is one eligible marked pair per component pair.
Its actual component-pair cost therefore rises from2 to33/16. No other
finite weight changes. The whole-component condition is essential.

The [producer](../scripts/fresh_order15_bd_basic_fractional_transfer_20261018.py)
and [certificate](fresh_order15_bd_basic_fractional_transfer_20261018.json)
bind the current [order14/contextual rule](fresh_order14_contextual_ambient59_frontier_20261018.md),
the [order15 B92/D92 predecessor](fresh_order15_bd_terminal_coupon_transfer_20261018.md),
and the universal [D14 near-fork theorem](fresh_order14_d_near_fork_half_transfer_20261018.md).
The D argument here uses a different, strict D4 end and energy1/2;
the old A3-end/energy3/4 argument is not assumed to apply unchanged.

## Canonical changes and the complete B domain

For affine B14 of order15, use spine0--1--...--12, doubled bond11--12,
and simple leaves13,14 at0. Delete the retained singleton(3,12). It has
distance9 and is invariant under the fork-leaf interchange. All earlier
zeros remain zero, and the resulting total is91.

For affine D14 of order15, use spine0--1--...--10 with leaves11,12 at0
and13,14 at10. Change precisely

    (2,13), (2,14), (8,11), (8,12)

from1 to15/16. They form one native automorphism orbit of distance9.
Its total debit is1/4, giving367/4. The earlier distance10 deletions,
including(1,13),(1,14), remain zero. The source reconstructs the complete
weighted tables, checks all native automorphisms, and retains every
recorded short8 and elementary-triangle unit.

The B singleton is basic, so a selected-map absence alone is insufficient.
The source reconstructs all536 eligible native embedding searches from
the complete high-overlap domain, including failures. There are exactly
215 successful patterns and548 maps, identical to the frozen predecessor;
none sends the marked pair to(3,12). It also freshly embeds each of the
74 complete protected native bases into the whole B15 target, not merely
their previously selected maps. None uses the changed pair. This covers
arbitrary supporting extensions of a protected base because its displayed
native diagram remains an induced subdiagram of any actual whole factor.

The inherited2614-context partition is1691 whole-order-at-most14 cases,
83 protected/adaptive cases,304 displayed-too-large cases and536 eligible
patterns. The83 consist of the74 protected bases and nine independently
paid adaptive finite cases. Low-core closures have order at most12;
the high-disjoint obstruction, complete low/protected routing and the
universal short-path routes are inherited inputs. The nine adaptive
outputs are finite and unchanged. Consequently neither a low closure,
an eligible selected payment nor a protected ideal payment can use this
new B singleton; the too-large branch cannot be a whole order15 factor.
The nonselected/current interfaces are detailed below. No new finite
coupon is needed for B.

## The actual D15 star and its complete profiles

Let an actual quadrilateral star be W+a+b+c+d, with opposite pairs
{a,c},{b,d}. Actual edge diagrams W plus one black root are strictly
positive definite; actual corner diagrams are positive semidefinite.
Use physical roots of square2. The inherited factor-fiber argument at
W+a+b supplies the whole varying affine factor F and a positive-definite
spectator diagram Z orthogonal to it, without a simple-cusp hypothesis.
It suffices to take a=2 and b=13 in F. Its positive kernel has coefficient1
on the spine and1/2 on each leaf.

If b.c, a.d or c.d is negative, a distinct actual corner has a retained
distance1 unit. If c,d share a nonzero spectator contact, a distance2
unit pays. Hence assume those products vanish and the nonnegative
spectator contact vectors x,y have disjoint supports. The opposite
products a.c and b.d remain unrestricted. Every nonzero native magnitude
is at least1; one greater than1 is at least sqrt2. This includes parallel
and dotted opposite products. Kernel heights are positive.

Deleting a leaves finite D4 on[1,0,11,12] and finite D10 on
[3,4,...,10,13,14]. Deleting b leaves finite D14, with long endpoint14.
In both large finite blocks, the inverse diagonals at the first two
long-arm vertices are1,2, and every other diagonal exceeds2. Thus a
positive-semidefinite one-root extension permits only no contact,
simple or doubled contact at the first vertex, or simple contact at the
second. Two contacts create a cycle in a connected diagram still
containing the fork, excluded by finite/affine classification; the
exceptional types are too small to contain D10. This is an all-native
classification, not a finite label cutoff.

A doubled-first or simple-second extension is affine with positive
kernel. In the actual PSD corner that kernel makes the entire component
orthogonal to every other corner root. Thus the corresponding near
profiles have no additional D4 or spectator contact. Their product with
the omitted opposite black a is still unrestricted.

The D4 inverse, in the stated order, is

    K = [[1,1,1/2,1/2], [1,2,1,1],
         [1/2,1,1,1/2], [1/2,1,1/2,1]].

Put t=-a.c, u=-1.c, v=-0.c, w=-11.c, z=-12.c and

    e=(u,v,w,z) K (u,v,w,z)^T,
    h0=t+u+v+(w+z)/2.

The strict actual edge contains D4+c, so e<2. A center contact already
has inverse diagonal2 and is impossible. Any two nonzero D4 contacts
have projection square at least3, also impossible. Each outer inverse
diagonal is1; therefore the only possible D4 contact is one simple
outer contact. In particular e is0 or1. This strict argument covers
arbitrary additional old and spectator contacts, whose principal
restrictions cannot weaken it.

Eliminating F minus a,b gives resistance R=10. The four near profiles
have height h, residual constant c0 and moment l as follows:

|Near profile|h|c0|l|
|---|---|---|---|
|Only D4/opposite contacts|h0|2-e|0|
|Simple contact at3|h0+1|1-e|1|
|Doubled contact at3|t+sqrt2|0|sqrt2|
|Simple contact at4|t+1|0|2|

For the far root put t=-b.d. Its complete profiles are:

|Far profile|k|d0|m|
|---|---|---|---|
|No contact to F minus b|t/2|2|0|
|Simple contact at14|(t+1)/2|1|1|
|Doubled contact at14|(t+sqrt2)/2|0|sqrt2|
|Simple contact at10|t/2+1|0|2|

Before eliminating Z, the new-root residual squares are c0+l²/10 and
d0+m²/10, and their cross projection is lm/10. The producer checks the
two complete inverse lists, all eight end identities, all16 cross
identities and the formal four-by-four determinant directly in the
physical norm-two Gram.

## Independent affine alternatives and native cost bounds

The last two near profiles have immediate distinct actual alternatives.
A doubled contact at3 gives whole affine B10 of order11, with marked
doubled-terminal/fork-leaf pair. That pair currently has weight zero;
apply the independently proved BC11-zero theorem at this actual corner
first. It supplies a full retained payment. A simple contact at4 gives
whole affine D10 of order11 with opposite fork-leaf marks, of current
weight1/2. That distinct half plus the selected15/16 is23/16. The source
checks both kernels, canonical native embeddings and current weights.
An old zero is never spent as a unit.

For each end set g=2l/h+c0/h². After those immediate near alternatives,
the complete near bound is g<=4. The unique extreme type is a simple
contact to old leaf11 or12, no opposite contact, h=1/2 and g=4. Every
other near profile has g<=3. Indeed, in the outer-only profile:

* no D4 contact forces t>=1 and gives g<=2;
* a simple contact at1 has h=t+1 and g<=1;
* a simple contact at11 or12 has h=t+1/2 and g=4 at t=0,
  whereas t>=1 gives g<=4/9.

For the simple-at3 profile, e is0 or1, hence0<=c0<=1 and h>=1,
so g<=3. These cover every allowed native opposite magnitude and every
strict D4 contact, not only endpoint substitutions in a polynomial.

At the far end, the isolated profile with t=1 and the simple-at14
profile with t=0 both have k=1/2 and g=8. Every other far profile has
g<=4: the isolated nonunit case has t>=sqrt2; the simple-at14 case with
t>=1 has g<=3; and the last two profiles have g<=4 for every t>=0.

Let J=G_Z^(-1). Elimination and the exact determinant identity give

    0 = 1-(g_c+g_d)/10 + E/10,
    E=(x/h-y/k)^T J (x/h-y/k) >= 0.

The spectator cross sign is part of the checked formal identity. If the
far end is nonextreme, g_c+g_d<=4+4=8, so the right side is at least1/5.
Thus the far end must be one of its two extreme types.

## Entire-support Cauchy and the energy1/2 boundary

Suppose first the near end is not extreme. Then k=1/2 and

    E=g_c+8-10=g_c-2<=1.

If y_i is nonzero, disjoint supports give x_i=0. Gram-metric coordinate
Cauchy and G_ii=2 imply

    E >= (x_i/h-y_i/k)²/G_ii = 2*y_i² >= 2,

a contradiction. Thus the **entire** far spectator vector is zero.
There is no assumption on the near vector or on the rank and shape of Z.
In the paying corner W+a+d the opposite root c is absent.

For the isolated far profile, d contacts only omitted b. It is a whole
A1 in this paying corner; F minus b is whole finite D14, with a at
canonical node9. The new coupon supplies1/16, in addition to the
selected15/16. For the simple-at14 profile, d extends the old D14 to
whole finite D15. Its marked pair is canonically(0,10), of retained
weight1. The current finite D15 reductions change other long pairs,
including(0,11),(0,12),(1,12), and do not remove this unit. The producer
checks the two literal whole-component roles and that exact current unit.

If both ends are extreme, h=k=1/2 and the determinant instead gives

    (x-y)^T J (x-y) = (4+8-10)/4 = 1/2.

Both separate spectator extensions are strictly positive definite,
being principal subdiagrams of actual edges. The complete arbitrary-rank
[below-one spectrum theorem](fresh_spectator_energy_below_one_20261017.md)
therefore applies. At energy1/2 its rational rank is r=1: there is exactly
one simple contact on one whole A1 component, and the other support is
empty. The two-adjacent-contact alternative cannot occur in A1, and the
irrational I2 alternatives have different energies. Every other
spectator component has zero contacts.

If the far support is empty, use the preceding far-corner payment.
If the near support is empty, c extends the old D4 by a simple contact
to leaf11 or12, giving whole D5 with c its unique long endpoint. The
other whole factor in W+b+c is D10, with b a short fork leaf. The
existing small-component/D10 coupon pays1/2 for any mark of a whole
component of order at most5 against a whole D10 short leaf. This gives
at least23/16 together with the selected corner. The other root d is
absent here, so its arbitrary contacts cannot enlarge either paying
factor. The exact D5-long/D10-shortleaf roles are reconstructed; no
new D5 coupon is being assumed.

More generally this proof allows every0<delta<=1/2: use selected weight
1-delta and the one new A1/D14 weight delta. The immediate D11 and old
D5/D10 halves still give at least1, the isolated far branch gives
(1-delta)+delta=1, and the remaining branches give independent full units.
The chosen delta=1/16 is the claimed rule. Other parameter values do not
inherit any finite-budget assertion from this packet.

## Simultaneous current payments and finite accounting

Apply the independent BC11-zero theorem first whenever such an actual
corner occurs. Its short1/2, small finite and retained D11-half outputs
avoid the new order15 roles and the new A1/D14 coupon. Next apply the
universal D15 proof whenever an actual newly reduced D15 corner occurs.
All of its payments above are independent of the B15 deletion. No new
ideal B15 pair is used in a D branch.

Otherwise the D ideal weights agree pointwise with the preceding rule.
Use the independently established current short8 interface, followed
by the complete inherited B selected/protected/low/disjoint routing.
The new B singleton cannot occur in any of its complete eligible or
protected maps. The universal cycle/radius-five tree replacements use
their protected short/basic paths, small whole factors and fixed-leaf
or fork outputs, not this degree2/doubled-terminal distance9 singleton.
The earlier D13/D14 immediate affine outputs have orders at most13;
the new D14 theorem uses order11, whole finite D14 or A4/D10 payments.
The current C/B projection patches have no new long order15 output.
The independent actual-corner small-ideal priorities and the later
sealed-factor transfers use their explicit retained finite roles.
Thus no earlier preferred witness is silently assumed to survive.

If a short8 finite alternative is changed by depth8 zero, its independent
short7 route applies; a changed finite7 alternative is covered by the
depth6 theorem's independent short6 or whole-A8/protected-leaf alternative.
These routes keep the new ideal distances9 unchanged as unused outputs.
The short6 cycle11 exception uses its independently proved second half
or finite short alternative. No new D15 proof branch needs a reduced
finite7/8 payment. The contextual A2/A13-middle, A2/D13-center and
conditioned A2/A2 changes are disjoint from every new output. The new
D5/D10 role is an existing broad small-component coupon, distinct from
all D5/D12 changes. Triangles retain their elementary short/protected
targets, checked in the complete weighted tables. Quadrilateral total1
and triangle total2 therefore persist for the simultaneous local rule.

The producer also reconstructs all216 current literal configurations
starting from the frozen contextual rank13 rule, whose D14 successor
has no literal change. It tests complete whole-factor native embeddings
and exact new coupon roles, and finds zero new ideal debits and zero
new coupon credits. All recorded corner lower bounds are nonnegative
and each face total remains at least1. This is supplementary interface
validation; the arbitrary-support D proof and full B routing give the
universal coverage.

The old actual A1/D14 cost is2, and the new unique node9 mark adds1/16.
Its total is33/16, explicitly including the other pre-existing roles.
No A2--A5/D16 coupon is restored, no same-component finite weight is
increased, and no finite budget is reused without this extra charge.
The ideal factor-fiber count is V/15 per pair, giving envelope367/60.
The completed unrestricted bound remains59. A complete finite count
and a new terminal-ideal comparison belong to a separate composition.

## Verification and scope

Generation session41241/chunk0c3a92 and separate fresh saved replay
session29179/chunke9816b both terminated exit0/PASS. The source checks complete B native-map
and protected-base domains, exact D4/large inverse lists, eight end and
16 cross identities, the formal determinant, native-bound endpoint
values, two immediate affine and three whole finite corners, current
weight/orbit/triangle interfaces and all216 literals. The complete
all-native extension classification, factor-fiber/actual-star framework
and arbitrary-rank spectator theorem are bound mathematical inputs;
they are not inferred solely from these finite checks.

Parent completed a full producer/report mathematical read. The strict
D4 contact reduction, complete near/far cost bounds, physical resistance
and spectator identity, entire-vector Cauchy argument, energy1/2 use,
whole finite roles, BC11 priority, complete B selected/protected domains,
current literal/short6 interfaces and exact1/16 cost pass. That review
did not execute calculations or repeat the historical census or the
inherited spectrum proof.

A second reviewer completed a full producer/report mathematical read
and reread the complete below-one spectrum source/report. The strict
D4 contact reduction, complete large-block/native-height bounds,
entire-vector Cauchy argument, energy1/2 whole-A1 classification, distinct
finite payments, BC11 priority, complete536+74 B support and current
short6/216 interfaces pass. The valid delta range and chosen1/16 cost
are correct, and the D5/D12 alternate outputs are disjoint. This review
did not execute arithmetic, maps, the spectrum checks or the census.

Source, certificate and report are final/frozen after both completed
mathematical reads, terminal generation/fresh replay and direct metadata
validation of120 immutable bindings, six local links and Python syntax.
No historical geometric census, global finite optimization, external
refereeing or proof-assistant formalization is asserted.

    python3 scripts/fresh_order15_bd_basic_fractional_transfer_20261018.py --verify
