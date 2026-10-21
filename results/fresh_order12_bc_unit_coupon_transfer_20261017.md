# Order12 B64 and C63 with whole A2/D11 unit coupons

The complete face-payment rule may remove the two remaining half-weights
of affine B11, and the two remaining nonbasic unit weights of affine C11.
Their order12 costs become64 and63 respectively. Upgrade the existing
half-coupon between a marked root of a whole A2 and a short fork leaf of
a separate whole finite D11 to a unit. The existing full A1/D11 short-leaf
coupon is retained. Every triangle still receives at least2 and every
quadrilateral at least1.

The order12 A/B/C/D costs are now60/64/63/62, so the normalized affine
rank11 envelope is16/3. The added finite cost is **2 per A2/D11 component
pair**, or4*a2*d11 in doubled units. This local theorem gives no finite
forest maximum or ambient dimension bound; the new coupon must be
included in any joint budget.

The [source](../scripts/fresh_order12_bc_unit_coupon_transfer_20261017.py)
and [certificate](fresh_order12_bc_unit_coupon_transfer_20261017.json)
reconstruct the complete selected support and exact local interfaces.
They retain the inherited all-label core classification and arbitrary
supporting-wall proofs stated below.

## Pair rules and complete support

For affine B11, number the spine0--...--9 with doubled bond8--9 and
short leaves10,11 attached to0. Remove its two half-weights at(9,10)
and(9,11), both distance10. The remaining64 pairs are exactly the basic
ideal5/feature5 set, all with weight1. This is the full leaf-swap orbit.

For affine C11, use the path0--...--11 with doubled bonds0--1 and10--11.
The endpoint pair(0,11) was already absent. Remove the reversal orbit
{(0,10),(1,11)}, both distance10. The remaining63 pairs are exactly the
basic ideal5/feature5 set, all with weight1. All short8 and elementary
triangle units remain for both types. No distance9 payment changes.

The checker re-enumerates all758 eligible connected high patterns into
both targets, including every pattern with no map, and compares every
successful map to the complete frozen support. The B search is exhaustive
injective induced-map recursion. The separate C algorithm tries every
contiguous path segment in both orientations, checking all labels and
nonedges. The affected support is:

|Type|Affected contexts|Affected induced maps|
|---|---:|---:|
|B|14|28|
|C|2|4|

The B contexts are

    98,100,119,121,164,167,171,1452,1728,1734,1745,2121,2127,2138.

The C contexts are164 and1452. Their abstract target maps are proper
extensions of selected finite B11 components. The next two sections
verify the actual geometry, not just a list of abstract marked images.

## B: exact whole finite alternatives replace the ideal halves

The [frozen B-half theorem](fresh_order12_affine_b_half_transfer_20261017.md)
already proves complete actual supporting-wall control for these14
contexts. The new checker reruns its literal inverse, inertia, contact
and corner reconstruction and compares every output to the frozen
certificate. Its arbitrary-support reasoning is essential here.

In eight radical contexts, the full indefinite core has signature(13,1,1).
Its radical is an actual linear relation: a null vector orthogonal to
a timelike vector in the Lorentz space must vanish. Every possible new
contact is in the checked nonzero support of that relation. Acuteness
forces a nonzero new contact vector to meet both signs, giving an unchanged
marked path of length at most3. In an otherwise unpaid face every extra
white is therefore orthogonal to the entire old core. The displayed
finite components remain whole after arbitrary supporting whites are
restored. There is no first contact from which an indirect tail could
start.

The other six cores are nonsingular Lorentzian of order14. With norm2
physical roots and saved gauge squares s_i, a native contact q at old
node i has projection square q^2*s_i*J_ii. The complete contact domain
splits into two ends; meeting both gives a retained short8 path. The
positive inverse blocks give a lower bound exceeding2 for every two-node
support on one end, including arbitrary larger contact magnitudes.
Exactly one singleton profile survives: a unit contact at an omitted
black root outside the selected component, with normalized projection2
and residual0. Labels at least4 are excluded.

In contexts164 and1452 no surviving contact can enter the selected B11.
It is a whole finite component and its marked pair retains a full unit.
These contexts cannot complete to an actual whole affine B11 or C11.

In the remaining four contexts there is either no such extra root or
exactly one. Its zero residual fixes its entire normal in the old
Lorentz span, so a second such normal is impossible under strict white
positive definiteness. Every core-orthogonal white is also orthogonal
to this forced root. Thus arbitrary further supporting whites cannot
attach through it or bridge to the displayed finite components. The
source reconstructs both actual possibilities, using the correct norm2
interpretation before applying the saved positive diagonal gauge.

There are16 finite-payment configurations: eight sealed radical cores
and both variants of four nonsingular cores. Their finite alternatives
are whole B12 units, whole A1/D11 short-leaf coupons, or whole A2/D11
short-leaf coupons. A1/D11 is already full by the
[universal D11 theorem](fresh_affine_d11_leaf_pair_transfer_20261017.md).
Making A2/D11 full allows the preferred ideal B11 contribution to vanish.
The complete finite totals are:

|New finite total|Configurations|
|---:|---:|
|1|8|
|2|8|

The checker uses the actual connected components and actual marked roles,
not a proper D11 subdiagram. Twelve configurations use the new A2 upgrade
and four use the existing A1 upgrade. All same-component units have
order at most13. The two already-finite B11 contexts separately retain
their full unit.

## C: the missing terminal is impossible in both contexts

Only contexts164 and1452 map a selected component to one of the removed
C pairs. The preceding argument already proves that their selected B11
component is whole and finite. A direct first-entry exclusion confirms
the four abstract C maps.

Each map omits one terminal of the12-node C path. The whole target would
require a new norm2 root with a label4 contact to an old marked node i.
The root is genuinely new: omitted old black roots are absent from this
corner, while every other old corner component is orthogonal to the
selected component. Its native squared contact magnitude is2. In every
map the correctly scaled bound is

    2*s_i*J_ii = 16 > 2.

Within one end, entrywise positivity means further contacts can only
increase this bound. Contacts on the other end give an unchanged short8
path. Thus no additional contacts, labels or indirect spectator branches
can repair the missing-terminal obstruction. No unspecified contact has
been assigned zero and no norm4 gauge vector has been mistaken for a
physical norm2 root.

## Every nonselected and later interface

The complete low-core closure has66 chosen cases and69 actual variants.
Its only affine B order12 images are basic retained pairs; it has no C
order12 output. These are closed actual components, not arbitrary
extensions of finite patterns. The disjoint high-core branch remains
excluded. Protected ideal outcomes have distance at most9 and remain
paid, while all finite protected outcomes are unchanged by this theorem.

The sole additional B provenance is the universal affine-D13 `inner4`
immediate alternative. The frozen B-half proof already excludes it for
all native opposite contact magnitudes and all spectator forests: its
end cost is at most2, while the opposite end costs at most6 or8, giving
strict determinant gaps1/9 or1/11. The nonnegative spectator correction
cannot cancel them. The checker reruns all three exact gap formulas.
No C immediate alternative exists. Universal D14 immediate factors have
order13 and are whole by their positive affine kernels.

The frozen C endpoint report audits every conclusion of the earlier
universal and local transfers. Its reasoning still applies because the
two newly removed pairs have distance10: cycle alternatives are finite
or short ideal paths; classical-tree alternatives have distance at most9
with persistent nearby features; all equality, partial-weight and
singleton alternatives have their stated finite components, separate
coupons or other whole affine factors. The cycle13/14 and longer-cycle
proofs therefore cannot introduce either distance10 C pair.

The later cycle12 and D11 transfers return finite alternatives or ideal
paths of length at most4 and5 respectively. The B14 projection refinement
introduces no replacement payment. The finite rank17 refinement changes
no ideal pair, and its finite deletions do not touch the components here.
Consequently there is no later nonselected B/C order12 obligation beyond
the support and immediate cases just treated.

The new large-D/D10 transfer changes only whole D_n with n>=20 and
separate small-component/D10 coupons. None of those halves or coupons
is used by the present B/C proof. Conversely, the large-D proof uses
short8/basic ideal paths and unchanged protected ideal units, so these
distance10 deletions do not debit it. Both local theorems compose, but
their finite coupons must be counted together.

## Finite cost and verification scope

A whole A2 has two marked positions and a whole D11 has two short fork
leaves. Four pair weights increase from1/2 to1, adding2 per component
pair. The fork-center coupons and every other old coupon are unchanged;
including its two fork-center halves, the full A2/D11 coupon total is5.
No coupon has been counted by a rank-only surrogate in the local proof.

Generationefbe5b and separate fresh replayd26928/b55e8a both terminated
exit0/PASS. The checker re-enumerates the complete758
pattern domain for both targets, reruns all14 literal B geometry
certificates and16 finite variants, checks the four C forced-terminal
projections, all basic/short8 sets and automorphisms, the low-domain
interfaces and three universal immediate gaps. The complete all-label
core census, finite/affine classification, factor-fiber theorem and
arbitrary-support sealing remain inherited frozen inputs.

A separate full source/report mathematical and interface read passed.
That review checked the inherited B-half geometry reconstruction, all16
finite variants with whole A2/D11 roles, both whole-finite-B11 C contexts
and scaled missing-terminal bounds, complete selected/low/immediate/later
interfaces, and compatibility with the large-D/D10 transfer without
consuming its coupons. It did not repeat the core census or numerical
geometry calculations.

Root completed a full source/report read and reread the frozen B-half
source/report and C endpoint report. That review passed all758 successful
and failed pattern cases, the16 finite variants with actual whole A1/A2
and D11 roles, scaled projection16, arbitrary radical and zero-residual
sealing without indirect bridges, all retained short/basic/low/immediate
and later interfaces, and independent composition with the large-D rule.
Root did not duplicate the geometry or census arithmetic. All37 direct
bindings, Python AST and four local report links pass. Source, certificate
and report are final/frozen.

No exact joint finite-forest maximum, ambient bound, polytope construction,
external refereeing or formal verification is asserted here.

    python3 scripts/fresh_order12_bc_unit_coupon_transfer_20261017.py
    python3 scripts/fresh_order12_bc_unit_coupon_transfer_20261017.py --verify
