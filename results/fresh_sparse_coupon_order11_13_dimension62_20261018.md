# Unrestricted finite-volume dimension bound62

The simultaneous payment rule below excludes dimension63. Together with
the completed higher-dimensional exclusion, it gives dimension at most62
for arbitrary finite-sided finite-volume hyperbolic Coxeter polytopes.
No arithmetic, compactness or simple-cusp hypothesis is imposed.
**Completed internally checked computer-assisted proof; final verification
and review scopes are recorded below. Source, data and report are frozen.**

At ambient63 its exact comparison is

    0 >= (1/2)*v_f + (139/46)*v_i +4960,

which is impossible. Here v_f and v_i count finite and ideal vertices.
The [source](../scripts/fresh_sparse_coupon_order11_13_dimension62_20261018.py)
and [certificate](fresh_sparse_coupon_order11_13_dimension62_20261018.json)
bind the completed local proofs and independently checked finite/ideal
arithmetic, and verify their precise simultaneous interfaces. They do not
rerun the historical geometric censuses, finite-forest dynamic programs
or2,204-profile scalar computation.

## The single simultaneous payment rule

Begin with the
[completed dimension63 rule](fresh_order11_order17_dimension63_20261018.md).
Make exactly these further ideal changes:

* whole affine C10 of order11: endpoint pair(0,10), weight1/2 to0, by the
  [C11 zero theorem](fresh_order11_c_endpoint_zero_transfer_20261018.md);
* whole affine B11 of order12: fork-center/doubled-terminal pair(0,9),
  weight1 to0, by the
  [B12 zero theorem](fresh_order12_b_center_terminal_transfer_20261018.md);
* whole affine D12 of order13: fork-center pair(0,8), weight1 to1/2, by the
  [D13 half theorem](fresh_order13_d_center_pair_half_transfer_20261018.md).

The order is the number of roots, one more than affine rank. The binder
reconstructs all four factor-family totals at each changed order:

|Order|A|B|C|D|Maximum|
|---:|---:|---:|---:|---:|---:|
|11|99/2|54|54|53|54|
|12|60|63|63|62|63|
|13|65|73|73|147/2|147/2|

All other ideal weights, including cycle11 halves and the order17 D
half-coupon improvement, retain their completed values. The new normalized
envelopes are54/11,21/4,147/26 at ranks10,11,12. Every other rank envelope
is checked against the preceding complete ideal calculation.

For finite corners, adopt the
[sparse whole-component theorem](fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018.md).
It gives A16=114, B16=104, D16=114, D17=125 by deleting six units, all
at marked distance at least10, and adding six typed component-pair
coupons. Retain every preceding coupon and the deep11 finite8 zeros.
Additionally include the B12 theorem's A2/B9 doubled-terminal half and
the D13 theorem's A2/D10-center and A3-middle/D10-center halves. No
payment is assigned merely because a component contains a named type
as a proper subdiagram: all these roles concern separate **whole** finite
components.

The nine added component-pair totals, in the order used by the budget,
are A4/A12:4, A4/A13:4, D5/A12:2, A2/A13:2, A2/A14:2, E6/A11:4,
A2/B9:1, A2/D10:1, A3/D10:1/2. The source compares every total with its
local role certificate and the final budget. The A3 coupon charges only
its unique middle. Existing coupons on the same type pairs have disjoint
marked roles and are retained. The enlarged E6/D11 total4 used in the
budget is only a conservative representative allowance: the actual
E6/D11 rule still costs2.

## Exhaustive actual-quadrilateral priority

Apply the following cases to the actual corners of each quadrilateral.
They do not refer to which corner an older algorithm selected. All
payments below use the final simultaneous rule.

**1. Some actual corner has the changed C11 endpoint pair.**
The independent C11 theorem analyzes every norm-two heavy-end profile,
including the affine alternatives. Its strict determinant gap at least1/5
forces a distinct actual adjacent-black or common-spectator path of
length1 or2. That corner retains a unit under every ideal and finite
change here. This branch consumes neither the deleted C11 half nor
another potentially changed affine half. It supplies total at least1.

**2. Otherwise some actual corner has the changed B12 center/terminal
pair.** The universal B12 theorem covers all supporting spectators and
native labels. Its possible outputs are a whole affine D11 factor of
order12 marked at a center and a short leaf, distance8; a whole finite
D12 center-to-long-endpoint distance9 unit; or two distinct finite
half-coupon corners, whole A2/B9 doubled-terminal and whole A2/D11
fork-center. The affine pair is still a unit: it is not an opposite-leaf
pair, and the new D13 change is in a different order. The finite D12
pair is unaffected by sparse rank16/17 distance-at-least10 deletions or
deep distance8 zeros. Both whole-component coupons are included in the
budget. The source checks the six literal boundary cases and the four
canonical affine center/leaf unit pairs. Total payment is at least1.

**3. Otherwise some actual corner has the changed D13 center pair.**
Keep its original half. The universal D13 theorem gives either a distinct
whole A2/D10-center half, a distinct whole A3-middle/D10-center half,
or a distinct retained unit. Its unit alternatives are whole order11
affine B10 center/terminal distance8, whole order11 affine D10
center/leaf distance7, or a path of length at most4. The two affine
center pairs are not the B/D11 terminal-leaf halves, and no C11 endpoint
or B12 center/terminal pair is used. All finite outputs are typed whole
D10 coupons or short paths; none uses a finite rank16/17 long unit.
The binder checks all four immediate affine unit roles and eight finite
coupon corners against the final weights. The half plus distinct half,
or the independent unit, gives total at least1. The local theorem applies
its bound inherited arbitrary-rank energy-one classification and proves
the new whole-component/corner consequence. Neither is inferred from
these literal records.

**4. No actual corner has any of the three changed ideal pairs.**
Apply the complete sparse finite **whole-face** theorem under the
baseline63 ideal weights. Its ideal weights now agree pointwise with
the final weights at every actual corner of this face. Its finite rule
is the final finite rule except for the B9/D10 coupon additions, which
can only increase payment. Therefore its total at least1 is inherited
directly. There is no requirement to predict which alternative corner
the sparse proof will select, and no circular appeal to a payment that
may have been removed elsewhere in the face.

These four cases are exhaustive and all handle the actual whole
components and actual corner marks. In particular, the sparse theorem's
old short8 alternatives need not be asserted to remain literal units
for every new ideal pair: any actual changed ideal corner was handled
first by its independent theorem.

## Scope of the inherited sparse finite theorem

The sparse theorem is a complete geometric input, not an inference from
its numerical forest bound. Its selected-support proof covers the
complete2,614-context high-overlap domain:394 eligible A16 patterns with
238 maps,394 D16 patterns with578 maps, and140 D17 patterns with168 maps.
Six proper D16 contexts are excluded by exact norm-two projections,
two-leaf Cauchy, or a negative fork residual; arbitrary additional
contacts are controlled before zeros are imposed. Ten radical contexts
give the new coupons, including the four obligations previously covered
by restored D17 units. Their zero radical and all-label contact domain
seal entire displayed components unless an actual short8 path occurs.

All216 inherited literal configurations are debited simultaneously;
the minimum remains1. This finite list is not treated as every
quadrilateral. The inherited low-core bound12, excluded high-disjoint
case, protected/adaptive leaf9 outputs, immediate and universal
cycle/tree branches, and later whole-component or typed-coupon
interfaces are retained separately in that theorem. In particular,
its two adaptive bound16 outputs have marked distance9, missing every
new deletion. Later baseline63 order11 outputs use components of order
at most14, protected A9/A18 endpoint8 or D18 leaf9, E8 path6, or typed
coupons; order17 uses no rank16/17 same-component unit.

If an inherited short8 path ends at a deep finite8 zero, the independent
[depth11 zero theorem](fresh_depth11_finite8_zero_transfer_20261018.md)
forces a short7 alternative. Under baseline63 weights its complete
small-tree, cycle11, cycle-tail and protected-feature interfaces pay.
The sparse deletions have distance at least10 and do not change those
outputs. This route is part of the sparse theorem under its stated
baseline ideal rule; the actual-corner priority above transfers it to
the final rule without reopening or circularly modifying the route.
The present author read the complete sparse source/report and its
inherited selected/literal exchange routines; the underlying all-label
censuses and published classification retain their inherited proof scope.

## Triangles and complete local payment

The [finite7/ideal5 triangle theorem](fresh_finite7_ideal5_triangle_payments_20261017.md)
uses finite distances at most7, ideal distances at most5, protected ideal
leaf distances6/7, or whole ideal factors of order at most9. The three
new ideal pairs have distances10,9,8 and orders11,12,13, respectively;
none is used by that triangle theorem. The finite deletions have distance
at least10, so triangle payment at least2 is unchanged. The previously
proved cycle11 triangle exclusion and all prior interfaces remain in
the baseline. The new coupons only add weight.

Thus every triangle receives at least2 and every quadrilateral at least1.
Weights are nonnegative, so faces with at least five sides need no
positive payment to bound5 minus their number of sides. Summing all
actual corner payments yields an upper bound on

    D5 = sum_(two-faces Q) (5 - number_of_sides(Q)).

## Exact finite and ideal vertex budgets

The frozen [joint finite budget](fresh_sparse_coupon_rank63_forest_budget_20261018.md)
includes all old and nine new coupons, with the exact marked-role
increments above. Its complete1,329 large multiplicity vectors and993
residual tables give maximum449 at rank63. Pointwise simultaneous
two-end type domination retains B9, E6 and A17 distinctly; its larger
E6/D11 allowance is an upper bound only. Every Bellman upper inequality
and an attaining predecessor are checked, and the actual formal forest
2D17+D29 has cost125+125+199=449. This proves the bound for all finite
vertex diagrams; no assertion that this witness occurs in a chamber
is needed. The curvature at a finite vertex is899/2, leaving surplus1/2.

At an ideal vertex, its Euclidean link is a product of simplices with
factor ranks summing to62. The
[independent ideal theorem](fresh_ambient63_ideal_margin_closure_20261018.md)
uses exactly the family envelopes above, together with all unchanged
rank envelopes. Its complete concentration domain has2,204 profiles:
arbitrary repeated ranks9 through17, at most one smaller positive rank,
and at most one rank at least18. Every normalized margin is positive;
the minimum is3505/201344 at(10,10,12,15,15).

The concentration moves preserve cost, decrease normalized curvature
and decrease the product V of factor orders. Since terminal margins
are positive, the absolute margin V times normalized margin also
decreases. Its exact minimum139/46 at the single rank62 factor is
therefore valid for every actual cusp profile, not just the terminal
list. The factor-fiber formula counts each marked varying-factor pair
with its exact multiplicity V/(r+1). Consequently the ideal payment C_v
satisfies

    kappa(L_v) - C_v >=139/46.

No cusp simplicity or arithmetic condition has been inserted. The
independent scalar theorem was explicitly conditional on the three
new envelopes; the local priority proof above now discharges those
geometric hypotheses, and the complete finite budget includes their
coupons.

## Global contradiction and higher dimensions

Let d=62. The
[shape-curvature theorem](fresh_cusp_link_shape_curvature_20261017.md)
gives the local curvature kappa(L)=(5G-dV)/2. Its extremum correction
can be strengthened for these links as follows. Every coefficient h_j
of a simplex or product-of-simplices H-polynomial is at least1. At each
of the two global extreme vertices,

    E-G = sum_(j>d/2) h_j*(2j-d)
        >= sum_(j>d/2)(2j-d) = floor(63²/4) =992.

The two extrema, each multiplied by5/2 in D5, thus contribute4960.
This corrects the same curvature lower bound once; it is not an extra
vertex payment or a second copy of a previously counted constant.

The local upper bound and global lower bound are therefore

    D5 <=449*v_f + sum_(ideal v) C_v,
    D5 >=(899/2)*v_f + sum_(ideal v) kappa(L_v) +4960.

Subtracting gives the impossible inequality at the start. It excludes
ambient dimension63 regardless of the numbers of finite or ideal
vertices. The completed unrestricted ceiling63 already excludes every
dimension at least64. Combining the two exclusions gives **dimension
at most62**. Compact and arithmetic-specific bounds are unchanged.
No new Coxeter polytope example is constructed by this argument.

## Verification and proof status

Source-bound generation4b3b23 and a separate fresh saved replaye48d4f
terminate exit0/PASS on12 complete
family interfaces, exact changed-pair identities, protected alternative
roles, all nine coupon increments, frozen budget/scalar bindings, finite
curvature and the strengthened extremum constant. Direct audit318d0f
passes42 source bindings,12 local links and Python syntax.

Root and two separate peer agents completed full source/report
mathematical composition reads. All three checked the exhaustive
actual-corner priority, pointwise agreement in the fourth branch,
protected unit and whole-component half outputs, all nine exact coupon
increments, conservative versus actual E6 allowance, triangle/deep11
interfaces, family envelopes, absolute ideal surplus, finite curvature,
single extremum correction and completed higher tail. Root separately
read the contributing local, scalar and budget proofs at their recorded
scopes. The first peer also read the frozen ideal, baseline63 and sparse
reports; the second completed full C11/B12 and sparse-budget source/report
reads in addition to its earlier D13, sparse-local and ideal reviews.
None of these composition reviewers reran the historical census,
Bellman tables or2,204 scalar rows. The inherited energy-one attribution
was clarified in the report; no mathematical or arithmetic change was
required. Source, data and report are final/frozen.

The constituent local and counting packets already have completed
generation/fresh replay and mathematical reviews at their recorded
scopes. This binder is a local-to-global composition check, not a new
independent enumeration of every historical geometric core. The
underlying finite/affine classification, all-label reduction and
all-rank concentration proofs remain explicit inherited mathematical
inputs. All review is internal; no external refereeing or formal
proof-assistant verification is claimed.

    python3 scripts/fresh_sparse_coupon_order11_13_dimension62_20261018.py
    python3 scripts/fresh_sparse_coupon_order11_13_dimension62_20261018.py --verify
