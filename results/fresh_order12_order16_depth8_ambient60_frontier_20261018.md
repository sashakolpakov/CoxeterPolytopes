# Ambient60 ideal closure and the remaining finite maximum409

The combined rule has **positive surplus at every ideal vertex**, and
its exact formal finite-forest maximum is **409**, against curvature405.
Every finite-sided, finite-volume real hyperbolic Coxeter polytope of
dimension60 therefore satisfies

    4 v_f ≥ 4500 + 54 v_i.

Here v_f and v_i count finite and ideal vertices. Consequently all-ideal
polytopes in this class have dimension at most59, using the completed
unrestricted dimension60 theorem for the higher-dimensional tail.
The unrestricted bound remains60. The formal maximum409 is attained
by 4A1+4D5+3D12; this is not a geometric realization of a vertex or polytope.

The [producer](../scripts/fresh_order12_order16_depth8_ambient60_frontier_20261018.py)
and [certificate](fresh_order12_order16_depth8_ambient60_frontier_20261018.json)
combine the frozen
[order13 complete rule](fresh_order13_sparse_ambient60_frontier_20261018.md),
[order12 B/C theorem](fresh_order12_bc_near_terminal_half_transfer_20261018.md),
[order16 B/D theorem](fresh_order16_bd_sparse_unit_transfer_20261018.md)
and [depth8 finite7 theorem](fresh_depth8_finite7_profile_transfer_20261018.md).
All previous same-component rules and whole-component coupons remain
present except for the explicitly stated new deletions below.

## The one simultaneous payment rule

At ideal order12 the family totals are A60/B62/C62/D62. The B11 affine
diagram halves its two fork-leaf/near-heavy-terminal pairs (8,10),(8,11);
the C11 affine path deletes (1,10). All three have distance9. The new
finite coupon is a half on the unique pair consisting of a whole A1 and
the penultimate node adjacent to the doubled terminal in a whole B11.
Thus its cost is exactly1/2 per A1/B11 component pair.

At ideal order16 the totals are A80/B102/C103/D103. The B15 affine diagram
deletes the three units (0,10),(10,14),(10,15); the D15 affine diagram
deletes the four halves (2,14),(2,15),(9,12),(9,13). Their distances are
at least10. The finite compensation upgrades the whole D5 unique long
endpoint against either short fork leaf of a whole D11 from a half to
a unit. The two eligible marked pairs add exactly1 per D5/D11 pair,
raising its total from15/2 to17/2. This does not change the separate
whole D5/A11 coupon or any finite D16/D17 same-component weight.

Finite distance7 pairs in whole A/B/D components receive zero whenever
both marks have feature depth at least8. Features are endpoints, branches
and endpoints of nonsimple bonds, as in the bound local theorem. The
previous threshold was10. The total numbers of deleted distance7 pairs
at rank n are

| Finite family | New depth8 count | Previous depth10 count |
| --- | --- | --- |
| A | max(0,n−23) | max(0,n−27) |
| B | max(0,n−24) | max(0,n−28) |
| D | max(0,n−25) | max(0,n−29) |

The binder subtracts only the difference of these counts from the
previous same-component costs. No coupon or ideal weight is changed by
this finite rule. In particular it affects none of the finite components
of ranks15–17 used by preceding compensation arguments.

## Exhaustive actual-corner composition

First consider an actual newly changed C12 corner. Its independent
universal near-terminal proof gives a distinct actual short1/2 unit.
Those units are unchanged by every present reduction, so this pays the
face without using its now-zero C12 pair.

Next consider an actual newly changed B12 corner. The independent B/C11
zero theorem has priority in its immediate whole B11 branch. Otherwise
the new B12 half is supplemented by a distinct retained D11 half, a
whole finite B12 unit, or the new whole A1/B11 half. The finite B12 unit
has marked distance10 and lies in a rank12 component, where none of the
current finite long-pair/depth rules changes it. The whole-component
proof allows arbitrary supporting roots: its Cauchy argument eliminates
the entire fork-side spectator vector, while the other root and its
arbitrary contacts are absent from the paying corner. The other new
ideal changes have orders16 and distances at least10, so none affects
these B12 outputs. The C12-first priority avoids recursion between the
two order12 changes.

After excluding every actual changed order12 corner, examine any actual
newly deleted finite distance7 corner. Its depth8 theorem proves either
an actual short6 alternative or a distinct retained whole A8 unit.
The latter has marked leaves. The former uses the independently
established current short6 interface: literal short units, the independent
cycle11 half transfer, or the universal cycle/tree replacement with
protected finite leaf/branch or whole small-component outputs. The
depth8 theorem explicitly checks these outputs, including the triangle
exception with its fixed marked leaf. None is a newly changed order12
distance9 or order16 distance10/11 pair; no new depth8 finite7 deletion
can remove a protected marked leaf or a whole A8 unit. This branch does
not call the selected large-affine exclusion back to establish its
short-path payment.

In the remaining case there is neither an actual new order12 corner
nor an actual new finite7 corner. On this whole face those two changes
agree pointwise with the preceding rule; their extra finite coupon can
only increase payment. Apply the order16 local theorem to the complete
preceding reduction, already including the order13 improvements. Its
independent short8 branch is established before its selected-factor
exclusions. The existing finite depth8 distance8 rule supplies short7
when necessary; the current finite7 rule supplies short6 or its retained
whole A8 unit. The earlier B/C11, cycle11, coupled order11, B12 center,
D13 and B/C13 priorities use only their stated smaller ideal roles or
separate finite components and remain independent.

In the no-short8 branch, the inherited low-core, high-disjoint,
protected/basic and adaptive routes are preserved. The order16 theorem
reconstructs all394 eligible patterns per family, including failed native
embeddings. Its changed B images are impossible by arbitrary-contact
physical projection or a saturated residual with no possible next white.
Its changed D images are excluded before imposing completion zeros or
seal a whole D5/D11 unit at a distinct corner. That new unit replaces the
deleted ideal half without consuming any finite same-component unit.
The complete routing and later interfaces are part of that local theorem,
not an inference from the current216 literal examples alone.

These cases exhaust actual faces. The order13 A4/B9 half and all older
compensations remain included. Triangle payment2 persists: the new ideal
pairs have distance at least9; every newly affected finite distance7
triangle output was checked in the depth8 theorem against its persistent
feature or fixed marked-leaf witness. Quadrilateral payment1 therefore
holds under one simultaneous rule with both new finite coupon costs.

## Complete current literal check

The binder reconstructs all216 inherited literal configurations from
their native Grams and starts with the final order13 cornerwise lower
payments. It maps every shared affine component of order12 or16 to the
new canonical targets and checks every finite corner for the D5/D11
upgrade. This complete literal finite domain is simply laced and its
largest whole finite component has order17. It therefore contains no
new whole B11 coupon and no newly deleted deep finite7 pair.

Exactly four configurations change: 1979,2372,2536,2540 in the inherited
D16 partial-weight package. Each loses one ideal half and gains the
distinct whole D5/D11 coupon half; each total stays1. There are no new
order12 debits and no B16 literal debit. All216 totals remain at least1,
and the four previously recorded D16 eighth credits are retained.
This exact saved-interface reconstruction supplements the universal
geometric composition above; it is not a replacement for it.

## Complete finite upper bound and actual attainment

Finite quantities in the producer use eighths. It checks simultaneous
same-component and cross-component dominance for all33,489 ordered type
pairs through rank60, including the exact new actual charges:

| Actual component pair | New actual total | Retained representative allowance |
| --- | --- | --- |
| D5/D11 | 17/2 | A5/D11:19/2 |
| A1/B11 | 1/2 | A1/D11:5/2 |
| A4/B9 | 1 | A4/B9:1 |

The last row is the preceding order13 charge, retained explicitly.
The first two new charges require no increase to the representative
interaction table. All earlier actual charges, including D5/A11 and
D5/D12, are still counted separately. Representative allowances are
upper bounds, not claims that their displayed representative forest
attains the actual value.

All1,044 fixed-component multiplicity vectors and792 interaction keys
are included. Fixed ranks are9 through17; the separate A17 representative
and residual D17 treatment are inherited unchanged. Residual ranks are
1 through8 and17 through60. Every actual forest is dominated by this
complete split, with both fixed-fixed and residual-residual cross costs
zero. The54,288 exact item identities verify every fixed interaction.

New residual tables are required because the large-rank same-component
items changed. Generation computes them by the item-first recurrence.
Replay checks all1,142,064 last-item upper inequalities and47,520 attaining
predecessors. Induction on residual rank gives both an upper bound and
attainment for each representative table, with no optimizer called in
replay. Twelve fixed-vector alternatives remain above405.

The largest representative value is409. The table displays
4A1+4A5+3D12 at this upper value, whose actual value is only358. This
does not make the bound unattainable: replacing the four A5 components
by D5 gives the actual formal forest

    4A1 + 4D5 + 3D12,

of total rank60 and exact cost409, checked directly with the full actual
cross rule. Thus409 is the exact maximum, not merely a conservative
upper estimate. A separate actual witness2D17+D26 now costs817/2=408.5;
it does not attain the maximum. No such formal forest is asserted to
occur as a vertex of a finite-volume hyperbolic Coxeter polytope.

## Exact ideal closure and its global consequence

The producer reconstructs the preceding exact integer convolutions and
scalar rows and verifies equality with the frozen order13 certificate.
The complete terminal domain has total rank59, at most one part from1
through8, and every other part at least9. Its2,113 profiles use the exact
family envelopes at every rank, with no constant high-rank substitution.
Only the established low9 concentration is used; arbitrary multiplicities
of ranks9 and above are retained.

For each row the new normalized cost saving is exactly

    count(rank11)/12 + count(rank15)/8.

These are total-cost savings1 at affine order12 and2 at affine order16,
divided by their true factor-fiber multiplicities12 and16. All ten
remaining negative profiles become positive; none is zero. The exact
normalized minimum is31/6000 at (9,9,9,9,9,14), and the exact absolute
minimum is54 at the singleton profile (59).

For an arbitrary actual cusp, inherited low9 concentration decreases
the normalized margin and does not increase the link-vertex product V.
Since every terminal normalized margin is now positive, multiplying
preserves the comparison: the actual absolute surplus is at least the
terminal absolute surplus, hence at least54. This positive-only product
argument was not available for earlier negative terminal rows. The
factor-fiber counting theorem converts these charges and surpluses to
the actual ideal vertices, rather than treating each factor as an
independent geometric vertex.

Write C(v) for each ideal vertex's charge and kappa(v) for its shape
curvature term. The simultaneous local payments and exact forest bound
give the global upper estimate

    D5 ≤ 409 v_f + sum_i C(i).

The inherited shape-curvature inequality in ambient60 gives

    D5 ≥ 405 v_f + sum_i kappa(i) + 4500.

Here the two extrema together contribute
5 floor(60²/4)=4500, counted once. Since kappa(i)−C(i)≥54, subtraction
gives4v_f≥4500+54v_i. In particular dimension60 requires at least1125
finite vertices, and at least1139 if an ideal vertex exists. An all-ideal
dimension60 polytope is impossible, and the completed unrestricted60
theorem excludes higher dimensions. Therefore the all-ideal class has
dimension at most59. The finite maximum409 exceeds405, so this argument
does not exclude dimension60 in general.

Restoring the older ideal17 center halves while reducing their finite
coupon is a separate possible successor. It is not applied in this
certificate; its finite count and simultaneous literal changes are not
silently assumed here.

## Verification scope

Generation session97211/chunkfa1358 and separate fresh replay
74572/chunk2b03bf both terminated exit0/PASS. Direct metadata checks
pass49 source bindings,6 local links and Python syntax. Root's full
producer/report mathematical composition read passed, including the
exhaustive actual-corner priority, exact representative domination and
actual attaining forest, simultaneous literal accounting and positive-only
absolute-surplus extension. Root independently checked the actual409
cost breakdown and a saved-margin update, but did not rerun the DP,
convolutions, literal computations or historical geometry. A second
independent full producer/report mathematical composition read passed:
it checked the same exhaustive priorities and protected outputs, both
new finite coupon allowances, the complete fixed/residual certificate,
actual attaining forest, exact ideal savings and positive-only product
comparison. That reviewer did not repeat geometric censuses, DP,
convolutions, literal or other arithmetic. Source, data and report are
final/frozen after these reviews and the direct metadata checks.
All bound local reports were frozen before this
generation. The producer checks source hashes, all finite type comparisons,
new complete Bellman certificates, reconstructed ideal convolutions and
all216 literal interfaces. The local arbitrary-support, all-label,
classification, concentration, factor-fiber and shape-curvature proofs
remain explicit inherited inputs. No historical geometric census is
rerun, and no external refereeing or proof-assistant validation is claimed.

    python3 scripts/fresh_order12_order16_depth8_ambient60_frontier_20261018.py --verify
