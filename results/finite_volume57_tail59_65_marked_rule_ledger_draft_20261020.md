# Exact marked rules for the ambient59–65 tail

Staging ledger,20 October2026. No imports, generators, scientific checks or
builds were run. This supplements the ambient66–97 ledger of batch_guard
and the intrinsic ambient58 appendix. It defines actual marked payments;
representative cost allowances are not used to define any corner weight.

Write `R_N` for the rule used to exclude ambient dimension N. The canonical
finite and affine labels are those of the intrinsic appendix: `A_r` path,
`B_r` path with doubled last bond, `D_r` spine0,...,r-3 with short leaves
r-2,r-1; affine order m means m nodes and rank m-1. An A endpoint set is
`E_r={0,r-1}` and `P_k(A_r)={k-1,r-k}`. A D long end is0. All components
in every cross-role or guard are **whole components** at the actual corner.

Each layer below is a mathematical assignment. A patch **sets** its stated
weight; it is not a debit subtracted again from a previously changed entry.
Other entries inherit the preceding layer. Thus these seven finite rules
are completely evaluable without running a scientific producer or reversing
the stronger, dimension58-specific rule. The layer order65→...→59 is a
definition of a finite table of functions, not the exposition of the
geometric proof.

## A. Base rule R66 and exact supplement lookup semantics

The complete R66 definition is in Section8.3–8.6 of
`finite_volume57_higher_tail_algorithmic_appendix_draft_20261020.md`
(held SHA256 `d8f775d30b7eba4e181c7083b0075b689248108c42ae7818123f9c8b5e9bb4bf`).
The finite part needed here can also be specified directly.

For rank at most16 every same-component pair is one. For rank at least18
keep the basic set: distances at most8, and distance9 at a whole-component
leaf. For `D_r,r>=20`, halve the distance9 short-fork-leaf pairs. For rank17
use the following positive sets (all unlisted pairs zero):

* A17: the basic set just stated, the ten pairs
  `(2,11),(2,12),(2,13),(2,14),(3,12),(3,13),(3,14),(4,13),(4,14),(5,14)`,
  and `(0,13),(3,16)`; every retained weight is one (cost114).
* B17: the basic set and `(2,11),(3,12),(4,13)` (cost105).
* D17: all pairs except
  `(0,10),(0,11),(1,10),(1,11),(1,12),(1,13),(1,14),(1,15),(1,16)`
  (cost127).

Its ordinary cross table is the broad half on every mark in a component
of rank<=5 against the center/short leaves of D11–D16, the middle of
A13/A15 or nodes5,11 of A17; also the half against D10 short leaves.
Additional roles are: A4-end/D13-node9 half; A2/B10 doubled terminal
half; A1/D11 and A1/D12 short-leaf upgrades to one; A2/D11 short-leaf
upgrades to one; A1/A11 middle one; A2/A11 middle half. Omitted roles
zero. There are no context guards. This is the literal table; replacing
B17 by D17 or E6 by a representative in a budget does not change it.

For R65 the affine-order12–17 base is read **directly**, with no inferred
cost-to-pair reconstruction, from

    fresh_deep_finite8_dimension64_20261018.json
      .short7_small_factor_interfaces

There is exactly one record for each `(order,kind)` in
`{12,...,17}×{A,B,C,D}`. Its `target` is the whole labelled graph, and
`weighted_pairs` is the weight table. The four arrays at each order have
the corroborating totals

|order|A|B|C|D|
|---:|---:|---:|---:|---:|
|12|60|64|63|62|
|13|65|73|73|74|
|14|70|84|84|82|
|15|75|94|91|94|
|16|80|107|103|106|
|17|85|111|108|119|

Every affine pair at orders<=11 is one for R65. Orders>=18 use distance
at most5 or distance6–9 with a mark within5 of a branch/nonsimple-bond
endpoint. This replaces the somewhat older R66 order15 C/D tables; it
is the explicit table adopted by the R65 composition, not an assumption
that the R66 envelope already had those values.

**Finite table schema.** A saved weighted-pair record is
`{pair:[i,j],weight:p/q}` with `0<=i<j<rank/order`. Fractions are exact.
An array asserted to be a complete table defines zero on absent pairs;
it may explicitly contain zero entries. An array named `retained_pairs`
defines an indicator of that set. A list named `changed_pairs` or
`half_pairs` is **only a patch**, not a complete table. Each locator below
specifies which of these interpretations is intended. Record selection
is by its type and node count (`kind`,`rank`, or length of `target`),
never by its position in an incidental JSON array. All target labels
must agree with the canonical labels, including nonedges.

## B. The all-rank depth table

A finite feature is an endpoint, branch or endpoint of a nonsimple bond
in the **whole** component; let delta be distance to that set. Beyond
the finite small-rank replacements in SectionC, the large-rank rule is
the R66 basic rule with these complete depth changes:

|N|distance6|distance7|distance8|
|---:|---|---|---|
|65|unchanged|unchanged|set1/2 if both deltas>=17|
|64,63,62|unchanged|unchanged|set0 if both deltas>=11|
|61|unchanged|unchanged|set0 if both deltas>=8|
|60|unchanged|set0 if both deltas>=6|set0 if both deltas>=8|
|59|set0 if both deltas>=6|set0 if both deltas>=6|set0 if both deltas>=8|

In particular **R59 is not the final R58 depth-six distance8 rule**.
No depth change affects the rank12–17 table replacements below.

## C. Finite same-component assignments

Start with the finite R66 rule in SectionA and apply the following layer
assignments. Use SectionD for all cross roles and SectionE for contexts.

### R65

Set the following units to zero:

```text
A16:(0,14),(0,15),(1,14),(1,15).
B16:(0,11),(0,12),(0,13),(0,14),(0,15),
    (1,12),(1,13),(1,14),(1,15),(2,13),(2,14),(2,15),
    (3,14),(3,15),(4,15),(5,15).
```

Set to one half:

    D16:(1,12),(1,13),(1,14),(1,15);
    D17:(2,13),(2,14),(2,15),(2,16).

The costs are A16/B16/D16=116/104/118 and A17/B17/D17=114/105/125.
These are exactly the tables in
`fresh_finite_d16_d17_long_pair_half_transfer_20261017.json`:
`A16_B16_records` uses **retained_pairs indicators**; `records`, selected
by rank16 or17, uses **weighted_pairs complete tables**. The D17 array
omits old zero entries, so absent means zero rather than one.

### R64

Set D16 `(0,11),(0,12)` to zero; restore D17 `(1,12),(1,13)` to one.
Other small finite entries remain R65. The costs become D16=116,D17=127.
Complete tables are `fresh_finite_d16_d17_two_unit_exchange_20261018.json`
`.rules`, selected by rank16/17, `.weighted_pairs`.

### R63

Set to zero:

    A16:(0,12),(3,15);
    D16:(0,10),(0,13);
    D17:(1,12),(1,13).

The full tables are
`fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018.json.rules`,
selected by kind and node count, `.weighted_pairs`.
The costs are A16/D16/D17=114/114/125. B16 remains104; rank15 remains
all units. The six sparse cross additions are in SectionD.

### R62

Set the following entries:

```text
A15:(0,13),(1,14) ->1/2; B15:(0,14)->0;
D15:(0,13),(0,14) ->1/2;
A16:(0,11),(4,15) ->1/2;
D16:(2,13),(2,14),(2,15) ->1/2;
A17:(1,12),(4,15) ->1/2;
D17:(0,12),(0,13) ->0;
    (0,14),(0,15),(0,16),(4,14),(4,15),(4,16) ->1/2.
```

The rank15 complete tables are
`fresh_finite_rank15_existing_half_transfer_20261018.json.rules[].weighted_pairs`.
The rank16 and A17 assignments are exactly the named patch lists
`A16_half_pairs,D16_half_pairs,A17_restored_half_pairs` in
`fresh_rank16_existing_payment_exchange_20261018.json`.
The full D17 table is
`fresh_finite_d17_five_unit_reduction_20261018.json.weighted_pairs`.
The costs are A/B/D15=104; A16=113,D16=225/2,A17=115,D17=120.

### R61

Replace the indicated finite tables by these complete saved tables:

|component|JSON stem (all have suffix `_20261018.json`)|field|
|---|---|---|
|A15|fresh_finite_a15_endpoint_coupon_zero|weighted_pairs|
|B15,D15|fresh_finite_b15_d15_existing_payment_transfer|rules selected by kind;weighted_pairs|
|A16|fresh_finite_a16_endpoint_sparse_half_transfer|weighted_pairs|
|D16,D17|fresh_finite_d16_fork_exchange|D16_weighted_pairs,D17_weighted_pairs|

A17 stays115. The other costs are A/B/D15=103; A16/D16=111;
B16=104; D17=119. For clarity the additional D17 change on the path to
this table halves `(3,13),(3,14),(3,15),(3,16)` and restores
`(1,15),(1,16)` to halves. It does not delete their final weights.
The new A16 halves are `(0,13),(2,15),(1,13),(2,14)`; the new D16
halves are `(0,14),(0,15),(2,12)`. These extra descriptions corroborate
the complete tables; they do not replace any inherited zeros or halves.

### R60

Replace A15/B15 by the complete kind-selected tables
`fresh_rank15_sparse_payment_exchange_20261018.json.rules[].weighted_pairs`.
Replace D15 by
`fresh_finite_d15_three_singleton_half_transfer_20261018.json.weighted_pairs`.
Their costs are101,97,201/2. There is no rank14 modification yet.

For A16/A17 start from the **complete** arrays `A16_weighted_pairs` and
`A17_weighted_pairs` in
`fresh_a16_a17_full_payment_exchange_20261018.json`, and set

    A16:(0,11),(4,15) ->1/8;
    A17:(1,12),(4,15) ->7/8.

For D16/D17 start from the **complete** arrays `D16_weighted_pairs` and
`D17_weighted_pairs` in `fresh_d16_d17_eighth_balance_20261018.json`, and set

    D16:(0,14),(0,15) ->7/8;
    D17:(1,15),(1,16) ->1/8.

These are the exact assignments of the center-reversion balance, not
the full arrays' earlier5/8 and3/8 values. The resulting costs are
A16=441/4,D16=443/4,A17=D17=467/4. Cross-D16 payments to ranks2–5
are removed separately in SectionD.

### R59

Replace B14/D14 by the kind-selected complete tables in
`fresh_finite_rank14_existing_payment_transfer_20261018.json.rules[].weighted_pairs`.
Replace A14 by
`fresh_finite_a14_near_endpoint_zero_transfer_20261018.json.weighted_pairs`.
Their costs are89,87,181/2: A14 zeros `(0,12),(1,13)`, B14 zeros
`(0,13),(1,13),(2,13),(3,13)`, D14 half `(0,11)`.

Replace A15 by
`fresh_finite_a15_three_coupon_zero_transfer_20261018.json.weighted_pairs`
and D15 by
`fresh_finite_d15_a4_a11_half_zero_transfer_20261018.json.weighted_pairs`.
Both cost100; B15 stays97. Replace D16 by
`fresh_finite_d16_alternate_half_zero_transfer_20261018.json.D16_weighted_pairs`
(cost439/4). A16 stays441/4. Replace A17 by
`fresh_a17_five_pair_alternate_transfer_20261018.json.weighted_pairs`
(cost451/4). Starting from R60's D17 table, make the disjoint assignments

    (2,12)->5/8; (2,14)->0; (3,15),(3,16)->0.

The D17 cost is919/8. The first two assignments are the exact patch
in `fresh_d17_existing_three_corner_transfer_20261018`; the last two
are `fresh_finite_d17_leaf_a16_half_transfer_20261018`. Each source also
has a complete table for its separately applied patch; those two separate
tables must **not** overwrite each other. D12,A13,B13,D13 remain all units.
No final R58 same-component presence guards are used in R59.

## D. Ordinary cross-role table by layer

Begin with R66's literal cross table in SectionA. At each following layer
retain all earlier cross roles, then **set** the rows below. Roles omitted
from the entire accumulated table have weight zero. The sets `E,P_k` are
defined above; the D short-leaf set has two nodes. A small component of
rank<=5 includes its non-A types and the full unbounded I2 family whenever
the broad R66 row names it.

|N|small marked role|large marked role|new weight|
|---:|---|---|---:|
|65|A1 or A2,any mark|D14 node10|1/2|
|64|E6 long ends|D11 short leaves|1/2|
|64|A1|B10 doubled terminal9|1/2|
|63|A4 endpoints|A12 endpoints|1|
|63|A4 endpoints|A13 P2|1|
|63|D5 long end0|A12 P2|1|
|63|A2 any mark|A13 endpoints|1/2|
|63|A2 any mark|A14 P2|1/2|
|63|E6 long ends|A11 P2|1|
|63|A2 any mark|B9 doubled terminal8|1/2|
|63|A2 any mark|D10 center7|1/2|
|63|A3 middle1|D10 center7|1/2|
|62|A2 any mark|D10 center7|1|
|62|A3 middle1|D10 center7|1|
|62|E6 long ends|D11 short leaves|1|
|61|A2 any mark|A13 endpoints|1|
|61|A1|A13 endpoints|1/2|
|61|A2 any mark|A14 P3|1/2|
|61|A4 endpoints|A12 P3|1/2|
|61|A1|B10 doubled terminal9|1|
|60|A2 any mark|D13 short leaves|1|
|60|A4 endpoints|D11 short leaves|1|
|60|A5 P2|D11 short leaves|1|
|60|D5 short leaves3,4|D12 center9 and leaves10,11|1|
|60|D5 long end0|D12 short leaves10,11|7/8|
|60|D5 long end0|A11 P2|1|
|60|A4 endpoints|B9 doubled terminal8|1/2|
|60|D5 long end0|D11 short leaves9,10|1|
|60|A1|B11 node9|1/2|
|60|every mark of any rank2–5 component|D16 center13 and leaves14,15|0|
|59|A2 any mark|A13 middle6|7/16|
|59|A2 any mark|D13 center10|7/16|
|59|A1|A13 middle6|9/16|
|59|A1|D13 node9|1/2|
|59|D5 long end0|D12 short leaves10,11|1/2|
|59|D5 long end0|A12 endpoints|3/8|
|59|A1|D5 long end0|3/8|
|59|D5 short leaves3,4|A12 P3|1|
|59|A1|A14 P3 and P4|1/2|
|59|A1|A15 P3|1/2|
|59|A1|A16 P4|1/2|
|59|A1|D14 node9|1/16|
|59|A1|A12 middle nodes5,6|1|
|59|A1|B12 node9|1/2|
|59|A1|D15 node11|1/16|
|59|A1|A12 endpoints|1/4|
|59|A1|A13 P2|1/4|
|59|A4 internal nodes1,2|A13 P2|1/2|
|59|A4 endpoints|A11 endpoints|1/2|

The D5-long/D12-center value remains its broad1/2 in R60 and R59;
its temporary3/4 is not adopted. The R60 D5-long/D12-leaf value7/8
is **replaced** by1/2 at R59, while the D5-short/D12 fork upgrades stay1.
Accordingly the full D5/D12 charges are45/4 at R60 and21/2 at R59.
The D5/A11 cost is2, D5/D11 cost17/2, and R59 D5/A12 cost27/4.
These marked distinctions are essential; a single type-pair total does
not say which long-end or short-leaf unit a proof may spend.

The R60 first three rows above are the actual upgrades in
`fresh_rank15_order15_ambient60_frontier_20261018`: they are not the
earlier conservative A2/D13 or A5/D11 allowances. R61 had allowed the
larger representative A2/D13 cost to dominate A2/A13, but its actual
A2/D13 short-leaf weights were still halves. R60 makes them units.
Similarly the R63 budget's E6/D11 allowance4 is not its actual charge2;
that geometric upgrade occurs at R62.

The R59 A1/D5 entry is unconditional. The final R58 rule's A11 or
A9+A2 guard, and its additional5/8 credit, are **not** tail rules.
No A1/A16 endpoint quarter is present in these layers. The A4/A12
third-node weight remains1/2 at R59, not its later R58 unit.

## E. Context guards: absent until R59

There are no context-dependent additions in R65,R64,R63,R62,R61 or R60.
At R59 add exactly the following two roles, and nothing else:

* For marks in two distinct whole A2 components, weight1/8 if an unmarked
  whole A10 and at least two unmarked whole A1 components are present.
* For any A2 mark versus the D5 long endpoint0, weight1/2 if unmarked
  whole A1 and A9 components are present.

Multiple witnesses do not multiply either payment. These ordinary roles
were zero before their guarded addition. Their total forest costs are

    (1/4) a2(a2-1)[a1>=2,a10>0]
      + a2 d5[a1>0,a9>0].

This is not the final ambient58 context polynomial. In particular the
first marked-pair value is1/8, not9/32, and there are no same-component
A13/A14/A15/D14/D15/A16/D16 restorations or rank15-presence guards here.

## F. Exact affine small-order assignments

Use R65's full order12–17 snapshot from SectionA. For each later layer,
keep all unchanged entries and apply the following replacements. All
orders10 remain entirely unit throughout59–65. All orders<=9 are also
entirely unit. Orders>=18 retain the common all-rank distance5/feature5
formula. Thus the following finite table is the whole ideal-rule ledger.

### R64

At order11 set A-cycle distance5 pairs to1/2; B pairs `(8,9),(8,10)`
to1/2; C pair `(0,10)` to1/2; D pairs `{7,8}×{9,10}` to1/2.
At D,order17 set `(3,15),(3,16),(9,13),(9,14)` to1/2.
These are the exact tables in `fresh_cycle11_half_pair_transfer`,
`fresh_order11_bcd_half_pair_transfer.factor_rules[].weighted_pairs`
and `fresh_order17_e6_half_coupon_transfer.weighted_pairs`, respectively.
All stems here and below have suffix `_20261018.json` unless17 is stated.

### R63

Set C,order11 `(0,10)` to zero; B,order12 `(0,9)` to zero;
D,order13 `(0,8)` to1/2. The complete C/D weighted tables are in
`fresh_order11_c_endpoint_zero_transfer` and
`fresh_order13_d_center_pair_half_transfer`. For B12 use the set
`fresh_order12_b_center_terminal_transfer.new_paid_pairs` as a unit
indicator table (not a `weighted_pairs` field that this file lacks).

### R62

Set D,order13 `(0,8)` to zero. Set B,order16 `(1,11),(2,12)` to zero.
Set D,order17 `(3,15),(3,16),(9,13),(9,14)` to zero. At order14 set
B `(9,12),(9,13)` and C `(1,12)` to zero. Complete tables are, respectively,

    fresh_order13_d_center_pair_zero_transfer.weighted_pairs;
    fresh_order16_affine_b_singleton_projection_transfer.weighted_pairs;
    fresh_order17_e6_unit_coupon_zero_transfer.weighted_pairs;
    fresh_order14_b_terminal_residual_transfer.weighted_pairs;
    fresh_order14_c_inner_endpoint_projection_transfer.weighted_pairs.

### R61

At B,order11 set `(8,9),(8,10)` to zero. At C,order11 set
`(0,9),(1,10)` to zero (its endpoint was already zero). At D,order16
set `(0,9),(2,11)` to1/2. These complete tables are
`fresh_order11_bc_zero_transfer.rules[].weighted_pairs` and
`fresh_order16_d_distance9_half_transfer.weighted_pairs`.

### R60

Replace the following whole affine tables by these complete snapshots:

|kind/order|JSON stem|field|
|---|---|---|
|B12,C12|fresh_order12_bc_near_terminal_half_transfer|rules by kind;weighted_pairs|
|B13,C13|fresh_order13_bc_sparse_terminal_transfer|records by kind;weighted_pairs|
|D13|fresh_order13_d_center_leaf_half_transfer|weighted_pairs,absent pairs zero|
|C14|fresh_order14_c_distance10_terminal_projection|weighted_pairs|
|B15,D15|fresh_order15_bd_terminal_coupon_transfer|rules by kind;weighted_pairs|
|B16,D16|fresh_order16_bd_sparse_unit_transfer|records by kind;weighted_pairs|
|D17|fresh_small_d16_coupon_and_ideal17_rebalance|local_proofs.weighted_pairs|

The last table restores the two center pairs `(0,10),(2,12)` to halves
and the four distance11 pairs `(2,15),(2,16),(10,13),(10,14)` to1/8.
It is the rebalance table, not the superseded quarter/three-quarter rule.

### R59

Replace these tables by the indicated complete snapshots:

|kind/order|JSON stem|field|
|---|---|---|
|B13,C13|fresh_order13_bc_inner_terminal_transfer|rules by kind;weighted_pairs|
|B14|fresh_order14_b_basic_distance9_transfer|weighted_pairs|
|C14|fresh_order14_c_basic_distance9_transfer|weighted_pairs|
|D14|fresh_order14_d_near_fork_half_transfer|weighted_pairs|
|B15,D15|fresh_order15_bd_basic_fractional_transfer|rules by kind;weighted_pairs|
|C16,D16|fresh_order16_cd_sparse_sixteenth_transfer|rules by kind;weighted_pairs|
|D17|fresh_order17_center_contextual_zero_transfer|local_proof.weighted_pairs|

Equivalently, the R59 ideal rule is exactly Section6 of the intrinsic
ambient58 appendix **except that every pair of every order10 factor
has weight one**. This equivalent description is a check against the
direct saved snapshots, not the definition used to recover them.

The following family maxima corroborate the tables; they do not define
the marked entries:

|N|m10|m11|m12|m13|m14|m15|m16|m17|
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|65|45|55|64|74|84|94|107|119|
|64|45|109/2|64|74|84|94|107|118|
|63|45|54|63|147/2|84|94|107|118|
|62|45|54|63|73|83|94|106|116|
|61|45|53|63|73|83|94|105|116|
|60|45|53|62|145/2|82|92|103|229/2|
|59|45|53|62|143/2|80|367/4|411/4|227/2|

## G. Why these are mathematical rule definitions

A reader evaluates any actual marked corner by selecting its ambient
layer N, its whole component type/labels and its two marked roles. The
finite large-rank formula, the complete small-rank table replacements,
the accumulated ordinary cross rows and the two explicitly Boolean guards
cover every finite type. The complete affine snapshots and common
all-rank formula cover every ideal type. A finite weighted-pair table
has a specified finite domain and exact default; no scalar maximizer,
upper representative or contextual budget relaxation is used as a lookup.
Thus the same function is the input to the geometric payment theorem
and the rank-N finite/cusp algorithms.

All finite snapshot files remain immutable. Each comes with its complete
target graph and pair domain or its explicitly bounded patch. Its SHA256
is specified by the corresponding frozen composition's `sources` mapping,
with transitive source bindings where the composition imports a local
theorem rather than the table itself. The exact global composition stems,
indexed by **excluded ambient dimension**, are:

```text
65 fresh_deep_finite8_dimension64_20261018
64 fresh_order11_order17_dimension63_20261018
63 fresh_sparse_coupon_order11_13_dimension62_20261018
62 fresh_rank15_rank16_dimension61_20261018
61 fresh_order11_order16_dimension60_20261018
60 fresh_depth6_dimension59_20261018
59 fresh_a15_d15_dimension58_20261018
```

This source manifest is provenance for the definitions just given, not
their substitute. In particular a successful finite budget does not prove
that an abstract embedding is an actual whole component. The local
payment proof for every layer still uses the actual positive/affine
projection lemma, complete core/overlap domain, arbitrary-support seals,
and independent actual-corner priorities. A new affine zero is handled
by its complete whole-face theorem first; otherwise old ideal weights
agree pointwise on the actual face, allowing the finite sparse theorem
to transfer. Deep finite corridors use the independent shorter-path
routes and protect every feature/leaf output. Finite half/zero exchanges
retain their stated distinct whole-component alternatives. Boolean guards
are proved by whole-component physical seals before their costs enter DP.
These are the distributed routing clauses from the current appendix,
with each layer's own numeric table; one cannot substitute the stronger
ambient58 restorations or coupons into a saved tail counting row.

The ledger is a source-level mathematical transcription. It does not
claim a new execution, new pair-table summation, new domain census or new
automated end-to-end verifier. The recorded local and composition results
belong to their respective frozen sources and certificates.
