# Rank 17 spectator extension: exact finite maximum 371.5

This separate scout extends the reviewed D16(3,13) deletion to a separate
whole A17 or D17 witness. The separate whole D_s, s >= 18, guard is retained.
The same actual-core argument proves the extension; the complete finite cost
is 743/2 = 371.5, attained by 2A17 + D23. Curvature remains 364. This does not
prove a dimension 56 bound. Every predecessor remains unchanged. The new rule
is a scout pending independent review.

Source: `scripts/fresh_d16_rank17_witness_extension_scout_20261018.py`.
Certificate: `results/fresh_d16_rank17_witness_extension_scout_20261018.json`.
Independent check: `scripts/check_d16_rank17_witness_extension_independent_20261018.py`.

## Exact extension and inherited actual geometry

Write r16 = a16 + d16, r17 = a17 + d17, and g for the presence of a whole
D_s with s >= 18. Counts always refer to actual whole components. The final
expanded rule deletes the marked whole D16(3,13) unit when g holds and a
separate whole A16, D16, A17 or D17 is present. At a different corner it adds
a unit to each A2/A13-fourth pair, with an unmarked whole A1 and the same
external witness and large D. This is four marked pairs, costing 4 per
A2/A13 component pair. Witness presence is Boolean. Existing A13/D13 coupons
and every other finite or ideal weight are retained.

Relative to the frozen completed D16 rule, the exact ordinary-cost increment
is

    g*[r17>0]*(-d16*[r16<2] +4*a2*a13*[a1>0]*[r16=0]).

The selected pair and its complete local domain are identical to those in
`fresh_d16_node3_rank16_large_d_guard_transfer_20261018.md`: 394 patterns,
578 maps, precisely 24 affected maps in fourteen contexts. Its six full and
eight proper core arguments do not assume the spectator has rank 16. They
prove whole-core orthogonality for every additional supporting white, using
the complete contact domain, radical relations, or the nonsingular inverse
projection inequalities. The two genuine proper completions are included.

In each surviving selected corner the displayed marked D16 uses sixteen
roots, and the displayed remainder has at most one A1. Every extra component
lies in the same orthogonal forest E at all four corners. A separate whole
A17 or D17 cannot lie in that remainder, and cannot straddle the core and E:
the arbitrary-support seal forbids its first connecting root. It therefore
lies wholly in E, as does the separate D_s >= 18. Both actual components
persist at the distinct whole A1+A2+A13 coupon corner. This proves the new
guard there; it is not an assumption about a literal or minimally displayed
forest.

The debit orbit is unchanged, so the completed output-retention audit applies
without a new affected role: no current literal, later consumer, newer guarded
seal, chosen-D12 output, protected/adaptive output, or ideal-priority output
calls the changed D16 pair. Its 216 literal configurations, 564 later corners,
newer core seals and the two chosen-D12 lists remain the frozen inputs.
Only additional spectator cases of the already proved selected route become
active. Short8 escape payments and all ideal payments are unchanged. Thus
the existing universal triangle/quadrilateral composition proof extends with
the guard-persistence argument above.

## Complete new positive family

A newly positive forest has whole A13, A1, A2, one A17 or D17 witness, and
D_s >= 18. These consume at least 51 ranks. Hence there is exactly one
witness and one A13, no rank16 component or other fixed rank9--17 factor,
and at most one rank5 factor. The residual rank is at most 6. This also
proves that its positive family is disjoint from the old rank16-coupon
family at rank57.

Use the frozen representative class with A13 and A17. The D17 witness is
covered separately by an extra base allowance of 1/4: its old ordinary base
is 113.5, while the A17 representative allowance is 113.25. Every actual
D17 cross charge with A13 or any rank <= 6 factor is exactly zero, checked
from the actual type function; the corresponding A17 representative charges
are nonnegative. Thus this pad dominates both witness choices. Negative
D17, D16 and A16 modifications are ignored.

The frozen representative item weights gain 4 at rank2 for the new coupon,
and 2 at rank4 for the already present guarded-D17 A4/A13 allowance. No
rank16 factor can occur, so no coupled-A16 positive allowance is needed.
Reserve all actual A1 and A2 multiplicities explicitly. For rank5 count
k=0,1, reserve D_s, where 18 <= s <= 24-5k. The remaining rank is
27-5k-s. Two small unbounded-item tables and all nine reserved-D rows give

    new positive-family cost <= 22917/64 = 358.078125 < 364.

The representative relaxation need not be attained by an actual forest.

## Complete finite cover

The original frozen rank57 certificate is 133MB. The producer streams only
its compact 38-class above-curvature summary and streams all hashes. Its
unlisted classes are at most 364. Exactly five listed classes have an old
maximum greater than 371.5. They contain no A13, no rank5 factor, and no
fixed ranks9--15:

| Fixed A17 count | Fixed rank16 count | Complete new maximum |
|---:|---:|---:|
| 0 | 0 | at most 370, inherited completed refinement |
| 0 | 1 | at most 371, inherited completed refinement |
| 0 | 2 | at most 371, inherited completed refinement |
| 1 | 0 | 1483/4 = 370.75, fresh refinement |
| 1 | 1 | 743/2 = 371.5, fresh refinement |

Here the fixed rank17 representative is **A17 only**:
`representative(A17)=A17`, whereas `representative(D17)=D17`. D17 is a
residual item, including in the fresh classes. The rank16 representative
fiber consists of all actual A16, B16 and D16 types. The fresh refinement
fixes A17 and each zero/one-factor choice in this complete rank16 fiber.
All other residual actual types whose representatives are frozen items are
included, with A1 multiplicity fixed explicitly.

For each A1 count, split into no large D and one reserved D_s >= 18. The
guarded table permits another residual large D whenever rank allows it.
All residual actual pair cross charges vanish, as explicitly asserted in
the source. After the fixed tests, the remaining context is additive in
residual components: A13 and all rank9--15 coupon factors are absent; there
cannot be two A16 factors; rank16 multiplicity is fixed; the A17 witness is
fixed; and large-D presence and A1 multiplicity are fixed. The old guarded
D17 decrement is consequently linear in residual D17 count. The new D16
debit is a fixed constant, and the new A2/A13 coupon is identically zero.
There is no residual-dependent presence test left. This algebraic argument
proves multiple-component additivity; checking single-item weights alone
would not suffice. Single-item equality for every permitted reserved D_s
is also checked, so a D18 table covers all actual large-D choices exactly.

Fresh rational Bellman inequalities, attaining predecessors and explicit
guard options certify both new class maxima. In the previous three classes
the new coupon is zero, so their frozen refined maxima can only decrease.

The five chronological modification families have complete upper bounds

    order17: 739/2; guarded D17: 365; coupled A16: 371;
    rank16-witness D16: 11513/32; rank17 extension: 22917/64.

All are below 371.5. Classify a forest by the last modification at which its
cost increases. That family's bound includes all earlier changes, and all
later steps for that forest are nonincreasing. If it never increases, it
decreases monotonically from the original rule and is covered by the
original classes plus the five refinements. This covers every finite
forest, including those absent from the original above-curvature summary.
Finally the actual formal forest 2A17+D23 has cost exactly 743/2 under the
new rule, proving the stated complete maximum. No realizability is inferred.

## Verification scope

The producer checks the new guard formula, the complete new-positive family,
the two new top classes, all actual rank16 fibers and reserved-D choices,
the inherited five-class cover, and source hashes. The independent checker
imports the frozen rule checker, not this producer, reconstructs the extended
cost directly, checks the saved Bellman certificates and cover arithmetic,
and checks residual-context additivity independently. Universal geometry is
the written extension of the already reviewed full-core theorem above.
Neither an ideal-profile recomputation nor a new dimension bound is claimed.

Generation `134636`, fresh replay `60026d`, and independent check `335fac`
all passed. The independent check verifies seven streamed hashes, 2472
coupon Bellman inequalities and 114 predecessors, 36579 top-class Bellman
inequalities and 1884 predecessors, 432 explicit guard options, and 34325
feasible residual pair-additivity identities. All processes exited; the
shared lock was released. No large certificate was loaded wholesale.
