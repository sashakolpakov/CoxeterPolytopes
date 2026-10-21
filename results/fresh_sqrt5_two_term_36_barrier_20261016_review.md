# Independent literal-graph review of the two-term barrier at 36

**PASS.** The [standard-library reviewer](../scripts/fresh_sqrt5_two_term_36_barrier_20261016_review.py)
and [certificate](fresh_sqrt5_two_term_36_barrier_20261016_review.json)
independently check the literal graphs and two exact inequalities in the
[restricted payment-model obstruction](fresh_sqrt5_marked_h4_remainder_free_bonus_20261016.md).
Generation and fresh saved replay pass. No producer code is imported.

The ten-node quadrilateral is reconstructed from all 45 pair labels.
Its four corners are A4+A4, A3+D5, A3+D5 and D4+D4. Explicit graph
isomorphisms identify the two marked D5 pairs with the same distance-three
orbit. With common external components consisting of a marked A1 and A12,
both payments are the same variable w. The other two corners receive zero
in this specified family, so its face obligation is 2w>=1.

A separate literal D5 recount finds eight baseline-paid pairs and two
distance-three pairs in that orbit. At the marked rank-21 vertex
A1+4D5 the weight is therefore 32+8w>=36. The exact dual multipliers 4
and 1 cancel w and reproduce the same lower bound. No H4 context change
affects this witness.

This checks an obstruction to the stated payment family, assuming its
baseline, discounts applying only from rank six, and absence of the
A4/D4 cross payments. It does not prove global geometric realizability of
the quadrilateral, audit the full dynamic program or all face coverage,
or exclude compact arithmetic dimension 21. Adding cross payments or
changing the D5 baseline can escape this obstruction.

```sh
python3 scripts/fresh_sqrt5_two_term_36_barrier_20261016_review.py --verify
```
