# Finite first-extension stopping census for the global dimension-21 proposal

Generation and a fresh saved replay pass. The [source](../scripts/fresh_sqrt5_global_b1_finite_first_20261016.py) and [exact data](fresh_sqrt5_global_b1_finite_first_20261016.json) classify all 2,464 finite first extensions of the potentially affected baseline-one minimal quadrilateral rows. Opposite-pair-preserving black permutations and all white permutations reduce them to 1,646 marked classes; every input has a checked explicit isomorphism map.

There are 196 admissible singular stopping classes and 1,450 nonsingular open classes. Each classification checks all four spherical corners, the exact field gauge, physical and conjugate inertia, integral cyclic products, and the signs of the root relation. The singular classes have exactly one relation, one physical negative direction and no conjugate negative direction. This is a complete first-extension classification, not a complete white-extension closure.

The payment scan sorts the singular classes by node count and index and stops at its first underpayment. Class 39, an eight-node extension of minimal class 3249, is the first failure. Its four corner forests are A1²+D4, A1²+A4, D4+I2(4), and B6. The B6 corner has the single baseline payment and one eligible adjacent-edge discount. With eight external A1 components and two I2(5) components, its minimum payment is approximately 0.77036334; the remaining three spectator ranks can be realized by A3. The exact profile is

`1 − delta_adjacent3 + (w14894 + w21728)/10`.

Both displayed compensation coefficients vanish in the [frozen restricted global proposal](fresh_sqrt5_global_discount_20261016.md). Only the payment scan stops early: all 196 singular classes and all 1,450 nonsingular classes are retained in the data. A [separate exact geometric witness](fresh_sqrt5_global_b6_discount_witness_20261016.md) checks that the first failure really supports a compact local quadrilateral.

The two compensation channels are A1 crossed with the A4 internal-root orbit and A1 crossed with a D4 leaf. Both use actual A1 count 10, capped I2(5) context 2 and no H3. The parameter metadata's A1 index 11 is the inherited shifted index, not an extra root.

No second or later attachment is enumerated here. An admissible singular Gram is not by itself asserted to realize a compact ambient polytope. The result identifies a missing constraint in the global proposal; it neither proves a new dimension bound nor affects the separate frozen compact arithmetic dimension-21 proof candidate.
