# Independent review of the golden balanced-cap construction

Reviewed the [producer](../scripts/fresh_h9_symmetric_balanced_caps_20261013.py), [standalone Fraction verifier](../scripts/fresh_h9_symmetric_balanced_caps_verify_20261013.py), and [literal input](../data/enumerations/fresh_h9_symmetric_balanced_caps_20261013.json). Independently reran the verifier: **PASS**, including rejection of a forged cap. No mathematical issue found.

The seed signature calculation is complete: the seven standard symmetric leaf directions are null, the seven standard antisymmetric directions have eigenvalue 4, the remaining symmetric two-by-two block has rank one and negative trace, and the antisymmetric block is positive definite. Thus the 18-by-18 Gram matrix has signature `(9,1,8)`. Conjugation gives `(10,0,8)`. The exact strict timelike center makes the Lorentz quotient a valid full-dimensional H9 cone.

All 70 cap formulas and their complete pairings are replayed. Each cap has norm 2, pairings 0 or minus phi with the leaf walls, and pairing minus one with each of the two center walls. The cap-cap pairing is `2-phi^2*h/4`, with Hamming distance `h` in `{2,4,6,8}`. Among finite dihedral labels, `2*cos(pi/m)` belongs to `Q(phi)` only for `m=2,3,5`, giving values `0,1,phi`; the three noncomplementary distances fail these possibilities or have positive pairing. Only complementary caps are compatible, and their pairing is dotted.

The exact escape for the selected complementary pair survives a small positive perturbation toward the timelike center, making all inequalities strict while retaining positive Lorentz square. The verifier separately checks its permutation transport for all 35 complementary pairs. Every compatible subset of the 70-cap dictionary extends to one of those pairs, so every such subset has infinite volume. This excludes only this finite dictionary; it makes no claim about additional mirrors or all caps breaking the symmetry.

Replay:

```sh
python scripts/fresh_h9_symmetric_balanced_caps_verify_20261013.py
```

The [PASS record](fresh_h9_symmetric_balanced_caps_verify_20261013.json) reports all 35 strict escapes and rejects the forged cap.
