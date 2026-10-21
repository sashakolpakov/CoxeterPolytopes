# Independent review of the opposite-gap cap

**PASS: this second cap repairs the targeted opposite direction and preserves all 45 preceding walls, but the 46-wall H9 chamber still has infinite volume.** The [independent checker](../scripts/fresh_h9_quartic_completion_20261013_opposite_independent.py) reconstructs the cap and verifies the [saved construction](../data/enumerations/fresh_h9_quartic_completion_20261013_opposite.json), without importing its producer. See the [PASS record](fresh_h9_quartic_completion_20261013_opposite_independent.json).

For the authoritative seed inverse `I`, exact arithmetic gives

```
I66=(4+2phi)/5, I99=(-18+11phi)/10, I69=0,
D9=(2-2I66)/I99=12+20phi.
```

The cap has negative seed pairings `-sqrt(2)` at zero-based index 6 and `-sqrt(D9)` at index 9, and zero elsewhere. Positive square roots are fixed by the exact ordered-field arithmetic. Its norm is exactly 2, its seed-6 angle is pi/4, its seed-9 relation is dotted, and the other eight seed relations are right angles. All 35 nonseed relations are independently checked to be dotted. The original strict timelike center remains strict.

The old feasible spacelike ray with `y6=2-phi`, `y9=1`, and other coordinates zero has norm `1-phi/2>0` and center pairing `-7/2`. All 45 preceding inequalities hold, with only old caps 67 and 310 vanishing among the nonseed walls. The new cap cuts this ray strictly. It also cuts exactly saved targets 1 and 30. Its `{6,9}` zero pattern differs from all 3,031 literal preceding entries, so its mirror is new within the archived search. No global novelty claim is made.

For the resulting 46-wall chamber, a different archived vector has positive exact Lorentz square, negative center pairing, and every one of the 46 inequalities strict. This proves infinite volume. The conclusion uses this exact certificate, not the numerical vertex counts that proposed it.

The checker adjoins the two radicals in the opposite order to the producer and explicitly converts the saved basis. Earlier caps are compared in exact composita with their positive radicals. Corrupted cap and escape certificates are rejected. During source review a concern about a reused zero-pattern filter was resolved by checking the current parameterized helper; the independent `{6,9}` census also passes with zero matching prior patterns. No unresolved issue remains.

```sh
python -m scripts.fresh_h9_quartic_completion_20261013_opposite_independent
```

This escape obstructs finite volume for the literal 46-wall chamber and its wall subsets. A subsequent [four-wall cusp obstruction](fresh_h9_fixed_core_cusp_obstruction_20261013.md) also excludes finite-volume completion by acute-compatible wall addition while retaining this core. Replacing inherited walls remains outside that theorem.
