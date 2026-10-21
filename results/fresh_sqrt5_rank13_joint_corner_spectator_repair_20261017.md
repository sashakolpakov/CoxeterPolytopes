# Batched spectator-family repair of the partial rank-13 weighting model

Adding all 682 certified local quadrilateral contexts preserves feasibility below 12. The exact maximum vertex payment is **64616136604/9999999999**, approximately 6.461613661. Every added context receives at least 10000000000/9999999999, strictly more than one.

The [source](../scripts/fresh_sqrt5_rank13_joint_corner_spectator_repair_20261017.py) and [saved model and certificates](fresh_sqrt5_rank13_joint_corner_spectator_repair_20261017.json) retain all predecessor rows and add the [complete finite spectator families](fresh_sqrt5_rank13_five_core_spectators_20261017.md) in one batch. Each of the 1,660 newly encountered whole vertex forests receives all 78 root pairs, grouped into exact unordered-pair orbits, with independent nonnegative payment variables and a full vertex budget.

The resulting model has 2,290 forests, 77,587 variables and 5,006 distinct face inequalities. Its 6,125 literal input representations include the two first-repair quadrilaterals both as individual certificates and as members of the complete spectator families; identical inequalities are deduplicated and both provenances are retained. No prior constraint or variable is removed or fixed.

The numerical minimum remains approximately 6.46161366006. A separate exact dual lower bound is 2326180916772/360000000037. Its 59 face rows involve 18 whole vertex types: 32 new non-ADE triangle rows, 18 earlier retained quadrilaterals and nine old triangle rows. None of the new spectator rows carries dual mass. The largest whole-vertex dual loads occur at A4+D5+H4, D5+2H4, A1+D4+2H4, D4+D5+H4 and A1+A4+2H4. This points toward missing stars at the double-H4 destination types, rather than further variants of the same spectator families.

Generation and a fresh saved reconstruction replay passed. Replay checks the exact primal against every face and every complete vertex budget, and checks the exact dual without calling the optimizer:

```sh
python scripts/fresh_sqrt5_rank13_joint_corner_spectator_repair_20261017.py
python scripts/fresh_sqrt5_rank13_joint_corner_spectator_repair_20261017.py --verify
```

This remains a partial local model. Of the 6,789 arithmetic-compatible finite forests of rank 13, 4,499 remain absent, including 1,614 whose old baseline exceeds 12. The new rule uses only 1,493 nonzero orbit variables and has substantial unconstrained face coverage. Complete triangle and quadrilateral domains at the changed types are still needed. The surviving slack is therefore not evidence of a completed dimension bound or compact ambient realizability.
