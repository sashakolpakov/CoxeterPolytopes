# First exact counterexample-guided repair of the partial rank-13 rule

Adding the two certified unpaid quadrilaterals preserves feasibility below the target 12. The new exact maximum vertex payment is again **21538712201/3333333333**, approximately 6.461613661. Each added face now receives 10000000000/9999999999, strictly more than one.

The [source](../scripts/fresh_sqrt5_rank13_joint_corner_repair1_20261017.py) and [saved model and certificates](fresh_sqrt5_rank13_joint_corner_repair1_20261017.json) retain every original constraint from the [first partial model](fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017.md), and add both [exact rank-13 local quadrilaterals](fresh_sqrt5_rank13_joint_corner_separation_20261017.md). The domain now has 4,326 distinct face constraints, with the same 630 whole vertex forests and 21,492 free nonnegative pair-orbit variables. No old coefficient is fixed and no constraint is removed.

The numerical minimum and exact dual lower bound are unchanged. Thus these two added faces repair the known underpayments without increasing the optimum of this partial constraint family. The rational primal and dual are checked against the full saved model during replay without rerunning the optimizer.

Generation and fresh saved replay passed:

```sh
python scripts/fresh_sqrt5_rank13_joint_corner_repair1_20261017.py
python scripts/fresh_sqrt5_rank13_joint_corner_repair1_20261017.py --verify
```

This is still an incomplete rule. The [already certified local core 7661](fresh_sqrt5_multiplicity_geometry_review_20261016.md), extended orthogonally by A1+A4 or A1+D4 to rank 13, gives further zero-payment faces. A separate successor will test complete finite spectator families of the available certified cores together, including previously absent whole vertex types. No dimension conclusion follows from this first repair.
