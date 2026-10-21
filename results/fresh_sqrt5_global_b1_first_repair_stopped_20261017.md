# Preserved unresolved first-stop repair

The first repair of the global dimension-21 payment proposal was stopped before it completed its first optimization iteration. This package preserves the exact initial cuts and the authoritative log. It proves neither feasibility nor infeasibility and supplies no dimension bound.

The rebuilt model has 32,554 variables and 70,156 distinct face-payment profiles. It includes the 996 triangle graphs, 37 earlier paid quadrilateral stops, 1,979 inherited baseline-zero quadrilateral rows, all 196 admissible finite singular first extensions, and all 14 admissible unpaid dotted first completions, with their common spectator contexts. The original proposal fails 52 of the 210 added stopping cases.

All 4,762 inherited literal vertex cuts were retained. The preserved state has no completed solve, no new oracle iteration, and no numerical weight vector. The task-owned process was explicitly stopped after approximately half an hour; its authoritative session returned exit status 143. The replacement using a bounded interior-point solve is a separate experiment and is outside this record.

The [checker](../scripts/fresh_sqrt5_global_b1_first_repair_stopped_20261017.py) validates the bound source hashes, the exact equality of the preserved and initial vertex cuts, and the completed constraint-build counters in the [preserved log](../data/enumerations/fresh_sqrt5_global_b1_first_repair_stopped_20261017.log). The [saved result](fresh_sqrt5_global_b1_first_repair_stopped_20261017.json) records these checks; the [preserved checkpoint](../data/enumerations/fresh_sqrt5_global_b1_first_repair_stopped_20261017.json) contains the literal cuts. It binds the earlier wrapper's completed literal recount and independent component assembly without repeating that arithmetic.

Generation and a fresh saved replay both pass:

```sh
python scripts/fresh_sqrt5_global_b1_first_repair_stopped_20261017.py
python scripts/fresh_sqrt5_global_b1_first_repair_stopped_20261017.py --verify
```

The deeper finite and dotted baseline-one attachment closures remain separate. This unresolved repair does not affect the frozen compact arithmetic dimension-21 upper-bound candidate.
