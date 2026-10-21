# Bounded repair: unresolved after five minutes

The replacement solver reached its configured 300-second limit without deciding whether the selected global dimension-21 model admits a budget of 40. This is an unresolved computation, not an infeasibility certificate.

The model and its 32,554 variables, 70,156 face profiles and 4,762 literal vertex cuts are unchanged from the [preserved original repair](fresh_sqrt5_global_b1_first_repair_stopped_20261017.md). The interior-point phase stopped making progress; the solver's automatic simplex cleanup then reached the time limit. The exact-recovery driver rejected this non-success status and exited before any numerical vector, rational certificate or new oracle cut was saved. No minimum-budget fallback ran.

The [preservation checker](../scripts/fresh_sqrt5_global_b1_first_repair_ipm_timeout_20261017.py) binds the [bounded wrapper](../scripts/fresh_sqrt5_global_b1_first_repair_ipm_20261017.py), [authoritative log](../data/enumerations/fresh_sqrt5_global_b1_first_repair_ipm_timeout_20261017.log), [checkpoint](../data/enumerations/fresh_sqrt5_global_b1_first_repair_ipm_timeout_20261017.json) and initial source manifests. It verifies that all initial cuts remain identical and that no completed optimization iteration or vector was recorded. Its [saved result](fresh_sqrt5_global_b1_first_repair_ipm_timeout_20261017.json) passes generation and a fresh saved replay.

```sh
python scripts/fresh_sqrt5_global_b1_first_repair_ipm_timeout_20261017.py
python scripts/fresh_sqrt5_global_b1_first_repair_ipm_timeout_20261017.py --verify
```

Complete finite and dotted attachment closures remain separate. Neither this timeout nor the earlier proposal establishes the attempted compact arithmetic bound 20. The older frozen bound-21 candidate is unchanged.
