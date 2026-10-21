# Quartic order-exchange obstruction: proof and reproduction

20 October2026. The literal lattices
`N_n=<-(phi+sqrt2-1)>+I_n` over `Z[sqrt2,phi]`, n>=8, have
noncofinite full integral reflection groups. This is a construction-family
exclusion, not a dimension ceiling, new example or exclusion of every
lattice in its quadratic-space commensurability class.

The [complete proof](../results/compact_biquadratic_unit_order_exchange_obstruction_20261020.md)
imports the [first edge](../results/compact_h9_biquadratic_unit_diagonal_first_edge_20261020.md)
and [second edge](../results/compact_h9_biquadratic_unit_second_edge_20261020.md).
Two different full rank-two faces have the same literal N_(n-2)
arrangement. Performing these operations in opposite orders gives two
actual full vertex germs. An explicit integral map between them
preserves the chamber and has infinite order. Its powers represent
infinitely many distinct reflection-group cosets.

Full [root](../results/review_compact_biquadratic_order_exchange_root_20261020.md)
and [outside-author](../results/review_biquadratic_unit_order_exchange_batch_guard_20261020.md)
mathematical reviews pass. The finite Coxeter classification, elementary
number-field foundations and Borel--Harish-Chandra arithmetic lattice
theorem are named external inputs. Actual edges, arbitrary-crossing-mirror
seals and chamber identification are written proofs, not script premises
quietly relabelled as external theorems.

## Exact arithmetic replay

From the repository root, run this single bounded job:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- scripts/check_quartic_order_exchange_period_20261020.py --verify
```

The wrapper enforces one proof process, one numerical thread, a59-second
wall limit and1GiB process-tree memory cap. Start gates are40GiB disk
and6GiB available memory; stop gates are8GiB and3GiB. Its process-list
watchdog must be permitted to run; do not bypass the watchdog.

The [fixed checker](../scripts/check_quartic_order_exchange_period_20261020.py)
recomputes both five-dimensional frames, both integral frame inverses,
the period and its inverse, the metric, oriented root images, future
points and the exact trace. It compares its entire result with the
[17,021-byte certificate](../results/quartic_order_exchange_period_20261020.json).
The active trace is `196phi+121+sqrt2*(138phi+86)>5`.
All input hashes, including checker and helper, are checked before and
after arithmetic. No root search, catalogue enumeration or indefinite
matrix elimination runs.

[Generation](../results/quartic_order_exchange_period_generation_20261020.json)
and [fresh replay](../results/quartic_order_exchange_period_replay_20261020.json)
pass:0.215 seconds total and56,311,808 bytes peak process-tree memory.
The [execution review](../results/review_quartic_order_exchange_execution_root_20261020.md)
records their exact scope. The [immutable manifest](../results/quartic_order_exchange_20261020_freeze.json)
pins the proof, reviews, checker, certificate, controls and conditional
logic sources. Replay assumes the imported complete geometry; it does
not prove geometric completeness.

## Lean and manuscript scope

The existing [ChamberPeriod kernel](../formal/CoxeterProofs/ChamberPeriod.lean)
has a recorded compiled seven-theorem axiom audit. It proves the abstract
trace/period/chamber-coset implications from explicit premises. This new
quartic geometric instance has not been instantiated in Lean, and the
matrix check is not a Lean proof. No new Lean or TeX job ran. No standalone
manuscript is added for this literal-family search exclusion.

Stop further edge, prefix, Pell, period and compactness searches in
these literal N_n for n>=8. N2 has a compact reflective quadrilateral;
no lower-dimensional cutoff, worldwide priority or field-wide exclusion
is claimed. Earlier OPEN first-edge notes are preserved historical stages.
