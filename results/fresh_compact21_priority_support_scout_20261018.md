# Bounded priority-support scout: exact candidate slightly above forty

18 October 2026. The new active/inactive charge rows improve the old
saved obstruction, but the first sparse support solve has not reached
the target forty. HiGHS reported numerical objective 40.00447916245.
Its rounded candidate was checked using exact fractions and has

  z=62506999/1562500=40.00447936.

All 2,805 selected inequalities pass exactly. This is an upper bound
for the selected support problem; it is not a lower-bound certificate,
full-model feasibility or a dimension theorem. A different exact dual
must be recovered before asserting that the support problem cannot
meet forty. The complete geometry theorem remains a real charge
improvement regardless of that later optimization outcome.

The rows are the 2,804 positive-support terms of the frozen A4-source
dual, with its D5+2E8 row replaced by active and inactive rows. There
are 6,121 coordinates, 17,387 nonzeros and 2,671 positive recovered
bonus values. The [127-KiB certificate](fresh_compact21_priority_support_scout_20261018.json)
stores only those positive values, the exact upper coordinate and
compact provenance/resource metadata. All other bonuses are zero.
Upward rounding to a 10^(-8) grid sufficed; the exact uniform face
scaling was one. The upper coordinate was then recomputed directly
from every selected forest row, and every selected inequality was
checked again.

## Resource guard and reproduction

The [scout](../scripts/fresh_compact21_priority_support_scout_20261018.py)
starts a separate numerical child only after a working memory-monitor
preflight. The driver holds the shared proof lock, applies the usual
40-GiB disk/6-GiB available-RAM start gate, limits numerical libraries
and HiGHS to one thread, monitors total process-tree RSS below one GiB,
and kills the process group on the 59-second wall cap, owner exit,
resource failure or excessive output. The child independently watches
its owner. Captured output is capped below 512 KiB. Only a successfully
verified final certificate is written atomically; its temporary is
removed immediately.

The first sandboxed launch stopped before solving because process
inspection was denied. A read-only escalated check confirmed no scout
process remained. The code then added a monitor preflight and cleanup
without the sandbox-denied process-group existence probe. The approved
bounded run used resource-monitor escalation so its live RSS watchdog
could operate. No failed output or temporary file was retained.

The successful run was session 80438, completion chunk `8105e9`, with
1.518 seconds guarded wall time. Combined process-tree peak RSS was
131,203,072 bytes; worker peak was 100,450,304 bytes. The start gate
reported 197.99 GiB free disk and 12.06 GiB available RAM. Captured
worker output was 61,623 bytes; the final certificate is 126,802 bytes.
The global lock was released and the numerical child reaped before
the next lane's calculation.

A fresh exact replay passed in chunk `a02692`. Replay streams the
same selected inequalities and source hashes; it does not import the
numerical libraries or rerun the optimizer:

```sh
python3 -B scripts/fresh_compact21_priority_support_scout_20261018.py --verify
```

There was no Lean process, giant JSON load or background numerical
worker. Full face-case separation and complete finite-forest verification
are not performed by this scout. They are deferred while this selected
candidate itself remains above forty. The immediate next step is to
retain and recover an exact support dual, then identify the forests
which require further geometry. Independent implementation review of
the numerical scout and guard remains separate from the existing
mathematical review of the charge rule.
