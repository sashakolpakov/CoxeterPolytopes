# Independent audit of source-1 compact CEGAR core 48

**Verdict:** `PASS`, with deliberately bounded scope.

The Record-0, `a=1`, source-1 compact CEGAR checkpoint now contains 48
lossless local compact-label cuts.  The new iteration 47 is independently
replayed using the classical Lanner tables and finite-diagram classifier from
the frozen source-3 hostile auditor, not the exploratory source-1 core
extractor.

Its local contradiction consists of the fixed `A4`, three additional `A5`
supports, and four rank-nine vertex restrictions.  The independent successive
label-state counts are

```text
108, 110, 16, 0.
```

The four stored facet complements are exactly the complements of the four
stored actual vertices, and the installed face cuts have the same projected
finite-diagram semantics.  Nine adjacent-transposition generators rebuild
the complete `S3 x S8` fixed-source group.  Generator-orbit BFS gives exactly
168 images, stabilizer order 1,440, and orbit hash

```text
81cd3f7b7c749f7ded551b2565ec6987122054d74411288346bc07afcd85e3a9.
```

All checks pass.  The 48 cuts have cumulative orbit-size sum 364,812; this is
not asserted to be a count of pairwise distinct clauses across different
cores.

The terminal solver status is still `OPTIMAL/core_cutoff`.  Therefore this is
**not** a fixed-source compact no-go theorem, does not address other sources,
and says nothing about finite-volume noncompact polytopes.  Only a no-timeout
terminal `INFEASIBLE` result after all installed full-source cuts would close
this source by this route.

Reproduction:

```bash
PYTHONPATH=src:. python3 \
  scripts/audit_h9_14_record0_a1_source1_checkpoint48_core.py --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_h9_14_record0_a1_source1_checkpoint48_core_audit.py
```

Machine-readable archive:
[`h9_14_record0_a1_source1_checkpoint48_core_audit.json`](../data/enumerations/h9_14_record0_a1_source1_checkpoint48_core_audit.json).
