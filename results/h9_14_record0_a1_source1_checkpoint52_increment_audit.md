# Independent audit of source-1 compact CEGAR cores 51--52

**Verdict:** `PASS`, with bounded nonterminal scope.

The Record-0, `a=1`, source-1 compact-label CEGAR has advanced from 50 to 52
complete local-core orbits.  The two appended cores were replayed with the
independent classical Lanner and finite-diagram code, and their complete
`S3 x S8` source orbits were rebuilt by generator BFS.

| core | local supports | local vertices | independent passing states | orbit | stabilizer |
|---:|---|---:|---|---:|---:|
| 51 | fixed `A4` + 3 `A5` | 2 | `108,110,48,0` | 10,080 | 24 |
| 52 | fixed `A4` + 5 `A5` | 2 | `108,110,62,103,10,0` | 10,080 | 24 |

The independently rebuilt orbit hashes are respectively

```text
c1fe586ce6ce9697c976a04f78bae6a6a6d774562633f48f9ae6386ceae48735
b89b8356f2372c2bb59300ad3292088c0ec80746204990636a6d15fec9694361
```

Every stored actual vertex has rank nine, every local facet complement is
exact, and each installed face cut has the same projected label semantics as
its actual vertices.  The first 50 iteration records are structurally
unchanged from the frozen checkpoint-50 input.

The cumulative sum of complete source-orbit sizes rises by 20,160, from
367,836 to 387,996.  This is a sum with multiplicity: it is **not** a claim
that clauses belonging to different core orbits are pairwise distinct.

The checkpoint status remains `OPTIMAL/core_cutoff`.  It is therefore not a
fixed-source compact no-go theorem.  The calculation concerns compact
labelings only; it neither checks nor claims anything about finite-volume
noncompact polytopes.

Reproduction:

```bash
PYTHONPATH=src:. python3 \
  scripts/audit_h9_14_record0_a1_source1_checkpoint52_increment.py \
  --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_h9_14_record0_a1_source1_checkpoint52_increment_audit.py
```

Machine-readable archive:
[`h9_14_record0_a1_source1_checkpoint52_increment_audit.json`](../data/enumerations/h9_14_record0_a1_source1_checkpoint52_increment_audit.json).
