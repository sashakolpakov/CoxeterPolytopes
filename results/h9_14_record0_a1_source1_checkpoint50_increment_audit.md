# Independent audit of source-1 compact CEGAR cores 49--50

**Verdict:** `PASS`, with bounded nonterminal scope.

The Record-0, `a=1`, source-1 compact-label CEGAR has advanced from 48 to 50
complete local-core orbits.  The two appended cores were replayed with the
independent classical Lanner and finite-diagram code and full `S3 x S8`
generator-orbit BFS.

| core | local supports | local vertices | passing states | orbit | stabilizer |
|---:|---|---:|---|---:|---:|
| 49 | fixed `A4` + 3 `A5` | 4 | `108,110,16,0` | 1,008 | 240 |
| 50 | fixed `A4` + 4 `A5` | 2 | `108,350,61,17,0` | 2,016 | 120 |

Their independently rebuilt orbit hashes are respectively

```text
b8c5321b73ab5716d056369fad95d4b5d62e832c044fe8b637d8feb3b37fd001
5538a6b0026f050a418636cc4543b64c75122a11d11390f4c2e70c8fab9db126.
```

Every stored actual vertex has rank nine, every stored facet complement is
exact, and each installed face cut has the same projected semantics as its
actual vertices.  The cumulative sum of full core-orbit sizes is now 367,836;
this is not asserted to count distinct clauses across different cores.

The checkpoint status remains `OPTIMAL/core_cutoff`.  It is therefore not a
fixed-source compact no-go theorem and makes no statement about other source
orbits or finite-volume noncompact polytopes.

Reproduction:

```bash
PYTHONPATH=src:. python3 \
  scripts/audit_h9_14_record0_a1_source1_checkpoint50_increment.py \
  --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_h9_14_record0_a1_source1_checkpoint50_increment_audit.py
```

Machine-readable archive:
[`h9_14_record0_a1_source1_checkpoint50_increment_audit.json`](../data/enumerations/h9_14_record0_a1_source1_checkpoint50_increment_audit.json).
