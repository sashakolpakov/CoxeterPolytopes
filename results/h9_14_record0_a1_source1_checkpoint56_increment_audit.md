# Independent audit of source-1 compact CEGAR cores 53--56

**Verdict:** `PASS`, with bounded nonterminal scope.

Frozen checkpoint SHA-256:

```text
a611447e70b7964c3643c5d837917427dd92c63b66d52ba50650457b40b93a46
```

The Record-0, `a=1`, source-1 compact-label CEGAR has advanced from 52 to
56 complete local-core orbits.  The four appended cores were replayed with
the independent classical Lanner and finite-diagram code, and their complete
`S3 x S8` source orbits were rebuilt by generator BFS.

| core | local supports | local vertices | independent passing states | orbit | stabilizer | solve time |
|---:|---|---:|---|---:|---:|---:|
| 53 | fixed `A4` + 3 `A5` | 6 | `108,98,16,0` | 504 | 480 | 146.091561 s |
| 54 | fixed `A4` + 2 `A5` | 6 | `108,110,0` | 504 | 480 | 625.765564 s |
| 55 | fixed `A4` + 3 `A5` | 2 | `108,256,110,0` | 2,016 | 120 | 1,401.682469 s |
| 56 | fixed `A4` + 4 `A5` | 2 | `108,350,122,21,0` | 10,080 | 24 | 253.580417 s |

The independently rebuilt new-orbit hashes are respectively

```text
5b4133e716b26897816a2297c3978996db85397e520eeb0bde435d5492c99f61
24b4b54878487513bdeeb43554989909f932c69f04cdd5c92b4ebb84f5d67de5
7adef54b708e265d991d69896b0f9a51be9d144a59a136dd63168695cb560149
e7cc7a1de3dff93fdedb77dd54785bfd7ad35f4a79d6bb46959f40b78d8449eb
```

Every stored local actual vertex has rank nine, every local facet complement
is exact, and each installed face cut has the same projected finite-diagram
semantics as its actual vertices.  The first 52 iteration records are
structurally unchanged from the frozen checkpoint-52 input; the separately
frozen checkpoint 53 is also an exact prefix of checkpoint 56.

## Exact global orbit census

Unlike the earlier incremental audits, this verifier rebuilds **all 56**
complete source orbits.  Their canonical representatives are pairwise
distinct, and their union has exactly **401,100 distinct clause signatures**.
Thus the increase from checkpoint 52 is exactly 13,104 distinct clauses, not
merely an orbit-size sum with possible cross-core multiplicity.  The global
hashes are

```text
canonical representatives  c0c0bb19fb3860ae4779c0b4459f1aa4031f471ca20535102d716296a25acd6a
distinct clause signatures  09556bd47f128544e250eaa10761e94b048ac84715e825cbacbb5c8887407fbb
```

This sharper clause census does **not** make the CEGAR search terminal.  The
checkpoint status is still `OPTIMAL/core_cutoff`: CP-SAT found another
incidence target before the requested 56-core bound.  Consequently this is
not a fixed-source compact no-go theorem.  The calculation concerns compact
labelings only and neither checks nor claims anything about finite-volume
noncompact polytopes.

Reproduction:

```bash
PYTHONPATH=src:. python3 \
  scripts/audit_h9_14_record0_a1_source1_checkpoint56_increment.py \
  --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_h9_14_record0_a1_source1_checkpoint56_increment_audit.py
```

Machine-readable archive:
[`h9_14_record0_a1_source1_checkpoint56_increment_audit.json`](../data/enumerations/h9_14_record0_a1_source1_checkpoint56_increment_audit.json).
Its SHA-256 is
`5f8c24c92dfa067af48fec06d7747d4059c8b2c30b409cbf9ad21e3c346d601c`.
