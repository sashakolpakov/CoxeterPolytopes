# H9/14 Record 0, `a=1`, source 1: compact CEGAR progress

Status date: 2026-08-08.

## Current bounded checkpoint

The alternate local-label-core CEGAR track for the fixed Record-0,
`a=1`, source orbit 1 has completed exactly 56 core iterations.  Its frozen
checkpoint is

`data/enumerations/h9_14_record0_a1_source1_label_core_cegar_checkpoint_56_more.json`

with SHA-256

```text
a611447e70b7964c3643c5d837917427dd92c63b66d52ba50650457b40b93a46
```

The exact status remains

```text
terminal_solver_status             OPTIMAL
stopped_reason                     core_cutoff
excluded local-core orbits         56
distinct orbit-image clauses       401100
fixed-source compact no-go         false
finite-volume noncompact checked   false
theorem certificate                false
```

The independent checkpoint-56 audit rebuilt all 56 complete source orbits
and verified that their 401,100 clause signatures are pairwise distinct.
This upgrades the previous orbit-size sum to an exact distinct-clause count.
`OPTIMAL/core_cutoff` only means that CP-SAT found the next incidence target
before the requested 56-core
bound was reached.  It is not an infeasibility result.

## Cores 53--56

Starting from the independently audited 52-core checkpoint, exact resumes
added iterations 52--55.  The previous 52 iteration records and fixed model
metadata are structurally unchanged.  The new exact metrics are:

| human core | representative digest | local supports | local vertices | independent passing states | orbit | stabilizer | branches / conflicts | wall time |
|---:|---|---|---:|---|---:|---:|---:|---:|
| 53 | `2b6671d5510370e353e50aa19e229e75097d968a0796aeddb6dbcb0f2a6790cf` | fixed `A4` + 3 `A5` | 6 | `108,98,16,0` | 504 | 480 | `793961 / 96087` | 146.091561 s |
| 54 | `10b118cba05ea78a4cfc67b1587950cc82863bea2e983dd8d09478379cb4e722` | fixed `A4` + 2 `A5` | 6 | `108,110,0` | 504 | 480 | `1930858 / 315027` | 625.765564 s |
| 55 | `f50a6e1bd82af47b227fc4eb70186e0bae30ddc779084adaa41fb9ef0626a6da` | fixed `A4` + 3 `A5` | 2 | `108,256,110,0` | 2,016 | 120 | `3842878 / 724416` | 1,401.682469 s |
| 56 | `09a56ee03c235d984358f4f368d013e805fe15acb4cb4a5892d72fc33ca54257` | fixed `A4` + 4 `A5` | 2 | `108,350,122,21,0` | 10,080 | 24 | `1252517 / 170032` | 253.580417 s |

All four complete targets again have 178 selected facet complements and 178
rank-nine actual vertices.  Their pairwise-disjoint full source orbits add
exactly 13,104 distinct clauses, raising the global count from 387,996 to
401,100.  The complete independent local/orbit replay and global
cross-core-disjointness census pass; see
[`h9_14_record0_a1_source1_checkpoint56_increment_audit.md`](h9_14_record0_a1_source1_checkpoint56_increment_audit.md).
This remains a bounded progress certificate for this exact CEGAR trajectory,
not a lower bound on the size of every possible certificate and not a
fixed-source compact no-go theorem.

## Cores 51--52

Starting from the independently audited 50-core checkpoint, the exact resume
added iterations 50 and 51.  The previous 50 iteration records, base-model
counts, and source-group order are structurally unchanged.  The new exact
metrics are:

| human core | representative digest | local supports | local vertices | independent passing states | orbit | stabilizer | branches / conflicts | wall time |
|---:|---|---|---:|---|---:|---:|---:|---:|
| 51 | `e0b91dc79bf9f4bfd7fc8db11d9a618a92d07480070d0d08e78a156d0f2a89a4` | fixed `A4` + 3 `A5` | 2 | `108,110,48,0` | 10,080 | 24 | `1251412 / 171427` | 339.801197 s |
| 52 | `e47a250ee632eec60f7f8d8690f2736888d7eb1abc6548d9f60125370d087bbc` | fixed `A4` + 5 `A5` | 2 | `108,110,62,103,10,0` | 10,080 | 24 | `1603138 / 246893` | 493.029114 s |

Both complete targets again have 178 selected facet complements and 178
rank-nine actual vertices.  Their full source orbits raise the cumulative
orbit-image count by 20,160, from 367,836 to 387,996.  The independent audit
rebuilds the two orbits with hashes

```text
c1fe586ce6ce9697c976a04f78bae6a6a6d774562633f48f9ae6386ceae48735
b89b8356f2372c2bb59300ad3292088c0ec80746204990636a6d15fec9694361
```

The bounded checkpoint and incremental audit pass exactly; see
[`h9_14_record0_a1_source1_checkpoint52_increment_audit.md`](h9_14_record0_a1_source1_checkpoint52_increment_audit.md).
This is still not a fixed-source compact no-go theorem.

## Cores 49--50

Starting from the independently audited 48-core checkpoint with SHA-256

```text
1d1f2ae760364730db8fc3aaaed78b2a3514237182174697186270d85c332931
```

the bounded resume added iterations 48 and 49.  The previous 48 iteration
records, base-model counts, and source-group order are unchanged exactly.
The new exact metrics are:

| human core | representative digest | local supports | local vertices | independent passing states | orbit | stabilizer | branches / conflicts | wall time |
|---:|---|---|---:|---|---:|---:|---:|---:|
| 49 | `030da2f5300465070a04f77c382eb641ac82e6dbf1bac2276c56fdbadd29bcae` | fixed `A4` + 3 `A5` | 4 | `108,110,16,0` | 1,008 | 240 | `179889 / 21602` | 40.141902 s |
| 50 | `14161bd66d1218225e6c55aab23b50c03cd9fb61dd407857ec7706509def0659` | fixed `A4` + 4 `A5` | 2 | `108,350,61,17,0` | 2,016 | 120 | `1524589 / 231472` | 285.209996 s |

Both complete targets again have 178 selected facet complements and 178
rank-nine actual vertices.  The two new orbits raise the cumulative
orbit-image count from 364,812 to 367,836.

## Independent exclusion check

Both new cores were replayed with the standalone Lanner tables and finite
Coxeter-diagram classifier from the no-trust source-3 auditor, rather than
with the exploratory core extractor.  It independently checked that

1. each support family consists of the fixed `A4={0,2,3,4}` and exactly the
   stored independent `A5` supports;
2. every stored facet complement is the complement of its rank-nine actual
   vertex;
3. the installed face restrictions have exactly the same projected
   finite-diagram semantics as those actual vertices;
4. both complete local label joins terminate with zero states; and
5. generator-orbit BFS for the fixed-source `S3 x S8` action gives exactly
   1,008 and 2,016 images, all inside the missing-five ground set.

All local and orbit checks passed.  The independently rebuilt orbit has
SHA-256 values

```text
b8c5321b73ab5716d056369fad95d4b5d62e832c044fe8b637d8feb3b37fd001
5538a6b0026f050a418636cc4543b64c75122a11d11390f4c2e70c8fab9db126
```

Thus cores 49 and 50 are valid lossless compact-label exclusions.  The
preceding core 48 remains independently frozen in
[`h9_14_record0_a1_source1_checkpoint48_core_audit.md`](h9_14_record0_a1_source1_checkpoint48_core_audit.md),
and the new incremental certificate is in
[`h9_14_record0_a1_source1_checkpoint50_increment_audit.md`](h9_14_record0_a1_source1_checkpoint50_increment_audit.md).
This does not turn the bounded search into a fixed-source no-go theorem.

## Reproduction and exact scope

The latest bounded resume was

```bash
caffeinate -i env PYTHONPATH=src:. python3 \
  scripts/explore_h9_14_record0_a1_label_core_cegar.py 1 \
  --max-cores 56 \
  --resume data/enumerations/h9_14_record0_a1_source1_label_core_cegar_checkpoint_53_more.json \
  --write data/enumerations/h9_14_record0_a1_source1_label_core_cegar_checkpoint_56_more.json \
  --summary-only
```

Only a terminal, no-timeout `INFEASIBLE` result after installing all complete
core orbits could close this fixed source by this route.  The present result
concerns compact Coxeter labelings for Record 0, `a=1`, source orbit 1 only.
It says nothing about finite-volume noncompact polytopes, other source
orbits, the whole Record-0 shadow, or the global H9/14 problem.
