# Record 1, `a=1`, source-orbit 0: terminal quotient engine

**Status: exact engine and bounded quotient theorem, 2026-08-07.  No full
seed class has been run.**  This is the corrected path to a complete incidence
census after the old `--max-models 10` cutoff was mistakenly described as an
exact nine-target enumeration.  It does not yet close this source orbit.

## Exact finite reduction

The fixed source consists of graph rows

```text
{0,13}, {1,13}, {2,13}
```

and the positive four-nonface `{0,3,4,5}`.  Its ground-set stabilizer is

```text
S2({1,2}) x S3({3,4,5}) x S7({6,...,12}),
```

of order 60,480.  Independent reconstruction of every vertex-, edge-, and
triple-link Macaulay interval gives:

| object | exact count |
|---|---:|
| local link constraints | 433 |
| independent five-sets | 1,497 |
| forced five-shadow | 9 |
| admissible additional-five candidates | 1,278 |
| base transversals of order at most four | 259 |
| base order-five covers | 634 |
| additional rows required | 19 |
| scheduled facets | 158 |
| scheduled active ridges | 711 |

The 1,278 candidate rows form exactly **21** orbits under the fixed-source
group.  Their orbit-size histogram is

```text
1^1 3^1 14^1 21^5 35^2 42^1 63^2 70^2 105^3 126^2 210^1.
```

This 21-class split is exact and is the durable restart boundary for the full
census.

## Why the class quotient is lossless

Each class fixes one representative candidate as a distinguished seed row.
Partial families in that class are quotiented only by the stabilizer of the
seed.  This restriction is mathematically necessary: pruning a search
conditioned on containing a seed by the full source group could map the seed
away and lose extensions that are not reachable in the kept branch.

Every final target contains at least one additional row.  Move that row to
the representative of its full-source candidate orbit; the resulting target
lies in the corresponding seed class.  Hence the union of all 21 complete
seed-class searches covers every full-source target orbit.  Final targets are
canonicalized under the full source group and deduplicated across classes.

Nauty canonicalizes the fixed colored bipartite incidence graph.  The memo
key stores the canonical row-to-ground incidence bits in six inline
`uint64_t` words.  Equality compares all 384 bits, so hash collisions affect
performance only.  At full depth the logical key uses 42 bytes.  The bounded
verification build reconstructs and compares every bit of nauty's dense
canonical graph.

## Exact search model

Within a seed class the search enforces:

1. all 433 local incidence intervals;
2. destruction of every base transversal of order at most four;
3. exactly 158 surviving order-five covers;
4. literal blocker dualization at every terminal candidate; and
5. a pure 158-member order-five blocker with exactly 711 active six-set
   ridges, each of degree two.

The least-current-destroyer branching is complete.  Partial-state quotienting
uses the distinguished-seed stabilizer, while target quotienting uses the full
source group.

## Independent bounded replay

The independent Python audit imports neither the C++ canonicalizer nor its
orbit keys.  It explicitly constructs all 60,480 source permutations and
compares equivalence relations rather than opaque key bytes.  It proves:

- all 1,278 singleton candidates give the same 21 orbits in the direct and
  nauty computations;
- for seed mask `31`, all 1,277 depth-two families give exactly 32 orbits in
  both computations under the order-20,160 seed stabilizer;
- a deterministic 58-family multirow pool, including the complete known
  nineteen-row target, gives 49 equivalence classes in both computations;
- the packed key round-trips every dense adjacency bit; and
- a real 64-node cutoff emits four non-proof checkpoints, a non-proof class
  terminal, and a non-proof run terminal.  It cannot be parsed as complete.

The bounded transcript digest is

```text
b5fc6064430986f9b0d152d2e052543f1ca66d6b9252c5837f1aa2ade9e5a530.
```

## Interruption and terminal policy

The engine writes flushed JSON Lines events:

```text
run_header
class_inventory
class_started
checkpoint ...
target ...
class_terminal
run_terminal
```

A checkpoint records progress but is not a resumable proof state.  A completed
`class_terminal` is the durable restart boundary.  An interrupted class must
be restarted from zero; already completed classes need not be repeated.

Only

```json
{"event":"run_terminal","terminal_status":"complete_exhaustive",
 "enumeration_complete":true,"whole_census_complete":true}
```

from an all-class run, or the corresponding complete output of the independent
log assembler, certifies the whole source quotient.  Targets, checkpoints,
cutoffs, inventory terminals, and selected-class terminals explicitly carry
false whole-census status.

The assembler accepts interrupted and cutoff logs but ignores them as proof
boundaries.  It independently rebuilds the class inventory and literally
rechecks every streamed target's local intervals, 158-facet blocker, and
711-ridge law.  It emits a complete terminal only after receiving all 21
complete class boundaries.

## Completion path

1. Compile the engine once without `RECORD1_VERIFY_PACKED_KEY` for production.
2. Run each of classes `0,...,20` separately, preserving its JSONL log.  This
   can proceed around interruptions at completed class boundaries.
3. Assemble all logs with `--require-complete`.  Any missing or cutoff class
   makes the assembler exit incomplete.
4. Deduplicate the streamed full-source target keys.  The already frozen
   60,480-member target orbit is compact-label impossible.  Any additional
   incidence orbit must receive its own exact label audit.
5. Only after every assembled target orbit is label-excluded may this source
   orbit be promoted to a compact no-go.

No full class was launched while the live Record-31 computation owned memory.

## Artifacts and replay

- `scripts/explore_h9_14_record1_a1_orbit0_quotient.cpp`
- `scripts/independent_h9_14_record1_a1_orbit0_quotient_audit.py`
- `scripts/assemble_h9_14_record1_a1_orbit0_quotient_logs.py`
- `data/enumerations/h9_14_record1_a1_orbit0_quotient_bounded_audit.json`
- `tests/test_h9_14_record1_a1_orbit0_quotient_engine.py`

```text
PYTHONPATH=src:. python3 \
  scripts/independent_h9_14_record1_a1_orbit0_quotient_audit.py --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_h9_14_record1_a1_orbit0_quotient_engine.py
```

Frozen SHA-256 digests:

```text
e6dfb0754fcf39cce13320e064bbba6c8be5576323ce0a63cd050e0bdef49316  engine
d43361d85ea317e6d1dd305a7db17f32164d83ac8577cc233152c6f9323a1bc1  bounded audit
753d944de8d6aa814b5091a8965f1110303a6057b4419fd3c4ccf05c2a447e8a  log assembler
76ab2dad7278deba46b9fe2917640ed6fbaf88e2e799a5012b3164ec0975a17a  audit JSON
1c249efee8dab9e8c943451d022c383d3d7ed03925429ccd64dc7a42deadd539  tests
```

## Scope

This is an exact bounded quotient theorem plus a production-ready enumerator.
It is not a completed target census, not a compact-label theorem beyond the
already frozen 60,480-member orbit, not a polytopality theorem, and not a
finite-volume statement.
