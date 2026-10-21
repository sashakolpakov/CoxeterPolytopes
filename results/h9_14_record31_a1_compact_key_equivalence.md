# Record 31 `a=1`: compact canonical-key equivalence

## Outcome

A separate memory-compact implementation of the production `--quotient`
engine now stores exactly the same nauty canonical graph in at most 22 bytes,
rather than storing its complete dense adjacency matrix.  The search tree,
canonicalization, pruning rules, target test, and corrected derived ridge
count 540 are unchanged.

The original source-orbit-3 run was neither stopped nor modified while this
implementation was built and tested.  After its dense representation became
operationally irrecoverable, that exact process was terminated cleanly and a
compact run was restarted from node zero; no state was spliced.

## Why the encoding is lossless

Every graph sent to nauty has a fixed colored bipartite form.  The first
color cell consists of the 14 ground vertices.  The remaining color cells
consist of the four graph-edge rows, the positive-A4 row, and the selected
A5 rows.  There are no ground-to-ground or row-to-row edges.

Nauty respects the supplied ordered color partition.  Consequently its
canonical graph is completely determined by the ordered incidence bits from
the canonical row vertices to the 14 canonical ground vertices.  The compact
key stores:

1. the number of row vertices; and
2. those row-to-ground bits in canonical order.

This reconstructs the original dense canonical graph injectively, so two
compact keys are equal exactly when the former dense keys are equal.

For the largest `a=1` partial family there are 12 row vertices.  The compact
payload therefore has

\[
1+\left\lceil\frac{12\cdot14}{8}\right\rceil=22
\]

bytes.  The old depth-six key for source orbit 3 used a 25-by-25 nauty dense
graph, or 200 bytes, before container and allocation overhead.

## Exact bounded equivalence replay

The dense and compact generators were compiled separately.  The compact
binary enabled an additional per-state assertion that decoded its packed key
and compared every recovered adjacency bit with the dense nauty graph.

Both binaries then completed the production quotient traversal for source
orbit 0, mask 60.  Their parsed JSON transcripts were exactly identical:

| count | value |
|---|---:|
| search nodes | 231,773 |
| depth-node counts | 1, 56, 420, 8,495, 195,250, 27,551 |
| canonical states | 1, 4, 47, 1,393, 51,035, 8,396 |
| ridge-valid target orbits | 1 |

Thus the compact serialization has both a structural injectivity argument
and an exact complete bounded-run replay.  The later compact source-orbit-3
production run also terminated normally at 469,990,577 nodes, but that
terminal result is logically separate from this serialization certificate
and still does not imply a full Record-31 schedule theorem.

## Dense-key memory failure mode and production policy

The original dense-key production run for source orbit 3, mask 8220, exposed
a representation-level memory failure mode.  At its 156-million-node
checkpoint it held 105,141,843 depth-six canonical states.  A depth-six key
is a 25-by-25 dense nauty graph and therefore has a 200-byte payload, so this
single set had more than 21 GB of key payload before string allocations,
hash-table nodes, and buckets.  System measurements showed severe compression
and swap churn.  A hash-table growth event then produced a long interval with
no node checkpoint.

This was **not** a mathematical failure.  The dense run eventually reached
283 million nodes but had covered only four of nine canonical depth-one
classes and 76 of 154 depth-two states.  An exact compact prefix diagnostic
established that the traversal was far from complete.  The dense process was
therefore terminated cleanly, its nonterminal checkpoint was frozen only as
diagnostic evidence, and the equivalence-certified compact binary was started
from zero.

The production cutover followed this policy:

1. preserve the dense run until its exact prefix coverage and memory state are
   measured;
2. retain its 283-million-node progress only as non-certificate diagnostic
   evidence;
3. restart the compact traversal from zero, with no state splicing;
4. require normal terminal JSON and independent literal validation; and
5. never infer the full Record-31 `a=1` schedule from this one source orbit.

The compact traversal closed normally at 469,990,577 nodes.  Its exact eleven
ridge-valid incidence orbits all pass topology and all fail independent
unrestricted compact-label joins; see
`h9_14_record31_a1_orbit3_corrected_quotient.md`.  This is a complete
source-orbit-3 result, not a complete six-source schedule census.

A further fixed-width variant stores each at-most-22-byte compact code in
three inline 64-bit words.  Its generator, exact hash/equality, syntax-only C++
smoke, and lightweight validator tests are frozen.  A complete 231,773-node
source-orbit-0 replay gives byte-identical dense, compact-string, and inline
JSON transcripts, while both compressed builds enable the dense-adjacency
round-trip assertion.  The inline generated-source SHA-256 is
`ec577df47018f1d03d48e8304c56a2541ef02180144937600569ec0a9a7e68c7`.
That exact inline source subsequently completed production source orbit 5 at
41,488,368 nodes.  A separate direct binary-matrix quotient completed the
same source at 6,042,612 nodes; their terminal pure-target ridge censuses are
identical and contain no `{2:540}` target.  This is evidence for the production
representation and a complete source-orbit-5 incidence no-go, not a complete
six-source schedule theorem; see
`h9_14_record31_a1_orbit5_corrected_quotients.md`.

## Artifacts

- Compact generator:
  [`generate_h9_14_record31_a1_compact_key_engine.py`](../scripts/generate_h9_14_record31_a1_compact_key_engine.py)
- Equivalence verifier:
  [`verify_h9_14_record31_a1_compact_key_equivalence.py`](../scripts/verify_h9_14_record31_a1_compact_key_equivalence.py)
- Frozen equivalence certificate:
  [`h9_14_record31_a1_compact_key_equivalence.json`](../data/enumerations/h9_14_record31_a1_compact_key_equivalence.json)
- Regression test:
  [`test_h9_14_record31_a1_compact_key.py`](../tests/test_h9_14_record31_a1_compact_key.py)
- Fixed-width inline-key generator:
  [`generate_h9_14_record31_a1_inline_key_engine.py`](../scripts/generate_h9_14_record31_a1_inline_key_engine.py)
- Prepared triple-equivalence verifier:
  [`verify_h9_14_record31_a1_inline_key_equivalence.py`](../scripts/verify_h9_14_record31_a1_inline_key_equivalence.py)
- Frozen triple-equivalence certificate:
  [`h9_14_record31_a1_inline_key_equivalence.json`](../data/enumerations/h9_14_record31_a1_inline_key_equivalence.json)

Replay with:

```bash
PYTHONPATH=src:. python3 -m scripts.verify_h9_14_record31_a1_compact_key_equivalence --check-archive
PYTHONPATH=src:. pytest -q tests/test_h9_14_record31_a1_compact_key.py
```
