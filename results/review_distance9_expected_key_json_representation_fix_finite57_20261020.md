# Distance-nine expected-key representation correction

Source/type review and authorized metadata self-tests only. No mathematical
worker was imported or run, and no mathematical replay coverage is claimed.
The author of the replay successor performed this review; it is not an
outside-author audit. Root independently reviewed the three-line repair.

## Actual failure and exact cause

The actual root-owned runtime ledger
`results/finite57_distance9_partitioned_replay_runtime_20261020.json` records
`distance9-domains-0-1-5`, FAILED after2.709s, worker SHA256
`126d5fe3588cfbb56009bff1f748dbab2154582e3df814feff50b32003cf023b`.
The traceback ends at line382, `ordered join domain`, through line92:
`json.loads(json.dumps(actual)) == expected`.

There are4442 projected join keys. Each actual key is constructed as a
tuple of two lists, then JSON-normalized to a list of two lists. Each
expected key was reconstructed as a tuple of two saved lists. Thus this
comparison necessarily failed for a nonempty domain regardless of the
key values: JSON lists do not equal Python tuples. The same construction
error occurred in the multiple-edge and overlap projected-key comparisons.
This diagnosis does not assert that an unexecuted numerical comparison
passed; root must replay the successor.

The repair changes only the expected-entry delimiters at lines382,394,408:

```python
[[row['left'], row['right']] for row in expected]
[[row['left'], row['right'], row['edges']] for row in expected]
[[row['left'], row['right'], row['match']] for row in expected]
```

The one-sided `equal` contract remains unchanged: canonicalized recomputed
output is compared against JSON-shaped saved/projected expected output.
No numeric conversion, value tolerance, sorting, filtering, domain bound,
witness check, counter, output order or coverage interval changed.

## Complete comparison inventory

There are20 call sites, excluding the definition of `equal`.

| Lines | Comparison | Expected representation |
|---|---|---|
|306|metadata constants|dict of values loaded by `saved`|
|330|isolation global scalars|dict of values loaded by `saved`|
|346|all overlap counters|saved JSON dict|
|358|complete profile domain/order|explicit JSON round-trip already present|
|362|boundary|saved JSON object|
|370|facet polynomial|saved JSON array|
|374|complete port census|saved JSON array|
|382|ordered join keys|**repaired nested lists**|
|386|join counters|dict of values loaded by `saved`|
|394|multi-edge subset keys/order|**repaired nested lists**|
|397|sorted inertia types|list of saved `ports` arrays; already correct|
|408|ordered overlap keys|**repaired nested lists**|
|436|single-join witnesses|raw selected JSON rows|
|453|multi-edge witnesses|raw selected JSON rows|
|468|exact inertia witnesses|raw selected JSON rows|
|503|overlap classifications|raw selected JSON rows|
|531|isolation core|saved JSON object|
|557|isolation cases|raw selected JSON rows|
|575|reduced profiles|raw selected JSON rows|
|578|dimension minima|raw selected JSON rows|

The other direct comparisons were also inspected: metadata's tuple of four
saved fields is compared with another tuple of four saved fields; context
ownership uses `(diagram,tuple(black))` on both sides; core-cache keys are
strings built from the same tuple-converted match; fixed source maps have
string keys; `old.INPUTS` and `stems` are both tuples. Scanner token/index,
set/counter and exact scalar comparisons introduce no analogous list/tuple
mismatch. No additional representation repair was identified.

Source correspondence was checked against immutable producer
`scripts/fresh_selective_distance9_dimension108_20261017.py`
(`f5cd7d0181cf9d5f96b10b6d4d6e7654caa3369b18ffa7473a16033f3124c81a`).
The successor reuses its context/port generator. Join ownership retains
the same ordered eligible-port product and distance test. Its supporting
white inverse is the same rationalized matrix inverse, cached per context
rather than per port; the exact rational tests `<1` and `>=1` are
complements. Sorted context groups, edge order and all subset sizes match
the producer. The overlap generator has the same ordered context product,
matching sizes and injectivity/Gram-consistency pruning, restricted only
by the already-reviewed left-context partition. The three fixes change
none of these operations. Saved fields are JSON lists, as illustrated by
the first `disjoint_obstructions.rows` entry of the frozen certificate.

## Source bindings and preservation

| Source | Predecessor SHA256 | Held successor SHA256 |
|---|---|---|
|worker|`126d5fe3588cfbb56009bff1f748dbab2154582e3df814feff50b32003cf023b`|`f5f0b47084135f8bb1886a37e0478d1f6b8c9c4327499f986a8ec9e173ecd843`|
|scanner fixture|`1f0bf6d952d536a0eb7316c088d8dceb23c52687c8e7fdc2536547b5371cf09b`|`f6fd1267462f1c95632998125aef419f939803d08aaaf3ba307feaa273701ff2`|
|coverage checker|`41f3c82b86ad05aa0877a34cf4718f92d5253e72fe669c083dd161a22c4342c6`|`456abc2bfe92262412106ebd6f94b07e48d01757cd1c76d8cef952c8bf3d3b0e`|
|resume coordinator|`1e3432ec083ce1fd6fce65f8b9ad50ca202e8680935eaa13c8cf6b97ebc4f6d1`|`cb776f5931e24f8541c9b8321123dde563c3dbb6cd1e3240cdceeb386e500541`|

The fixture/coverage edits replace only their worker pin; the coordinator
edit replaces only its coverage-checker pin. All four predecessors were
saved before editing in `manuscripts/source_archives/`, using the original
stem followed by `_before_json_key_lists_20261020.py.gz`; decompressed
bytes were verified identical to the original and the hashes above.
The frozen producer/certificate, all runtime rows, and the prior Boolean
repair receipt remain unchanged. That prior receipt is
`results/review_distance9_symbolic_boolean_count_fix_finite57_20261020.md`,
SHA256 `3e9b97f2a0f09072e45deb5cc0dd4b74967f98b0ec6e59e550b9aad73d9f736f`.

Authorized metadata-only runs actually completed: coverage `--self-test`
(synthetic oracle fixtures, no replay credit), coordinator `--self-test`,
and coordinator `--plan` against the current runtime ledger. All exited0.
The plan credits zero current-source units, ignores all6 historical rows
(2 stale successful workers and4 non-successful records), and proposes
`metadata [0,1)` first. The complete2196-unit domain remains unverified
under the new worker. The ledger's before/after SHA256 was unchanged:
`28e8d04e47a57b348ce802b8c0d22aa99aec82d84dd417ebd66ae85743c1db95`.
Only root may run the resource-gated successor after source review.
