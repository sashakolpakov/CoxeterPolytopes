# Distance-nine JSON-key correction: independent source review

**PASS for the limited representation correction and dependent pin updates.** This is source review, not a mathematical replay or scanner-fixture execution. I previously reviewed the D13 projection helper and related receipts; I did not author these three corrections. No worker was imported or run, and no canonical source was edited during this review.

I decompressed the four `manuscripts/source_archives/*_before_json_key_lists_20261020.py.gz` predecessors, verified their decoded hashes, and compared their complete text with the current sources:

| Source under `scripts/` | Decoded predecessor SHA-256 | Current SHA-256 |
|---|---|---|
| `check_selective_distance9_batches_20261020.py` | `126d5fe3588cfbb56009bff1f748dbab2154582e3df814feff50b32003cf023b` | `f5f0b47084135f8bb1886a37e0478d1f6b8c9c4327499f986a8ec9e173ecd843` |
| `check_distance9_array_scanner_fixtures_20261020.py` | `1f0bf6d952d536a0eb7316c088d8dceb23c52687c8e7fdc2536547b5371cf09b` | `f6fd1267462f1c95632998125aef419f939803d08aaaf3ba307feaa273701ff2` |
| `check_partitioned_replay_coverage_20261020.py` | `41f3c82b86ad05aa0877a34cf4718f92d5253e72fe669c083dd161a22c4342c6` | `456abc2bfe92262412106ebd6f94b07e48d01757cd1c76d8cef952c8bf3d3b0e` |
| `resume_distance9_replay_20261020.py` | `1e3432ec083ce1fd6fce65f8b9ad50ca202e8680935eaa13c8cf6b97ebc4f6d1` | `cb776f5931e24f8541c9b8321123dde563c3dbb6cd1e3240cdceeb386e500541` |

The worker has exactly three one-line changes, at current lines382,394,408. Each changes a reconstructed expected key from a tuple to a list: ordered join `(left,right)`, multiple-edge `(left,right,edges)`, and overlap `(left,right,match)`. The unchanged `equal` function compares `json.loads(json.dumps(actual))` with the expected value. The predecessor therefore compared normalized list entries against tuple entries and necessarily rejected these nonempty domains even when their values agreed. The repaired expected lists match the representation of the saved JSON fields.

The original immutable producer `scripts/fresh_selective_distance9_dimension108_20261017.py`, SHA-256 `f5cd7d0181cf9d5f96b10b6d4d6e7654caa3369b18ffa7473a16033f3124c81a`, explicitly normalizes its entire result through the same JSON round-trip before verification. I compared its ordered port-product joins, sorted context-group/subset loops and ordered overlap/matching loops with the successor's corresponding domain generation. The correction changes no key value, ordering, multiplicity, sorting, filter, range, Gram calculation, exact inequality, witness, counter, normalization function or selected-field path. In particular, the4442 joins,1332 multiple-edge subsets,20718 consistent overlaps,85 inertia types and all phase bounds remain unchanged. This is equality in the intended JSON representation, not a tolerance or mathematical weakening.

The fixture and coverage checker each change only their worker-hash literal; the coordinator changes only its coverage-checker hash. All three new pins equal the actual corresponding source hashes. The producer, certificate, D13 helper, bounded wrapper and manuscript driver pins are unchanged. The explicit disabled-assertion rejection remains in the worker. No parser or coverage logic changed in this delta.

I read the author's receipt `results/review_distance9_expected_key_json_representation_fix_finite57_20261020.md`, SHA-256 `114e38930209f4bc23f965b7745d7ca000efe0d107859ad2b3625b3519b1e185`. Its description of the representation defect and exact four-file delta matches the inspected sources. Its metadata/self-test execution statements are not new executions by this reviewer. Old fixture results and old-source mathematical receipts do not establish current-source execution or complete partition coverage. Fresh bounded execution and actual coverage validation remain root-owned obligations; no verification status is advanced by this report.
