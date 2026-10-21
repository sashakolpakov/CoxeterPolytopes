# Source audit of the bounded D13 replay partition

20 October 2026. **PASS for source equivalence of the complete 19-batch
conjunction, with the effective assertion guard in the successor below.**
I found no omitted original output, input-status obligation, inherited
pin obligation, or mathematical recomputation. This is a source review,
not a statement that the 19 batches have all executed successfully.
No checker, project import, solver, Lean or TeX job was run. Small bounded
certificate metadata and file sizes were inspected; no scientific values
were recomputed. I previously helped assess the resource/projection
approach, so this is internal project review, not external refereeing.

## Exact sources

Final reviewed worker:
`scripts/check_d13_transfer_batches_20261020.py`, SHA256
`be8a42d1d3687204af4d559c6bd8a28779dbfac68aeccdd95d41a7172645d555`.

The full comparison was initially made against predecessor
`4ba0d5d7163953bd0f5c965a7126f62fd62a3d9708fb7ee2614ee159114f6f2a`.
Its byte-preserved payload is in
`manuscripts/source_archives/check_d13_transfer_batches_20261020_before_assertion_guard_20261020.py.gz`.
I checked the decompressed payload hash and that the final source differs
by exactly one replacement: the ineffective `assert __debug__` is replaced
by `if not __debug__: raise RuntimeError(...)`. This rejects optimized
worker execution before any mathematical definitions are imported or
replayed. No mathematical branch, selector or pin changed in that delta.

Immutable comparison boundary:

* Original producer `scripts/fresh_d13_three_rank15_context_transfer_20261018.py`,
  SHA256 `a95ea679854ba907ead8640964a258c3baeb3e998052a3e2e3ff20f230cce486`.
* Original certificate `results/fresh_d13_three_rank15_context_transfer_20261018.json`,
  SHA256 `2c7d615cd38123bef20cffad68de8ab30baa91826c147a6bb084892590e8ca86`.
  Its size is 1,116,554 bytes, below the worker's bounded read cap.
* Pinned wrapper `scripts/run_bounded_research_check_20261018.py`, SHA256
  `4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132`.

## The union covers the original run

The worker requires the exact 13 top-level keys and nine
`current_interfaces` keys of the immutable certificate. Every one has
an explicit batch owner:

| Batch or group | Original obligation replayed |
| --- | --- |
| `pins-0` through `pins-3` | Every one of the 161 input statuses and every inherited `sources` map |
| `metadata` | Status, all 484 outer source pins, ambient dimension, changed pairs, scope, and the four constant interface fields |
| `support` | Complete selected support and all selected proofs |
| `rank15` | Rank-15 contact classification |
| `one-sided` | Both one-sided whole models |
| `quiet` | All quiet repairs |
| `literal` | Complete literal corner rows |
| `later-0` through `later-3` | All later corner rows from the 22 original finite source entries |
| `newer` | The same eight newer proof collections in their original order |
| `variants` | All 24 D14 variants |
| `d12` | Both full positive-output walks and their D13 rows |
| `universal` | The complete original universal-output inventory |
| `finite` | The original final finite-count composition |

`chunks` checks that concatenating its four slices recovers the original
sequence and that the sequence is duplicate-free. For 161 stems its
slice sizes are 41, 41, 41, 38; for the 22 finite source entries they are
6, 6, 6, 4. Thus no source group is silently dropped. The metadata branch
also requires all 564 saved later rows to belong to those finite sources;
each group is compared with its complete saved sublist. The immutable
certificate fixes their original order. No assertion of completeness is
inferred merely from a count.

The 484-path outer domain is reconstructed as the original producer plus
the script/JSON/Markdown triples of the original 161 stems and compared
as a set with the certificate's source map. Every batch checks all these
hashes before and after its work. The four pin batches additionally
project every original input's status and optional inherited source map.
Optional absence of `sources` has exactly the original `get(...,{})`
meaning. An inherited pin outside the outer map is digested when first
encountered and checked again at the end of that pin batch. Repeated
inherited pins cannot conceal conflicting expected hashes.

All 19 batches are required. The worker itself prints that requirement;
it does not turn one successful named batch into a complete replay.
An execution summary must therefore bind all names, the same final
worker hash, the same frozen certificate and successful bounded runs.

## Projection does not narrow a mathematical domain

I compared the selectors with the actual bodies of the original routines,
including `T.literal_inputs`, `B.corners` and `S.W.walk`.

`support` loads **all** isolation cases and cores, then calls the unchanged
original `support` and `selected_proof`. It does not preload only the 33
successful affected cases. `rank15` and `one-sided` call the original
functions directly. `quiet` supplies both complete fields the original
function reads, including the independently named affine-B target.

The literal projection supplies exactly these eight original object
paths: the specified isolation core; the whole `exceptional_D15` and
`order17` objects; D13 and D14 equality rows; the separate D17 and D16
partial-transfer `cases`; and the full `two_A1_A3_spectators.rows`.
These cover every access of `T.literal_inputs`, including its inherited
finite/component assertions and all corner rows. The unchanged
`long_images` then checks the same marked whole-component role.

Each later source loads the complete top-level proof fields selected by
the original `FINITE` dictionary. Calling the original recursive
`corners` on each field with the same initial path preserves inherited
Gram/white data, dictionary/list traversal, multiplicity and path strings.
The eight newer collections retain the same physical-seal fallback and
all corners. The variant batch retains all `selected_proofs` and all
their variants, rather than selecting only saved hits.

For D12, both documents are deliberately decoded whole. The original
recursive `walk` therefore has exactly its original domain and root
paths; its 280/52 positive-domain tests and D13 filtering are unchanged.
The universal batch supplies all three fields actually read by
`universal_roles`. The final batch supplies the complete frozen eight
obstructions and exact prior maximum read by the original `finite`.

These latter calls retain their original scope: the final composition
recomputes the affected old obstructions and exact increasing families;
it is **not** a new unrestricted forest DP. Similarly, the original
universal inventory explicitly inherits written arbitrary-support
classification/routing results. Partitioning that inventory does not
prove those geometric theorems afresh. No stronger mathematical
completeness claim is justified by this worker.

`equal` uses the original encoder followed by the same JSON normalization
as the old command-line verification. Native tuples and exact-number
encodings consequently receive the same comparison semantics. The split
does not replace rational/inertia/embedding calculations with saved
booleans: it calls the original mathematical functions, or reproduces
the original interface loops with unchanged native helpers.

## Parser and resource semantics

The read-only mmap parser builds a trie of exact object paths, rejects
overlapping selectors and duplicate keys in visited objects, and treats
quoted strings as whole lexical tokens while scanning skipped composites.
Braces or brackets inside strings cannot change the nesting stack.
Every requested leaf is decoded by `json.loads`; required paths and the
end of the root document are checked. A root-empty selector decodes the
whole document, as needed for the two D12 walks.

The parser intentionally does not validate every skipped scalar's JSON
grammar. This is sound at the declared boundary: the entire input file
is pinned before and after each batch to the immutable original JSON,
and skipped fields do not enter that batch's calculation. This is an
input projection, not a general validator for arbitrary mutable JSON.
I did not rely on the separately reported parser-fixture PASS as proof
of these source properties, and did not execute or re-audit that fixture
checker here.

The outer certificate has a bounded snapshot and final byte comparison.
Scientific file hashes stream in fixed-size chunks. Unselected composites
are scanned without decoding or copying their arrays. The 161 current
JSON input sizes are below the 40 MiB file cap (maximum 32,888,237 bytes).
The decoded leaf cap is a raw-byte cap, not a bound on all Python/SymPy
allocations; the wrapper's actual memory/time watchdog remains necessary.
Source review does not establish that every mathematical batch fits its
59-second/1-GiB execution limits.

The final effective optimization guard resolves the only concrete
correction found in this audit. The parent-lock check, parent-liveness
checks, worker/input post-pins, and replay-only/no-certificate-write
behavior remain intact. No source-equivalence gap remains identified;
fresh successful execution of the full conjunction is a separate task.
