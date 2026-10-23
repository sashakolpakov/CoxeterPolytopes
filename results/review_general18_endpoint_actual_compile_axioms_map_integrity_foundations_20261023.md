# General18 endpoint kernel: actual runtime, map and integrity review

23 October2026. **FULL actual-evidence/source-binding integrity PASS.**
I read both complete saved tool transcripts, both successful ledgers, the
entire2296-byte actual axiom payload, frozen runtime manifest and current
map. I checked all20 direct file bindings, current compiler/input/object
hashes, the exact map archive roundtrip/delta and owned-work absence.
No compiler, Lean, worker, mathematical import, enumeration, test, build
or Git operation was run; the checks here read/hash existing evidence only.
No bound source, map, receipt or previous review was modified.

The root authored the NEW Lean code. I independently reviewed that code
before execution in the source review below; I authored the underlying
written A18 analytic proof and contributed inherited endpoint foundations.
This is different-author runtime/code-binding review, not an independent
mathematical self-review of that written proof. Its separate outside/root
acceptances remain the mathematical authorities for that premise.

## Exact immutable authority chain

| Object | Bytes | SHA256 |
|---|---:|---|
|results/lean_general18_endpoint_refinement_verification_20261023/manifest.json|12055|9576ff6555c0eca49326af4b43e87f5f94606ed90f878417873948aaf318edd6|
|formal/CoxeterProofs/GeneralEighteenEndpointRefinement.lean|10668|f6630f0f3275295c71b1bf8cc49e8476dc022177fb918de006ac4e9316d831c5|
|docs/general_eighteen_endpoint_refinement_logic_20261023.md|6906|0166a9a14c57514954b994a5a154477f83070cdd84a378b36e1c1baa35cd5deb|
|results/review_general_eighteen_endpoint_refinement_full_source_foundations_20261023.md|10865|c456ef88ce7208b7e61666e5bcd2d456f60fdc45bc8d427ddee970eb31d40eb4|
|results/review_general_eighteen_endpoint_refinement_precompile_root_gate_20261023.md|2942|61a2176e5263fc73d403fe40348a43049b2204c1bc51495bbaf9dfa356dc421b|

Every one of the manifest's20 path/byte/SHA256 bindings matches current
bytes, including all written-source/review/acceptance authorities,
the corrected A17-I source0087, old d45f source, both transcripts/ledgers
and the pre-execution map archive/manifest. No stale successor, omitted
current hash or source replacement was found. The source remains exactly
the f663 bytes reviewed BEFORE these two jobs.

## Separate actual jobs and input identity

The successful compile ledger is1648bytes, SHA
e42d95b7ab0ac6bc77742895f842340699dc63168dd3880beb258d469954454a.
Its complete JSON equals the manifest's `compile` record. It reports
`KERNEL_COMPILED` for the endpoint module in0.74s,2048MiB/one thread,
with101.68GiB free disk and7.79GiB available RAM at admission.

The SEPARATE audit ledger is1567bytes, SHA
b3db1d74d37c648e9699e341205c1464e2b130af70691b0d36e880c6084dac28.
Its complete JSON equals the manifest's `separate_axiom_audit` record.
It reports `KERNEL_AXIOM_AUDIT` for the SAME source, dependency and compiler
inputs in0.51s, with101.68GiB free disk and7.76GiB available RAM.
These are admission observations and elapsed times, not a peak-RSS claim.

Both saved operational ledgers are byte-identical to their immutable
copies in the verification folder. Both exact compiler/input snapshots
match current bytes: the formal lean-toolchain/lakefile/lake-manifest and
all FOUR named compiler binary/shared-library files. These files were
hashed, never executed. Both use the same recorded empty override
environment and the unchanged2048MiB/one-thread settings.

The compile tool transcript is1851bytes, SHA
a041fcdcd1ea0080b976e8540bcffe72134e22caa7561c73f19583e4e1263e7d.
Its JSON ledger matches the successful compile record and its terminal
`PASS_SEQUENTIAL_MODULE_CHECK` lists ONLY the endpoint module as newly
capped, with no legacy uncapped dependencies. The audit tool transcript
is4017bytes, SHA
db3f1bd2fdfefa2e0bc6a9b89756f7a2b1ff84c9e3dfae8ddf471461722ecf76.
Its terminal PASS has an EMPTY newly-capped list and no legacy uncapped
module. It retains the raw actual payload followed by its matching audit
ledger. Thus the audit is separate, rather than compilation prints being
silently reused as a claimed second job.

## All eleven actual declarations and warnings

The complete actual payload is2296bytes, SHA
9e24915e908b562301781d4a77b80546fc77ecacbcf24a0d96ab5cd3838e77fe.
It is an EXACT byte prefix of the actual audit tool transcript. BOTH
successful ledgers report this same byte count/hash. Parsing the actual
quoted declaration names and multiline axiom lists gives EXACTLY the
eleven requested source names, in order, with no missing/extra/duplicate
request and no mismatch against the manifest's declaration summary:

| Declaration under CoxeterProofs.GeneralEighteenEndpointRefinement | Actual axioms |
|---|---|
|exact_five_cases|none|
|a17_all_spectators_exclusion|none|
|classified_five_case_refinement|propext|
|same_vertex_five_case_refinement|propext|
|Arithmetic.all_angle_integer_seals|none|
|Arithmetic.machin_integer_seals|none|
|Arithmetic.radical_integer_seals|none|
|Arithmetic.mutual_three_integer_signs|none|
|Arithmetic.mutual_four_endrange_signs|none|
|Arithmetic.all_mutual_four_gap_seals|none|
|Arithmetic.consecutive_gap_excludes_root|propext, Quot.sound|

There are eight axiom-free outputs, two using only propext, and one using
propext/Quot.sound. Neither sorryAx nor Classical.choice appears. The raw
payload retains EXACTLY five unnecessarySimpa warnings, at source lines
72,79,86,93,100. No zero-warning claim, warning suppression, source polish
or unnecessary successful-job rerun is inferred or needed.

## Dependency, object, archive and cleanup

The imported source remains d45f3a3652e82b8fb908e2b58857cbf2f18b8b1acb255c5c2309296ab2656194.
Its current capped object SHA35ae12d6380bfae171e8a58af83a89f8124620cf11c4f6775871ba79e174da41
matches both job inputs and the prior successful six-group runtime manifest.
The transcripts show it was reused without rerunning its compile or audit.
The signal-safe/base runner bytes still match their earlier bound pins
c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d and
97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c.

The new endpoint .olean is currently759480bytes, SHA
061b448e3a8e705def17aeb83e5999acd42528b019a2de3eb220fa632ff55f13,
matching the compile record and cache-only declaration. It is under the
explicit /formal/.lake/ ignore rule; .olean is also globally ignored.
The owned formal/.proof-batches/work directory is ABSENT. This audit did
not recreate work or delete evidence. An ignored object is a local cache,
not a promised public source asset.

The archive manifest b744ef3f66a2369b7bdad2727dde111f8bd1d9102c81354719ec55a98afcdc97
binds the2806-byte gzip cfe1d23027dc33a9eaea12bbcc49f35702578478f14c104b5276899c53598059.
Decompressing its saved bytes gives EXACTLY6124bytes and SHA
39810d34cb886504cee320643f86344208a082294047991c05ba324e6adedd5e,
the preparation map bound by the full source review/root precompile gate.
I compared that complete predecessor with the current6906-byte map. ONLY
the opening execution-status paragraph and the final execution/publication
section changed. The geometric premise list, integer meanings, formal
coverage, external bridges and fallback boundary are unchanged. The new
paragraphs agree with the actual two ledgers and all eleven raw outputs.

## Exact accepted formal scope

Actual compilation and a separate actual standard-only axiom audit are
established for this conditional kernel. It derives the five forest cases
at the SAME supplied selected vertex from nine inherited plus two new
EXPLICIT written premises. It also proves the thirty Taylor integer seals,
finite Machin/radical/sign certificates and discrete gap logic reviewed
before execution. Explicit theorem parameters are not discharged by a
standard-only printed axiom list.

The code has no real trigonometric or hyperbolic instance. Machin identity,
analytic remainder/monotonicity/sign transport, actual walls/Gram/incidence,
own/native-label completeness and all geometric exclusions remain written
application inputs. The written A18 theorem covers ONLY two actual endpoint
opposites that are both isolated WITHIN whole A18 and have finite own
contact. The kernel neither excludes whole A18 nor proves global compact17
or any new example. Ceilings57/18/16 and the six frozen papers stay unchanged.
The cyclotomic fallback remains UNEXECUTED and unused. No integrity or scope
correction is requested; the manifest and final map are suitable for the
stated bounded conditional verification claim.

