# Lorentz kernel: actual verification and final map review

21 October2026. **PASS for the bounded source/output/receipt bindings and
the stated conditional scope.** This is metadata/output inspection, not
a rerun of Lean, a compiler, tests or mathematical jobs.

The immutable manifest
`results/lean_lorentz_obstruction_verification_20261021/manifest.json`
matches SHA256
`63d8bf415a344c3551956d98b862673d66290861c76ff51dbb37e353a917db86`
(4617 bytes). Final Lean source5c0f8e57 remains unchanged (4765 bytes).
Its object matches the compile ledger at
`3f42079ce3022d40f5bcb11a65235c2363c87027a249ffe3b7b9522f391ad4b7`
(257696 bytes). The module imports Std only; the empty local dependency
map is not a separately checked transitive Std closure.

I read the actual1434-byte child output, SHA256
`479036b2bdb8a7be1df2ecc32d69628e8a7a01a921458fed8108b4f6bf2322d0`.
It has exactly the seven requested theorem names, each once. The first
six list only propext, Classical.choice and Quot.sound; the final Nat
rank contradiction lists only propext and Quot.sound. No successful
output has sorryAx. One unused-section-variable warning in the Schur
identity is retained. These are seven actual declarations, not seven
new geometric theorems or counts inferred from source commands.

Both immutable ledger bindings match:

* compile: `59060c1edd68b0654bee12de98ff3d93934c46c3a7154b5b538d24577fcebe34`,
  1522 bytes, terminal KERNEL_COMPILED,0.52s;
* separate audit: `9c6dec43417023a8a3cebe2113834cb90b152157467212ba41280480dce3e9d3`,
  1441 bytes, terminal KERNEL_AXIOM_AUDIT,0.51s.

They bind the same final source, compiler/toolchain inputs and actual
child-output digest, with2048MiB/one thread. I also checked both outer
stdout hashes against the manifest. The audit outer output's first1434
bytes equal the preserved child output exactly; its remaining ledger
and sequential-wrapper records report PASS. Initial free disk was
134.56GiB, available RAM7.54/7.57GiB respectively. These are admission
samples, not measured peak RSS. No successful dependency rerun is claimed.

The two failed predecessor sources and failed child/outer outputs remain
bound by the separate preservation receipts in the manifest. Their
recorded hashes match. Both have ZERO verification credit, regardless
of recovery axiom prints. The final static review chain is the complete
seven-statement review70a73b90, explicit-transitivity repair reviewd5a728d1,
and final source5c0/correct-name review920d2655. The last explicitly
corrects my earlier misqualified transitivity-API suggestion. The
actual terminal evidence above supplies elaboration credit that none of
those static reviews claimed.

At root's authorization I updated ONLY the logic map's execution status:
`docs/lorentz_obstruction_kernel_logic_20261021.md`, final SHA256
`47a68353036a316203e088ff42663d7e19920d78b5d96790a86bfde02be820ab`
(5061 bytes). Its prepared predecessor9a51f5b4 was byte-compared with its
gzip decode before editing; archive SHA256
`da61d8f6bcd65765d6775f1b3972ca0526a3917332a06ee7cb275a3d92815d5b`
is preserved at
`results/source_archive/lorentz_obstruction_logic_before_actual_verification_20261021.md.gz`.
Receipt `results/lorentz_obstruction_logic_actual_status_preservation_20261021.json`
matches `6d56422cc7218454392991c64ea0d7adb1fc2c66ff6566113e1605e7b0d6194b`.
The earlier a8242e69 map and its squared-norm wording correction also
remain preserved. The immutable verification manifest was not edited.

The final map accurately distinguishes internal ring/order/rank logic
from finite-order averaging, actual Cauchy realization, original-lattice
and full-chamber facts, covolume/index arguments, complete white support
and Gram-rank interpretation. Printed axioms do not discharge these
theorem parameters. No polytope construction, host exclusion, finite57
charge improvement or dimension theorem follows from this kernel alone.
