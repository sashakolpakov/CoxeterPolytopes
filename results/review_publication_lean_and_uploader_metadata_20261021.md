# Publication Lean evidence and serial uploader: source/metadata review

21 October 2026. **PASS at the conditional-kernel and trusted-publication-input boundaries below.** I read the saved ledgers/stdout, compared their byte bindings and current small sources/objects, and read the complete uploader and watchdog. I ran no Lean, mathematical worker, test, Git command, network request, uploader or document build. This review is metadata/source inspection, not a reproduction of the recorded executions or a claim of remote availability.

## 1. Two actual conditional Lean verifications

The immutable authority is `results/lean_publication_completion_20261021/manifest.json`, **7162B**, SHA256 `593f5f4aa2dbae6e3c790cb07208099c9c497a966cec88587fae08e5dd2b9419`. All eight referenced ledger/stdout files match its stated bytes and SHA256. The sources and resulting objects also match current bytes:

| Module | Source SHA256 | Object SHA256 |
|---|---|---|
| RegularIdeal24CellHeightKernel | `633ae6fef8c5d0f0d0cafb2363d99f38e6922405d61c25c0e3a569e7bceaecec` | `29d8602bbd4d0f1285c5366ed6fde5d6186f8550fff941b0492b81b19b442dbf` |
| ClassicalDimensionComposition | `c035f62de642c123ea168a61871bc6d3156113343b3f5301a877844f0062e102` | `2c9fc7aac4a8136d0a42037dd0c751ff2571dcdcd03d4bea0b5b5802ae8961d0` |

The four individual records, under the manifest directory, are:

| Run | Ledger: bytes / SHA256 | Full captured stdout: bytes / SHA256 |
|---|---|---|
| Regular24 compile |1871 / `ffff22a543bb020c8a93c53d270c0cff1e0d3346997a50c7c37d81a273f3623c`|2071 / `5766abeb00dda04f953f14d4dd196e4a4c9aada4993fa2d9eba90f4c6e05aae5`|
| Regular24 separate axiom audit |1790 / `10f1114de18ca2c226f8cb7c8598adb47051be0357fbe5a5e14b551b3f7511fd`|3784 / `3f4277854e86b2c0b8d6f82e34547f447d19653228d6a1bece60d89b019fafb7`|
| Classical compile |1526 / `23880fb4e944a908f596c7cdf37b1900dffc8389ef218607e06c1e3fcb0f0db8`|1725 / `575b4c52b168f4682c1b3b56ec81b51e4416d4b9366c4d794cb8a752251d4083`|
| Classical separate axiom audit |1445 / `611b2445c02736e5a2f63fd7b57bfc216fada1882bcadd9b11134c1d2ec47578`|2078 / `35785938a3bc7baeee251141469cbe7d915897564013389616c3c8583de3ad5a`|

I parsed exactly the five requested Regular24 declaration names and four Classical names from the actual axiom output, and matched them to the source print commands and manifest lists. Regular24 uses only `propext`, `Classical.choice`, `Quot.sound`; Classical uses only `propext`, `Quot.sound`. No `sorryAx` appears. Regular24 retains two unused-section-variable warnings; Classical has none.

The actual child-output portions have1840B/SHA `cd4c3a37207ceda207291104de50b63a57dd4e56ea3232dc7b1ceab7c4f3c23a` and479B/SHA `2aeb78620e8d3c268dfd878665a79762332f1a609985426c66fb3610890ac6b3`, respectively. Hashing those exact prefixes of the audit captures reproduces BOTH corresponding compile/audit ledgers' child-output digests. These child digests correctly differ from the full captures containing runner JSON.

The three recorded Regular24 dependency objects also match: RightAngledSharpDimensions `f2bd3858fee23605e95212d21784afcc8f6a66f148e11f524ea3cbeb077d7f91`; RightAngledCubicScalarKernel `8d06c91e92439233c0469782bd382f5cb87fc7056488800e40298bf4e00c781c`; RightAngledProjectionHeightKernel `b2e265c8445767735f7ea149cfcf8344df074229bf62bafe158b0b482d7ef2a3`. Classical's local dependency map is empty. Current toolchain/lake configuration hashes match the ledgers; this review did not independently rehash the installed compiler libraries or transitive Std closure.

All four records are successful compile/audit records with2048MiB/one thread. Recorded admission RAM is8.89–8.94GiB and disk132.42–132.43GiB, above unchanged40/6GiB start gates. Compile summaries name only the new module; audits name no newly compiled module, no legacy uncapped target, and capped dependency provenance. The receipt records strictly serial execution with no successful dependency reruns. Elapsed times are1.13/.76s for Regular24 and.34/.34s for Classical.

**Scope:** nine actual conditional internal declarations, not nine new geometric theorems. Regular24's metric-model transport and actual-cell interpretation remain inputs. Classical's forest, field, geometric, script and written branch interfaces remain premises; only its stated composition and integer arithmetic are formalized. Printed foundational axioms do not enumerate these theorem parameters. The six paper editions and Lean sources are unchanged.

## 2. Serial release uploader and watchdog

Reviewed final `tools/upload_public_data.py`, SHA256 **`fbb9d14925f97ebc86759645725f08ceab6bb2dac8ca1d1e57555b5a4274abf1`**, and `tools/run_publication_io.py`, SHA256 **`88d221c48b4a32282661d2551d1a7721b8d89bdcc34581d553fd89fa27831e43`**. I did not author the uploader or watchdog. My review found and root corrected a genuine tag guard defect before execution: checking the release's `target_commitish` does not establish where its actual existing tag points. Final code resolves `repos/{repo}/commits/{tag}` and compares that SHA to local HEAD, retaining draft rejection. Reversing only that text correction reproduces predecessor SHA `3b622edde8a0580cf992c2d33b52d1db913a580ee0a1ebc770a34e4ad9a34335`.

The uploader fixes the manifest repository/tag to `sashakolpakov/CoxeterPolytopes` / `research-data-2026-10-21`, requires its publication-I/O guard, and requires an existing published release at the actual current commit. It handles the33 manifest archives (2,178,662,957B) serially. Each local archive is size/SHA checked using1MiB reads before credit or upload; an existing remote asset is rechecked, not trusted from a previous progress report. Remote size, GitHub SHA256 digest and `uploaded` state must match. Duplicate assets, incomplete listings and mismatches fail. The upload command has no `--clobber` or deletion; a conflicting asset is not replaced. Local progress reports are owned temporary publication records, not immutable mathematical receipts.

Removing only the watchdog's uploader allowlist line reproduces reviewed predecessor SHA `d683e0577b428e768f860586cca923337d0d0c9981c1c1cd5d199623ceff858d`. Monitoring remains sampled256MiB owned-process-group RSS, one-second polling and publication floors3GiB RAM/20GiB disk. These are publication-I/O controls, not lowered mathematical proof admission gates or an allocation cap excluding subsecond peaks.

**Scope:** trusted frozen manifest/local archives and no malicious concurrent filesystem or release mutation. This source review certifies no completed upload, public restore, remote release availability, redistribution right, mathematical enumeration, or end-to-end geometric formalization. Actual publication and its remote verification remain root's separate actions.
