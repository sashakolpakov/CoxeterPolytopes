# Six-paper selected verifier modes: source audit

20 October2026. **No actual mathematical coverage hole from the ignored `expected_status` field was found for the current selected argv.** This is an internal source/control-flow audit, not a checker execution or a new mathematical proof audit. No scientific module was imported, no Lean/TeX/checker job was run, and no canonical file was edited. The compact20 and ideal4 entrypoints were also read by the two other agents; their findings were checked against the selected arguments. I previously contributed several sources and scope documents, so this is not external refereeing.

Bound catalogue: `manuscripts/catalogue.json`, SHA-256 `9e8199f39a1b7dbfcbabea8b9a3b201a05aa5b05fe7ab5ee8937c6a076e06a52`. The six entries are `finite57`, `compact20`, `ideal_ra4_volume_and_structure`, `right_angled_sharp`, `compact_farthest_bounds` and `compact_ra4`; the last has role `primary_paper` rather than `primary`. Together their full selections name48 distinct Python checker sources. Their aggregate source binding is `ef6b0cab8e8fcb93a52034b28038db907fc9110eafc4b89b4ad44e3dcfdd5cbb`: SHA-256 of the UTF-8 concatenation of `file_sha256 + "  " + repository_path + "\n"`, with those48 paths sorted lexicographically and deduplicated from the catalogue argv. This binds the inspected source snapshot, not a runtime result.

Shared sources:

| Path under `scripts/` | SHA-256 |
|---|---|
| `verify_manuscript.py` | `f6a0b39f981c41f3c7d75b184ee3d274e309616bcfa72b693a3b018bfc7394a1` |
| `run_lean_batches_20261018.py` | `97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c` |
| `run_bounded_research_check_20261018.py` | `4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132` |
| `build_illustrated_ra_manuscript_20261018.py` | `27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15` |

## Scientific argv and narrower intentional scopes

`verify_manuscript.py:46–110` builds the commands from the catalogue; it does not propagate or enforce `expected_status`. Lines139–155 classify a job from its exit code. Consequently the status field alone is no assurance against arbitrary future command changes. For this snapshot, the selected scientific commands have the following actual paths:

- **Finite57:** every selected producer receives `--verify`, performs its `run`/certificate checks and compares the result with saved output. For example, `fresh_d13_context_dimension57_20261018.py:151–155` calls `run()` before equality and success output. The historical tails' optional `--scalar-only` paths are not selected: see `fresh_ideal5_piecewise_dimension77_20261017.py:130–145` and `fresh_distance8_ideal_dimension97_20261017.py:116–124`. The saved-table path in `fresh_d17_a16_coupled_ambient58_frontier_20261018.py:118–148` still checks types/weights, every transition, tight predecessors and unreachable states; its within-call table cache is not a skip-verification mode. All19 selected D13 parts require `--verify --batch`, validate pins, invoke `replay`, and only then print their batch success (`check_d13_transfer_batches_20261020.py:376–414`). The parser-contract check deliberately tests projection semantics, not scientific facts (`check_d13_json_projection_contract_20261020.py:130–187`); it is separate from those19 parts.
- **Compact20:** the selected primal reconstruction's `--verify` path always reconstructs and verifies before comparing saved output (`fresh_compact21_self_hub_exact40_reconstruction_20261018.py:126–135`). Its exact maximum is asserted at most40 at53–54. The independent `--candidate` checker has no verification/generation switch: it checks the candidate pins,729 contexts, native faces and forests before printing (`check_compact21_self_hub_primal_independent_20261018.py:298–339`). A precise limitation is that its `primal_at_most_40` field is printed, not separately asserted there; it is a general exact-primal evaluator. The current other two checks enforce40, including `check_compact_arithmetic20_composition_20261018.py:99–106`. Thus this nuance does not permit a false complete compact20 selection PASS.
- **Ideal4:** both selected producer `--verify` paths call their full `build()` before saved equality (`fresh_ideal_ra4_double24_volume_20261018.py:172–178`, `fresh_ideal_ra4_equality_structure_20261018.py:215–223`). The mode-free saved checkers always validate before success (`check_ideal_ra4_second_volume_saved_20261018.py:162–196`, `check_ideal_ra4_equality_structure_saved_20261018.py:12–146`). The path-root checker always runs its finite calculation (`check_ideal_ra4_octahedral_path_roots_20261018.py:85–93`).
- **SharpRA:** the mode-free24-cell checker always constructs/checks its model (`check_ideal_right_angled_dimension4_model_20261020.py:23–96`). The SOS worker's required mutually exclusive modes select `--verify`, so `check()` and exact saved equality occur before its `PASS` (`check_right_angled_sharp_dimension_sos_20261020.py:128–150`). No write/stdout-only branch is selected.
- **Compact23:** the explicit `--manifest ...dimension23_packet...` avoids the verifier's older dimension24 default. `verify_compact_farthest_vertex_packet_20261020.py:14–62` checks packet bytes and prints `SOURCE_INTEGRITY_PASS` with `mathematical_or_lean_verification=False`. This is intentionally a metadata-tier integrity check, not mathematical verification. It enters the full plan, while core omits that tier (`verify_manuscript.py:48–52`).
- **CompactRA4:** `verification=[]` and `lean.modules=[]`. Its PDF-only plan therefore verifies typesetting only. Without a PDF or any declared job, the driver refuses to claim a pass (`verify_manuscript.py:202–203`).

The bounded wrapper invokes the selected script as `__main__` and propagates failure, rather than substituting a planner (`run_bounded_research_check_20261018.py:25–38,76–100`). This assessment concerns the current assertion-enabled Python commands and does not certify arbitrary altered interpreter environments or future argv.

## Exact legacy Lean qualification

Eight selected targets currently have existing objects, no per-module capped compilation ledger, no axiom printers, and objects not older than their sources:

```text
finite57: GuardedRank58, Concentration58, ContextCover59, ShapeCurvatureArithmetic
compact20: CompactHVector, CompactFieldBudget, DimensionComposition
ideal4: IdealRA4Gap
```

All names are under `CoxeterProofs`. These commands **can return zero without fresh elaboration**: `run_lean_batches_20261018.py:167–178` recognizes `LEGACY_CACHED_UNCAPPED`, and329–340 reuses it while explicitly reporting the legacy list and `all_dependencies_have_capped_ledger=false`. This is not a new capped compilation, a fresh source/object hash binding for those legacy objects, or an axiom audit. It must retain that qualification. It is separate from the five newly verified source/object-bound kernels and their44 parsed actual axiom prints.

For targets containing actual printer commands, the driver adds `--audit-axioms` (`verify_manuscript.py:79–89`). The runner always calls `batch(..., audit=True)` even when compilation is cached (`run_lean_batches_20261018.py:336–337`), so their audit is not a check-only/cache no-op. Nevertheless neither exit zero nor `KERNEL_AXIOM_AUDIT` parses/approves axiom sets automatically; the separate actual-output summaries remain necessary. The current five-kernel receipts supply that distinct validation. No broader end-to-end geometry or fresh legacy verification follows from the orchestration status.

## PDF and planning qualification

All six current entries use the default PDF builder. Its `--check-only` mode really invokes TeX in a temporary directory, requires successful output, checks the source is unchanged and hashes the resulting PDF (`build_illustrated_ra_manuscript_20261018.py:81–139`). It suppresses replacement of published artifacts at140–142; it is not a cache or no-build mode. Its success explicitly says `mathematical_verification_claimed=False` at137.

The driver's own `--list` and `--plan` return zero with `INVENTORY_ONLY` and `PLAN_ONLY_NOT_VERIFIED` respectively (`verify_manuscript.py:188–201`). No current selected checker argv calls those paths as a scientific verification job. Their separately saved metadata plans are not runtime proofs. No frozen driver, packet, script, certificate or manuscript was changed by this audit.
