# Reproducibility and Lean companion: dimension bounds 57/17/15

This is the companion to the single consolidated
[dimension-bound article](../../arxiv/hyperbolic_coxeter_dimension_bounds/main.pdf).
Use it inside a full checkout of this repository. It connects the article's
ideas and named propositions to complete written proofs, exact finite
certificates, recorded executions and conditional Lean kernels. It creates
no additional manuscript family. The four arithmetic certificate replays
below require no bulk-data download.

The [machine-readable manifest](manifest.json) gives exact source and evidence
hashes, all thirteen catalogue items, the four explicit arithmetic input
snapshots, and the selected Lean source closure. The
[proof catalogue](../../arxiv/hyperbolic_coxeter_dimension_bounds/ADDENDA.md)
states the mathematical domains in full. The complete
[551-page supplement](../../manuscripts/addenda/dimension_bounds_20261030/main.pdf)
and its [source](../../manuscripts/addenda/dimension_bounds_20261030/main.tex)
are preserved byte for byte from the existing proof edition at
[`36eeaddf`](https://github.com/sashakolpakov/CoxeterPolytopes/tree/36eeaddf646eff25d4990c9ab730a20ada21f80f).

## Statement-to-proof map

Each item has its article labels, complete-proof labels, checker/input/receipt
links and applicable Lean modules in `manifest.json`. A written proof has no
numerical execution requirement unless a finite output is explicitly identified.

| Item | Mathematical obligation | Complete proof labels and checking route |
|---|---|---|
| FV1 | Intrinsic triangle and quadrilateral payments | `fv:thm:local`, `fv:app:routing`, `fv:app:weights`; current-rule application conjunction |
| FV2 | Complete necessary spherical-forest maximum 377 | `fv:prop:finite`, `fv:app:finite-contexts`; finite endpoint, producer-replacement and Bellman records |
| FV3 | All 2552 terminal profiles and cusp surplus 36 | `fv:prop:terminal`, `fv:lem:concentration`, `fv:lem:low8`; exact endpoint/low8 records and concentration kernel |
| FV4 | Independent exclusions in every higher dimension | `fv:app:tail`, `fv:app:narrow-tail`; recorded bounded-tail checks, scalar kernel, written unbounded tail |
| C1 | Necessary rank-18 cover and proper-row exclusions | `compact:sec:components`, `compact:sec:forests`; written same-vector component and row proofs |
| C2 | Whole D18 and the separate higher-rank tail | `compact:sec:d18-setup`, `compact:sec:composition`, `compact:app:tail:completeconjunction`; written cases plus two compact-dimension kernels |
| A0 | Actual-field and uniformly funded component estimates | `arith:app:arithmeticchain`, golden-field/packing and summed-foot appendices; written proofs and authenticated literature input |
| A1 | Empty quadratic A17 endpoint domain | `arith:prop:pair-output`(i), `arith:sec:algorithms`; exact quadratic certificate replay |
| A2 | Empty cubic A17 finite-anchor domains | `arith:prop:pair-output`(ii), `arith:sec:algorithms`; exact cubic certificate replay |
| A3 | Two necessary quadratic B/D17 pairs | `arith:prop:pair-output`(iii); source-bound batches, saved-record combine and distinct B17/D17 geometric consumers |
| A4 | One necessary quadratic whole-B16 pair | `arith:prop:pair-output`(iv); exact quadratic certificate replay and written geometric consumer |
| A5 | Complete rank-17/rank-16 reductions to whole D16 | `arith:prop:f17`, `arith:prop:f16`, `arith:prop:d16reduce`, `arith:lem:D16initial`; written covers and actual-field reduction |
| A6 | Manual dyadic exclusion and arithmetic bound 15 | `arith:sec:odd-inert`, `arith:sec:ramified`, `arith:sec:split`, `arith:prop:noD16`; written proof and dyadic kernel |

Necessary forests and surviving algebraic pairs are not asserted to be
realized polytopes. Geometric domain inclusions and consumers are supplied
by the complete supplement; a computational certificate does not replace them.

## Check this companion's metadata

Run from the repository root. This checks selected file hashes, Git index
membership, article and complete-proof labels, navigation links, and the
literal input closure of the four arithmetic snapshots. It imports no scientific
worker and runs no enumeration, Lean compiler or TeX process.

The guarded entry points require Python 3.11 or later on macOS, with
`vm_stat`, `ps`, `sysctl` and permission to inspect the owned process tree.
Create the ignored runtime directory once in a fresh checkout:

```sh
mkdir -p formal/.proof-batches
```

```sh
python3 -B scripts/run_proof_batch_20261024.py check \
  --ledger results/local_dimension_companion_metadata_watchdog.json -- \
  scripts/check_dimension_bounds_companion_metadata_20261104.py \
  --output results/local_dimension_companion_metadata.json
```

Choose fresh output names. `PASS_COMPANION_METADATA_ONLY` establishes those
metadata checks; it is not a new mathematical verification. Current article
and guide files are navigation links, while the preserved supplement and
selected historical evidence have exact byte bindings.

## Replay the four finite arithmetic outputs

The four workers use Python's standard library only. Their original source,
literal pinned inputs, saved parts, successful execution ledgers and stdout
are public repository files. No scientific worker is modified by this companion.

| Selection | Recorded domain/output | What its combine replay recomputes |
|---|---|---|
| `A1` | 25 fields, 53 models, 159 ordered pairs; empty | Every field/model/pair, including all rejections, using the original exact algorithm |
| `A2` | 2432 coefficient triples, nine strip models; empty | Every record in both full integer-ring coefficient boxes, using the original exact algorithm |
| `A3` | 36 fields, 120 pairs; two survivors | All field records, counters and surviving pairs; **the 118 rejected pair tests retain their original batch evidence** |
| `A4` | 43 fields, 273 pairs; one survivor | Every field and pair, including rejected pairs, using the original exact algorithm |

The runtime manifests for
[A1](../../results/a17_two_isolated_quadratic_endpoints_runtime_manifest_20261022.json),
[A2](../../results/a17_finite_cubic_endpoint_runtime_manifest_20261022.json),
[A3](../../results/d17_one_isolated_quadratic_fork_runtime_manifest_20261022.json)
and [A4](../../results/b16_isolated_quadratic_pair_census_verification_20261023.json)
retain their distinct original scopes. Same-code replay is not an independent
implementation. The four small A4 certificates published with this companion
are the exact existing outputs bound by the 23 October receipt, not new runs.

Inspect a snapshot without running mathematics:

```sh
python3 -B scripts/run_proof_batch_20261024.py check \
  --ledger results/local_A1_snapshot_plan.json -- \
  scripts/replay_dimension_bounds_arithmetic_20261104.py A1 --plan
```

Use `A2`, `A3` or `A4` with a new ledger for the other plans. For an actual
replay, remove `--plan` and supply fresh certificate and receipt names:

```sh
python3 -B scripts/run_proof_batch_20261024.py check \
  --ledger results/local_A1_replay_watchdog.json -- \
  scripts/replay_dimension_bounds_arithmetic_20261104.py A1 \
  --output results/local_A1_replay.jsonl \
  --receipt results/local_A1_replay.json
```

Run one selection at a time. The launcher copies only the named small inputs
into an owned temporary directory, executes the byte-identical frozen worker,
preserves the fresh result and removes its own temporary files. The current
outer watchdog and shared exclusive lock govern the entire process tree.
The recorded runner environment is macOS. The resource/process monitors use
macOS APIs and need permission to inspect their owned process tree; the
standard-library arithmetic workers alone do not make those monitors portable.
Install Python 3.11 or later and the pinned Lean toolchain before using their respective
commands. Follow the current repository resource policy on that host.

For the established direct bounded-runner route, use the same worker arguments
after this prefix (with a new ledger):

```sh
python3 -B scripts/run_bounded_research_check_20261018.py \
  --parent-lock --ledger results/local_A1_direct_watchdog.json -- \
  scripts/replay_dimension_bounds_arithmetic_20261104.py A1 --plan
```

This is the existing monitored runner, with the current policy and shared
lock; it is not the archived runner staged below. Execution approval and
working process monitoring remain required in restricted environments.

Three original workers hash a historical Lean-runner source as an input. Its
25 October resource-policy update changed that live file. The launcher stages
the [exact archived preimage](../../results/source_archive/before_user_percentage_resource_policy_20261025/run_lean_batches_20261018.py)
at the expected **temporary input path**. The worker only hashes these bytes:
the archived runner is never executed or imported, and the live resource policy
is never replaced. The manifest records this mapping explicitly.

## Lean source, premises and actual audits

The environment is pinned by [lean-toolchain](../../formal/lean-toolchain)
to `leanprover/lean4:v4.32.2` and by
[lakefile.toml](../../formal/lakefile.toml). The selected sources and their
local Lean imports are listed in the manifest. `Std` belongs to the pinned
Lean distribution. Broader historical Python checks retain their separate
[requirements file](../requirements-proof-20261018.txt).

| Argument | Selected actual evidence | What remains a mathematical input |
|---|---|---|
| Finite-volume bound 57 | Nine inherited selected roots; source/receipt entries and [logic map](../../docs/publication_logic_formalization_20261020.md) distinguish their scopes | Actual face/cusp geometry, geometric interpretation of finite domains, numerical/certificate premises and analytic comparison |
| Compact bound 17 | [Two modules, separate compilation and actual audits, ten declarations](../../results/compact17_targeted_lean_completion_20261025/manifest.json) | Analytic Taylor/value/derivative envelopes, geometric inventories, funded references, branch exclusions and the prior global bound 18 |
| Earlier arithmetic bound 16 | [Three separately audited declarations](../../results/lean_arithmetic16_composition_verification_20261022/manifest.json) | Actual arithmetic and field reductions, covers and finite outputs at that recorded interface |
| Arithmetic bound 15 | [Fourteen separately audited declarations](../../results/lean_quadratic_dyadic_triple_verification_20261023/manifest.json) | Actual number fields/valuations, norm and square conditions, geometric exclusions; the final integer deduction assumes the prior bound 16 and exclusion of dimension 16 |

These are conditional internal kernels. Their ordinary theorem parameters
are explained in the [compact-bound premise map](../../docs/compact_dimension17_conditional_logic_map_20261024.md)
and [dyadic premise map](../../docs/quadratic_dyadic_triple_kernel_logic_20261023.md).
`#print axioms` does not list those ordinary parameters. The package supplies
written proofs for the surrounding mathematics; it does not claim an end-to-end
Lean formalization of hyperbolic polytopes or number fields.

Some frozen sources still say “uncompiled” or “unexecuted” in preparation-time
comments. Subsequent exact-source execution receipts above establish their
later status. Neither those comments nor the recorded evidence is rewritten.

If reproducing a selected Lean module, compile it and then run a **separate**
actual axiom audit, with distinct new ledgers:

```sh
python3 -B scripts/run_proof_batch_20261024.py lean-compile \
  CoxeterProofs.CompactDimensionSeventeenComposition \
  --ledger results/local_compact17_compile.json
python3 -B scripts/run_proof_batch_20261024.py lean-audit \
  CoxeterProofs.CompactDimensionSeventeenComposition \
  --ledger results/local_compact17_axiom_audit.json
```

The other two direct `Std` targets are
`CoxeterProofs.D18UniformDottedHeadScalarKernel` and
`CoxeterProofs.QuadraticDyadicTripleKernel`. Successful unchanged modules need
no rerun merely to update this guide. Follow the [formal guide](../../formal/README.md)
for the inherited finite-volume targets and their dependency evidence.

## Finite-volume bound 57 and the bulk-input boundary

The [finite57 reproduction guide](../../docs/finite_volume57_reproduce.md)
documents the existing current-rule conjunction: 259 actual scoped units,
64 shards, 28,028 maps, 12,518 nonunit applications and 83 prefixes. The
companion binds its input manifests and selected direct receipt/source links.
This metadata inventory does not establish the complete dynamic input closure
of that reader or replay all native enumerations. Its explicit inherited
geometric and finite-data premises remain in force.

The primary `manuscripts/catalogue.json` deliberately keeps historical edition
selectors: `compact_farthest_bounds` selects 23 and `compact20` selects 20.
They are not commands for the current 17/15 arguments. This companion maps the
current claims without changing those historical selections.

The separate [bulk-data manifest](../public-data-manifest.json) describes
936 files in 33 prepared archives. **No bulk asset has been uploaded.** Its
restoration mechanism does not currently provide a public full regeneration
package. The small A1–A4 snapshots here are explicit public inputs, separate
from that collection. Older source records also refer to superseded guide
bytes or retired PDF locations; the manifest lists those unresolved historical
links as provenance, not prerequisites for these small arithmetic replays.

All actual jobs use the existing serial runners: scientific work is limited
to 59 seconds, 1 GiB and one thread, and Lean to one 2048 MiB module followed
by a separate audit. Normal admission requires 2% free disk and 15% available
RAM; an active job stops below 1% disk or 6% RAM. Plans, metadata validation, restored
bytes and typeset PDFs are kept separate from mathematical execution evidence.

## Preparation checks

The [4 November preparation receipt](../../release/dimension_bounds/ideas_companion_20261104.json)
records successful metadata validation of all 13 proposition entries and
523 public evidence files, plus all four snapshot plans. These checks ran
under the current bounded monitor and launched no mathematical worker or
Lean compiler. Original mathematical execution remains in the dated records
linked above.
