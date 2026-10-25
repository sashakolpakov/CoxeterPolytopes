# Reproduce the six papers

The [six-paper reading table](../manuscripts/README.md) links the completed
TeX editions and matching PDFs. The frozen compact PDFs prove bounds of 23 and 20;
reviewed revisions in those families prove 17 and 15, respectively, and now
have matching verification PDFs. Their length and layout still need review.
These are distinct proof records. Neither written review nor selected
computation is an end-to-end geometric Lean proof.

## Inspect the selected checks

Run from the repository root. Start with plans, which launch no proof child:

```sh
python3 -B scripts/verify_manuscript.py --list
python3 -B scripts/verify_manuscript.py finite57 --mode full --pdf --plan
python3 -B scripts/verify_manuscript.py right_angled_sharp --mode core --pdf --plan
```

The [catalogue](../manuscripts/catalogue.json) defines the selected commands;
each paper's map states their mathematical scope and the inherited inputs.
Use the corresponding ID with `--mode core --pdf --plan` to inspect another
paper. A plan or a catalogue entry is not execution evidence.

| Paper | Verification ID | Proof and reproduction map |
| --- | --- | --- |
| General finite volume 57 | `finite57` | [Forests, cusps and payments](../docs/finite_volume57_lemma_map.md); [bounded reproduction](../docs/finite_volume57_reproduce.md) |
| Compact arithmetic 20 | `compact20` | [Fields and exact angle rules](../docs/compact_arithmetic20_proof_map.md) |
| General compact 23 | `compact_farthest_bounds` | [Farthest vertex and component energies](../docs/compact_farthest_vertex_bounds_proof_map_20261020.md) |
| Sharp right-angled 8/ideal 4 | `right_angled_sharp` | [Projection, heights and codes](../docs/right_angled_sharp_dimensions_publication_proof_map_20261021.md) |
| Ideal RA4 volume and structure | `ideal_ra4_volume_and_structure` | [Volume, equality and restricted chains](../docs/ideal_ra4_volume_and_structure_proof_map.md) |
| Compact RA4 volume constraints | `compact_ra4` | [Capacity and excess](../docs/compact_ra4_volume_constraints_reproduce_20261019.md) |

**Keep finite57's full command in plan mode.** Its historical146-job selection
contains136 scripts, nine conditional Lean roots and one PDF build. All 25
distinct bounded replacement units have actual PASS receipts and are
[integrated into the primary catalogue](../results/finite57_complete_loader_catalogue_integration_20261022.json).
The111 other command matches retain their inherited execution scope;
the core11 selection is unchanged. There is no combined146-job runtime
PASS or complete transitive regeneration. The
[current-rule evidence guide](../docs/finite_volume57_reproduce.md)
records the successful 64-shard/259-unit application conjunction and its
explicit geometric and finite-data premises. This is a separate verification
scope, with inherited inputs rather than complete transitive regeneration.
Earlier failures and resource deferrals remain preserved.

The [verification guide](../docs/manuscript_verification.md) explains the
`scripts`, `logic`, `core` and `full` selections. The smaller finite57 core
does not replay the 136 scripts. SharpRA has recorded exact/model and selected
conditional Lean checks; historical combined runs keep their original source
bindings. Compact RA4's current selection checks typesetting only.

## Inputs and bounded execution

The [public-data guide](PUBLIC_DATA.md) explains manifest-bound release assets.
Inspect their plan without network requests or writes:

```sh
python3 -B scripts/fetch_public_data.py --plan
```

Once assets are available, `--restore` restores all missing manifest-listed
inputs sequentially. Repeated `--path` selects exact paths;
`--from-directory release-assets` uses downloaded archives. Byte restoration
does not verify mathematics, and these instructions do not claim an upload.
Recorded proof-package versions are in
[requirements-proof-20261018.txt](requirements-proof-20261018.txt); Lean's
toolchain is pinned in [formal/lean-toolchain](../formal/lean-toolchain).
Do not launch a clean whole-project build as a reproduction shortcut.

All scientific, Lean and TeX checks share one exclusive lock and run serially.
The user's25October override measures percentages of total filesystem
capacity and physical RAM: admit at **2% free disk and 15% available RAM**,
stop active jobs below **1% disk or6% RAM**, and start cleanup below 2% disk.
Gate on exact byte ratios and record totals, percentages and the policy source
hash. The coordinator uses the existing bounded research runner
`scripts/run_bounded_research_check_20261018.py` for scientific/TeX children
(59 seconds, 1 GiB, one thread) and
`scripts/run_lean_batches_signal_safe_20261018.py` for one Lean module
(2,048 MiB, one thread), with a separate axiom audit. Keep monitoring enabled;
resource refusal is not permission to lower a floor. Use new ledger paths,
preserve failed/deferred receipts, and remove only owned unused intermediates
after a batch. Retain sources, certificates, PDFs and evidence.

Lean verifies selected internal implications with explicit geometric,
classification and script premises. The
[assumption register](../docs/manuscript_assumptions.md) identifies the boundary;
printed axioms do not enumerate theorem parameters. Regular24 now has
**bounded compile PASS and five actual axiom outputs**; ClassicalDimensionComposition
has **bounded compile PASS and four actual outputs**, each from a separate
successful audit. The
[completion receipt](../results/lean_publication_completion_20261021/manifest.json)
binds unchanged sources and the four serial runs, with no successful dependency
reruns. Only standard axioms occur and no `sorryAx` appears. The coordinate,
geometry, field and script boundaries remain explicit in their maps. These
two selected compilation/audit prerequisites are complete; final publication
checks and remote upload remain separate actions.

## Build documentation only

The [Sphinx site](../docs/site/index.rst) uses pinned Sphinx 9.1.0 and its
built-in Alabaster theme. A documentation build runs no scientific, Lean or
TeX checks:

```sh
python3 -m pip install -r repro/requirements-docs.txt
python3 -m sphinx -W --keep-going -b html docs/site docs/site/_build/html
```

The former thirteen-document aggregate, atlas regeneration and raw historical
commands remain byte-preserved in the
[archived README](../results/source_archive/repro_readme_before_publication_20261021/README.md.gz)
and its [preservation manifest](../results/source_archive/repro_readme_before_publication_20261021/manifest.json).
They are historical instructions, not the current six-paper recipe.
