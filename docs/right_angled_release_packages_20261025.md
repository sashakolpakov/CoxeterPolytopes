# The two right-angled release packages

25 October2026. These are two existing paper families. The assembled
[current manifest](../results/right_angled_release_packages_current_20261025.json)
binds the unchanged paper editions and current documentation
to their builds, written audits, unchanged computational sources and
prior actual formal-verification records. Assembly is metadata work;
no successful mathematical check or Lean module was rerun.
The [predecessor assembly review](../results/review_right_angled_release_assembly_source_only_20261025.md)
accepts the reused evidence at its stated boundaries. The current assembly
refreshes documentation and retrieval pins while preserving every execution
record. The explicit legacy CertificateBounds source entry supplies no new
source-to-object verification credit.

| Paper | Current edition | Written result |
|---|---|---|
| [Sharp right-angled dimensions](../manuscripts/right_angled_sharp_dimension_bounds_publication.tex) | 5-page PDF | Finite-volume dimension at most8, entirely ideal dimension at most4; both internally sharp. |
| [Ideal RA4 volume and structure](../manuscripts/ideal_ra4_volume_and_structure_publication.tex) | 17-page PDF | Published24-cell minimum; internal next-volume and unique-double theorem; restricted octahedral-chain decomposition. |

The current source/PDF hashes match their separate successful zero-warning
build records. Neither historical combined core ledger has the current
TeX hash. Reusing its unchanged computational and formal components does
not turn it into a fresh combined run of the current paper.

## Written proof and exact computations

The [sharp-dimension proof map](right_angled_sharp_dimensions_publication_proof_map_20261021.md)
separates actual facet projection and face inheritance from the cubic
spherical-code argument and the entirely ideal branch. Its two exact
workers check the specified24-cell coordinates and universal polynomial
identities. The latter expands36 rank-two and216 rank-three tensor terms
over rational polynomials in twelve independent variables. This is an
algebraic identity check, not a search for spherical codes. The coordinate
worker does not claim complete vertex enumeration of arbitrary polytopes.

The [ideal-volume proof map](ideal_ra4_volume_and_structure_proof_map.md)
keeps the defect gap, equality facet spectrum, actual cuts, recognition
and gluing as distinct written steps. Its fixed-model algorithms are
specified and justified in the publication source:

- Enumerate10,626 four-wall subsets of the24-cell, solve nonsingular
  rational systems, retain feasible solutions and remove duplicates.
  Every vertex has four independent active walls, which proves completeness.
- Recover faces by closure under intersection with facet vertex sets.
  There are finitely many subsets, and every face is an intersection
  of its containing facets. The double's complete vertex set is first
  justified by its convex reflection-union description.
- Enumerate364 wall triples for the rational cuboctahedral model;
  enumerate the fifteen six-vertex perfect matchings for the six-cusp
  graph; solve the six cube equations in eight unknowns exactly.
- Check all276 unordered D4 root pairs and the specified reflected
  octahedral incidences. This finite support does not prove the universal
  surface, cutting or metric recognition arguments.

The reused core has four actual exact-check records: saved double faces,
equality-model replay, an independent equality-data checker,
and D4 roots. The full catalogue additionally offers the historical
double-model producer replay. No complete four-polytope census is assumed.

I compared the complete16-to17-page TeX difference during assembly. It
adds the fixed-model algorithm descriptions and their termination and
completeness arguments, plus one spacing correction. The theorem
statements, universal gap/equality/cutting proof and model-check inputs
are unchanged. The [fresh complete source audit](../results/review_ideal_ra4_fresh_audit_root_acceptance_20261024.md)
accepts the current edition at the declared external/finite interfaces.

## Formal verification and assumed inputs

| Selected kernel | Actual audited declarations | Remaining inputs |
|---|---:|---|
| RightAngledSharpDimensions |12|Actual face/point interpretation, higher-face inheritance and named published facts.|
| RightAngledCubicScalarKernel |10|Actual tensor/code realization and the stated contraction interfaces.|
| RightAngledProjectionHeightKernel |8|Actual walls and ridges, positive normalizers, future-sheet interpretation and real geometric realization.|
| RegularIdeal24CellHeightKernel |5|The regular24-cell geometric instance and named scalar/coordinate hypotheses.|
| ManuscriptIdeal4 |8|Published minimum and classifications, defect/equality stage hypotheses, actual cuts, recognition, chain and gluing inputs.|

Each listed module has successful recorded compilation and a separate
actual printed-axiom audit. The manifest identifies their exact sources
and ledgers. It does not infer assumptions from the foundational axiom
list: ordinary theorem parameters are exposed in the
[assumption register](manuscript_assumptions.md), the
[projection logic map](right_angled_projection_height_logic_20261020.md)
and [formal guide](../formal/README.md).

The sharp paper's catalogue core selects the first three sharp kernels,
with30 audited declarations. RegularIdeal24CellHeightKernel is a separately
recorded supplement with five declarations; the core command below does
not select it. The release manifest includes that supplement's existing
successful compilation and separate actual axiom audit explicitly.

The ideal paper's historical core also reused compiled IdealRA4Gap and
CertificateBounds dependencies. That core explicitly reports legacy
uncapped dependency caches; it supplies no fresh capped compilation or
new separate axiom audit for those dependencies. Their status is not
upgraded by this assembly. Their current sources are pinned for retrieval;
no current-source-to-old-dependency-object link is inferred from those pins.
No end-to-end geometric Lean proof is claimed.

The precise RA6 minimum of27 facets is directly authenticated in
Alexandrov's accessible statement, attributed there to Dufour; Dufour's
original full text was not obtained. The ideal paper's Figure7 incidence
transcription was independently compared with the author's original
image, and the published classification passage was read separately.
These provenance findings are in the
[primary-input review](../results/review_ra_release_primary_inputs_root_20261024.md).
The third-party image is not part of the release payload.

The [predecessor manifest](../results/right_angled_release_packages_20261025.json)
retains its original documentation snapshots and execution scopes. The current
assembly refreshes those document pins and directly records the already used
model producer, descending-height proof and pinned Lean configuration. It
introduces no new mathematical or formal-verification execution.

## Cleanup and subsequent replay

The [completed cleanup](../results/tex_lean_repro_only_cleanup_20261025.json)
removed loose compiler products after verification. Sources, final PDFs,
certificates, compilation/audit records and scientific inputs remain.
Source-bound selected objects may be compiled again when a new replay is
needed; compilation must precede its separate actual axiom audit.
The local preservation archive retains the exact legacy objects and their
trace/hash provenance. It is not a new source-to-legacy-object verification,
and a source checkout alone does not supply that local archive or the
unuploaded scientific assets. Restore at the recorded provenance boundary
or use the documented bounded source compilation as appropriate. The runner
rejects objects older than their sources/imports; never change timestamps to
force stale-cache acceptance.

## Running a necessary replay

Use the guarded public entry from the repository root. Planning runs no
mathematical child:

```sh
python3 -B scripts/run_proof_batch_20261024.py manuscript right_angled_sharp --mode core --pdf --plan
python3 -B scripts/run_proof_batch_20261024.py manuscript ideal_ra4_volume_and_structure --mode core --pdf --plan
```

To execute a necessary replay, omit `--plan` and supply a fresh `--ledger`
path in `results/`. Execution remains serial under the shared lock and
requires at least 2% free filesystem capacity and 15% available physical RAM.
Active jobs stop below 1% disk or 6% RAM. Gates use exact byte ratios and
record totals, percentages and the policy source hash. Compile and actual-axiom
audit are separate stages; the temporary PDF build preserves the existing
edition. Successful unchanged stages need no rerun merely to refresh a date.

If a substantive change requires replaying the separate regular24-cell
supplement, use its two named stages with new ledger paths:

```sh
python3 -B scripts/run_proof_batch_20261024.py lean-compile CoxeterProofs.RegularIdeal24CellHeightKernel --ledger results/regular24_needed_compile_new.json
python3 -B scripts/run_proof_batch_20261024.py lean-audit CoxeterProofs.RegularIdeal24CellHeightKernel --ledger results/regular24_needed_axioms_new.json
```

No such replay is necessary for this unchanged supplement in the present
assembly. The second stage audits the matching compiled target and does
not silently compile a missing one.

This package assembly is separate from the final strict documentation
build, actual Git-index audit and publication push. No bulk scientific
archive has been uploaded.
