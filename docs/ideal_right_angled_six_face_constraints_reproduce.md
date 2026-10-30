# Ideal right-angled5/6: proof and reproduction map

For current reading, see the [current sharp right-angled manuscript](../arxiv/sharp_dimension_bounds_8_and_4/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

The current12-page manuscript
(supporting argument) is ONE
paper with four vector figures and shared geometric arguments. It proves
necessary conditions, not existence or a dimension exclusion.

|Entirely ideal class|Facets at least|Vertices at least|
|---|---:|---:|
|Dimension5|63|211|
|Dimension6|275|6106|

It also proves272V>=6039F in dimension6. In mixed ambient dimensions11/12,
the ideal five-face fraction is<4/47 and the ideal six-face fraction is
<15/248, using the explicitly stated published NK/Dufour inputs.
The [full theorem/dependency map](ideal_right_angled_six_face_constraints_preparation_20261019.md)
identifies the shared lemmas and all four figures.

## Current verification

The complete written arguments have root and independent mathematical
reviews. The new eleven-vertex equality exclusion is a geometric proof,
not a bounded enumeration of eleven-vertex polyhedra. Its internally
proved small-facet bridge cites the unified ideal4 paper and Figure7
extension; only the stated three-dimensional classification is published
input. Arbitrary four-polytopes are not assumed to decompose into24-cells.

The [bounded build ledger](../results/ideal56_figure_spacing_build_20261019.json)
records12 pages, two passes, zero warnings,1.073s and86,228,992 bytes
peak process-tree RSS. It started with187.21GiB disk and6.11GiB available
RAM. Temporary files were removed. Source SHA256:
`782cea54d88970b70e1d6692934f48f5727a46b262bf6b49e342bd9cd68493ad`;
PDF SHA256:
`a89a13677f3d7256afe5c2837552a669c5dc11e5967119324c6881b7a582d376`.
All four figures were visually inspected; a spacing correction in Figure4
was rebuilt and rechecked. Temporary raster previews were removed.
Typesetting success is separate from mathematical or kernel verification.

The historical carrier, transport and sharp-reciprocal exact checkers
remain unchanged. Their small frame and sign identities support the
universal written arguments; an88-frame list is not a proof for every
cusp width. Their old final counts are preserved as predecessor results,
not substituted for the current63/211 and275/6106.

The new [conditional Lean source](../formal/CoxeterProofs/RightAngledElevenVertexCounts.lean)
and [assumption handoff](right_angled_eleven_vertex_counts_lean_preparation_20261019.md)
derive the current integer counts and cleared fractions from explicit
geometric/count structures. The [current compilation and separate axiom audit](right_angled_eleven_vertex_counts_verification_20261019.md)
pass with no warnings or placeholders. The first resource deferral and
corrected elaboration error are documented separately from that success.
The preceding62/206 and269/5875 implications have their separate successful
[Lean execution histories](right_angled_cusp_capacity_verification_20261019.md).
The geometric statements remain written internal obligations, not newly
introduced published axioms.

The [current composed verification ledger](../results/ideal56_current_composed_verification_20261019.json)
binds the three successful exact replays, current conditional Lean/axiom
audit, separate final PDF build and visual review. The attempted single
five-job replay passed its first four jobs, then deferred the check-only
PDF job at the memory gate. Its original aggregate ledger still records
that deferral. The composed record instead checks that the separately
successful PDF build has the identical final source and canonical PDF;
it does not relabel a deferred execution as successful.

## Push-button replay

From the repository root, inspect the exact job list without computation:


Replay the three declared exact checkers, current conditional Lean logic
and a fresh PDF compilation serially, preserving the canonical PDF:


Choose a fresh ledger filename for another replay. The command refuses
to overwrite an existing ledger. Each scientific/TeX job takes the shared
lock and retains its resource gates and caps; a resource deferral is not
a successful verification. The resulting record states each job's scope
and the remaining geometric gaps. The script results may be taken as
explicit premises for logical verification; this does not formalize the
enumerator implementations or the underlying hyperbolic geometry.

For only the latest one-module logic and its printed-axiom audit:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.RightAngledElevenVertexCounts --audit-axioms
```

A fresh check-only typesetting command is:


Check-only mode prints the build record and publishes neither a new PDF
nor the unused ledger path. The old five-page source and PDF are preserved
in byte-verified gzip archives under`manuscripts/source_archives`; the
[historical cusp-package ledger](ideal6_cusp_package_audit_20261018.json)
retains its original provenance. No frozen proof packet was replaced.
