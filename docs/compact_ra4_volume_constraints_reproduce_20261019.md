# Compact right-angled four-polytopes: proof and reproduction map

For current reading, see the [complete compact RA4 constraint source](../manuscripts/compact_ra4_volume_constraints_publication.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

19 October 2026. One current paper:
[compact_ra4_volume_constraints.tex](../manuscripts/compact_ra4_volume_constraints.tex),
SHA256 `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`.
The consolidated W>=17 source has
[root](../results/review_compact_ra4_volume_constraints_W17_root_20261019.md)
and [independent full manuscript reviews](../results/review_compact_ra4_volume_constraints_W17_finite57_20261019.md).
All four vector figures are source-native. The
current14-page PDF
has a successful [source-bound build ledger](../results/compact_ra4_volume_constraints_W17_build_20261019.json):
zero warnings,1.310 seconds,under86MB process-tree RSS,with temporary
files removed. Its SHA256 is
`920edadbbd7a973d70931e6ced9b8d0d5cfa12922d3b5dc0def3dc19ed2e8bf2`.
All four figures have a [visual-review ledger](../results/compact_ra4_w17_figure_visual_review_20261019.json);
disposable raster previews were removed after inspection.

The paper supersedes the capacity note in the same manuscript family.
The seven-page predecessor
and its [source/PDF freeze](../results/compact_ra4_facet_capacity_20261018_freeze.json)
are preserved. That PDF proves the older capacity theorem, not the
successor results below.

## Results and dependencies

For a compact right-angled four-polytope, let f and v count facets and
vertices. Put `W=sum_R(k_R-5)=v-5f`, summed over polygonal ridges. A
facet with V vertices has excess `d=(V-20)/2`, and `sum_F d(F)=2W`.

|Result|Proof location|Verification boundary|
|---|---|---|
|A three-dimensional k-gon requires at least 2k+2 facets; equality has the Loebell face lattice|Section 3|Complete dual-disc and closing-ear argument in the paper; no numerical census.|
|Four-dimensional capacity f>=2V+5|Section 4|Complete perpendicular-edge and incidence proof inherited from the frozen predecessor.|
|W>=ceil(3(V-20)/4)|Section 5|Sum excesses over distinct neighboring facets using the three-dimensional capacity and exact excess identity.|
|Facet sizes at most 46 below the 120-cell volume, 46 at that volume, and 96 below the double volume|Section 5|The volume identity and preceding inequalities; these are necessary local cutoffs, not volume minima.|
|W=0 forces the regular right-angled 120-cell|Section 6|Explicit three-dimensional dodecahedron metric input followed by full-facet propagation.|
|A degree-two positive facet forces W>=20|Section7|Full L6 metric and actual common-base proof; two disjoint rows each supply five extra facet-excess units.|
|W=0 or W>=17|Sections7–9|One complete polygon and graph proof excludes every0<W<=16 together,using the degree-two lemma,the universal five-neighbor-link exclusion and the ten-vertex complement lemma; no graph census is assumed.|

The volume identity is

```
Vol(P) = pi^2*(f+W+16)/12.
```

It follows from the stated reflection-orbifold Euler and
Chern–Gauss–Bonnet inputs and elementary face counts. The paper does not
prove `f+W>=120`, identify the second compact
volume, or improve the compact right-angled dimension ceiling of four.
The preceding W>=15 and W>=16 source/PDF versions are preserved
supporting evidence,not additional current papers.

## Explicit external and internal inputs

External inputs are the named acute-chamber background, elementary
Lorentz and convex-polytope facts, orbifold Gauss–Bonnet, the compact
three-dimensional twelve-facet metric dodecahedron classification, and
Andreev's marked metric uniqueness. The metric references and their
precise uses are checked in the full manuscript review. No unpublished
four-dimensional minimum theorem is used.

The perpendicular-face construction, dual-disc capacity, neighboring
excess estimate, four-edge closure, eligible-neighbor condition and
small-complement graph proof are internal written arguments. They have
written mathematical reviews. They are not geometric Lean theorems or
script-enumeration assumptions.

## Push-button build

From the repository root, run one serial batch, using a NEW ledger path:

```sh
python3 -B scripts/build_illustrated_ra_manuscript_with_ledger_20261019.py manuscripts/compact_ra4_volume_constraints.tex --check-only --ledger results/compact_ra4_volume_constraints_rebuild_20261019.json
```

This entrypoint enforces the shared lock, 40 GiB disk and 6 GiB available
memory launch gates, 8/3 GiB abort gates, 59-second limit and 1 GiB
process-tree limit. It removes its own TeX intermediates. Once a
PDF exists,use `--check-only` as above and a new unused ledger
argument to check a rebuild without replacing it; check-only reports JSON to
stdout. No unrestricted Lean aggregate build belongs to this paper's
reproduction path.

The machine-readable [catalogue](../manuscripts/catalogue.json) maps the
same paper to `compact_ra4`. It lists no numerical or Lean check for
the universal geometry. A successful TeX build certifies rendering and
source identity; the written proof review remains a separate layer.

## Preserved supporting evidence

* [Preparation and full proof map](compact_ra4_volume_constraints_preparation_20261019.md).
* [Neighboring excess and zero-excess recognition](../results/compact_ra4_neighbor_excess_and_zero_excess_rigidity_20261019.md),
  with [root review](../results/review_compact_ra4_neighbor_excess_root_20261019.md).
* [Positive-facet gap proof](../results/compact_ra4_positive_facet_link_gap_fifteen_20261019.md),
  with [root](../results/review_compact_ra4_positive_facet_gap_fifteen_root_20261019.md)
  and [independent](../results/review_compact_ra4_positive_facet_link_gap_fifteen_finite57_20261019.md)
  complete mathematical reviews.
* [Fifteen-excess exclusion](../results/compact_ra4_fifteen_excess_equality_exclusion_research_20261019.md),
  with [root](../results/review_compact_ra4_fifteen_excess_exclusion_root_20261019.md)
  and [independent](../results/review_compact_ra4_fifteen_excess_equality_exclusion_finite57_20261019.md)
  complete reviews.
* [Sixteen-excess exclusion](../results/compact_ra4_sixteen_excess_equality_exclusion_research_20261019.md),
  with [root review](../results/review_compact_ra4_sixteen_excess_exclusion_root_20261019.md),
  and [degree-two excess20](../results/compact_ra4_degree_two_excess_twenty_20261019.md),
  with [root review](../results/review_compact_ra4_degree_two_twenty_root_20261019.md).
  The current source combines these without repeating staged global proofs.
* [Primary-source comparison](../results/compact_ra4_primary_literature_comparison_20261019.md):
  the perpendicular-opposite mechanism and excess identity have explicit
  published antecedents; the comparison does not claim worldwide priority.

The current supporting source contains the resulting constraint proofs.
Earlier draft production is recorded only by the dated execution receipts;
those drafts are not separate reading or reproduction targets.
