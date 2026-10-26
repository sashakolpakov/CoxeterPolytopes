# Right-angled results and diagram-scope audit

26 October 2026. Author/reviewer `/root/ra_release_completion_audit`.
Scope: the mathematical claims used in the two new Sphinx chapters
`docs/site/right_angled.rst` and `docs/site/volume_problems.rst`, and the
caption contracts and source semantics of their two root-authored SVGs. I authored these
chapters, so this note is not their independent editorial acceptance.
It is a source-bound audit of the existing mathematical arguments and
interfaces, written by other authors. A separate author reviews the new
chapters; root owns acceptance and every actual build/deployment.

## Sources and accepted results

| Unchanged accepted source | SHA256 | PDF scope retained |
|---|---|---|
| `manuscripts/right_angled_sharp_dimension_bounds_publication.tex` | `0c29a7869fa46e4e1a542d8ae6f4dfc7e616136c7a1aad3984d941a95cf89722` | Five pages: finite-volume RA dimension at most8; entirely ideal at most4; both internally sharp at named published inputs. |
| `manuscripts/ideal_ra4_volume_and_structure_publication.tex` | `e0c2ca24cec24ffcd736e201d3c7e06d75b28ef7be8093decd321ae3ac54bda1` | Seventeen pages: published minimum; internal next volume/unique double; restricted Oct/D chains. |
| `manuscripts/compact_ra4_volume_constraints.tex` | `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f` | Fourteen pages: capacity/excess/recognition constraints, not unrestricted minimality or a second-volume theorem. |

I inspected the complete sharp and ideal proof chains during the preceding
source-only review, rechecked their unchanged identities and relevant
interfaces here, and read the compact capacity, metric, graph and final
W>=17 argument. The pages give full-article page counts, not invented
minimum essential-page counts. The finite/ideal/compact defect definitions
are kept separate.

## Mathematical claim checks, beyond acceptance labels

**Sharp dimensions.** Lemma `lem:projection` checks the future sheet,
unit norm and every other actual wall, including non-neighbors. Maximum
unit-normal height is used for the next height bound, not as a substitute
for feasibility. The strict finite descent, actual intrinsic ridge walls
and face inheritance lead to one distinct code direction per actual facet.
The tensor-square identity first gives N<=27; the four-vector negative
square excludes equality and is essential to obtaining26. The ideal branch
uses Musin for cardinality and Kolpakov for metric facet equality, then the
opposite-normal height contradiction. It assumes no kissing-code
uniqueness or arithmeticity. No missing all-wall or strict-endpoint step
was found in the stated written argument.

**Ideal volume and structure.** The defect identity and actual
opposite-cusp injection give W>=12 for a noncell. At equality, the
facet-spectrum reduction is proved before the chain theorem. The
seven-positive-facet branch uses actual embedded medial recovery and
cuboctahedral cusp widths; finite model data alone do not classify an
arbitrary remaining facet. Under the explicit metric Oct/D hypothesis,
the compactified internal triangles are proved to be closed surfaces.
The Euler bound does not assume spherical topology. Equal unit-normal
widths construct a Lorentzian hyperplane; convex-hull and actual-wall
containments prove it is an actual octahedral section. Outside facets
cannot touch its finite face/diagonal interiors. Cutting creates no finite
vertices and reduces the total D-count by eight. Whole-facet propagation,
the D4 common-neighbor obstruction, disjoint interiors and facet-symmetry
extension supply recognition, opposition and no-twist placement. No missing
recognition or outside-wall check was identified at the stated inputs.

**Compact capacity and excess.** The extra-facet count retains the
singleton rule: e_i=1 consumes a distinct extra that cannot serve another
neighbor, while remaining extras have capacity two. This gives f>=2V+5,
not merely a weaker incidence estimate. The neighboring-excess sum is
compatible with the exact factor-two ridge count. The W>=17 argument
distinguishes weighted nonpentagonal graph G from full positive-adjacency
graph H. Its five-neighbor exclusion uses marked hexagon/pentagon metric
data; the weaker abstract graph conditions alone admit C5[K2]. The
ten-vertex complement lemma closes all 0<W<=16. W=0 recognition requires
every facet to be a metric dodecahedron. No unsupported graph census,
single-facet propagation or minimum-volume conclusion was found in these
stated theorems.

## Evidence actually available and remaining interfaces

The current RA assembly/root acceptance are
`results/right_angled_release_packages_current_20261025.json` and
`results/right_angled_current_release_root_acceptance_20261025.json`.
The independent source/evidence review is
`results/review_right_angled_current_release_documentation_20261025.md`.
Their current-paper build receipts are
`results/sharp_ra_audit_edition_readability_build_20261021.json` and
`results/ideal_ra4_volume_and_structure_publication_build_20261018.json`.
The historical combined cores bind predecessor TeX; the chapters do not
promote them to combined current-paper executions.

The sharp exact workers cover the displayed 24-cell and universal rational
tensor identities. Its actual selected Lean audits cover12+10+8
declarations, plus five in the separate regular-cell supplement. The ideal
core has four exact model/combinatorial jobs and eight actually audited
ManuscriptIdeal4 declarations. `GapPartition`, `WrittenEqualityStages` and
`OctDGeometry` expose actual partition, equality, cutting, recognition,
path and gluing inputs. The stronger coordinate sharp kernels still need
actual wall/code realization and the identified contraction interfaces.
Printed foundational axioms do not enumerate ordinary theorem parameters.

Legacy IdealRA4Gap/CertificateBounds objects receive no new source-to-old-
object or separate-audit credit. A complete clean-checkout replay is not
claimed. Sources/configuration and local archived object bytes have
different retrieval/provenance boundaries. Successful unchanged kernels
were not rerun for this task.

The compact14-page paper has the actual zero-warning receipt
`results/compact_ra4_volume_constraints_W17_build_20261019.json`, four-figure
visual receipt `results/compact_ra4_w17_figure_visual_review_20261019.json`,
and complete written reviews
`results/review_compact_ra4_volume_constraints_W17_root_20261019.md` and
`results/review_compact_ra4_volume_constraints_W17_finite57_20261019.md`.
Its catalogue does not certify the universal geometry by a numerical
polytope census or geometric Lean proof. The new chapter says so.

Named literature inputs remain external. The exact RA6 statement was
authenticated in Alexandrov, with Dufour attribution; Dufour's original
full text was not obtained. Figure7's incidence transcription was checked
against the author's original image, separately from the external
small-facet completeness passage. The image is not redistributed. The
compact metric dodecahedron and Andreev-uniqueness inputs are stated
separately from the internal capacity and graph arguments. The site uses
the primary anchors already identified in these sources and provenance
reviews; it adds no new literature or worldwide-priority claim.

## Diagram contracts and findings

The accepted TeX sources contain1,5 and4 native figures respectively.
No accepted TeX or PDF was changed. I read the root-authored SVG sources
and aligned the two chapters' alt text and captions with their actual
flowcharts, rather than describing every arrow as geometric:

- `_static/right-angled-projection.svg`, SHA256
  `d4847e30791dd665234a4d0ad8daf7dec0accac5ac63348b11b370ba51283877`,
  now says "Finite interior point; maximum-height facet". The shared stage
  is the normalized all-wall-checked projection. Later arrows are logical
  dependencies in the separate general/ideal branches, not polytope edges.
  The text gives the strict9->8->7->6 height sequence; the SVG is a branch
  flowchart rather than a metric picture of that sequence.
- `_static/ideal24-chain.svg`, SHA256
  `edc9bea743cb32831d04b91fee75849f4c32095e8d348068f2e62e3d687a7dd0`,
  explicitly labels the entry from the arbitrary-polytope equality
  reduction to the Oct/D class "at equality: 8π²/3". Upper boxes describe
  restricted hypotheses and proof steps; only the bottom cell blocks and
  horizontal arrows show the abstract cell-adjacency path. Intermediate
  cuts are opposite and gluing has no independent twist. No existence
  for every chain length or decomposition/volume quantization for
  arbitrary ideal polytopes is implied.

The current reviewed chapter hashes are
`3322501fbcf56821bc279bb48f4d914a7050e199c5747a774efa69be0015e1c4`
for `right_angled.rst` and
`251369fe86d33b2cd20b18d182cb79c11090fd6bb807e0c5cf95b62d60b108b1`
for `volume_problems.rst`. SVG/caption source agreement does not certify
rendered layout or browser math rendering.

The six new TikZ figures in the three other working dimension manuscripts
have a separate independent source/byte-restoration review in
`results/review_illustrated_manuscripts_diagram_semantics_20261026.md`.
Actual builds and rendered inspection remain root-owned evidence.

**No material mathematical gap was identified within the advertised
theorems at their declared inputs.** Important limits remain explicit:
unrestricted compact120-cell minimality, a compact second-volume theorem,
general ideal24-cell decomposition, external classifications and actual
geometric Lean instantiation are not established by this documentation.
These open or external statements are not silently treated as proved.

One operational documentation finding is retained: the dated compact14
reproduction guide still prints the historical40/6GiB admission and8/3GiB
abort prose. The new chapter points readers to the current percentage
policy instead of copying it. The old guide was not edited here.
The complete independent chapter review is
`results/review_sphinx_right_angled_volume_chapters_20261026.md`.
Its three findings have been corrected: a missing backslash in the height
recurrence's spacing command, wording that conflated eight printed
declarations with eight separate audit processes, and an incomplete
explanatory sentence. The single actual ManuscriptIdeal4 audit ledger and
its eight printed declarations are now distinguished. No theorem,
formula, mathematical scope or existing execution record was changed.

This note launches no mathematical checks. No worker was imported; no
enumeration, test, Lean/TeX compilation, PDF/image rendering, publication
push or deployment was performed by this reviewer. Source inspection,
documentation authoring and metadata reads supply no new runtime PASS.
