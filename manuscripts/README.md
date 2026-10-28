# The six manuscript families

Read the **[paper guide](https://sashakolpakov.github.io/CoxeterPolytopes/papers.html)**
for the mathematical methods, assumptions, figures, page counts and evidence.
The 28 October right-angled submission editions add introductions, precise
reproduction links and presentation improvements within the existing paper
families. The other illustrated copies preserve their reviewed bases. Earlier
editions remain available as provenance.

| Family | Current PDF | TeX |
|---|---|---|
| General finite volume, dimension ≤ 57 | [Illustrated PDF](build/finite_volume57_illustrated.pdf) | [Source](finite_volume57_illustrated.tex) |
| General compact, dimension ≤ 17 | [Illustrated PDF](build/compact_farthest_vertex_bounds_illustrated.pdf) | [Source](compact_farthest_vertex_bounds_illustrated.tex) |
| Compact arithmetic, dimension ≤ 15 | [Illustrated PDF](build/compact_arithmetic_dimension_bound_illustrated.pdf) | [Source](compact_arithmetic_dimension_bound_illustrated.tex) |
| Sharp right-angled dimensions 8/4 | [7-page PDF, 2 figures](build/arxiv_sharp_dimension_bounds_8_and_4_20261028_layout_final.pdf) | [Submission source](revisions/right_angled_sharp_dimension_bounds_arxiv_20261028.tex) |
| Ideal RA4 volume and structure | [19-page PDF, 5 figures](build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_historical_context.pdf) | [Submission source](revisions/ideal_ra4_volume_and_structure_arxiv_20261028.tex) |
| Compact RA4 volume constraints | [PDF](build/compact_ra4_volume_constraints.pdf) | [Source](compact_ra4_volume_constraints.tex) |

The long compact papers contain necessary case proofs. Their approximately
19/38-page main expositions are not standalone proofs; the arithmetic bound
of 15 also uses the general compact bound of 17. No minimum essential length has been
established. Compact RA4 supplies complete constraint theorems and partial
progress toward the compact volume problems. The external
[Guo–Zheng preprint](https://arxiv.org/abs/2610.33053) now establishes the
compact 120-cell minimum; the compact second-volume problem remains open here.

Both submission PDFs compiled from their isolated single-file sources with
zero warnings. Their figures and bibliographies are inline; both name Alexander
Kolpakov and include an AI-use acknowledgment. The
[edition and reproduction guide](../docs/arxiv_submission_editions_20261028.md)
records exact evidence and the two suite IDs. No new mathematical checks or
Lean modules were executed for these presentation revisions.

The frozen 21 October editions state finite volume 57, compact 23, arithmetic 20,
sharp right-angled 8/4 and the two RA4 volume results. Their source/PDF hashes
remain in the [delivery record](../results/six_paper_final_tex_delivery_20261021.json).
The [result-to-proof map](../docs/result_manuscript_correspondence.md) and
[catalogue](catalogue.json) distinguish these editions and their selected
checks. A catalogue ID does not automatically replay a later illustrated copy.
Its `--pdf` option still builds the preserved 5/17-page right-angled publication
editions, rather than the new 7/19-page submission copies linked above.

For serial plans and resource rules, use [repro/README.md](../repro/README.md).
Written proof, source audit, exact computation, conditional Lean and PDF
compilation retain separate verification records.
