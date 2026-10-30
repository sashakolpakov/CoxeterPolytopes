# Manuscript organization and proof methods

Updated 30 October 2026. There are three arXiv submission manuscripts:
the combined dimension-bounds article, the sharp right-angled dimension article,
and the second ideal four-volume article. They consolidate the research
families without turning each numerical improvement into a paper.
The [submission guide](arxiv_submission_editions_20261028.md) lists their
standalone sources and the three current PDFs.

## Common methods and case divisions

| Submission | Common argument | Where the arguments separate |
|---|---|---|
| [Dimension bounds 57 / 17 / 15](../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex) | A common introduction presents the three theorems, their external inputs and proof organization. | The finite-volume proof uses angle payments and incidence curvature. The compact proof uses one original farthest vector and component energies. The arithmetic proof adds field restrictions and exact algebraic exclusions to that geometry, using the compact theorem internally. |
| [Sharp right-angled dimensions 8 / 4](../arxiv/sharp_dimension_bounds_8_and_4/main.tex) | Projection from a finite interior point to a facet of maximum height. | The general case uses a cubic spherical-code bound; the ideal case uses a cardinality bound and a metric-equality exclusion. |
| [Second ideal four-volume](../arxiv/ideal_four_dimensional_volume_and_structure/main.tex) | Defect counting and cusp geometry. | The volume gap and equality facet spectrum come first. A separately proved cutting and recognition theorem applies to the restricted octahedral facet spectrum, then identifies the equality double. |

The dimension article places all three main arguments before their complete
case appendices. Its 71-page opening exposition does not replace those proofs.
The compact and arithmetic arguments use the same original farthest vector
throughout each component-energy estimate. Every auxiliary point must satisfy
every original wall inequality. Arithmetic restrictions enter where integral
cycles and positivity at nonphysical embeddings constrain the roots; they are
not hypotheses of the general bounds.

The ideal-volume paper proves the facet hypothesis before invoking its chain
theorem. It does not claim that arbitrary ideal four-polytopes decompose into
24-cells. The compact RA4 constraint source remains supporting work rather
than a fourth arXiv submission.

## How computations are presented

A finite computation enters when the geometric proof has reduced to a precisely
stated finite obligation. The manuscript specifies:

1. The objects, equivalence relation, input domain and required output.
2. The algorithm and the invariant that justifies every pruning rule.
3. A termination bound and a proof that every necessary case is covered.
4. The exact conclusion obtained from the run and the proposition that uses it.

The [computational proof contract](computational_proof_contract.md) makes those
requirements explicit. A necessary forest list is not a list of realized
polytopes. Finite-output assumptions, inherited classifications and geometric
interpretation remain visible where the final theorem uses them.

The composition of Kedlaya, Kolpakov, Poonen and Rubinstein’s
[*Space vectors forming rational angles*](https://math.mit.edu/~poonen/papers/space_vectors.pdf)
provides a useful editorial model: an early proof strategy, an algorithm
followed by its correctness and termination theorem, and a separate account
of exact acceptance and completeness. The mathematical reductions here differ.

## Reading papers and auditing records

The papers contain definitions, complete arguments, mathematical algorithms
and explanatory figures. Source hashes, execution logs, raw enumeration rows
and earlier failed attempts belong in the repository’s verification records.
Selected Lean modules formalize conditional deductions; their geometric,
literature and finite-data premises remain explicit.

The [result-to-manuscript map](result_manuscript_correspondence.md) identifies
the proof homes, detailed maps and verification boundaries. The
[reproduction guide](../repro/README.md) distinguishes mathematical checks
from typesetting a current submission.
