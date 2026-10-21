# The exact six-simplex minimal-diagram frontier

This is the necessary-domain certificate for the
[six-simplex finite-edge proof](fresh_finite_edge_cusp_tradeoff_28over3_20261017.md).
It reuses the frozen26 high-rank quasi-Lannér types and does not repeat
their146,262-attachment census.

Seven black roots and zero through three internal white roots give539
minimal-diagram choices. Exactly175 have at least two affine black
deletions,290 have one, and74 have none. Affine-kernel orthogonality and
the complete scalar-threshold lists prove that the first465 cases have
at most three ideal vertices when at least two vertices are finite.
The ideal-count envelopes are140 copies of `{2}`,35 of `{3}`,247 of
`{1,2}`,36 of `{1,3}`, and seven of `{1,7}`. No easy case permits exactly
five ideal vertices.

The74 remaining black blocks comprise four affine rank-six blocks and70
positive-definite rank-seven blocks. Their complete direct-white domains
use a strengthened choice: choose the indefinite subdiagram L with
minimum order in the entire actual star. For an affine deletion F=L-h,
where h is white, any other white u touching F makes F+u indefinite.
Its proper subdiagrams have smaller order and are PSD, so F+u belongs to
the same-order quasi-Lannér catalogue. Enumerating every affine deletion
and every label-preserving graph isomorphism gives the complete list of
contact profiles on the literal F. Root norms are uniformly normalized
to two; edge labels and their normalized products are preserved.

The [checker](../scripts/fresh_six_simplex_ideal_frontier_20261017.py)
and [JSON](fresh_six_simplex_ideal_frontier_20261017.json) record all539
choices, every easy threshold and tie, all profile provenance maps,
all literal projection-cost matrices, and every mandatory-root pair test.
The finite list concerns whites directly touching F only. Whites right
to F may have arbitrary rank and graph; their analytic Schur correction
is handled in the linked proof. No finite-rank extrapolation or full
six-simplex conclusion is inferred from this domain computation alone.

Generation and fresh saved replay pass. The companion theorem is now a
complete internally checked lemma, with a separate
[independent review](review_finite_edge_28over3_20261017.md). That review
checks the all-size proof and independently recomputes the remaining
small exact systems; it does not repeat this539-case catalogue and
profile-isomorphism enumeration. The frozen source/data retain their
historical status wording. No predecessor was edited. **FINAL/FROZEN.**
