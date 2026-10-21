# Independent review of the two E8 cusp completions

**PASS.** Source and geometric proof review found no issue in the [79-wall construction](fresh_h20_e8_cusp_completion_20261014.md). A fresh primary saved replay and a [separate Fraction checker](../scripts/fresh_h20_e8_cusp_completion_20261014_review.py) both pass; its [certificate](fresh_h20_e8_cusp_completion_20261014_review.json) binds the frozen source and input data.

The independent checker imports no producer arithmetic. It uses squarefree-radical multiplication in the basis 1,√2,√79,√158 and outward rational square-root bounds at 160 bits, whereas the primary uses nested exact squaring signs.

The following were independently reconstructed or checked:

- The E8 inverse, its norm-two negative fundamental weight, and the exact formulas r_i=−ω+8ν_i.
- All 79 norms, the 155 new relations (39 right, two π/3, 114 dotted), the common strict center, and every escape inequality.
- Complete incidence at both selected null rays, each positive full-support affine kernel, a positive-definite principal minor in every component, and the exact positive root sum proportional to the corresponding null ray.
- All 70 new-cap pairings at the 35 preserved H12 finite vertices.

The two links have affine block sizes 4+9+9 and 12+9. Their mutually orthogonal positive-semidefinite quotient blocks have total rank 19. Each positive kernel forces any link recession vector to be orthogonal to its whole block; together the blocks span the nineteen-dimensional Euclidean quotient. The recession cone is therefore zero. This is a bounded-link proof using positive kernels and full rank, not a rank-only inference.

The actual-facet argument was source-reviewed: every new and retained facet receives the exact timelike projection of the common strict center, and every nonincident inequality is checked strictly. The independently certified 77-wall geometry and H12 source completeness remain explicit dependencies; this review does not repeat their complete pair or vertex censuses.

Both new roots project onto positive multiples of future null rays on the H12 face, so they do not change its intrinsic chamber. Their pairings with all 35 finite face vertices are strictly negative. The old exact spacelike escape also survives all 79 walls, with square

\[
\frac{269344645863312077}{40000000000000000}>0.
\]

Thus the two specified ambient cusp links are repaired, while the full chamber still has infinite volume. No complete ambient cusp or vertex census is claimed.

    python scripts/fresh_h20_e8_cusp_completion_20261014.py --verify
    python scripts/fresh_h20_e8_cusp_completion_20261014_review.py --verify
