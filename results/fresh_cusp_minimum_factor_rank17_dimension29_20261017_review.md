# Independent review of the rank17 cusp-factor dimension29 composition

**FINAL/FROZEN: mathematical/source review, independent exact generation and fresh saved replay all pass.** This review concerns the restricted statement that every irreducible affine factor at every cusp has rank at least17. It makes no unrestricted dimension29 claim.

The [review checker](../scripts/fresh_cusp_minimum_factor_rank17_dimension29_20261017_review.py) uses only independent `Fraction` polynomial arithmetic, with no producer arithmetic imports or symbolic algebra library. Its [certificate](fresh_cusp_minimum_factor_rank17_dimension29_20261017_review.json) binds the final composition source, data and report, and both affine-carrier proofs. The geometric review below is a written proof/source audit; it does not rerun the inherited large enumerations or independently rederive Vinberg's published finite-diagram tables.

## Geometric hypotheses and exhaustive payments

An actual connected affine subdiagram has a positive null root combination. Its product with a strict chamber center is negative, so that combination is nonzero. Acuteness makes the null point feasible. Strict convexity of the Klein ball makes it an ideal vertex. In the Euclidean product-of-simplices link, an irreducible affine component has only spherical proper principal subdiagrams. Thus the connected affine subdiagram is an entire cusp factor. The rank17 hypothesis consequently forbids every connected affine diagram of order at most17; it is stronger than a bound on the number of cusp factors.

At a two-face ideal corner, the link edge varies in one simplex factor. Both marked side roots lie in that factor, and the other factors have one node omitted and are spherical spectators. Every marked pair in a factor of order r+1 has exactly V/(r+1) such fibers. Therefore the all-rank factor cost3(r+1) gives a total ideal-vertex cost3sV. This count does not assume that the ideal vertex is simple. The finite-vertex cost remains3n.

For a triangle, every indefinite subdiagram of its actual full star contains all three black roots: omitting one black leaves a positive-semidefinite corner. A minimal indefinite connected subdiagram is Lannér or quasi-Lannér. The published quasi-Lannér order bound10 supplies a proper affine subdiagram of order at most9, prohibited here. The remaining Lannér core has order at most5. The audited finite triangle cases provide the required payments. If a contributing corner is ideal, a connected path through at most four remaining core nodes has distance at most3 within its varying factor. The high-bond H3/H4 case cannot instead lie in a classical affine varying factor. Extra cusp factors do not change that role assignment. No compact-ambient theorem is being applied to a noncompact polytope.

For quadrilaterals with an ideal corner, the varying factor has order at least18, so the separately reviewed large-factor payment proof applies. Its spherical spectators give a positive-semidefinite Schur correction with the required sign. The exceptional E7/E8 roles have their own attachment checks; they are not inferred from monotonicity of arbitrary finite component weights.

For a quadrilateral with four finite corners, the actual face is compact. Its four corner diagrams are positive definite, the opposite-side diagrams are indefinite, and its full star has rank n+1 and corank one. These local facts remain valid when other vertices of the ambient polytope are nonsimple and ideal.

The two carrier exclusions apply to every connected affine subdiagram of H−U or H−V, rather than merely to a private diagram that happens to be wholly affine. This distinction is essential. Every such affine carrier contains the appropriate opposite black pair. The cycle and tree cases exhaust all possible large affine carriers, so both complements are affine-free. Their subdiagrams U,V are affine-free as well. The subsequent finite-quadrilateral transfer uses only its audited carriers of at most17 nodes, or positivity of diagrams orthogonal to an actual indefinite subdiagram, together with the retained published finite classifications and determinant arguments. These hypotheses are now all available. The nonzero determinant coefficient follows from the actual full-star dependence, not a hypothetical smaller face star.

Thus triangles receive at least2, quadrilaterals at least1, and larger polygons at least5−k by nonnegativity. The branch composition is exhaustive.

## Independent hostile read of the carrier arguments

In the [cycle proof](fresh_affine_cycle_carrier_exclusion_20261017.md), the two old arcs between any pair of contacts force that pair to be the two marked roots: otherwise one cycle survives a spherical deletion. This also excludes three contacts. A sole contact has label3/4 by the long finite deletion, and creates a forbidden small affine leaf witness. A shortest connecting white path consequently supplies a two-contact node and an exterior predecessor. The remaining terminal B-path cannot be extended at its doubled end in the actual finite corner. The proof never assumes the cycle exhausts private V.

In the [tree proof](fresh_affine_tree_carrier_exclusion_20261017.md), the pair-path argument again bounds support by two. The actual connecting path justifies the high-label exclusion by finite-corner classification and payment. Grounded tree resistance gives

    q(s)=t_x² nu_x² s+t_y² nu_y²(L−s)<2.

Both coefficients are at least1/4, so L<8. The complete endpoint-role positions give precisely the stated11 surviving rows; the range is derived, not a large-rank sample. The fork-switch continuation is a small affine doubled star. The two-step terminal continuation gives the full affine B/D diagram inside a finite corner, where its large order is irrelevant. Only the terminal triangle remains.

Once that terminal connector exists, the long finite B/D deletion has only its long-end extension available. A second neighbor at that end gives an affine positive kernel or a negative perturbation, excluding every additional surviving contact. A remaining sole contact at the removed terminal is ruled out by its small affine witness or the high-bond finite-corner payment. The opposite black roots are already orthogonal to the carrier, so that last additional node is white and the cited white path belongs to the same finite corner. The unique connector then makes the whole carrier exactly private V, allowing the frozen private-tree exclusion. No private-exhaustion assumption is used prematurely.

These checks include the source's role range and graph conventions, but not a duplicate908-case arithmetic run. The independent review found no remaining gap in either carrier argument or their use in this composition.

## Independent all-rank rational checks

The checker derives all four breakpoint formulas for the increment of the mean absolute value of a centered discrete uniform. With N negative terms,

    sum_{l=0}^r |y+l−r/2|=(r+1)y+N(r−2y−N+1).

Substituting the appropriate N for integer and half-integer y and for both parities of r gives the four polynomial identities after clearing denominators. The zero boundary cases can be checked directly; the displayed formula is not being interpreted as a count with negative N. Piecewise linearity and zero tails extend the maxima from these breakpoints to every real y.

For r=17+2t or18+2t, the maximum increment is bounded by5/19, with exact slack

    t / [38(19+2t)] >= 0.

Averaging over the other independent uniforms preserves this bound. Thus increasing one factor rank raises normalized curvature by at least

    3/4−(5/2)(5/19)=7/76>0.

For s≥2 factors the minimum is attained when all ranks are17. Cauchy–Schwarz and variance323s/12 give normalized curvature strictly larger than3s, because

    48[(39s/4)²−(25/4)(323s/12)] = s(4563s−8075),
    4563(2+t)−8075 = 1051+4563t > 0.

Both compared quantities are positive before squaring. For a single factor, and for the finite-vertex cost, the parity identities are

    8[beta_(30+2t)−3(30+2t)] = 4t²+60t,
    8[beta_(31+2t)−3(31+2t)] = 4t²+64t+36.

These coefficientwise identities prove the unbounded scalar comparisons. Eighteen literal discrete-uniform ranks corroborate the formulas but do not supply their all-rank scope. No ambient163 cutoff is used.

Finally, the shape inequality has the positive additive term5(n−1). Consequently equality in the local cost coefficient at n=30 still yields a strict global contradiction. Combining the complete payments with the shape inequality excludes every n≥30 under the stated minimum-factor-rank hypothesis.

The final producer hashes bound by the certificate are PY `77c6692ab69f5dfa332ebc15bf886c406cd74b30edd2e2cf5c13c45f31aa35bd`, JSON `197f5620927cfbed9d974a4ee6d590395b1248bded3d8289d1d846766f136e14`, and MD `cd378191653bffd4a3a340d8b4de920ec997b3b05e55dd6c24799fa255a99816`. Generation and a fresh process running `--verify` both terminate with exit0/PASS. Syntax, direct hashes and local report links pass. Both carrier proofs were frozen before binding.

This is internal mathematical review plus an independent small exact scalar certificate. It is distinct from external refereeing, formal verification, a new example, or a publication-priority assessment.
