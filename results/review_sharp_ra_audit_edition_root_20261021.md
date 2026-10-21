# Sharp right-angled audit edition: root review

21 October 2026. **PASS at the stated published and geometric inputs.**
Final source: `manuscripts/right_angled_sharp_dimension_bounds_publication.tex`,
SHA256 `0c29a7869fa46e4e1a542d8ae6f4dfc7e616136c7a1aad3984d941a95cf89722`.
The root read the preceding complete source and checked the exact two
successor hunks against the byte-preserved `0d607b48…` source. This is an
internal mathematical/readability review, not external refereeing.

The first hunk explains a previously stated standard inheritance fact.
At a finite relative-interior point of a ridge of the selected face,
the incident original walls are orthogonal, so their induced walls are
orthogonal. A face's vertices in the closed Klein polytope are vertices
of the original polytope. Entire ideality follows, and a finite
triangulation into simplices with finite or ideal vertices proves
finite volume in the dimensions used. No step descends to a polygon.
Here the manuscript's phrase "finite or ideal simplices" includes
simplices with both kinds of vertices. The second hunk explicitly
defines `A_j`, already used as the maximum height after j projections.
Neither changes the statements, assumptions, formulas or case splits.

The projection has future norm minus one, is strict at every unselected
wall, and lies on the actual selected facet. Intrinsic ridge normals
remain unit and orthogonal to the selected normal. The recurrence gives
the strict three-step height bound and distinct radial directions with
pairwise product below one quarter. Direct contraction of the displayed
rank-two and rank-three tensors yields the cubic positive-square sum.
The bound N<=27 is strict at equality: four vectors with all products
minus one half would have squared sum minus two. The published minimum
of 27 facets for a finite-volume right-angled six-polytope therefore
excludes dimension nine and, by face inheritance, every higher dimension.

For the ideal branch, the same projection gives heights below one.
Musin supplies only cardinality; Kolpakov supplies hyperbolic metric
24-cell rigidity. Opposite spatial normal pairs then contradict those
strict heights. No uniqueness assertion for a kissing code is used.
The external inputs remain exactly those identified in the text and
its existing published-input audits; they were not newly retrieved here.

The final source-bound build is recorded separately in
`results/sharp_ra_audit_edition_readability_build_20261021.json`:
five pages, zero warnings, PDF SHA256
`609a9fc025b43ae5c899b1b37b1051e592de6a8142f3ffd275fcd9dca360d443`.
The build establishes typesetting only. The unchanged figure keeps its
earlier visual-review scope; no new all-page visual inspection is claimed.
Prior exact checks and selected conditional Lean records retain their
original bindings. No fresh script replay, combined core PASS or
end-to-end geometric Lean verification follows from this editorial edit.
