# Independent review of the cap-parameter exclusion

The [parameter diameter proof](fresh_h9_h8_facet_cover_20261014_parameter_gate.md) and [whole-line exclusion](fresh_h9_h8_facet_cover_20261014_line_exclusion.md) pass source and mathematical review. No domain or endpoint gap was found.

The [separate checker](../scripts/fresh_h9_h8_facet_cover_20261014_line_review.py) and [review data](fresh_h9_h8_facet_cover_20261014_line_review.json) import no producer arithmetic. They use a separate 192-bit rational interval evaluator of the physical quartic coordinates, and SymPy only for two general symbolic derivative identities. Generation and saved replay pass. The audit independently verifies both signs bracketing the future-null crossing, positivity of all 24 derivative endpoint values, and all six strict scalar margins at the lower and split parameters.

The quadratic cap norm is strictly increasing on the segment: its leading coefficient is negative while its derivative at one is positive. The checked null bracket therefore isolates its only crossing there. Every old product remains negative, and the exact positive linear numerator makes each normalized magnitude γ decrease wherever the cap is all-dotted. This also decreases each projected norm κ². For distinct old walls, the center-distance formula has negative partial derivatives in the γ variables because their mutual products are nonpositive; ball radii increase as well.

The parameter checker computes the cosh of each ball-diameter bound using algebraic operations and outward intervals. Its strict lower-end comparison excludes every earlier all-dotted cap; its upper-end center-distance comparison excludes every later one. Convexity of hyperbolic distance makes the complete compact source diameter equal to the maximum over its 41 vertices. The finite-volume argument only needs the necessary closed-ball coverage, not a sufficiency claim at the boundary.

The whole-line split is valid at both ends. Before and including `T=0.9967505`, wall 7 lies strictly in the forbidden norm interval `(a,2a)`. At and after T, four distinct old centers have projected norm below `4a`; the source gap lemma forces them onto only three vertex rays. The center-distance formula proves distinctness uniformly, including originally right old-wall pairs. Both arguments include T, and the outer exact diameter exclusions cover the rest of the all-dotted segment.

Fresh saved primary replays also pass. Those replays reuse the producers and are not an independent repetition of the 820 source-distance comparisons. The source contact-gap census is an explicit dependency, not repeated by this review. Its convex-decomposition proof was separately read by the frontier agent: fixing an exact zero pattern prevents an exceptional vertex from entering another pattern; reverse-Lorentz superadditivity and the future recession cone give the minimum-at-a-vertex argument.

The conclusion remains specific to the fixed eleven-facet compact H8 source, perpendicular attachment, all twelve retained old walls, and the specified H–E cap segment. It is not an exclusion of compact H9 polytopes or of all lower-dimensional source attachments.

```sh
python scripts/fresh_h9_h8_facet_cover_20261014_line_review.py --verify
```
