# Review of the balanced all-radius radial gap

Source and proof review found no issue. The separate [review checker](../scripts/fresh_h20_balanced_half_20261014_review.py) and [result](fresh_h20_balanced_half_20261014_review.json) pass; they independently reconstruct the excluded-angle minimal polynomials with SymPy and check the rational radial constants and secondary half-probe inequalities. This is not a second run of the 48,216 transported power evaluations.

The [primary proof](fresh_h20_balanced_half_20261014.md) uses a valid necessary angle relaxation: an omitted cosine is zero, `1/2`, or at least `1/sqrt(2)`. Its prior continuous interval and exceptional-radius list therefore cover all finite labels and dotted/parallel relations. For the all-right and mixed A6/A7 rows, the additional short angle intervals leave the explicitly listed cosines; the nonintegral monic minimal polynomials exclude them because `2cos(pi/m)` is always an algebraic integer. This restriction concerns individual angles and does not assume that the hyperbolic reflection group is arithmetic.

The transported power formula has the correct unit/norm-two factor of two. The inverse of the retained simply-laced finite Gram is entrywise nonnegative, all target distances are positive, and each endpoint profile has spatial cost at most one. Consequently the ball power decreases with radius; evaluating it at a necessary upper bound is safe. The archived minimum `164/25` at `2z/5` also matches the direct own-vertex all-right computation `196(3/5)^2-8^2`.

The secondary strict-half-cutter statement retains the zero lower bounds separately. The five A7 profiles require radius greater than `49/8`; their remaining cosine interval contains only π/4. The mixed-pair numerator and the unequal-A7 direction-difference bounds correctly place their normalized cosines strictly between zero and `1/2`. Thus the strict-cut incompatibility does not silently treat equality as a strict cut.

Finally, `X=(2z/5,3111/200,1)` has square `1/4`, with all seed inequalities strict. Endpoint ball powers are at least `164/25`; the independently established [interior-center central-contact theorem](fresh_h20_balanced_20261014_bridge.md) gives power at least `9/25` for every ball in that separate family. Both exceed `1/4`, so the same strict future spacelike vector survives their entire finite union. This proves infinite volume for that union. It does not cover arbitrary boundary-center caps or a replaced central facet.

Replay the small independent arithmetic review with:

```sh
python3 scripts/fresh_h20_balanced_half_20261014_review.py --verify
```
