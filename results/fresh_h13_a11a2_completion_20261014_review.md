# Independent review of the complete H13 source

The [separate checker](../scripts/fresh_h13_a11a2_completion_20261014_review.py)
and [certificate](fresh_h13_a11a2_completion_20261014_review.json) pass.
The source has 15 actual facets, 44 finite vertices, and two bounded cusp
links of quotient rank 12. All 105 wall pairs have the required Coxeter
angles: 91 right, 13 at pi/3, and one at pi/4.

This review uses a different complete vertex enumeration from the producer.
In the thirteen nonnegative corner-gap coordinates there are only two
additional inequalities. A vertex can therefore have at most two positive
coordinates. The reviewer tests the origin, all thirteen coordinate axes,
and all 78 coordinate pairs. The resulting 46 points agree exactly with
the producer's census of all thirteen-wall active sets. The two positive
inequalities bound every coordinate; all vertices lie in the closed future
Lorentz cone. Positive cusp kernel vectors and quotient rank 12 certify
finite volume at the two ideal points.

The reviewer also reconstructs the rational Gram matrix, every strict facet
witness, and all pairings of the normalized transfer in Q(sqrt2). It imports
no producer code. SymPy exact rational matrices provide a separate matrix
implementation. The 210-wall projection census is inherited rather than
repeated; intrinsic source completeness does not depend on that census.
No compatibility with an ambient H20 chamber is asserted by this review.

Generation and saved replay passed:

`python3 scripts/fresh_h13_a11a2_completion_20261014_review.py --verify`
