# Independent review of the 860-facet π/6 construction

The construction and its stated geometric scope pass mathematical source review and a separate exact arithmetic replay. The [review checker](../scripts/fresh_h22_multiquadratic_20261013_review.py) and [saved audit](fresh_h22_multiquadratic_20261013_review.json) bind the frozen [producer](../scripts/fresh_h22_multiquadratic_20261013.py), [certificate](fresh_h22_multiquadratic_20261013.json), and rational source data by SHA-256. The main [construction report](fresh_h22_multiquadratic_20261013.md) makes no finite-volume or new-dimension claim; no mathematical overstatement was found.

The separate checker imports no producer arithmetic, interval, or projection implementation. It reconstructs the rational 812-wall carrier from its literal inputs, then represents the new coordinates as rational linear combinations of `sqrt(n)` for distinct positive squarefree integers n. Multiplication reduces products by exact integer gcd identities. Exact zero and norm identities are coefficient equalities; strict signs use rational square-root enclosures obtained by integer square roots. No floating-point value enters a verdict. The previously certified pair gate among the 812 rational parents is inherited rather than duplicated.

## Construction and angles

For each original recipe, orthogonality of u and v gives

`(-(q sqrt(3)/2)u + epsilon sqrt(q/v² - 3q²u²/(4v²))v)² = q`.

The finite parent pairings follow from exact identities, not interval overlap. Contacts proportional to `sqrt(3)` have normalized cosine squared `3/4`, and the retained right and π/4 contacts remain exact. Every other pairing is checked afresh. The rejected type has cosine squared `48/49`; its corresponding `2cos(2π/m)=94/49` cannot be a rational algebraic integer, so it is forbidden for every integer Coxeter label.

The independent replay verifies all 48 signed cap norms and strict center signs. It checks the 19,488 positive-height cap–parent pairs directly: 504 right, 26 of label 6, two of label 4, and 18,956 dotted. The negative-height checks follow from the exact height-reflection symmetry of the rational parent set. All 1,128 cap–cap pairs are checked independently: 1,120 dotted, six right, and two parallel.

Positive norm, a common strict timelike point, and nonpositive cross pairings make every wall an actual facet. Specifically, projecting that point orthogonally to one normal remains timelike and is strict on every other inequality. This justifies the count of 860 facets.

## The compact quadrilateral is an actual face

The exact frame identities are

`Z²=-2`, `u²=2`, `v²=1`, `(u,v)=1`, `(Z,u)=(Z,v)=0`.

The reviewed cancellation `(Z,cap17)=0` is `sqrt(93)(3/2-3/2)=0`. It establishes the required contact with the substituted π/6 cap, rather than inheriting orthogonality from the old π/7 configuration. The two new norm-two formulas are mutually perpendicular by

`-2*209*241 + 241*418 = 0`.

All four frame points are reconstructed independently and tested against all 860 inequalities. Their squared norms are −2, −28,800, −843,292,800, and −29,281; their complete active sets each contain exactly 22 walls. Positive definiteness follows from principal blocks of the checked 21-wall positive-definite parent Gram matrix, together with mutually orthogonal positive-norm added roots. Thus these are actual finite corners, not merely feasible timelike points.

The 20 common walls define H2. Inverting the released-wall dual Gram `[[2,1],[1,1]]` gives `[[1,-1],[-1,2]]`, so the old corner has intrinsic angle π/4. The other three angles are right. All remaining inequalities contain the four finite corners and hence their compact convex quadrilateral. This certifies the actual compact face and, by reflection, its second copy. It does not establish compactness of the H22 chamber.

## Field and volume conclusions

The independent reconstruction confirms the stated multiquadratic coordinate-field upper bound and the cyclic product `208+64sqrt(3)`. Consequently the cyclic-product field contains `Q(sqrt(3))`; no assertion about the full minimal field is needed.

In the normalization with diagonal entries two, the retained 42 short rational roots have a rational principal Gram block whose all-ones value is −126. A nontrivial embedding of the cyclic-product field, extended to the Gram field, leaves this indefinite block fixed. It therefore cannot make the conjugate Gram matrix positive semidefinite, as required for quasi-arithmeticity. This obstructs a quasi-arithmetic finite-volume completion retaining those reflections. The [separate field audit](fresh_h22_extension_feasibility_20261013.md) explains the distinct pseudo-arithmetic restriction. A multiquadratic coordinate field removes the earlier cubic obstruction; it does not itself prove pseudo-arithmeticity.

The saved integer escape is independently checked to have norm 1,586, center pairing −814, and strict negative pairing with all 860 normals. It is future spacelike, so it gives an open unbounded end and proves infinite volume. The recorded old null cusp also survives every added cap. These conclusions agree with the producer and its report.

During replay, an equality-test failure was traced to tuple versus JSON-list serialization of the four vertex frame coefficients. The producer changed those coefficients to lists. This was a metadata fix; the independent norm, angle, and geometry checks had already passed.

Replay the separate checker with:

```sh
python3 scripts/fresh_h22_multiquadratic_20261013_review.py --verify
```
