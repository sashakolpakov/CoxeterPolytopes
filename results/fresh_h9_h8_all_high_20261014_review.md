# Independent review of the all-high attachment obstruction

The [separate reviewer](../scripts/fresh_h9_h8_all_high_20261014_review.py) reproduces the [three-pair obstruction](fresh_h9_h8_all_high_20261014.md) with a different exact algebraic-field implementation and 192-bit rational intervals. It imports neither the producer nor its quartic Fraction and 128-bit arithmetic. The [certificate](fresh_h9_h8_all_high_20261014_review.json) verifies all three strict cap-norm contradictions.

The checker reconstructs source vertices 32 and 33 from their literal incident walls, verifies their chamber inequalities and exact distance, and checks all 66 old-wall pairings. It derives the exact threshold identity

`cosh(4 asinh(sqrt(2/(sqrt(5)-1)))) = 17+8sqrt(5)`.

Thus any two balls covering these source vertices must belong to one of the three old dotted pairs. For each pair it checks the contact bound forced by coverage, literal orthogonality to a positive support, exact positive LDL pivots, and entrywise nonnegativity of the support Gram inverse. The resulting lower bounds on the cap square are respectively greater than 5.0096, 1737.9235 and 3.2160, each contradicting its normalization to two. Center-distance monotonicity and the Lorentz projection formula are independently checked symbolically.

This review excludes the entire all-high branch for coverage by the twelve frozen old-wall balls. Together with the separately reviewed low-center obstructions it closes the raw perpendicular attachment of this particular source. Further nonperpendicular covering walls, changed old walls and other source chambers remain outside that conclusion.

```sh
python3 scripts/fresh_h9_h8_all_high_20261014_review.py --verify
```
