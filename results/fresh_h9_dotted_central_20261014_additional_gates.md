# Two complementary block obstructions in case5

These two all-label gates use the same literal four-root block and four-path supports as the [dense/high gate](fresh_h9_dotted_central_20261014_gate.md). Every support is explicitly contained in a saved actual vertex of case5. They are local Lorentzian negative-index obstructions; neither claims that the entire incidence is excluded.

All roots below have square1. Let A and B be disjoint dihedral blocks, and let x and y be leaves whose own missing triples are A∪{x} and B∪{y}. The saved supports include A∪B, B∪{x}, A∪{y}, and each {A_i,x,y,B_j}. Inverse entries of a positive dihedral Gram matrix are nonnegative.

## A strong contact on the second leaf

Assume the A label is at least4, the B label at least6, both x–A contacts are nonright, and one specified y–B_j contact has label at least4. The [checker](../scripts/fresh_h9_dotted_central_20261014_strong_contact_gate.py) verifies **88 placements and88 condition sets**, including the additional actual support A∪{y,B_j}.

The label at least6 in B forces every A–B cross contact and every x–B contact to be right. If y had a nonright contact with A, its projection square onto A would be at least `(1/4)/(1−1/2)=1/2`. The specified y–B_j contact contributes at least another1/2. This contradicts positivity of the actual support A∪{y,B_j}. Thus y–A is right as well.

The four-path determinant is at most `3/8−δ²`, where δ is the x–y cosine magnitude. A label at least4 fails, so δ²≤1/4. The x projection square is at least `1+√2/2`, and the y projection square is at least2. Their residual squares are negative; their residual determinant is at least

`√2/2−1/4 > 0`.

This produces two negative directions after eliminating the positive rank-four block.

## Dense leaves on both blocks

Assume instead that the A label is at least6, the B label at least4, and both own contacts of each leaf are nonright. The [second checker](../scripts/fresh_h9_dotted_central_20261014_dense_pair_gate.py) verifies **44 placements and44 condition sets**, including the two actual triangles {x,y,B_j}.

The high A label forces all A–B and y–A contacts right. Write H for the own projection square of x onto A and Y for that of y onto B. The elementary bounds

`cos(π/6)>5/6`, `cos(π/4)>2/3`

give H>3 and Y>3/2, respectively.

If x–B is right, the original four-path bound gives δ²≤1/2. The residual determinant is greater than `(3−1)(3/2−1)−1/2=1/2`.

Otherwise positivity of B∪{x} forces the B label to be4 or5, and x has exactly one foreign contact, of label3. Two contacts would make a cycle; a foreign label at least4 gives a nonpositive three-node path. The actual triangle formed by x, y and the touched B root forces x–y right. Let f and v be the B projections of x and y, and set F=f². Since `cos²(π/5)<2/3`,

`F = 1/[4sin²(π/b)] < 3/4`.

The residual determinant is

`(H−1+F)(Y−1)−(f,v)²`.

Cauchy's inequality in the positive B plane gives `(f,v)²≤FY`. Consequently this determinant is at least `(H−1)(Y−1)−F`, which is strictly greater than1/4. Both residual diagonal entries are negative, again contradicting Lorentzian index1.

The three trigonometric comparisons reduce to the exact rational square inequalities `3/4>25/36`, `1/2>4/9`, and `5<49/9`. No numerical angle estimate or finite cutoff is used. All other contacts and all dotted distances remain free.

Both sources and JSON certificates are frozen; generation and saved replay pass:

```sh
python scripts/fresh_h9_dotted_central_20261014_strong_contact_gate.py --verify
python scripts/fresh_h9_dotted_central_20261014_dense_pair_gate.py --verify
```

The [independent review](fresh_h9_dotted_central_20261014_extra_gates_review.md) reconstructs all444 literal vertices, the88 and44 condition sets, and1,100 support occurrences. It separately checks the path, triangle and residual determinant formulas with standard-library exact arithmetic. Both primary saved replays and its independent saved replay pass; no issue was found.
