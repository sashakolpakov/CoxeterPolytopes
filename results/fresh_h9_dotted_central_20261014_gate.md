# A dense/high dihedral block obstruction in tree case5

The [exact gate](fresh_h9_dotted_central_20261014_gate.json) certifies **44 literal placements, giving24 distinct label conditions**, in the saved case5 incidence. Four conditions apply to its previously interrupted search path. This excludes a family of label assignments, not the whole incidence. All other labels and dotted distances remain unrestricted.

Normalize all roots to square1. Let A=I2(a), with a≥6, and B=I2(b), with b≥7, be disjoint blocks. Let x have nonright contacts with both roots of A, and let y have at least one nonright contact with B. The latter is necessary because B together with y is a missing Lannér triple. Require the following literal subsets to be finite faces: A∪B, B∪{x}, A∪{y}, and every {x,y,A_i,B_j}. The checker supplies an actual containing vertex for each support; it makes no incidence-automorphism assumption.

A finite connected component of rank at least3 cannot contain a label at least6: a three-node path would have determinant at most 1−3/4−1/4=0. Hence all A–B cross contacts, x–B contacts and y–A contacts are right. The four block roots have a positive definite direct-sum Gram matrix.

Choose a nonright y–B contact. On each corresponding four-node path A_i–x–y–B_j, write the three cosine magnitudes as α,δ,β. Its determinant is

`(1−α²)(1−β²)−δ²`.

Since α,β≥1/2, positivity gives δ²<9/16. A label at least5 violates this, because cos²(π/5)=(3+√5)/8>9/16. Thus x–y has label at most4 and δ²≤1/2.

The projection square of x onto A is at least

`1/[2(1−cos(π/a))] ≥ 2+√3`.

The projection square of y onto B is strictly greater than5/4. Indeed cos²(π/7)>4/5, so even one cosine of magnitude1/2 contributes more than `1/[4(1−4/5)]=5/4`. The checker verifies the Chebyshev identity

`T7(z)−1=(z−1)(8z³+4z²−4z−1)²`.

The cubic is negative at3/5 and positive at1; its largest root is cos(2π/7), proving the strict bound.

After eliminating the positive rank-four block, x and y have negative residual squares of magnitudes at least1+√3 and strictly greater than1/4. Their residual mutual product is their original product, because their block projections are disjoint. The residual determinant is therefore strictly greater than

`(1+√3)/4−1/2=(√3−1)/4>0`.

This negative definite two-dimensional residual contradicts Lorentzian negative index1. The proof covers every integer label in the indicated ranges, including the unbounded categories. The source and certificate are frozen; generation and saved replay pass.

```sh
python scripts/fresh_h9_dotted_central_20261014_gate.py --verify
```
