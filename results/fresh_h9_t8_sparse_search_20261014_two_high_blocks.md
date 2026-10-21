# Two high blocks force two negative directions in T8

The literal case-8 incidence has no compact H9 Coxeter realization under just these four conditions on finite integer labels:

```text
(1,2) >= 7,  (4,5) >= 7,  (7,8) >= 4,  (12,13) >= 4.
```

Every other label and dotted distance is unrestricted. This excludes the [latest saved mandatory assignment](fresh_h9_t8_sparse_search_20261014_next_case_08.json). It does not exclude the whole incidence.

The [standalone Fraction checker](../scripts/fresh_h9_t8_sparse_search_20261014_two_high_blocks.py) and [certificate](fresh_h9_t8_sparse_search_20261014_two_high_blocks.json) pass generation and saved replay. The [independent review](fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.md) separately rebuilds the literal supports, positive-null matrices and all-label cost argument; generation and saved replay also pass. There is no new polynomial parameter bank: positive null vectors, a uniform Lannér projection-cost bound and a two-by-two negative-definite residual suffice.

## All contacts used in the argument follow from the four conditions

Normalize all wall normals to square one. The actual eight-wall support

```text
S = [1,2,4,5,7,8,12,13]
```

splits into four orthogonal rank-two blocks. Twenty cross-pair witnesses use a three-node path whose internal label is at least seven and whose tested cross contact has label at least three. Lowering the internal label to six gives a positive null vector `(sqrt(3)/2,1,1/2)`; the actual internal label makes its quadratic value negative. The remaining four cross-pair witnesses use the affine `4-3-4` path. Every witness support is contained in an explicitly checked actual vertex, so these failures force every cross-block contact right.

The blocks themselves are positive definite, since their labels are finite. Hence S is positive definite of rank eight. The checker also verifies that every pair among `S + {11,14}` is a meeting pair. No parallel or dotted relation is inserted into this positive-support calculation.

The three literal missing triples are

```text
{1,2,14}, {4,5,14}, {11,12,13}.
```

Compactness makes each of their Gram diagrams Lannér. Thus wall 14 has a nonright contact to each high block, and wall 11 has a nonright own contact to the last block.

For wall 11, retain an active own root from `{12,13}` and delete its partner from S. Both choices give actual eight-wall subsets after adding 11. A nonright contact to either high block costs more than one in the Schur projection, because

`cos(pi/a)^2 / sin(pi/m)^2 >= 1/[4 sin(pi/m)^2] > 1`

for `a>=3`, `m>=7`. The positive own-contact cost makes this incompatible with a positive-definite vertex. Wall 11 is therefore right to both high blocks.

For wall 14, choose one active root from each high block, remove their two partners from S, and add 14. All four resulting seven-wall supports are actual faces. The two retained contacts cost at least `1/4+1/4`. Any nonright contact into either remaining full block costs at least `1/2`, because its internal label is at least four. The total reaches the forbidden Schur value one. Consequently wall 14 is right to roots 7, 8, 12 and 13.

Finally, for each choice of active high roots `i` and `j`, the four-wall set `{11,14,i,j}` is an actual face. Testing label at least four on `11–14` and lowering the two active contacts to label three gives an affine star with positive null vector `(sqrt(2)/2,1,1/2,1/2)` in order `(11,14,i,j)`. Larger contacts only decrease its quadratic value. Thus the mutual label of 11 and 14 is two or three.

These arguments check every possible active-root orientation; no symmetry of the incidence is assumed.

## A uniform projection cost greater than five quarters

Let `c=cos(pi/n)` be the cosine on the last block, `n>=4`, and let `x,y` be wall 11's own cosine magnitudes. Its projection cost on that block is

`U = (x^2+2cxy+y^2)/(1-c^2)`.

The missing triple forces `U>1`. The discreteness of finite Coxeter labels strengthens this to `U>5/4` uniformly:

| Own contacts | Necessary minimum | Lower bound for U |
|---|---|---|
| Both nonright | `x,y>=1/2` | `1+sqrt(2)/2 > 5/4` |
| Single, `n=4` | Own label at least 5 | `(3+sqrt(5))/4 > 5/4` |
| Single, `n=5` | Own label at least 4 | `1+sqrt(5)/5 > 5/4` |
| Single, `n=6` | Own label at least 4 | `2 > 5/4` |
| Single, `n>=7` | Own label at least 3 | `1/[4(1-c^2)] > 5/4` |

For `n=4`, an own label at most four would make the triple positive semidefinite; for `n=5`, own label three is positive definite; for `n=6`, own label three is affine. Those options are forbidden for a missing triple of a compact chamber.

The last row uses `cos(pi/7)^2>4/5`. The checker independently regenerates the Chebyshev polynomial `U6`; after setting its argument squared to `z`, the polynomial is `64z^3-80z^2+24z-1`. It has value `-29/125` at `4/5` and is strictly increasing above `3/4`. Since `cos(pi/7)^2>3/4`, its root lies above `4/5`.

Thus the residual of wall 11 after projection off S has square strictly less than `-1/4`; any additional projection onto the `(7,8)` block only decreases it further.

## The negative plane

Each high block contributes more than one to wall 14's projection cost. Its residual therefore has square strictly less than `-1`.

Wall 11's support projection uses only the two low blocks, whereas wall 14's uses only the two high blocks. Their projections are orthogonal. Hence their residual mutual pairing equals their original normalized pairing: either zero or `-1/2`.

The two residual diagonal entries are strictly negative, and their determinant is strictly greater than

`(1/4)*1-(1/2)^2 = 0`.

Their residual Gram is negative definite. Together with the positive rank-eight support, this gives a principal Gram with two negative directions, impossible in Lorentz signature `(9,1)`. This proves the all-label branch obstruction without evaluating any dotted distance or fixing any outside leaf contact.

```sh
python scripts/fresh_h9_t8_sparse_search_20261014_two_high_blocks.py --verify
```

The certificate binds the original incidence and the displayed saved assignment. Its new local witness arithmetic is independent of the prior polynomial producers. No later continuation is performed by this checker.
