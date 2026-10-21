# Independent review of the two additional case5 gates

**PASS.** No mathematical issue was found in the [strong-contact and two-dense-leaf arguments](fresh_h9_dotted_central_20261014_additional_gates.md). Both are sufficient local all-label obstructions. Neither excludes all label assignments of case5 by itself.

The separate [checker](../scripts/fresh_h9_dotted_central_20261014_extra_gates_review.py) and [certificate](fresh_h9_dotted_central_20261014_extra_gates_review.json) use only standard-library rational arithmetic and a small independently implemented polynomial ring; they import no producer arithmetic. Source hashes bind both gate scripts/data and the literal incidence source, but not the prose report.

## Literal premises

The checker directly verifies all 444 saved rational vertex points of case5 against its 18 defining inequalities. Every point has precisely the declared nine active facets, whose normals have rank nine. It then independently reconstructs the 44 ordered block/leaf placements from the incidence, including the missing own triples and all required finite supports.

This reproduces the 88 strong-contact placements and 44 two-dense-leaf placements, with their distinct condition sets. All 1,100 listed containing-vertex occurrences, involving 89 distinct saved vertices, are checked directly. The condition sets and their witness-map references agree with the saved gate data. The completeness of the original incidence census is inherited; this review does not re-enumerate the original polytope from scratch.

## Schur bounds and unrestricted labels

The independently expanded three-node determinant is `1−a²−b²−c²−2abc`. Since every nonright finite Coxeter contact has cosine at least `1/2`, an edge of label at least 6 cannot meet another nonright edge in a positive-definite three-node Gram. This justifies the forced-right cross-block contacts without a label cutoff.

The four-node path determinant is independently recovered as

`(1−α²)(1−β²)−δ²`.

For the strong-contact gate, positivity of the literal `A∪{y,B_j}` support forces `y–A` right: a nonright A contact contributes at least `1/2` to the projection square, as does the specified strong `y–B_j` contact. The four-path determinant then forces the middle label to be at most 3. Eliminating the positive rank-four block gives negative residual diagonals and determinant at least `√2/2−1/4>0`.

For the dense-pair gate, the minimum dense projection square is independently obtained from the two-by-two inverse as `1/[2(1−cos θ)]`. The indicated block labels give the strict bounds `H>3` and `Y>3/2`. If there is no foreign contact, the four-path bound gives a residual determinant greater than `1/2`.

If a foreign contact exists, positivity forces the B label to be 4 or 5 and exactly one foreign contact, of label 3. The literal touched-root triangle forces `x–y` right. The foreign square satisfies `F<3/4`. The checker independently verifies the identity

`det(residual)−[(H−1)(Y−1)−F] = FY−K²`,

where `K` is the inner product of the two B-plane projections. Positive-metric Cauchy makes the right side nonnegative, giving a strict determinant bound greater than `1/4`. Both residual diagonals are negative. Sylvester inertia therefore supplies two negative directions, impossible in a Lorentzian Gram matrix.

All strict radical comparisons used in these bounds are rechecked by exact rational squaring. Unlisted labels and dotted distances remain free. The proof uses no incidence automorphism and does not depend on a truncated angle list.

Generation and saved replay of the independent checker pass. Fresh saved replays of both primary gates also pass. No full compact-H9 search traversal is repeated or certified by this note.

```sh
python scripts/fresh_h9_dotted_central_20261014_extra_gates_review.py --verify
```
