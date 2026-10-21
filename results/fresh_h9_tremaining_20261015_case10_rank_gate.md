# T10: exact exclusion of the first necessary-label survivor

The first new T10 traversal stopped after 9,008 visits at a complete mandatory-label assignment. This certificate excludes **that entire assignment family**, including every finite integer in its `>=7` category, every omitted meeting label and every dotted distance. T10 itself remains open; the ten pending DFS frames are preserved and have not yet been resumed.

The source and saved data are [the checker](../scripts/fresh_h9_tremaining_20261015_case10_rank_gate.py) and [its certificate](fresh_h9_tremaining_20261015_case10_rank_gate.json). Generation and saved replay pass. Arithmetic uses the separately written Fraction field implementation from the earlier H20 review; this is primary arithmetic reuse, not a claim of independent review of this new gate.

The fixed positive support consists of the four root pairs `(3,4)`, `(6,7)`, `(9,10)`, `(12,13)`, with labels `4,4,m,4`, where `m>=7`. The [separate high-label gate](fresh_h9_tremaining_20261015_case10_high_bound.json), independently checked by [the root review](fresh_h9_tremaining_20261015_case10_high_bound_review.json), proves all cross-pair contacts right and gives the further needed right contacts. It excludes every `m>=12`: the residual plane of normals 11 and 16 has determinant at least `sqrt(6)/2-1>0`, with both diagonal entries negative.

For the five remaining labels 7–11, use residual normals 11, 16 and 2. Each of 11 and 2 is right to the high block. Four literal eight-node faces, obtained by deleting one own root from the positive support and inserting the corresponding leaf, prove that each leaf has at most one foreign contact, necessarily label 3. The retained own root costs `1/4`; a foreign label-3 contact costs at least `1/2`; two such contacts or one label at least 4 contradict positive definiteness. A contact to the high block costs at least 1 even after relaxing `m` to 6. Two additional literal vertex failures bound the 11–2 and 16–2 labels by 3. These give precisely a necessary relaxation of `5*5*2*2=100` residual cases.

Write `d=1/(4 sin²(pi/m))`. Only the residual square of normal 16 depends on `d`, so each three-by-three residual determinant is `A*d+B` over `Q(sqrt(2))`. All 100 slopes `A` are nonzero. Since an eight-dimensional positive span in the ten-dimensional Lorentz space leaves a two-dimensional complement, this determinant must vanish. Consequently `d=-B/A`.

Every actual integer label must satisfy the elementary Chebyshev identity

`T_m(1-1/(2*d))=1`.

The checker evaluates all 500 necessary identities exactly. Only four cases at `m=8` remain. Their residual `(11,2)` planes have negative diagonals and positive determinants, respectively `1/2+sqrt(2)/4` or `5/8+sqrt(2)/2`. Each therefore contributes two negative directions, which no Lorentzian Gram matrix admits. The Chebyshev test is only a necessary test; it does not assume a field-degree classification or omit any real branch.

The rank argument is specific to H9. The two-negative-direction conclusions themselves are dimension independent. All actual-face and contact guards are checked for this literal T10 incidence; no transfer from a different tree type or unproved global symmetry is used.
