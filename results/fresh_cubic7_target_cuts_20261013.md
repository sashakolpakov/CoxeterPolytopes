# Targeted cubic-seven H9 cuts: symmetry obstruction, replacements, and scalar mirrors

No compact or finite-volume H9 example was found. New reflection mirrors cut successive certified escapes, but all tested compatible cones still have exact spacelike witnesses. The bounded searches below explore different support patterns in one host, rather than uniformly increasing an ambient Vinberg cutoff. A separate exact projection identity rules out one entire family of compatible cuts at every height.

The host is

\[
G=\langle-2\tau\rangle\oplus B^4\oplus\langle1\rangle,
\quad B=\begin{pmatrix}2&-C\\-C&2\end{pmatrix},
\]

where `τ=2cos(2π/7)`, `τ³+τ²−2τ−1=0`, `C=τ²+τ−1`, `D=2−τ²`, and `β=2−τ`. The identities `C²=2+τ`, `D²=2−C`, and `1/D=1+τ` are exact. The form has signature (9,1) at the chosen embedding and is positive definite at the other two embeddings. The starting point is the [25-wall height-30 prefix](fresh_cubic7_dihedral_h9_2theta_h30_20261013.json), whose zero-height chamber is `I2(14)^4+A1`. Coordinates are time, four ordered binary blocks, and a final scalar.

## A geometric obstruction tells us which symmetry must be broken

The original spacelike escape is

```
v = (1,0,0,0,0,1,1,(-10+10τ+8τ²)/7,(24+4τ−8τ²)/7,0),
v² = (114−44τ−38τ²)/7 > 0.
```

Let `P` average only coordinates 7 and 8. This is orthogonal projection onto the hyperplane of original wall 7, the reflection swapping those coordinates. With the zero-based wall indices of the original file, the exact identity is

\[
P(v)=(\tau-1)w_{13}+(2\tau-2)w_4+(\tau^2-1)w_5.
\]

All coefficients are positive. All three supporting walls lie in `w7⊥`. For any real new normal `r` with equal coordinates 7 and 8, and nonpositive pairing with every retained original wall, this implies `(r,v)=(r,P(v))≤0`.

Thus **no full-dimensional finite-wall completion using only such added normals can have finite volume**, regardless of their other coordinates, field, height, or reflection integrality. Any cut of this particular escape must break that single block-swap symmetry. This is not a nonreflectivity theorem for the full host and does not apply after discarding its supporting original walls.

The [certificate](fresh_cubic7_block_symmetry_obstruction_20261013.json) and [Fraction-only checker](../scripts/fresh_cubic7_block_symmetry_obstruction_20261013.py) verify the projection, positivity, original inequalities, and escape norm. A deliberately changed coefficient is rejected.

Before the stronger obstruction, a complete height-100 search with all four blocks equal in their two coordinates found 18 individually compatible new mirrors and no cut of this escape. Its [artifact](fresh_cubic7_equal_block_cuts_h100_20261013.json) is retained as a calibration. Raising that particular bound cannot solve the volume problem.

## Exact bounded searches that break the symmetry

The full primitive norm classes, up to square units, are `{1,2,β,2β}`. In the scalar-zero slice evenness reduces these to `{2,2β}`. Norm classes 4 and 4β would be imprimitive. All candidate roots must pass integrality of the reflection coefficients, every retained-wall inequality, and the finite-angle gate.

This last gate is essential: integral acute reflection pairs in this field can have angle `2π/7`, which is not a Coxeter angle. For a finite-angle pair the allowed values of `4cos²` in this host are exactly

```
0, 1, 2, 3, 2+τ, 2+C.
```

Values at least 4 represent parallel or disjoint mirrors. A regression explicitly rejects an integral norm-two pair with value `τ²` and accepts the corresponding `π/7` example. Mere acuteness was not used as a substitute for this test in the reported dictionaries.

The fast [thin alphabet enumerator](../scripts/fresh_cubic7_target_alphabet_20261013.py) uses all three real embeddings. Certified rational root intervals and the exact trace-dual basis bound each integral coefficient; for each pair of the last two coefficients it intersects the three permitted intervals for the first. Every survivor passes exact field-sign checks. Four comparisons against the original complete coefficient-box enumerator pass.

| Support and retained base | Height | New dictionary roots | Tested final wall counts | Result |
|---|---:|---:|---:|---|
| One asymmetric block + one equal block, scalar zero; original 25 | 1000 | 48 | 73 | Exact infinite-volume escape |
| Two arbitrary blocks, scalar zero; retain the previous 73 | 1000 | 204 | 158, including separate Pell cap | Exact infinite-volume escape |
| Two arbitrary blocks, scalar zero; rebuild from original 25 | 1000 | 258 | 122, including separate Pell cap | Exact infinite-volume escape |
| One equal block + scalar; original-25 dictionary | 1000 | 25 | 182 or 146 after completion/replacement | Exact infinite-volume escape |
| One full block + one equal block + scalar; original-25 dictionary | 1000 | 663 | 742 or 706 after completion/replacement | Exact infinite-volume escape |

These counts describe dictionaries of individually compatible mirrors, not new polytopes. Distinct dictionary roots can conflict, and the selected compatible subsets are not asserted to have maximum cardinality. The bounded support searches are not complete ambient height searches.

The first asymmetric search found no new norm-two compatible roots and 48 norm-2β roots, all mutually compatible. One of the latter cuts the original escape. It has time `2+7τ+3τ²`, equal coordinates of the same value in block 2, and block 3 coordinates `(6+14τ+6τ²,4+16τ+8τ²)`. Its height is `61+171τ+76τ²` and its escape pairing is `−8+4τ+2τ²>0`. All remaining blocks and the scalar vanish. The resulting 73-wall cone has another exact spacelike escape. See the [dictionary and certificate](fresh_cubic7_asymmetric_cuts_h1000_20261013.json).

Allowing two full blocks finds useful earlier mirrors that conflict with some of the 48 trial caps. Rebuilding from the original 25 walls retains 36 of those caps, selects 96 dictionary mirrors altogether, and gives 122 walls after a separate scalar Pell cap. Keeping all 73 earlier walls instead gives 158 walls. This replacement is necessary to avoid treating premature trial mirrors as compulsory walls of every possible completion. The [two-block report](fresh_cubic7_two_block_report_20261013.md) gives both full bounded enumerations and their common escape.

## Scalar cuts and rejection of a premature Pell cap

For one equal block `(a/D,a/D)` and scalar `z`, the norm equation is

\[
2a^2+z^2-2\tau x^2=q.
\]

This integral slice is saturated. Primitive β-divisible norms are impossible: reflection integrality would force `x,a,z` all divisible by β, contradicting norm valuation one. The complete normalized primitive pool is therefore `{1,2}`. All four placements of the equal block are enumerated through height 1000.

There are 25 new original-25-compatible mirrors: four of norm 1 and 21 of norm 2. Seven cut the common escape of the 158- and 122-wall cones. Twenty-four fit both retained cones; they give 182 and 146 walls. A third run rebuilds directly from the original 25 using the union of the two-block and scalar dictionaries. It selects a pure-scalar norm-two mirror of height `5+14τ+6τ²` and correctly rejects the older norm-one Pell cap, which conflicts with the selected walls. This alternative also has 146 walls.

All three choices have another exact escape, with norm

\[
(-2885-26821\tau+26557\tau^2)/5041>0.
\]

Its support uses one equal block, one asymmetric block, and the scalar. The [scalar artifact](fresh_cubic7_scalar_h1000_20261013.json) saves all three cases, the complete declared dictionary, precise cross-compatibility flags, and rejection of the old Pell cap. The Pell cap has height `648+1816τ+808τ²`, above 1000, and was never counted as part of a height-1000 dictionary.

## A full block, an equal block, and the scalar together

The next search allows exactly that more general support, including zero and equal free blocks. It restores the full primitive norm pool `{1,2,β,2β}` and solves

\[
E(y_1,y_2)+2a^2+z^2-2\tau x^2=q
\]

by an exact energy join over all twelve ordered block placements. The same height-1000 bounds apply to all three embeddings; this is a new support family rather than a larger cutoff. There are 663 original-25-compatible dictionary mirrors: 32 of norm 1, 559 of norm 2, none of norm β, and 72 of norm 2β. Seventy-six cut the previous common escape.

Each of the three preceding cones admits 560 additional mirrors, all mutually Coxeter-compatible. They give one 742-wall and two 706-wall cones. Rebuilding from the original 25 using the union with the previous two-block dictionary gives a fourth 706-wall cone. The combined dictionary has 855 distinct mirrors; no claim is made that it is a full ambient height-1000 dictionary.

All four cones have an exact common escape with norm

\[
(53964-19648\tau-18749\tau^2)/5041>0.
\]

The [mixed scalar artifact](fresh_cubic7_mixed_scalar_h1000_20261013.json) records every root, source index, compatibility flag, selected wall, interior point, and escape; the [dedicated report](fresh_cubic7_mixed_scalar_report_20261013.md) gives reproduction commands and the complete census. The producer completed in approximately 398 seconds, including all-pair final angle checks and geometry. The larger wall counts still do not establish a finite-volume example.

The final check is stronger than failure of four selections: that same escape satisfies **all 855 mirrors in the combined dictionary, all 25 original walls, and the separate Pell cap**. A common strict timelike center is also certified. Hence every compatible subset of this entire finite bank still has infinite volume. The union itself need not be Coxeter; it is used only as an overinclusive system of inequalities. The [common-escape certificate](fresh_cubic7_mixed_scalar_common_escape_20261013.json) and [independent replay](fresh_cubic7_common_escape_independent_audit_20261013.json) establish this bounded-dictionary result. They do not exclude larger heights or support patterns outside the searched union.

## Independent verification and limits

The [independent target verifier](../scripts/fresh_cubic7_target_independent_verify_20261013.py) imports only a separate standard-library/Fraction field implementation, with no producer arithmetic or numerical geometry. It reconstructs the Gram form; verifies every saved root, norm, height, reflection coefficient and support; checks all retained and final Coxeter angles; verifies selected-wall provenance and scalar replacement flags; and replays each strict timelike interior and spacelike witness. Every witness is also perturbed to a strictly feasible spacelike vector, proving an open infinite-volume end. Forged escape norms are rejected.

The independent replay includes the equal-block, asymmetric, retained two-block, replacement two-block, three scalar, and four mixed scalar cases. Its [record](fresh_cubic7_target_independent_audit_20261013.json) identifies each completed check. This verifier does not independently repeat the entire positive-height candidate enumeration; the stated completeness is specific to the producers' embedding bounds and support patterns. A separate [completeness review](fresh_cubic7_two_block_completeness_audit_20261013.md) checks those bounds and primitive norm reductions, independently brute-forces six small alphabet cases, and compares the full host enumerator with a separate two-block enumeration through height 30. All 20 roots agree. It also reviews the scalar joins and replacement provenance.

```sh
PYTHONPATH=src:. python scripts/fresh_cubic7_target_independent_verify_20261013.py \
  results/fresh_cubic7_equal_block_cuts_h100_20261013.json \
  results/fresh_cubic7_asymmetric_cuts_h1000_20261013.json \
  results/fresh_cubic7_two_block_h1000_20261013.json \
  results/fresh_cubic7_two_block_replacement_h1000_20261013.json \
  results/fresh_cubic7_scalar_h1000_20261013.json \
  results/fresh_cubic7_mixed_scalar_h1000_20261013.json \
  --output results/fresh_cubic7_target_independent_audit_20261013.json
python scripts/fresh_cubic7_block_symmetry_obstruction_20261013.py --verify
python scripts/fresh_cubic7_target_independent_verify_20261013.py \
  results/fresh_cubic7_mixed_scalar_common_escape_20261013.json \
  --output results/fresh_cubic7_common_escape_independent_audit_20261013.json
PYTHONPATH=src:. python -m pytest -q tests/test_fresh_cubic7_target_cuts_20261013.py
```

The whole host remains unresolved. Neither wall counts nor repeated successful cuts of individual witnesses establish finite volume; a complete positive volume certificate would be required for any claimed new-dimensional example.
