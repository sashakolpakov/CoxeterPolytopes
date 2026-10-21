# Independent review of the three-large-block obstruction

**PASS.** The four predicates in the [primary theorem](fresh_h9_star_three_large_blocks_20261015.md) follow from their literal compact-face premises. No own-contact labels or noncore labels need to be fixed.

The separate [checker](../scripts/fresh_h9_star_three_large_blocks_20261015_review.py) and [certificate](fresh_h9_star_three_large_blocks_20261015_review.json) import no producer code. They reconstruct the complete cross/foreign pair coverage of all four placements, check all 96 path occurrences against saved actual vertices, check the 12 missing triples and four actual leaf triangles, and independently verify the positive kernel and Schur determinant using exact SymPy algebra. The prior primal realization and old parameter banks are not replayed.

The normalized 6–3 path has positive null vector ((\sqrt3,2,1)). A hypothetical nonright trial contact gives an actual-face Gram entrywise no larger than that matrix. Its quadratic form on the same positive vector is therefore nonpositive, contradicting the positive definiteness required at a compact vertex. Equality at label 6 causes no exception. This forces every claimed cross-root and foreign-leaf contact to be right.

The resulting six-root Gram is positive definite as an orthogonal sum of three finite dihedral blocks. Those six roots need not themselves be contained in one saved vertex: positive definiteness follows algebraically from the verified pair relations. Each own triple is minimally nonfacial. Its proper pairs are spherical; the whole triple cannot be spherical, because acute spherical wall sets meet in an actual face, and cannot be parabolic in a compact chamber. Consequently its determinant and its Schur residual square are strictly negative.

Projecting the three leaves off the six-root positive subspace therefore gives three timelike residual vectors. Their projections onto distinct own blocks are orthogonal. If an ambient leaf pair were right, its two timelike residuals would be orthogonal, giving two negative directions in a Lorentzian space. Thus every leaf pair is nonright. All three pairs belong to the verified actual leaf triangle, so their finite Coxeter labels are at least 3. The triangle Gram evaluated on the all-ones vector is at most (3-2(1/2+1/2+1/2)=0), contradicting positive definiteness.

This proves at most two of the four indicated internal labels can be at least 6 in a compact realization of the saved incidence. It does not rule out that entire incidence, and it supplies no occurrence assertion for arbitrary high-dimensional polytopes.

Generation and saved replay pass:

```sh
python scripts/fresh_h9_star_three_large_blocks_20261015_review.py --verify
```
