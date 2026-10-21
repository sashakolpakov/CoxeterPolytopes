# Independent review of the monotone cut-contact family

**PASS.** No compact H9 Coxeter realization exists in the stated family: the four internal B2 labels remain 4, each positive leaf has its two own contacts labeled 3, and both own contacts of every cut facet are nonright. This is not an exclusion of the star incidence. Own right contacts, changed positive-leaf labels and changed B2 labels remain outside the argument.

The [independent checker](../scripts/fresh_h9_star_cube_metric_20261014_cut_profiles_review.py) and [certificate](fresh_h9_star_cube_metric_20261014_cut_profiles_review.json) import neither SymPy nor producer arithmetic. They use Fraction arithmetic in the basis `1,sqrt(2),sqrt(5),sqrt(10)`, independently of the producer's algebraic-number field. Every primary determinant is converted from its primitive-element encoding and compared exactly.

The review verifies all eight new support guards in both own-root orientations. If a cut has an own label at least 4, the actual eight-root spherical support has own projection cost at least `1/2`; one foreign label-3 contact would add at least `1/2`. Positive definiteness therefore forbids every foreign nonright contact. A `(3,3)` cut retains the seven earlier patterns. The actual three-node paths consisting of a positive leaf, one own B2 root and its cut have determinant `1-1/4-3/4=0` at label 6; monotonicity excludes all larger labels. Thus cut own labels reduce exactly to 3, 4 or 5 within this family.

The 30 earlier forced-right and anchor-angle guards are inherited from the [separate first-star review](fresh_h9_star_cube_metric_20261014_review.md). Every changed contact increases its cosine magnitude, so none can restore positive definiteness to a failed relaxed vertex Gram. No unproved finite-label cutoff is used.

The independent census reproduces all **81 profiles and 6,300 residual cases**, with 765 distinct determinant values. There are 105 distinct nondegenerate anchor cases. Exactly 6,296 residual determinants are nonzero, giving rank at least eleven after restoring the positive rank-eight block.

The four zero determinants are not discarded. They have `C0=(3,3)`, `C3=(4,4)`, with the two foreign contacts choosing one root from each of their B2 blocks. Each forces a four-cycle with four label-3 bonds and two right diagonals, on one of the actual supports:

| Four-root support | Actual vertex index |
|---|---:|
| 3,11,12,14 | 162 |
| 3,11,13,14 | 163 |
| 4,11,12,14 | 189 |
| 4,11,13,14 | 190 |

Its normalized Gram annihilates the positive all-ones vector. This affine A3 principal block cannot lie in a positive-definite Gram at a compact vertex. The checker reconstructs all four incidences and null relations separately and compares the primary witnesses.

The final result uses both the rank obstruction and these four spherical-vertex obstructions. It would be incorrect to describe the rank test alone as exhaustive exclusion, or to extend the vertex argument to arbitrary noncompact realizations. The labels of the other two cuts may vary within the stated nonright family because they are unnecessary for this principal-block contradiction.

Generation and saved replay pass. Source and data hashes are bound; prior frozen reports and arithmetic files were not modified.

```sh
python3 scripts/fresh_h9_star_cube_metric_20261014_cut_profiles_review.py --verify
```
