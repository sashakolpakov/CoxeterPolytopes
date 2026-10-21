# Independent review of the complete fixed-B2 cut-profile branch

**PASS.** The [primary theorem](fresh_h9_star_cube_metric_20261014_all_cut_profiles.md) excludes compact H9 Coxeter realizations of star case 0 when the four internal B2 edges have label 4 and every positive leaf has own contacts `(3,3)`. The cut contacts are unrestricted within this statement. It does not exclude changed B2 labels, changed positive-leaf contacts or the entire star incidence.

The [review script](../scripts/fresh_h9_star_cube_metric_20261014_all_cut_profiles_review.py) and [certificate](fresh_h9_star_cube_metric_20261014_all_cut_profiles_review.json) use the preceding independent Fraction field engine in the basis `1,sqrt(2),sqrt(5),sqrt(10)`. They import no producer or SymPy arithmetic. Signs use rational square-root enclosures obtained by integer square roots; determinant zero tests are exact coefficient comparisons.

The new right-contact profiles require fresh proofs, and those checks pass:

- All **24 cross-B2 forced-right witnesses** are recomputed using only the fixed B2 and positive-leaf labels, with cut contacts relaxed to right.
- All **44 positive-leaf/cut-pair witnesses** are checked on their literal three-root faces. Before imposing a label bound, a Lanner cut triangle has some own label at least 3. Testing a nonright positive-leaf/cut contact gives an all-ones quadratic value at most `1-2cos(pi/m) <= 0`, so that pair must be right.
- Both orientations of all **eight own-root path/support guards** are verified. With the preceding right contact established, an own label 6 makes the actual `3–6` path determinant zero; every larger label also fails. The Lanner condition then leaves exactly the eleven ordered profiles, including `(2,5)` and `(5,2)`.
- The inverse matrices of the actual eight-root supports justify all foreign-contact budgets. The `(3,3)` profile has seven possibilities. Every other allowed profile has an own label at least 4 and therefore permits no foreign nonright contact.
- All **242 fresh anchor-label upper witnesses** are independently recomputed by ordinary exact row elimination and leading determinants. Each uses only the two specified cut profiles; every unused cut is relaxed to right. No old witness is transferred monotonically through a decreased contact.

The independent census checks every one of the **8,092 residual determinants** and compares it with the primary algebraic-field encoding. There are 119 distinct nondegenerate anchor parameter cases. Exactly **8,088** residual determinants are nonzero, forcing full principal rank at least eleven in an ambient dimension of ten.

The remaining four cases are precisely the earlier `(3,3)/(4,4)` cases. Each has an actual four-root affine A3 cycle with a positive all-ones null vector. The review reconstructs these four supports and also independently verifies the primary full-vertex failure minors. Thus they cannot occur at a compact vertex. No rank-zero case is silently dropped.

The proof order avoids circularity: first force the positive-leaf/cut pair right without any upper label bound; next exclude labels at least 6; then enumerate the eleven finite profiles and their residuals. The final obstruction combines rank and compact-vertex positivity. It is not a claim about arbitrary noncompact realizations.

Generation and saved replay pass. The first-star and 6,300-case review files remain frozen. The new positive-leaf `(3,4)` assignment recorded by the continuation is outside this theorem and remains a separate metric problem.

```sh
python3 scripts/fresh_h9_star_cube_metric_20261014_all_cut_profiles_review.py --verify
```
