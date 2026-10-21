# Independent review of the E7 face in the 87-wall chamber

**PASS.** The exact face and the all-label Schur restriction are sound. A separate recursive `Fraction` implementation replays all 87 projections, the 67 future-timelike nonfacets, the actual finite apex, the rational intrinsic form, and the unique nonright profile. No producer arithmetic is imported, and signs are decided by exact recursive quadratic comparisons rather than numerical intervals.

The frozen inputs are [the E7-face certificate](fresh_h20_e7_face_20261014.json), [its producer](../scripts/fresh_h20_e7_face_20261014.py), and [the 87-wall source](fresh_h20_e8_interior_caps_20261014_exact.json). The [independent checker](../scripts/fresh_h20_e7_face_20261014_review.py) and [review certificate](fresh_h20_e7_face_20261014_review.json) bind those sources by SHA-256. Generation and fresh saved-result replay passed:

```sh
python scripts/fresh_h20_e7_face_20261014_review.py --verify
```

The seven walls `[12,13,14,15,19,33,67]` have the norm-two E7 Gram matrix: determinant two and arm lengths one, two, three. SymPy independently reconstructs its inverse. All 49 entries are positive, and its diagonal is

`[7/2, 12, 6, 2, 15/2, 3/2, 4]`.

Projecting the saved strict center to their orthogonal complement gives a future timelike relative-interior point, strictly inside every other inequality. Thus this is an actual codimension-seven face. The only positive-square projected walls are `[0,1,2,3,4,5,6,7,8,9,10,18,83]`; they were already orthogonal to E7. Their Gram matrix is `A11 + A2`, positive definite with determinant 36. The seven E7 projections vanish, and all remaining 67 projected normals are strictly future timelike. Those 67 inequalities are redundant throughout the future hyperbolic face. The thirteen independent remaining normals each give an actual facet: a sufficiently small perturbation of the common apex can make any one pairing zero and the other twelve negative.

The normalized timelike vector `T` is independently checked from `V = s T`, where

`s = (-3204 - 153 sqrt(79) + 5537 sqrt(11) + 356 sqrt(869))/24 > 0`.

Its square is `-2`, so the intrinsic form is `diag(A11,A2,-2)`, with determinant `-72` and signature `(13,1)`. The embedding basis has coefficients in `Q(sqrt(79),sqrt(11))`. All 87 saved reduced projection coordinate vectors are independently reconstructed and mapped back to the ambient vectors; this does not assert that every reduced coefficient is rational or establish an arithmetic lattice.

The face is currently an infinite-volume H13 cone. The saved vector `Y` is orthogonal to the apex `V`, has square 145, and pairs to `-1` with each of the thirteen face normals. Hence the path `V+tY` reaches the future light cone with all thirteen inequalities strict. The 67 future-timelike inequalities are also strict there. A neighborhood of that ideal point is therefore uncovered, proving infinite volume. The review rechecks the apex against all 87 walls; it does not merely infer actuality from a positive-definite subdiagram.

For a new acute-compatible normal `r` of square `q>0`, write `c_i=-(r,s_i)>=0` for the seven E7 pairings. Its projected square is

`q_F = q - c^T B^{-1} c`.

An actual hyperbolic facet of the E7 face requires `q_F>0`. Every nonzero Coxeter pairing with a norm-two E7 normal has `c_i^2 >= q/2`; this covers all integer finite labels and also parallel or dotted pairs. Positivity of the inverse implies that any contact at an index with inverse diagonal at least two already costs at least `q`. The only possible contact is therefore wall 33, whose inverse diagonal is `3/2`. There, label three costs `3q/4` and leaves `q_F=q/4`; every label at least four costs at least `3q/2` and is impossible. Thus the complete necessary list is: all seven right, or exactly one π/3 contact at wall 33 and six right contacts. This is an all-label argument, not a label cutoff.

If `omega_33` is the checked E7 fundamental weight, then `omega_33^2=3/2`. An intrinsic normal `y` with square `q_F>0` has the two necessary lifting forms

`r=y`, or `r=y-sqrt(2q_F) omega_33`.

The second has square `4q_F` and normalized cosine `1/2` at wall 33. These formulas do **not** guarantee compatibility with the other ambient walls; every proposed lift still needs those checks.

The face cannot disappear in a nonempty full-dimensional acute-compatible extension retaining its seven walls. Project that extension's own strict timelike point to the E7 complement. Because `B^{-1}` is positive, the projected point is the strict point plus a positive combination of the seven outward normals, so every other acute inequality stays strict. A finite-volume finite-facet extension would consequently need a finite-volume completion of this H13 face. This is a necessary lower-dimensional construction problem, not a proof of global noncompletion. The review does not repeat the complete 87-wall pair classifications or ambient vertex/edge census.
