# Review of the changed compact-face H9 seed

**PASS.** I source-reviewed the [construction](../scripts/fresh_h9_compactification_20261013.py), read its [scope and proof](fresh_h9_compactification_20261013.md), and independently invoked its saved-result replay:

```sh
python scripts/fresh_h9_compactification_20261013.py --verify
```

It passed with 810 proposals: 60 nonpositive scales, 438 forbidden old angles, and 312 individually compatible caps. The [review record](fresh_h9_compactification_20261013_review.json) binds the frozen source, output, and imported scalar implementation. This is mathematical/source review and a reused-checker replay, not a second implementation of all root arithmetic.

The compact H5 assertion has a global boundedness certificate. A strictly positive combination of all seven normals gives a timelike dual center, so its affine slice is bounded. Enumerating every five-wall support then gives exactly ten timelike extreme rays. The explicit ten PD9 lifts and strict facet projections justify the stated H9 seed geometry.

The field argument is valid. For `K=Q(sqrt(5))` and `L=K(sqrt(R))`, normality would require `Norm(R)=16/11` to be a square in K. A positive rational number is a square in K only when it or its quotient by five is a rational square; neither condition holds. Thus L is a nonnormal quartic with unique quadratic subfield K. Its abelian subfields are therefore contained in K, giving exactly the finite squared-angle labels 2, 3, 4, 5, 6, and 10 used for these profiles. This bound does not apply to arbitrary caps in larger fields.

The end-wall Cartan cycle introduces the quadratic extension over K. The conjugation negating its square root fixes the indefinite five-wall block, so the stated Vinberg semidefiniteness obstruction to quasi-arithmeticity is justified. This conclusion is conditional on retaining these facets in any proposed finite-volume completion.

Finally, the common strict spacelike escape is checked against every accepted cap by exact signs and squared comparisons. All 312 caps are strict at `(h,0)` because their old components are positive multiples of future timelike vertices. A sufficiently small inward A4 perturbation gives a common strict timelike point for this finite bank. The escape therefore proves that the whole bank, and every subset, has infinite volume. This closes only the specified vertex-centered family; it does not exclude completing the changed seed by other normals.
