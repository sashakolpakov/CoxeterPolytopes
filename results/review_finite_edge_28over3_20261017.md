# Independent internal review of the sharp finite-edge coefficient

The [six-simplex proof and coefficient 28/3](fresh_finite_edge_cusp_tradeoff_28over3_20261017.md)
pass a separate mathematical/source review and the independent exact checks
described below. No gap is identified in this lemma. This is a complete
internally checked incidence theorem, not an ambient dimension bound,
external refereeing or a proof-assistant formalization.

The mathematical review checked the minimum-order choice of the indefinite
subdiagram and why it makes the direct-white profile list exhaustive. It
also checked the treatment of arbitrarily many indirect whites through a
positive-definite Stieltjes block, the entrywise PSD relaxation and proved
multiplicity bounds. Thus the finite profile enumeration is not a hidden
cutoff on the supporting white rank.

For the remaining models, the rank-one singularity argument, affine linear
determinant implications, nonnegative-kernel correction argument and grounded
threshold calculations are sound. Grounding at a nonzero kernel coordinate
is justified by subtracting a multiple of that kernel from an arbitrary
vector; the rank-one correction preserves the kernel because its coordinate
is orthogonal to it. This gives the full scalar PSD interval, not just a test
at a sampled correction.

The [independent checker](../scripts/review_finite_edge_28over3_20261017.py)
uses Python `Fraction` and elementary Gaussian elimination. It imports
neither the producer's arithmetic engine nor SymPy. From the frozen profile
cost matrices it reconstructs and checks:

- All 111 single-profile maximum/tie envelopes.
- All 651 five-ideal determinant systems: 628 inconsistent systems, 21
  forcing at least six ideal vertices, and two forcing a positive offdiagonal
  entry of K, namely 1/4.
- All three remaining triple-profile systems, including their PSD residuals,
  nonnegative kernels, grounded thresholds and limiting corrected matrices.
- The rational rank-ten sharpness Gram, its inverse, a strict timelike center,
  seven finite vertices and three proper affine cusp links.
- Every saved coefficient case and the attained value `2(35+63)/21 = 28/3`.

The sharpness geometry was reviewed separately: the ten vertex rays are
future causal and generate the whole simplicial cone, while the three
connected affine cusp links are bounded. It is the classical finite-volume
H9 Coxeter simplex, not a new construction. The inherited finite-ended
simplex-block theorem then supplies the global summation of triangle costs.

Generation and a fresh saved replay both passed. The
[exact review certificate](review_finite_edge_28over3_20261017.json)
binds both frozen producer packages and records the 651 independent outcomes.
Syntax, direct hashes and local links passed. The earlier 539-case catalogue
and affine-profile isomorphism enumeration are reviewed as inputs; this
checker does not independently repeat that entire census. The arithmetic
and mathematical scopes are deliberately distinct. **FINAL/FROZEN.**

```sh
python3 scripts/review_finite_edge_28over3_20261017.py --verify
```
