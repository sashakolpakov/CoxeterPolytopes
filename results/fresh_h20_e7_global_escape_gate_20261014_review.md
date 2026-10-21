# Independent arithmetic review of the E7-face escape gate

**2026-10-14 — generation and saved replay PASS.** The [separate checker](../scripts/fresh_h20_e7_global_escape_gate_20261014_review.py) and [certificate](fresh_h20_e7_global_escape_gate_20261014_review.json) reconstruct the [primary gate](fresh_h20_e7_global_escape_gate_20261014.md) from the physical frozen 87-wall normals, the escape and the E7 face data. No producer or project arithmetic is imported. Sparse multiplication in the eight square-root monomials and outward 192-bit rational enclosures replace the primary checker's recursive quadratic sign algorithm.

The reconstruction verifies the original escape's positive norm, future orientation and all 87 strict signs. It recomputes the thirteen face weights, the full rank-20 finite-corner Gram and inverse, every reduced old-wall coefficient, and the target coefficients a, B and w. Every stored coefficient in the primary dual, elimination and convexity certificates agrees exactly with the independent reconstruction.

For the right lift, all three rational multipliers are positive and all fourteen dual residual coefficients are strictly negative. For the nonright lift, the signs of the two old-wall timelike coefficients are checked before eliminating t. The combined inequality forces d12 below one; after d12=0 the wall-16 inequality forces d0 above one. Thus the two discrete Coxeter gaps are applied in the correct direction.

The independent rank-one positive-definiteness gate, positive denominator and eleven positive gradient numerators certify the stated constrained minimizer. The exact quantity `(4/3)D−2C²` is positive, so the minimum projected norm is strictly greater than 1/2. This contradicts the required norm of the unique nonright E7 lift under an escape-cutting assumption. The same argument excludes zero pairing with the escape: replace the strict upper bound on t by a weak one, and the final positive norm margin still yields a contradiction.

The automatic time-orientation argument is also sound. The full corner inverse is entrywise nonnegative, T is future with norm −2, and the twenty corner roots span its positive orthogonal complement. A valid proposed common strict interior point then forces t>0 for any added acute normal. Retaining the old finite apex is therefore a consequence; it is not an extra assumed constraint on candidate walls.

The proof covers every magnitude in the necessary Coxeter relaxation, hence all finite labels and the parallel/dotted cases. It concerns only walls inducing actual facets on this selected E7 face while the old 87 facets are retained. Walls outside that class and replacements of old facets remain available. No old full pair/corner census was repeated, and no unrestricted H20 completion obstruction is claimed.

```sh
python3 scripts/fresh_h20_e7_global_escape_gate_20261014_review.py --verify
```
