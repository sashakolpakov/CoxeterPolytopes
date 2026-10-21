# Independent review of the sixteen star metric orientations

**PASS.** The [primary family theorem](fresh_h9_star_cube_metric_20261014_orientations.md) excludes compact H9 Coxeter realizations of star case 0 with the four fixed B2 blocks and, in each original positive-leaf/cut-leaf pair, at least one leaf having own contacts `(3,3)`. The opposite leaf has arbitrary allowable Coxeter contacts. This does not exclude the full star incidence, changed B2 labels, or a block in which neither leaf has own contacts `(3,3)`.

The [independent checker](../scripts/fresh_h9_star_cube_metric_20261014_orientations_review.py) and [certificate](fresh_h9_star_cube_metric_20261014_orientations_review.json) import no producer or SymPy arithmetic. They reuse the frozen independent Fraction field engine in the basis `1,sqrt(2),sqrt(5),sqrt(10)` and its exact rational-enclosure sign tests. Every distinct residual determinant is recomputed; the review does not merely trust cached primary values or reuse the preceding determinant table.

The domain checks pass in their necessary logical order:

- The **24 cross-B2 right-contact witnesses** use only the four B2 bonds, with every leaf contact relaxed to right. The review also reconstructs each actual four-root `4–3–4` path and its positive null vector. Thus this preliminary gate is independent of the leaf orientation.
- All **16 eight-root budget supports** are checked against the literal vertex incidence, covering both retained-own-root choices of all eight leaves. Their inverses are reconstructed exactly. A `(3,3)` leaf permits at most one foreign contact, necessarily label 3; every other admissible profile permits none.
- The own positive/cut pair is forced right before imposing a label bound. A free leaf's missing triangle has a nonright own contact, while its partner's contact is 3. Testing a nonright mutual contact gives a nonpositive all-ones quadratic on an actual three-root face. The review checks **704 orientation/profile instances** and the generic inequality. The eight actual `3–6` paths then exclude every own label at least 6. The Lanner inequality leaves exactly the eleven ordered profiles, including `(2,5)` and `(5,2)`.
- All **2,112 upper-contact witnesses** are freshly reconstructed for their literal orientations and profiles. Ordinary exact row elimination verifies the saved nonpositive leading minors and their pivot identities. Unused free-leaf contacts are relaxed to right. Restoring any omitted acute contacts cannot repair positive definiteness, by monotonicity of the largest eigenvalue of the nonnegative matrix in the normalized Gram expression `I-A`.

No positive-leaf/cut-leaf swap is assumed to be a polytope automorphism. The rank calculation keeps the literal three normals `P3=11`, `C3=17`, and `C0=14`. When the orientation fixes `C0`, its profile is fixed to `(3,3)` rather than replacing it by a different facet. Its two tested anchor pairs are actual meeting pairs, with freshly verified bounds to labels 2 or 3.

The separate census recomputes **12,852 distinct residual determinants**, comparing them with every one of the **91,392 saved occurrences**. The 189 distinct two-anchor residual matrices are nondegenerate. There are **91,328 nonzero three-dimensional residual determinants**; together with the PD rank-eight block, each would force rank at least eleven in the ten-dimensional ambient vector space.

Exactly **64 occurrences** have zero residual determinant. Every one has an independently reconstructed nonpositive principal minor inside an actual vertex, so none is accepted from rank alone. In 32 occurrences the saved vertex also contains an explicit four-root affine cycle with a positive all-ones null vector; the remaining occurrences are closed by their full actual-vertex failure minors. The review makes no unnecessary diagram-type assertion for those remaining cases.

Generation and saved-result replay pass. The certificate binds the final primary source/data and recursively checks their source manifests. The older first-state, monotone-cut, and complete-cut reviews remain frozen. This review does not audit the separate continuation state or its subsequent single-assignment metric gate.

```sh
python3 scripts/fresh_h9_star_cube_metric_20261014_orientations_review.py --verify
```
