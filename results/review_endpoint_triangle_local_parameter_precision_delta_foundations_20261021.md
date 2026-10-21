# Endpoint triangle: final local-parameter precision delta

21 October2026. **Static exact-delta PASS.** No execution, compilation,
test, import or edit to a bound source/map/evidence file was performed.

Current source `formal/CoxeterProofs/EndpointTriangleNormKernel.lean` is
SHA256 `ce2ad70b636880d68e0e99f345f82b8db3c7bd7197636ea9c93c0d23d98df1a4`.
The complete diff against the preserved predecessor is exactly ONE line:

```diff
-  grind only [scaled]
+  grind only
```

The explicit local `scaled` equality immediately before it remains in
scope. All ten theorem statements, assumptions, the universal cancellation,
rewrites, other proof bodies and ten requested printers are byte-unchanged.
Removing its redundant explicit tactic argument changes no mathematical
input or derived identity. The exact ring reduction in the previous
addendum still applies: the rewritten goal difference is
`4*h*(17*A²-8*U-17*b²+17)`, supplied as zero by that local equality.

This binds the original complete source/map review
`results/review_endpoint_triangle_norm_kernel_source_foundations_20261021.md`
(SHA `d8318de0c4b89cf40685b073cf1f9fef0e92e700c86b1e0dce78914fa4bc63af`)
and explicit-cancellation addendum
`results/review_endpoint_triangle_explicit_cancellation_delta_foundations_20261021.md`
(SHA `d92a00ae36c3852bde3d299a75cda710faff652e3b5f2985476c6a27ae4be290`)
to the current source without reopening unchanged proofs. Application
geometry, number-field norm interpretation and the separately reviewed
ordinary-A1 extension retain exactly the same external boundary.

The predecessor is preserved at
`results/source_archive/endpoint_triangle_norm_before_redundant_parameter_fix_20261021.lean.txt`,
SHA `92c9467a3b0ee96d2f41633235a48ab7d0325f30e8c76e8a35d8cbbb90231d16`.
I read the full second failure: it reports only the redundant `scaled`
parameter error, the retained unused-section-variable warning, and recovery
axiom output. Preservation record
`results/endpoint_triangle_norm_redundant_parameter_failure_preservation_20261021.json`
has SHA `5d2781b1eff51753b67ebf4183201f02a120672c62d413acc216a0ae4a260f64`;
its failure text has SHA
`452a2c1dd8ec43e528c4eb36e65196bf999bcb98ba4c34a49d5f1f476ac0c72d`.
Those bindings match the current files. Both failed batches and their
recovery `sorryAx` prints receive ZERO verification credit.

No correction requested. This review approves the unchanged logic at the
final source identity; it does not claim a successful retry or actual
separate axiom audit. Root owns that future evidence.
