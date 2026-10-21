# Endpoint triangle: explicit-cancellation successor review

21 October2026. **Static exact-delta/logic PASS.** No source correction
requested. No compilation, test, worker import, mathematical execution or
file mutation beyond this new review was performed.

Reviewed current source:
`formal/CoxeterProofs/EndpointTriangleNormKernel.lean`, SHA256
`92c9467a3b0ee96d2f41633235a48ab7d0325f30e8c76e8a35d8cbbb90231d16`.
The original source is preserved byte-exact at
`results/source_archive/endpoint_triangle_norm_before_explicit_cancellation_20261021.lean.txt`,
SHA `e187298414441c7d65013928dd06e4c00c9cc72fb58474962c8c6dbbe4ef8a95`.
Its independent full source/map review remains unchanged at
`results/review_endpoint_triangle_norm_kernel_source_foundations_20261021.md`,
SHA `d8318de0c4b89cf40685b073cf1f9fef0e92e700c86b1e0dce78914fa4bc63af`.

The complete unified diff has ONE hunk: only the tactic body of
`cycle_expansion` changed. Every theorem statement, parameter, import,
other proof and all ten requested axiom-printer commands are unchanged.
The map is also unchanged at SHA
`c9beffd99c3d4b963573dad65fcfb28637b4d85de215c187a2e6b6886fef4a75`.

I independently checked the replacement, not merely its hash. The new
local `cancellation` is the universal commutative-ring identity

    4(Aa-br)²+8br(Aa-br)+4b²
      =4(A²a²-b²r²+b²).

Rewriting lambda in the supplied cycle identity, applying that identity,
then rewriting a²=h and r²=1+h leaves

    rho=4(A²h-b²(1+h)+b²).

Rewriting rho and v=4h in the goal makes its left side minus right side
equal to

    4h*(17A²-8U-17b²+17).

This is exactly zero by the supplied time equation multiplied by4h,
which the new `congrArg` supplies explicitly. Thus the final
`grind only [scaled]` has only that ring consequence to normalize;
it needs no new assumption, division, positivity premise or desired
cycle conclusion as an input. The replacement exposes the same
mathematical derivation checked in the original full source review.

The initial failure is honestly preserved in
`results/endpoint_triangle_norm_initial_elaboration_failure_20261021.txt`,
SHA `ffbdc29e8ef0fa7e76a80a5751c2a1fda69e997239e8d110673c5a0d92ef0184`.
The preservation record
`results/endpoint_triangle_norm_initial_failure_preservation_20261021.json`,
SHA `dfe7283a1576bb9af9ef3864cd05049d199edce3f4d483874636072e0ffe2f7c`,
binds that failure and exact predecessor. I read the failure: it reports
the original `grind` failure at `cycle_expansion` and an error-recovery
`sorryAx` output. It is NOT a verified declaration or a successful
separate axiom audit; no declaration receives credit from that batch.

This addendum transfers the unchanged full source/assumption review to
the repaired tactic source, with ordinary-A1 and all geometric/number-field
application boundaries still external. It does not claim elaboration has
now succeeded or predict an actual axiom list. A fresh guarded compile
and separate actual audit remain root-owned, pending execution evidence.
