# SummedActualFootKernel: source-review successor

21 October 2026. **Source delta PASS.** No Lean or other proof job ran in
this review; no Lean source was edited.

The preserved gzip
`formal/source_archives/SummedActualFootKernel_before_sum_linearity_elaboration_fix_20261021.lean.gz`
has SHA256 `1675908400e23d99f8768a5859bdc83d40a736e87d0db8df0a2eb469695fd6e3`.
I independently decompressed and hashed it: its14247 bytes are exactly
the previously reviewed source
`f4dcb86bc2fe13d310fa22d030ae6b5f4f99cbb30376e0127cb7873601e3d3ec`.

The current `formal/CoxeterProofs/SummedActualFootKernel.lean` has14269
bytes, SHA256
`41d183ca03e0e405ce3d86ee838bdc3859c03f23ca19bb2cf9cd5239da050c28`.
A complete byte-derived text diff has exactly one changed line, at293:

```lean
-    exact entries
+    simpa only [sumN_add] using entries
```

This applies finite-sum addition normalization when matching the proved
entry identity to its goal. It changes no definition, theorem statement,
hypothesis, coefficient, finite range, conclusion, import or one of the
ten axiom-print requests. The previous mathematical source PASS in
`results/review_summed_actual_foot_kernel_source_20261021.md`, SHA256
`504c12b717cc69cae69624fb144ec13f129d69beecb0556d876b1f87c8fac629`,
therefore extends to this successor unchanged.

The earlier compilation failure and the actual successor compilation/
axiom-output receipts are separate execution evidence. This note does
not independently audit those executions. The kernel remains conditional
finite-coordinate algebra, not a formal whole-A20 or dimension theorem.
