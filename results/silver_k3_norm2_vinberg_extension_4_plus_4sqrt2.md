# 142 new silver walls at `x=4+4sqrt(2)`

The intrinsic norm-two Vinberg search is exhausted exactly through
`x=4+4sqrt(2)`.  This shell adds 142 walls to the certified 1,076-wall
prefix, producing an acute 1,218-wall prefix.  It is not a complete chamber,
and the shell certificate alone makes no compactness or finite-volume claim.

The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_extension_4_plus_4sqrt2.json`, with
canonical digest

```text
2caccab1462719095839350612a2278f52e01d66ac24ea25cd0c4748f5d10f8f
```

The independent block replay digest is

```text
34cd4e858e2bc4181b06bbf20c1d445ec2b08fb83c0a4d0d39585712147ccc42
```

## Exact shell census

For `s=sqrt(2)`, the finite-weight target is

\[
2+(1+s)(4+4s)^2=114+80s.
\]

The alphabet has 21 elements, with identity/conjugate square-root uppers
`(22,2)` and coefficient bound 24.  The production Schur transcript is

```text
visited: 1,21,441,1533,3192,4746,35049,122157,95718,36141
pruned:  0,0,368,1381,2966,3077,29232,117599,93997,27184
```

There are 957 exact weights.  Of these, 491 reconstruct nonintegrally; among
the 466 integral roots, 324 are blocked and 142 survive.  Blocker-set sizes
are `{1:255,2:69}`.

The independent replay exhausts raw component boxes
`(194481,194481,21)`, retains `(226,61,21)` bounded vectors, and uses 4,746
norm-join probes.  It reproduces the complete production accepted-root and
rejection-witness sets.  The frozen set digests are

```text
all reconstructed roots: 270095d6680c128cafcee53d25a025e71d319a31f9e71cc2a8e23284960d3652
rejection witnesses:     0261491d7dba61e2595b23b8daec9a0734ed7b71df0e4f2da9519254d3cea860
accepted-root set:        a172ba63ff68152d5fde252b6c128d3ad41705155b52a3c017978716a7b692d3
production root order:   22e9d1f86340ba1e4dba8995b0cfa8eeabf0ca5246eeca88481a6b6d8318ee80
```

All 10,011 distinct same-height pairings are strictly negative, so the
entire cohort enters simultaneously.  The full prefix inertias are
`(9,1,1208)` and `(10,0,1208)`.

The targeted `A2+F4+B2` edge retains finite endpoints `(33,622)` and gains no
new endpoint: all 142 new extensions have negative determinant, identity
inertia `(8,1,0)`, and conjugate inertia `(9,0,0)`.  This local statement is
not a global stopping condition.

The next admissible coordinate at this checkpoint was `x=6+3sqrt(2)`.  The
separate production-order
[1,218-wall endpoint/parabolic diagnostic](silver_k3_norm2_vinberg_1218_global_diagnostic.md)
has since been completed, and the next
[73-wall shell](silver_k3_norm2_vinberg_extension_6_plus_3sqrt2.md) produces a
1,291-wall prefix whose own
[global diagnostic](silver_k3_norm2_vinberg_1291_global_diagnostic.md) also
fails both stopping tests.  These remain finite-prefix failures, not an
eventual infinite-volume theorem; the status is
`undecidable_incomplete_simple_roots`.
