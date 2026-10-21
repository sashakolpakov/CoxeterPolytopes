# 73 new silver walls at `x=6+3sqrt(2)`

The intrinsic norm-two Vinberg search is exhausted exactly through
`x=6+3sqrt(2)`.  This shell adds 73 walls to the certified 1,218-wall
prefix, producing an acute 1,291-wall prefix.  It is not a complete chamber,
and the shell certificate alone makes no compactness or finite-volume claim.

The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_extension_6_plus_3sqrt2.json`, with
canonical digest

```text
01c8c722b0db00f757d67789d82674f2e9724ec2e81e4b7e720870e1700a67f5
```

The independent block replay digest is

```text
da8031eeceb64650441027f81afd7321f6c910fdbc81e00965fe4d2b171a6295
```

## Exact shell census

For `s=sqrt(2)`, the finite-weight target is

\[
2+(1+s)(6+3s)^2=128+90s.
\]

The alphabet has 20 elements, with identity/conjugate square-root uppers
`(23,2)` and coefficient bound 25.  The production Schur transcript is

```text
visited: 1,20,400,1380,2780,3920,28820,98760,76020,27780
pruned:  0,0,331,1241,2584,2479,23882,94959,74631,20465
```

There are 697 exact weights.  Of these, 425 reconstruct nonintegrally; among
the 272 integral roots, 199 are blocked and 73 survive.  Blocker-set sizes
are `{1:127,2:72}`.  The manifest stores all 199 rejected roots with their
complete blocker-index sets and all 73 accepted roots literally.

The independent replay exhausts raw component boxes `(160000,160000,20)`,
retains `(196,51,20)` bounded vectors, and uses 3,920 norm-join probes.  It
reproduces the complete production accepted-root and rejection-witness sets.
The frozen set digests are

```text
all reconstructed roots: 33bf391225bf203da4b3b4e15e1fddf1b61b663476d70c027df80b28930d8e74
rejection witnesses:     f0b09f02437f8a624c9fa125d0d4e888020c49bee2f6339114fe9e90bca95729
accepted-root set:        266f2683d83a9f70b6b464fa9d12e2fddc43863efb73533a73cceacbe7324333
production root order:   5ef1a4b9924f1f672101c4c2992252946d75d6f4e9a052593e8cd0ce18d9ea7d
```

All 2,628 distinct same-height pairings are strictly negative, so the entire
cohort enters simultaneously.  The full prefix inertias are `(9,1,1281)`
and `(10,0,1281)`.

The next admissible coordinate is `x=5+4sqrt(2)`.  The separate production-
order [1,291-wall endpoint/parabolic diagnostic](silver_k3_norm2_vinberg_1291_global_diagnostic.md)
is now complete.  It finds 43,891 elliptic rank-eight edge diagrams and no
connected affine diagram or admissible ideal vertex; 28,898 rank-eight edge
diagrams are endpoint-deficient.  Neither compact nor finite-volume stopping
fires, so no compact or finite-volume example is certified.  This is only a
finite-prefix failure, not an eventual infinite-volume theorem, and the
status remains `undecidable_incomplete_simple_roots`.
