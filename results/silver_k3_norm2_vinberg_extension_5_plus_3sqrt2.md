# 350 new silver walls at `x=5+3sqrt(2)`

The intrinsic norm-two Vinberg search on the designed silver K3 carrier is
exhausted exactly through

\[
x=5+3\sqrt2.
\]

This shell adds 350 walls to the certified 726-wall prefix, producing an
acute 1,076-wall prefix.  This is not a complete chamber: no global Vinberg
endpoint/parabolic stopping gate is part of this shell certificate, and no
compactness or finite-volume conclusion follows from a bounded prefix.

The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_extension_5_plus_3sqrt2.json`; its
canonical digest is

```text
1de8fe6d72f76951410950c34ce6b7aee223f12dcdb3e78ca8f481c17c12e189
```

The independent block-factorized replay digest is

```text
1a6f3890736e6bb389ad4f55f87ad21d69ae5616c04d6db2bbcb31d3e6b70d87
```

## Exact production enumeration

Write `s=sqrt(2)`.  The finite-weight equation is

\[
z^T K^{-1}z=2+(1+s)(5+3s)^2=105+73s.
\]

The exact two-embedding bounds give identity and conjugate square-root
uppers 21 and 2, coefficient bound 23, and a 27-element alphabet.  The
canonical left-to-right Schur recursion gives

```text
visited by depth:       1,27,729,3456,8775,15120,146556,639441,553014,191916
Schur-pruned by depth:  0,0,601,3131,8215,9692,122873,618959,545906,142182
exact weights:          861
nonintegral roots:      0
blocked integral roots: 511
accepted roots:         350
```

Of the 511 blocked roots, 323 have one complete lower-wall blocker and 188
have two.  The certificate freezes both the full blocker-set digest and a
second digest including every positive blocker pairing.

## Independent exhaustive replay

The replay exhausts the orthogonal `B4`, `F4`, and `A1` component boxes
directly.  Its raw boxes have sizes `531441,531441,27`; exact bounds retain
`560,126,27` vectors, and the norm join uses 15,120 outer probes.  It obtains
the same 861 integral roots, the same 511 complete rejection sets, and the
same 350 accepted-root set.  The exact set digests are

```text
all reconstructed roots: 2cdc2b15a5beabaf2f1b06e2db61e8b2f6c60181004f31035fdd585a7530aaf4
rejection witnesses:     d0b303829d3b6deaa4b159697b18ebf8a14dfb4ff6f6b72c54e1f877a04889a4
accepted roots:           9b9afbf1c5b6e48edcf56fe06c9351d2e6476db29c155af6923ad4bb77469c94
```

The production ordering of the accepted cohort has digest
`56d394220b332d5c50b591f460ab348dadb7aee89cd3389191c87a87444a8bab`.
The manifest archives all 1,076 prefix roots and all 350 shell roots
literally in production order for downstream diagnostics.

## Acuteness, inertia, and the targeted edge

Every accepted root pairs nonpositively with every prior wall.  The 61,075
distinct same-height pairings split into 60,759 negative and 316 zero values,
so all 350 roots enter simultaneously.  The 1,076-wall Gram inertias are

```text
identity embedding:  (9,1,1066)
conjugate embedding: (10,0,1066)
```

The previously targeted `A2+F4+B2` edge retains finite endpoints `(33,622)`.
All 350 new edge extensions have negative determinant, identity inertia
`(8,1,0)`, and conjugate inertia `(9,0,0)`; this shell adds no endpoint to
that already locally closed edge.  This one-edge statement is not a global
Vinberg condition.

The next admissible control coordinate is `x=4+4sqrt(2)`.  Until exhaustive
simple-root search and a separate global endpoint/parabolic gate say
otherwise, the stopping status remains
`undecidable_incomplete_simple_roots`.
