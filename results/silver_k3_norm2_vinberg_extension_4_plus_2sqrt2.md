# One hundred twenty-one new silver walls at `x=4+2sqrt(2)`

The intrinsic norm-two Vinberg search on the designed silver K3 carrier is
now exhausted through

\[
x=4+2\sqrt2.
\]

This shell adds exactly 121 walls relative to the frozen 153-wall prefix,
giving a 274-wall acute prefix.  None of the 121 walls supplies the missing
second finite endpoint of the targeted `A2+F4+B2` edge.

The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_extension_4_plus_2sqrt2.json`; its
canonical digest is

```text
fb9d687f2e536cc00688d365db95935ee7effb70733d2c82c597c6712c465ff5
```

The block-factorized replay has digest

```text
fc05a1563fdd4e905850bb5ccc756bfb3c2d18e516cf5b0a85f13a3e83718726
```

## Exact shell enumeration

Write `s=sqrt(2)`.  The finite-weight equation is

\[
z^T K^{-1}z=2+(1+s)(4+2s)^2=58+40s.
\]

The exact bounds at both real embeddings give the 19-element alphabet

```text
0, 1, s, 1+s, 2+s, 3+s, 2+2s, 3+2s, 4+2s, 3+3s,
4+3s, 5+3s, 4+4s, 5+4s, 6+4s, 7+4s, 6+5s, 7+5s, 8+5s.
```

The production left-to-right recursion gives

```text
visited by depth:         1,19,361,1178,2261,3211,22363,73568,56012,21299
Schur-pruned by depth:    0,0,299,1059,2092,2034,18491,70620,54891,15783
exact weight solutions:  653
nonintegral roots:       330
integral roots:          323
rejected by old walls:   202
accepted at this shell:  121
```

Among the 202 rejected integral roots, 160 have one explicit lower-wall
blocker and 42 have two.  The manifest records every blocker index and
positive pairing through a frozen witness digest, as well as the complete
blocking-wall occurrence histogram.

## Independent block replay

The hostile replay independently derives the two-embedding coefficient box,
then exhausts the orthogonal `B4`, `F4`, and `A1` component boxes and performs
an exact field-norm join.  It shares neither the production coordinate order
nor its leading-Schur pruning tree.  The raw boxes have sizes
`130321,130321,19`; exact bounds retain `169,46,19` vectors, and the keyed join
makes 3,211 outer probes.  It obtains the same 653 weights, 330 nonintegral
reconstructions, 202 rejected integral roots with the same complete blocker
sets, and the same 121 accepted roots.

The replay also independently reconstructs the full prefix audit: all
accepted-to-old and accepted-to-accepted pairings are nonpositive, the first
ten old walls remain a rank-ten witness, the prefix inertias are
`(9,1,264)` and `(10,0,264)`, the targeted finite endpoint set is exactly
wall 33, and the next admissible coordinate is `3+3sqrt(2)`.

## Simultaneous admission and prefix inertia

Every accepted root pairs nonpositively with all 153 lower-height walls.  Of
the `binomial(121,2)=7260` distinct same-height pairings, exactly 7,180 are
negative and 80 are zero.  Thus no root in this cohort cuts off another, and
all 121 enter independently of candidate order.

Here “simple” means relative to the exhausted lower-height norm-two prefix and
the complete same-height cohort.  It does not mean that the full chamber wall
list has terminated.

The 274-wall prefix Gram matrix is acute and has exact inertias

```text
identity embedding:   (9,1,264)
conjugate embedding:  (10,0,264).
```

The principal Gram matrix on walls 0 through 9 has inertias `(9,1,0)` and
`(10,0,0)`, so it spans the rank-ten carrier; factorization through that
carrier forces the other 264 directions to be null.

## The targeted edge remains open

The diagnostic elliptic edge is still

```text
E = (1,2,4,5,6,7,9,12),   type A2+F4+B2.
```

For every new wall `r153,...,r273`, the principal Gram matrix on
`E union {r_i}` has identity inertia `(8,1,0)`, conjugate inertia `(9,0,0)`,
and negative identity-embedding determinant.  Hence none is an elliptic
rank-nine endpoint.  Wall 33 remains the only certified finite endpoint, and
the second endpoint status is

```text
not_found_through_this_shell.
```

An ideal endpoint cannot be decided from this bounded prefix.  Vinberg's
criterion requires the full componentwise-affine zero-neighbourhood in the
complete simple-wall system.

## Literal status

The next admissible control coordinate is

\[
x=3+3\sqrt2.
\]

The formal status remains

```text
undecidable_incomplete_simple_roots.
```

There is no proof here of wall completeness, no global audit of every
elliptic edge and affine neighbourhood, and therefore no conclusion of
compactness, finite volume, or infinite volume.  Primitive reflective vectors
of norms other than two also remain outside this enumeration.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_silver_k3_vinberg_4_plus_2sqrt2.py
PYTHONPATH=src:. python3 scripts/independent_silver_k3_vinberg_4_plus_2sqrt2.py
PYTHONPATH=src:. pytest -q tests/test_silver_k3_vinberg_4_plus_2sqrt2.py
```
