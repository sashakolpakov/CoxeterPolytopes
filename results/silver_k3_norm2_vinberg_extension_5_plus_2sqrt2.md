# Exact empty silver shell at `x=5+2sqrt(2)`

The intrinsic norm-two Vinberg search on the designed silver K3 carrier is
now exhausted through

\[
x=5+2\sqrt2.
\]

This shell has one integral norm-equation solution.  Its reconstructed root
is blocked by a lower-height wall, so the shell contributes no simple wall
and the certified acute prefix remains at 333 walls.  The targeted
`A2+F4+B2` edge still has only wall 33 as a certified finite endpoint.

The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_extension_5_plus_2sqrt2.json`; its
canonical digest is

```text
e2283e6330c93632035febe1b84ef1d499a7b9507931008940b956f1470e10d5
```

The block-factorized replay has digest

```text
b4b101544f6eb548ab43cc649275665a87f42a61ab88987ffdcea9638ce3b76a
```

## Exact shell enumeration

Write `s=sqrt(2)`.  The finite-weight equation is

\[
z^T K^{-1}z=2+(1+s)(5+2s)^2=75+53s.
\]

The exact bounds at both real embeddings give identity and conjugate square-
root uppers 18 and 1, coefficient bound 19, and the four-element alphabet

```text
0, 3+2s, 4+3s, 7+5s.
```

The production left-to-right recursion gives

```text
visited by depth:         1,4,16,16,16,12,32,44,24,12
Schur-pruned by depth:    0,0,12,12,13,4,21,38,21,4
exact weight solutions:  1
nonintegral roots:       0
integral roots:          1
rejected by old walls:   1
accepted at this shell:  0
```

The rejected root has exactly one positive lower-wall pairing, to wall 9.  Production
records the complete blocker tuple; its witness digest is reproduced exactly
by the independent search.

## Independent block replay

The hostile replay independently derives the coefficient box, exhausts the
orthogonal `B4`, `F4`, and `A1` component boxes, and performs an exact field-
norm join.  It shares neither the production coordinate order nor its
leading-Schur pruning tree.  The raw boxes have sizes `256,256,4`; exact
bounds retain `3,1,4` vectors, and the keyed join makes 12 outer probes.

It obtains the same unique integral root and the same complete one-wall
blocker tuple.  The reconstructed-root and blocker-witness digests are

```text
997eb239f9c0498f8c8cc2483c2dc6a164e184b22f69895cf0e6f92777e5d781
8330422302f1de59f9d6bebd43171577cd69a69793414753c9203d5ef42e4d94.
```

The replay also reconstructs the unchanged prefix inertias, endpoint census,
and next admissible coordinate without using the production Gram audit.

## Empty cohort and unchanged prefix

There are no accepted same-height roots, hence no same-height pairings and no
new-to-old acute inequalities to add.  The 333-wall prefix remains acute with
exact inertias

```text
identity embedding:   (9,1,323)
conjugate embedding:  (10,0,323).
```

The principal Gram matrix on walls 0 through 9 still has inertias `(9,1,0)`
and `(10,0,0)`, so it spans the rank-ten carrier; factorization through that
carrier forces the other 323 directions to be null.

## The targeted edge remains open

The diagnostic elliptic edge remains

```text
E = (1,2,4,5,6,7,9,12),   type A2+F4+B2.
```

Because this shell contributes no wall, it contributes no finite endpoint.
The independently recomputed finite endpoint set remains exactly wall 33, and
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
x=4+3\sqrt2.
\]

The formal status remains

```text
undecidable_incomplete_simple_roots.
```

“Empty shell” is a bounded norm-two statement, not a stopping theorem.  There
is no proof here of wall completeness, no global audit of every elliptic edge
and affine neighbourhood, and therefore no conclusion of compactness, finite
volume, or infinite volume.  Primitive reflective vectors of norms other
than two also remain outside this enumeration.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_silver_k3_vinberg_5_plus_2sqrt2.py
PYTHONPATH=src:. python3 scripts/independent_silver_k3_vinberg_5_plus_2sqrt2.py
PYTHONPATH=src:. pytest -q tests/test_silver_k3_vinberg_5_plus_2sqrt2.py
```
