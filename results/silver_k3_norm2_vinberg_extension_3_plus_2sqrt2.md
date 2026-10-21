# Sixty new silver walls at `x=3+2sqrt(2)`

The intrinsic norm-two Vinberg search on the designed silver K3 carrier has
now been exhausted through

\[
x=3+2\sqrt2.
\]

This shell adds exactly 60 walls relative to the frozen 93-wall prefix, giving
a 153-wall acute prefix.  None of the 60 walls supplies the missing second
finite endpoint of the targeted `A2+F4+B2` edge.

The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_extension_3_plus_2sqrt2.json`; its
canonical digest is

```text
b10548b35ce263769c858b54a572b02c6ab148f508136ce14946f1f6d9016372
```

The block-factorized replay has digest

```text
f45bcd9e02db07db2fb0c7a0bb9fe0363b7490803d45da5bc697e8fc8ae65600
```

## Exact shell enumeration

Write `s=sqrt(2)`.  The finite-weight equation is

\[
z^T K^{-1}z=2+(1+s)(3+2s)^2=43+29s.
\]

The exact bounds at both real embeddings give the 18-element alphabet

```text
0, 1, s, 1+s, 2+s, 1+2s, 3+s, 2+2s, 3+2s,
4+2s, 3+3s, 4+3s, 5+3s, 4+4s, 6+3s, 5+4s,
6+4s, 7+4s.
```

The production left-to-right recursion gives

```text
visited by depth:         1,18,324,1080,2016,2772,20142,66078,47934,16290
Schur-pruned by depth:    0,0,264,968,1862,1653,16471,63415,47029,11339
exact weight solutions:  138
nonintegral roots:       0
integral roots:          138
rejected by old walls:   78
accepted at this shell:  60
```

Among the 78 rejected roots, 45 have one explicit lower-wall blocker and 33
have two.  The frozen manifest records the blocker-wall histogram and a
digest of every root, blocker index, and positive pairing.

The hostile replay independently rebuilds the alphabet, exhausts the
`B4`, `F4`, and `A1` coordinate boxes directly, and performs an exact
field-norm join.  It has neither the production coordinate order nor its
leading-Schur pruning tree.  The raw component boxes have sizes
`104976,104976,18`; exact two-embedding norm bounds retain `154,38,18`
vectors, and the keyed join makes 2,772 outer probes.  It recovers the same
138 weights.  The reconstructed-root and rejection-witness digests are frozen
in the manifest.

## Simultaneous admission

Every accepted root pairs nonpositively with all 93 lower-height walls.  Of
the `binomial(60,2)=1770` distinct same-height pairings, exactly 1,702 are
negative and 68 are zero.  Thus no root in this cohort cuts off another, and
all 60 enter independently of candidate order.  The accepted-root and full
same-height-pairing digests are frozen in the manifest.

Here “simple” means relative to the exhausted lower-height norm-two prefix and
the complete same-height cohort.  It does not mean that the full chamber wall
list has terminated.

The 153-wall prefix Gram matrix is acute.  Its exact inertias are

```text
identity embedding:   (9,1,143)
conjugate embedding:  (10,0,143).
```

These inertias use a compact exact rank witness: the principal Gram matrix on
walls 0 through 9 has inertias `(9,1,0)` and `(10,0,0)`.  It therefore spans
the rank-ten carrier, while factorization of the full wall Gram through that
carrier forces the other 143 directions to be null.

## The targeted edge remains open

The diagnostic elliptic edge is still

```text
E = (1,2,4,5,6,7,9,12),   type A2+F4+B2.
```

For every new wall `r93,...,r152`, the principal Gram matrix on
`E union {r_i}` has identity inertia `(8,1,0)` and conjugate inertia
`(9,0,0)`.  Hence none is an elliptic rank-nine endpoint.  The exact
determinants and their multiplicities are frozen in the manifest.

Wall 33 therefore remains the only certified finite endpoint, and the second
endpoint status remains

```text
not_found_through_this_shell.
```

An ideal endpoint cannot be decided from this bounded prefix.  Vinberg's
criterion requires the full componentwise-affine zero-neighbourhood in the
complete simple-wall system.

## Literal status

The next admissible control coordinate is

\[
x=4+2\sqrt2.
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
PYTHONPATH=src:. python3 scripts/certify_silver_k3_vinberg_3_plus_2sqrt2.py
PYTHONPATH=src:. python3 scripts/independent_silver_k3_vinberg_3_plus_2sqrt2.py
PYTHONPATH=src:. pytest -q tests/test_silver_k3_vinberg_3_plus_2sqrt2.py
```
