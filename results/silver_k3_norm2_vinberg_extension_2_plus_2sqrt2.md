# Fifty-six new silver walls at `x=2+2sqrt(2)`

The intrinsic norm-two Vinberg search on the designed silver K3 carrier has
now been exhausted through

\[
x=2+2\sqrt2.
\]

This shell adds exactly 56 walls relative to the frozen 37-wall prefix, giving
a 93-wall acute prefix.  None of the 56 walls supplies the missing second
finite endpoint of the targeted `A2+F4+B2` edge.

The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_extension_2_plus_2sqrt2.json`; its
canonical digest is

```text
4cbc5d390cccf160cac882918c229900b3b5540331d12db82d1686ab2beb6e49
```

The reversed-coordinate replay has digest

```text
a31ecd05c1d816380ba20f78db78b352d0e617e9e9e2eae6749abc17d9a524f3
```

## Exact shell enumeration

Write `s=sqrt(2)`.  The finite-weight equation is

\[
z^T K^{-1}z=2+(1+s)(2+2s)^2=30+20s.
\]

The exact bounds at both real embeddings give the 16-element alphabet

```text
0, 1, s, 1+s, 2+s, 1+2s, 3+s, 2+2s,
3+2s, 4+2s, 3+3s, 4+3s, 5+3s, 4+4s, 6+3s, 5+4s.
```

The production left-to-right recursion gives

```text
visited by depth:         1,16,256,656,1072,1392,7984,22016,15840,6368
Schur-pruned by depth:    0,0,215,589,985,893,6608,21026,15442,4706
exact weight solutions:  277
nonintegral roots:       129
integral roots:          148
rejected by old walls:   92
accepted at this shell:  56
```

Among the 92 rejected integral roots, 77 have one explicit lower-wall blocker
and 15 have two.  The frozen manifest records the blocker-wall histogram and
a digest of every root, blocker index, and positive pairing.

The hostile replay fixes the nine weight coordinates in the opposite order
and rebuilds its own leading Schur complements.  It obtains the same 277
weight solutions, 129 nonintegral reconstructions, 92 rejected integral roots,
and 56 accepted roots.  Its distinct traversal is

```text
visited by depth:         1,16,256,2336,8816,6832,2592,14160,37264,67264
Schur-pruned by depth:    0,0,110,1785,8389,6670,1707,11831,33060,65602
```

The reconstructed-root digest is
`2927acfebaef149f8df0aff541d3090119ba1de4efab093aaf22f2a0734b63c3`,
and the independently sorted rejection-witness digest is
`cbb2a20806b09c9d09ae3f627c4141e7ed7eecd6564e91ce06d9cd7d9fbe9028`.

## Simultaneous admission

Every accepted root pairs nonpositively with all 37 lower-height walls.  Of
the `binomial(56,2)=1540` distinct same-height pairings, exactly 1,482 are
negative and 58 are zero.  Thus no root in this cohort cuts off another, and
all 56 enter independently of candidate order.  The accepted-root digest is

```text
f87130615c11d7f550f5901875144bff46f2c0df83b7dd614f24ad9add4851aa
```

Here “simple” means relative to the exhausted lower-height norm-two prefix and
the complete same-height cohort.  It does not mean that the full chamber wall
list has terminated.

The complete 93-wall Gram matrix is acute, has digest

```text
8b034920fab8305dd34240968a62b4d6ea2178584585bf2ebf5eef1cb0f054cd
```

and exact inertias

```text
identity embedding:   (9,1,83)
conjugate embedding:  (10,0,83).
```

## The targeted edge remains open

The diagnostic elliptic edge is still

```text
E = (1,2,4,5,6,7,9,12),   type A2+F4+B2.
```

For every new wall `r37,...,r92`, the principal Gram matrix on
`E union {r_i}` has identity inertia `(8,1,0)` and conjugate inertia
`(9,0,0)`.  Hence none is an elliptic rank-nine endpoint.  The exact
identity-embedding determinants, grouped with multiplicity, are

```text
3*(-2172-1536s)   1*(-1632-1156s)   2*(-1392-986s)
2*(-1380-978s)    1*(-1080-766s)    3*(-1060-750s)
1*(-1000-708s)    1*(-888-632s)     2*(-876-620s)
4*(-816-578s)     4*(-804-570s)     3*(-444-316s)
1*(-420-302s)     3*(-420-300s)     1*(-408-294s)
3*(-324-230s)     3*(-312-222s)     4*(-180-130s)
2*(-160-114s)     3*(-108-80s)      2*(-100-72s)
2*(-96-74s)       1*(-96-68s)       2*(-84-66s)
1*(-60-48s)       1*(-4-6s).
```

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
x=3+2\sqrt2.
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
PYTHONPATH=src:. python3 scripts/certify_silver_k3_vinberg_2_plus_2sqrt2.py
PYTHONPATH=src:. python3 scripts/independent_silver_k3_vinberg_2_plus_2sqrt2.py
PYTHONPATH=src:. pytest -q tests/test_silver_k3_vinberg_2_plus_2sqrt2.py
```
