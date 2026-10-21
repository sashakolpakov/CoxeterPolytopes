# The silver shell `x=1+2sqrt(2)` is exactly empty

This certificate continues the intrinsic norm-two silver K3 Vinberg search
through the first admissible control coordinate after `x=2+sqrt(2)`.  It is a
bounded negative result: the shell contributes no new simple wall, so it does
not supply the missing second endpoint of the targeted rank-eight edge.

The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_extension_1_plus_2sqrt2.json`; its
canonical digest is

```text
e644e073dacd4c691821409920c723f8c2edebc2ee7c9d751cb4b2059d0950c2
```

## Exact coordinate gate

For a norm-two root with control coordinate `x`, the positive conjugate
equation gives

\[
(x')^2\leq 2(1+\sqrt2).
\]

Simultaneously imposing the identity bound through `1+2sqrt(2)` leaves
exactly

\[
1,\quad \sqrt2,\quad 2,\quad 1+\sqrt2,\quad
2+\sqrt2,\quad 1+2\sqrt2.
\]

The first five levels were exhausted by the previous certificates.  The next
admissible value after the present shell is `3+sqrt(2)`.

## Exact shell enumeration

At

\[
x=1+2\sqrt2
\]

the nine finite-simple-root weights satisfy

\[
z^T K^{-1}z=2+(1+\sqrt2)x^2=19+13\sqrt2.
\]

The two embeddings leave the seven-element alphabet

```text
0, 1, 1+sqrt(2), 2+sqrt(2), 2+2sqrt(2),
3+2sqrt(2), 4+3sqrt(2).
```

Exact Schur-complement recursion gives

```text
visited by depth:         1,7,49,77,84,70,266,469,273,112
Schur-pruned by depth:    0,0,38,65,74,32,199,430,257,58
exact weight solutions:  4
nonintegral roots:       0
rejected by old walls:   4
accepted simple roots:   0
```

All four reconstructed norm-two roots fail the inequality against the same
previous simple wall, global index `9`.  Their exact positive pairings with
that wall are respectively

\[
1+\sqrt2,\quad \sqrt2,\quad \sqrt2,\quad \sqrt2.
\]

The manifest records every rejected root and its blocker index.  Hence this
is an exhaustive empty shell, not a timeout or a numerical failure.

## Endpoint status and literal completion criterion

The targeted edge remains

```text
E = (1,2,4,5,6,7,9,12),   type A2+F4+B2.
```

Among the 34 accepted walls through `x=1+2sqrt(2)`, wall `33` remains its
only certified finite endpoint.  No wall is added at the new shell, so no
second finite endpoint and no new ideal-endpoint data arise there.

The edge is complete only after the complete simple-root list is known and
one of the following supplies a distinct second endpoint:

1. a further accepted wall `r` for which the principal Gram matrix on
   `E union {r}` is positive definite of rank nine; or
2. a full componentwise-affine zero-neighbourhood incident to `E`, of total
   parabolic rank eight, counted once as an ideal endpoint.

Even closing this one edge would not finish the chamber.  Vinberg's literal
finite-volume test must then verify exactly two finite-or-ideal endpoints for
every elliptic rank-eight edge and full neighbourhoods for every connected
affine subdiagram.  Compactness additionally forbids all affine endpoints.

The formal status therefore remains

```text
undecidable_incomplete_simple_roots
```

No conclusion about compactness, finite volume, or infinite volume is drawn
from the empty shell.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_silver_k3_vinberg_1_plus_2sqrt2.py
PYTHONPATH=src:. pytest -q tests/test_silver_k3_vinberg_1_plus_2sqrt2.py
```
