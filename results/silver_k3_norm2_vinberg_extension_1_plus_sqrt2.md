# First intrinsic silver K3 Vinberg extension: `x=1+sqrt(2)`

> **Subsequent extension.**  An exact two-embedding Schur census reduces the
> unique zero-endpoint edge to eleven possible incidence patterns.  The next
> shell `x=2+sqrt(2)` adds eighteen walls, exactly one of which supplies a
> finite endpoint.  See
> [`silver_k3_norm2_vinberg_edge_extension.md`](silver_k3_norm2_vinberg_edge_extension.md).

This is the smallest exact continuation of the intrinsic norm-two search in
[`silver_k3_norm2_vinberg_prefix.md`](silver_k3_norm2_vinberg_prefix.md).
It exhausts the first admissible control-coordinate level beyond `x=2`; it
does **not** claim that the simple-root list is complete.

The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_extension_1_plus_sqrt2.json`.  Its
canonical certificate digest is

```text
f77614309c85b85366c829aacdfd8f03a028866b293ca5325b8646983f72cf81
```

## Why this is the next level

Write `x=a+b sqrt(2)`.  The norm-two equation at the positive conjugate
embedding gives

\[
(x')^2\leq 2(1+\sqrt2).
\]

An exact enumeration in `Z[sqrt(2)]`, simultaneously imposing `0<x<=3`,
leaves precisely

\[
1,\quad \sqrt2,\quad 2,\quad 1+\sqrt2.
\]

The first three levels are the previous certificate.  Thus
`x=1+sqrt(2)` is the unique unseen level through identity bound three.

## Exact shell enumeration

For fixed `x`, the nine nonnegative finite-simple-root pairings `z` obey

\[
z^T K^{-1}z=2+(1+\sqrt2)x^2.
\]

At `x=1+sqrt(2)`, the two real embeddings give an eight-element coordinate
alphabet.  Exact Schur-complement pruning gives the following replay data:

```text
visited by depth:       1,8,64,104,120,112,456,904,512,184
pruned by depth:        0,0,51,89,106,55,343,840,489,103
exact weight solutions: 6
nonintegral roots:      0
rejected by old walls:  2
accepted simple roots:  4
```

The four new roots are

```text
(1+sqrt2, 2+sqrt2, 2+2sqrt2, 4+3sqrt2, 4+3sqrt2,
 0, 0, 0, 0, 1+sqrt2)

(1+sqrt2, 2+sqrt2, 2+2sqrt2, 4+3sqrt2, 4+3sqrt2,
 2+sqrt2, 4+2sqrt2, 3+3sqrt2, 2+2sqrt2, 0)

(1+sqrt2, 2+sqrt2, 2+2sqrt2, 4+3sqrt2, 4+3sqrt2,
 2+2sqrt2, 3+3sqrt2, 4+2sqrt2, 2+sqrt2, 0)

(1+sqrt2, 2+2sqrt2, 4+2sqrt2, 5+3sqrt2, 4+3sqrt2,
 2+sqrt2, 3+2sqrt2, 3+2sqrt2, 2+sqrt2, 0)
```

Together with the previous twelve walls they form a 16-wall acute prefix.
Its Gram inertias are

\[
(9,1,6)\quad\hbox{and}\quad(10,0,6)
\]

at the identity and conjugate embeddings.

## Exact diagram diagnostic

Before an exact inertia test, the verifier uses only the following necessary
principal-minor condition.  Since all roots have square two and the Gram is
acute, every positive-semidefinite principal submatrix has pairings at least
`-2`, and every positive-definite one has pairings strictly greater than
`-2`.  This is an exact pruning lemma, not a floating-point classification;
all surviving subsets are still checked by exact silver-field inertia.

The complete subset audit gives

```text
elliptic rank-8 edge diagrams:       267
elliptic rank-9 finite vertices:      41
finite endpoint histogram:             1 edge with 0 endpoints
                                      163 edges with 1 endpoint
                                      103 edges with 2 endpoints
connected affine subdiagrams:           0
```

The unique zero-finite-endpoint edge is the wall-index set

```text
(1,2,4,5,6,7,9,12)
```

and its rank-eight Gram determinant is `6`.  It is a concrete target for a
later shell, not an obstruction: unseen walls can still supply finite or
ideal endpoints.  Likewise, the increase from 48 open edges in the 12-wall
prefix to 164 in the 16-wall prefix is not paradoxical.  Adding facets creates
new spherical rank-eight subsets, so open-edge counts of incomplete prefixes
need not decrease.

The formal status remains

```text
undecidable_incomplete_simple_roots
```

In particular, neither compactness, finite covolume, nor infinite covolume
follows from this bounded prefix.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_silver_k3_vinberg_extension.py
PYTHONPATH=src:. pytest -q tests/test_silver_k3_vinberg_extension.py
```
