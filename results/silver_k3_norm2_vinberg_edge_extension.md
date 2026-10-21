# The first missing silver edge endpoint appears at `x=2+sqrt(2)`

> **Subsequent extension.**  The next admissible shell
> `x=1+2sqrt(2)` has four exact integral norm-equation solutions, all blocked
> by wall 9, and therefore adds no simple root.  The missing second endpoint
> remains open; the next admissible coordinate is `x=3+sqrt(2)`.  See
> [`silver_k3_norm2_vinberg_extension_1_plus_2sqrt2.md`](silver_k3_norm2_vinberg_extension_1_plus_2sqrt2.md).

This certificate continues the intrinsic norm-two silver K3 search through
the next admissible control coordinate after `x=1+sqrt(2)`.  It combines a
finite, two-embedding Schur-complement census for one open rank-eight edge
with an exact enumeration of the next Vinberg shell.

The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_edge_extension.json`; its canonical
digest is

```text
9538f7d9845b0f31088280438be3843603f9dbbca7da2da450b9cf4c5e314b01
```

## The targeted edge

In the 16-wall prefix, the unique elliptic rank-eight subset having no
finite endpoint was

```text
E = (1,2,4,5,6,7,9,12).
```

Its diagram and determinant are

\[
E=A_2\sqcup F_4\sqcup B_2,
\qquad \det G_E=6.
\]

Let `r` be any future norm-two simple root which extends `E` to an elliptic
rank-nine diagram, and put `z_i=(r,e_i)` for `e_i` in `E`.  The identity
Gram is acute, while the conjugate carrier Gram is positive definite.
Consequently every two-by-two principal minor gives

\[
-2<z_i\leq0,
\qquad |z_i'|<2.
\]

The only elements of `Z[sqrt(2)]` satisfying both inequalities are

\[
z_i\in\{0,-1,-\sqrt2\}.
\]

The rank-nine Schur complement must also satisfy, at both embeddings,

\[
z^T G_E^{-1}z<2.
\]

Checking the resulting `3^8` vectors exactly leaves only eleven patterns.
Writing `0`, `3`, and `4` for the pairings `0`, `-1`, and `-sqrt(2)`, they
are

```text
00000000  00000003  00000030
03000000  03000003  03000030  04000000
30000000  30000003  30000030  40000000
```

In particular, an elliptic endpoint wall cannot meet the `F4` component.
This finite pattern theorem applies to every future norm-two root, not only
to the next shell.

## The next exact shell

The admissible coordinates through the next value are

\[
1,\ \sqrt2,\ 2,\ 1+\sqrt2,\ 2+\sqrt2,
\]

and the following admissible value is `1+2sqrt(2)`.  At
`x=2+sqrt(2)`, exact positive-definite enumeration gives

```text
weight alphabet size:    11
visited by depth:         1,11,121,264,374,396,1936,4532,3047,1199
Schur-pruned by depth:    0,0,97,230,338,220,1524,4255,2938,798
exact weight solutions:  72
nonintegral roots:        43
rejected by earlier walls:11
accepted simple roots:    18
```

Thus the intrinsic norm-two prefix is exact through `x=2+sqrt(2)` and has
34 walls.  It remains acute.  Its two field Gram inertias are

\[
(9,1,24)\quad\hbox{and}\quad(10,0,24).
\]

Exactly one new wall extends the targeted edge elliptically.  It is global
wall index `33`, with lattice coordinates

```text
(2+sqrt2,
 4+2sqrt2, 4+4sqrt2, 6+5sqrt2, 6+4sqrt2,
 0, 0, 0, 0, 1+sqrt2).
```

Its pairings with `E` are

```text
(-sqrt2,0,0,0,0,0,0,0),
```

and the resulting rank-nine Gram has determinant `4` and inertia `(9,0,0)`
at both embeddings.  It is therefore a genuine finite endpoint of this
edge.

## What this does and does not decide

The former zero-endpoint edge now has one certified finite endpoint.  A
second finite or ideal endpoint has not yet been found or excluded.

No global all-subset Vinberg audit is claimed for the incomplete 34-wall
prefix.  In particular, the targeted success does not prove that all other
edges close, that the wall list terminates, or that the eventual chamber is
compact or finite-volume.  The formal status remains

```text
undecidable_incomplete_simple_roots
```

The `x=1+2sqrt(2)` shell is now exhausted and empty.  The exact next target is
`x=3+sqrt(2)`, together with the search for a second finite endpoint pattern
among the same eleven possibilities and, after wall completeness, the full
ideal-endpoint/parabolic-neighbourhood audit.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_silver_k3_vinberg_edge_extension.py
PYTHONPATH=src:. pytest -q tests/test_silver_k3_vinberg_edge_extension.py
```
