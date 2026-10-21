# The complete local silver Weyl star is non-Coxeter

## Outcome

The first silver-Alexeev obstruction is not an accident of one side of one
ambient face.  Its bad ridge lies in an ambient spherical face of type

\[
A_1+A_2.
\]

The complete local star has twelve ambient Weyl chambers.  Exactly six of
those chambers cut the fixed silver eigenspace in a full-dimensional
hyperbolic-nine sector.  **All six sectors fail the Coxeter angle gate at the
common ridge.**  Thus the original framing and all five alternative
full-dimensional framings in this local star are excluded exactly.

This remains a local theorem.  Weyl framings outside this twelve-chamber star
are not excluded.

## The spherical face and its finite local group

Use the length-108 framing and the exact timelike ridge witness from
`silver_alexeev_h9_folding_obstruction.json`.  It pairs to zero precisely with
Alexeev walls `2`, `4`, and `5`, and strictly positively with the other 47
walls.  Their ambient Gram matrix is

\[
\begin{pmatrix}
-2&0&0\\
0&-2&1\\
0&1&-2
\end{pmatrix},
\]

so the local standard parabolic subgroup is finite:

\[
W_J\cong W(A_1)\times W(A_2),\qquad |W_J|=12,
\quad J=\{2,4,5\}.
\]

The certificate enumerates all twelve integral matrices in this group.  The
common witness is fixed by every one of them.  Consequently, for every
chamber `uC` in the star it is still zero on the three transformed walls and
strict on every other transformed wall.  No floating-point chamber or facet
decision enters the proof.

## Exact sector test

Let `E` be the positive `sqrt(2)` eigenspace and let `n_2,n_4,n_5` be the
three projected transverse normals for a local chamber.  They span a
negative-definite plane.  Their exact `3 x 3` Gram matrix has rank two, and
the certificate computes its one-dimensional kernel

\[
k_2n_2+k_4n_4+k_5n_5=0
\]

over `Q(sqrt(2))`.

The strict Gordan alternative now gives a purely exact decision:

- a full transverse sector exists iff the three coefficients contain both
  signs;
- if one coefficient is zero, the other two walls project to the same
  oriented facet;
- if no coefficient is zero, the uniquely signed normal is redundant and
  the other two normals bound the sector.

The twelve shortest lexicographic words below use
`u=s_(i_1)...s_(i_k)` and the transformed chamber `uC`.

| word | full `H^9` sector | projected facet groups | `cos^2(theta)` |
|---|---:|---|---|
| `[]` | yes | `{2,5}`, `{4}` | `1/4 + sqrt(2)/8` |
| `[2]` | no | -- | -- |
| `[4]` | yes | `{4}`, `{5}` | `3/8 - sqrt(2)/4` |
| `[5]` | no | -- | -- |
| `[2,4]` | no | -- | -- |
| `[2,5]` | yes | `{2,5}`, `{4}` | `1/4 + sqrt(2)/8` |
| `[4,5]` | yes | `{2,4}`, `{5}` | `1/4 + sqrt(2)/8` |
| `[5,4]` | no | -- | -- |
| `[2,4,5]` | no | -- | -- |
| `[2,5,4]` | yes | `{4}`, `{5}` | `3/8 - sqrt(2)/4` |
| `[4,5,4]` | no | -- | -- |
| `[2,4,5,4]` | yes | `{2,4}`, `{5}` | `1/4 + sqrt(2)/8` |

The six negative decisions are also global for the intersection with their
respective ambient chambers, not merely failures of an infinitesimal test.
If a strict point existed elsewhere in the convex slice, the segment from
the common face point to that strict point would produce a strict transverse
sector locally, contradicting the exact Gordan certificate.

## Both angles are non-Coxeter

Four full sectors have

\[
\cos^2\theta=\frac14+\frac{\sqrt2}{8},
\qquad
\frac14<\cos^2\theta<\frac12,
\]

strictly between the `m=3` and `m=4` Coxeter values.  The other two have

\[
\cos^2\theta=\frac38-\frac{\sqrt2}{4},
\qquad
0<\cos^2\theta<\frac14,
\]

strictly between the `m=2` and `m=3` values.  Therefore none is `pi/m` for an
integer `m >= 2`.

## Scope and volume warning

The exact conclusion is:

> Every full-dimensional silver `H^9` slice among the twelve Alexeev Weyl
> chambers incident to this `A1+A2` face is non-Coxeter.

This does not enumerate the global Weyl orbit, change the silver action, or
construct the chamber of the intrinsic norm-two silver reflection group.
Those directions remain open.

Vinberg's finite-volume criterion cannot be applied to any of the six local
slices: the prerequisite that the boundary be an acute Coxeter polyhedron
already fails.  The theorem makes no finite- or infinite-volume claim for a
different silver chamber.

## Reproduction

The frozen certificate is

- `data/candidates/silver_alexeev_h9_local_star.json`.

Its SHA-256 certificate digest is

```text
fe055505befad84e1cb0e9956b54fbdeb88184e645bbcc423a14e42ca655225e
```

Replay all twelve chambers with

```sh
PYTHONPATH=src:. python3 scripts/certify_silver_alexeev_h9_local_star.py
```

The base chamber and ridge data come from the preceding
[`silver_alexeev_h9_folding_obstruction`](silver_alexeev_h9_folding_obstruction.md)
certificate, reconstructed from Alexeev's published `(20,2,1)` chamber in
[Reflective hyperbolic 2-elementary lattices](https://arxiv.org/abs/2209.09110).
