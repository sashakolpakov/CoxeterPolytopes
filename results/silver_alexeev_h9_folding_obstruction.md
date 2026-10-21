# The natural silver slice of Alexeev's `(20,2,1)` chamber is not Coxeter

## Outcome

The silver divisor-lattice design reaches a genuine singular K3 carrier, but
the most direct attempt to inherit a Coxeter chamber from that surface fails
an exact angle test.

Start with the one-4 Lannér cycle

\[
C_5(4,3,3,3,3),
\]

saturate its trace lattice by its unique nontrivial even
`sqrt(2)`-stable index-two class, and add the positive silver block
`F4 + A1`.  The integral trace lattice is

\[
II_{1,9}\oplus E_8(-1)\oplus A_1(-1)^2.
\]

This is the unique even 2-elementary K3 lattice with invariants `(20,2,1)`,
the lattice of Alexeev's known 50-wall finite-volume noncompact reflection
chamber in `H^19`.

The construction supplies an explicit integral self-adjoint operator `tau`
with

\[
\tau^2=2,
\qquad
(\tau x,y)=(x,\tau y),
\qquad
(\tau x,\tau y)=2(x,y).
\]

Its `sqrt(2)`-eigenspace has signature `(1,9)`.  It is therefore a natural
candidate `H^9` section of the algebraic reflection chamber.

It is **not** a Coxeter section.

## Exact Weyl framing

Alexeev's displayed chamber does not meet the initially marked positive
eigenspace.  Consequently it would be invalid simply to project the
published diagram.  A regular positive eigenvector was first moved into the
50-wall chamber.  The resulting Weyl word has length 108.  Every step is an
integral reflection in one of the published `(-2)` or `(-4)` roots, and the
stored word verifies exactly that

\[
w^T G w=G,
\qquad
\tau'=w\tau w^{-1},
\qquad
(\tau')^2=2,
\qquad
(\tau')^TG=G\tau'.
\]

This chamber framing is a necessary gate: before framing, the eigenspace and
the selected chamber have empty intersection.

## A non-Coxeter ridge

Let `r_2` and `r_4` be walls 2 and 4 in zero-based indexing.  Four times the
Gram matrix of their positive-eigenvalue projections is

\[
\begin{pmatrix}
-4+2\sqrt2 & \sqrt2\\
\sqrt2 & -4
\end{pmatrix}.
\]

Hence their hyperplanes meet with

\[
\cos^2\theta
=\frac{(\sqrt2)^2}{(-4+2\sqrt2)(-4)}
=\frac14+\frac{\sqrt2}{8}.
\]

This is strictly between the consecutive Coxeter values

\[
\cos^2(\pi/3)=\frac14,
\qquad
\cos^2(\pi/4)=\frac12,
\]

so `theta` is not `pi/m` for any integer `m`.

The certificate also proves that this is a real ridge of the restricted
polyhedron rather than an irrelevant intersection outside it.  In a fixed
ten-column basis of the `sqrt(2)`-eigenspace, the stored point has
coefficients

\[
\begin{aligned}
(&-69-49\sqrt2,\;-26-19\sqrt2,\;-96-68\sqrt2,\\
 &-74-53\sqrt2,\;-48-35\sqrt2,\;121+84\sqrt2,\\
 &142+100\sqrt2,\;90+63\sqrt2,\;138+97\sqrt2,
   -1-\tfrac12\sqrt2).
\end{aligned}
\]

Its square is

\[
2238+1601\sqrt2>0.
\]

It pairs to zero exactly with walls 2, 4, and 5, and strictly positively
with the other 47 published walls.  Walls 2 and 5 have the same projected
mirror, so after removing that duplicate, walls 2 and 4 are two independent
facets meeting along the displayed hyperbolic ridge.

## What this proves—and what it does not

This closes one particularly rigid construction:

> The naturally marked silver action, framed by the explicit regular
> eigenvector above, does not turn Alexeev's 50-wall `(20,2,1)` chamber into a
> Coxeter polytope in `H^9`.

It does not exclude every Weyl framing, a different silver action on the same
K3 lattice, or the chamber of the norm-two silver reflection group itself.
Those are separate finite searches.  It also does not weaken Alexeev's
finite-volume theorem for the ambient `H^19` chamber.

A subsequent exact local-star audit strengthens the first clause: all six
full-dimensional silver slices in the complete twelve-chamber `A1+A2` Weyl
star at this ridge are non-Coxeter.  See
[`silver_alexeev_h9_local_star.md`](silver_alexeev_h9_local_star.md).  Weyl
framings outside that finite local star remain open.

Most importantly for volume claims, Vinberg's edge/cusp criterion is **not
applicable** to this 16-facet restricted boundary: the acute Coxeter-Gram
prerequisite already fails.  Declaring it finite-volume from inherited cusp
data would therefore be invalid.

## Reproduction

The self-contained frozen certificate is

- `data/candidates/silver_alexeev_h9_folding_obstruction.json`.

Its SHA-256 certificate digest is

```text
d4a8502d51ea7b33225af999825d31a18506f0ac58ee86d96f4db1069daacb23
```

Replay it with

```sh
PYTHONPATH=src:. python3 scripts/certify_silver_alexeev_h9_folding_obstruction.py
```

The exploratory reconstruction and exact Weyl framing are in
`scripts/probe_silver_alexeev_h9_slice.py`.  Rebuilding from the primary
ancillary file additionally verifies its published SHA-256 digest.  The
ambient matrix is from Alexeev,
[Reflective hyperbolic 2-elementary lattices, arXiv:2209.09110v4](https://arxiv.org/abs/2209.09110).
