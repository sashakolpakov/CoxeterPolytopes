# Intrinsic silver K3 Vinberg prefix through `x <= 2`

> **Subsequent extension.**  The first unseen admissible level is
> `x=1+sqrt(2)`.  It adds four walls, producing an exact 16-wall prefix which
> still does not satisfy a Vinberg stopping criterion.  See
> [`silver_k3_norm2_vinberg_extension_1_plus_sqrt2.md`](silver_k3_norm2_vinberg_extension_1_plus_sqrt2.md).

This is an exact bounded result for the Alexeev-aligned silver K3 carrier.
It studies the intrinsic norm-two reflection system over

\[
O=\mathbf Z[\sqrt2],
\]

not a restriction of Alexeev's 50-wall chamber.  The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_prefix.json` and is independently
recomputed by `scripts/certify_silver_k3_vinberg_prefix.py`.  Its canonical
certificate digest is

```text
43364d42b68527be4192f905e7c905ee255e4d88569eb1e25e7af9b7518a4dc6
```

## The free saturated field lattice

Put

\[
y={\sqrt2\over2}(e_1+e_2+e_4).
\]

Then `e0,e1,e2,e3,y` is a free `O`-basis of the saturated one-4 `C5`
block.  In this basis its Gram matrix is

\[
\begin{pmatrix}
2&-1&0&0&-\sqrt2\\
-1&2&-\sqrt2&0&-1+\sqrt2\\
0&-\sqrt2&2&-1&-1+\sqrt2\\
0&0&-1&2&-\sqrt2\\
-\sqrt2&-1+\sqrt2&-1+\sqrt2&-\sqrt2&3-\sqrt2
\end{pmatrix}.
\]

Its determinant is `-(1+sqrt(2))`, a unit.  Orthogonally adjoining the
field blocks `F4+A1` gives a rank-ten form with signatures

\[
(9,1)\quad\hbox{and}\quad(10,0)
\]

at the identity and conjugate embeddings.

The integral control vector

\[
p=(\sqrt2-1,\sqrt2-1,2-\sqrt2,1,1,0,0,0,0,0)
\]

satisfies

\[
p^2=1-\sqrt2,\qquad (p,r)=-x,
\]

where `x` is the first coordinate of `r`.  For a norm-two root,

\[
\delta(r)={1+\sqrt2\over2}x^2.
\]

## Height zero is exhaustive

The positive-definite trace lattices perpendicular to `O p` split into
blocks with short-vector counts

| block | trace `(-2)` vectors | field norm-two vectors |
|---|---:|---:|
| saturated `C5` complement | 240 | 32 |
| `F4` trace block | 240 | 48 |
| `A1` trace block | 4 | 2 |

An exact short-vector enumeration, followed independently by reflection
closure, gives the same 82 field roots.  Their nine simple walls have type

\[
B_4+F_4+A_1.
\]

Thus the height-zero part is proved complete; it is not merely the visible
set inherited from the designed diagram.

## Exhaustive positive-height shells

The conjugate norm equation gives

\[
(x')^2\le 2(1+\sqrt2).
\]

Together with `0<x<=2`, this leaves exactly

\[
x=1,\ \sqrt2,\ 2.
\]

For fixed `x`, let `z` be the nine pairings with the finite simple roots.
The enumeration solves exactly

\[
z^T K^{-1}z=2+(1+\sqrt2)x^2,
\]

at both real embeddings, reconstructs the lattice vector, and then applies
all previously accepted wall inequalities.

| `x` | weight alphabet | exact weight solutions | nonintegral | blocked by earlier walls | accepted |
|---|---:|---:|---:|---:|---:|
| `1` | 4 | 1 | 0 | 0 | 1 |
| `sqrt(2)` | 5 | 5 | 3 | 1 | 1 |
| `2` | 4 | 7 | 2 | 4 | 1 |

The three accepted roots, after the nine height-zero walls, are

```text
(1,
 2, 2sqrt2, 2+2sqrt2, 2+2sqrt2,
 0, 0, 0, 0, 0)

(sqrt2,
 sqrt2, 2, 2+sqrt2, 2+sqrt2,
 2+sqrt2, 3+2sqrt2, 3+2sqrt2, 2+sqrt2, 0)

(2,
 2, 2sqrt2, 2+2sqrt2, 2+2sqrt2,
 0, 0, 0, 0, 1+sqrt2)
```

All twelve walls have norm two and are pairwise acute.  Their Gram inertia
is `(9,1,2)` at the identity embedding and `(10,0,2)` at the conjugate
embedding.  The first new wall completes the intended one-4 `C5` Lanner
cycle.  The second completes another one-4 `C5` cycle on the `F4` side and
connects the two cycles by a dotted edge.  The third connects the `A1`
factor.

## Vinberg stopping test: not yet

The exact rank-eight/rank-nine diagnostic for this 12-wall prefix gives

```text
elliptic rank-8 edge diagrams:       69
elliptic rank-9 finite vertices:     10
finite endpoint histogram:           48 edges with 1, 21 edges with 2
connected affine subdiagrams:         0
```

Therefore 48 edge diagrams still have an open end.  Neither compact nor
finite-volume Vinberg stopping is available at this bound.  The formal
status is

```text
undecidable_incomplete_simple_roots
```

The endpoint counts are a stopping diagnostic only.  They are not being
used to infer the volume of the eventual chamber.

## The saturation divisor is not an intrinsic field reflection

The new `(-2)` divisor introduced by saturation is the fifth basis vector
`y`.  Its field norm is

\[
h(y,y)=3-\sqrt2,
\]

and its trace square is `-2`.  However, its reflection coefficient on `e0`
is

\[
{2h(e_0,y)\over h(y,y)}
=-{4+6\sqrt2\over7}\notin O.
\]

Consequently this divisor mirror does not define an integral reflection of
the silver `O`-lattice.  This separates the K3 root-closure warning from the
intrinsic norm-two reflection search: the divisor is a genuine K3
`(-2)`-class, but not a suitable silver reflection root.

## Scope

Proved here:

- the free saturated `O`-lattice presentation;
- complete height-zero norm-two roots;
- every norm-two Vinberg level through `x=2`;
- the exact 12-wall acute prefix and its failed stopping diagnostic;
- nonintegrality of the saturation-divisor reflection.

Not proved here:

- completeness of the eventual norm-two chamber;
- compactness or finite covolume;
- exhaustion of primitive reflective vectors of norms other than two.
