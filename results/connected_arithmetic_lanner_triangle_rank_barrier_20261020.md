# Compact Lanner exchange and the unbounded connected-star barrier

20 October 2026. Own exact paper calculation; no scientific execution.

Compactness supplies a short global exchange path, but connectedness does
not bound local continuation rank. Both points have elementary proofs.

## 1. Every facet has a Lanner exchange path of length at most five

Let `P` be a compact hyperbolic Coxeter polytope and fix any Lanner
principal diagram `L`. Then every facet node `u` belongs to a Lanner
diagram `L_u`, and its graph distance from `L` is at most five. Edges in
this statement mean all nonorthogonal contacts, including dotted bonds.

Indeed a compact Coxeter polytope is simple, and its dual nerve is a
simplicial sphere. Its faces are the elliptic principal subsets; its
minimal nonfaces are the Lanner subsets, of size at most five. If a node
`u` belonged to no minimal nonface, adjoining `u` to any face not already
containing it would still give a face. The nerve would be a cone at `u`,
contrary to its sphere homology. This proves the asserted `L_u`.

Two disjoint Lanner diagrams cannot be mutually orthogonal: their direct
sum would have two negative directions in the ambient Lorentz space.
If `L_u` meets `L`, follow a path inside connected `L_u` to the
intersection, using at most four edges. Otherwise some node of `L_u` has
a nonorthogonal contact with `L`; reaching it inside `L_u` and then
taking that contact uses at most five edges. This does not assert that
the intervening nodes are simultaneous white nodes of one actual face.

## 2. A genuine compact triangular face with unbounded connected star

The following is a countermodel to a proposed **local** rank reduction,
not a compact example or a dimension-bound counterexample. Even the
complete connected star of an actual compact triangular face can satisfy
integral reflection and nonphysical-positive-definiteness conditions in
arbitrarily high dimension. Thus the earlier rank-13 calculation for a
single-4 cycle with a pendant chain cannot be extended to all connected
Lanner continuations.

For any integer `m >= 1`, put `n = m + 2` and take the basis
`a,b,c,w_1,...,w_m` over `K = Q(sqrt(3))`. Give every basis root square two.
The only nonzero off-diagonal products are

\[
 (a,b)=-\sqrt3,\qquad (a,c)=(b,c)=-1,\qquad
 (c,w_1)=(w_i,w_{i+1})=-1.
\]

The three black roots form the compact Lanner triangle `(3,3,6)`. The
white support `W = {w_1,...,w_m}` is `A_m` and is attached to `c` by a single
simple bond. The entire support graph is connected.

### Exact signatures and arithmetic data

Eliminating `W` leaves the black Gram

\[
 \begin{pmatrix}
 2&-\sqrt3&-1\\
 -\sqrt3&2&-1\\
 -1&-1&t
 \end{pmatrix},\qquad t=\frac{m+2}{m+1}\in(1,2).
\]

Its physical determinant is `t - 4 - 2 sqrt(3) < 0`; the leading two-node
block is positive definite. Consequently the full Gram has signature
`(m+2,1)`. Under the other embedding of `K`, its determinant is
`t - 4 + 2 sqrt(3) > 2 sqrt(3) - 3 > 0`; the first two leading minors
are again positive. The full conjugate Gram is therefore positive
definite for every `m`.

Every Gram entry belongs to `Z[sqrt(3)]`, and every basis-root reflection
has formula `x -> x - (x,r)r`, so preserves this integral lattice. The
length-three cyclic product is `-sqrt(3)`, and hence the cyclic field is
exactly `K`. These are genuine admissible integral reflection data; no
bound on field degree or reflection integrality repairs the local rank
argument.

### The triangular face is actual and compact

Deleting `c` gives the finite diagram `I_2(6) + A_m`. Deleting either `a`
or `b` gives `A_(m+2)`. Thus all three complete proposed corner diagrams
are positive definite, regardless of `m`.

The physical Gram is `2I - A` with `A` nonnegative and irreducible.
Perron--Frobenius supplies a positive-coefficient timelike vector having
strictly negative product with every wall. This is an interior point of
the acute halfspace chamber. Projection to any positive-definite wall
subset remains timelike and feasible for every other wall: the inverse
of its positive-definite Gram with nonpositive off-diagonal entries is
entrywise nonnegative. In particular the three rank-`n` subsets just
listed give actual finite vertices, strict on their omitted black wall.
Within `W^perp`, their convex hull is the compact hyperbolic triangle
bounded by `a,b,c`. All walls of this chamber have an actual relative
interior, by the same projection argument for a single wall.

### The global compactness gap remains essential

The chamber is **infinite volume**. Deleting the terminal white root
`w_m` leaves the same Lorentz construction with a chain of length `m-1`
(the bare `(3,3,6)` triangle when `m=1`). Thus the dual vector `v`
satisfying `(v,w_m)=-1` and zero products with every other wall is
spacelike: its square is the ratio of that negative principal determinant
to the negative full determinant. Joining an interior timelike vector to
a sufficiently large positive multiple of `v` crosses the light cone
while all wall inequalities remain strict. An open ideal boundary patch
therefore survives, which proves infinite volume.

No compact completion by further walls is asserted. The model shows
precisely what the moving-core route must add: **global compactness or
coverage of the other escaping directions**, beyond connectedness,
admissibility, actual compactness of one small face, and its three full
finite corner diagrams. It does not invalidate the separately scoped
single-4 rank-13 calculation.
