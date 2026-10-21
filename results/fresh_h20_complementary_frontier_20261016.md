# The135-wall H20 frontier: a persistent null ray and an exhausted old bank

The new strict escape E1 survives **all847 archived exact nonright-S61
recipes**, even when their mutual incompatibilities are ignored. Thus
switching among old-bank subsets alongside cap133 cannot remove this gap.
The complementary wall also creates a useful exact local target: a
persistent ideal point with a4-dimensional unbounded direction space in
its horospherical link.

The [checker](../scripts/fresh_h20_complementary_frontier_20261016.py) and
[result](fresh_h20_complementary_frontier_20261016.json) pass generation
and fresh saved replay. Exact geometric conclusions and the bounded
numerical compatibility diagnostic are explicitly separated. No frozen
predecessor is modified; no additional wall is asserted here.

## An exact common escape for the entire old recipe bank

Use the vector from the [135-wall certificate](fresh_h20_complementary_cap_20261016.md):

```
E1=x-e/100-(49/50)*f-omega17
   -(1/100000)*sum_{j!=17}(omega_j),
E1²=1842000049/4000000000>0.
```

The old [bank gate](fresh_h20_nonright_s61_20261014.md) independently of
its numerical discovery coordinates defines every recipe by named exact
support contacts, two exact anchor equations, norm2 and an isolated
quadratic branch. It saves rigorous rational intervals for all847
residual frames. This new certificate evaluates E1 against those exact
intervals, reconstructing each support coefficient from its label3,4,5
or6 using rational square-root enclosures. Every pairing is strictly
negative. The largest upper bound is below−0.20769, at recipe677.

The same old strict center works for the bank and cap133, and E1 is
already strict for all135 actual walls. Consequently the relaxed
intersection of all847 recipe halfspaces, the133 original walls and
cap133 has infinite volume. Every compatible bank subset has the same
escape. This remains a finite archived-bank theorem; it does not exclude
new profiles, new anchor equations or other added walls.

## Bounded compatibility diagnostic

Using exact-recipe interval midpoints and tolerance10^-10,765 of the847
recipes pass the numerical mutual-angle screen with cap133. Among them,
476 cut E0, but none cut E1. Their induced compatibility graph has3978
edges and maximal-clique counts:

| Size | Number |
|---|---:|
|1|173|
|2|3975|
|3|1|

The sole triple is `[100,805,810]`. Proposal44, now exactly certified as
cap134, has no distinct compatible old-bank partner in this screen.
Thus the bank yields no candidate to add directly to the current135-wall
branch. Alternative branches can retain more old-bank walls: for example,
proposal468 cuts E0 and has67 compatible neighbors after the cap133 gate.
Its contacts are π/3 at support roots5 and87, with anchor labels4 and6
at65 and95. None of those alternatives defeats the exact E1 bank gate.

These graph counts are **numerical diagnostics**, not exact angle
acceptance or exclusion theorems. The exact common-escape conclusion does
not depend on this graph or a finite angle-search cutoff.

## Complete incidence at a new future null ray

Let r denote cap134. Its exact contacts with old roots0 and1 give the
affine diagram

```
(r,0,1): [[2,-√3,0],[-√3,2,-1],[0,-1,2]],
positive kernel (√3,2,1).
```

Therefore

```
nu=√3*r+2*n0+n1
```

is null. It is nonzero and future-oriented because its pairing with the
common future center is a strictly negative positive combination of
three strict wall inequalities. For any retained wall outside these
three, acuteness makes every summand in its pairing with nu nonpositive.
Equality is possible exactly when that wall is right to all three.

The complete incident set is

```
{134,0,1} union {3,4,6,7,8,9,10,12,13,14,15,17,31}.
```

The checker reconstructs all needed old physical pairings exactly; the
only other right neighbors of r, namely2,61,87,95, have a negative
pairing with0 or1 and are excluded. All remaining old walls have negative
pairing with r itself, so cannot be incident. This includes cap133.

The incident diagram is

```
affine G2 + A2[3,4] + A5[6,...,10] + A4[12,...,15]
          + A1[17] + A1[31].
```

Its Euclidean normal rank is2+13=15 in the19-dimensional horospherical
link. Hence the link has a4-dimensional lineality space and is unbounded.
This is not a proper bounded cusp, and no complete ambient ideal-point
census is claimed.

## A concrete necessary completion condition

Any further coherent acute wall w retaining the affine triple satisfies

```
w.nu=√3*(w.r)+2*(w.n0)+(w.n1)<=0.
```

Thus **no such wall can cut nu off**. If w passes through nu, all three
summands must vanish, so w must be right to the entire affineG2 component.
Walls strictly negative at nu do not affect a sufficiently small local
horosphere. A finite-volume completion must therefore add incident walls
whose projections positively span the4 currently free directions.
At least5 additional incident walls are necessary; at most4 homogeneous
halfspace normals cannot bound a4-dimensional Euclidean space.

This gives a finite geometric target for the next construction: complete
the remaining five spherical paths A2,A5,A4,A1,A1 inside the17-dimensional
factor orthogonal to affineG2. It is a local cusp-completion problem,
separate from cutting the surviving spacelike escape E1.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_complementary_frontier_20261016.py --verify
```
