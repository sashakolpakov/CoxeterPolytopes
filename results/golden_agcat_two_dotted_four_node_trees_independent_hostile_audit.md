# Independent hostile audit: golden two-dotted four-node trees

This is a clean-room audit: its certifier imports no project module and
rebuilds `Z[phi]`, the `4 < x < 62^2`, `0 < x' < 4` dotted alphabet, exact
real-place signs, and all counts from integer pairs `(a,b)`.

The original blanket rejection of the disjoint-dotted path was false.  Its
two indefinite edge blocks are coupled by the finite middle edge, so their
negative indices do not add.  For example, with finite label `3`, endpoint
products

```text
x = 3 + phi,   z = 2 + 2 phi
```

give full determinant `-2 phi`: it is negative at the selected embedding
and positive at the conjugate embedding.  Exact leading-principal-minor
signs give selected inertia `(3,1,0)` and conjugate inertia `(4,0,0)`.

After that correction, the audit passes.  It independently obtains 6,869
dotted values and the same conjugate-order SHA-256, and brute `S_4` orbit
enumeration on small alphabets proves the roles used by the closed count:
stars and disjoint paths use unordered dotted pairs, while adjacent paths
use ordered remote/adjacent pairs.  Thus the raw count is

```text
5 * (2*C(6869+1,2) + 6869^2) = 471,865,955.
```

The exact retained rows `(star, disjoint, adjacent)` for finite labels
`3,4,5,6,10` are respectively

```text
(6,634,798,  7, 17,689,960)
(2,948,447,  5, 11,792,207)
(9,650,285, 44, 21,334,861)
(  737,201,  2,  5,895,415)
(5,052,550, 23, 15,436,911).
```

Therefore exactly **97,172,716** packets survive: 25,023,281 stars, 81
disjoint-dotted paths, and 72,149,354 adjacent-dotted paths.  All comparisons
are strict and exact.  The conjugate-PD formulas follow from Schur complements
and Sylvester's criterion; the selected disjoint path has index one exactly
when its full determinant is negative.

The integral Cartan/symmetrizer claim also passes.  Orient each tree edge
product `c` as `(-1,-c)`.  With the product of all edge products at a root,
divide recursively by `c` along oriented edges.  This gives an integral,
totally positive diagonal symmetrizer, with no cycle-consistency condition.

Strict scope: this audit covers connected four-node trees with exactly two
dotted edges.  It does not enumerate cyclic two-dotted packets.
