# Golden AGCat: exactly-two-dotted four-node trees

This certificate closes the complete connected four-node **tree** tranche
with exactly two dotted edges.  It uses the full 6,869-value dotted golden
alphabet and all five nonorthogonal finite Coxeter products.

Modulo graph automorphisms there are 471,865,955 canonical assignments before
the two-place Gram gate.  Exact conjugate-order counting retains 97,172,716:

- 25,023,281 stars;
- 81 paths with disjoint dotted edges (a thin determinant-sign locus);
- 72,149,354 paths with adjacent dotted edges.

The cycle-squareclass condition is vacuous on a tree, and an integral Cartan
gauge is obtained by orienting every edge product `c` as `(-1,-c)`.  At the
selected embedding the retained star and adjacent-path cases are automatically
Lorentzian.  At the conjugate embedding their positive-definiteness conditions
reduce exactly to

```text
star:          x' + y' + z' < 4,
adjacent path: (4-x')(4-z') - 4y' > 0.
```

The counts use exact arithmetic, monotone two-pointer joins, and an exact
two-dimensional dominance join, with no floating point.  Their abundance
shows that local AGCat compatibility does
not itself prove golden nonexistence.  Cyclic two-dotted packets, packets with
more dotted edges, global ranks 23--30, orders/conductors/stable lattices, and
complete Vinberg hosts remain open.  This certificate is neither a compact
example nor a global nonexistence theorem.
