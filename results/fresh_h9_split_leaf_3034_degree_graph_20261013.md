# Degree restriction tightens the exact cap compatibility graph

The second graph retains **3,192,114 possible edges** among the 3,040
archived cap IDs, representing 3,034 distinct mirrors. It certifies
**1,427,166 nonedges**, including **96 additional exclusions** beyond the
preserved [all-label supergraph](fresh_h9_split_leaf_3034_safe_graph_20261013.md).
Every rejected pair was replayed with general interval multiplication.
This is a compatibility supergraph, not a finite-volume conclusion.

The [graph data](../data/enumerations/fresh_h9_split_leaf_3034_degree_graph_20261013.json)
uses the same `adjacency_hex` format: `int(adjacency_hex[i],16)` is the
bitset of possible neighbors of zero-based archive ID `i`.
`additional_nonedges` lists the 96 newly excluded pairs explicitly.
The original graph is unchanged, and replay checks that the new graph
is a subgraph of it.

The [script](../scripts/fresh_h9_split_leaf_3034_degree_graph_20261013.py)
and [replay result](fresh_h9_split_leaf_3034_degree_graph_20261013_verify.json)
are standard-library-only. The source normals are anchored by their file
SHA-256 digest. Their norm-two, seed compatibility, inertia, and strict
center properties were separately established by the
[independent normal audit](fresh_h9_split_leaf_3040_independent_20261013.md).

## Why only fourteen finite labels are possible

Let `K=Q(phi)`, where `phi=(1+sqrt(5))/2`. Every archived normal and its
covector have coefficients in `K(sqrt(D))`, for a literal `D` in `K`.
The pairing of two normals therefore belongs to

```
F = K(sqrt(D),sqrt(E)).
```

Each adjoining step has degree one or two. Consequently `[F:Q]` is 2,
4, or 8, and the degree of any individual pairing divides 8 by the
tower law. This remains true for zero or square radicands and for equal
quadratic extensions; no radicand-independence assumption is made.

Both cap norms are two, so a finite Coxeter angle requires the pairing
`-2*cos(pi/m)` for an integer `m>=2`. Put `z=exp(pi*i/m)`, a primitive
`2m`-th root of unity. The cyclotomic field has degree `totient(2m)`;
its real subfield `Q(z+z^-1)` has index two, since `z` satisfies
`X^2-(z+z^-1)X+1` and is nonreal. Thus

```
[Q(2*cos(pi/m)):Q] = totient(2m)/2.
```

That degree must divide 8, so `totient(2m)` belongs to `{2,4,8,16}`.
Every odd prime divisor `p` of `2m` then has `p-1` a power of two at
most 16, leaving only `p=3,5,17`. Each such prime has exponent at most
one: a larger exponent would make the totient divisible by the odd
prime. Write

```
2m = 2^a * 3^b * 5^c * 17^d,
a>=1, b,c,d in {0,1}.
totient(2m) = 2^(a-1+b+2c+4d).
```

The exponent must be between 1 and 4. Exhausting this finite exponent
budget, rather than imposing a cutoff on `m`, gives exactly

```
m in {2,3,4,5,6,8,10,12,15,16,17,20,24,30}.
```

The script independently enumerates and archives each exponent tuple
and corresponding real cyclotomic degree. The degree test is necessary,
not sufficient; in particular label 17 is retained conservatively.
No unbounded high-label tail remains necessary.

## Exact interval gate and replay

The coefficient and angle enclosures are unchanged from the earlier
graph. Coefficients are enclosed at scale `10^24`, after integer-square
radical bounds at scale `10^80`. Rigorous finite-angle enclosures come
from Machin's rational arctangent identity and an outward-rounded cosine
Taylor series. No floating-point trigonometry is used.

A pair is rejected only when its entire closed pairing interval is
strictly positive, or when it lies above `-2` and misses every enclosure
of the fourteen possible finite values. Every interval that could
contain a parallel/dotted pairing, an allowed finite angle, or an
unresolved equality remains an edge.

| Classification | Pairs |
|---|---:|
| Retained possible edges | 3,192,114 |
| Rejected: strictly positive | 509,018 |
| Rejected: finite-angle gap | 918,148 |
| Total | 4,619,280 |

Replay checks all 787 distinct coefficient enclosures using a separate
exact Fraction field-sign implementation. It then rechecks **all
1,427,166 nonedges**, including the 96 additions, using general
four-product interval arithmetic rather than the construction's
sign-optimized multiplication. It also checks symmetry, the adjacency
digest, all counts, exact containment in the earlier graph, retention of
every allowed-angle enclosure, and rejection of every excluded label's
enclosure through 120. The latter is a regression check; completeness
comes from the exponent-budget proof, not that finite test range.

Both construction and replay passed, as did Python compilation.

```sh
python scripts/fresh_h9_split_leaf_3034_degree_graph_20261013.py
python scripts/fresh_h9_split_leaf_3034_degree_graph_20261013.py --verify
```

The adjacency digest is
`85233677a8d57e40445e6d33a3bf0c4eb9a64e163c31b05a5bbb5aa410332391`.

The six duplicate profile IDs remain for compatibility with earlier
selection data. A geometric facet set can use one representative per
mirror. A retained clique still requires exact compatibility and a
volume test. An unsatisfiability proof using this graph and independently
certified necessary escape-cut clauses could exclude the literal bank;
this graph alone does not. No claim concerns mirrors outside the bank
or arbitrary completions of the seed.
