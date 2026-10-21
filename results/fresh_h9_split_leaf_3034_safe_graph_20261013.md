# A safe compatibility supergraph for the larger cap bank

The complete graph retains **3,192,210 possible edges** and certifies
**1,427,070 nonedges** among the 3,040 archived cap IDs. Every genuine
Coxeter-compatible pair of distinct cap mirrors remains an edge. Some
retained pairs may still be incompatible; this is deliberately a
supergraph suitable for proving an exclusion from a subsequent exact
unsatisfiability certificate.

The [graph data](../data/enumerations/fresh_h9_split_leaf_3034_safe_graph_20261013.json)
contains `adjacency_hex`: converting entry `i` with `int(value,16)` gives
the bitset of possible neighbors of zero-based profile ID `i`. It retains
all original IDs, including the six duplicate pairs identified in the
[normal audit](fresh_h9_split_leaf_3040_independent_20261013.md). A geometric
facet selection can use one representative per mirror; there are 3,034
distinct mirrors. The source exact-normal certificate is anchored by its
SHA-256 digest.

The [constructor and replay](../scripts/fresh_h9_split_leaf_3034_safe_graph_20261013.py)
use only the standard library. The separate
[verification result](fresh_h9_split_leaf_3034_safe_graph_20261013_verify.json)
records the independent field-sign check of all coefficient enclosures and
the general interval-product replay of every nonedge.

## Exact coefficient and angle bounds

Normal coefficients and covector entries are enclosed by closed intervals
whose endpoints are integers divided by `10^24`. Construction first bounds
`phi` at scale `10^80`, using an integer square root and checking opposite
signs of `x^2-x-1` at the endpoints. It then bounds each positive radical
`sqrt(D)` by integer-square inequalities. Signed interval arithmetic gives
the coefficient enclosures, including the degenerate `D=0` cases.

The replay checks these enclosures using a different method: the
independent Fraction field-sign implementation from the normal audit tests
both exact inequalities at every distinct coefficient interval. Thus the
integer radical-bound construction is not trusted without an algebraic
cross-check.

For each pair, the sixteen interval products bound the actual pairing.
Both cap norms equal two, so a finite Coxeter angle requires pairing
`-2*cos(pi/m)` for an integer `m>=2`.

All labels from 2 through 120 have rigorous closed enclosures. The pi
interval comes from the rational alternating series in Machin's identity,

```
pi = 16*atan(1/5) - 4*atan(1/239).
```

The arctangent sums use 60 and 18 terms, respectively, with the next term
bounds retained. The tangent addition identity is checked rationally.
For each `pi/m` interval, a 48-term fixed-integer cosine calculation at
scale `10^70` includes the alternating tail. Every rounding is outward;
no floating-point trigonometric value is used.

Labels above 120 need no algebraic-degree assumption. Their values lie
between `-2` and `-2*cos(pi/120)`, so the **entire interval** is retained,
using the rigorous upper enclosure of the latter endpoint. All dotted
and parallel pairings are retained as well.

## Certified rejection rules

A pair is rejected only if its whole closed pairing interval is:

1. strictly positive; or
2. above the retained high-label tail and disjoint from every finite-angle
   enclosure for labels 2 through 119.

Intervals touching an allowed enclosure are retained. These rules cover
every possible finite label, as well as parallel and dotted relations.
They certify the following counts:

| Classification | Pairs |
|---|---:|
| Retained possible edges | 3,192,210 |
| Rejected: strictly positive | 509,018 |
| Rejected: finite-angle gap | 918,052 |
| Total pairs | 4,619,280 |

The graph retains 732 more edges than the earlier numerical scout graph.
That difference is not a count of newly certified compatible pairs.

Construction exploits the already checked nonpositive covector intervals
to compute each product with two endpoint multiplications. Replay uses
general four-endpoint interval multiplication for every rejected pair and
checks graph symmetry, the full adjacency digest, and all rejection
counts. This independently checks the optimization in the graph producer.

## Replay and limits

```sh
python scripts/fresh_h9_split_leaf_3034_safe_graph_20261013.py
python scripts/fresh_h9_split_leaf_3034_safe_graph_20261013.py --verify
```

The graph alone proves neither a compatible chamber nor infinite volume.
A retained clique still needs exact pair compatibility and a volume test.
Conversely, unsatisfiability using this supergraph and separately certified
necessary escape-cut clauses can exclude the entire bank. Numerical
clauses alone would not suffice for that conclusion.

The prerequisite normal audit establishes norm two, the H9 seed inertia,
seed compatibility, and a common strict timelike center. This graph adds
only conservative cap-pair exclusions. It makes no completeness claim for
caps outside the literal bank and no exclusion of arbitrary completions
of the seed.
