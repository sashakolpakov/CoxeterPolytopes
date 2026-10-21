# Global Vinberg audit of the 1,076-wall silver prefix

No compact or finite-volume hyperbolic Coxeter polytope is certified by the
1,076-wall prefix.  The exact global audit finds 38,564 elliptic rank-eight
edges, but 25,406 still lack exactly two finite-or-ideal endpoints.  Neither
Vinberg stopping flag fires.

This is a negative result for the current incomplete prefix, not an
infinite-volume theorem for the eventual chamber.  Later simple walls may
close an edge or remove the apparent escape represented by it.

The frozen diagnostic is
`data/candidates/silver_k3_norm2_vinberg_1076_global_diagnostic.json`, with
canonical digest

```text
16bb5292d7c5bd2ed7790989031d38125ac38b85a1fa2f3761bc51043cab890c
```

Its complete off-diagonal pairing histogram has digest

```text
0e8c890c2694456414aed1954ea058ef5f3a9780f11177cd3ccd4c9fdf1636f6
```

The histogram has 163 distinct exact silver values, contains all
`binomial(1076,2)=578350` pairings, and has no invalid finite-angle pair.

## Exact endpoint and parabolic census

```text
elliptic rank-eight edges:    38,564
elliptic rank-nine vertices:   5,352
connected affine diagrams:         0
admissible ideal vertices:          0
invalid finite-angle pairs:         0
```

The finite and finite-plus-ideal endpoint histograms coincide:

| endpoints | rank-eight edges |
| ---: | ---: |
| 0 | 3,554 |
| 1 | 21,852 |
| 2 | 13,158 |

Thus 13,158 current edges satisfy the two-endpoint condition and 25,406 do
not.  The exact connected-affine growth transcript, recording the number of
candidate diagrams passed to the next order, is

```text
570133,16193,28460,32822,23870,14597,6887,1827,6,0.
```

No affine diagram survives, so there is no ideal endpoint and no nontrivial
full zero-neighbourhood to accept or reject.

The hereditary positive-definite search uses the literal production wall
order and visits

```text
1069,9832,35327,70512,91253,85977,63734,38564
```

positive prefixes at orders one through eight.  A lexicographically ordered
independent root set gives slightly different intermediate visit counts, as
expected for an order-dependent backtracking trace, but the invariant edge
set and all endpoint totals agree exactly.

## Logical status

```text
compact stopping available:        no
finite-volume stopping available:  no
final chamber infinite volume:     not proved
```

The full parabolic-neighbourhood rule is used: an ideal endpoint would be
counted only when the entire zero-neighbourhood decomposes into affine
components of total Coxeter rank eight.  The absence of one in this bounded
prefix does not rule out later cusps.  Likewise, later walls may close or cut
off current open edges, and reflective norms other than two remain
unexhausted.  This diagnostic is neither a new example nor a non-existence
theorem.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_silver_k3_1076_global_vinberg.py
PYTHONPATH=src:. python3 scripts/certify_silver_k3_1076_global_vinberg.py --fresh
PYTHONPATH=src:. pytest -q tests/test_silver_k3_1076_global_vinberg.py
```
