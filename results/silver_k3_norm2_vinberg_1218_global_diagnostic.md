# Global Vinberg audit of the 1,218-wall silver prefix

No compact or finite-volume hyperbolic Coxeter polytope is certified by the
1,218-wall prefix.  The exact global audit finds 42,297 elliptic rank-eight
edges, but 27,854 still lack exactly two finite-or-ideal endpoints.  Neither
Vinberg stopping flag fires, so this prefix is not a new genuinely hyperbolic
example.

This is a negative result for the current incomplete prefix, not an
infinite-volume theorem for the eventual chamber.  Later simple walls may
close an edge or remove the apparent escape represented by it.

The frozen diagnostic is
`data/candidates/silver_k3_norm2_vinberg_1218_global_diagnostic.json`, with
canonical digest

```text
cfb9272463b66ba8295a57b79a80bea4a311f64b6873a987e3098d4a5f1b0dea
```

Its complete off-diagonal pairing histogram has digest

```text
29989ed7befef001484c5d68731b208fbac47538d0e801d48760200e0c7ca99d
```

The histogram has 173 distinct exact silver values, contains all
`binomial(1218,2)=741153` pairings, and has no invalid finite-angle pair.

## Exact endpoint and parabolic census

```text
elliptic rank-eight edges:    42,297
elliptic rank-nine vertices:   5,875
connected affine diagrams:         0
admissible ideal vertices:          0
invalid finite-angle pairs:         0
```

The finite and finite-plus-ideal endpoint histograms coincide:

| endpoints | rank-eight edges |
| ---: | ---: |
| 0 | 3,865 |
| 1 | 23,989 |
| 2 | 14,443 |

Thus 14,443 current edges satisfy the two-endpoint condition and 27,854 do
not.  The exact connected-affine growth transcript, recording the number of
candidate diagrams passed to the next order, is

```text
731828,17502,30365,34816,25215,15458,7310,1949,6,0.
```

No affine diagram survives, so there is no ideal endpoint and no nontrivial
full zero-neighbourhood to accept or reject.

The hereditary positive-definite search uses the literal production wall
order and visits

```text
1211,11054,39647,78884,101295,94467,68937,42297
```

positive prefixes at orders one through eight.

## Logical status

```text
compact stopping available:        no
finite-volume stopping available:  no
new genuinely hyperbolic example:  no
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
PYTHONPATH=src:. python3 scripts/certify_silver_k3_1218_global_vinberg.py
PYTHONPATH=src:. python3 scripts/certify_silver_k3_1218_global_vinberg.py --fresh
PYTHONPATH=src:. pytest -q tests/test_silver_k3_1218_global_vinberg.py
```
