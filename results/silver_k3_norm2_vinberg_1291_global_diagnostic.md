# Global Vinberg audit of the 1,291-wall silver prefix

The `x=6+3sqrt(2)` shell extends the wall list to this 1,291-wall prefix.  No
compact or finite-volume example is certified by it.  The exact global audit
finds 43,891 elliptic rank-eight edge diagrams, of which 28,898 are
endpoint-deficient.  There is no connected affine diagram or admissible ideal vertex.
Neither Vinberg stopping flag fires, so this prefix is not a new genuinely
hyperbolic example.

As a standalone certificate this is only a finite-prefix failure: later
simple walls could close an edge or remove an apparent escape.  That logical
limitation of the endpoint census remains unchanged.

> **Superseding complete-chamber result.**  The eventual intrinsic norm-two
> chamber is now independently proved infinite-volume by two exact
> complete-chamber symmetries whose product has infinite order.  See the
> [infinite-stabilizer theorem](silver_k3_norm2_infinite_stabilizer.md).  The
> result does not infer infinite volume from the deficient edges recorded
> here; it uses an equality-inclusive global separator bound and arithmetic
> reflection-index argument.

The frozen diagnostic is
`data/candidates/silver_k3_norm2_vinberg_1291_global_diagnostic.json`, with
canonical digest

```text
cbe9681f0a15f13881d32b19946021cdd07f6321adba2af5bdc7f3ab776f53e3
```

Its complete off-diagonal pairing histogram has digest

```text
f541b395c3ebd41329e298eddd5638b7d07af01c6739ed5d0be3381e93eba141
```

The histogram has 185 distinct exact silver values, contains all
`binomial(1291,2)=832695` pairings, and has no invalid finite-angle pair.

## Exact endpoint and parabolic census

```text
elliptic rank-eight edges:    43,891
elliptic rank-nine vertices:   6,091
connected affine diagrams:         0
admissible ideal vertices:          0
invalid finite-angle pairs:         0
```

The finite and finite-plus-ideal endpoint histograms coincide:

| endpoints | rank-eight edges |
| ---: | ---: |
| 0 | 4,065 |
| 1 | 24,833 |
| 2 | 14,993 |

Thus 14,993 current edges satisfy the two-endpoint condition and 28,898 do
not.  The exact connected-affine growth transcript is

```text
822876,18363,31504,35929,25962,15884,7468,1985,6,0.
```

No affine diagram survives, so there is no ideal endpoint and no nontrivial
full zero-neighbourhood to accept or reject.  The hereditary
positive-definite search uses literal production wall order and visits

```text
1284,11744,42112,83534,106799,98835,71690,43891
```

positive prefixes at orders one through eight.

## Logical status

```text
compact stopping available:        no
finite-volume stopping available:  no
new genuinely hyperbolic example:  no
final chamber infinite volume:     not proved by this prefix census
complete norm-two chamber volume:  infinite (separate stabilizer theorem)
```

The full parabolic-neighbourhood rule is used.  The absence of an ideal
vertex in this bounded prefix does not rule out later cusps.  Later walls may
close or cut off current open edges, and reflective norms other than two
remain unexhausted.  This diagnostic is neither a new example nor a
non-existence theorem by itself.  The separate stabilizer theorem supersedes
only its complete norm-two volume status; finite-sidedness and reflective
norms other than two remain open.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_silver_k3_1291_global_vinberg.py
PYTHONPATH=src:. python3 scripts/certify_silver_k3_1291_global_vinberg.py --fresh
PYTHONPATH=src:. pytest -q tests/test_silver_k3_1291_global_vinberg.py
```
