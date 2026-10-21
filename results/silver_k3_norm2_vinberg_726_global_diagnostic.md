# Global Vinberg audit of the 726-wall silver prefix

No new compact or finite-volume hyperbolic Coxeter polytope is certified by
the 726-wall prefix.  The exact global audit finds many closed edges, but
16,760 of its 25,764 elliptic rank-eight edges still do not have exactly two
finite-or-ideal endpoints.  Thus neither Vinberg stopping flag fires.

This is a negative result for the current prefix, not an infinite-volume
theorem for the eventual chamber.  Later simple walls may close an edge or
cut off the apparent escape represented by it.

## Exact census

The frozen prefix consists of 726 norm-two walls in the silver rank-ten
carrier of hyperbolic dimension nine.  Hereditary positive-definite
backtracking gives

```text
elliptic rank-eight edges:       25,764
elliptic rank-nine vertices:      3,604
connected affine subdiagrams:         0
admissible ideal vertices:            0
invalid finite-angle pairs:           0
```

The complete endpoint histogram is

| finite + ideal endpoints | rank-eight edges |
| ---: | ---: |
| 0 | 2,332 |
| 1 | 14,428 |
| 2 | 9,004 |

Because there are no affine subdiagrams in this prefix, the ideal-endpoint
contribution is zero and the combined histogram equals the finite-endpoint
histogram.  Exactly 9,004 edges currently satisfy the two-endpoint
condition.  The other 16,760 do not.

The positive-definite search visits, at orders one through eight,

```text
719, 6,429, 23,154, 46,687, 61,224, 57,746, 43,123, 25,764.
```

These are exact search counts, not samples.

## Finite-volume criterion used

For a finite-volume Coxeter polytope in hyperbolic dimension nine, every
elliptic rank-eight edge diagram must have exactly two endpoints.  An
endpoint may be either an elliptic rank-nine finite vertex or a full
parabolic ideal vertex.

The ideal case is not inferred merely from finding a parabolic subset.  For
every connected affine diagram `J`, the audit forms its entire zero
neighbourhood

```text
N(J) = J union {walls orthogonal to every wall of J}
```

and requires `N(J)` to decompose completely into affine components of total
Coxeter rank eight.  If there are `c` components, this is the size condition
`|N(J)| = 8+c` together with the exact parabolic-rank test.  The resulting
full neighbourhood is counted once as an ideal vertex.  An incident edge is
obtained by omitting exactly one node from each affine component.  This is
the reducible-cusp incidence rule; it prevents counting components or
arbitrary parabolic subsets as separate endpoints.

The implementation avoids scanning all eight-element subsets of 726 walls.
It grows only hereditary positive-definite diagrams using compatibility
bitsets, but retains exact arithmetic for every residual classification.
Connected affine diagrams are grown from connected elliptic diagrams, and
their full zero-neighbourhoods are obtained by bitset intersection.  Small
regressions compare the sparse implementation with the generic exact
Vinberg verifier, including an ideal triangle, a contaminated
zero-neighbourhood, and a reducible cusp.

## Logical status

```text
compact stopping available:        no
finite-volume stopping available:  no
final chamber infinite volume:     not proved
```

The failure has only prefix scope.  In particular:

- the absence of affine diagrams now does not rule out cusps after later
  walls enter;
- an edge with zero or one current endpoint can acquire another endpoint or
  be removed from the final chamber by later walls;
- reflective vectors of norms other than two have not been exhausted.

Accordingly, this audit supplies neither a new example nor a non-existence
theorem.  It does supply an exact quantitative target for the next dual
step: 16,760 current edge defects must be resolved before this carrier can
yield a finite-volume example through the norm-two Vinberg track.

## Reproduction

The frozen literal transcript is
`data/candidates/silver_k3_norm2_vinberg_726_global_diagnostic.json`.

```bash
PYTHONPATH=src:. python3 scripts/certify_silver_k3_726_global_vinberg.py
PYTHONPATH=src:. python3 scripts/certify_silver_k3_726_global_vinberg.py --fresh
PYTHONPATH=src:. pytest -q tests/test_silver_sparse_vinberg_volume.py tests/test_silver_k3_726_global_vinberg.py
```
