# Independent replay of the 726-wall silver diagnostic

The independent replay confirms the negative prefix result: the present
726-wall list is neither a compact example nor a certified finite-volume
example.  It has 16,760 elliptic rank-eight edges whose current
finite-or-ideal endpoint count is not two.

This remains a statement about an incomplete simple-root prefix.  It does
**not** prove that the eventual chamber has infinite volume: later walls may
close a current edge or remove it from the final chamber.

## Independent exact checks

The replay imports neither the production sparse Vinberg verifier nor either
production 726-wall diagnose/certify script.  Starting from the 726 literal
roots, it independently rebuilds the rank-ten wall Gram matrix.  Its canonical
Gram digest is exactly the manifest value.

The independent diagram search finds the following positive principal-subset
counts:

| order | elliptic subsets |
| ---: | ---: |
| 1 | 726 |
| 2 | 9,098 |
| 3 | 42,057 |
| 4 | 103,017 |
| 5 | 151,756 |
| 6 | 140,105 |
| 7 | 79,842 |
| 8 | 25,764 |
| 9 | 3,604 |

The 25,764 rank-eight tuples and 3,604 rank-nine tuples agree with the
production transcript as complete ordered literal sets, not merely as
counts.  The replay then bypasses the finite-type shortcut and applies direct
exact `silver_inertia` elimination to every one of those 29,368 archived
subgrams.  All are positive definite.

Deleting one node from every rank-nine vertex reconstructs the production
endpoint rows exactly.  The double count is

```text
finite edge-vertex incidences = 32,436 = 9 * 3,604.
```

The resulting finite-endpoint histogram is

| endpoints | rank-eight edges |
| ---: | ---: |
| 0 | 2,332 |
| 1 | 14,428 |
| 2 | 9,004 |

## Affine completeness

Every compatible two-node diagram has pairing `0`, `-1`, or `-sqrt(2)`.
Finite prefixes are therefore decided by the exact finite Coxeter
classification.  The replay grows every possible proper elliptic prefix
through order eleven.  This bound is sufficient: an irreducible affine
diagram on `k` nodes has rank `k-1`, while the ambient wall Gram has rank at
most ten.

Rank-two affine pairs are scanned separately.  At higher orders, every proper
subdiagram of an irreducible affine diagram is elliptic, so compatibility and
elliptic-prefix pruning cannot discard one.

For each connected non-finite boundary diagram, nonzero determinant is
certified in at least one of

```text
F_5[s]/(s^2-2),  F_13[s]/(s^2-2).
```

Both are fields because two is a quadratic nonresidue modulo 5 and 13.  A
nonzero reduction is a rigorous witness that the determinant in
`Z[sqrt(2)]` is nonzero.  The replay records 39,043 such witnesses; no
candidate survives to the direct exact-inertia fallback.  Together with the
rank-two scan, this gives zero connected affine diagrams through order eleven.
Consequently there are no ideal endpoints in this prefix, and Vinberg's full
zero-neighbourhood test is vacuous here.

For additional diagnostics only, numerical eigenvalues had a smallest
accepted minimum eigenvalue of `0.010956209263` and a closest rejected minimum
eigenvalue of `-0.193993277713`, with no sign disagreements.  These floating
values play no role in any exact decision.

## Frozen result and reproduction

The compact replay certificate is
`data/candidates/silver_k3_norm2_vinberg_726_independent_replay.json`.
Its certificate digest is
`8e7e0c39cbf6ea2248efdfd7224def31ef199e9cbeb61a95c6f92420bc959a65`.

```bash
PYTHONPATH=src:. python3 scripts/certify_silver_k3_726_independent_replay.py --verify
PYTHONPATH=src:. python3 scripts/certify_silver_k3_726_independent_replay.py --write
PYTHONPATH=src:. pytest -q tests/test_silver_k3_726_independent_replay.py
```

The `--write` replay is the full calculation, including all 29,368 direct
exact inertia checks.  `--verify` and the focused regression test validate the
compact archive, its linked literal transcript, both literal-set digests, and
the 32,436-incidence identity.
