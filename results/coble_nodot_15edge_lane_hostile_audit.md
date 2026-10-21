# Hostile audit of the 15-edge no-parallel Coble lane

**Status (2026-08-09): PASS.  No mathematical or scope defect found.**

This audit is independent of the production lane's enumeration, cone,
overlattice, Vinberg, and graph-isomorphism helpers.  Its executable and
archive are
[`audit_coble_nodot_15edge_lane_hostile.py`](../scripts/audit_coble_nodot_15edge_lane_hostile.py)
and
[`coble_nodot_15edge_lane_hostile_audit.json`](../data/enumerations/coble_nodot_15edge_lane_hostile_audit.json).

## Independent denominator reconstruction

A direct call to nauty `geng`, followed by exact integer rank and rational
congruence rather than the production numerical spectral prefilter, gives

```text
5,741 connected 12-vertex, 15-edge, maximum-degree-three supports
rank 12: 5,588
rank 11:   151
rank 10:     2
```

Both rank-ten matrices have exact inertia `(9,1,2)`.  Their graph6 strings
are exactly

```text
K??E@_KK_YWC
K?AA@_ghECHG.
```

Their Smith diagonals are respectively `1^8,2,2,0,0` and
`1^10,0,0`, so the saturated radical quotients have determinants `-4` and
`-1`.  Direct binary half-coset enumeration—not the general production
discriminant-group routine—finds exactly two isotropic order-two classes in
the determinant-four case and none in the unimodular case.  Hence the exact
branch count is

```text
2 root Grams -> 2+1 = 3 integral E10 branches.
```

This exact rank histogram also rules out the main possible numerical-filter
failure: there is no additional rank-ten support close enough to zero to
have been discarded by floating-point inertia.

## Independent cone and return checks

Direct `cddexec_gmp` calls with a separate parser give the following counts
per `E10` branch:

| quotient | branch | finite | ideal | hyperideal |
|---|---:|---:|---:|---:|
| `det=-4` | 0 | 38 | 17 | 0 |
| `det=-4` | 1 | 38 | 17 | 0 |
| `det=-1` | 0 | 18 | 8 | 24 |

All twelve input inequalities are irredundant in every branch.  The complete
24-ray hyperideal set in the rejected branch has independent primitive-list
SHA-256

```text
0c5fb0b79f03096f5d51b55824c62295d47fb9be751aec0932077d85a61897cc.
```

Thus its eight ideal rays do not imply finite volume: the 24 spacelike rays
give the claimed infinite-volume obstruction.

Nauty canonicalization identifies only the determinant-four root graph with
the archived Type-II graph.  For each of its two integral extensions, the
audit independently lifts all 17 cusps, adds the characteristic boundary,
and obtains

```text
30 facets; 1,333 finite + 207 ideal + 0 hyperideal rays.
```

It then constructs an explicit 30-facet permutation and a marked integral
unimodular isometry to the frozen Type-II lattice and walls.  The two
isometry determinants are `-1` and `+1`; all 900 wall-Gram entries agree,
all walls map, and both wall stocks span `I_(1,10)` with index one.  The
permutation hashes are

```text
6a3649b1560f23e34fc9483233b2abaab1ead2cd50e86c40ccbf089f4d50484a
623928b5df5e1f492f37d1f80bafbc2814b046e93c1723d89e9aaf2b61514735.
```

This is stronger than an unweighted graph recognition.  Exact marked
identity transports the independently frozen Type-II Vinberg census, while
the production lane also recomputes that census directly.

## Scope verdict

The production report is conservative and its L4 disclaimer is correct.
The result is a genuine scoped, surface-free inverse calibration because the
lane is declared before the survivors and every finite-volume survivor
returns through the already exhaustive `raw_30` compiler.  It is not strict
L4 for the full nodot class because it omits:

- every ordinary support of maximum degree above three;
- every 15-edge diagram with pairing-two edges;
- all other edge counts not covered by the separate lanes.

Moreover, the determinant-one branch is eliminated by cdd before algebraic
geometry enters.  This is therefore not an AG-derived nonexistence theorem
for abstract polytopes.  The sole returned full chamber is the already known
Type-II `H10/30` chamber, so there is no new abstract-polytope claim.  The
split-K3/Coble existence theorem is imported only after literal full-Gram
identity with `raw_30`, exactly as the production archive states.

## Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_nodot_15edge_lane.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_nodot_15edge_lane.py
PYTHONPATH=src:. python3 scripts/audit_coble_nodot_15edge_lane_hostile.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_nodot_15edge_lane_hostile_audit.py
```
