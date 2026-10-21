# TOOT: Jacobian--Kummer `H12` five-fork `B=3` census

Status: **PASS/FREEZE; exact positive `SURFACE -> POLYTOPE` theorem**
(2026-08-10). Worldwide priority remains open.

## Positive result

The degree-three five-vertex fork over the frozen 112-wall
Jacobian--Kummer doubling stock has exactly **72,758** ordinary Coxeter types
with three `B`-coloured edges. They are pairwise nonisomorphic,
finite-volume noncompact polytopes in `H^12`, of reflection index and volume
multiplier five. Every row has exactly **779 facets**.

Write `(a,d|b,c)` for the central long edge, its continuation, and the
unordered short-arm pair. The exact layer is

| edge-class pattern | raw proper colourings | exact seed-action orbits |
|---|---:|---:|
| `BB|AB` | 13,862,400 | 36,247 |
| `BA|BB` | 6,858,240 | 18,055 |
| `AB|BB` | 7,004,160 | 18,456 |
| **total** | **27,724,800** | **72,758** |

For each row, an explicit orbit cover agrees independently with the Burnside
fixed-configuration sum for the exact 384-element seed action. Thus this is a
lossless census of the complete `B=3` fork layer, not a sample.

Together with the frozen `B<=2` result, the fork now supplies exactly
**92,977** pairwise nonisomorphic finite-volume noncompact `H^12` polytopes:

| `B` occurrences | exact fork types | facets |
|---:|---:|---:|
| 0 | 116 | 752 |
| 1 | 2,042 | 761 |
| 2 | 18,061 | 770 |
| 3 | 72,758 | 779 |
| **total** | **92,977** | **752--779** |

The exact remaining fork denominator is **106,385**, precisely the `B=4`
shard. No `B=4` matrix-census claim is made.

## Exact cross-shape increment

Every `B=3` fork row was compared with every completed `P5` and `K1,4` row
having 779 facets. The complete and disjoint predecessor ledger is

| predecessor source | rows |
|---|---:|
| `P5 ABBB` | 36,247 |
| `P5 BABB` | 36,648 |
| star `B=3` | 6,117 |
| **total** | **79,012** |

Digest filtering found 2,618 predecessor candidates. Literal replay produced
2,690 equal-degree pairs; exact colour refinement separated all 2,690. The
unresolved cross-shape pair count is zero. Hence `B=3` is an exact **+72,758**
increment over the completed `P5` and star sets, not merely a different
marked-tree presentation of old matrices.

One clean exact union now contains

- all 201,747 complete `P5` rows;
- all 15,102 complete star rows with `B>=3`;
- all 92,977 fork rows with `B<=3`.

These **309,826** ordinary Coxeter matrices are pairwise nonisomorphic. The
remaining 1,770 star rows with `B<=2` are deliberately omitted because their
separation from `P5` at the same facet headers remains open.

## Exact matrix and finite-volume gates

For fork colours `(a,d|b,c)`, the five chamber transformations are

```text
1, R_a, R_b, R_c, R_a R_d.
```

The production certifier reconstructs all 72,758 boundaries in the rank-13
Picard lattice, verifies the literal marked adjacency, deletes exactly four
internal mirrors, and checks on every row:

- primitive integral reflection roots;
- a strict timelike summed interior, with minimum strict pairing five;
- the fixed `W1--W60` Smith-`1^13` full-lattice core;
- full rank 13;
- only ordinary labels `m2,m3,m4,m6,infinity`;
- the exact facet formula `752+9B`.

Within the fork layer, degree profiles have 1,968 collision groups containing
4,043 members. All 4,043 were replayed literally and every literal collision
was separated by canonical exact colour refinement. Across all internal and
cross-shape equality candidates, 6,569 current rows carry a refinement
certificate.

Allcock's finite-subtree theorem supplies a five-copy reflection fundamental
polytope for every proper colouring. Therefore every row is finite-volume of
index/volume multiplier five. Finite index in the nonuniform seed lattice
also proves every row noncompact. This conclusion does not come from an
incomplete ideal-face enumeration.

## Hostile post-run acceptance

The launch-time production source and its four direct exact imports are
SHA-pinned, as are every input and both output archives. A separate terminal
acceptance certifier then performs two independent checks without rebuilding
the 72,758 boundaries:

1. It streams the archived colours and degree digests from all 79,012
   predecessors, regenerates their orbit order, and proves row-by-row equality
   with the production prior iterator. This removes the positional
   colour/digest association weakness identified by the static hostile audit.
2. It independently rehashes and parses all 72,758 compact rows. This corrects
   a diagnostic-only production undercount: the exact number of refined rows
   with discrete canonical refinement is **6,446**, not the raw summary's
   3,975. The type count and nonisomorphism theorem are unchanged.

The acceptance replay has its own deterministic transcript hash, unique
`(source,index)` checks, exact source split `36,247+36,648+6,117`, exact
three-`B` check, and archive-equality mode. The frozen acceptance payload has
SHA-256

```text
4a385121b4f0f14484e258223d62e15fae8b78b42c31ef85c0fda71380176d38
```

## Surface carrier and claim boundary

Every row retains the same weak full-Picard K3 carrier

```text
U(2) + D5(-1)^2 + <-4>
```

on the pinned seven-dimensional family. This is a genuine
`SURFACE -> POLYTOPE` construction, but no full-nef, effectivity, or surface
automorphism statement is made. Worldwide novelty and priority are not
claimed.

The following remain explicitly open:

- the 106,385-row fork `B=4` matrix census;
- star versus `P5` separation for the omitted 1,770 `B<=2` star rows;
- full-nef/action promotion and worldwide-priority audit.

## Resource certificate and reproduction

The exact production command was

```bash
python3 scripts/run_with_rss_limit.py --memory-gib 1.0 --timeout 7200 -- \
  python3 scripts/certify_jacobian_kummer_h12_five_vertex_fork_streaming_census.py \
  --min-b 3 --max-b 3
```

It terminated normally in 4,443.676 seconds with sampled peak RSS
**0.457 GiB**. The archive-equality hostile replay terminated normally in
217.919 seconds with peak **0.216 GiB**.

Verification:

```bash
python3 scripts/certify_jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance.py \
  --check-archive
python3 scripts/certify_jacobian_kummer_h12_five_vertex_fork_b_le_3_complete_census.py \
  --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_jacobian_kummer_h12_five_vertex_fork_b3_complete_census.py \
  tests/test_jacobian_kummer_h12_five_vertex_fork_b_le_2_complete_census.py
```

The new and `B<=2` focused tests pass **8/8**. Six additional focused
`P5 ABBB`, `P5 BABB`, and star predecessor tests pass **6/6**.

Frozen artifacts:

- [`B=3` production summary](../data/enumerations/jacobian_kummer_h12_five_vertex_fork_b3_to_b3_streaming_census.json)
- [`B=3` compact transcript](../data/enumerations/jacobian_kummer_h12_five_vertex_fork_b3_to_b3_streaming_census_rows.jsonl.gz)
- [`hostile terminal acceptance certifier`](../scripts/certify_jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance.py)
- [`hostile terminal acceptance archive`](../data/enumerations/jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance.json)
- [`combined B<=3 theorem certifier`](../scripts/certify_jacobian_kummer_h12_five_vertex_fork_b_le_3_complete_census.py)
- [`combined B<=3 theorem archive`](../data/enumerations/jacobian_kummer_h12_five_vertex_fork_b_le_3_complete_census.json)
- [`focused tests`](../tests/test_jacobian_kummer_h12_five_vertex_fork_b3_complete_census.py)
