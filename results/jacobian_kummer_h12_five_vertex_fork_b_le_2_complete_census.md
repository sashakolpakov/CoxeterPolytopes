# TOOT: Jacobian--Kummer `H12` five-fork `B<=2` census

Status: **PASS/FREEZE; exact positive `SURFACE -> POLYTOPE` theorem**
(2026-08-10). Worldwide priority remains open.

## Positive result

The degree-three five-vertex fork over the frozen 112-wall
Jacobian--Kummer doubling stock has exactly **20,219** ordinary Coxeter types
with at most two `B`-coloured edges. They are pairwise nonisomorphic,
finite-volume noncompact polytopes in `H^12`, of reflection index and volume
multiplier five:

| `B` occurrences | exact types | facets |
|---:|---:|---:|
| 0 | 116 | 752 |
| 1 | 2,042 | 761 |
| 2 | 18,061 | 770 |
| **total** | **20,219** | **752--770** |

Every fork row was also compared with every completed `P5` and `K1,4` row
having the same facet count. There are no cross-shape isomorphisms. Thus this
is an exact **+20,219** increment over both completed shape sets, not merely a
new marked-tree presentation of old matrices.

One especially clean exact union now contains:

- all 201,747 complete `P5` rows;
- all 15,102 complete star rows with `B>=3`;
- all 20,219 fork rows with `B<=2`.

These **237,068** matrices are pairwise nonisomorphic. The smaller-header
star--`P5` comparison remains open, so no 238,838-type all-three-shape union
is claimed.

## Exact full-fork denominator and terminating remainder

The aggregate 199,362-orbit fork denominator is now refined into all twelve
invariant edge-role patterns. Write `(a,d|bc)` for the central long edge,
its extension, and the unordered short-arm pair.

| `B` count | pattern orbit counts | total |
|---:|---|---:|
| 0 | `AA|AA`: 116 | 116 |
| 1 | `BA|AA`: 556; `AB|AA`: 479; `AA|AB`: 1,007 | 2,042 |
| 2 | `BB|AA`: 2,941; `BA|AB`: 6,220; `AB|AB`: 5,896; `AA|BB`: 3,004 | 18,061 |
| 3 | `BB|AB`: 36,247; `BA|BB`: 18,055; `AB|BB`: 18,456 | 72,758 |
| 4 | `BB|BB`: 106,385 | 106,385 |
| **all** |  | **199,362** |

For every row of the table, an explicit orbit cover agrees with an independent
Burnside fixed-configuration sum for the exact 384-element seed action. The
raw proper-colouring counts sum to 75,897,360, recovering the earlier
aggregate denominator.

The complete fork matrix census therefore has an exact terminating remainder:
**179,143 rows**, precisely the `B=3` and `B=4` shards. This is a matrix-census
remainder, not an UNKNOWN finite-volume gate: Allcock's finite-subtree theorem
already supplies the polytope and finite-volume conclusion for every proper
colouring.

## Direct Coxeter and nonisomorphism gates

For fork colours `(a,d|b,c)`, the five chamber transformations are

```text
1, R_a, R_b, R_c, R_a R_d.
```

The certifier verifies the literal marked adjacency, removes exactly four
internal mirrors, reconstructs every exterior normal in the rank-13 Picard
lattice, and checks on every row:

- primitive integral reflection roots;
- a strict timelike summed interior;
- the fixed `W1--W60` Smith-`1^13` full-lattice core;
- full rank 13;
- only ordinary Coxeter labels `m2,m3,m4,m6,infinity`;
- the exact facet formula `752+9B`.

Within the 20,219 fork rows, exact degree profiles leave 317 collision groups
with 637 members. Every member was rebuilt from its four colours and all
literal degree collisions were separated by canonical colour refinement.

Against the 20,042 same-facet completed `P5`/star rows in the `B=2` layer,
digest filtering finds 291 prior candidates. Literal replay produces 299
equal-degree pairs; exact refinement separates all 299. The `B<=1` layer has
no equal-degree cross-shape candidate among its 2,421 comparison rows.
Consequently the unresolved pair count is zero. Different digests are used
only in the logically safe direction: different digest implies different
literal bytes; equality candidates are never trusted without replay.

## Finite volume and surface carrier

The frozen 172-facet seed is directly certified finite-volume noncompact, and
all 112 selected colours are pairwise disjoint doubling walls. A proper
five-vertex subtree therefore gives an index-five reflection fundamental
polytope and volume multiplier five. Finite index in the nonuniform seed
lattice preserves nonuniformity.

Every row retains the same weak full-Picard K3 carrier

```text
U(2) + D5(-1)^2 + <-4>
```

on the pinned seven-dimensional family. No full-nef, effectivity, or surface
automorphism claim is made.

## Resource certificate and reproduction

The successful `B=2` production run was guarded by a hard 1 GiB ceiling:

```text
python3 scripts/run_with_rss_limit.py --memory-gib 1.0 --timeout 900 -- \
  python3 scripts/certify_jacobian_kummer_h12_five_vertex_fork_streaming_census.py \
  --min-b 2 --max-b 2
```

It completed in 594.390 seconds with sampled peak RSS **0.332 GiB**. The
`B<=1` rebuild completed in 52.011 seconds with peak **0.231 GiB**. Boundaries
are processed one at a time; only compact row digests and exact replay inputs
are retained.

Focused verification:

```bash
PYTHONPATH=src:. pytest -q \
  tests/test_jacobian_kummer_h12_five_vertex_fork_b_le_2_complete_census.py
python3 scripts/certify_jacobian_kummer_h12_five_vertex_fork_b_le_2_complete_census.py \
  --check-archive
```

Frozen artifacts:

- [`fork streaming certifier`](../scripts/certify_jacobian_kummer_h12_five_vertex_fork_streaming_census.py)
- [`combined theorem certifier`](../scripts/certify_jacobian_kummer_h12_five_vertex_fork_b_le_2_complete_census.py)
- [`twelve-pattern denominator`](../data/enumerations/jacobian_kummer_h12_five_vertex_fork_pattern_denominator.json)
- [`combined B<=2 theorem archive`](../data/enumerations/jacobian_kummer_h12_five_vertex_fork_b_le_2_complete_census.json)
- [`focused tests`](../tests/test_jacobian_kummer_h12_five_vertex_fork_b_le_2_complete_census.py)

The `B=3` and `B=4` production runs, complete fork census, full-nef/action
promotion, and worldwide-priority audit remain open.
