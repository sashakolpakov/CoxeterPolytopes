# TOOT: complete Jacobian--Kummer `H12` five-vertex tree census

Status: **PASS/FREEZE; exact positive `SURFACE -> POLYTOPE` theorem**
(2026-08-10).  Worldwide priority remains open.

## Positive theorem

Every proper five-vertex tree colouring over the frozen 112-wall
Jacobian--Kummer doubling stock has now been matrix-classified.  The three
unlabelled tree shapes give exactly

| tree shape | exact ordinary Coxeter types |
|---|---:|
| path `P5` | 201,747 |
| star `K1,4` | 16,872 |
| degree-three fork | 199,362 |
| **complete union** | **417,981** |

All **417,981** ordinary Coxeter matrices are pairwise nonisomorphic.  Every
one defines a finite-volume noncompact polytope in `H^12`, with reflection
index and volume multiplier five over the seed chamber.  The aggregate facet
histogram is

```text
752: 277;  761: 4,302;  770: 38,103;
779: 151,770;  788: 223,529.
```

This is a complete theorem for the fixed Jacobian--Kummer five-vertex
doubling-stock denominator.  It is not a classification of all hyperbolic
Coxeter polytopes in `H^12`.

## The closed fork denominator

The last open layer was the all-`B` fork.  Its complete result is **106,385**
pairwise nonisomorphic 788-facet types:

| `B` occurrences | exact fork types | facets |
|---:|---:|---:|
| 0 | 116 | 752 |
| 1 | 2,042 | 761 |
| 2 | 18,061 | 770 |
| 3 | 72,758 | 779 |
| 4 | 106,385 | 788 |
| **total** | **199,362** | **752--788** |

For `B=4`, the single `BB|BB` class contains 40,720,800 raw proper
colourings.  The explicit seed-action orbit cover has 106,385 rows and its
orbit-size sum is exactly 40,720,800; the independent Burnside fixed sum is
40,851,840 for the 384-element action, again giving 106,385 orbits.  Thus the
layer is a lossless quotient, not a search sample.

Every `B=4` boundary was reconstructed in the rank-13 divisor lattice.  All
primitive-integral reflection, strict-interior, fixed Smith-`1^13` core,
rank, signed-int64, and ordinary-Coxeter gates pass.  Each boundary deletes
exactly the four marked internal mirrors and has exactly 788 surviving
facets.

## Exact nonisomorphism and cross-shape separation

Inside `B=4`, degree-profile filtering produces 4,559 collision groups with
9,410 members.  Every member was replayed as literal degree-profile bytes;
canonical exact colour refinement separates every literal collision.

The complete same-facet predecessor set is

| source | rows |
|---|---:|
| `P5 BBBB` | 108,159 |
| star `B=4` | 8,985 |
| **total** | **117,144** |

Digest filtering finds 5,186 predecessor candidates.  Literal replay yields
5,536 equal-degree pairs, and exact refinement separates all 5,536.  The
unresolved cross-shape pair count is zero.

Across the complete fork, the three production tranches contain 6,844
literal degree-collision groups with 14,090 members.  The complete
same-header fork-versus-`P5`/star audits contain 8,525 literal equal-degree
pairs; all are refinement-separated.  Different `B` layers have different
facet counts, so they are automatically mutually nonisomorphic.  The
independent complete star--`P5` theorem supplies the remaining cross-shape
edge.  These facts prove the exact 417,981-type union.

## Hostile terminal acceptance

The production source, its exact imports, the RSS watchdog, every input, and
both `B=4` output archives are SHA-pinned at launch.  A separate terminal
acceptance certifier then performs two checks without rebuilding the 106,385
boundaries:

1. It streams colours and degree digests together from all 117,144 archived
   predecessors, regenerates their orbit order, and proves exact row-by-row
   equality with the production prior iterator.  The source split is exactly
   `108,159+8,985`, every source/index key is unique, and every row has four
   `B` colours.  This eliminates positional colour/digest association as a
   provenance weakness.
2. It independently parses and rehashes all 106,385 compact output rows.  The
   final number of refined rows with discrete canonical form is **13,945**,
   not the production summary's diagnostic-only count 9,300; the latter did
   not count cross-only refinements.  The type count and nonisomorphism
   theorem are unchanged.

The acceptance certifier was then rerun in archive-equality mode.  Both
passes terminate green.  The frozen acceptance archive has SHA-256

```text
25ac2f56bd33b41b2d40e04d6659addc6653226cc3c96dab0e034cb6d71cd2e1
```

The combined complete-three-shape theorem archive has SHA-256

```text
d8545b484ef56038043b653964f35f734a6d64d5e4cfb61756ecda9d04b337a5
```

## Finite-volume gate and surface carrier

Allcock's finite-subtree theorem makes every proper colouring a five-copy
reflection fundamental polytope.  Hence every row has finite volume and
index/volume multiplier five.  Finite index in the nonuniform seed
reflection lattice proves noncompactness; no inference is made from an
incomplete ideal-face enumeration.

Every row retains the same weak full-Picard K3 carrier

```text
U(2) + D5(-1)^2 + <-4>
```

on the pinned seven-dimensional K3 family.  This is genuine
`SURFACE -> POLYTOPE` leverage and a large positive factory, but it does not
identify the constructed chambers with full nef cones and does not prove
effectivity or a surface automorphism action.  Worldwide novelty and priority
are not claimed.

## Resources and reproduction

The terminal `B=4` production command was

```bash
python3 scripts/run_with_rss_limit.py --memory-gib 1.0 --timeout 14400 -- \
  python3 scripts/certify_jacobian_kummer_h12_five_vertex_fork_streaming_census.py \
  --min-b 4 --max-b 4
```

It terminated normally in 11,549.659 seconds with sampled direct peak RSS
0.483 GiB.  The first hostile acceptance pass took 399.360 seconds at 0.294
GiB; the independent archive-equality replay took 401.747 seconds at 0.327
GiB.

Verification:

```bash
python3 scripts/certify_jacobian_kummer_h12_five_vertex_fork_b4_terminal_acceptance.py \
  --check-archive
python3 scripts/certify_jacobian_kummer_h12_five_vertex_three_shape_complete_census.py \
  --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_jacobian_kummer_h12_five_vertex_complete_three_shape_census.py \
  tests/test_jacobian_kummer_h12_five_vertex_fork_b3_complete_census.py \
  tests/test_jacobian_kummer_h12_star_p5_b_le_2_exact_separation.py
```

The focused regression set passes **12/12**.

Frozen artifacts:

- [`B=4` production summary](../data/enumerations/jacobian_kummer_h12_five_vertex_fork_b4_to_b4_streaming_census.json)
- [`B=4` compact transcript](../data/enumerations/jacobian_kummer_h12_five_vertex_fork_b4_to_b4_streaming_census_rows.jsonl.gz)
- [`B=4` hostile terminal acceptance certifier](../scripts/certify_jacobian_kummer_h12_five_vertex_fork_b4_terminal_acceptance.py)
- [`B=4` hostile terminal acceptance archive](../data/enumerations/jacobian_kummer_h12_five_vertex_fork_b4_terminal_acceptance.json)
- [`complete three-shape theorem certifier`](../scripts/certify_jacobian_kummer_h12_five_vertex_three_shape_complete_census.py)
- [`complete three-shape theorem archive`](../data/enumerations/jacobian_kummer_h12_five_vertex_three_shape_complete_census.json)
- [`focused tests`](../tests/test_jacobian_kummer_h12_five_vertex_complete_three_shape_census.py)
