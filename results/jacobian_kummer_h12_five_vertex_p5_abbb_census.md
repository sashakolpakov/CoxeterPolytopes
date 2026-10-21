# TOOT candidate: complete Jacobian–Kummer H12 `P5` `ABBB` census

Status: **PASS/FREEZE; corrected metadata, full final-source archive replay,
and focused tests green** (2026-08-08).

## Target positive result

The next exact binary class-word shard after the frozen `B<=2` theorem is

```text
ABBB (identified with BBBA by path reversal).
```

Its complete denominator is

```text
13,862,400 raw unrooted proper colorings
       -> 36,247 seed-and-path-reversal orbits.
```

Every resulting path has three `B`-edge occurrences and hence 779 facets by

```text
f(P5) = 172 + 145*4 + 9*3 = 779.
```

The production run rebuilds every orbit representative over the exact
rank-13 Picard lattice, rather than extrapolating from the single earlier
`B=3` pilot.  All 36,247 representatives pass, and exact ordinary-matrix
invariants prove that they are pairwise nonisomorphic within this shard.

## Exact denominator and separation from the frozen tranche

For the canonical nonpalindromic class word `ABBB`, the raw count is

```text
16 * 96 * 95 * 95 = 13,862,400.
```

Path reversal identifies this word with `BBBA`; fixing the `ABBB`
orientation therefore leaves the 384-element seed action as the exact
quotient group.  Explicit orbit cover and Burnside independently give
36,247 orbits.

The frozen `B<=2` layer has facet counts 752, 761, and 770.  Every `ABBB`
boundary has 779 facets, so the two completed tranches are exactly disjoint
as ordinary Coxeter matrices without relying on hashes.  The combined exact
completed path count is therefore

```text
20,693 + 36,247 = 56,940.
```

The two remaining binary shards will then be:

| remaining class word | raw unrooted colorings | exact orbit denominator | facets |
|---|---:|---:|---:|
| `BABB` | 14,008,320 | 36,648 | 779 |
| `BBBB` | 41,154,000 | 108,159 | 788 |
| **remaining** | **55,162,320** | **144,807** |  |

Because `BABB` has the same 779-facet header as `ABBB`, its eventual merge
will require an exact cross-shard ordinary-matrix comparison.  `BBBB` is
header-separated at 788 facets.

## Direct gates and finite volume

The separate certifier imports the now-immutable `B<=2` boundary engine and
pins both its builder and archive SHA-256 digests.  For every representative
it constructs

```text
1, r_a, r_a r_b, r_a r_b r_c, r_a r_b r_c r_d
```

and verifies the marked path, exactly four internal mirrors, a strict
timelike interior, primitive integral reflective boundary normals, the
Smith-`1^13` Picard core, and all ordinary Coxeter relations.  Exact degree
profiles and canonical refinement are used to separate unmarked ordinary
matrices; color orbits are not assumed to be Coxeter types.

Of the 36,247 types, 35,846 use four distinct colors and 401 repeat one
nonadjacent `B` color.  The largest chamber-transform entry observed is
301,368, the largest boundary coordinate is 462,070, and the largest
absolute boundary pairing is 476,144; every dynamic pre-multiplication bound
remains below `2^63`.

## Pairwise ordinary-Coxeter nonisomorphism

The exact multiset of five ordinary-label degrees at every facet gives
35,202 initial invariant classes.  The complete residue consists of 1,023
collision groups, 2,068 members, and maximum group size four.  Every one of
those 2,068 ordinary matrices is rebuilt and canonical color refinement is
discrete in every case.  The resulting exact canonical words are distinct
inside each initial class:

```text
36,247 color orbits = 36,247 pairwise nonisomorphic ordinary matrices.
```

Finite volume is transported by Allcock's finite-subtree theorem from the
separately certified nonuniform 172-facet seed.  Thus every successful row is
a finite-volume noncompact `H^12` Coxeter polytope of reflection index and
volume multiplier five.  No finite-volume claim is inferred merely from the
divisor lattice or rational polyhedrality.

## Carrier and novelty boundary

All examples remain in the weak full-Picard K3 carrier

```text
U(2) + D5(-1)^2 + <-4>.
```

No full-nef, effectivity, or surface-automorphism realization is claimed.
The original `B=3` pilot is recovered inside this shard.  The theorem gives
36,246 further `ABBB` types beyond it and 56,936 types beyond all five
original path pilots after combining with the frozen tranche.

Comparison against the 16,872-type star archive is only a collision-resistant
SHA-256 screen.  Both the exact degree-profile and combined-type screens have
zero matches, but this is not promoted to a pure exact cross-shape theorem.
Worldwide priority remains open.

## Reproduction

```bash
python3 scripts/certify_jacobian_kummer_h12_five_vertex_p5_abbb_census.py --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_jacobian_kummer_h12_five_vertex_p5_abbb_census.py
```

Final production verification:

```text
archive_check=PASS
all_checks_pass=true
ABBB_types=36247
combined_completed_P5_types=56940
metadata_semantics=PASS
focused tests: 6 passed in 12.11s
builder SHA-256: 6217c4902b8a66bbb4cef57e4a2edc188da723d6ceffcc2cad82cd8d172eb18e
archive SHA-256: b535bbe5d1639358e5a336ce606a4d7eae8dedd0b55835f9ee057321ac394299
```

Artifacts:

- [`certify_jacobian_kummer_h12_five_vertex_p5_abbb_census.py`](../scripts/certify_jacobian_kummer_h12_five_vertex_p5_abbb_census.py)
- [`jacobian_kummer_h12_five_vertex_p5_abbb_census.json`](../data/enumerations/jacobian_kummer_h12_five_vertex_p5_abbb_census.json)
- [`test_jacobian_kummer_h12_five_vertex_p5_abbb_census.py`](../tests/test_jacobian_kummer_h12_five_vertex_p5_abbb_census.py)
