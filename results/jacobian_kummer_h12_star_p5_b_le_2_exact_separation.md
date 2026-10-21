# TOOT: complete low-`B` star--`P5` exact separation in `H12`

Status: **PASS/FREEZE; exact positive cross-shape theorem** (2026-08-10).
Worldwide priority remains open.

## Result

The previously omitted 1,770 five-vertex star rows with `B<=2` are exactly
separated from all 20,693 `P5` rows with the same three facet headers:

| `B` | star rows | `P5` rows | facets |
|---:|---:|---:|---:|
| 0 | 19 | 142 | 752 |
| 1 | 169 | 2,091 | 761 |
| 2 | 1,582 | 18,460 | 770 |
| **total** | **1,770** | **20,693** |  |

Every one of the 22,463 boundaries was rebuilt from its archived colour word.
For every facet, the certifier forms its five ordinary-label degrees, sorts
the complete facet-profile multiset, and serializes the actual bytes.  The
1,770 star byte keys are retained and every rebuilt `P5` key is tested by
literal byte membership.  The exact result is

```text
literal_cross_shape_degree_match_pair_count = 0.
```

Thus no collision-resistance assumption enters the cross-shape theorem.
SHA-256 only pins source provenance and checks each rebuilt row against its
frozen archive entry.

Together with the already frozen high-`B` comparisons, all **218,619**
complete star and complete `P5` matrices are pairwise nonisomorphic.  The
fork `B<=3` theorem is already separated from both shapes at every matching
facet header, so the current exact three-shape union is

```text
201,747 P5 + 16,872 star + 92,977 fork B<=3 = 311,596.
```

Every row is a finite-volume noncompact `H12` polytope of reflection
index/volume multiplier five.  Only the exact 106,385-row fork `B=4` matrix
census remains in the five-vertex denominator.

## Provenance and resource bound

The production archive pins the frozen literal-byte star--`ABBB` comparator,
the streamed archive-colour reader, the P5/star boundary builders, the exact
colour-refinement builder, both source archives, the complete-P5 and high-B
separation theorems, and the fork `B<=3` theorem.  It also records its own
builder hash.

The guarded run terminated normally in 787.549 seconds with sampled peak RSS
**0.163 GiB**, below its 0.6 GiB ceiling.  It retained only the 1,770 star
keys and one boundary/label matrix at a time.

Frozen hashes:

```text
builder  5c13641fca1eef719c917bd1c1b6b62a97176463b61c6a3c1056b2e25af7112e
archive  b73e85b9cd3fa61fc727724dca966e906af33feadcaa5b656c758f6640dee15f
```

Focused verification:

```bash
PYTHONPATH=src:. pytest -q \
  tests/test_jacobian_kummer_h12_star_p5_b_le_2_exact_separation.py
```

Artifacts:

- [`exact comparator`](../scripts/certify_jacobian_kummer_h12_star_p5_b_le_2_exact_separation.py)
- [`exact archive`](../data/enumerations/jacobian_kummer_h12_star_p5_b_le_2_exact_separation.json)
- [`focused tests`](../tests/test_jacobian_kummer_h12_star_p5_b_le_2_exact_separation.py)

The common carrier remains the weak full-Picard K3 divisor lattice
`U(2)+D5(-1)^2+<-4>`.  No full-nef, effectivity, surface-action, or worldwide
novelty claim is made.
