# TOOT: complete Jacobian–Kummer H12 `P5`, `B<=2` census

Status: **PASS/FREEZE; full production archive replay and focused tests
green** (2026-08-08).

## Positive result

The complete union of the seven five-vertex path shards having at most two
`B`-colored edges contains

```text
20,693 seed-and-path-reversal color orbits.
```

Every representative is rebuilt over the exact rank-13 Picard lattice of the
172-facet Jacobian–Kummer seed.  Exact ordinary-matrix invariants prove that
all 20,693 are pairwise nonisomorphic within this completed layer.  The facet
histogram is

```text
752^142, 761^2091, 770^18460.
```

## Exact disjoint denominator

The frozen 384-element seed action preserves the split

```text
A = W61,...,W76       (16 colors),
B = W77,...,W172      (96 colors).
```

Thus a binary edge-class word, taken modulo reversal, defines an invariant
and pairwise disjoint shard.  The complete production denominator is:

| class word modulo reversal | `B` occurrences | raw unrooted proper colorings | seed-and-reversal orbits | facets |
|---|---:|---:|---:|---:|
| `AAAA` | 0 | 27,000 | 142 | 752 |
| `AAAB` | 1 | 345,600 | 1,007 | 761 |
| `AABA` | 1 | 368,640 | 1,084 | 761 |
| `AABB` | 2 | 2,188,800 | 5,819 | 770 |
| `ABAB` | 2 | 2,359,296 | 6,304 | 770 |
| `ABBA` | 2 | 1,167,360 | 3,256 | 770 |
| `BAAB` | 2 | 1,105,920 | 3,081 | 770 |
| **total** |  | **7,562,616** | **20,693** |  |

For every row two computations agree: an explicit orbit cover with exact
orbit-size sum, and Burnside's lemma on the same 384 permutations, including
the reversal-fixed contribution for palindromic class words.

At the freeze date of this tranche, the three remaining shards formed the
following resumable boundary:

| remaining class word | raw unrooted proper colorings | exact orbit denominator | facets |
|---|---:|---:|---:|
| `ABBB` | 13,862,400 | 36,247 | 779 |
| `BABB` | 14,008,320 | 36,648 | 779 |
| `BBBB` | 41,154,000 | 108,159 | 788 |
| **remaining** | **69,024,720** | **181,054** |  |

Subsequent update: `ABBB` is now fully certified as 36,247 pairwise
nonisomorphic types.  It is no longer open.  The current remaining
color-orbit denominator is therefore `BABB+BBBB=144,807`; this is an orbit
denominator, not yet a certified type count.  See
[`jacobian_kummer_h12_five_vertex_p5_abbb_census.md`](jacobian_kummer_h12_five_vertex_p5_abbb_census.md).

Together the ten rows recover the independent all-`P5` denominator
`201,747`.  No row can overlap another under the seed action because that
action preserves `A` and `B`.

## Direct hyperbolic and lattice gates

For a proper edge word `(a,b,c,d)`, the engine constructs the five chamber
transformations

```text
1, r_a, r_a r_b, r_a r_b r_c, r_a r_b r_c r_d.
```

This allows nonadjacent repeated colors; it does not inherit the four-distinct
color restriction of the first five pilots.  For every representative the
engine checks the actual marked `P5`, exactly four internal mirrors, a strict
timelike interior, primitive integral reflective boundary normals, the fixed
Smith-`1^13` core `W1,...,W60`, and every off-diagonal ordinary Coxeter
relation.  The facet identity is checked independently on every boundary:

```text
f(P5) = 172 + 145*4 + 9B = 752 + 9B.
```

The census contains 19,677 four-distinct-color types, 1,005 types with one
nonadjacent repetition, and 11 with two repetitions.

Signed 64-bit arithmetic is used only after global transform bounds and
dynamic root/interior Gram-product bounds lie below `2^63`.  The largest
chamber-transform entry observed is 79,836; the largest boundary coordinate
is 131,522; and the largest absolute boundary pairing is 130,916.

## Pairwise ordinary-Coxeter nonisomorphism

Seed-color orbits are not assumed to be unmarked Coxeter types.  The exact
multiset of the five ordinary-label degrees at every facet gives 20,029
initial invariant classes.  Its complete residue consists of 629 collision
groups, 1,293 members in total, and maximum group size four.

Every residual matrix is rebuilt and passed through canonically numbered
one-dimensional color refinement.  The resulting invariant separates every
member inside its initial exact class:

```text
20,693 color orbits = 20,693 pairwise nonisomorphic ordinary matrices.
```

Of the 1,293 refined matrices, 1,165 acquire a discrete canonical coloring.
For the remaining 128, the distinct stable quotient inventories suffice
inside their respective initial degree-profile classes; no assertion that a
nondiscrete stable quotient is a complete canonical form is needed.

## Why these are finite-volume noncompact polytopes

Finite volume is **not** inferred from the divisor lattice or from rational
polyhedrality.  The pinned seed has already passed the literal Vinberg and
independent exact-cone finite-volume gates.  Its 112 selected walls have been
checked directly to be doubling walls and pairwise parallel or
ultraparallel.  Allcock's finite-subtree theorem therefore makes every
five-vertex path a Coxeter fundamental polytope of reflection index and
volume multiplier five.  Finite index in the nonuniform seed lattice
preserves nonuniformity, hence every example is finite-volume noncompact.

## Divisor-lattice carrier and claim boundary

Every path remains in the same full Néron–Severi lattice

```text
U(2) + D5(-1)^2 + <-4>
```

of the pinned seven-dimensional K3 family.  The fixed Smith core excludes a
proper-sublattice artifact.  This is a weak full-Picard divisor-reflection
carrier: the result does not claim that a glued chamber is a full nef cone,
that every desired class is effective, or that a surface automorphism
realizes the subgroup.

Three earlier pilots (`B=0,1,2`) are recovered in exact seed-action/path-
reversal orbits.  The completed layer therefore contributes 20,690 types
beyond all five earlier `P5` pilots; the remaining `B=3,4` pilots have
different facet counts.  Both the exact degree-profile SHA-256 screen and the
combined type-invariant SHA-256 screen have zero matches against all 16,872
archived star types.  This is a collision-resistant archival comparison, not
promoted to a pure exact cross-shape nonisomorphism theorem.  Worldwide
priority remains open: Allcock's exponential dimension-12 abundance is
classical, and no publication-complete literature comparison has yet been
made for these diagrams.

## Reproduction

The parameterized certifier can replay the production layer, any remaining
binary shard, or the full ten-shard census:

```bash
python3 scripts/certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py --check-archive
python3 scripts/certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py \
  --pattern ABBB --output data/enumerations/jacobian_kummer_h12_five_vertex_p5_abbb_census.json
python3 scripts/certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py \
  --all-patterns --output data/enumerations/jacobian_kummer_h12_five_vertex_p5_complete_census.json
PYTHONPATH=src:. pytest -q \
  tests/test_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py
```

Fresh verification:

```text
archive_check=PASS
P5_types=20693
type_facet_histogram={"752":142,"761":2091,"770":18460}
focused tests: 6 passed in 15.48s
builder SHA-256: d7c5c02803911f8d205a63c59027d8d9584cce998de5efa0a7acb134d56d00e4
archive SHA-256: 32943e8531bdaa1874af41d45635ea792c9337694adc4f20dfc152aac3d7713a
```

Artifacts:

- [`certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py`](../scripts/certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py)
- [`jacobian_kummer_h12_five_vertex_p5_b_le_2_census.json`](../data/enumerations/jacobian_kummer_h12_five_vertex_p5_b_le_2_census.json)
- [`test_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py`](../tests/test_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py)
