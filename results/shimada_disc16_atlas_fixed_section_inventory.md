# Shimada discriminant-16 atlas fixed-section denominator

## Exact result

Shimada's pinned `disc16.txt` atlas contains 4,539 published induced chambers.
Exactly 98 have archived lattice/chamber stabilizer `Aut(D)` of order greater
than one.  Enumerating every subgroup of each of those 98 archived groups and
quotienting by inner conjugacy gives

```text
346 subgroups
257 subgroup conjugacy classes
 98 identity classes, one per chamber
159 proper subgroup classes.
```

Exactly 107 **proper** subgroup classes have fixed rank at least 11, hence a
fixed hyperbolic section above `H9`:

| fixed hyperbolic dimension | proper archived-`Aut(D)` classes |
|---:|---:|
| `H10` | 4 |
| `H11` | 101 |
| `H14` | 2 |
| **total** | **107** |

These 107 classes occur across 94 chambers.  The archived JSON field
`high_dimensional_subgroup_class_count = 205` is also correct but includes the
98 identity classes, whose fixed section is the unchanged ambient `H19`
chamber.  The proper fishing denominator is therefore `205 - 98 = 107`, not
205.

All six proper casts whose dimension is not `H11` are concentrated in two
chambers (indices are the atlas's zero-based indices):

| chamber | walls | archived `Aut(D)` order | proper non-`H11` classes |
|---:|---:|---:|---|
| 4 | 158 | 48 | one `H14`; three `H10` |
| 5 (`D5`) | 165 | 20 | one `H14`; one `H10` |

The remaining 101 casts are in `H11`.  The inventory does not identify casts
across different chambers up to isometry, so 107 is a class count in the
archived chamber-by-chamber denominator, not a claim of 107 distinct abstract
Coxeter polytopes.

## Downstream closure status

The denominator certificate itself remains an inventory, but its entire
107-cast queue has now been processed.  Three casts map exactly to the
previously closed `D5` stock.  Of the remaining 104, 100 fail the exact
Coxeter-angle gate and four `D0004` casts survive.  [Exact normalized
squared-Gram isomorphisms](../data/enumerations/shimada_disc16_d4_duplicate_survivor_identifications.json)
identify `D0004/C001` with the already certified
`D5 H14/60` chamber and both `C003` and `C004` with the certified
`D5 H10/30` chamber.  The distinct order-four `C007` cast is a proved
27-wall finite-volume noncompact `H10` chamber; see its
[exact lineage theorem](shimada_disc16_d4_c7_h10_27_lineage.md).

Consequently the archived lane is now closed downstream as well: 101 failed
cast instances and six finite-volume cast instances representing three
abstract chamber types.  The separate full weighted-wall audit adds 575
fresh rational-symmetry casts, all of which fail the exact angle gate; see the
[575-cast no-go](shimada_disc16_fresh_weighted_high_section_no_go.md).

## Scope boundary

This certificate is complete for exactly the following input:

> the integral lattice/chamber stabilizer matrices `Aut(D)` archived by
> Shimada for all 4,539 chambers in `disc16.txt`.

It is **not** the full weighted-wall graph automorphism census.  In particular:

- it does not assume that archived `Aut(D)` equals the full exact symmetry
  group of each weighted wall configuration;
- the subsequently completed full weighted-wall graph computation is a
  separate certificate and adds 575 subgroup casts;
- the two denominators must not be added together or described as one closed
  full-symmetry census;
- no wall-orbit reduction, Coxeter-angle gate, irredundancy proof, or Vinberg
  volume test is performed here.

Thus the conclusion of this inventory artifact itself is:

> The archived-stabilizer denominator is closed and contains 107 proper casts
> above `H9`.

The downstream paragraph above records the later completed angle and volume
work without retroactively expanding this artifact's proof scope.

## Stabilizer-order audit

Across all 4,539 chambers the archived stabilizer-order histogram is

```text
order  1: 4441 chambers
order  2:   84
order  3:    4
order  4:    2
order  6:    5
order  8:    1
order 20:    1
order 48:    1.
```

For every nontrivial archived group the certificate constructs its faithful
left regular permutation action, enumerates all subgroups, computes inner
conjugacy classes, and then computes fixed ranks from the original integral
20-by-20 lattice matrices.  The regular action is only an exact subgroup-
enumeration device; it is not used to approximate the lattice fixed space.

## Artifacts and replay

- generator:
  [`inventory_shimada_disc16_atlas_fixed_sections.py`](../scripts/inventory_shimada_disc16_atlas_fixed_sections.py)
- certificate:
  [`shimada_disc16_atlas_fixed_section_inventory.json`](../data/enumerations/shimada_disc16_atlas_fixed_section_inventory.json)
- regression test:
  [`test_shimada_disc16_atlas_fixed_section_inventory.py`](../tests/test_shimada_disc16_atlas_fixed_section_inventory.py)

Replay from the repository root:

```bash
PYTHONPATH=src:. python3 scripts/inventory_shimada_disc16_atlas_fixed_sections.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_shimada_disc16_atlas_fixed_section_inventory.py
```

The full archive replay takes materially longer than the JSON regression test
because it re-enumerates every subgroup of all 98 archived nontrivial
stabilizers.
