# Full weighted-wall symmetry denominator for the discriminant-16 atlas

## Result

The exact symmetry denominator is now complete for every published induced
chamber in Shimada's pinned `disc16.txt` atlas.  There are 4,539 chambers.
Their full oriented weighted-wall groups have order distribution

| group order | chambers |
|---:|---:|
| 1 | 3,866 |
| 2 | 659 |
| 3 | 4 |
| 4 | 2 |
| 6 | 5 |
| 8 | 1 |
| 20 | 1 |
| 48 | 1 |

Thus 673 chambers have nontrivial weighted-wall symmetry.  In 3,964 chambers
the exact weighted-wall group equals Shimada's archived integral
stabilizer `Aut(D)`.  In the remaining 575 chambers it is strictly larger.

The enlargement is exceptionally rigid:

- every one of the 575 archived stabilizers is trivial;
- every full weighted-wall group is `C2`;
- every extra involution lifts uniquely to an exact rational Lorentzian
  isometry of the complete oriented wall configuration;
- none of the 575 lifted primal actions is integral on the published
  Neron--Severi lattice.

After enumerating every subgroup and dividing by inner conjugacy, the full
weighted-wall denominator has 1,496 subgroups in 1,407 conjugacy classes.
Exactly 682 **proper** classes have fixed rank at least 11:

| hyperbolic dimension | proper classes |
|---:|---:|
| 10 | 4 |
| 11 | 101 |
| 13 | 1 |
| 14 | 2 |
| 18 | 574 |
| **total** | **682** |

The first 107 rows are precisely the previously frozen archived-`Aut(D)`
queue.  The full weighted-wall audit therefore adds exactly 575 casts:
574 in `H18` and one in `H13`.

This artifact is a denominator theorem, not by itself a polytope theorem.
Its downstream exact screen is now complete, however: every one of the 575
fresh casts has a forbidden pair among independently certified genuine
facets.  There are zero Coxeter survivors and hence zero volume jobs in this
strict-extension branch; see the
[complete 575-cast no-go](shimada_disc16_fresh_weighted_high_section_no_go.md).

## Why the extra symmetries are exact

For each chamber, all primitive oriented wall covectors are read from the
official atlas.  Their full rational Gram matrix supplies an exact coloured
graph: wall norms colour vertices and every nonzero pairing colours an edge;
zero pairings are its nonedges.  Exact colour refinement is followed by a
complete graph-automorphism enumeration.

The walls span the nondegenerate rank-20 dual space.  Consequently a
pairing-preserving wall permutation has at most one linear lift.  The
certificate selects 20 independent walls, solves for that lift over `Q`, and
then verifies all of the following entrywise:

1. every published wall is sent to its prescribed oriented image;
2. the dual Lorentzian Gram matrix is preserved;
3. the conjugate primal action preserves the Neron--Severi Gram matrix; and
4. the full finite permutation group and every subgroup are closed exactly.

For groups equal to archived `Aut(D)`, the already frozen integral matrix
inventory is reused byte-for-byte after exact equality of wall permutations
is proved.  Fresh subgroup enumeration is performed only for the 575 strict
extensions.  This avoids conflating a surface/lattice automorphism with a
rational symmetry of the divisor-wall configuration.

The accelerated colour-refinement enumerator has also been compared
set-for-set with the earlier generic exact graph matcher on both a small
order-six archived group and a newly detected rational `C2` extension.

## Algebraic meaning and limitation

Although the extra involutions do not preserve the integral Neron--Severi
lattice, they preserve the complete divisor-wall configuration over `Q`.
Their fixed spaces are therefore rational Lorentzian sections, and their
intersections with the Neron--Severi lattice are saturated integral lattices.
This keeps the construction anchored in an algebraic divisor stock without
pretending that every folding is an automorphism of the surface.

The nonintegrality matters for lineage: a positive section may have a
reflection-stable mirror lattice different from its saturated divisor-fixed
lattice.  Those two lattices must be separated explicitly before any
arithmetic or commensurability claim.

## Frozen artifacts and replay

- [complete certificate](../data/enumerations/shimada_disc16_weighted_wall_symmetry_inventory.json)
- [certifier](../scripts/inventory_shimada_disc16_weighted_wall_symmetries.py)
- [focused tests](../tests/test_shimada_disc16_weighted_wall_symmetries.py)
- [archived integral-stabilizer denominator](shimada_disc16_atlas_fixed_section_inventory.md)

The complete reconstruction is deliberately expensive but deterministic:

```sh
PYTHONPATH=src:. python3 -u scripts/inventory_shimada_disc16_weighted_wall_symmetries.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_shimada_disc16_weighted_wall_symmetries.py
```

The primary `disc16.txt` member has SHA-256

```text
a19b797d6c314bc85e1ac4991fe3b190fe29a18a4746103bc6fa9b05e26a7db6
```

and the frozen weighted-wall certificate has SHA-256

```text
e2b1da1b0fbb8bac659dcb77bd59a640df17d95996e066c5065808ded9160fb8
```
