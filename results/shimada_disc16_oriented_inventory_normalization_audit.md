# Discriminant-16 orientation-first inventory audit

## Result

The complete discriminant-16 weighted-wall fishing denominator is unchanged
when every wall is oriented toward Shimada's published chamber **before** the
exact Gram-automorphism group is computed.

The orientation-first reconstruction covers all 4,539 pinned chambers.  It
resolves 105 wall orbits whose published Weyl pairing is zero by using the
archived opposite adjacent chamber.  The resulting exact census is

| item | exact count |
|---|---:|
| chambers | 4,539 |
| chambers with nontrivial oriented weighted wall group | 673 |
| strict extensions beyond archived `Aut(D)` | 575 |
| subgroups in nontrivial groups | 1,496 |
| subgroup conjugacy classes | 1,407 |
| proper classes with fixed rank at least 11 | 682 |

The high-dimensional split is

```text
H18:574, H14:2, H13:1, H11:101, H10:4.
```

This is exactly the earlier complete raw-mirror weighted inventory.  More
strongly:

- all 4,539 chamber wall counts, archived group orders, full group orders,
  and strict-extension flags agree;
- for every one of the 673 nontrivial groups, the complete multiset of
  subgroup-class order, fixed rank, conjugacy-class size, and element-order
  data agrees; and
- for all 575 strict extensions, the full wall permutation and exact rational
  dual action of every group element agrees entrywise.

Thus no new oriented-wall action or high-rank subgroup class was omitted by
the earlier split screens.  The 107 archived proper classes are precisely
the already closed archived-`Aut(D)` lane, and the 575 strict-extension
classes are precisely the already closed fresh lane.  A second facet/angle
screen would duplicate identical actions and inequalities.

The complete primal-action denominator histogram is

```text
integral:4735, denominator 2:574, denominator 7:1.
```

The 574 denominator-two actions are the primitive `(-4)`, divisibility-one
`H18` family.  The unique denominator-seven action is the previously screened
`H13` involution at chamber 600.  No additional denominator or dimension
appears.

This audit is a denominator-normalization theorem only.  The no-go decisions
come from the already frozen exact strict-facet Coxeter screens; it does not
reinterpret a non-Coxeter cast as an infinite-volume polytope.

## Reproduction

```bash
PYTHONPATH=.:src python3 scripts/inventory_shimada_oriented_weighted_wall_symmetries.py --discriminant 16
PYTHONPATH=.:src python3 scripts/certify_shimada_disc16_oriented_inventory_normalization.py
PYTHONPATH=.:src pytest -q tests/test_shimada_disc16_oriented_inventory_normalization.py
```

Frozen artifacts:

- `data/enumerations/shimada_disc16_oriented_weighted_wall_symmetry_inventory.json`
- `data/enumerations/shimada_disc16_oriented_weighted_wall_chambers/`
- `data/enumerations/shimada_disc16_oriented_inventory_normalization_audit.json`

