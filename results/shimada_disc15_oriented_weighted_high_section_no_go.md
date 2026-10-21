# Determinant-15 oriented-wall fixed-section no-go

## Result

The complete pinned Shimada determinant-15 atlas contains no Coxeter
polytope in the following high-dimensional fixed-section denominator.

For every one of the 2,051 published induced chambers (D), orient every
published wall toward the chamber side, let (W_{\mathrm{wall}}(D)) be the
full exact automorphism group of the resulting oriented weighted wall Gram
matrix, and take every inner-conjugacy class of proper subgroups
(H\leq W_{\mathrm{wall}}(D)) for which

\[
 \operatorname{rank}\bigl(NS\cap \operatorname{Fix}_{\mathbb Q}(H)\bigr)
 \geq 11.
\]

There are exactly 236 such classes.  In every class, two independently
certified hyperbolic facets of the folded cone have a forbidden normalized
Gram product.  Hence none of the 236 folded cones is a hyperbolic Coxeter
polytope.  In particular there is no finite-volume or compact example in
this denominator; the volume gate is never reached.

This is a denominator theorem, not a universal non-existence theorem for
determinant 15 and not a classification of all divisor constructions on all
K3 surfaces.

## Exact denominator

The source is `/private/tmp/threediscs/disc15.txt`, pinned by SHA-256

```text
f2fdcf80bd980c363fbdb68957a85878284bdb5e7d75c4dac92604aefefc4081
```

Its rank-20 Gram determinant is (-15).  The complete chamber census is:

| item | exact count |
|---|---:|
| published chambers | 2,051 |
| chambers with nontrivial archived (Aut(D)) | 53 |
| chambers with nontrivial full oriented wall group | 211 |
| strict extensions beyond archived (Aut(D)) | 171 |
| weighted-group subgroup classes, all ranks | 564 |
| proper classes with fixed rank at least 11 | 236 |

The full oriented group order histogram is

```text
1:1840, 2:196, 3:2, 4:6, 8:2, 12:1, 16:1, 24:2, 48:1.
```

The high-dimensional queue separates into three exact arithmetic families:

| hyperbolic dimension | classes | subgroup | relation to archived (Aut(D)) | primal action | determinant of saturated fixed divisor lattice |
|---:|---:|---:|---|---|---:|
| 18 | 144 | (C_2) | new wall symmetry | denominator 2 | (+60) |
| 13 | 33 | (C_2) | new wall symmetry | integral | (-320) |
| 11 | 59 | (C_2) | contained in archived (Aut(D)) | integral | (-3840) |

The archived stabilizers alone have 60 proper H11 classes.  Passing to inner
conjugacy in the larger oriented wall groups fuses one archived distinction,
leaving the 59 H11 classes in the full denominator.

## Orientation is part of the theorem

Shimada's wall rows specify mirrors, not a globally consistent choice of
half-space.  The inventory orients each entire wall orbit toward the
published projected Weyl vector.  Across all 2,051 chambers every orbit has
a uniform sign; no orbit has mixed Weyl-side signs.  There are 114 orbits on
which the published Weyl pairing vanishes.  For these, the explicitly
archived adjacent chamber fixes the opposite side (94 `_backto`, 18 `_isom`,
2 `_newcham`).

Only after this orientation is fixed is the complete exact weighted Gram
automorphism group enumerated.  Every archived (Aut(D)) action is checked
to preserve the oriented wall rows exactly and to lie in the resulting full
wall group.

## Divisor lattice, rational lift, and mirror lattice

For a subgroup (H), the lifted wall action may be rational on (NS\otimes
\mathbb Q).  The divisor lattice used by the screen is not a guessed
eigenlattice and not the lattice spanned by mirrors.  It is the saturated
intersection

\[
 L_H = NS\cap \operatorname{Fix}_{\mathbb Q}(H),
\]

computed as an integral Smith kernel after clearing the exact rational
denominators.  The complete rational constraint row space is reduced to an
independent exact basis before Smith decomposition and the resulting lattice
is verified against every subgroup action.

Restricted mirror directions are then computed in (L_H^\vee\).  Their
reflections are tested individually on (L_H); reflection integrality is
never inferred from a wall permutation.  None of the 236 saturated divisor
lattices is stable under every retained mirror reflection.  Thus the
algebro-geometric divisor lattice and a possible reflection-stable mirror
lattice remain distinct in every case.

## Exact Coxeter obstruction

Every published wall is folded by the subgroup.  Conic redundancy is solved
over (\mathbb Q), with each discarded wall sum expressed using the final
retained conic facets.  A retained spacelike mirror is called an actual
hyperbolic facet only after an exact timelike point is exhibited on that
mirror and strictly inside every other half-space.

The final exact pair gate found 5,260 forbidden pairs.  Every one of the 236
cases has at least one; 12 cases have exactly one.  There is no zero bin.
All facets were resolved in 149 cases.  In the remaining 87 cases some
additional facet-causality questions remain unresolved, but this does not
weaken the no-go: the forbidden pair is between two facets already proved
individually by strict exact witnesses.

The most common forbidden rational normalized products are (1/16) (1,469
occurrences), (5/8) (534), (1/6) (507), (3/128) (506), and (5/16)
(452).  These counts are diagnostics; one forbidden pair is already
decisive.

## Primal-dual status

This cast moved the theorem boundary: it exhausts a substantially larger
algebro-geometric fishing denominator than archived surface automorphisms,
including 177 genuinely new high-dimensional wall-symmetry classes (144 H18
and 33 H13).  It did not move the example boundary: no Coxeter candidate
survived, so no Vinberg finite-volume test was appropriate.

The next independent cast is the full pinned determinant-11 atlas.  Existing
directed RM/eigen-nef constructions do not replace that exhaustive
oriented-wall subgroup census.

## Reproduction

```bash
PYTHONPATH=.:src python3 scripts/inventory_shimada_oriented_weighted_wall_symmetries.py --discriminant 15
PYTHONPATH=.:src python3 scripts/screen_shimada_oriented_weighted_high_sections.py --discriminant 15
PYTHONPATH=.:src python3 scripts/screen_shimada_oriented_weighted_high_sections.py --discriminant 15 --audit-archive
PYTHONPATH=.:src pytest -q tests/test_shimada_disc15_oriented_weighted_high_sections.py
```

Frozen primary artifacts:

- `data/enumerations/shimada_disc15_oriented_weighted_wall_symmetry_inventory.json`
- `data/enumerations/shimada_disc15_oriented_weighted_wall_chambers/`
- `data/enumerations/shimada_disc15_oriented_weighted_high_section_census.json`
- `data/enumerations/shimada_disc15_oriented_weighted_high_sections/`

