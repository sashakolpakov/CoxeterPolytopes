# Determinant-11 oriented-wall fixed-section no-go

## Result

The complete pinned Shimada determinant-11 atlas contains no Coxeter
polytope in the full high-dimensional oriented-wall fixed-section
denominator.

For each of the 1,098 published induced chambers (D), let
(W_{\mathrm{wall}}(D)) be the full exact automorphism group of its inward
oriented weighted wall Gram matrix.  Enumerate every inner-conjugacy class of
proper subgroups (H\leq W_{\mathrm{wall}}(D)) satisfying

\[
 \operatorname{rank}\bigl(NS\cap \operatorname{Fix}_{\mathbb Q}(H)\bigr)
 \geq 11.
\]

There are exactly 294 classes.  Every class has two independently certified
hyperbolic facets with a forbidden normalized Gram product.  Hence none of
the folded cones is a hyperbolic Coxeter polytope, and none advances to the
compactness/finite-volume gate.

This is a complete theorem for the stated published-chamber and
oriented-wall-symmetry denominator.  It is not a universal non-existence
theorem for determinant 11 or for all divisor constructions on K3 surfaces.

## Pinned atlas and complete group census

The source `/private/tmp/threediscs/disc11.txt` has SHA-256

```text
fa39e5b59f4a4b2b1b7fe7632e5737d1734ee3eb7e7023ae26bf2b83387afb68
```

and a rank-20 Gram matrix of determinant (-11).

| item | exact count |
|---|---:|
| published chambers | 1,098 |
| chambers with nontrivial archived (Aut(D)) | 71 |
| chambers with nontrivial full oriented wall group | 284 |
| strict extensions beyond archived (Aut(D)) | 216 |
| weighted-group subgroup classes, all ranks | 603 |
| proper classes with fixed rank at least 11 | 294 |

The full oriented wall-group order histogram is

```text
1:814, 2:270, 3:2, 4:5, 6:5, 8:1, 20:1.
```

The archived stabilizers alone contribute 69 H11 classes.  The full exact
wall denominator contributes 225 further high-dimensional classes.

## Arithmetic families

All 294 classes fall into the following exact families:

| dimension | exact subgroup classes | relation to archived (Aut(D)) | primal lift | determinant of saturated fixed divisor lattice |
|---:|---|---|---|---:|
| H18 | (216\times C_2) | new | denominator 2 | (+44) |
| H17 | (C_2,C_4) | new | denominator 2 | (-176) |
| H12 | (2\times C_2) | new | denominator 2 | (+2816) |
| H11 | (69\times C_2) | archived | integral | (-2816) |
| H11 | (2\times C_4) | new | denominator 2 | (-2816) |
| H10 | (2\times C_4,C_8) | new | denominator 2 | (+5632) |

The H17 and H10 families arise in chamber 654, where the archived order-2
stabilizer extends to a full oriented wall group of order 8.  All seven
high-rank subgroup classes of that chamber have every facet resolved; their
forbidden-pair counts are 6, 6, 3, 6, 3, 3, and 3.

## Exact orientation and action separation

Unlike determinant 15, no determinant-11 wall orbit lies on the published
projected Weyl vector.  Every wall orbit has one strict uniform Weyl-side
sign.  The half-spaces are oriented before the weighted Gram automorphism
group is computed, and every archived (Aut(D)) action is checked to
preserve those oriented rows exactly.

The wall-group lift is allowed to be rational.  The divisor lattice is then
computed independently as the saturated Smith kernel

\[
 L_H=NS\cap\operatorname{Fix}_{\mathbb Q}(H).
\]

Denominators are cleared only for the fixed-space equations; rational action
matrices are never silently treated as integral surface automorphisms.
Likewise, the restricted mirror lattice is not identified with (L_H).
Every active mirror reflection is tested on (L_H), and none of the 294
fixed divisor lattices is stable under all active reflections.

## Exact Coxeter obstruction

Every published wall is folded.  Conic redundancy is replayed over
(\mathbb Q), and an inequality counts as a hyperbolic facet only after an
exact timelike point is exhibited on it and strictly inside every other
half-space.  The pair gate then uses exact normalized Gram products.

Across the 294 cases it certifies 3,557 forbidden facet pairs.  Three cases
have exactly one forbidden pair and no case has zero.  All facets are
resolved in 136 cases.  The other 158 retain some unresolved *additional*
facet-causality questions, but the no-go conclusion is rigorous: in each
case at least one forbidden pair joins two facets already established by
independent strict exact witnesses.

The two most frequent forbidden products are (4/81) (807 occurrences) and
(11/36) (669).  These frequencies are diagnostic only; a single forbidden
pair is decisive.

## Primal-dual status and lineage

This closes an exhaustive chamber-fixed subgroup denominator that is
different from the earlier directed determinant-11 RM/eigen-nef searches.
The larger wall groups reveal H18, H17, H12, and H10 arithmetic sections that
are invisible in the archived surface-automorphism denominator.  None is a
Coxeter polytope.

Together with the determinant-15 census, this moves the theorem boundary by
530 exact high-dimensional no-go classes.  It does not move the example
boundary and does not lower a universal Vinberg non-existence dimension
bound, because the quantifier is over these pinned algebraic wall stocks and
their complete oriented symmetry groups.

## Reproduction

```bash
PYTHONPATH=.:src python3 scripts/inventory_shimada_oriented_weighted_wall_symmetries.py --discriminant 11
PYTHONPATH=.:src python3 scripts/screen_shimada_oriented_weighted_high_sections.py --discriminant 11
PYTHONPATH=.:src python3 scripts/screen_shimada_oriented_weighted_high_sections.py --discriminant 11 --audit-archive
PYTHONPATH=.:src pytest -q tests/test_shimada_disc11_oriented_weighted_high_sections.py
```

Frozen primary artifacts:

- `data/enumerations/shimada_disc11_oriented_weighted_wall_symmetry_inventory.json`
- `data/enumerations/shimada_disc11_oriented_weighted_wall_chambers/`
- `data/enumerations/shimada_disc11_oriented_weighted_high_section_census.json`
- `data/enumerations/shimada_disc11_oriented_weighted_high_sections/`

