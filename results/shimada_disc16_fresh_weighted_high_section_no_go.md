# Complete no-go census for the 575 fresh `disc16` weighted-wall sections

## The result

The full exact weighted-wall symmetry inventory of Shimada's 4,539
discriminant-16 chambers has 575 strict extensions beyond the archived
chamber stabilizers.  Their complete proper fixed-section stock above `H9`
consists of one involution section in `H13` and 574 in `H18`.

> **Every one of the 575 fresh sections is not a hyperbolic Coxeter
> polyhedron.**

For each section, at least one pair of walls is first proved to consist of
genuine hyperbolic facets by exact strict timelike witnesses, and then shown
to have a forbidden normalized product or the wrong acute orientation.
Thus this is an exact angle obstruction; no volume inference is involved.
There are no zero-bad candidates and hence no volume jobs in this branch.

This closes the *fresh strict-extension* denominator.  It does not erase the
positive chambers found in the separate archived-`Aut(D)` branch, including
the proved 27-wall and 30-wall examples in `H10`.

## Exact denominator

The parent inventory exhausts exact automorphisms of the complete signed
weighted wall Gram graph of every published `disc16` chamber.  All 575 strict
extensions occur over chambers with trivial archived stabilizer and full
weighted-wall group `C2`.  Consequently there is one proper subgroup class
per chamber:

| dimension | classes |
|---:|---:|
| `H13` | 1 |
| `H18` | 574 |
| **total** | **575** |

The original chamber wall counts, and the resulting involution-orbit counts,
are:

| published walls | cases | typical orbit walls |
|---:|---:|---:|
| 20 | 221 | 19; the exceptional `H13` case has 14 |
| 21 | 294 | 20 |
| 22 | 59 | 21 |
| 23 | 1 | 22 |

All orbit normals span the complete rational fixed space.

## Orientation correction

Shimada's `wallorbit` vectors record mirrors, but `_isom` rows are not
uniformly oriented as chamber half-spaces.  Using their raw signs can produce
a formally correct weighted graph and an incorrect chamber angle test.

The screen therefore orients every wall toward the projected published Weyl
point.  When that point lies on a recorded wall, the side is chosen opposite
the adjacent chamber's Weyl point; an exact dyadic perturbation then supplies
a strict interior point.  The weighted involution is required to preserve
these orientation signs before any folding is accepted.

This issue occurs already in the unique `H13` case: 17 raw rows have positive
Weyl pairing and three must be reversed.  The final obstruction is unchanged
in count but now belongs to the correctly oriented chamber.

## Divisor lattice versus mirror lattice

The new weighted-graph involutions lift to rational Lorentz isometries, not
generally integral Néron--Severi isometries.  The screen therefore keeps two
objects separate.

The algebraic object is the saturated divisor intersection

```text
NS(X) intersection Fix_Q(H),
```

computed as an integral Smith kernel after clearing the rational action's
denominators.  Hodge index gives its exact Lorentz signature.  Restricted
wall covectors are then represented in its rational Lorentz space.  For each
mirror, the corresponding reflection is separately tested for integrality on
the saturated divisor lattice.  No artifact calls a rational weighted-graph
lift a surface automorphism, and no reflection-stable mirror lattice is
silently identified with the divisor intersection.

In fact, none of the 575 saturated fixed divisor lattices is stable under all
retained mirror reflections.  This is not needed for the geometric angle
no-go, but it is important arithmetic lineage information.

## Facet certification

Every conic redundancy is an exact nonnegative rational combination of the
final conic facets.  A retained wall enters the angle gate only after a point
of positive Lorentz square is found which is zero on that wall and strictly
positive on every other folded inequality.

The ordinary orthogonal Weyl projection certifies most facets.  Two exact
fallbacks prevent numerical conditioning from creating false candidates:

- for a simplicial cone, exact slack coordinates perturb a timelike face ray
  into the relative interior;
- for a preliminary zero-bad nonsimplicial case, all face extreme rays are
  enumerated from the rank conditions, and an exact ray sum or dyadic
  perturbation gives the witness.

The latter step eliminated every temporary zero-bad case.  Some negative
cases retain unclassified additional conic walls, but their conclusion is
still complete: the archived forbidden pair lies among two independently
certified genuine facets, so any remaining wall status cannot repair the
Coxeter failure.  All facets are resolved in 225 cases; the other 350 already
have a decisive certified obstruction.

## Angle census

There are 5,524 forbidden certified facet pairs in total.  Every case has at
least one; the range is 1 to 23.  Twenty-nine cases attain the minimum of one
bad pair.  The complete case histogram is frozen in the manifest.

The unique `H13` case is `W0600_C001`.  Its 20 walls fold to 14 strict facets
in a saturated rank-14 divisor lattice of determinant `-7168`; clearing the
rational fixed-space equations requires denominator 7.  Its five forbidden
normalized squared products are

```text
2/5, 4/13, 9/65, 1/45, 1/117.
```

Thus the only fresh opportunity below `H18` is already killed at the exact
pair gate.

## Frozen artifacts and replay

- [full weighted-wall parent inventory](../data/enumerations/shimada_disc16_weighted_wall_symmetry_inventory.json)
- [575-case exact manifest](../data/enumerations/shimada_disc16_fresh_weighted_high_section_census.json)
- [per-case certificate directory](../data/enumerations/shimada_disc16_fresh_weighted_high_sections)
- [resumable screen and archive auditor](../scripts/screen_shimada_disc16_fresh_weighted_high_sections.py)
- [focused regression tests](../tests/test_shimada_disc16_fresh_weighted_high_sections.py)

The all-case audit recomputes every decisive exact Gram/pair gate and fully
replays the divisor, mirror, redundancy, and witness arithmetic on selected
orientation/fallback representatives:

```sh
PYTHONPATH=src:. python3 scripts/screen_shimada_disc16_fresh_weighted_high_sections.py --audit-archive
PYTHONPATH=src:. python3 -m pytest -q tests/test_shimada_disc16_fresh_weighted_high_sections.py
```
