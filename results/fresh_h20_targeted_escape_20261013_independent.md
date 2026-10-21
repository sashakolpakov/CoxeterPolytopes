# Standalone audits of the 86th, 87th and 88th H20 walls

**PASS.** The new wall cuts the former 85-wall chamber's strict escape and
creates a genuine finite vertex. A different exact strict escape proves
that the resulting 86-facet chamber still has infinite volume.

The [standalone checker](../scripts/fresh_h20_targeted_escape_20261013_independent.py)
imports no construction or project arithmetic module. It implements Fraction
arithmetic in the basis 1, √6, √10, √60 and uses exact comparisons of squares
for signs. Relations with the three earlier radical caps require one further
quadratic extension, handled separately. Its [saved result](fresh_h20_targeted_escape_20261013_independent.json)
checks all 85 new pairs: nineteen right, three π/4 and 63 dotted.

Let v be the saved rejected polar 27, whose last coordinate is −1 and whose
square is 1. The new normal retains its first twenty coordinates and changes
the last to −4+√10. Its square is 76−24√10>0; its pairing with the strict
center is −110. Its pairing with the former strict escape is
1225/128−3√10>0, so that escape is removed. All 107 inherited cusp rays remain
strictly inside the new halfspace.

Writing t=4−√10 and v₀ for the first twenty coordinates, the point
x=(3t v₀,−2) has square −6(76−24√10)<0. It is incident to the new wall and
nineteen independent old walls, and strict on every other wall. The checker
verifies the incident rank, orientation, exact equalities, and strict signs.
The existing 85-wall certificate remains a hash-anchored dependency.

The new [strict escape](../data/enumerations/fresh_h20_targeted_escape_20261013_volume.json)
has square 11485671335440549/6400000000000000000>0 and center pairing
−1829587279/160000000<0. All 86 halfspace inequalities are independently
checked strictly. Together with the strict timelike center, this gives an
open set at infinity and proves infinite volume.

Forged normal, norm, and old-escape-cut values are rejected. Reproduce with:

```sh
python scripts/fresh_h20_targeted_escape_20261013_independent.py --verify
```

The audited [new-wall archive](../data/enumerations/fresh_h20_targeted_escape_20261013_exact.json)
is frozen. Later extensions have separate artifacts. No finite-volume or
global novelty claim follows from this checkpoint.

## Independent replay of the 87th wall

The [next checker](../scripts/fresh_h20_targeted_escape_20261013_second_independent.py)
also passes. It uses a separate parameterized Fraction implementation of
Q(√6, √(268/3+36√6)), followed by one extra radical for each cross pair.
It imports no producer or project arithmetic implementation.

The [87th normal](../data/enumerations/fresh_h20_targeted_escape_20261013_second.json)
has square 2 and passes all 86 previous-wall relations: eighteen right,
one π/3, and 67 dotted. It cuts the 86-wall escape, preserves the strict
center and all 107 inherited cusp rays, and retains every old facet.
The [independent result](fresh_h20_targeted_escape_20261013_second_independent.json)
also rejects forged normal, norm and cut values.

A different strict escape has square 73993/8388608>0 and center pairing
−28919/2048. All 87 strict inequalities pass the separate arithmetic,
so this stage also has infinite volume. Reproduce with:

```sh
python scripts/fresh_h20_targeted_escape_20261013_second_independent.py --verify
```

## Independent replay of the 88th wall

The [third checker](../scripts/fresh_h20_targeted_escape_20261013_third_independent.py)
passes both its initial run and deterministic replay. It reuses the preceding
independent Fraction arithmetic, imports no producer, and checks all 87 new
pairs: eighteen right angles, two of angle π/3, and 67 dotted relations.
The cap has norm two, preserves the strict center and all 107 inherited
cusp rays, and cuts an exact strict point of the preceding chamber.

The actual finite vertex lies on the eighteen saved short walls, old wall
42, and new wall 87. Its twenty incident normals have a rational
positive-definite Gram matrix, checked by exact congruence pivots. The
saved vertex is future timelike and strict on all other walls. The
[audit result](fresh_h20_targeted_escape_20261013_third_independent.json)
also rejects forged normal, norm and target-cut values.

The final strict escape has square
478875483478515721/32000000000000000000>0 and center pairing
−48815755199/4000000000. Every one of its 88 wall inequalities is strictly
negative in the independent arithmetic. Thus the chamber still has
infinite volume. This audit does not certify the numerical 1,593-cap bank.

```sh
python scripts/fresh_h20_targeted_escape_20261013_third_independent.py --verify
```
