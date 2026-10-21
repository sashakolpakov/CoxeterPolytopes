# Silver-field H24 `E8^3` triage through H10

This is an exact low-cost screen of

\[
\langle-a\rangle\oplus E_8^3
\quad\text{over}\quad\mathbf Z[\sqrt2]
\]

for `a=sqrt(2)` and `a=2sqrt(2)`.  Both seeds have nonzero
positive-height Vinberg walls, but both exhibit the same decomposable
product-of-simplices first chamber.  This is a negative triage result, not a
new finite-volume polytope.

## Admissible coefficients and mirror norms

An arithmetic hyperbolic form must have `a>0` at the distinguished embedding
and `a'<0` at the conjugate embedding.  Therefore `sqrt(2)` and
`2sqrt(2)` are admissible.  The proposed `(1+sqrt(2))*sqrt(2)=2+sqrt(2)` is
positive at both embeddings, so its negative line stays negative after
conjugation; it is not admissible.

Write `pi=sqrt(2)`, so `(2)=(pi)^2`, and put
`I(r)=(a*x,G_E8^3*y)`.  Reflection integrality gives `(q(r))|2I(r)`.
Away from `(pi)`, unimodularity and primitivity make the pairing ideal a
unit.  At `(pi)`, evenness of `E8^3` supplies the lower valuation bound that
is missed by a bare ideal-divisibility argument:

- for `a=pi`, a primitive root has `v_pi(q)=1` when `x` is a unit, and
  `v_pi(q)=2` when a spatial coordinate is a unit;
- for `a=pi^3=2sqrt(2)`, it has `v_pi(q)=3` when `x` is a unit and the
  spatial vector is divisible by `pi`, and `v_pi(q)=2` otherwise.

The fundamental unit `epsilon=1+sqrt(2)` has signs `(+,-)` and
`O*=+/-epsilon^Z`; consequently every totally-positive unit is a square.
The complete totally-positive norm representatives, modulo root rescaling by
unit squares, are therefore:

| coefficient | complete canonical root norms |
|:--|:--|
| `sqrt(2)` | `2`, `2+sqrt(2)=epsilon*pi` |
| `2sqrt(2)` | `2`, `4+2sqrt(2)=epsilon*pi^3` |

## Exact H10 wall screen

The exact height-zero dominant census has three norm-2 orbits for each seed,
one per `E8`; these expand to exactly 720 roots.  Positive levels use complete
Minkowski rectangles, dominant E8 Dynkin weights, integral inverse-Cartan
coordinates, exact two-embedding Schur bounds, and simultaneous equal-height
packet selection.

| coefficient | controls | candidates | accepted positive walls | total walls |
|:--|--:|--:|--:|--:|
| `sqrt(2)` | 5 | 6 | 3 | 27 |
| `2sqrt(2)` | 3 | 6 | 3 | 27 |

For both seeds, the accepted height is `1+sqrt(2)`.  There is one accepted
wall supported on each E8 block.  Their norms are respectively `2+sqrt(2)`
and `4+2sqrt(2)`.  The other three H10 candidates are rejected by lower
walls.

## Exact first-packet geometry

In dominant Dynkin-weight coordinates, each accepted wall is a scalar
multiple of the E8 highest root and caps only its own block.  At the
`t=1` section, every block is an 8-simplex.  The complete cone section is
therefore their Cartesian product, with exactly `9^3=729` projective rays.
Exact norm classification gives the same result for both coefficients:

| finite | cusps | hyperideal | finite volume |
|--:|--:|--:|:--:|
| 7 | 0 | 722 | no |

This confirms the decomposable `E8^3` pathology and motivates the requested
irreducible `D24+` Niemeier lane.  A current hyperideal ray may still be cut
by a later wall, so this is not a nonreflectivity theorem.

The `D24+` lane cannot reuse three independent E8 dominant enumerators.  Its
height-zero root system is `D24`, while the ambient unimodular lattice is the
index-two spinor glue extension.  A correct implementation must enumerate
dominant D24 weights in **both** lattice cosets and derive the reflection
pairing ideal in the glued integral basis; treating D24's Cartan root lattice
as if it were already unimodular would give the wrong norm classes.  The
natural first seed is `a=2sqrt(2)`, using the same admissibility gate, but its
glue-aware backend remains the next bounded task rather than an inferred
result here.

## Reproduction

- Backend: `src/coxeter_polytopes/search/silver_e8_cubes.py`
- Generator: `scripts/search_silver_e8_cubes_h24.py`
- Tests: `tests/test_silver_e8_cubes.py` (4 passed)
- Archive: `data/enumerations/silver_e8_cubes_h24_h10.json`
- SHA-256:
  `998f7009f3133f255c5437439c07f30129c89e000f7dfa51b697a7a09976aca1`

```console
PYTHONPATH=src python3 scripts/search_silver_e8_cubes_h24.py --height 10
PYTHONPATH=src pytest -q tests/test_silver_e8_cubes.py
```
