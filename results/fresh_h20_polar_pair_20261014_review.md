# Independent review of the rational polar pair and its new vertices

**PASS.** The [independent checker](../scripts/fresh_h20_polar_pair_20261014_review.py) and [certificate](fresh_h20_polar_pair_20261014_review.json) validate the two polar walls, the complete finite-vertex census on them, and the surviving strict escape. Generation and saved replay pass. No mathematical issue was found in the [primary construction](fresh_h20_polar_pair_20261014.md) or the [separate census](fresh_h20_polar_pair_20261014_census.json).

The review imports the earlier independent sparse-radical implementation, not the polar producer's arithmetic. It reconstructs all inherited normals needed for the new calculations. The frozen 98-wall all-pair and old-vertex certificates remain explicit dependencies; their full censuses are not repeated.

## New walls and geometry

The literal rational normals have squares 2 and 4 and are mutually right. All **197** new pairs are checked independently:

| Relation | Count |
|---|---:|
| Right | 53 |
| `pi/3` | 4 |
| `pi/4` | 3 |
| Parallel | 12 |
| Strictly dotted | 125 |

Strict dotted inequalities are certified by outward rational intervals; finite and parallel equalities are checked in the exact field. The complete finite-neighbor lists contain 36 and 25 roots respectively, including the other new wall. Every old finite neighbor lies in the rational pool used by the census, so that pool is exhaustive for vertices incident to either new wall.

The common center is reconstructed as the sum of the first 95 physical normals. Its old norm and signs come from the bound independent whole-face certificate. The three later old-wall pairings follow by exact summation of their already certified pairings; the two new pairings are independently evaluated. All 100 signs are strict.

For each new wall, the checker reconstructs its facet point `T-(T·r/r²)r`, its negative square, and all 100 pairings from the projection identity. Thus both new facets are actual; the common-center acute-extension argument also preserves all old facets. The two normals themselves are feasible spacelike old rays. Their old zero sets have independently checked rational rank 20, confirming the polar extreme-ray interpretation.

The three inherited bounded cusps pair strictly negatively with both new normals. Their existing local links therefore remain unchanged. Every old finite vertex is also strictly retained by the acute-extension lemma.

## Completeness of the new-vertex census

The independent enumeration uses a different route from the census producer. Starting from the verified full finite-neighbor graph, a direct fixed-size clique recursion enumerates all **20,810** distinct 20-root candidates containing a new wall. It does not use the producer's maximal-clique list.

Positive definiteness is then tested by graph structure and exact leaf elimination:

- A connected support graph containing a cycle cannot be positive definite. Normalize its diagonal entries to 1. Every nonzero edge has cosine at least `1/2`, so the vector equal to 1 on a cycle and 0 elsewhere has nonpositive square. Extra chords only decrease it.
- On a tree, eliminate a leaf with diagonal `d` and edge pairing `g`, replacing its neighbor's diagonal by `d_neighbor-g²/d`. Every eliminated pivot and the final remaining diagonal must be positive. This is exact rational Schur elimination without fill-in.

These tests independently agree with every one of the producer's **17,171** cached component decisions: 14,716 cycle obstructions, 1,977 nonpositive leaf pivots, 299 nonpositive final pivots, and 179 positive tree components. No finite Coxeter classification catalogue or numerical eigenvalue threshold is needed.

Exactly **16** positive 20-root supports remain. The saved primitive rational point for each is independently evaluated against all 100 physical normals, giving **1,600 exact inequalities**. Every point is future timelike, has precisely its 20 stated zero pairings, and is strict at the other 80 walls. The points are distinct.

There are 14 vertices on wall 98 and 14 on wall 99; 12 lie on both. This is a complete census on the union of the two new facets. Any finite vertex there contains a positive-definite 20-root support in the enumerated finite graph; the strict signs exclude additional incident facets. It is not a complete census of the old part of the chamber.

The inherited lower bound of 118 finite vertices is strictly preserved. The 16 vertices on the new facets are new, so the final chamber has **at least 134 finite vertices**.

## Volume and the surviving H2 face

The old rational escape still has square

`20181174579467302741 / 10^18 > 0`,

and pairs strictly negatively with both new walls. Its inherited 98 strict inequalities are unchanged. The resulting 100-wall chamber therefore remains infinite-volume.

Both new walls are strictly dotted to wall 75. Since that wall has norm 2 and inverse support diagonal 5 in S75, the nonnegative support inverse gives projected square strictly less than `-9q` for a new normal of square `q`. Projection of the common strict center fixes the future orientation. Thus both polar additions leave the entire S75 hyperbolic plane strict, consistent with the [complete S75 projection gate](fresh_h20_s75_projection_gate_20261014.md). This last statement is a Schur-bound consequence of the checked pair data, not a new full projection census.

```sh
python3 scripts/fresh_h20_polar_pair_20261014_review.py --verify
```

The review binds the final primary and census source/data plus its arithmetic dependencies. It does not claim finite volume, a complete cusp census, or a full census away from the new facets.
