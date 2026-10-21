# Five remote geodesic lines in the 27-wall chamber

The complete edge census finds **446 actual geodesic edges** of the frozen [27-wall H9 chamber](fresh_h9_makarov_joint_20261013.md). Besides the [441 edges incident to finite vertices](fresh_h9_joint_frontier_20261013.md), there are exactly **five full geodesic lines with no finite or ideal vertex at either end**. Their ten open ends provide new construction targets that the finite-corner profile bank did not cover.

The [exact certificate](fresh_h9_joint_frontier_20261013_all_edges.json) exhausts all 1,105 eight-cliques of the finite-pair graph. Exactly 446 have positive-definite Gram; the other 659 are rejected. Their existing endpoint classification is:

| Edge type | Count |
| --- | ---: |
| Two finite endpoints | 144 |
| One finite endpoint and one open end | 297 |
| Two open ends, no finite endpoint | 5 |
| Ideal vertex extensions of these supports | 0 |

## Why this census includes all edges

Every positive-definite eight-wall support is an actual codimension-eight face. Project the strict center perpendicular to its normals. The inverse of its acute positive-definite Gram is entrywise nonnegative, so the projected point remains strict on every other inequality. Its norm remains negative. Thus the support defines a genuine hyperbolic geodesic edge, not merely a formal positive-definite subdiagram.

Conversely, the incident normals at an interior point of an edge form the simple rank-eight spherical chamber normal to that edge. Hence every edge occurs in this eight-clique census. For every support, the checker examines all remaining existing walls: a finite endpoint requires a positive ninth Schur complement; an ideal endpoint requires zero. All positive cases agree with the previously complete 65-corner list. No zero case occurs.

For each of the five newly found lines, the archive constructs a strict relative-interior timelike point `p` and an orthogonal spacelike direction `d` in its Lorentzian two-plane. It supplies two rational parameters with feasible spacelike points `p+t*d` and `p-t*d`. Both points have exactly the eight support equalities and are strict on all other walls. Convexity certifies the whole intervening hyperbolic line and strictness at both ideal ends. These proofs resolve the relevant geometry without trusting the numerical scout's exterior incidence sets.

## New endpoint profiles

All IDs refer to the current 27-wall source. The new line IDs extend the earlier 0–440 edge list.

| Line | Eight incident walls | Finite spherical profiles | Isolated dihedral families |
| --- | --- | ---: | ---: |
| 441 | 0,2,3,4,6,7,9,25 | 26 | 2 |
| 442 | 0,2,3,4,6,8,9,25 | 26 | 2 |
| 443 | 0,2,3,4,6,9,25,26 | 21 | 2 |
| 444 | 0,2,4,6,7,9,25,26 | 35 | 4 |
| 445 | 0,2,4,6,8,9,25,26 | 35 | 4 |

There are **143 additional finite profile records and 14 symbolic `I2(m)` families**, with every integer `m>=3` retained. The forest, valence and exact Schur-complement argument is the same complete local spherical classification used for the finite-corner frontier. It applies independently to both ends of each line. All walls, including polar wall 6, are normalized to norm two.

Line **443** is the only one of these five whose finite profile list allows a nonright contact to a current cap. A new normal may have a single label-three contact to cap 26, all other seven support contacts right, with residual norm `1`. Alternatively it may contact cap 25 at label three, with residual norm `1-sqrt(5)/5`. On the other four lines every finite profile is right to its current caps. The complete data were sent directly to the joint cap search, which is now targeting these remote ends.

These are necessary local profiles, not certified completing walls. In particular, the data do not assert that a new normal can cut an open end while remaining compatible with all 27 existing walls. New affine endpoint profiles are outside this spherical list.

## Replay and independence

The [producer](../scripts/fresh_h9_joint_frontier_20261013_all_edges.py) reuses the frozen exact arithmetic but leaves the prior finite-frontier files unchanged. The [separate checker](../scripts/fresh_h9_joint_frontier_20261013_all_edges_independent.py) imports only the earlier independent `Fraction`/rational-interval implementation. It independently exhausts the eight-cliques, verifies every endpoint classification, replays all ten spacelike continuations, and re-enumerates the 143 profiles and 14 tails. Its [audit result](fresh_h9_joint_frontier_20261013_all_edges_independent.json) binds the source chain.

```sh
python3 scripts/fresh_h9_joint_frontier_20261013_all_edges.py --verify
python3 scripts/fresh_h9_joint_frontier_20261013_all_edges_independent.py --verify
```

Both implementations pass. The chamber remains infinite-volume; this is an exact new frontier for further construction, not a new compact or finite-volume example.
