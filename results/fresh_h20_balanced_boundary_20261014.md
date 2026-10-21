# Balanced H20: coherent contact chains and a whole-family escape

There is an exact **1,206-facet Coxeter chamber** in the unit-inradius `A6 × A6 × A7` model. Seven inner caps connect the central facet to two complete 392-cap outer layers. The second layer cuts all 392 previously obstructing two-fifths probes, while all 392 original finite endpoints remain actual. **The chamber still has infinite volume:** one exact strict spacelike point survives even the larger, not necessarily compatible, 1,176-cap three-level bank.

The construction retains the 415-wall baseline consisting of the 22 seed walls, the central radius-two cap, and all 392 radius-`sqrt(28)` vertex caps. It **omits the three-cap motif** in the earlier [418-wall report](fresh_h20_balanced_20261014_bridge.md). No compatibility with that motif or with the separate integer-gap construction is asserted.

The [exact checker](../scripts/fresh_h20_balanced_boundary_20261014_chain.py) and [certificate](fresh_h20_balanced_boundary_20261014_chain.json) use `Fraction` arithmetic in `Q(sqrt(65))`. Generation and saved replay passed:

```sh
python3 scripts/fresh_h20_balanced_boundary_20261014_chain.py --verify
```

The baseline [exact archive](fresh_h20_balanced_20261014_bridge_exact.json), all of its declared source hashes, the shared [coordinate module](../scripts/fresh_h20_balanced_20261014.py), and the new producer are hash-bound. Baseline geometry is a dependency; the new checker is not a second implementation of the old pair census. A separate [rational-matrix audit](fresh_h20_balanced_bridge_review_20261014.json) already checked that predecessor's geometry.

## Exact construction

Use spatial Gram `B = diag(A6/2, A6/2, A7/2)` and Lorentz coordinates `(x,D,U)` with square `x²−2DU`. Write `z_abc` for the 392 product vertices, indexed by `0≤a,b≤6`, `0≤c≤7`. Every `z` has square `196`. The least nonzero vertex distance squared is `84`.

A sphere with center `c` and radius `R` has outward normal

```text
r(c,R) = (c, (c²−R²)/2, 1),       r² = R².
```

Set

```text
k = (33−sqrt(65))/32,
k²−(33/16)k+1 = 0,       0<k<1,
C_j(z): center (1−k^j)z, radius 2k^j.
```

Select `C_1(z_aaa)` for `a=0,...,6`, and select **every** `C_2(z)` and `C_3(z)`. These are 791 new caps. The seven first-level directions have minimum squared separation `280`.

Every selected cap has 19 seed contacts at angle `pi/3`; its other three seed pairs are strictly dotted. The first-level caps meet the central sphere at `pi/3`; levels two and three are strictly dotted to it. All selected caps are strictly dotted to all 392 old vertex caps.

For two chain caps, put `s_i=k^i`, `g_i=1−s_i`, and `d=|z−w|²`. Their center distance squared is exactly

```text
196(s_i−s_j)² + g_i g_j d.
```

The coefficient of `d` is positive. Consequently the following eight exact gates cover all new mutual pairs, with no numerical distance tolerance:

| levels | vertex-distance condition | relation |
|---|---:|---|
| 1,1 | distinct selected directions, `d≥280` | dotted |
| 1,2 | same direction | `pi/3` |
| 1,2 | different directions, `d≥84` | dotted |
| 1,3 | `d≥0` | dotted |
| 2,2 | different directions, `d≥84` | dotted |
| 2,3 | same direction | `pi/3` |
| 2,3 | different directions, `d≥84` | dotted |
| 3,3 | different directions, `d≥84` | dotted |

There are exactly 640,710 new pairs: 15,435 at `pi/3` and 625,275 dotted. Including the baseline gives 726,615 pairs: 16,263 at `pi/3`, 7,657 right, 392 parallel, and 702,303 dotted.

The point `(0,200,1)` has square `−400` and is strict on every wall. Projecting it orthogonally to any one positive-norm wall remains timelike and stays strict on all other inequalities, because all distinct wall pairings are nonpositive. This proves that all 1,206 inequalities define actual facets.

Every old endpoint `(z,112,1)` has square `−28`. Its old 19 seed contacts and its own vertex cap remain a positive-definite rank-20 incidence set, and every added cap is strictly negative there. Thus all 392 old finite endpoints survive. This is not a complete census of new vertices or cusps.

For each old vertex, `C_2(z)` has negative ball power at `(2/5)z`. It therefore cuts the previous strict spacelike lift of square `1/4`, and also removes that ideal-boundary probe. This checks a full orbit of old gaps; it does not prove a global covering.

## Exact escape for all three complete layers

Allow all 392 first-level caps as well, without imposing their mutual compatibility. This gives a relaxed bank of 1,176 caps. Take

```text
x = (3/5)z_000,
E = (x, 14087/400, 1),       E² = 1/8.
```

For level `j`, the least ball power over all product vertices occurs at `z_000` and equals

```text
P_j = 196(k^j−2/5)² − 4k^(2j).
```

Indeed the power at another vertex adds `(3/5)(1−k^j)|z_000−z|²≥0`. Exact field signs give `P_j>1/8` for `j=1,2,3`; their approximate values are `25.7696`, `6.94869`, and `0.156626`. The central power is `1664/25`, the minimum old vertex-cap power is `84/25`, and the minimum seed gap is `2/5`. Thus `E` is strict on the baseline and the entire relaxed bank. Its pairing with the strict center is `−94087/400<0`.

Together with a strict timelike point, this future spacelike point gives an open ideal-boundary patch and hence infinite volume. In particular, **no finite selection from these three fixed chain levels can complete the retained baseline to finite volume**. This is not an obstruction to different endpoint radii, changed chain profiles, or additional cap families.

The next same-ray level is incompatible with the retained endpoint layer: `C_4(z)` has positive normal pairing with the old cap at `z`. Continuing the same chain therefore requires an explicit replacement of that old facet.

The chain branch also lies outside the quasi-arithmetic route, conditionally on a finite-volume completion retaining its facets. A first-level cap and an omitted A6 seed wall have Cartan two-cycle

```text
(2353+273sqrt(65))/512.
```

Hence the cyclic field contains `Q(sqrt(65))`. The retained rational seed/central Gram contains an indefinite restriction `[[0,−7], [−7,4]]`; it remains indefinite under every field embedding. An embedding moving `sqrt(65)` therefore contradicts Vinberg's semidefiniteness condition for quasi-arithmeticity. This uses the finite-volume criterion stated in [Bogachev–Yorov, Theorem 2.2](https://nvbogachev.github.io/Quasi_arithmetic_Coxeter_prisms-2024.pdf); it is not a claim that the present infinite-volume chamber is a lattice.

## Separate bounded boundary-center scout

Before constructing the chains, the [numerical scout](../scripts/fresh_h20_balanced_boundary_20261014.py) tested 64,829 profiles on one old 19-wall vertex support. At most three prescribed nonright seed contacts were allowed, with labels `3,4,5,6`; the other retained contacts were right. The three omitted seed contacts could be finite at labels `2,...,60`, parallel, or dotted. Central labels were `2,...,16` or parallel. Both quadratic radius branches were considered, restricting to `R>1.8`, which is necessary for a central-contact ball to reach a point of norm `5.6`.

The [saved numerical result](fresh_h20_balanced_boundary_20261014.json) has 1,037,264 profile/central-label trials and zero surviving proposals. This is a bounded floating-point result, not an exact all-label or all-support exclusion. It motivated the change to chains whose outer caps are dotted to the central facet. The earlier [all-radius endpoint review](fresh_h20_balanced_half_20261014_review.md) and radius-five central-contact bound do not exclude those outer caps.
