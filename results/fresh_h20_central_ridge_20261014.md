# Finite vertices on every central golden ridge, and one finite edge

Starting from the weighted 39-wall core, an explicit 288-cap construction creates a finite vertex on each of its sixteen central π/5 ridges. One further cap closes an actual edge on a selected ridge while retaining all preceding walls. The final chamber has **328 actual facets, exactly 17 finite vertices, and one finite edge**. A saved strict spacelike point proves that its volume remains infinite.

This construction uses only the 39-wall core. The preceding binary and outer-code banks are not constraints here. It establishes the requested finite-corner construction; it does not claim a global completion.

The [independent review](fresh_h20_central_ridge_20261014_review.md) passed with a different field basis and no producer imports, including all 4,608 vertex/cap checks and the complete graph counts. Both primary saved replays passed.

## Explicit orthogonal clusters

Use physical Euclidean spatial coordinates, with ten unit inradii, one inradius √2 and eight inradii φ. The Lorentz form is `x²−2DU`. A vertical side has unit normal `(±e_j,b_j,0)`, and the central normal is `a=(0,−2,1)`, of square four.

Fix a signed golden side `b=(σe_g,φ,0)`. Put

```text
L = 3−φ,       R = 1/3,
s = sqrt(R²+18L),
α = (s−R)/18,       0<α<1/3.
```

For each of the eighteen coordinates `i≠g`, define

```text
c_i = σφ e_g + α 1_rest + R e_i,
r_i = (c_i,2,1).
```

The identity `18α²+2Rα=L` gives

```text
r_i² = R²,       r_i·r_j = 0 for i≠j,
r_i·a = r_i·b = 0.
```

All other core walls are strictly dotted to these caps. Indeed, every nontarget center coordinate lies between zero and `α+R<2/3`, while every inradius is at least one. Its smaller side gap is therefore greater than `R`. The opposite target gap is `2φ>R`.

The finite corner is explicit:

```text
P = (σφ e_g + (L/s) 1_rest, 2, 1),
P² = −LR²/(R²+18L) < 0.
```

Its physical incident Gram is `[[4,−φ],[−φ,1]] ⊕ R²I₁₈`, positive definite because the first block has determinant L. Every other core inequality is strict.

Repeat for all eight golden coordinates and both signs. Different targets have center differences exceeding `φ−2/3>5/6` in two coordinates, so their squared distances exceed `25/18>4R²`; their caps are strictly dotted. Opposite signs at the same target are also dotted. The [327-wall checker](../scripts/fresh_h20_central_ridge_20261014.py) and [certificate](fresh_h20_central_ridge_20261014.json) verify the exact formulas, all 41,328 cap pairs, and all sixteen points against every cap. There are 2,448 right cap pairs and 38,880 dotted cap pairs; the 11,232 new/core pairs split into 576 right and 10,656 dotted.

The arithmetic field is `Q(φ,α)`. The rational norm of `R²+18L` is `132031/81`, whose numerator is nonsquare; the quadratic representation and sign decisions have an explicit nonsquare premise. The original heavy side length √2 is needed only in the core and its uniform lower-bound comparisons.

The timelike point `(0,100,1)` is strict for all walls. Acute pairings and projection to any individual wall prove that every wall is an actual facet.

## Complete vertex and edge counts

Every cap has exactly nineteen other finite neighbors: its seventeen cluster companions and a,b. Thus any finite vertex incident to it must be the explicit twenty-wall corner above. A finite vertex using only old core walls would require the central wall and one side of every coordinate; its Gram Schur complement is `4−(20+8φ)<0`. There are consequently exactly sixteen finite vertices in the 327-wall chamber.

The same graph argument completely counts the positive-definite rank-nineteen subsets, and hence its hyperbolic edges. Cap-containing subsets give twenty edges at each of the sixteen corners, all open in the other direction. Any possible further old wall has a dotted pair with a retained root. No further wall can supply a second finite vertex or a bounded-link ideal vertex. These rays do reach the ideal boundary; their limiting links are unbounded. The certificate's “without second finite or ideal endpoint” field refers to this absence of a terminating vertex, not absence of an ideal-boundary limit.

Without a cap, a subset containing the central wall and eighteen vertical walls has negative Schur complement, since each side weight is at least one. The remaining subsets choose one of each nineteen opposite vertical side pairs: exactly `2¹⁹=524,288` vertical edges. They join the common infinity cusp, with its bounded box link, to distinct Euclidean box vertices with unbounded links. The latter points lie strictly outside every cap: `|z|>5`, `|c|<7/3` and `R=1/3` give ball power greater than seven. Thus no finite-only walk would detect this additional edge family.

The complete 327-wall edge count is **524,608**: 320 finite-ended open edges and 524,288 seed-only vertical edges.

## One exact edge closure

Choose the positive target coordinate `g=11` and release cluster cap `r₀`, with unit-coordinate axis zero. Its nineteen retained incident walls are a,b and the other seventeen cluster caps. Set

```text
κ = L/(R²+18L),       t = 4,
A = 2t/(1−κt²) = (120824−1152φ)/1711,
B = (1+κt²)/(1−κt²) = (28495−288φ)/1711,
r_new = A P − B r₀.
```

Since `B²−κA²=1`, the new normal has square `R²`. It is right to all nineteen retained edge normals and strictly dotted to every other old wall. Its last coordinate `A−B` is positive. The [edge checker](../scripts/fresh_h20_central_ridge_20261014_edge.py) and [certificate](fresh_h20_central_ridge_20261014_edge.json) check all 327 pairings exactly: nineteen right and 308 dotted.

The new endpoint is

```text
Q = P − λr₀,       λ=Aκ/B,
Q² = −κR²/B² < 0.
```

All required incidences and every other old inequality are checked directly. The seventeen retained caps, the new cap and a,b again give `I₂(5)+18A₁`. All sixteen old corners survive strictly behind the new wall.

The new cap has exactly nineteen finite neighbors, so Q is its only possible finite vertex. Exactly nineteen new rank-nineteen subsets contain it. The final counts are therefore:

| Quantity | 328-wall chamber |
|---|---:|
| Actual facets | 328 |
| Finite vertices | 17 |
| Finite edges | 1 |
| Open finite-ended edges | 338 |
| Seed-only vertical edges | 524,288 |
| All hyperbolic edges | 524,627 |

The increase in open edges is explicit; this is a local finite-edge construction, not progress toward a verified bounded chamber.

## Remaining escape and stopping rule

Let x have all ten unit coordinates `99/100` and all heavy/golden coordinates zero. Then

```text
X = (x,8801/2000,1),       X²=1.
```

Every core inequality is strict. For all 288 cluster caps, the pairing is strictly less than `−5541/2000`; the edge certificate also verifies a negative pairing with the last cap. This future spacelike escape proves infinite volume at both stages.

There is also a structural limit to iterating this particular construction. Its induced central π/5 face remains right-angled. The [noncrystallographic face gate](fresh_noncrystallographic_face_gate_20261014.md) excludes finite-volume completion by continuing only within that right-angled face class. General nonright extensions or replacement of the golden ridge are outside that obstruction. No further local cap iteration is claimed.

Both saved replays passed:

```sh
python scripts/fresh_h20_central_ridge_20261014.py --verify
python scripts/fresh_h20_central_ridge_20261014_edge.py --verify
```
