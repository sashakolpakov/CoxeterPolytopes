# Independent review of the transported H22 cusp repair

**PASS, with an infinite-volume conclusion.** I source-reviewed the orbit construction and exact pair gates, freshly reran the complete primary verifier, and independently reconstructed the representative link using a separate standard-library Fraction implementation. The completed configuration has 11,688,030 actual facets and483,840 specified bounded cusp links, but an exact strict spacelike point proves infinite volume. A further local argument identifies an uncapped ideal direction in every newly added facet.

## What was independently reconstructed

The [separate checker](../scripts/fresh_h22_transported_cusp_completion_20261013_local_independent.py) imports no producer, project arithmetic, NumPy or SymPy. It reads the explicit rational proposed normals and raw original form, representing the added coordinate as its coefficient of sqrt(11). It verifies:

* ambient signature(22,1), strict center square−14, all 13 new norms and all 2,730 pairings with the original210 walls;
* all 465 pairings among the 31 incident normals:440 right,19 pi/3, and6 parallel;
* exactly the prescribed17 incident original walls, the existing repair wall, and13 new walls;
* normal-span rank 22, Gram inertia(21,0,10), and the full-support positive relation`sum(r_i)=66*delta`.

The positive relation makes every recession direction orthogonal to every quotient normal. Their rank 21 spans the full21-dimensional Euclidean quotient, so the link is bounded. This checks the positive-dependence condition that a rank calculation alone would miss.

The [local certificate](fresh_h22_transported_cusp_completion_20261013_local_independent.json) and its saved replay both passed.

## Review of the global transfer

The primary verifier checks the explicit241,920-element orbit arrays under every original symmetry generator, together with a depth-decreasing tree reaching every frame from the representative. Each generator preserves the original form, strict center, and original wall set. The same generators preserve the inherited complete banks, including the four replacement incidence walls and the 560-wall binary family. Consequently representative pair gates transfer to every transported frame.

The height signs are handled correctly: all new positive-height coordinates are positive; the opposite-sign pairing is smaller than the same-sign pairing. Nonlocal new–new gates explicitly check both signs. Old radical families use exact comparisons of`A+B*sqrt(d)` against the dotted threshold. Native contractions have stated integer overflow bounds, and the final inequalities use Python integers or Fractions. The two original center walls are intentionally replaced by the four incidence walls; they are not missing retained constraints.

Every nonlocal pair is strictly dotted. Within one cusp the local gates are right or parallel for new–new relations. The strict timelike center and nonpositive pairings prove each inequality is an actual facet by orthogonal projection; positive-proportional duplicate normals would violate the pair gates. The specified bounded-link count is not a census of all ideal ends or a finite-volume criterion.

Fresh primary replay passed:

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h22_transported_cusp_completion_20261013_verify.py --verify
```

## Exact global escape and structural obstruction

The [escape verifier](../scripts/fresh_h22_transported_cusp_completion_20261013_escape.py) replays the previous5,398,110 inequalities and checks all 6,289,920 new inequalities through13 exact full-family maxima. For each family it records integers`A<0`, `B>=0`, and`A^2>11*B^2`, proving both height signs strictly negative. Its rational point has square59356314529/250000000000>0 and future orientation. Its saved replay passed. This full-array check reuses the inherited exact checker; it is not described as a second independent full-orbit enumeration.

The separate [facet diagnostic](../scripts/fresh_h22_transported_cusp_completion_20261013_facet_obstruction.py) constructs an explicit future null direction in each representative new facet using only Fraction arithmetic. Every other local inequality is strict there. The primary dotted-pair certificate implies all nonlocal projected normals are future timelike, so they too are strict. A small motion beyond that null direction yields a spacelike point, and a further motion toward the strict center makes every wall strict.

This reveals a geometric limitation: all hyperplanes meeting one of the new facets pass through its same cusp delta. No finite vertex lies on that facet. A finite acute-compatible completion retaining the facet must eventually add a wall that meets it at finite hyperbolic points **without containingdelta**, or must replace the facet. Parallel and dotted hyperplanes cannot bound its intrinsic hyperbolic interior. This is a condition on the incidence pattern, not a global H22 nonexistence theorem or a prohibition on arbitrary new wall families.

[Combined review script](../scripts/fresh_h22_transported_cusp_completion_20261013_review.py) · [review certificate](fresh_h22_transported_cusp_completion_20261013_review.json) · [facet directions](fresh_h22_transported_cusp_completion_20261013_facet_obstruction.json) · [exact escape](fresh_h22_transported_cusp_completion_20261013_escape.json).
