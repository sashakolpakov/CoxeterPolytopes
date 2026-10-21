# Review of the finite polar-bank obstruction

The frozen [primary verifier](../scripts/fresh_h9_makarov_polar_closure_20261013_verify.py) and its [certificate](fresh_h9_makarov_polar_closure_20261013_verify.json) pass a fresh saved-result replay. Source and proof review found no mathematical issue. This is a finite-bank obstruction: the 1,564 literal polars cannot complete a finite-volume Coxeter chamber, even when all 27 old walls are optional. It is not a classification of every possible polar or new wall.

The verifier reconstructs each rank-nine nullspace direction, checks its positive norm and orientation, and normalizes it to norm two. The 21 future spacelike witnesses satisfy every old inequality strictly. A possible cutter includes every uncertain or zero pairing, so interval uncertainty enlarges the candidate set. The pair gate excludes only certified positive pairings or strict gaps between allowed finite-angle intervals. Its retained tail includes every label above 120, parallel pairs, dotted pairs, and undecided comparisons. The resulting compatibility supergraph therefore contains every genuinely compatible selection.

The 31,308 exclusions and exhaustive 35-state cover search prove that no selection cuts all witnesses. Deleting redundant inequalities from any proposed polytope gives a subset of actual facets; their acute Coxeter pair relations satisfy this same relaxed graph. Such deletion cannot evade the cover obstruction. A surviving spacelike witness and a strict timelike interior point produce an open feasible region beyond the light cone, hence infinite volume.

A [separate standard-library checker](../scripts/fresh_h9_makarov_polar_closure_20261013_review.py) rechecks all 32,844 saved witness/polar interval decisions, all 31,308 saved conflict intervals, and the Boolean cover with an independently written symmetric-graph DFS. It also returns UNSAT, in 33 states. Its [result](fresh_h9_makarov_polar_closure_20261013_review.json) binds the frozen certificate. This separate check uses the saved interval data; it is not a second implementation of root reconstruction or rigorous interval arithmetic.

Both commands pass:

```sh
python3 scripts/fresh_h9_makarov_polar_closure_20261013_verify.py --verify
python3 scripts/fresh_h9_makarov_polar_closure_20261013_review.py --verify
```

Frozen primary source SHA-256: `1800b2183205cc7ab26a1c4f6546ac97bd612b2e91db22cc2402a1447e320a78`.
Frozen certificate SHA-256: `61d6714320cbff6d6d0cfd96f41deb1933d7e37f58776f833c7de51778e47317`.
Frozen recipes SHA-256: `c6cc16781bac8272205f3a5146b339f2568fb5453de9789ba3016c05256c8cbb`.
