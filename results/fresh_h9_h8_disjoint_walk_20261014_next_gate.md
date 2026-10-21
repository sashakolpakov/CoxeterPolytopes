# All current finite-endpoint extensions miss the quartet's escape

Every compatible finite-endpoint cap on the nineteen-facet chamber's 127 open finite-corner germs misses the same literal strict escape. This follows by combining the [original 62-germ theorem](fresh_h9_h8_disjoint_walk_20261014_universal_gate.md) with the new [71-germ certificate](fresh_h9_h8_disjoint_walk_20261014_next_gate.json). The original theorem covers the 56 old germs that remain open; the new theorem covers all 71 germs at the nine new finite vertices.

The [new checker](../scripts/fresh_h9_h8_disjoint_walk_20261014_next_gate.py) reconstructs the quartet's exact geometry in `Q(√5,√(6+16√5/5))`. Its outward interval embedding uses that radical explicitly. It verifies the normalized support Gram matrices, inverses, future Lorentz frames and all monotonicity hypotheses before the interval proof begins.

The new census has 1,065 finite spherical profiles and fifty isolated dihedral tails. Every real endpoint parameter is covered. Each unbounded dihedral family is enlarged to the continuous interval `√3≤2cos(π/m)≤2`. All nonsupport angle labels are retained by the same conservative pairing relaxation used in the original theorem. The new proof terminates in 3,269 boxes. The exact local profile census and the frozen interval engine are inherited; this is not a separate reimplementation of either.

The conclusion is a limit on the next finite-endpoint step. It does not exclude adding intermediate walls, reaching further supports after additional steps, adding walls without these finite endpoints, or changing the retained chamber. Bounded ideal endpoints on the new supports are also outside this certificate. In particular, this does not prove that no finite-volume completion exists.

```sh
python scripts/fresh_h9_h8_disjoint_walk_20261014_next_gate.py --verify
```
