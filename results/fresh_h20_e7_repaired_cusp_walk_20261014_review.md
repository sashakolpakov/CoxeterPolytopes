# Independent review of the repaired cusp's complete edge star

**PASS.** The [separate checker](../scripts/fresh_h20_e7_repaired_cusp_walk_20261014_review.py) and [certificate](fresh_h20_e7_repaired_cusp_walk_20261014_review.json) reproduce the complete classification of 1,728 edge germs in the [cusp census](fresh_h20_e7_repaired_cusp_walk_20261014.md): 26 finite endpoints, 16 ideal contact points with unbounded links, and 1,686 open opposite ends. Twenty finite vertices are new relative to the previous 96, giving at least 116 verified finite vertices.

The review imports no producer arithmetic. It uses the independent sparse squarefree field implementation, exact relative-norm comparisons for the bridge extension, and outward rational interval bounds.

It reconstructs all 25 component projections directly from the physical normals and rational inverse Gram matrices. The six components have positive affine kernels; every omitted-root submatrix is positive definite, and different components are orthogonal. Consequently choosing one omission per component gives exactly all 1,728 supports of rank 19. Every saved projection, squared norm and all 95 wall pairings agree exactly with the independent reconstruction.

All 26 finite and 16 ideal endpoint coordinates are checked against the full 95-wall chamber. Their norms, future orientations, full active sets, distinctness and incidence with the starting cusp agree. Each shares precisely the stated 19 independent walls with that cusp; every other inequality is strict at at least one endpoint. Thus the open segment is an actual edge with no intervening wall. Six finite endpoints match the previous list; twenty do not.

For an open germ, let `W` be the strict center minus the six chosen projections and let `ν` be the starting cusp. The other future null ray in the edge plane is

```text
Z = (−2〈H,ν〉)W + 〈W,W〉ν.
```

The review proves directly that all 76 inequalities outside the support are strict at `Z`. It uses 128-bit outward rational enclosures derived from the exact field coefficients: **128,136 strict signs**, with zero exact-sign fallbacks required. The support products vanish by the reconstructed projection identities. This proves the 1,686 open ends independently of the producer's threshold maximization. The numerical ordering used to propose a winning wall is not part of this independent argument; first-contact winners beyond the light cone are not separately re-maximized.

Each of the 16 other ideal contacts has a saved positive-norm recession vector. The review checks its norm, orthogonality to the null ray, and nonpositive products with every incident normal. Modulo the null direction it gives a nonzero Euclidean recession in the horospherical link, proving that link unbounded. Horospherical rank 19 alone would not establish boundedness.

Generation and saved replay pass. This is a complete star at one specified bounded cusp, not a complete vertex/edge census of the 95-wall chamber. The inherited strict escape and the explicitly open ends continue to prove infinite volume; no new finite-volume example is claimed.
